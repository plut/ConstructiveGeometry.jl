module ConstructiveGeometry
# Module imports««1
# using Printf
using LinearAlgebra
using StaticArrays
using FixedPointNumbers
using Logging
using FastClosures
# using DataStructures

import Polyhedra # for convex hull
import GLPK

# using MiniQhull
import Rotations
import Colors: Colors, Colorant
import Clipper


import Base: show, print
import Base: union, intersect, setdiff, copy, isempty, merge
import Base: *, +, -, ∈, inv, sign, iszero

include("Shapes.jl")
using .Shapes
# include("SpatialSorting.jl")
# include("TriangleIntersections.jl")
include("CornerTables.jl")
using .CornerTables
include("scad.jl")

# Types««1
# Angle types««2
# to keep it simple, angles are just Float64 (in degrees).
# Just in case we might change this in the future, we define two types:
# Angle and AnyAngle (= input type).
const Angle = Float64
const AnyAngle = Real
const ° = 1.
@inline radians(x::Angle) = π/180*x
@inline radians(x::AnyAngle) = radians(Angle(x))
@inline degrees(x::Angle) = x
@inline degrees(x::AnyAngle) = degrees(Angle(x))

# General tools««1
norm²(v) = sum(v .* v)
distance²(p, q) = norm²(p-q)

struct Consecutives{T,V} <: AbstractVector{T}
	parent::V
end
consecutives(v::AbstractVector{T}) where{T} =
	Consecutives{T,typeof(v)}(v)
Base.getindex(c::Consecutives, i::Integer) =
	(c.parent[i], c.parent[mod1(i+1, length(c.parent))])
Base.size(c::Consecutives) = size(c.parent)
# uniquereindex #««
# """"
#      uniquereindex(a)
# 
# Returns `(new_a, idx)`, where
#  - `new_a` is `unique(sorted(a))`;
#  - `new_a[idx]` is same as `a`.
# """
# function uniquereindex(a)
# 	perm=sortperm(a)
# 	idx = similar(perm)
# 	@inbounds t = perm[1]
# 	@inbounds idx[t] = 1
# 	@inbounds l = a[t]
# 	new_a = sizehint!([l], length(a))
# 	i = 1
# 	@inbounds while i < length(perm)
# 		i+= 1
# 		t = perm[i]
# 		if a[t] ≠ l
# 			l = a[t]
# 			push!(new_a, l)
# 		end
# 		idx[t] = length(new_a)
# 	end
# 	return (new_a, idx)
# end
# #»»
# small determinants««2

# 2-dimensional determinant, useful for computing orientation
# @inline det2(v1, v2) = v1[1]*v2[2] - v1[2]*v2[1]
@inline det2(pt1, pt2, pt3) = det2(pt2-pt1, pt3-pt1)

# 3-dimensional determinant
# @inline det3(v1::Vec{3}, v2::Vec{3}, v3::Vec{3}) = det([v1 v2 v3])
# @inline det3(p1::Point{3}, p2::Point{3}, p3::Point{3}, p4::Point{3}) =
# 	det3(p2-p1, p3-p1, p4-p1)
# @inline det3(v1, v2, v3) = det([v1 v2 v3])
# @inline det3(p1, p2, p3, p4) = det3(p2-p1, p3-p1, p4-p1)
# Rows view««2
struct RowsView{T,M<:AbstractMatrix{T}} <:
		AbstractVector{SubArray{T,1,M,Tuple{T,Base.Slice{Base.OneTo{T}}},true}}
	source::M
	RowsView(m::AbstractMatrix) = new{eltype(m), typeof(m)}(m)
end
Base.size(r::RowsView) = (size(r.source,1),)
Base.getindex(r::RowsView, i::Integer) = view(r.source, i, :)

# Reflection matrices««2
"""
		Reflection

A type containing compressed information for an orthogonal reflection
matrix. Inspired by the `Diagonal` type.
"""
struct Reflection{T,V<:AbstractVector} <: AbstractMatrix{T}
# The reflection formula is r(x) = x - 2 (x⋅v/‖v‖²) v
# hence we store v and (v/‖v‖²) separately (this allows exact
# computations where needed):
	axis::V # v
	proj::Transpose{T,V} # v/‖v‖²
end
StaticArrays.similar_type(::Vector, T::Type) = Vector{T} # piracy!
@inline Reflection(v::AbstractVector) =
	Reflection{typeof(v[1]/1)}(v)
@inline Reflection{T}(v::AbstractVector) where{T} =
	Reflection{T, similar_type(v, T)}(v, transpose(v*inv(dot(v,v))))
# normed case:
@inline Reflection(k::Val{:normed}, v::AbstractVector) =
	Reflection{eltype(v)}(k, v)
@inline Reflection{T}(::Val{:normed}, v::AbstractVector) where{T} =
	Reflection{T, typeof(v)}(v, transpose(v))

# allow this to behave as a StaticMatrix when needed:
@inline Size(::Type{Reflection{T,V}}) where{T,V} = (Size(V,1), Size(V,1))
@inline Size(r::Reflection) = Size(typeof(r))
@inline Base.size(r::Reflection) = (size(r.axis,1), size(r.axis,1),)
@inline det(r::Reflection) = -one(eltype(r))
@inline tr(r::Reflection) = size(r,1)-2*one(eltype(r))

# @inline getindex(r::Reflection, i::Int) =
# 	getindex(R, fld1(i,D), mod1(i,D))
@inline Base.getindex(r::Reflection, i::Int, j::Int) =
	(i == j) - 2*r.axis[i]*r.proj[j]
# we need two methods here for disambiguation:
@inline Base.:*(r::Reflection, v::AbstractVector) = v - 2*r.axis*(r.proj*v)
@inline Base.:*(r::Reflection, v::AbstractMatrix) = v - 2*r.axis*(r.proj*v)

# Affine transformations««2
abstract type AbstractAffineMap end
struct AffineMap{A,B} <: AbstractAffineMap
	a::A
	b::B
end
struct LinearMap{A} <: AbstractAffineMap
	a::A
end
struct TranslationMap{B} <: AbstractAffineMap
	b::B
end
@inline AffineMap(a; center=Val(false)) = _affine_map_center(a, center)
@inline _affine_map_center(a, ::Val{false}) = LinearMap(a)
@inline _affine_map_center(a, c) = AffineMap(a, a*c-c)
@inline linear(f::AbstractAffineMap, x) = f.a*x
@inline linear(f::TranslationMap, x) = x
@inline translate(f::AbstractAffineMap, x) = x+f.b
@inline translate(f::LinearMap, x) = x
@inline affine(f::AbstractAffineMap, x) = translate(f, linear(f, x))
@inline (f::AbstractAffineMap)(x) = affine(f,x)

@inline signdet(a::Number) = sign(a)
@inline signdet(a::AbstractMatrix) = sign(det(a))
@inline signdet(f::AbstractAffineMap) = signdet(f.a)
@inline signdet(f::TranslationMap) = +1

# @inline (f::Affine)(p::Point) = Point(f.a * coordinates(p) + f.b)
# function (f::Affine)(p::Polygon)
# 	if sign(f) > 0
# 		return Polygon(f.(vertices(p)))
# 	else
# 		return Polygon(reverse(f.(vertices(p))))
# 	end
# end
# (f::Affine)(p::nolygonXor) = PolygonXor(f.(paths(p)))

# 	[apply(f, p) for p in points]

# @inline Base.sign(f::Affine{<:Number}) = sign(f.a)
# @inline Base.sign(f::Affine{<:AbstractMatrix}) = sign(det(f.a))

# # I/O: ««3
# # OpenSCAD only uses 4×4 matrices for transformations;
# # we pad the matrix to this size if needed:
# function scad_parameters(io::IO, f::Affine)
# 	m = [ mat33(f.a) vec3(f.b); 0 0 0 1 ]
# 	print(io, "[")
# 	join(io, map(i->Float64.(view(m,i,:)),1:size(m,1)),",")
# 	print(io, "]")
# end
# 
# @inline mat33(a::AbstractMatrix) = [ get(a, (i,j), i==j) for i=1:3, j=1:3 ]
# @inline mat33(a::Diagonal) = Diagonal(vec3(diag(a)))
# @inline mat33(a::Real) = SDiagonal(a,a,a)
# @inline mat33(::Val{b}) where{b} = mat33(b)
# 
# @inline vec3(b::AbstractVector) = SVector{3}(get(b, i, 0) for i in 1:3)
# @inline vec3(::Val{b}) where{b} = SA[b, b, b]


#————————————————————— Tools —————————————————————————————— ««1

#»»1
# Parameters««1
# Definitions««2
# Accuracy is the absolute deviation allowed.
# Default value is 2.0 (from OpenSCAD `$fs`), interpreted as 2mm.
#
# Precision is the relative deviation allowed.
# Default value is 0.02 (1-cos(180°/`$fa`)).
# FIXME: explain why .005 works better

_DEFAULT_PARAMETERS = (accuracy = 0.1, precision = .005, symmetry = 1,
	type = Float64, ε = 1e-8)

@inline get_parameter(parameters, name) =
	get(parameters, name, _DEFAULT_PARAMETERS[name])

# Circles««2

round(m::Integer, x::Integer, ::typeof(RoundUp)) = x + m - mod1(x, m)
"""
    sides(radius, parameters)

Returns the number of sides used to draw a circle (arc) of given angle.
The base value `n` is given by the minimum of:
 - accuracy: each sagitta (s= r(1-cos 2π/n)) must not be smaller
 than `accuracy`, or n = π √(r/2 accuracy);
 - precision: s/r = 1-cos(2π/n)  not smaller than precision,
 or n = π /√(2*precision).
"""
function sides(r::Real, parameters)
  ε = max(get_parameter(parameters,:precision),
		get_parameter(parameters,:accuracy)/r)
	base = ceil(Int, π/√(2*ε))
	# a circle always has at least 4 sides
	return round(get_parameter(parameters,:symmetry), max(4, base), RoundUp)
end


"""
    unit_n_gon(T::Type, n::Int)
		unit_n_gon(r, parameters::NamedTuple)

Returns the vertices of a regular n-gon inscribed in the unit circle
as points with coordinates of type `T`, while avoiding using too many
trigonometric computations.
"""
function unit_n_gon(T::Type{<:Real}, n::Int)
	ω = cis(2π/n) # exp(2iπ/n)
	z = Vector{Complex{T}}(undef, n)
	z[1] = one(T)
	# TODO: use 2-fold, 4-fold symmetry if present
	# n=3: 2..2
	# n=4: 2..2 (+ point 3 = -1)
	# n=5: 2..3
	# n=6: 2..3 (+ point 4 = -1)
	# points with y>0:
	for i in 2:(n+1)>>1
		@inbounds z[i] = z[i-1]*ω; z[i] /= abs(z[i]) # for radius stability
		# z[n] = conj(z[2]), etc.
		@inbounds z[n+2-i] = conj(z[i])
	end
	if iseven(n)
		@inbounds z[n>>1+1] = -1
	end
	reinterpret(SVector{2,T}, z)
end
@inline unit_n_gon(r, parameters::NamedTuple) =
	r*unit_n_gon(get_parameter(parameters, :type), sides(r, parameters))

# Spheres««2
"""
    sphere_vertices(r::Real, parameters::NamedTuple)

Returns the number `n` of points on a sphere according to these
parameters.

This produces n points on the sphere, hence 2n-4 triangular faces
(genus 0). Average area of a triangular face is 4πr²/(2n-4)=2πr²/(n-2),
hence square of edge length is d²≈ (8π/√3) r²/(n-2).
(unit equilateral triangle area: A=√3d²/4, i.e. d²=(4/√3)A).

Sagitta is given by
s/r = 1-√{1-d²/4r²}
≈ 1-(1-d²/8 r²)
≈ 1-(1-(π/√3)/(n-2))
≈ (π/√3)/(n-2).
Hence n ≈ 2 + (π/√3)/(precision).

"""
function sphere_vertices(r::Real, parameters::NamedTuple = _DEFAULT_PARAMETERS)
  ε = max(get_parameter(parameters,:precision), get_parameter(parameters,:accuracy)/r)
	base = 2 + ceil(Int, (π/√3)/ε)
	# a sphere always has at least 6 vertices
	return max(6, base)
end

const golden_angle = 2π/MathConstants.φ
"""
    fibonacci_sphere_points(T::Type{<:Real}, n::Int)

Returns a set of `n` well-spaced points, of type `SVector{3,T}`, on the unit
sphere.

TODO: use ideas from
http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/

to optimize for volume of convex hull.
"""
function fibonacci_sphere_points(T::Type{<:Real}, n::Int)

	v = Vector{SVector{3,T}}(undef, n)
	for i in eachindex(v)
		θ = i*T(golden_angle)
		z = (n+1-2i)/T(n)
		ρ = √(1-z^2)
		(s,c) = sincos(θ)
		@inbounds v[i] = SA[c*ρ, s*ρ, z]
	end
	return v
end
@inline fibonacci_sphere_points(r::Real, parameters...) =
	r*fibonacci_sphere_points(real_type(r),
		sphere_vertices(r, parameters...))

#————————————————————— Definitions of objects —————————————————————————————— ««1

#»»1
# Abstract type««1
"""
    AbstractGeometry{D}

A `D`-dimensional geometric object.
Interface:
 - `children`
 - `coordtype`: coordinate type (if applicable).
"""

abstract type AbstractGeometry{D} end
@inline Base.ndims(::AbstractGeometry{D}) where{D} = D
@inline children(::AbstractGeometry) = AbstractGeometry[]

abstract type AbstractGeometryCoord{D,T} <: AbstractGeometry{D} end
@inline coordtype(::AbstractGeometryCoord{D,T}) where{D,T} = T

# 2d primitives««1

# Ortho: square or cube (orthotope) with a corner at zero««2
"""
    Ortho{D,T}

An orthotope (Cartesian product of intervals), with lower corner at zero,
dimension `D`, and coordinate type `T`.
"""
struct Ortho{D,T} <: AbstractGeometryCoord{D,T}
	size::SVector{D,T}
end

@inline Ortho{D}(a::AbstractVector) where{D} = Ortho{D,eltype(a)}(a)
@inline Ortho{D}(a::Real...) where{D} = Ortho{D}(SVector(a))
@inline Ortho(a::Real...) = Ortho(SVector{length(a)}(a...))
Square = Ortho{2}

@inline scad_info(s::Square) = (:square, (size=s.size,))

function vertices(s::Square)
	# in trigonometric order:
	z = zero(coordtype(s)); v = s.size
	return [
		SA[z   ,z   ],
		SA[v[1],z   ],
		SA[v[1],v[2]],
		SA[z   ,v[2]]]
end
mesh(s::Square, parameters) =
	PolygonXor{get_parameter(parameters, :type)}(vertices(s))

# Ball: circle or sphere, centered at zero««2
"""
    Ball{D,T}

A `D`-dimensional ball, with radius of type `T`, centered at the origin.
"""
struct Ball{D,T} <: AbstractGeometryCoord{D,T}
	radius::T
end
@inline Ball{D}(a::Real) where{D} = Ball{D,typeof(a)}(a)

Circle = Ball{2}
@inline scad_info(s::Circle) = (:circle, (r=s.radius,))

# Polygon««2
"""
    Polygon{T}

A simple, closed polygon enclosed by the given vertices.
`T` is the coordinate type.
"""
struct Polygon{T} <: AbstractGeometryCoord{2,T}
	points::Vector{SVector{2,T}}
end
@inline Polygon(points::AbstractVector{<:AbstractVector{<:Real}}) =
	Polygon{eltype(eltype(points))}(points)

@inline scad_info(s::Polygon) = (:polygon, (points=s.points,))
@inline mesh(s::Polygon, parameters) =
	PolygonXor{get_parameter(parameters, :type)}(s.points)

# Draw ««2
"""
    draw(path, width; kwargs...)
    ends=:round|:square|:butt|:loop
    join=:round|:miter|:square
"""
@inline draw(args...; kwargs...) = Draw(args...; kwargs...)
struct Draw{T} <: AbstractGeometryCoord{2,T}
	path::Vector{SVector{2,T}}
	width::Float64
  ends::Symbol
	join::Symbol
	miter_limit::Float64
end
Draw(path, width; ends=:round, join=:round, miter_limit=2.) =
	Draw{coordtype(path)}(path, width, ends, join, miter_limit)

"""
    draw(path, width; kwargs)
    ends = :loop|:butt|:square|:round
		join = :square|:round|:miter
    miter_limit = 2.0

Draws a path of given width.
"""
draw(path, width; kwargs...) = Draw(path, width; kwargs...)

# 3d primitives««1

Cube = Ortho{3}
@inline scad_info(s::Cube) = (:cube, (size=s.size,))
function mesh(s::Cube, parameters)
	pts = vertices(s, parameters)
	return CornerTable(pts, [ # 12 triangular faces:
	 (6, 5, 7), (7, 8, 6), (7, 3, 4), (4, 8, 7),
	 (4, 2, 6), (6, 8, 4), (5, 1, 3), (3, 7, 5),
	 (2, 1, 5), (5, 6, 2), (3, 1, 2), (2, 4, 3),
	])
end

Sphere = Ball{3}
@inline scad_info(s::Sphere) = (:sphere, (r=s.radius,))
# FIXME: mesh = convex hull

# Cylinder is an extrusion

# Surface««2
"""
    Surface{T}

A surface delimited by the given faces.
"""
struct Surface{T} <: AbstractGeometryCoord{3,T}
	points::Vector{SVector{3,T}}
	faces::Vector{NTuple{3,Int}}
end
@inline Surface(points, faces) = Surface{eltype(eltype(points))}(points, faces)
@inline Surface(points, faces::AbstractVector{<:AbstractVector}) =
# FIXME: in this case, triangulate!
	Surface(points, [(f...,) for f in faces])

@inline scad_info(s::Surface) =
	(:surface, (points=s.points, faces = [ f .- 1 for f in s.faces ]))
@inline mesh(s::Surface, parameters) = CornerTable(s.points, s.faces)


# Constructive geometry operations««1
# https://www.usenix.org/legacy/event/usenix05/tech/freenix/full_papers/kirsch/kirsch.pdf
# Type definition««2
"""
		ConstructedSolid{D,S}

A type representing CSG operations on solids. `D` is the dimension and
`S` is a symbol representing the operation (union, intersection etc.)
"""
struct ConstructedSolid{S,V,D} <: AbstractGeometry{D}
	children::V # Vector{<:AbstractGeometry}, or tuple etc.
	# passing a vector or tuple:
	@inline ConstructedSolid{S,V,D}(v::V) where{S,V,D} = new{S,V,D}(v)
end

@inline children(s::ConstructedSolid) = s.children
@inline scad_info(s::ConstructedSolid{S}) where{S} = (S, ())

constructed_solid_type(s::Symbol, T = Vector{<:AbstractGeometry}) =
	ConstructedSolid{s,T}

# Generic code for associativity etc.««2
# make operators associative; see definition of + in operators.jl
for op in (:union, :intersect, :minkowski, :hull)
	Q=QuoteNode(op)
	# union, intersection, minkowski are trivial on single objects:
	op != :hull &&  @eval ($op)(a::AbstractGeometry) = a
	@eval begin
	# all of these are associative:
	# we leave out the binary case, which will be defined on a case-by-case
	# basis depending on the operators (see below).
#		($op)(a::Geometry, b::Geometry) =
#			ConstructedSolid{$Q}([unroll(a, Val($Q)); unroll(b, Val($Q))])
	($op)(a::AbstractGeometry, b::AbstractGeometry, c::AbstractGeometry, x...) =
		Base.afoldl($op, ($op)(($op)(a,b),c), x...)
	end
end

"""
		unroll(x::AbstractGeometry, Val(sym1), Val(sym2)...)

Returns either `[x]` or, if `x` is a `ConstructedSolid` matching one of the
symbols `sym1`, `sym2`..., `children(x)`. (This helps reduce nesting).
"""
@inline unroll(s::AbstractGeometry, ::Val, tail...) = unroll(s, tail...)
@inline unroll(s::AbstractGeometry) = s
@inline unroll(s::ConstructedSolid{D, S}, ::Val{S}, tail...) where{D, S} =
	children(s)
@inline unroll2(s::AbstractGeometry, t::AbstractGeometry, tail...) =
	[unroll(s, tail...); unroll(t, tail...)]

# Union««2
CSGUnion = constructed_solid_type(:union)
# this calls union of CornerTable or PolygonXor as needed:
@inline mesh(s::CSGUnion{2}, parameters) =
	union(mesh(x, parameters) for x in children(s))
@inline mesh(s::CSGUnion{3}, parameters) =
	CornerTables.combine([mesh(x, parameters) for x in children(s)], 1,
		parameters.ε)
@inline union(a1::AbstractGeometry{D}, a2::AbstractGeometry{D}) where{D} =
	CSGUnion{D}(unroll2(a1, a2, Val(:union)))

# Intersection««2
CSGInter = constructed_solid_type(:intersection)
@inline mesh(s::CSGInter, parameters) =
	intersect(mesh(x, parameters) for x in children(s))
# FIXME allow intersection of 2d and 3d primitives (as 2d result)
@inline intersect(a1::AbstractGeometry, a2::AbstractGeometry) =
	CSGInter{minimum(embeddim.((a1,a2)))}(unroll2(a1, a2, Val(:intersection)))

# Difference««2
# this is a binary operator:
CSGDiff = constructed_solid_type(:difference,
	Tuple{<:AbstractGeometry,<:AbstractGeometry})
@inline mesh(s::CSGDiff, parameters) =
	setdiff(mesh(s.children[1], parameters), mesh(s.children[2], parameters))
@inline setdiff(x::AbstractGeometry{D}, y::AbstractGeometry{D}) where{D} =
	CSGDiff{D}((x,y))
# added interface: setdiff([x...], [y...])
@inline setdiff(x::AbstractVector{<:AbstractGeometry},
                y::AbstractVector{<:AbstractGeometry}) =
	difference(union(x...), union(y...))
Base.:\(x::AbstractGeometry, y::AbstractGeometry) = setdiff(x, y)

# Complement««2
# TODO
CSGComplement = constructed_solid_type(:complement, Tuple{<:AbstractGeometry})

# Empty unions and intersects««2
"""
    EmptyUnion

A convenience type representing the union of nothing.
This is removed whenever it is `union()`-ed with anything else.
"""
struct EmptyUnion end
"""
    EmptyIntersect

A convenience type representing the intersection of nothing.
This is removed whenever it is `intersect()`-ed with anything else.
"""
struct EmptyIntersect end

union() = EmptyUnion()
intersect() = EmptyIntersect()
Base.show(io::IO, ::EmptyUnion) = print(io, "union()")
Base.show(io::IO, ::EmptyIntersect) = print(io, "intersect())")

macro define_neutral(op, what, result)
	quote
	@inline $(esc(op))(neutral, absorb::$what) = $result
	@inline $(esc(op))(absorb::$what, neutral) = $result
	@inline $(esc(op))(x::$what, ::$what) = x
	end
end
function minkowski end
function hull end
@define_neutral union EmptyUnion neutral
@define_neutral union EmptyIntersect  absorb
@define_neutral intersect EmptyUnion absorb
@define_neutral intersect EmptyIntersect  neutral
@define_neutral minkowski EmptyUnion absorb
@define_neutral minkowski EmptyIntersect  absorb
@define_neutral hull EmptyUnion neutral
@define_neutral hull EmptyIntersect  absorb

# Convex hull««2
# TODO
CSGHull = constructed_solid_type(:hull)
"""
    hull(s::AbstractGeometry...)

Represents the convex hull of given solids.
"""
@inline hull(s::AbstractGeometry...) =
	CSGHull{maximum(embeddim.(s))}(
		[unroll(t, Val.((:hull, :union))...) for t in s])


# Minkowski sum and difference««2
# TODO
CSGMinkowski = constructed_solid_type(:minkowski)
"""
    minkowski(s::AbstractGeometry...)

Represents the Minkowski sum of given solids.
"""
@inline minkowski(a1::AbstractGeometry, a2::AbstractGeometry) =
	CSGMinkowski{maximum(embeddim.((a1,a2)))}(unroll2(a1, a2, Val(:minkowski)))

# Transformations««1
# AbstractTransform type««2
"""
    AbstractTransform{S,D}

Represents a solid of dimension `D` obtained via a transformation with
name `S` (a symbol).

This type defines functions allowing to chain transforms; these are used by
`multmatrix`, `color` etc. operations (see below).

The minimal job left to concrete types (see e.g. `AffineTransform` as an
example) is to define a type and a constructor:
    Frobnicate = Transform{:frobnicate}
		frobnicate(x::real, s...) = Frobnicate((x=x,), s...)
"""
abstract type AbstractTransform{D} <: AbstractGeometry{D} end

# AffineTransform««2
struct AffineTransform{D,A<:AbstractAffineMap} <: AbstractTransform{D}
	f::A
	child::AbstractGeometry{D}
end
# AffineTransform(f, child) constructor is defined

# FIXME
function mesh(s::AffineTransform{3}, parameters)
	m = mesh(s.child, parameters)
	m.points .= s.f.(m.points)
	return m
end
@inline scad_info(s::AffineTransform) = (:multmatrix, (m=mat44(s.f),))

"""
    mult_matrix(a, [center=c], solid...)
    mult_matrix(a, b, solid...)
    mult_matrix(a, b) * solid

Represents the affine operation `x -> a*x + b`.

# Extended help
!!! note "Types of `mult_matrix` parameters"

    The precise type of parameters `a` and `b` is not specified.
    Usually, `a` will be a matrix and `b` a vector, but this is left open
    on purpose; for instance, `a` can be a scalar (for a scaling)
    and `b` can be `Val(false)` for a linear operation. Any types so that
    `a * Vector + b` is defined will be accepted.

    Conversion to a matrix will be done when converting to OpenSCAD
    format.

!!! note "Matrix multiplication"

    Chained `mult_matrix` operations will be combined into a single
    operation when possible. This saves time: multiple
    (3 × n) matrix multiplications are replaced by
    (3 × 3) multiplications, followed by a single (3 × n).
"""
@inline mult_matrix(a, s...; kwargs...) =
	AffineTransform(Affine(a; kwargs...), s...)
@inline mult_matrix(a, b::Union{AbstractVector,Val{false}}, s...) =
	AffineTransform(Affine(a, b), s...)
@inline parameters(s::AffineTransform) = (m=s.data,)

# # these two functions are now enough to pre-compose all affine transforms
# # *before* applying them to objects:
# @inline assoc(::Curry{:multmatrix}, ::Curry{:multmatrix}) = :left
# @inline function compose(c1::Curry{:multmatrix}, c2::Curry{:multmatrix})
# 	(f1, f2) = (extract(c1), extract(c2))
# 	mult_matrix(f1.a*f2.a, f1.a*f2.b + f1.b)
# end
# 
# # Translation, scaling, rotation, mirror««2
# # FIXME change this '1' to a compile-time constant?
# """
#     translate(v, s...)
#     translate(v) * s
# 
# Translates solids `s...` by vector `v`.
# """
# @inline translate(v::AbstractVector, s...) = mult_matrix(1, v, s...)
# """
#     scale(a, s...; center=0)
#     scale(a; center=0) * s
# Scales solids `s` by factor `a`. If `center` is given then this will be
# the invariant point.
# 
# `a` may also be a vector, in which case coordinates will be multiplied by
# the associated diagonal matrix.
# """
# @inline scale(a::Real, s...; kwargs...) = mult_matrix(a, s...; kwargs...)
# @inline scale(a::AbstractVector, s...; kwargs...) =
# 	mult_matrix(Diagonal(a), s...; kwargs...)
# """
#     mirror(v, s...; center=0)
#     mirror(v; center=0) * s
# 
# Reflection with axis given by the hyperplane normal to `v`.
# If `center` is given, then the affine hyperplane through this point will
# be used.
# """
# @inline mirror(v::AbstractVector, s...; kwargs...) =
# 	mult_matrix(Reflection(v), s...; kwargs...)
# 
# @inline rotation(θ::AnyAngle; axis=SA[0,0,1], kwargs...) =
# 	real_type(θ,axis...).(Rotations.AngleAxis(radians(θ), axis...))
# @inline rotation(θ::AbstractVector{<:AnyAngle}; kwargs...) =
# 	real_type(θ,axis...).(Rotations.RotZYX(radians.(θ); kwargs...))
# 
# """
#     rotate(θ, {center=center}, {solid...})
#     rotate(θ, axis=axis, {center=center}, {solid...})
# 
# Rotation around the Z-axis (in trigonometric direction, i.e.
# counter-clockwise).
# """
# @inline rotate(θ, s...; kwargs...) = mult_matrix(rotation(θ; kwargs...), s...)
# """
#     rotate((θ,φ,ψ), {center=center}, {solid...})
# 
# Rotation given by Euler angles (ZYX; same ordering as OpenSCAD).
# """
# @inline rotate(θ::Real, φ::Real, ψ::Real, s...; kwargs...) =
# 	mult_matrix(rotation((θ,φ,ψ); kwargs...), s...)
# 
# 
# struct Transform{S,D,D1,X} <: AbstractGeometry{D}
# 	data::X
# 	child::AbstractGeometry{D1}
# 	Transform{S,D,D1}(data, child::AbstractGeometry{D1}) where{S,D,D1} =
# 		new{S,D,D1,typeof(data)}(data, child)
# 	Transform{S,D}(data, child::AbstractGeometry) where{S,D} =
# 		Transform{S,D,embeddim(child)}(data, child)
# 	# Default case: D = D1
# 	Transform{S}(data, child::AbstractGeometry) where{S} =
# 		Transform{S,embeddim(child)}(data, child)
# end
# # more constructors, including unary curryfied constructor:
# @inline (T::Type{<:Transform{S}})(f, s1::AbstractGeometry,
# 		s2::AbstractGeometry, tail::AbstractGeometry...) where{S} =
# 		T(f, union(s, s2, tail...))
# @inline (T::Type{<:Transform{S}})(f, s::Vector{<:AbstractGeometry}) where{S} =
# 	T(f, s...)
# @inline (T::Type{<:Transform{S}})(f) where{S} = Curry{S}((s...)->T(f, s...))
# @inline (T::Type{<:Transform})(f, ::typeof(extract)) = f
# 
# # default values for I/O:
# # (parameters in `data` are assumed to be stored in a NamedTuple).
# # @inline children(f::Transform) = [f.child]
# # @inline scad_name(f::Transform{S}) where{S} = S
# # @inline parameters(f::Transform) = f.data
# 
# # Curry««2
# """
#     Curry{S}
# 
# A structure representing partially-evaluated functions.
# This allows chaining transformations by overloading the multiplication
# operator: each factor in such a 'product', except the last one,
# is a `Curry` object.
# 
# `S` is a datum indicating the type of transformation performed by the
# function. It is used to compose functions when possible.
# 
# # Examples
# ```jldoctest
# julia> add(a)=Curry(x->x+a)
# julia> add(1)*add(2)*4
# 7
# ```
# """
# struct Curry{S}
#   f # either a Function or a Type...
# end
# # poor man's associative functor...
# 
# # fall-back case:
# @inline Base.:*(f::Curry) = f
# # binary rules:
# @inline Base.:*(f::Curry, g::Curry) = compose(f, g)
# @inline Base.:*(f::Curry, x) = f.f(x)
# # ternary rule for associativity: we use the `assoc` type trait to
# # decide whether to associate left or right.
# @inline Base.:*(f::Curry, g::Curry, args...) =
# 	_comp(Val(assoc(f,g)), f, g, args...)
# @inline _comp(::Val{:left} , f, g, args...) = *(compose(f, g), args...)
# @inline _comp(::Val{:right}, f, g, args...) = *(f, *(g, args...))
# 
# # default values for the traits: transforms are right-associative and
# # composition is trivial.
# @inline assoc(::Curry, ::Curry) = :right
# @inline compose(f::Curry, g::Curry) = Curry{:∘}(f.f ∘ g.f)
# 
# # We can extract the `f` value from the above in the following way:
# """
#     extract(c::Curry)
# 
# Given a `Curry` object with function `s -> Transform{...}(f, s)`,
# recovers the parameter `f`.
# """
# function extract end
# @inline extract(c::Curry) = c.f(extract)
# 
# # SetParameters««2
# SetParameters = Transform{:parameters}
# """
#     set_parameters(;accuracy, precision, symmetry) * solid...
# 
# A transformation which passes down the specified parameter values to its
# child. Roughly similar to setting `\$fs` and `\$fa` in OpenSCAD.
# """
# @inline set_parameters(s...; parameters...) =
# 	SetParameters(parameters.data, s...)
# 
# # Color««2
# Color = Transform{:color}
# 
# """
#     color(c::Colorant, s...)
#     color(c::AbstractString, s...)
#     color(c::AbstractString, α::Real, s...)
#     color(c) * s...
# 
# Colors objects `s...` in the given color.
# """
# @inline color(c::Colorant, s...) = Color((color=c,), s...)
# @inline color(c::AbstractString, s...) =
# 	color(parse(Colorant, c), s...)
# @inline color(c::AbstractString, a::Real, s...) =
# 	color(Colors.coloralpha(parse(Colorant, c), a), s...)
# 
# # Linear extrusion««2
# LinearExtrude = Transform{:linear_extrude,3,2}
# """
#     linear_extrude(h, s...)
#     linear_extrude(h) * s...
# 
# Linear extrusion to height `h`.
# """
# @inline linear_extrude(h, scale::AbstractVector, s...; center=false)=
# 	LinearExtrude((height=h, scale=scale, center=center,), s...)
# @inline linear_extrude(h, scale::Real, s...; kwargs...) =
# 	linear_extrude(h, SA[scale, scale], s...; kwargs...)
# @inline linear_extrude(h, s...; kwargs...) =
# 	linear_extrude(h, 1, s...; kwargs...)
# 
# # Rotational extrusion««2
# """
#     rotate_extrude([angle = 360°], solid...)
#     rotate_extrude([angle = 360°]) * solid
# 
# Similar to OpenSCAD's `rotate_extrude` primitive.
# """
# @inline rotate_extrude(s...) = rotate_extrude(360, s...)
# @inline rotate_extrude(angle::Real, s...) =
# 	RotateExtrude((angle=angle,), s...)
# RotateExtrude = Transform{:rotate_extrude,3,2}
# # Offset
# """
#     offset(r, solid...; kwargs...)
#     offset(r; kwargs...) * solid
# 
# Offsets by given radius.
# 
#     ends=:round|:square|:butt|:loop
#     join=:round|:miter|:square
# """
# @inline offset(r::Real, s...; join=:round, miter_limit=2.) =
# 	Offset((r=r, join=join, miter_limit=miter_limit), s...)
# Offset = Transform{:offset}
# Syntactic sugar: operators««2
@inline +(v::AbstractVector, x::AbstractGeometry) = translate(v,x)
@inline +(x::AbstractGeometry, v::AbstractVector) = translate(v,x)

# this purposely does not define a method for -(x::AbstractGeometry).
@inline Base.:-(x::AbstractGeometry, y::AbstractGeometry, tail...) =
	difference(x, [y, tail...])
@inline Base.:-(x::AbstractGeometry{D}) where{D} = difference(intersect(), x)
@inline Base.:-(x::AbstractVector{<:AbstractGeometry},
                y::AbstractVector{<:AbstractGeometry}) = difference(x, y)

# @inline *(f::AbstractAffineMap, x::AbstractGeometry) = mult_matrix(f, x)
@inline *(s::Union{Real,AbstractVector}, x::AbstractGeometry) = scale(s,x)

⋃ = Base.union
⋂ = Base.intersect
# Offset ««2

#————————————————————— Meshing (2d) —————————————————————————————— ««1

#»»1
# Generic code for 2d and 3d meshing««1
mesh(s::AbstractGeometry) = mesh(s, _DEFAULT_PARAMETERS)
# “thickness” of points, edges etc. for computing intersections:
# Transformations««2
# @inline mesh(s::SetParameters, parameters) =
# 	mesh(s.child, merge(parameters, s.data))
# Generic case (e.g. `color`): do nothing
# @inline mesh(s::Transform, parameters) = mesh(s.child, parameters)

# Polyhedra interface««1
@inline polyhedra_lib(T::Type{<:Real}) =
	Polyhedra.DefaultLibrary{T}(GLPK.Optimizer)

# fixing an oversight in Polyhedra.jl: it has multiplications but no
# divisions
@inline Base.:(/)(h::Polyhedra.HyperPlane, α::Real) =
	Polyhedra.HyperPlane(h.a / α, h.β / α)

# converts path to matrix with points as rows:
@inline poly_vrep(points::Path) = vcat(transpose.(Vector.(points))...)
@inline poly_vrep(points::Matrix) = points
@inline poly_eltype(points::Path) = eltype(eltype(points))
@inline poly_eltype(points::Matrix) = eltype(points)
"""
    vpoly(points...)

Returns a `Polyhedra.polyhedron` in vrep from a list of points.
"""
@inline function vpoly(points; lib=true)
	PH = Polyhedra
	if lib
		return PH.polyhedron(PH.vrep(poly_vrep(points)),
			polyhedra_lib(poly_eltype(points)))
	else
		return PH.polyhedron(PH.vrep(poly_vrep(points)))
	end
end

# HRepElement is the supertype of HalfSpace and HyperPlane
@inline direction(h::Polyhedra.HRepElement) = h.a
@inline function normalize(h::Polyhedra.HRepElement)
	n = norm(direction(h))
	(n ≠ 0) ? (h / n) : h
end
@inline (h::Polyhedra.HRepElement)(p) = h.a ⋅ coordinates(p) - h.β
@inline ∈(p, h::Polyhedra.HyperPlane) = iszero(h(p))
@inline Base.convert(T::Type{<:Polyhedra.HRepElement},
		h::Polyhedra.HRepElement) = T(h.a, h.β)

# Intersections (2d)««1

"""
    inter(path, hyperplane::Polyhedra.HyperPlane)

intersection of simplex and hyperplane
"""
function inter(simplex, hyperplane::Polyhedra.HyperPlane)
	n = length(simplex)
	s = [hyperplane(p) for p in simplex]
	newpath = similar(simplex, n); c = 0
	for i in 1:n
		if s[i] == 0
			newpath[c+= 1] = simplex[i]
		end
		for j in 1:i-1
			# separating these two cases avoids some painful `-0.` expressions:
			if s[i] < 0 < s[j]
				newpath[c+= 1] = (s[j]*simplex[i]-s[i]*simplex[j])/(s[j]-s[i])
			elseif s[j] < 0 < s[i]
				newpath[c+= 1] = (s[i]*simplex[j]-s[j]*simplex[i])/(s[i]-s[j])
			end
		end
	end
	return newpath[1:c]
end
"""
    inter(path, halfplane)

Computes intersection of a (planar) convex closed loop and the half-plane [h≥0].
The intersection is returned as a vector of points.
"""
function inter(path::Path{2}, halfplane::Polyhedra.HalfSpace)
	s = [halfplane(p) for p in path]
	boundary = convert(Polyhedra.HyperPlane, halfplane)
	n = length(path)
	# we know that we add at most 1 new point (cutting a corner).
	newpath = similar(path, 0); sizehint!(newpath, n+1)
	for i in eachindex(path)
		j = mod1(i+1, n)
		(si, sj) = (s[i], s[j])
		(si >= 0) &&  push!(newpath, path[i])
		if si*sj < 0
		# whiskers would generate two new points; we remove the second one
			newpoint = inter(path[[i,j]], boundary)[1]
			(c==0|| newpath[c] != newpoint) && push!(newpath, newpoint)
		end
	end
	return newpath
end
# @inline inter(path::AnyPath, h::Polyhedra.HRepElement,
# 		t::Polyhedra.HRepElement...) =
# 	inter(inter(path, h), t...)

"""
    line(p1=>p2)
"""
# XXX
function line(p12::Pair{<:StaticVector{2}})
	(x1, y1) = coordinates(p12[1])
	(x2, y2) = coordinates(p12[2])
	a = SA[y1-y2, x2-x1]
	b = y1*x2 - x1*y2
	return Polyhedra.HyperPlane(a, b)
end
"""
    halfplane(p1=>p2, p3)

Returns the half-plane through (p1, p2) such that h(p3) > 0.
"""
function halfplane(p12::Pair{<:StaticVector{2}}, p3::StaticVector{2})
	l = line(p12)
# 	(x1, y1) = p12[1]
# 	(x2, y2) = p12[2]
# 	a = SA[y1-y2, x2-x1]
# 	b = y1*x2 - x1*y2
	s = sign(l.a ⋅ coordinates(p3) - l.β)
	return Polyhedra.HalfSpace(s*l.a, s*l.β)
end

#Convex hull (3d)««1
"""
    convex_hull(x::Geometry{3}...)

Returns the convex hull of the union of all the given solids, as a
pair `(points, faces)`. `faces` is a list of triangles.
"""
@inline convex_hull(x::AbstractGeometry{3}) =
	convex_hull(vcat([vertices(y) for y in x]...))

"""
    convex_hull(vector of 3d points)

Returns the convex hull of these points, as a pair `(points, faces)`.
All the faces are triangles.
"""
function convex_hull(p::AbstractVector{<:StaticVector{3,T}}) where{T}
	M = hcat(Vector.(coordinates.(p))...)
	PH = Polyhedra
	poly = PH.polyhedron(PH.vrep(transpose(M)), polyhedra_lib(T))
	R = PH.removevredundancy!(poly)
	V = SVector{3,T}.(collect(PH.points(poly)))

	triangles = SVector{3,Int}[]
	for i in PH.eachindex(PH.halfspaces(poly)) # index of halfspace
		h = PH.get(poly, i)
		pts = PH.incidentpointindices(poly, i) # vector of indices of points
		for t in triangulate_face(
				[SVector{3}(PH.get(poly, j)) for j in pts];
				direction = h.a,
				map = [j.value for j in pts],
				convex = Val(true))
			(a,b,c) = (V[j] for j in t)
			k = det([b-a c-a h.a])
			push!(triangles, (k > 0) ? t : SA[t[1], t[3], t[2]])
		end
	end
	return (points=V, faces=triangles)
end

# # 2d Minkowski sum««1
# # Convolution of polygons««2
# # http://acg.cs.tau.ac.il/tau-members-area/general%20publications/m.sc.-theses/thesis-lienchapter.pdf
# """
#     circularcmp(v1, v2, v3, [Val(:offset)])
# 
# Circular comparison predicate; returns true iff directions of vectors
# `v1`, `v2`, `v3` are arranged in a trigonometric ordering along the unit
# circle.
# 
# If `Val(:offset)` is passed then `v1`, `v3` are infinitesimally rotated
# in the positive direction compared to `v2`.
# """
# function circularcmp(v1, v2, v3)
# 	d1 = v2[1]*v3[2] ≥ v2[2]*v3[1]
# 	d2 = v3[1]*v1[2] ≥ v3[2]*v1[1]
# 	d3 = v1[1]*v2[2] ≥ v1[2]*v2[1]
# 	return (d1+d2+d3) ≥ 2
# end
# function circularcmp(v1, v2, v3, ::Val{:offset})
# 	d1 = v2[1]*v3[2] > v2[2]*v3[1]
# 	d2 = v3[1]*v1[2] ≥ v3[2]*v1[1]
# 	d3 = v1[1]*v2[2] ≥ v1[2]*v2[1]
# 	return (d1+d2+d3) ≥ 2
# end
# 
# # TODO: 3d Minkowski««2
# 
# # 2d meshing««1
# # Primitive objects««2
# @inline vertices(c::Circle, parameters) = unit_n_gon(c.radius, parameters)
# 
# mesh(s::PolygonXor, parameters) = s
# mesh(s::Circle, parameters) = PolygonXor{get_parameter(parameters, :type)}(vertices(s, parameters))
# 
# # Transforms««2
# # Reduction of CSG operations««2
# function mesh(s::CSGHull{2}, parameters)
# 	l = [mesh(x, parameters) for x in children(s)]
# 	return PolygonXor(convex_hull([vertices.(l)...;]))
# end
# 
# function mesh(s::CSGMinkowski{2}, parameters)
# 	l = [mesh(x, parameters) for x in children(s)]
# 	global G = minkowski(l[1], l[2])
# # 	return PolygonXor(reduce((p,q)->minkowski(p,q), l)...)
# end
# 
# function mesh(s::AffineTransform{2}, parameters)
# 	g = mesh(s.child, parameters)
# 	f = s.data
# 	b = sign(f)
# 	@assert b ≠ 0 "Only invertible linear transforms are supported (for now)"
# 	return f(g) # reversal (if b < 0) is done polygon-by-polygon there
# end
# # Set-wise operations:
# # # Minkowski sum:
# # function (U::Type{<:PolyUnion})(s::ConstructedSolid{2,:minkowski},
# # 	parameters)::U
# # 	reduce((a,b)->U(minkowski(a.poly, b.poly)),
# # 		_convert(U, s.children, parameters))
# # end
# # function _combine2(::Val{:minkowski}, a::PolyUnion{T}, b::PolyUnion{T}) where{T}
# # 	# not implemented in Clipper.jl...
# # end
# 
# # Offset and draw««2
# function mesh(s::Offset, parameters)
# 	T = coordtype(s)
# 	m = mesh(s.child, parameters)
# 	ε = max(get_parameter(parameters,:accuracy), get_parameter(parameters,:precision) * s.data.r)
# 	return PolygonXor(offset(vertices.(paths(m)), s.data.r;
# 		join = s.data.join,
# 		ends = :fill,
# 		miter_limit = s.data.miter_limit,
# 		precision = ε)...)
# end
# 
# # function mesh(s::Draw, parameters)
# # 	r = one_half(s.width)
# # 	ε = max(get_parameter(parameters,:accuracy), get_parameter(parameters,:precision) * r)
# # 	p = offset([s.path], r; join=s.join, ends=s.ends, miter_limit = s.miter_limit)
# # 	return PolygonXor(p...)
# # end
# 
# # """
# #     draw(path, width; kwargs...)
# # 
# #     ends=:round|:square|:butt|:closed
# #     join=:round|:miter|:square
# # """
# # function draw(path::Path{2,T}, width::Real;
# # 		ends::Symbol = :round, join::Symbol = :round,
# # 		miter_limit::Float64 = 2.0, precision::Real = 0.2) where{T}
# # 	CT = clipper_type(T)
# # 	RT = clipper_rettype(T)
# # 	c = ClipperOffset(miter_limit, clipper_float(CT, precision))
# # 	println("join=$join, round=$round")
# # 	Clipper.add_path!(c, clipper_path(path),
# # 		JoinType(Val(join)), EndType(Val(ends)))
# # 	println("$(clipper_type(T)) $(CT(1.)); prec=$(Float64(CT(precision)))")
# # 	ret = clipper_unpath.(RT, Clipper.execute(c, clipper_float(CT, width)/2))
# # 	return PolyUnion(ret)
# # end
# # 
# # # Convex hull««2
# # # """
# # # 		convex_hull(x::Geometry{2}...)
# # # 
# # # Returns the convex hull of the union of all the given solids, as a
# # # `PolyUnion` structure.
# # # """
# # @inline convex_hull(x::Geometry{2}...) =
# # 	convex_hull(PolyUnion(union(x...)))
# # 
# # @inline convex_hull(u::PolyUnion) = convex_hull(Vec{2}.(vertices(u)))
# # 
# 
#————————————————————— Meshing (3d) —————————————————————————————— ««1
#=
mesh(geom, parameters = _DEFAULT_PARAMETERS) returns a Surface
mesh(x::Surface, parameters...) = x
=#
# function mesh(s::AffineTransform{3}, parameters)
# 	g = mesh(s.child, parameters)
# 	b = sign(s.data)
# 	@assert b ≠ 0 "Only invertible linear transforms are supported (for now)"
# 	if b > 0
# 		return (typeof(g))(s.data.(vertices(g)), faces(g))
# 	else
# 		return (typeof(g))(s.data.(vertices(g)), reverse.(faces(g)))
# 	end
# end
#»»1
# Converting 3d objects to Surfaces««1
# Primitive objects««2
mesh(s::Surface, parameters) = s

function vertices(s::Cube, parameters)
	z = zero(coordtype(s))
	v = s.size
	return [
		SA[z   ,z   ,z   ],
		SA[z   ,z   ,v[3]],
		SA[z   ,v[2],z   ],
		SA[z   ,v[2],v[3]],
		SA[v[1],z   ,z   ],
		SA[v[1],z   ,v[3]],
		SA[v[1],v[2],z   ],
		SA[v[1],v[2],v[3]],
	]
end
# function vertices(c::Cylinder, parameters)
# 	p1 = unit_n_gon(c.r1, parameters)
# 	p2 = unit_n_gon(c.r2, parameters)
# 	return vcat([ c.origin + [ p; 0 ] for p in p1],
# 	            [ c.origin + [p; c.height ] for p in p2 ])
# end
@inline vertices(s::Sphere, parameters) =
	[ s.center + p for p in fibonacci_sphere_points(s.radius, parameters) ]

# All of these are convex, so we use the lazy approach and just take
# convex hull of all the points.
# function mesh(s::Union{Cylinder, Sphere}, parameters)
# 	p = vertices(s, parameters)
# 	(pts, faces) = convex_hull(p)
# 	return triangulate(pts, faces)
# end
# CSG operations««2
function mesh(s::CSGHull{3}, parameters)
	l = [mesh(x, parameters) for x in children(s)]
	(pts, faces) = convex_hull(vcat(vertices.(l)...))
	return Surface(pts, faces)
end
#————————————————————— Extra tools —————————————————————————————— ««1
#»»1
# OpenSCAD output
# Return top-level objects from included file««1

# FIXME: replace Main by caller module?
# FIXME: add some blob to represent function arguments
"""
		ConstructiveGeometry.include(file::AbstractString, f::Function)

Reads given `file` and returns the union of all top-level `Geometry`
objects (except the results of assignments) found in the file.

```
#### Example: contents of file `example.jl`
C=ConstructiveGeometry.cube(1)
S=ConstructiveGeometry.square(1)
ConstructiveGeometry.circle(3)
S

julia> ConstructiveGeometry.include("example.jl")
union() {
 circle(radius=3.0);
 square(size=[1.0, 1.0], center=false);
}
```

"""
function include(file::AbstractString)
	global toplevel_objs = AbstractGeometry[]
	Base.include(x->expr_filter(obj_filter, x), Main, file)
	return union(toplevel_objs...)
end
# # TODO: somehow attach a comment indicating the origin of these objects
# # last_linenumber holds the last LineNumberNode value encountered before
# # printing this object; we use this to insert relevant comments in the
"""
    obj_filter(x)

Appends `x` to the global list of returned objects if `x` is a `AbstractGeometry`.
"""
@inline obj_filter(x) = x
@inline obj_filter(x::AbstractGeometry) =
	(global toplevel_objs; push!(toplevel_objs, x); return x)

"""
		expr_filter(E)

Read top-level expression `E` and decides if a ConstructiveGeometry object is returned.

The function `expr_filter` transforms top-level expressions by wrapping
each non-assignment expression inside a call to function `obj_filter`,
which can then dispatch on the type of the expression.
"""
# Numeric values, LineNumber expressions etc. will never be useful to us:
expr_filter(f::Function, e::Any) = e
# expr_filter(e::LineNumberNode) = (global last_linenumber = e)
# A Symbol might be an AbstractGeometry variable name, so we add it:
expr_filter(f::Function, e::Symbol) = :($f($e))

# we add any top-level expressions, except assignments:
expr_filter(f::Function, e::Expr) = expr_filter(f, Val(e.head), e)
expr_filter(f::Function, ::Val, e::Expr) = :($f($e))
expr_filter(f::Function, ::Val{:(=)}, e::Expr) = e

# if several expressions are semicolon-chained, we pass this down
expr_filter(f::Function, ::Val{:toplevel}, x::Expr) =
	Expr(:toplevel, expr_filter.(f, x.args)...)

# # Attachments««1
# # Anchor system««2
# """
# 		find_anchor(x::AbstractGeometry, name)
# 
# Returns the anchor (an affine rotation) found for the given `name` for
# the solid `x`.
# 
# Some types of `name` that can be used include:
#  - a symbol: either one of the six standard directions (`:left`, `:right`,
# 	 `:front`, `:back`, `:top`, `:bottom`, `:center`)
# 	 or a custom-defined label (TODO);
# 	 (*Note:* for 2-dimensional solids, `:bottom` is equivalent to `:front`
# 	 and `:top` is equivalent to `:back`)
#  - a list (tuple) of symbols, which is interpreted as the sum of the
# 	 corresponding directions;
#  - for standard convex bodies, a vector of the same dimension as `x` is
# 	 normalized to a point at the boundary of `x` (see below);
#  - a way to designate a point at the boundary of `x` (see below).
# """
# @inline find_anchor(x::AbstractGeometry, labels::NTuple{N,Symbol}) where{N} =
# 	find_anchor(x, sum([ labeled_anchor(x, l) for l in labels]))
# @inline function find_anchor(x::AbstractGeometry, label::Symbol)
# 	y = labeled_anchor(x, label)
# 	y isa Missing && error("No anchor named '$label' found in $(scad_name(x))")
# 	return find_anchor(x, y)
# end
# 
# default_positions = (
# 	left  = SA[-1,0,0],
# 	right = SA[+1,0,0],
# 	front = SA[0,-1,0],
# 	back  = SA[0,+1,0],
# 	bot	  = SA[0,0,-1],
# 	bottom= SA[0,0,-1],
# 	top	  = SA[0,0,+1],
# 	center= SA[0,0,0],
# )
# @inline labeled_anchor(x::Geometry{3}, label::Symbol) =
# 	get(default_positions, label, missing)
# @inline labeled_anchor(x::Geometry{2}, label::Symbol) =
# 	_labeled_anchor_3to2(get(default_positions, label, missing))
# @inline _labeled_anchor_3to2(::Missing) = missing
# @inline _labeled_anchor_3to2(v::StaticVector{3}) =
# # for 2d objects, allow (:left..:right, :bot..:top) as anchor names:
# 	(v[1] == 0 && v[2] == 0 && v[3] != 0) ? SA[v[1], v[3]] : SA[v[1], v[2]]
# # Define named anchors ««2
# # first column is translation, second (if existing) rotation,
# # spin is angle in 2d
# struct AnchorData{D,T,R}
# 	origin::Vec{D,T}
# 	direction::R
# 	spin::T
# 	@inline AnchorData{3,T}(o::AnyVec{3}, r::AnyVec{3}, s::Real) where{T} =
# 		new{3,T,Vec{3,T}}(o, r, s)
# 	@inline AnchorData{2,T}(o::AnyVec{2}, s::Real) where{T} =
# 		new{2,T,Nothing}(o, nothing, s)
# end
# 
# @inline AnchorData{3,T}(v::AnyVec{3}) where{T} =
# 	AnchorData{3,T}(v,SA[0,0,1], zero(T))
# @inline AnchorData{3,T}(data::Tuple{<:AnyVec{3},<:AnyVec{3}}) where{T} =
# 	AnchorData{3,T}(data[1], data[2], zero(T))
# @inline AnchorData{3,T}(data::Tuple{<:AnyVec{3},<:AnyVec{3},<:Real}) where{T} =
# 	AnchorData{3,T}(data[1], data[2], T(radians(data[3])))
# @inline text(x::AnchorData{3}) =
# 	"origin=$(Float64.(x.origin)) direction=$(Float64.(x.direction)) spin=$(Float64(x.spin))"
# 
# @inline AnchorData{2,T}(v::AnyVec{2}) where{T} =
# 	AnchorData{2,T}(v, zero(T))
# @inline AnchorData{2,T}(data::Tuple{<:AnyVec{2},<:Real}) where{T} =
# 	AnchorData{2,T}(data[1], T(radians(data[2])))
# @inline text(x::AnchorData{2}) =
# 	"origin=$(Float64.(x.origin)) angle=$(Float64(x.spin))"
# 
# @inline affine(x::AnchorData) = AffineMap(linear(x), x.origin)
# @inline linear(x::AnchorData{3}) =
# 	rotation_between(SA[0,0,1], x.direction) * RotZ(x.spin)
# @inline rotation(x::AnchorData{2}) = Angle2d(x.spin)
# 
# 
# """
# 		struct NamedAnchors
# 
# Wraps an object, adding symbolic anchors to it.
# """
# struct NamedAnchors{D,T} <: Geometry{D,T}
# 	child::Geometry{D,T}
# 	anchors::Dict{Symbol, AnchorData{D,T}}
# end
# """
# 		named_anchors(x, label=>anchor...)
# 
# Adds symbolic anchors to an object. `anchor` may be either
# 
#  - a vector: the associated anchor is a translation.
#  - a pair (origin, direction): the associated anchor is the affine
# 	 rotation to given direction.
#  - a triple (origin, direction, spin).
#  - (in 2d) a pair (origin, angle).
# """
# @inline named_anchors(x::Geometry{D,T},
# 	a::Pair{Symbol}...) where{D,T} =
# 	NamedAnchors{D,T}(x, Dict(k => AnchorData{D,T}(v) for (k,v) in a))
# 
# function scad(io::IO, x::NamedAnchors, spaces::AbstractString = "")
# 	println(io, spaces, "// Object with named anchors:")
# 	for (label, anchor) in x.anchors
# 		println(io, spaces, "// $label: $(text(anchor))")
# 	end
# 	scad(io, x.child, spaces)
# end
# 
# function labeled_anchor(x::NamedAnchors, label::Symbol)
# 	y = get(x.anchors, label, missing)
# 	if y isa Missing
# 		return get(default_anchors, label, missing)
# 	end
# end
# # Anchors for convex bodies ««2
# """
# 		find_anchor(x::Square, position::SVector{2})
# 		find_anchor(x::Circle, position::SVector{2})
# 		find_anchor(x::Cube, position::SVector{3})
# 		find_anchor(x::Cylinder, position::SVector{3})
# 		find_anchor(x::Sphere, position::SVector{3})
# 
# For a convex body, the anchor corresponding to `position` has its origin
# at the surface fo the body and maps the unit vector `[0,1]` (resp.
# `[0,0,1]` in dimension 3) to the normal vector at this position.
# 
# If `position` is zero, then the translation to the center of the body is
# returned.
# 
# The rotation is computed using `Rotations.rotation_between`. This returns
# a well-defined result even for the rotation mapping `[0,0,1]` to
# `[0,0,-1]`.
# """
# function find_anchor(x::Ortho{D}, pos::Vec{D,<:Real}) where{D}
# 	center = ~x.center*one_half(x.size) # the center of the square/cube
# 	if iszero(pos) return Translation(center) end
# 	maxc = findmax(abs.(pos))[2] # max abs value of a coordinate
# 	v = sum(pos[abs.(pos) .= maxc]), # sum of coordinates with max abs value
# 
# 	# we don't need to normalize v since `rotation_between` does it for us:
# 	AffineMap(rotation_between(SA[0,0,1], v), center + v .* x.size)
# end
# 
# function find_anchor(x::Circle, pos::Vec{2,<:Real})
# 	if iszero(pos) return Translation(SA[0,0]) end
# 	p1 = pos / sqrt(pos[1]^2+pos[2]^2)
# 	AffineMap(SA[p1[2] p1[1]; -p1[1] p1[2]], radius(x)*p1)
# end
# 
# function find_anchor(x::Sphere, pos::Vec{3,<:Real})
# 	if iszero(pos) return Translation(SA[0,0,0]) end
# 	return AffineMap(rotation_between(SA[0,0,1], pos), pos*x.radius / norm(pos))
# end
# 
# function find_anchor(x::Cylinder, pos::Vec{3,<:Real})
# 	center = ~x.center*one_half(x.h)
# 	if iszero(pos) return Translation(center) end
# 	r = sqrt(pos[1]*pos[1]+pos[2]*pos[2]) # convert to cylindrical coords
# 	if pos[3]*x.r2 > r # top face: normalize to pos[3]==1
# 		return Translation(SA[pos[1]/pos[3], pos[2]/pos[3], center+one_half(x.h)])
# 	elseif pos[3]*x.r1 < -r # bottom face: normalize to pos[3]==-1
# 		return AffineMap(rotation_between(SA[0,0,1], SA[0,0,-1]),
# 			SA[-pos[1]/pos[3], -pos[2]/pos[3], center-one_half(x.h)])
# 	end
# 	# the line equation is 2r = (r2-r1) z + (r1+r2)
# 	r3 = one_half(x.r1+x.r2+pos[3]*(x.r2-x.r1)) # radius at given z
# 	# in cyl coordinates, the contact point is (r=r3, z=pos[3]*h/2+center)
# 	p = SA[pos[1]*r3/r, pos[2]*r3/r, center + one_half(x.h)*pos[3]]
# 	# normal vector is 2 dr = (r2-r1) dz
# 	n = SA[2*pos[1]/r, 2*pos[2]/r, x.r2-x.r1]
# 	AffineMap(rotation_between(SA[0,0,1], n), p)
# end
# 
# # Coordinates on circle & sphere ««2
# """
# 		find_anchor(x::Circle, angle::Real)
# Returns anchor at point `(-sin(angle), cos(angle))` (trig. orientation,
# starting at top of circle; angle in **degrees**) with outwards normal vector
# (the start at top guarantees that angle 0 preserves upwards-pointing
# vector).
# """
# function find_anchor(x::Circle{T}, angle::Real) where{T}
# # 	a = T(radians(angle))
# 	(s, c) = T.(sincosd(a))
# 	AffineMap(SA[c -s;s c], SA[-radius(x)*s, radius(x)*c])
# end
# """
# 		find_anchor(x::Sphere, (latitude, longitude))
# 
# Returns normal vector to sphere at this position (angles in **degrees**).
# """
# function find_anchor(x::Sphere{T}, angle::AnyVec{2,<:Real}) where{T}
# # 	(lat, lon) = T.(radians.(angle))
# 	(s1, c1) = T.(sincosd(lat))
# 	(s2, c2) = T.(sincosd(lon))
# 	r = x.radius
# 	AffineMap(SA[s1*c2 -s2 c1*c2;s1*s2 c2 c1*s2; -c1 0 s1],
# 						SA[r*c1*c2, r*c1*s2, r*s1])
# end
# 
# 
# # attach ««2
# """
# 		attach(parent, {:label => child}...)
# 
# Moves (by rotations) all children so that their anchor matches the
# anchors of `parent` defined by the given labels.
# """
# function attach(parent::Geometry, list::Pair{<:Any,<:Geometry}...)
# 	union(parent, [ attach_at(parent, pos, child) for (pos,child) in list]...)
# end
# function attach_at(parent::Geometry{D}, label, child::Geometry) where{D}
# 	m = find_anchor(parent, label)
# 	mult_matrix(m, child)
# end
# 
# # anchor ««2
# """
# 		half_turn(q::UnitQuaternion)
# 
# Returns the half-turn rotation with same axis as `q`.
# """
# @inline half_turn(q::Rotations.UnitQuaternion) =
# 	Rotations.UnitQuaternion(0, q.x, q.y, q.z) # the constructor normalizes this
# """
# 		half_turn_complement(q::Rotation{3})
# 
# Returns the unique rotation `r` such that `qr=rq` is a half-turn.
# """
# @inline half_turn_complement(q::Rotations.Rotation{3}) =
# 	inv(q)*half_turn(q)
# 
# """
# 		anchor(solid, label)
# 
# Translates the solid so that the anchor with name `label` is at origin
# and the corresponding anchor vector points inward.
# """
# function anchor(x::Geometry, label)
# # Ax + B maps [0,0,0] to anchor point p, and e_z=[0,0,1] to anchor vec v
# # namely: A e_z = v, B = p
# # we want to map p to [0,0,0] and -v to e_z
# # i.e. A' v = - e_z and A' p + B' = 0, or B' = -A' p = -A' B
# #
# 
# 	m = find_anchor(x, label)
# 	a1= half_turn_complement(linear(m))
# # ⚠ -a1*m is interpreted as (-a1)*m, and a1 is a quaternion ⇒ -a1≡a1 (as
# # a rotation)
# 	mult_matrix(AffineMap(a1, a1*-translation(m)), x)
# end
# 
# # # Annotations ««1
# # abstract type AbstractAnnotation{D} end
# # 
# # struct DimensionArrow{D,T} <: AbstractAnnotation{D}
# # 	center::Vec{D,T}
# # 	vec::Vec{D,T}
# # 	label::AbstractString
# # 	fontsize::Float64
# # 	offset::Float64
# # end
# # """
# #     DiameterArrow{X}
# # 
# # Indicates that a diameter arrow should be drawn for the given object. The
# # parameter `X` is a type indicating which type of arrow should be drawn.
# # 
# #  - `Circle`: parametrized by center (`Vec{2}`) and radius (scalar),
# #  and possibly preferred orientation (vector if non-zero).
# # 
# #  - `Sphere`: parametrized by center (`Vec{3}`) and radius (scalar),
# #   and possibly preferred orientation.
# # 
# #  - `Cylinder`: shows a circle in 3d space, parametrized by center (`Vec{3}`), normal vector (non-zero), radius (scalar), and preferred orientation (vector; should be in the circle plane).
# # """
# # struct DiameterArrow{X<:Geometry,T,D,N} <: AbstractAnnotation{D}
# # 	center::Vec{D,T}
# # 	radius::T
# # 	orientation::Vec{D,T}
# # 	normal::N
# # 	# inner constructors enforce the correct type for N
# # 	DiameterArrow{Circle,T}(center, radius, orientation) where{T} =
# # 		new{Circle, T, 2, Nothing}(center, radius, orientation, nothing)
# # 	DiameterArrow{Sphere,T}(center, radius, orientation) where{T} =
# # 		new{Sphere, T, 3, Nothing}(center, radius, orientation, nothing)
# # 	DiameterArrow{Cylinder,T}(center, radius, orientation, normal) where{T} =
# # 		new{Cylinder, T, 3, Vec{3,T}}(center, radius, orientation, normal)
# # end
# # 
# # struct Annotate{D,T} <: Geometry{D,T}
# # 	annotations::Vector{<:AbstractAnnotation{D}}
# # 	child::Geometry{D,T}
# # end
# # # the offset is just a hint; we let the visualizer take care of using
# # # this
# # 
# #
# Exports ««1
export square, circle, cube, sphere, cylinder, polygon, surface
export offset, draw
export mult_matrix, translate, scale, rotate, mirror
export linear_extrude, rotate_extrude, path_extrude
export color, set_parameters
export mesh
export difference, ⋃, ⋂, offset, hull, minkowski
# don't export include, of course
# »»1
# function explain(s::AbstractSurface, io::IO = stdout; scale=1,
# 		offset=[0.,0.,0.], name=:m )
# 	println(io, """
# module $name(pos=$offset, c="gray", s=$scale) {
# translate($scale*pos) {
# """)
# 	for (i, p) in pairs(ConstructiveGeometry.vertices(s))
# 		println(io, """
# translate(s*$(Vector{Float64}(coordinates(p)))) {
# 	color("red") sphere(1);
# 	color("black", .8) linear_extrude(1) text("$i", size=5);
# }
# """)
# 	end
# 	println(io, "color(c, .7) polyhedron([")
# 	b = false
# 	for p in ConstructiveGeometry.vertices(s)
# 		print(io, b ? "," : ""); b = true
# 		print(io, " s*",Vector{Float64}(coordinates(p)))
# 	end
# 	println(io, "],[")
# 	b = false
# 	for f in ConstructiveGeometry.faces(s)
# 		print(io, b ? "," : ""); b = true
# 		println(io, " ", Vector{Int}(f) .- 1, " //", Vector{Int}(f))
# 	end
# 	println(io, "]); } }\n$name();")
# end
# @inline explain(s::AbstractSurface, f::AbstractString; kwargs...) = begin
# 	println("writing a surface with $(nvertices(s)) points to $f")
# 	open(f, "w") do io explain(s, io; kwargs...) end
# end
end #««1 module
# »»1
# vim: fdm=marker fmr=««,»» noet:
