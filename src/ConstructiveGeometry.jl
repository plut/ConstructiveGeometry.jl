module ConstructiveGeometry
#————————————————————— Tools —————————————————————————————— ««1

#»»1
# Module imports««1
# using Printf
using LinearAlgebra
using StaticArrays
using FixedPointNumbers
using Logging
using FastClosures
# using DataStructures

import Rotations
import Colors: Colors, Colorant

import Base: show, print
import Base: union, intersect, setdiff, copy, isempty, merge
import Base: *, +, -, ∈, inv, sign, iszero

# sub-packages:
include("ConvexHull.jl")
using .ConvexHull
include("Shapes.jl")
using .Shapes

# include("SpatialSorting.jl")
# include("TriangleIntersections.jl")
include("CornerTables.jl")
using .CornerTables
@inline corner_table(points, faces, c::Colorant) =
	CornerTable(points, faces, Iterators.repeated(c))
include("scad.jl")

# General tools««1
struct Consecutives{T,V} <: AbstractVector{T}
	parent::V
end
consecutives(v::AbstractVector{T}) where{T} =
	Consecutives{T,typeof(v)}(v)
Base.getindex(c::Consecutives, i::Integer) =
	(c.parent[i], c.parent[mod1(i+1, length(c.parent))])
Base.size(c::Consecutives) = size(c.parent)

@inline one_half(x::Real) = x/2
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

# Parameters««1
# Definitions««2
# Accuracy is the absolute deviation allowed.
# Default value is 2.0 (from OpenSCAD `$fs`), interpreted as 2mm.
#
# Precision is the relative deviation allowed.
# Default value is 0.02 (1-cos(180°/`$fa`)).
# FIXME: explain why .005 works better
#
# ε = “thickness” of points, edges etc. for computing intersections:

const _DEFAULT_PARAMETERS = (
	accuracy = 0.1, precision = .005, symmetry = 1,
	type = Float64, ε = 1/65536,
	color = Colors.RGBA{N0f8}(.5,.5,.5,1.), # gray
)

@inline get_parameter(parameters, name) =
	get(parameters, name, _DEFAULT_PARAMETERS[name])

#————————————————————— Geometry —————————————————————————————— ««1

#»»1
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

# Affine transformations««1
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
@inline AffineMap(a; center=nothing) = _affine_map_center(a, center)
@inline _affine_map_center(a, ::Nothing) = LinearMap(a)
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
@inline det(::Rotations.Rotation{D,T}) where{D,T} = one(T)

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


# Reflection matrices««1
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

@inline Base.getindex(r::Reflection, i::Int, j::Int) =
	(i == j) - 2*r.axis[i]*r.proj[j]
# we need two methods here for disambiguation:
@inline Base.:*(r::Reflection, v::AbstractVector) = v - 2*r.axis*(r.proj*v)
@inline Base.:*(r::Reflection, v::AbstractMatrix) = v - 2*r.axis*(r.proj*v)

# Circles and spheres««1
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
function unit_n_gon(r::T, n::Int) where{T<:Real}
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
	reinterpret(SVector{2,T}, r*z)
end
@inline unit_n_gon(r, parameters::NamedTuple)= unit_n_gon(r,sides(r,parameters))

# Spheres««2
"""
    sphere_nvertices(r::Real, parameters::NamedTuple)

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
function sphere_nvertices(r::Real, parameters)
	ε = max(get_parameter(parameters,:precision),
		get_parameter(parameters,:accuracy)/r)
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
function fibonacci_sphere_points(r::T, n::Int) where{T<:Real}
	v = Vector{SVector{3,T}}(undef, n)
	for i in eachindex(v)
		θ = i*T(golden_angle)
		z = (n+1-2i)/T(n)
		ρ = √(1-z^2)
		(s,c) = sincos(θ)
		@inbounds v[i] = SA[r*c*ρ, r*s*ρ, r*z]
	end
	return v
end
@inline fibonacci_sphere_points(r::Real, parameters::NamedTuple) =
	fibonacci_sphere_points(r, sphere_nvertices(r, parameters))
# tools for rotate extrusion««2
"""
    _rotate_extrude(point, data, parameters)

Extrudes a single point, returning a vector of 3d points
(x,y) ↦ (x cosθ, x sinθ, y).
"""
function _rotate_extrude(p::StaticVector{2,T}, angle, parameters) where{T}
	@assert p[1] ≥ 0
	# special case: point is on the y-axis; returns a single point:
	iszero(p[1]) && return [SA[p[1], p[1], p[2]]]
	n = cld(sides(p[1], parameters) * angle, 360)

	ω = Complex{T}(cosd(angle/n), sind(angle/n))
	z = Vector{Complex{T}}(undef, n+1)
	z[1] = one(T)
	for i in 2:n
		@inbounds z[i] = z[i-1]*ω; z[i]/= abs(z[i])
	end
	# close the loop:
	z[n+1] = Complex{T}(cosd(angle), sind(angle))
	return [SA[p[1]*real(u), p[1]*imag(u), p[2]] for u in z]
end
"""
    ladder_triangles(n1, n2, start1, start2)

Given two integers m, n, triangulate as a ladder between these integers;
example: ladder(5, 4, a, b)

    a──a+1──a+2──a+3──a+4
    │ ╱ ╲  ╱   ╲ ╱  ╲  │
    b────b+1────b+2───b+3

Returns the m+n-2 triangles (a,b,a+1), (b,b+1,a+1), (a+1,b+1,a+2)…
"""
function ladder_triangles(n1, n2, start1, start2)
	p1 = p2 = 1
	triangles = NTuple{3,Int}[]
	while true
		(p1 == n1) && (p2 == n2) && break
		# put up to scale (n1-1)*(n2-1):
		# front = (i,j) corresponds to ((n2-1)*i, (n1-1)*j)
		c1 = (p1)*(n2-1)
		c2 = (p2)*(n1-1)
		if c1 < c2 || ((c1 == c2) && (n1 <= n2))
			push!(triangles, (p1+start1-1, p2+start2-1, p1+start1))
			p1+= 1
		else
			push!(triangles, (p1+start1-1, p2+start2-1, p2+start2))
			p2+= 1
		end
	end
	return triangles
end

#————————————————————— Objects and meshing —————————————————————————————— ««1

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

"""
    Mesh{T}

A structure storing the parameters for the (abstract object)->(mesh) conversion.
The coordinate type determines the Julia type of the returned object;
therefore, it is stored as a type parameter of `Mesh`.
"""
struct Mesh{T<:Real,C}
	parameters::NamedTuple
	color::C
	@inline Mesh{T,C}(parameters) where{T,C} =
		new{T,C}(parameters, parameters.color)
	@inline Mesh{T}(parameters) where{T} =
		Mesh{T,typeof(parameters.color)}(parameters)
	@inline Mesh(parameters::NamedTuple) = Mesh{parameters.type}(parameters)
end
@inline coordtype(::Mesh{T}) where{T} = T
@inline mesh(s::AbstractGeometry; kwargs...) =
	Mesh(merge(_DEFAULT_PARAMETERS, kwargs.data))(s)
@inline get_parameter(g::Mesh, name) = get_parameter(g.parameters, name)

@inline corner_table(g::Mesh{T}, points, faces) where{T} =
	CornerTable{Int,SVector{3,T},typeof(g.color)}(
		points, faces, Iterators.repeated(g.color))
@inline polygon_xor(::Mesh{T}, args...) where{T} = PolygonXor{T}(args...)
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

@inline square_vertices(u, v) = [ SA[0,0], SA[u,0], SA[u,v], SA[0,v]]

@inline (::Mesh{T})(s::Square) where{T} =
	PolygonXor{T}(square_vertices(T(s.size[1]), T(s.size[2])))

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
@inline (g::Mesh{T})(s::Circle) where{T} =
	PolygonXor{T}(unit_n_gon(T(s.radius), g.parameters))

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
@inline (::Mesh{T})(s::Polygon) where{T} = PolygonXor{T}(s.points)

# Draw ««2
struct Draw{T} <: AbstractGeometryCoord{2,T}
	path::Vector{SVector{2,T}}
	width::Float64
	ends::Symbol
	join::Symbol
	miter_limit::Float64
end
Draw(path, width; ends=:round, join=:round, miter_limit=2.) =
	Draw{coordtype(path)}(path, width, ends, join, miter_limit)

function (g::Mesh{T})(s::Draw) where{T}
	r = one_half(T(s.width))
	ε = max(get_parameter(g,:accuracy), get_parameter(g,:precision) * r)
	return PolygonXor{T}(offset([s.path], r;
		join=s.join, ends=s.ends, miter_limit = s.miter_limit)...)
end

# 3d primitives««1
Cube = Ortho{3}
@inline scad_info(s::Cube) = (:cube, (size=s.size,))

@inline cube_vertices(u, v, w) = [
		SA[0,0,0], SA[0,0,w], SA[0,v,0], SA[0,v,w],
		SA[u,0,0], SA[u,0,w], SA[u,v,0], SA[u,v,w]]

(g::Mesh{T})(s::Cube) where{T} =
	corner_table(g, cube_vertices(T(s.size[1]), T(s.size[2]), T(s.size[3])),
	[ # 12 triangular faces:
	 (6, 5, 7), (7, 8, 6), (7, 3, 4), (4, 8, 7),
	 (4, 2, 6), (6, 8, 4), (5, 1, 3), (3, 7, 5),
	 (2, 1, 5), (5, 6, 2), (3, 1, 2), (2, 4, 3),
	])

Sphere = Ball{3}
@inline scad_info(s::Sphere) = (:sphere, (r=s.radius,))

function (g::Mesh{T})(s::Sphere) where{T}
	plist = fibonacci_sphere_points(T(s.radius), g.parameters)
	(pts, faces) = convex_hull(plist)
	return corner_table(g, pts, faces)
end


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
	Surface(points, [(f...,) for f in faces])

@inline scad_info(s::Surface) =
	(:surface, (points=s.points, faces = [ f .- 1 for f in s.faces ]))
@inline (g::Mesh)(s::Surface) = corner_table(g, s.points, s.faces)

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

# Union««2
struct CSGUnion{D} <: AbstractGeometry{D}
	children::Vector{<:AbstractGeometry{D}}
end
@inline (g::Mesh)(s::CSGUnion{2}) = Shapes.clip(:union, g.(s.children)...)
@inline (g::Mesh)(s::CSGUnion{3}) =
	CornerTables.combine(g.(s.children), 1, g.parameters.ε)

# Intersection««2
struct CSGInter{D} <: AbstractGeometry{D}
	children::Vector{<:AbstractGeometry{D}}
end
@inline (g::Mesh)(s::CSGInter{2}) =
	Shapes.clip(:intersection, g.(s.children)...)
@inline (g::Mesh)(s::CSGInter{3}) =
	CornerTables.combine(g.(s.children), length(s.children), g.parameters.ε)

# Difference««2
struct CSGDiff{D} <: AbstractGeometry{D} # this is a binary operator:
	children::Tuple{<:AbstractGeometry{D},<:AbstractGeometry{D}}
end
@inline (g::Mesh)(s::CSGDiff{2}) =
	Shapes.clip(:difference, g(s.children[1]), g(s.children[2]))
@inline (g::Mesh)(s::CSGDiff{3}) =
	CornerTables.combine([g(s.children[1]), reverse(g(s.children[2]))],
		2, g.parameters.ε)

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

# Convex hull (TODO)««2
struct CSGHull2
	children::Vector{AbstractGeometry{2}}
end

struct CSGHull3
	children::Vector{AbstractGeometry}
end

CSGHull = constructed_solid_type(:hull)

function (g::Mesh{T})(s::CSGHull2) where{T}
	v = reduce(vcat, Shapes.vertices(g(x)) for x in children(s))
	return PolygonXor{T}(convex_hull(v))
end

@inline vertices3(m::CornerTable) = CornerTables.points(m)
@inline vertices3(m::PolygonXor) =
	(SA[p[1], p[2], 0] for p in Shapes.vertices(m))
function (g::Mesh{T})(s::CSGHull{3}) where{T}
	v = SVector{3,T}[]
	for x in children(s); push!(v, vertices3(g(x))...); end
	(p, f) = convex_hull(v)
	return corner_table(g, p, f)
end
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

# Minkowski sum and difference (TODO)««2
CSGMinkowski = constructed_solid_type(:minkowski)
# Transformations««1
abstract type AbstractTransform{D} <: AbstractGeometry{D} end
# AffineTransform««2
struct AffineTransform{D,A<:AbstractAffineMap} <: AbstractTransform{D}
	f::A
	child::AbstractGeometry{D}
end
# AffineTransform(f, child) constructor is defined

function (g::Mesh)(s::AffineTransform{3})
	# FIXME what to do if signdet(s.f) == 0 ?
	m = g(s.child)
	m.points .= s.f.(m.points)
	return signdet(s.f) > 0 ? m : reverse(m)
end

function (g::Mesh)(s::AffineTransform{2})
	m = g(s.child)
	d = signdet(s.f)
	d ≠ 0 || error("Only invertible linear transforms are supported (for now)")
	for p in m.paths
		p .= s.f.(p)
		d < 0 && reverse!(p)
	end
	return m
end

# Project and cut (TODO)««2
struct Project <: AbstractGeometry{2}
	child::AbstractGeometry{3}
end

function (g::Mesh)(s::Project)
	m = g(s.child)
	p = CornerTables.points(m)
	f = CornerTables.faces(m)
	return PolygonXor([SA[x[1], x[2]] for x in p], f)
end

struct Cut <: AbstractGeometry{2}
	child::AbstractGeometry{3}
end

function (g::Mesh)(s::Cut)
	(pts, seg) = CornerTables.plane_cut(g(s.child), get_parameter(g, :ε))
	return Shapes.glue_segments(pts, seg)
end
# SetParameters««2
# (including colors)
struct SetParameters{D} <: AbstractTransform{D}
	parameters
	child::AbstractGeometry{D}
	@inline SetParameters(parameters, child::AbstractGeometry{D}) where{D} =
		new{D}(parameters, child)
end

@inline (g::Mesh)(s::SetParameters) =
	(typeof(g))(merge(g.parameters, s.parameters))(s.child)
# Linear extrusion««2
struct LinearExtrude{T} <: AbstractTransform{3}
	height::T
	child::AbstractGeometry{2}
end

function (g::Mesh)(s::LinearExtrude)
	m = g(s.child)
	@assert m isa PolygonXor
	pts2 = Shapes.vertices(m)
	tri = Shapes.triangulate(m)
	peri = Shapes.perimeters(m)
	# perimeters are oriented ↺, holes ↻

	n = length(pts2)
	pts3 = vcat([[SA[p..., z] for p in pts2] for z in [0, s.height]]...)
	# for a perimeter p=p1, p2, p3... outward: ↺
	# with top vertices q1, q2, q3...
	# faces = [p1,p2,q1], [p2,q2,q1], [p2,p3,q2],  [p2,q3,q2]...
	#  - bottom: identical to tri
	#  - top: reverse of tri + n
	#  - sides:
	faces = [ tri;
	  [ reverse(f) .+ n for f in tri ];
	  vcat([[(i,j,i+n) for (i,j) in consecutives(p) ] for p in peri]...);
	  vcat([[(j,j+n,i+n) for (i,j) in consecutives(p) ] for p in peri]...);
	]
	return corner_table(g, pts3, faces)
end

# Cone««2
# Note: this is *not* a convex hull - a generic cone is not convex
# (the base shape may be nonconvex, or even have holes).
struct Cone{T} <: AbstractTransform{3}
	apex::SVector{3,T}
	child::AbstractGeometry{2}
end

function (g::Mesh{T})(s::Cone) where{T}
	m = g(s.child)
	@assert m isa PolygonXor
	pts2 = Shapes.vertices(m)
	tri = Shapes.triangulate(m)
	peri = Shapes.perimeters(m)
	n = length(pts2)
	pts3 = [SVector{3,T}[p[1],p[2], 0] for p in pts2]
	push!(pts3, SVector{3,T}(s.apex))
	faces = [ tri;
		vcat([[(i,j,n+1) for (i,j) in consecutives(p)] for p in peri]...);]
	return corner_table(g, pts3, faces)
end
# Rotate extrusion««2
struct RotateExtrude{T} <: AbstractTransform{3}
	angle::T
	child::AbstractGeometry{2}
end

function (g::Mesh{T})(s::RotateExtrude) where{T}
	# right half of child:
	m0 = g(s.child)::PolygonXor{T}
	m = intersect(Shapes.HalfPlane(SA[1,0],0), m0)
	pts2 = Shapes.vertices(m)
	tri = Shapes.triangulate(m)
	peri = Shapes.perimeters(m) # oriented ↺
	n = length(pts2)
	
	pts3 = _rotate_extrude(pts2[1], s.angle, g.parameters)
	firstindex = [1]
	arclength = Int[length(pts3)]
	@debug "newpoints[$(pts2[1])] = $(length(pts3))"
	for p in pts2[2:end]
		push!(firstindex, length(pts3)+1)
		newpoints = _rotate_extrude(p, s.angle, g.parameters)
		@debug "newpoints[$p] = $(length(newpoints))"
		push!(arclength, length(newpoints))
		pts3 = vcat(pts3, newpoints)
	end
	@debug "pts3: $pts3"
	@debug "firstindex: $firstindex"
	@debug "arclength: $arclength"
	# point i ∈ polygonxor: firstindex[i]:(firstindex[i]-1+arclength[i])
	triangles = vcat(
		[ firstindex[t] for t in tri ],
		[ firstindex[t] .+ arclength[t] .- 1 for t in reverse.(tri) ]
	)
	@debug "triangles: $(Vector.(triangles))"
	for l in peri
		@debug "triangulating perimeter l=$l««"
		for (i1, p1) in pairs(l)
			i2 = mod1(i1+1, length(l)); p2 = l[i2]
			# point p1: firstindex[p1] .. firstindex[p1]-1+arclength[p1]
			# point p2: etc.
			@debug "triangulating between points $p1 and $p2:"
			@debug "   $(firstindex[p1])..$(firstindex[p1]-1+arclength[p1])"
			@debug "   $(firstindex[p2])..$(firstindex[p2]-1+arclength[p2])"
			nt = SVector.(ladder_triangles(
				arclength[p2], arclength[p1],
				firstindex[p2], firstindex[p1],
				))
			push!(triangles, nt...)
			@debug "new triangles: $(Vector.(nt))"
		end
		@debug "(end perimeter)»»"
	end
	@debug "triangles = $triangles"
	return corner_table(g, pts3, triangles)
end


# Offset««2

struct Offset <: AbstractTransform{2}
	radius::Float64
	ends::Symbol
	join::Symbol
	miter_limit::Float64
	child::AbstractGeometry{2}
end

function (g::Mesh)(s::Offset)
	m = g(s.child)
	ε = max(get_parameter(g,:accuracy), get_parameter(g,:precision) * s.radius)
	return PolygonXor(offset(vertices.(paths(m)), s.radius;
	join = s.join, ends = s.ends, miter_limit = s.miter_limit, precision = ε)...)
end

#————————————————————— Front-end —————————————————————————————— ««1

# Constructors««1
# Squares and cubes««2
# TODO: add syntactic sugar (origin, center, etc.) here:
"""
    square(size; origin, center=false)

An axis-parallel square or rectangle  with given `size`
(scalar or vector of length 2).
"""
@inline square(a::Real) = Square(a,a)
@inline square(a::AbstractVector) = Square(a...)

@inline cube(a::Real) = Cube(a,a,a)
@inline cube(a::AbstractVector) = Cube(a...)

# Circles and spheres««2
@inline circle(a::Real) = Circle(a)
@inline sphere(a::Real) = Sphere(a)

# TODO: cylinders and cones««2
"""
    cylinder(h, r1, r2 [, center=false])
    cylinder(h, (r1, r2) [, center=false])
    cylinder(h, r [, center=false])

**Warning:** `cylinder(h,r)` is interpreted as `cylinder(h,r,r)`,
not `(h,r,0)` as in OpenSCAD.
"""
@inline cylinder(h::Real, r::Real) = linear_extrude(h)*circle(r)

"""
    cone(h, r)
    cone(apex, r)

Circular cone.
"""
@inline cone(v::AbstractVector, r::Real) = Cone(v, circle(r))
@inline cone(h::Real, r::Real) = cone(SA[0,0,v], r)

# Polygon««2
"""
    polygon(path)
"""
@inline polygon(path) = Polygon(path)
# Draw««2
"""
    draw(path, width; kwargs)
    ends = :loop|:butt|:square|:round
    join = :square|:round|:miter
    miter_limit = 2.0

Draws a path of given width.
"""
draw(path, width; kwargs...) = Draw(path, width; kwargs...)


# Surface««2
function surface(points, faces)
	# FIXME: triangulate!
	return Surface(points, faces)
end

# CSG operations««1
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

# Unrolling: prevent nesting of similar constructions
"""
		unroll(x::AbstractGeometry, Val(sym1), Val(sym2)...)

Returns either `[x]` or, if `x` is a `ConstructedSolid` matching one of the
symbols `sym1`, `sym2`..., `children(x)`.
"""
@inline unroll(s::AbstractGeometry, ::Val, tail...) = unroll(s, tail...)
@inline unroll(s::AbstractGeometry) = s
@inline unroll(s::ConstructedSolid{D, S}, ::Val{S}, tail...) where{D, S} =
	children(s)
@inline unroll2(s::AbstractGeometry, t::AbstractGeometry, tail...) =
	[unroll(s, tail...); unroll(t, tail...)]

# Booleans««2
@inline union(a::AbstractGeometry{D}, b::AbstractGeometry{D}) where{D} =
	CSGUnion{D}(unroll2(a, b, Val(:union)))

@inline intersect(a::AbstractGeometry{D}, b::AbstractGeometry{D}) where{D} =
	CSGInter{D}(unroll2(a, b, Val(:intersection)))
@inline Base.intersect(a::AbstractGeometry{3}, b::AbstractGeometry{2}) =
	intersect(cut(a), b)
@inline Base.intersect(a::AbstractGeometry{2}, b::AbstractGeometry{3}) =
	intersect(a, cut(b))

@inline setdiff(x::AbstractGeometry{D}, y::AbstractGeometry{D}) where{D} =
	CSGDiff{D}((x,y))
@inline setdiff(a::AbstractGeometry{2}, a::AbstractGeometry{3}) =
	setdiff(a, cut(b))

# added interface: setdiff([x...], [y...])
@inline setdiff(x::AbstractVector{<:AbstractGeometry},
                y::AbstractVector{<:AbstractGeometry}) =
	setdiff(union(x...), union(y...))
# Convex hull««2
"""
    hull(s::AbstractGeometry...)
		hull(s::AbstractGeometry{2} | StaticVector{2}...)

Represents the convex hull of given solids (and, possibly, points).
"""
@inline hull(s::AbstractGeometry{2}...) =
	CSGHull{2}([unroll(t, Val.((:hull, :union))...) for t in s])

function hull(s::Union{AbstractGeometry{2},StaticVector{2,<:Real}}...)
	l = AbstractGeometry{2}[unroll(t, Val.((:hull, :union))...)
		for t in s if t isa AbstractGeometry]
	v = filter(x->x isa AbstractVector, s)
	push!(l, polygon([v...]))
	return CSGHull{2}(l)
end

function hull(s::Union{AbstractGeometry,AbstractVector}...)
	l = AbstractGeometry[]
	T = Bool; v = SVector{3,T}[]
	for x in s
		if x isa AbstractVector
			y = (length(x) == 3) ? SVector{3}(x) :
			    (length(x) == 2) ? SVector{3}(x[1], x[2], 0) :
					error("bad point dimension: $(length(x))")
			T = promote_type(T, eltype(y))
			v = push!(SVector{3,T}.(v), y)
		else
			push!(l, x)
		end
	end
	!isempty(v) && push!(l, surface(v, [(i,i,i) for i in 1:length(v)]))
	return CSGHull{3}(l)
end
# Minkowski««2
"""
    minkowski(s::AbstractGeometry...)

Represents the Minkowski sum of given solids.
"""
@inline minkowski(a1::AbstractGeometry, a2::AbstractGeometry) =
	CSGMinkowski{maximum(embeddim.((a1,a2)))}(unroll2(a1, a2, Val(:minkowski)))


# Transformations««1
# Transform type««2
# An operator F may be called in either full or currified form:
# F(parameters, object) - applied
# F(parameters) - abstract transform (allows chaining)
"""
    Transform{S,F}

This type defines functions allowing to chain transforms; these are used by
`multmatrix`, `color` etc. operations (see below).

The minimal job left to concrete types (see e.g. `AffineTransform` as an
example) is to define a type and a constructor:
    Frobnicate = Transform{:frobnicate}
		frobnicate(x::real, s...) = Frobnicate((x=x,), s...)
"""
struct Transform{S,F}
	f::F
	@inline Transform{S}(f) where{S} = new{S,typeof(f)}(f)
end
@inline operator(f, x) = Transform{f}(@closure y -> f(x..., y))
@inline operator(f, x, s) = operator(f, x)*s
@inline operator(f, x, s::AbstractVector) = operator(f, x, s...)
@inline operator(f, x, s...) = operator(f, x, union(s...))
@inline Base.:*(u::Transform, s) = u.f(s)
@inline (u::Transform)(s) = u*s
@inline Base.:*(u::Transform, v::Transform, s...)= _comp(Val(assoc(u,v)),u,v,s...)
@inline _comp(::Val{:left}, u, v, s...) = *(compose(u,v), s...)
@inline _comp(::Val{:right}, u, v, s...) = *(u, *(v, s...))
@inline assoc(::Transform, ::Transform) = :right
@inline Base.:*(u::Transform, v::Transform) = compose(u, v)
@inline compose(u::Transform, v::Transform) = Transform{typeof(∘)}(u.f∘v.f)
# extrusions««2
"""
    linear_extrude(h, s...)
    linear_extrude(h) * s...

Linear extrusion to height `h`.
"""
@inline linear_extrude(height, s...; center=false) =
	operator(_linear_extrude, (height, center), s...)
@inline function _linear_extrude(height, center, s)
	m = LinearExtrude(height, s)
	center ? translate([0,0,-one_half(height)], m) : m
end

"""
    cone(h, shape)
    cone(h)*shape
    cone(apex, shape)
    cone(apex)*shape

Cone with arbitrary base.
"""
@inline cone(v::AbstractVector, s...) = operator(Cone, (v,), s...)
@inline cone(h::Real, s...) = cone(SA[0,0,h], s...)

"""
    rotate_extrude([angle = 360°], solid...)
    rotate_extrude([angle = 360°]) * solid

Similar to OpenSCAD's `rotate_extrude` primitive.
"""
@inline rotate_extrude(angle::Real, s...) =
	operator(RotateExtrude, (angle,), s...)
@inline rotate_extrude(s...) = rotate_extrude(360, s...)

# offset««2
"""
    offset(r, solid...; kwargs...)
    offset(r; kwargs...) * solid

Offsets by given radius.

    ends=:round|:square|:butt|:loop
    join=:round|:miter|:square
"""
@inline offset(r::Real, s...; join=:round, miter_limit=2.) =
	Transform(Offset,((r=r, join=join, miter_limit=miter_limit),), s...)

# set_parameters««2
"""
    set_parameters(;accuracy, precision, symmetry, ε, type) * solid...

A transformation which passes down the specified parameter values to its
child. Roughly similar to setting `\$fs` and `\$fa` in OpenSCAD.
"""
@inline set_parameters(s...; parameters...) =
	operator(SetParameters, (parameters.data,), s...)

"""
    color(c::Colorant, s...)
    color(c::AbstractString, s...)
    color(c::AbstractString, α::Real, s...)
    color(c) * s...

Colors objects `s...` in the given color.
"""
@inline color(c::Colorant, s...) = operator(SetParameters, ((color=c,),), s...)
@inline color(c::AbstractString, s...) = color(parse(Colorant, c), s...)
@inline color(c::AbstractString, a::Real, s...) =
	color(Colors.coloralpha(parse(Colorant, c), a), s...)

# Affine transformations««1
# mult_matrix««2
"""
    mult_matrix(a, [center=c], solid...)
    mult_matrix(a, b, solid...)
    mult_matrix(a, b) * solid

Represents the affine operation `x -> a*x + b`.

# Extended help
!!! note "Types of `mult_matrix` parameters"

    The precise type of parameters `a` and `b` is not specified.
    Usually, `a` will be a matrix and `b` a vector, but this is left open
    on purpose; for instance, `a` can be a scalar (for a scaling).
    Any types so that `a * Vector + b` is defined will be accepted.

    Conversion to a matrix will be done when converting to OpenSCAD
    format.

!!! note "Matrix multiplication"

    Chained `mult_matrix` operations will be combined into a single
    operation when possible. This saves time: multiple
    (3 × n) matrix multiplications are replaced by
    (3 × 3) multiplications, followed by a single (3 × n).
"""
@inline mult_matrix(a, s...; center=nothing) =
	operator(AffineTransform, (_affine_map_center(a, center),), s...)
# @inline affine(a, b::AbstractVector, s...) =
# 	operator(AffineTransform, (AffineMap(a, b),), s...)
# translate ««2
"""
    translate(v, s...)
    translate(v) * s

Translates solids `s...` by vector `v`.
"""
@inline translate(a,s...)= operator(AffineTransform,(TranslationMap(a),), s...)

# scale««2
# `scale` is really the same as one of the forms of `mult_matrix`:
"""
    scale(a, s...; center=0)
    scale(a; center=0) * s
Scales solids `s` by factor `a`. If `center` is given then this will be
the invariant point.

`a` may also be a vector, in which case coordinates will be multiplied by
the associated diagonal matrix.
"""
@inline scale(a,s...; kwargs...) = mult_matrix(scaling(a), s...; kwargs...)
@inline scaling(a) = a
@inline scaling(a::AbstractVector) = Diagonal(a)

# mirror««2
"""
    mirror(v, s...; center=0)
    mirror(v; center=0) * s

Reflection with axis given by the hyperplane normal to `v`.
If `center` is given, then the affine hyperplane through this point will
be used.
"""
@inline mirror(v::AbstractVector, s...; kwargs...) =
	mult_matrix(Reflection(v), s...; kwargs...)
# rotate««2
# FIXME: add a method Angle2d*StaticVector{3}
@inline rotation(θ::Real, ::Nothing) = Rotations.Angle2d(θ)
@inline rotation(θ::Real, axis) = Rotations.AngleAxis(θ, axis)
"""
    rotate(θ, {center=center}, {solid...})
    rotate(θ, axis=axis, {center=center}, {solid...})

Rotation around the Z-axis (in trigonometric direction, i.e.
counter-clockwise).
"""
@inline rotate(θ::Real, s...; axis=nothing, kwargs...) =
	mult_matrix(rotation(float(θ), axis), s...; kwargs...)
"""
    rotate((θ,φ,ψ), {center=center}, {solid...})

Rotation given by Euler angles (ZYX; same ordering as OpenSCAD).
"""
@inline rotate(angles, s...; kwargs...) =
	mult_matrix(Rotations.RotZYX(radians.(angles)...), s...; kwargs...)

# project and cut««2
@inline project(s...) = operator(Project, (), s...)
@inline cut(s...) = operator(Cut, (), s...)

# overloading Julia operators««1
# backslash replaces U+2216 ∖ SET MINUS, which is not an allowed Julia operator
@inline Base.:\(x::AbstractGeometry, y::AbstractGeometry) = setdiff(x, y)
@inline Base.:-(x::AbstractGeometry, y::AbstractGeometry,
	tail::AbstractGeometry...) = setdiff(x, union(y, tail...))
# this purposely does not define a method for -(x::AbstractGeometry).
# (complement could be defined as ~x)
@inline Base.:-(x::AbstractGeometry{D}) where{D} = difference(intersect(), x)
@inline Base.:-(x::AbstractVector{<:AbstractGeometry},
                y::AbstractVector{<:AbstractGeometry}) = difference(x, y)

@inline Base.:+(v::AbstractVector, x::AbstractGeometry) = translate(v, x)
@inline Base.:+(x::AbstractGeometry, v::AbstractVector) = translate(v, x)
@inline Base.:-(x::AbstractGeometry, v::AbstractVector) = translate(-v, x)

@inline Base.:*(c::Real, x::AbstractGeometry) = scale(c, x)
@inline Base.:*(c::AbstractVector, x::AbstractGeometry) = scale(c, x)

⋃ = Base.union
⋂ = Base.intersect
#————————————————————— I/O —————————————————————————————— ««1

# File inclusion««1

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

# STL ««1
#
# SVG ««1

#————————————————————— Code to rewrite: —————————————————————————————— ««1
#————————————————————— Meshing (2d) —————————————————————————————— ««1

#»»1
# # 2d meshing««1
# # Offset and draw««2
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
#————————————————————— Extra tools —————————————————————————————— ««1
#»»1
# OpenSCAD output
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
