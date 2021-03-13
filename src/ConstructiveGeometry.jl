module ConstructiveGeometry
# Module imports««1
# using Printf
using LinearAlgebra
using StaticArrays
using FixedPointNumbers
using SparseArrays
using Logging

import Polyhedra # for convex hull
import GLPK
module AbstractMeshes
	using Meshes
# 	import Meshes: Geometry, Primitive, Point, Vec
# 	import Meshes: Box
	HyperSphere = Meshes.Sphere
	export Geometry, Point, Vec, Segment, Triangle, HyperSphere
	export coordtype, coordinates, embeddim, vertices
end; using .AbstractMeshes
import .AbstractMeshes: vertices, coordtype

# using MiniQhull
import Rotations
import Colors: Colors, Colorant
import Clipper
import LightGraphs

import Base: show, print
import Base: length, getindex, size, iterate, keys, eltype, convert
import Base: union, intersect, setdiff, copy, isempty, merge
import Base: *, +, -, ∈, inv, sign, iszero

#————————————————————— Ideal objects —————————————————————————————— ««1
#»»1
# Types««1
# Numeric types ««2

"""
    ConstructiveGeometry._FIXED
The type used whenever a fixed-precision real number is needed (e.g.
when interfacing with `Clipper.jl`).
"""
const _FIXED = Fixed{Int64,16}
"""
    ConstructiveGeometry._REAL

The default type used for computing coordinates (i.e. the type to which
integers are converted). Defaults to a `Fixed` type (fixed-precision
real), but this module should work just as well with floats.
"""
# const _REAL = _FIXED
const _REAL = Float64
"""
    real_type(T)

The type of real numbers to which objects of type `T` are promoted.
"""
@inline real_type(T::Type{<:Real}) = T
@inline real_type(::Type{<:Integer}) = _REAL
@inline real_type(a::Type{<:Real}, b::Type{<:Real}...) =
	real_type(promote_type(a, b...))
@inline real_type(a::Real...) = real_type(typeof.(a)...)
@inline to_real(x::T) where{T} = convert(real_type(T), x)

# divide by two without losing type information
@inline one_half(x::Real) = x/2
@inline one_half(x::T) where{T<:Fixed} = reinterpret(T, reinterpret(x)>>1)
@inline one_half(x::Integer=1) = one_half(to_real(x))
@inline one_half(x::AbstractArray) = one_half.(x)

# by default, sqrt(::Fixed) is a Float. We do not want this.
# This is (disappointingly...) faster than implementing a custom
# integer square-root function:
@inline Base.sqrt(a::T) where{T<:FixedPoint} = T(Base.sqrt(Float64(a)))
@inline Base.rem2pi(a::T, r) where{T<:FixedPoint} =
	T(Base.rem2pi(Float64(a), r))

# Dict of lists
DictOfLists{A,B} = Dict{A,Vector{B}}
@inline function listpush!(d::DictOfLists, (key, val)) # listpush!(l, k=>v)
	!haskey(d, key) && (d[key] = [])
	push!(d[key], val)
end
	
# Array indices ««2
"""
    cyclindex(i, a, [count=1])

Returns the index following `i` in the indices of array `a`, looping past
the end of `a` if needed.
If `count` is provided, advance by `count` steps instead.
"""
@inline cyclindex(i::Int, a, count...) = _cyclindex(i, eachindex(a), count...)
@inline _cyclindex(i::Int, a::Base.OneTo) = mod(i, a.stop)+1
@inline _cyclindex(i::Int, a::Base.OneTo, count::Int) =
	mod1(i+count, a.stop)

# Static vectors and paths««2
# functions missing in Meshes.jl««
(::Type{Point{N}})(a::AbstractVector) where{N} = Point{N,eltype(a)}(a)
Base.getindex(p::Point, i::Integer) = getindex(coordinates(p), i)
Base.getindex(p::Point, i::AbstractVector) = Point(getindex(coordinates(p), i))
Base.iszero(p::Point) = iszero(coordinates(p))
#»»
coordtype(v::AbstractVector{<:Number}) = eltype(v)
@inline coordinates(p::Point) = p.coords
@inline function barycenter(points::Point...)
	T = promote_type(coordtype.(points)...)
	return Point(sum(points) / convert(real_type(T), length(points)))
end
# """
#     Vec{D,T}
# The type used for representing `D`-dimensional vectors (or points) with
# coordinates of type `T`. An alias for `SVector{D,T}`.
# """
# const Vec{D,T} = SVector{D,T} # type alias
# this comes with the following constructors:
# Vec{2}(SA[1,2])		Vec(SA[1,2])
# Vec{2}((1,2))			Vec((1,2))
# Vec{2}(1,2)				Vec(1,2)
# Vec{2}([1,2])			- (needs explicit size) -

# XXX: we might need an “embedding” function (maybe not a constructor though):
# Vec{3,T}(x::AnyVec{2,T}, fill = zero(T)) where{T} = 
#		Vec{3,T}(get(x, i, fill) for i in 1:3)
#
# we need to extend get to tuples for the next function.
# this will fail if wrong index, as intended.
# Base.get(t::NTuple{D,T}, i, fill) where{D,T} = t[i]
# Vec{D,T}(v::AnyVec{D,T}, fill = zero(T)) where{D,T} =
#		Vec{D,T}(get(v, i, fill) for i in 1:D)

# for the path type, three possibilities:
#  (1) HybridArray{Tuple{D,StaticArrays.Dynamic()}, T}
#  (2) Vector{SVector{D,T}}
#  (3) a trivial wrapper for (1) or (2)
#
# (3) is probably lots of code for nothing (unless we someday need to
# distinguish between different possible uses for a Vector{SVector{D,T}},
# which is possible). (1) would be nice, except for the oversight in
# HybridArrays that makes a matrix product (SMatrix * HybridMatrix) a
# (fully-dynamic) Matrix. So we went for (2).
#
# XXX We might someday need a path type which uses something else than
# Vector (e.g. a generator, for computed paths).

const Path{D,T<:Real} = Vector{Point{D,T}}
# constructors:
# the following functions are wrappers, which should allow us to
# rewrite `Path` objects as matrices instead (instead of `length` and
# `getindex`, return one matrix dimension and one column:)
# See also `apply` below.

# This represents plausible user input for a `Path` object:
const AnyList{T} = Union{AbstractVector{<:T},NTuple{N,<:T} where{N}}
const AnyPath{D,T<:Real} = AnyList{<:Point{D,T}}

norm²(v::Vec) = sum(v .* v)
distance²(p::Point, q::Point) = norm²(p-q)


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
isunique(array) = length(unique(array)) == length(array)
struct Consecutives{T,V} <: AbstractVector{T}
	parent::V
end
consecutives(v::AbstractVector{T}) where{T} =
	Consecutives{T,typeof(v)}(v)
Base.getindex(c::Consecutives, i::Integer) =
	(c.parent[i], c.parent[mod1(i+1, length(c.parent))])
Base.size(c::Consecutives) = size(c.parent)
# findextrema««2
"""
    findextrema(itr; lt=isless)

Like `findmin`, except that
 - it returns both extrema, as a `NamedTuple`;
 - it accepts an `lt` parameter, like `sort`.
"""
function findextrema(itr; lt=isless)
  p = pairs(itr); y = iterate(p)
  if y == nothing
    throw(ArgumentError("collection must be non-empty"))
  end
  (mi, m), s = y; (Mi, M) = (mi, m); i = mi
  while true
    y = iterate(p, s)
    y == nothing && break
    (ai, a), s = y
    if lt(a, m) m = a; mi = ai; end
    if lt(M, a) M = a; Mi = ai; end
  end
  return (min=(m, mi), max=(M, Mi))
end

# small determinants««2

# 2-dimensional determinant, useful for computing orientation
@inline det2(v1, v2) = v1[1]*v2[2] - v1[2]*v2[1]
@inline det2(pt1, pt2, pt3) = det2(pt2-pt1, pt3-pt1)

# 3-dimensional determinant
@inline det3(v1::Vec{3}, v2::Vec{3}, v3::Vec{3}) = det([v1 v2 v3])
@inline det3(p1::Point{3}, p2::Point{3}, p3::Point{3}, p4::Point{3}) =
	det3(p2-p1, p3-p1, p4-p1)
# @inline det3(v1, v2, v3) = det([v1 v2 v3])
# @inline det3(p1, p2, p3, p4) = det3(p2-p1, p3-p1, p4-p1)
# Rows view««2
struct ViewRows{T,M<:AbstractMatrix{T}} <:
		AbstractVector{SubArray{T,1,M,Tuple{T,Base.Slice{Base.OneTo{T}}},true}}
	source::M
	ViewRows(m::AbstractMatrix) = new{eltype(m), typeof(m)}(m)
end
Base.size(r::ViewRows) = (size(r.source,1),)
Base.getindex(r::ViewRows, i::Integer) = view(r.source, i, :)

# Primitive solids««1
# Base type««2
"""
    PrimitiveSolid{S,D,T,X}

A type used to factorize code for the various primitive solids.
`D` and `T` are as usual the dimension and coordinate type.
`S` is a symbol indicating the type of solid (`:square`, etc.);
`X` is a NamedTuple type holding the parameters for this solid
(e.g. `@NamedTuple{radius::T}` for a sphere).

From the info in `X` we can derive a default constructor
(using keyword arguments matching the `NamedTuple`) and a
`show` method.

The coordinate type `T` is guessed from the values of the tuple using the
`_infer_type(H;kwargs...)` method.
"""
struct PrimitiveSolid{S,D,T,X} <: Geometry{D,T}
	parameters::X
	@inline PrimitiveSolid{S,D,T,X}(;kwargs...) where{S,D,T,X} =
		new{S,D,T,X}(kwargs.data)
end
@inline (H::Type{<:PrimitiveSolid{S,D,T}})(;kwargs...) where{S,D,T} =
	H{typeof(kwargs.data)}(kwargs.data)
@inline (H::Type{<:PrimitiveSolid{S,D}})(;kwargs...) where{S,D} =
	H{_infer_type(H;kwargs...)}(;kwargs...)

@inline scad_name(::PrimitiveSolid{S}) where{S} = S
@inline parameters(s::PrimitiveSolid) = getfield(s, :parameters)
@inline Base.getproperty(s::PrimitiveSolid, name::Symbol) =
	Base.getproperty(parameters(s), name)

# Glue for types imported from Meshes.jl: Square, Cube, Circle, Sphere««2
# --mesh
Ortho = AbstractMeshes.Box
@inline Ortho{D}(p1::AbstractVector, p2::AbstractVector) where{D} =
	Ortho{D,promote_type(eltype.((p1,p2))...)}(Point{D}(p1), Point{D}(p2))
@inline Ortho{D}(v::AbstractVector; origin=zero(v), center=false) where{D} =
	let p1 = center ? origin - one_half(v) : origin
	Ortho{D}(p1, p1 + v)
	end
@inline width(b::Ortho) = b.max - b.min

"""
    square(size; origin, center=false)

An axis-parallel square or rectangle  with given `size`
(scalar or vector of length 2).
"""
square(args...; kwargs...) = Square(args...; kwargs...)
square(a::Number, b::Number; kwargs...) = square(SA[a,b]; kwargs...)
square(a::Number; kwargs...) = square(a, a; kwargs...)
Square = Ortho{2}

"""
    cube(size; origin, center=false)

A cube or parallelepiped  with given `size`
(scalar or vector of length 3).
"""
cube(args...; kwargs...) = Cube(args...; kwargs...)
cube(a::Number, b::Number, c::Number; kwargs...) = cube(SA[a,b,c]; kwargs...)
cube(a::Number; kwargs...) = cube(a, a, a; kwargs...)
Cube = Ortho{3}

@inline HyperSphere{D}(r::T) where{D,T} =
	HyperSphere(Point(zero(SVector{D,T})), r)

"""
    circle(radius)

A circle. Discretization is done via the `accuracy` and `precision`
parameters.
"""
circle(args...; kwargs...) = Circle(args...; kwargs...)
Circle = HyperSphere{2}

"""
    sphere(radius)

A sphere. Discretization is done via the `accuracy` and `precision`
parameters.
"""
sphere(args...; kwargs...) = Sphere(args...; kwargs...)
Sphere = HyperSphere{3}

# Cylinder ««2
"""
    cylinder(h, r1, r2 [, center=false])
    cylinder(h, (r1, r2) [, center=false])
    cylinder(h, r [, center=false])

**Warning:** `cylinder(h,r)` is interpreted as `cylinder(h,r,r)`,
not `(h,r,0)` as in OpenSCAD.
"""
@inline cylinder(args...; kwargs...) = Cylinder(args...; kwargs...)
struct Cylinder{T} <: AbstractMeshes.Primitive{3,T}
	origin::Point{3,T}
	height::T
	r1::T
	r2::T
end
@inline Cylinder{T}(h, r1, r2; center=false, origin=zero(Point{3,T})) where{T} =
	Cylinder{T}(center ? origin-SA[0,0,one_half(h)] : origin, h, r1, r2)
@inline Cylinder(args...; kwargs...) =
	Cylinder{real_type(args...)}(args...; kwargs...)
@inline Cylinder(h, r; kwargs...) = Cylinder(h, r, r; kwargs...)


# Polygon ««2
"""
    polygon{T}
    polygon([point1, point2, ...])
    polygon(point1, point2, ...)

A simple, closed polygon enclosed by the given vertices.
"""
@inline polygon(args...; kwargs...) = Polygon(args...; kwargs...)
struct Polygon{T} <: Geometry{2,T}
	points::Vector{Point{2,T}}
	@inline Polygon(points::AbstractVector{<:Point{2}}) =
		new{real_type(coordtype.(points)...)}(points)
	@inline Polygon(points::AbstractVector{<:AbstractVector{<:Real}}) =
		Polygon(Point.(points))
	@inline Polygon(points::AbstractVector{<:Real}...) = Polygon([points...])
	# points in rows, for easy notation:
	@inline Polygon(m::AbstractMatrix{<:Real}) =
		Polygon([m[i,:] for i in 1:size(m,1)])
end
# Polygon{T} = PrimitiveSolid{:polygon,2,T,
# 	@NamedTuple{points::Path{2,T}}}
# @inline _infer_type(::Type{<:Polygon}; points) = real_type(eltype.(points)...)
# @inline (T::Type{<:Polygon})(points::AnyPath{2}) = T(points=points)
# @inline (T::Type{<:Polygon})(points::AnyVec{2,<:Real}...) = T([points...])

@inline vertices(p::Polygon) = p.points
@inline nvertices(p) = length(vertices(p))
# # Region (as explicit union of polygons-with-holes)««2
# """
#     Cheese
# 
# A polygonal area with polygonal holes.
# All these polygons are assumed to be simple closed loops, in direct order.
# """
# struct Cheese{T} <: Geometry{2,T}
#   exterior::Polygon{T}
#   holes::Vector{Polygon{T}}
# end
# 
# Cheese(p::Polygon) = Cheese{coordtype(p)}(p, [])
# Cheese(p::Polygon, holes::Polygon...) =
# 	Cheese{real_type(coordtype(p), coordtype.(holes)...)}(p, [holes...])
# 
# @inline scad_name(::Cheese) = :polygon
# function scad_parameters(p::Cheese)
# 	firstindex = similar(p.holes, Int)
# 	if length(p.holes) > 1 # else length(firstindex) == 0...
# 		firstindex[1] = nvertices(p.exterior)
# 		for i in 1:length(p.holes)-1
# 			firstindex[i+1] = firstindex[i] + nvertices(p.holes[i])
# 		end
# 	end
# 	points = vcat(vertices(p.exterior), vertices.(p.holes)...)
# 	paths = [[0:nvertices(p.exterior)-1;],
# 		[[firstindex[i] : firstindex[i]+nvertices(p.holes[i])-1;]
# 			for i in eachindex(p.holes)]...]
# 	return (points=points, paths=paths)
# end
# 
# """
#     Region
# 
# A distinct union of polygonal areas with polygonal holes.
# 
# Should be fillable with parity rule.
# """
# struct Region{T} <: Geometry{2,T}
#   children::Vector{Cheese{T}}
# end
# Region(p::Cheese) = Region{coordtype(p)}([p])
# Region(p::Polygon) = Region(Cheese(p))
# Region(v::Vector{<:Point{2}}) = Region(Polygon(v))
# 
# children(r::Region) = r.children
# @inline scad_name(::Region) = :union
# @inline scad_parameters(::Region) = NamedTuple()
# 
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

macro define_neutral(op, what, result)
	quote
	@inline $(esc(op))(neutral, absorb::$what) = $result
	@inline $(esc(op))(absorb::$what, neutral) = $result
	@inline $(esc(op))(x::$what, ::$what) = x
	end
end
union() = EmptyUnion()
intersect() = EmptyIntersect()
Base.show(io::IO, ::EmptyUnion) = print(io, "union()")
Base.show(io::IO, ::EmptyIntersect) = print(io, "intersect())")

# these are necessary for the following macros:
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

# # Somewhat reduce type I/O clutter««1
# Base.show(io::IO, ::Type{_FIXED}) = print(io, "_FIXED")
# Mesh types in 2d and 3d««1
# Draw ««2
"""
    draw(path, width; kwargs...)
    ends=:round|:square|:butt|:loop
    join=:round|:miter|:square
"""
@inline draw(args...; kwargs...) = Draw(args...; kwargs...)
struct Draw{T} <: Geometry{2,T}
	path::Vector{Point{2,T}}
	width::Float64
  ends::Symbol
	join::Symbol
	miter_limit::Float64
end
Draw(path::Union{AnyPath{2},AbstractVector{<:AbstractVector{<:Real}}}, width;
		ends=:round, join=:round, miter_limit=2.) =
	Draw{real_type(coordtype.(path)...)}(path, width, ends, join, miter_limit)

"""
    draw(path, width; kwargs)
    ends = :loop|:butt|:square|:round
		join = :square|:round|:miter
    miter_limit = 2.0

Draws a path of given width.
"""
draw(path, width; kwargs...) = Draw(path, width; kwargs...)

# PolygonXor««2
"""
    PolygonXor(polygon...)
    PolygonXor([points], [points]...)

Exclusive union of several polygons.
(This may be used to represent either the disjoint union of several polygons,
or a polygon with holes, etc.)
"""
struct PolygonXor{T} <: Geometry{2,T}
	paths::Vector{Polygon{T}}
end
PolygonXor(p::Polygon...) = PolygonXor{real_type(coordtype.(p)...)}([p...])
PolygonXor(points::AbstractVector{<:Point{2}}...) =
	PolygonXor(Polygon.(points)...)
PolygonXor(points::AbstractVector{<:AbstractVector{<:Real}}...) =
	PolygonXor([Point{2}.(p) for p in points]...)

@inline paths(p::PolygonXor) = p.paths
@inline vertices(p::PolygonXor) = [vertices.(paths(p))...;]
"""
    perimeters(::PolygonXor)

Returns a list of perimeters and holes for this region.
Perimeters are oriented ↺ and holes ↻.
"""
function perimeters(p::PolygonXor)
	firstindex = zeros(Int, length(p.paths))
	firstindex[1] = 0
	for i in 1:length(p.paths)-1
		firstindex[i+1] = firstindex[i] + nvertices(p.paths[i])
	end
	return [[firstindex[i]+1 : firstindex[i] + nvertices(p.paths[i]);]
		for i in eachindex(p.paths)]
end

scad_name(::PolygonXor) = :polygon
scad_parameters(p::PolygonXor) =
	(points = vertices(p),
	paths = [ f .- 1 for f in perimeters(p) ])

# Surface««2
# TODO also store some face-neighbour information
# this can be updated (not fully recomputed) for subtriangulation
"""
    AbstractSurface{T}

Abstract supertype for all triangulated surfaces.

### Parameters
 - `T` is the coordinate type

### Interface
 - `vertices(s)`: `AbstractDictionary` of indices
 - `faces(s)`: `AbstractDictionary` of faces

### TODO
Implement various concrete subtypes, e.g. surface with AABB tree,
or with extra incidence information.
"""
AbstractSurface = AbstractMeshes.Mesh{3}

@inline vertices(s::AbstractSurface) = _vertices(surface(s))
@inline faces(s::AbstractSurface) = _faces(surface(s))
@inline nvertices(s::AbstractSurface) = length(vertices(s))
@inline nfaces(s::AbstractSurface) = length(faces(s))


"""
    surface([points...], [faces...])

Encodes information about a surface.
"""
@inline surface(args...; kwargs...) = Surface(args...; kwargs...)
struct Surface{T} <: AbstractSurface{T}
	vertices::Vector{Point{3,T}}
	faces::Vector{SVector{3,Int}}
	Surface{T}(vertices, faces) where{T} = new(vertices, faces)
end

TriangulatedSurface = Surface # temporary alias

(T::Type{<:Surface})(points, faces) =
	T(Point{3}.(points), SVector{3,Int}.(faces))

Surface(points::AbstractVector{<:Point{3}},
		faces::AbstractVector{<:SVector{3,<:Integer}}, args...) =
	Surface{real_type(coordtype.(points)...)}(points, faces, args...)
Surface(points, faces) = triangulate(points, faces)

@inline surface(s::Surface) = s
@inline _vertices(s::Surface) = s.vertices
@inline _faces(s::Surface) = s.faces

# a lot of functions operating on 'Surface' values are defined later in
# the meshing part of this file.

#ConstructedSolid««1
# https://www.usenix.org/legacy/event/usenix05/tech/freenix/full_papers/kirsch/kirsch.pdf
"""
		ConstructedSolid{D,S}

A type representing CSG operations on solids. `D` is the dimension and
`S` is a symbol representing the operation (union, intersection etc.)
"""
struct ConstructedSolid{S,V,D,T} <: Geometry{D,T}
	children::V # Vector{<:AbstractGeometry}, or tuple etc.
	# passing a vector or tuple:
	@inline ConstructedSolid{S,V,D,T}(v::V) where{S,V,D,T} = new{S,V,D,T}(v)
end
@inline ConstructedSolid{S,V,D}(s) where{S,V,D} =
	ConstructedSolid{S,V,D,real_type(coordtype.(s)...)}(s)
@inline children(s::ConstructedSolid) = s.children
@inline scad_name(::ConstructedSolid{S}) where{S} = S

ConstructedSolid(s::Symbol, T = Vector{<:Geometry}) = ConstructedSolid{s,T}
CSGUnion = ConstructedSolid(:union)
CSGInter = ConstructedSolid(:intersection)
CSGDiff = ConstructedSolid(:difference,Tuple{<:Geometry,<:Geometry})
CSGComplement = ConstructedSolid(:complement,Tuple{<:Geometry})
CSGHull = ConstructedSolid(:hull)
CSGMinkowski = ConstructedSolid(:minkowski)
# @inline scad_name(::ConstructedSolid{D, :intersection}) where{D} = :intersection

# make operators associative; see definition of + in operators.jl
for op in (:union, :intersect, :minkowski, :hull)
	Q=QuoteNode(op)
	# union, intersection, minkowski are trivial on single objects:
	op != :hull &&  @eval ($op)(a::Geometry) = a
	@eval begin
	# all of these are associative:
	# we leave out the binary case, which will be defined on a case-by-case
	# basis depending on the operators (see below).
#		($op)(a::Geometry, b::Geometry) =
#			ConstructedSolid{$Q}([unroll(a, Val($Q)); unroll(b, Val($Q))])
	($op)(a::Geometry, b::Geometry, c::Geometry, x...) =
		Base.afoldl($op, ($op)(($op)(a,b),c), x...)
	end
end

"""
    union(s::Geometry...)
    s1 ∪ s2

Represents the union of given solids.
"""
@inline union(a1::Geometry, a2::Geometry) =
	CSGUnion{maximum(embeddim.((a1,a2)))}(unroll2(a1, a2, Val(:union)))
"""
    intersect(s::Geometry...)
    s1 ∩ s2

Represents the intersection of given solids.
"""
@inline intersect(a1::Geometry, a2::Geometry) =
	CSGInter{minimum(embeddim.((a1,a2)))}(unroll2(a1, a2, Val(:intersection)))
"""
    minkowski(s::Geometry...)

Represents the Minkowski sum of given solids.
"""
@inline minkowski(a1::Geometry, a2::Geometry) =
	CSGMinkowski{maximum(embeddim.((a1,a2)))}(unroll2(a1, a2, Val(:minkowski)))
"""
    hull(s::Geometry...)

Represents the convex hull of given solids.
"""
@inline hull(s::Geometry...) =
	CSGHull{maximum(embeddim.(s))}(
		[unroll(t, Val.((:hull, :union))...) for t in s])

"""
		unroll(x::Geometry, Val(sym1), Val(sym2)...)

Returns either `[x]` or, if `x` is a `ConstructedSolid` matching one of the
symbols `sym1`, `sym2`..., `children(x)`. (This helps reduce nesting).
"""
@inline unroll(s::Geometry, ::Val, tail...) = unroll(s, tail...)
@inline unroll(s::Geometry) = s
@inline unroll(s::ConstructedSolid{D, S}, ::Val{S}, tail...) where{D, S} =
	children(s)
@inline unroll2(s::Geometry, t::Geometry, tail...) =
	[unroll(s, tail...); unroll(t, tail...)]

# minus operator is binary:
"""
    ConstructiveGeometry.difference(s1, s2)

Represents the difference `s1 ∖ s2`.
"""
@inline difference(x::Geometry, y::Geometry) =
	CSGDiff{embeddim(x)}((x, y))
Base.:\(x::Geometry, y::Geometry) = difference(x, y)
# 		[unroll(x, Val(:difference)); unroll.(y, Val(:union))...])
# added interface: difference([x...], [y...])
@inline difference(x::AbstractVector{<:Geometry},
				y::AbstractVector{<:Geometry}) =
	difference(union(x...), union(y...))
# 	ConstructedSolid{foldr(max,embeddim.(x);init=2),:difference}(union(x...), y...)

@inline complement(x::Geometry{D}) where{D} =
	CSGComplement{embeddim(x)}((x,))

# General transforms««1
# Curry««2
"""
    Curry{S}

A structure representing partially-evaluated functions.
This allows chaining transformations by overloading the multiplication
operator: each factor in such a 'product', except the last one,
is a `Curry` object.

`S` is a datum indicating the type of transformation performed by the
function. It is used to compose functions when possible.

# Examples
```jldoctest
julia> add(a)=Curry(x->x+a)
julia> add(1)*add(2)*4
7
```
"""
struct Curry{S}
  f # either a Function or a Type...
end
# poor man's associative functor...

# fall-back case:
@inline Base.:*(f::Curry) = f
# binary rules:
@inline Base.:*(f::Curry, g::Curry) = compose(f, g)
@inline Base.:*(f::Curry, x) = f.f(x)
# ternary rule for associativity: we use the `assoc` type trait to
# decide whether to associate left or right.
@inline Base.:*(f::Curry, g::Curry, args...) =
	_comp(Val(assoc(f,g)), f, g, args...)
@inline _comp(::Val{:left} , f, g, args...) = *(compose(f, g), args...)
@inline _comp(::Val{:right}, f, g, args...) = *(f, *(g, args...))

# default values for the traits: transforms are right-associative and
# composition is trivial.
@inline assoc(::Curry, ::Curry) = :right
@inline compose(f::Curry, g::Curry) = Curry{:∘}(f.f ∘ g.f)

# Transform type««2
"""
    Transform{S,D,T,X}

Represents a solid of dimension `D` obtained via a transformation with
name `S` (a symbol).

This type defines functions allowing to chain transforms; these are used by
`multmatrix`, `color` etc. operations (see below).

The minimal job left to concrete types (see e.g. `AffineTransform` as an
example) is to define a type and a constructor:
    Frobnicate = Transform{:frobnicate}
		frobnicate(x::real, s...) = Frobnicate((x=x,), s...)
"""
struct Transform{S,D,D1,T,X} <: Geometry{D,T}
	data::X
	child::Geometry{D1,T}
	Transform{S,D,D1}(data, child::Geometry{D1}) where{S,D,D1} =
		new{S,D,D1,coordtype(child),typeof(data)}(data, child)
	Transform{S,D}(data, child::Geometry) where{S,D} =
		Transform{S,D,embeddim(child)}(data, child)
	# Default case: D = D1
	Transform{S}(data, child::Geometry) where{S} =
		Transform{S,embeddim(child)}(data, child)
end
# more constructors, including unary curryfied constructor:
@inline (T::Type{<:Transform{S}})(f, s1::Geometry,
		s2::Geometry, tail::Geometry...) where{S} =
		T(f, union(s, s2, tail...))
@inline (T::Type{<:Transform{S}})(f, s::Vector{<:Geometry}) where{S} =
	T(f, s...)
@inline (T::Type{<:Transform{S}})(f) where{S} = Curry{S}((s...)->T(f, s...))
# We can extract the `f` value from the above in the following way:
"""
    extract(c::Curry)

Given a `Curry` object with function `s -> Transform{...}(f, s)`,
recovers the parameter `f`.
"""
function extract end
@inline (T::Type{<:Transform})(f, ::typeof(extract)) = f
@inline extract(c::Curry) = c.f(extract)

# default values for I/O:
# (parameters in `data` are assumed to be stored in a NamedTuple).
@inline children(f::Transform) = [f.child]
@inline scad_name(f::Transform{S}) where{S} = S
@inline parameters(f::Transform) = f.data

# SetParameters««2
SetParameters = Transform{:parameters}
"""
    set_parameters(;accuracy, precision, symmetry) * solid...

A transformation which passes down the specified parameter values to its
child. Roughly similar to setting `\$fs` and `\$fa` in OpenSCAD.
"""
@inline set_parameters(s...; parameters...) =
	SetParameters(parameters.data, s...)

# Color««2
Color = Transform{:color}

"""
    color(c::Colorant, s...)
    color(c::AbstractString, s...)
    color(c::AbstractString, α::Real, s...)
    color(c) * s...

Colors objects `s...` in the given color.
"""
@inline color(c::Colorant, s...) = Color((color=c,), s...)
@inline color(c::AbstractString, s...) =
	color(parse(Colorant, c), s...)
@inline color(c::AbstractString, a::Real, s...) =
	color(Colors.coloralpha(parse(Colorant, c), a), s...)

# Linear extrusion««2
LinearExtrude = Transform{:linear_extrude,3,2}
"""
    linear_extrude(h, s...)
    linear_extrude(h) * s...

Linear extrusion to height `h`.
"""
@inline linear_extrude(h, scale::AbstractVector, s...; center=false)=
	LinearExtrude((height=h, scale=scale, center=center,), s...)
@inline linear_extrude(h, scale::Real, s...; kwargs...) =
	linear_extrude(h, SA[scale, scale], s...; kwargs...)
@inline linear_extrude(h, s...; kwargs...) =
	linear_extrude(h, 1, s...; kwargs...)

# Rotational extrusion««2
"""
    rotate_extrude([angle = 360°], solid...)
    rotate_extrude([angle = 360°]) * solid

Similar to OpenSCAD's `rotate_extrude` primitive.
"""
@inline rotate_extrude(s...) = rotate_extrude(360, s...)
@inline rotate_extrude(angle::Real, s...) =
	RotateExtrude((angle=angle,), s...)
RotateExtrude = Transform{:rotate_extrude,3,2}
# Offset
"""
    offset(r, solid...; kwargs...)
    offset(r; kwargs...) * solid

Offsets by given radius.

    ends=:round|:square|:butt|:loop
    join=:round|:miter|:square
"""
@inline offset(r::Real, s...; join=:round, miter_limit=2.) =
	Offset((r=r, join=join, miter_limit=miter_limit), s...)
Offset = Transform{:offset}
@inline scad_parameters(io::IO, s::Offset) =
	scad_parameters(io, s, Val(parameters(s).join), parameters(s))
@inline scad_parameters(io::IO, ::Offset, ::Val{:round}, param) =
	scad_parameters(io, (r=param.r,))
@inline scad_parameters(io::IO, ::Offset, ::Val{:miter}, param) =
	scad_parameters(io, (delta=param.r, chamfer=false,))
@inline scad_parameters(io::IO, ::Offset, ::Val{:square}, param) =
	scad_parameters(io, (delta=param.r, chamfer=true,))

# Affine transforms««1
# Affine type««2
# type and constructors««3
"""
    Affine(a, b)
		Affine(a)
		Affine(a, center=c)

A structure representing an affine transformation `x -> a*x + b`.
This is purposely kept as generic as possible.

As a special case, `a == Val(true)` corresponds to translations,
while `b == Val(false)` corresponds to linear maps. (The `+` and `*`
operators are overloaded correspondingly).
"""
struct Affine{A,B}
	a::A
	b::B
end
# default constructors: Affine{A,B}(::A,::B), Affine(::A,::B)
# @inline Affine(a; center) = Affine(a, a*center - center)
@inline Affine(a;center=Val(false)) = Affine(a, Val(false))

@inline (f::Affine)(p::Point) = Point(f.a * coordinates(p) + f.b)
function (f::Affine)(p::Polygon)
	if sign(f) > 0
		return Polygon(f.(vertices(p)))
	else
		return Polygon(reverse(f.(vertices(p))))
	end
end
(f::Affine)(p::PolygonXor) = PolygonXor(f.(paths(p)))

# @inline apply(f::Affine, p::Point) = Point(f.a * v.coords + f.b)
# @inline apply(f::Affine, points::AbstractVector{<:Point}) =
# 	[apply(f, p) for p in points]

@inline sign(f::Affine{<:Number}) = sign(f.a)
@inline sign(f::Affine{<:AbstractMatrix}) = sign(det(f.a))

# neutral elements: ««3
# this could in principle be defined for Val{T} where{T}, but we try
# to pirate a minimum number of functions in Base.
@inline Base.:*(::Val{true}, v) = v
@inline Base.:*(a, v::Val{false}) = v
@inline Base.:+(v, ::Val{false}) = v
@inline Base.:-(v, ::Val{false}) = v

# I/O: ««3
# OpenSCAD only uses 4×4 matrices for transformations;
# we pad the matrix to this size if needed:
function scad_parameters(io::IO, f::Affine)
	m = [ mat33(f.a) vec3(f.b); 0 0 0 1 ]
	print(io, "[")
	join(io, map(i->Float64.(view(m,i,:)),1:size(m,1)),",")
	print(io, "]")
end

@inline mat33(a::AbstractMatrix) = [ get(a, (i,j), i==j) for i=1:3, j=1:3 ]
@inline mat33(a::Diagonal) = Diagonal(vec3(diag(a)))
@inline mat33(a::Real) = SDiagonal(a,a,a)
@inline mat33(::Val{b}) where{b} = mat33(b)

@inline vec3(b::AbstractVector) = SVector{3}(get(b, i, 0) for i in 1:3)
@inline vec3(::Val{b}) where{b} = SA[b, b, b]

# Reflections««2
"""
		Reflection

A type containing compressed information for an orthogonal reflection
matrix. Inspired by the `Diagonal` type.
"""
struct Reflection{D,T,V<:AbstractVector} <: AbstractMatrix{T}
	axis::V
	@inline Reflection{D,T,V}(v::V) where {D,T,V<:AbstractVector{T}} =
		new{D,T,V}(v)
end
# we use only static vectors in Reflections:
@inline Reflection(v::AbstractVector) = Reflection(SVector{length(v)}(v))
# FIXME: add some method where we know that v is normed already
@inline Reflection(v::Vec{D}) where{D} = let u = v/norm(v)
	Reflection{D,typeof(u[1]),typeof(u)}(u)
end
# Mat{3}(R::Reflection{2}) = Reflection([R.axis;SA[0]])

@inline size(R::Reflection{D}) where{D} = (D,D)
@inline getindex(R::Reflection{D}, i::Int) where{D} =
	getindex(R, fld1(i,D), mod1(i,D))
@inline getindex(R::Reflection, i::Int, j::Int) = let a = R.axis
	(i == j) - 2*a[i]*a[j]
end
function Matrix{T}(R::Reflection{D,T}) where{D,T}
	# R(x) = x - 2<x, axis>/‖axis‖² · axis
	a = R.axis
	I(D) - 2*a*a'
end

# AffineTransform««2
AffineTransform = Transform{:multmatrix}

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

# these two functions are now enough to pre-compose all affine transforms
# *before* applying them to objects:
@inline assoc(::Curry{:multmatrix}, ::Curry{:multmatrix}) = :left
@inline function compose(c1::Curry{:multmatrix}, c2::Curry{:multmatrix})
	(f1, f2) = (extract(c1), extract(c2))
	mult_matrix(f1.a*f2.a, f1.a*f2.b + f1.b)
end

# Translation, scaling, rotation, mirror««2
# FIXME change this '1' to a compile-time constant?
"""
    translate(v, s...)
    translate(v) * s

Translates solids `s...` by vector `v`.
"""
@inline translate(v::AbstractVector, s...) = mult_matrix(1, v, s...)
"""
    scale(a, s...; center=0)
    scale(a; center=0) * s
Scales solids `s` by factor `a`. If `center` is given then this will be
the invariant point.

`a` may also be a vector, in which case coordinates will be multiplied by
the associated diagonal matrix.
"""
@inline scale(a::Real, s...; kwargs...) = mult_matrix(a, s...; kwargs...)
@inline scale(a::AbstractVector, s...; kwargs...) =
	mult_matrix(Diagonal(a), s...; kwargs...)
"""
    mirror(v, s...; center=0)
    mirror(v; center=0) * s

Reflection with axis given by the hyperplane normal to `v`.
If `center` is given, then the affine hyperplane through this point will
be used.
"""
@inline mirror(v::AbstractVector, s...; kwargs...) =
	mult_matrix(Reflection(v), s...; kwargs...)

@inline rotation(θ::AnyAngle; axis=SA[0,0,1], kwargs...) =
	real_type(θ,axis...).(Rotations.AngleAxis(radians(θ), axis...))
@inline rotation(θ::AnyList{<:AnyAngle}; kwargs...) =
	real_type(θ,axis...).(Rotations.RotZYX(radians.(θ); kwargs...))

"""
    rotate(θ, {center=center}, {solid...})
    rotate(θ, axis=axis, {center=center}, {solid...})

Rotation around the Z-axis (in trigonometric direction, i.e.
counter-clockwise).
"""
@inline rotate(θ, s...; kwargs...) = mult_matrix(rotation(θ; kwargs...), s...)
"""
    rotate((θ,φ,ψ), {center=center}, {solid...})

Rotation given by Euler angles (ZYX; same ordering as OpenSCAD).
"""
@inline rotate(θ::Real, φ::Real, ψ::Real, s...; kwargs...) =
	mult_matrix(rotation((θ,φ,ψ); kwargs...), s...)

# Operators««2
@inline +(v::AbstractVector, x::Geometry) = translate(v,x)
@inline +(x::Geometry, v::AbstractVector) = translate(v,x)

# this purposely does not define a method for -(x::Geometry).
@inline Base.:-(x::Geometry, y::Geometry, tail...) =
	difference(x, [y, tail...])
@inline Base.:-(x::Geometry{D}) where{D} = difference(intersect(), x)
@inline Base.:-(x::AbstractVector{<:Geometry},
                y::AbstractVector{<:Geometry}) = difference(x, y)

# @inline *(f::AbstractAffineMap, x::Geometry) = mult_matrix(f, x)
@inline *(s::Union{Real,AbstractVector}, x::Geometry) = scale(s,x)

⋃ = Base.union
⋂ = Base.intersect

# OpenSCAD output ««1
@inline indent(io::IO) = print(io, " "^get(io, :indent, 0))
@inline add_indent(io::IO, n=1) =
	IOContext(io, :indent => get(io, :indent, 0)+1)
@inline function Base.show(io::IO, l::Geometry...)
	for s in l
		indent(io); scad(io, s)
	end
end

@inline children(::Geometry) = nothing
@inline parameters(::Geometry) = NamedTuple()

"""
    scad(filename::AbstractString, s::Geometry...)
    scad(io::IO, s::Geometry)

Prints an OpenSCAD-like representation of the given solid(s).

## The various `scad_*` functions

    `scad_name(s)`
Returns, in printable form (e.g. `Symbol` or `String`), the OpenSCAD name
of this object.

    `scad_parameters(s)`
Returns a `NamedTuple` representing parameters of this object.

    `scad_transform(s)`
Possible transformation prepended to the object.

    `to_scad(x)`
Represents value `x` (number, array, etc.) in OpenSCAD format.

"""
function scad(io::IO, s::Geometry)
	indent(io)
	print(io, scad_transform(s))
	print(io, scad_name(s), "(")
	f = true;
	for (k, v) in pairs(scad_parameters(s))
		if f f = false; else print(io, ", "); end
		print(io, k, "=", to_scad(v))
	end
	print(io, ")")
	if children(s) isa Nothing
		println(io, ";")
	else
		scad_children(io, s)
	end
end

@inline scad(filename::AbstractString, s::Geometry...) =
	open(filename, "w") do f scad(f, s...) end

function scad_children(io::IO, s::Geometry)
	print(io, " {\n")
	io2 = add_indent(io)
	for c in children(s)
		scad(io2, c)
	end
	indent(io); print(io, "}")
end

@inline scad_name(::Square) = :square
@inline scad_name(::Cube) = :cube
@inline scad_name(::Circle) = :circle
@inline scad_name(::Sphere) = :sphere
@inline scad_name(::Cylinder) = :cylinder
@inline scad_name(::Polygon) = :polygon


@inline scad_parameters(s::Geometry) = parameters(s)
@inline scad_parameters(s::Ortho) = (size=Vector{Float64}(width(s)),)
@inline scad_parameters(s::HyperSphere) = (r=s.radius,)
@inline scad_parameters(s::Cylinder) = (h=s.height, r1=s.r1, r2=s.r2,)
@inline scad_parameters(p::Polygon) = (points=p.points,)

@inline scad_transform(s::Geometry) = ""
@inline scad_transform(s::Ortho) = scad_origin(minimum(s))
@inline scad_transform(s::HyperSphere) = scad_origin(s.center)
@inline scad_transform(s::Cylinder) = scad_origin(s.origin)

@inline scad_origin(p::Point) = scad_origin(coordinates(p))
@inline scad_origin(p::Vec) =
	iszero(p) ? "" : string("translate(", Vector{Float64}(p), ")")

@inline to_scad(x) = x
@inline to_scad(p::Point) = Vector{Float64}(coordinates(p))
@inline to_scad(p::Vec) = Vector{Float64}(p)
# kill any SVector{...} appearing in output:
@inline to_scad(v::AbstractVector) = Vector(to_scad.(v))
@inline to_scad(c::Colorant) = round.(Float64.([
	Colors.red(c), Colors.green(c), Colors.blue(c), Colors.alpha(c)]), digits=3)

# special case: Surface, with annotations for points
function scad(io::IO, s::AbstractSurface)
	println(io, "polyhedron(points=[ // ", nvertices(s), " points:")
	for (i,p) in pairs(vertices(s))
		indent(io)
		print(io, " ", Vector{Float64}(coordinates(p)))
		if i < nvertices(s) print(io, ","); end
		println(io, " // ", i)
	end
	println(io, "], faces=[ // ", nfaces(s), " faces:")
	for (i,f) in pairs(faces(s))
		indent(io)
		print(io, " ", Vector(f .- 1))
		if i < nfaces(s) print(io, ","); end
		println(io, " // ", i, "=", Vector(f))
	end
	indent(io); println(io, "] );")
end

function scad(io::IO, s::SetParameters)
	indent(io); println(io, "{ // SetParameters")
	for (i,j) in pairs(s.data)
		indent(io); println(io, "// ", i, "=", j, ";")
	end
	scad(io, s.child)
	indent(io); println(io, "} // SetParameters")
end

#————————————————————— Meshing (2d) —————————————————————————————— ««1

include("SpatialSorting.jl")

#»»1
# Generic code for 2d and 3d meshing««1
mesh(s::Geometry) = mesh(s, _DEFAULT_PARAMETERS)
# “thickness” of points, edges etc. for computing intersections:
const _THICKNESS = 1e-8

# Transformations««2
function mesh(s::AffineTransform{3}, parameters)
	g = mesh(s.child, parameters)
	b = sign(s.data)
	@assert b ≠ 0 "Only invertible linear transforms are supported (for now)"
	if b > 0
		return (typeof(g))(s.data.(vertices(g)), faces(g))
	else
		return (typeof(g))(s.data.(vertices(g)), reverse.(faces(g)))
	end
end
@inline mesh(s::SetParameters, parameters) =
	mesh(s.child, merge(parameters, s.data))
# Generic case (e.g. `color`): do nothing
@inline mesh(s::Transform, parameters) = mesh(s.child, parameters)

# Converting circles to polygons««1
# Accuracy is the absolute deviation allowed.
# Default value is 2.0 (from OpenSCAD `$fs`), interpreted as 2mm.
#
# Precision is the relative deviation allowed.
# Default value is 0.02 (1-cos(180°/`$fa`)).
# FIXME: explain why .005 works better

_DEFAULT_PARAMETERS = (accuracy = 0.1, precision = .005, symmetry = 1)

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
	reinterpret(Vec{2,T}, z)
end
@inline unit_n_gon(r, parameters::NamedTuple) =
	r*unit_n_gon(real_type(r), sides(r, parameters))

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

Returns a set of `n` well-spaced points, of type `Vec{3,T}`, on the unit
sphere.

TODO: use ideas from
http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/

to optimize for volume of convex hull.
"""
function fibonacci_sphere_points(T::Type{<:Real}, n::Int)

	v = Vector{Vec{3,T}}(undef, n)
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

# Clipper.jl interface: clip, offset, simplify««1
# This is the only section in this file which contains code directly
# related to `Clipper.jl`. The entry points to this section are the
# functions `clip` and `offset` defined below.
# Types««2
# default number of bits for Clipper types
const _CLIPPER_BITS = FixedPointNumbers.nbitsfrac(_FIXED)
# this must be a 64-bit type, even if _FIXED is modified:
const _CLIPPER_FIXED = Fixed{Int64,_CLIPPER_BITS}
# constants
const _CLIPPER_ENUM = (#««
	clip=(
		union       =Clipper.ClipTypeUnion,
		intersection=Clipper.ClipTypeIntersection,
		difference  =Clipper.ClipTypeDifference,
		xor         =Clipper.ClipTypeXor,
	),
	# the enum order in `Clipper.jl` is wrong, we fix it here:
	ends=(
		fill  = Clipper.EndType(0), #Clipper.EndTypeClosedPolygon,
		loop  = Clipper.EndType(1), #Clipper.EndTypeClosedLine,
		butt  = Clipper.EndType(2), #Clipper.EndTypeOpenButt,
		square= Clipper.EndType(3), #Clipper.EndTypeOpenSquare,
		round = Clipper.EndType(4), #Clipper.EndTypeOpenRound,
	),
	join=(
		square=Clipper.JoinTypeSquare,
		round =Clipper.JoinTypeRound,
		miter =Clipper.JoinTypeMiter,
	),
	fill=(
		nonzero=Clipper.PolyFillTypeNonZero,
		evenodd=Clipper.PolyFillTypeEvenOdd,
		positive=Clipper.PolyFillTypePositive,
	),
)#»»

"""
    to_clipper(OriginalType, ...)
		from_clipper(OriginalType, ...)

Converts stuff (numbers, vectors, paths...) to and from `Clipper.jl` types.
"""
@inline to_clipper(::Type{<:Real}) = _CLIPPER_FIXED
@inline to_clipper(T::Type{<:FixedPoint{<:Int64}}) = T
@inline to_clipper(T::Type, x::Real) = reinterpret(convert(to_clipper(T), x))
@inline to_clipper(T::Type, v::Point{2,<:Real}) =
	Clipper.IntPoint(to_clipper(T, v[1]), to_clipper(T, v[2]))
@inline to_clipper(T, p::AbstractVector{<:Point{2}}) =
	[to_clipper(T, v) for v in p]
@inline to_clipper(T, p::Vector{<:AbstractVector{<:Point{2}}}) =
	[to_clipper(T, v) for v in p]

# special case: if the data is already compatible, we just wrap it
@inline to_clipper(::Type{T}, p::Path{2,T}) where{T<:FixedPoint{<:Int64}} =
	GC.@preserve p unsafe_wrap(Array,
		pointer(reinterpret(Clipper.IntPoint, p)),
		length(p))
# FIXME: in the general case, check stride of array

# this is a special case; the `delta` parameter wants a Float64,
# and we *might* use a different type for `delta` than for coordinates:
@inline to_clipper_float(T::Type{<:Real}, x)::Float64 =
	x*FixedPointNumbers.rawone(to_clipper(T))

# numbers...
@inline from_clipper(T::Type{<:Real}, x::Int64) =
	convert(T, reinterpret(to_clipper(T), x))
# points...
@inline from_clipper(T::Type{<:Real}, p::Clipper.IntPoint) =
	Point{2}(SA[from_clipper(T, p.X), from_clipper(T, p.Y)])
# paths...
@inline from_clipper(T::Type{<:Real}, p::Vector{Clipper.IntPoint}) =
	[ from_clipper(T, v) for v in p ]
@inline from_clipper(T::Type{<:Fixed{Int64}}, p::Vector{Clipper.IntPoint}) =
	reinterpret(Vec{2,T}, p)
# vectors of paths...
@inline from_clipper(T, polys::Vector{Vector{Clipper.IntPoint}}) =
	[ from_clipper(T, p) for p in polys ]
# Wrappers for Clipper calls««2
# We wrap all Clipper objects in a NamedTuple with the original type
struct Marked{T,X}
	data::X
	@inline Marked{T}(x::X) where{T,X} = new{T,X}(x)
end

@inline ClipperClip(T::Type) = Marked{T}(Clipper.Clip())
@inline ClipperOffset(T::Type, miterLimit::Real, roundPrecision::Real) =
	Marked{T}(Clipper.ClipperOffset(Float64(miterLimit),
		to_clipper_float(T, roundPrecision)))

@inline add_path!(c::Marked{T}, path::AnyPath{2}, args...) where{T} =
	Clipper.add_path!(c.data, to_clipper(T, path), args...)
@inline add_paths!(c::Marked{T}, paths::Vector{<:AnyPath{2}}, args...) where{T}=
	Clipper.add_paths!(c.data, [ to_clipper(T, p) for p in paths], args...)

@inline execute(c::Marked{T,Clipper.Clip}, args...) where{T} =
	from_clipper(T, Clipper.execute(c.data, args...)[2])
@inline execute(c::Marked{T,Clipper.ClipperOffset}, delta::Real) where{T} =
	from_clipper(T, Clipper.execute(c.data, to_clipper_float(T, delta)))

# Calls on Path values««2
function clip(op::Symbol,
		v1::AbstractVector{Path{2,T}},
		v2::AbstractVector{Path{2,T}};
		fill = :evenodd)::Vector{Path{2,T}} where {T}
	c = ClipperClip(T)
	add_paths!(c, v1, Clipper.PolyTypeSubject, true) # closed=true
	add_paths!(c, v2, Clipper.PolyTypeClip, true)

	f = _CLIPPER_ENUM.fill[fill]
	return execute(c, _CLIPPER_ENUM.clip[op], f, f)
end
function offset(v::AbstractVector{Path{2,T}}, r::Real;
		join = :round,
		ends = :fill,
		miter_limit = 2.,
		precision = 0.2
		)::Vector{Path{2,T}} where{T}
	c = ClipperOffset(T, miter_limit, precision)
	add_paths!(c, v, _CLIPPER_ENUM.join[join], _CLIPPER_ENUM.ends[ends])
	execute(c, r)
end
function offset(v::AbstractVector{Path{2,T}}, r::AbstractVector{<:Real};
		join = :round,
		ends = :fill,
		miter_limit = 2.,
		precision = 0.2
		)::Vector{Vector{Path{2,T}}} where{T}
	# “Simultaneously” computes offset for several offset values.
	# Used by path_extrude() below.
	c = ClipperOffset(T, miter_limit, precision)
	add_paths!(c, v, _CLIPPER_ENUM.join[join], _CLIPPER_ENUM.ends[ends])
	[ execute(c, ρ) for ρ in r]
end
@inline function simplify(p::Vector{<:AnyPath{2,T}}; fill=:nonzero) where{T}
	return from_clipper(T,
		Clipper.simplify_polygons(to_clipper(T, p), _CLIPPER_ENUM.fill[fill]))
end
"""
    orientation(p::Path{2})

Returns `true` iff p is a direct loop (i.e. if area >= 0).
"""
@inline function orientation(p::Path{2,T}) where{T}
	return Clipper.orientation(to_clipper(T, p))
end
"""
    pointinpolygon(pt::Vec{2}, p::Path{2})

Returns 1 if point is in the interior, -1 on boundary, and 0 outside the
given polygon.

Polygon is assumed not self-intersecting.
"""
@inline function point_in_polygon(point::Point{2,T},
	path::AnyPath{2,T}) where{T}
	return Clipper.pointinpolygon(to_clipper(T, point), to_clipper(T, path))
end
@inline point_in_polygon(point::Point{2}, path::Polygon) =
	point_in_polygon(point, vertices(path))
# Polyhedra interface and intersections««1
# Polyhedra types««2
@inline polyhedra_lib(T::Type{<:Real}) =
	Polyhedra.DefaultLibrary{T}(GLPK.Optimizer)

# fixing an oversight in Polyhedra.jl: it has multiplications but no
# divisions
@inline Base.:(/)(h::Polyhedra.HyperPlane, α::Real) =
	Polyhedra.HyperPlane(h.a / α, h.β / α)

# converts path to matrix with points as rows:
@inline poly_vrep(points::AnyPath) = vcat(transpose.(Vector.(points))...)
@inline poly_vrep(points::Matrix) = points
@inline poly_eltype(points::AnyPath) = eltype(eltype(points))
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
@inline (h::Polyhedra.HRepElement)(p::Point) = h.a ⋅ coordinates(p) - h.β
@inline ∈(p::Point, h::Polyhedra.HyperPlane) = iszero(h(p))
@inline convert(T::Type{<:Polyhedra.HRepElement}, h::Polyhedra.HRepElement) =
	T(h.a, h.β)

# Intersections (2d)««2
# we need a couple of functions in the particular case of simplexes.
# `Polyhedra.jl` is a bit slow for these simple cases, so we write them
# here:

standardize(x::Float64) = (x == -0.0) ? 0.0 : x
"""
    inter(path, hyperplane::Polyhedra.HyperPlane)

intersection of simplex and hyperplane
"""
function inter(path::AnyPath, hyperplane::Polyhedra.HyperPlane)
	n = length(path)
	s = [hyperplane(p) for p in path]
	newpath = similar(path, n); c = 0
	for i in 1:n
		if s[i] == 0
			newpath[c+= 1] = path[i]
		end
		for j in 1:i-1
			if s[i] * s[j] < 0
				newpath[c+= 1] = standardize.((s[j]*path[i]-s[i]*path[j])/(s[j]-s[i]))
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
function inter(path::AnyPath, halfplane::Polyhedra.HalfSpace)
	s = [halfplane(p) for p in path]
	boundary = convert(Polyhedra.HyperPlane, halfplane)
	n = length(path)
	# we know that we add at most 1 new point (cutting a corner).
	newpath = similar(path, n+1); c = 0
	for i in eachindex(path)
		j = mod1(i+1, n)
		(si, sj) = (s[i], s[j])
		if si >= 0   newpath[c+=1] = path[i]; end
		if si*sj < 0
		# whiskers would generate two new points; we remove the second one
			newpoint = inter(path[[i,j]], boundary)[1]
			if (c==0|| newpath[c] != newpoint) newpath[c+=1] = newpoint; end
		end
	end
	return newpath[1:c]
end
@inline inter(path::AnyPath, h::Polyhedra.HRepElement,
		t::Polyhedra.HRepElement...) =
	inter(inter(path, h), t...)

"""
    line(p1=>p2)
"""
# XXX
function line(p12::Pair{<:Point{2}})
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
function halfplane(p12::Pair{<:Point{2}}, p3::Point{2})
	l = line(p12)
# 	(x1, y1) = p12[1]
# 	(x2, y2) = p12[2]
# 	a = SA[y1-y2, x2-x1]
# 	b = y1*x2 - x1*y2
	s = sign(l.a ⋅ coordinates(p3) - l.β)
	return Polyhedra.HalfSpace(s*l.a, s*l.β)
end

"""
    hrep(pt1, pt2, pt3)

Returns the triple of half-planes delimiting the interior of the given
triangle.
"""
function hrep(p1::Point{2}, p2::Point{2}, p3::Point{2})
	return (halfplane(p1=>p2, p3),
	        halfplane(p2=>p3, p1),
	        halfplane(p3=>p1, p2))
end

function line_inter(s1::Segment{2}, s2::Segment{2})
	((x1,y1), (x2,y2)) = coordinates.(vertices(s1))
	((x3,y3), (x4,y4)) = coordinates.(vertices(s2))
	d=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
	iszero(d) && return nothing
	a = x1*y2-y1*x2; b = x3*y4-y3*x4
	d1 = inv(to_real(d))
	return Point(standardize.(d1 .* SA[a*(x3-x4)-b*(x1-x2), a*(y3-y4)-b*(y1-y2)]))
end
function inter(s1::Segment{2}, s2::Segment{2}; thickness = 0)
	c = line_inter(s1, s2)
	(c == nothing) && return nothing
	((x1,y1), (x2,y2)) = coordinates.(vertices(s1))
	((x3,y3), (x4,y4)) = coordinates.(vertices(s2))
	(c[1] > x1 + thickness && c[1] > x2 + thickness) && return nothing
	(c[1] > x3 + thickness && c[1] > x4 + thickness) && return nothing
	(c[2] > y1 + thickness && c[2] > y2 + thickness) && return nothing
	(c[2] > y3 + thickness && c[2] > y4 + thickness) && return nothing
	(c[1] < x1 - thickness && c[1] < x2 - thickness) && return nothing
	(c[1] < x3 - thickness && c[1] < x4 - thickness) && return nothing
	(c[2] < y1 - thickness && c[2] < y2 - thickness) && return nothing
	(c[2] < y3 - thickness && c[2] < y4 - thickness) && return nothing
# 	if c ∉ boundingbox(s1) return nothing; end
# 	if c ∉ boundingbox(s2) return nothing; end
	return c
end
#Convex hull««1
# 2d convex hull ««2

# """
#     convex_hull([vector of 2d points])
# 
# Returns the convex hull (as a vector of 2d points).
# """
# function convex_hull(points::Union{AbstractVector{<:Vec{2}},
# 	AbstractMatrix{<:Number}})
# # M is a matrix with the points as *columns*, hence the transpose
# # below:
# 	PH = Polyhedra
# 	poly = vpoly(points)
# 	PH.removevredundancy!(poly)
# 	return Vec{2}.(collect(PH.points(poly)))
# end

"""
    convex_hull_list(points)

Returns the convex hull of the points, as a list of indexes (in direct
order, starting at a reproducible index in the list of points).
"""
function convex_hull_list(points)
  # Uses Graham scan
  # In practice this is faster than using `Polyhedra.jl`.
#   println("points=$points, length=$(length(points))")
  i0 = findextrema(points;
    lt=(p,q)->(p[2]<q[2])|| (p[2]==q[2] && p[1]>q[1])).min[2]
  @inline detp2(i,j,k) = det2(points[[i,j,k]]...)
	# 1024 is an arbitrary value for detecting “aligned” points (i.e. up to
	# representation errors), which should be fast for both Float and Fixed
	# division
  @inline function are_aligned(i,j,k)
    v1 = points[j]-points[i]
    v2 = points[k]-points[j]
    d = det2(v1, v2)
    c = v1 ⋅ v2
    return abs(d) < abs(c)/1024
  end
  scan = sort(filter(!isequal(i0), eachindex(points)),
    lt=(p,q)->detp2(i0,p,q) > 0)
#   println("i0=$i0, scan=$scan")
  stack = [i0, scan[1]]
  for h in scan[2:end]
#     println("scanning: $stack + $h")
    v1 = points[stack[end]] - points[stack[end-1]]
    v2 = points[h] - points[stack[end]]
    s = det2(v1, v2)
    c = v1 ⋅ v2
    if abs(s) < abs(c)/1024 && c < 0 # points are aligned and backwards
			# here we know that we can insert at (end)
			# look for an insertion point i:
			i = length(stack)
			while i > 2
# 				println(" try to insert at $i")
				v1 = points[stack[i]] - points[stack[i-1]]
				v2 = points[h] - points[stack[i]]
				s = det2(v1, v2)
				c = v1 ⋅ v2
				if s < -1e-3*abs(c)
# 					println(" break at $i")
					break
				end
				i -= 1
			end
# 			println("  inserting at $i")
			insert!(stack, i, h)
			continue
# 			println("  now stack=$stack")
    end
    while detp2(last(stack), h, stack[end-1]) < 0
      pop!(stack)
    end
    push!(stack, h)
  end
  return stack
end

"""
    convex_hull([vector of 2d points])

Returns the convex hull (as a vector of 2d points, ordered in direct
order).
"""
@inline convex_hull(points::AbstractVector{<:Point{2}}) =
	points[convex_hull_list(points)]


# this is the version using MiniQhull: #««
# convex_hull(points::AbstractVector{<:AnyVec(2)}) =
#		# Delaunay triangulation:
#		let T = MiniQhull.delaunay([p[i] for i in 1:2, p in points]),
#				N = length(points),
#				M = zeros(Bool, N, N) # ⚠ memory O(N²)
#		# mark all edges:
#		for (a,b,c) in eachcol(T)
#			b1 = points[b] - points[a]
#			c1 = points[c] - points[a]
#			d = b1[1]*c1[2] - b1[2]*c1[1] # determinant === orientation of triangle
#			if d < 0
#				M[a,b] = M[b,c] = M[c,a] = true
#			else
#				M[b,a] = M[c,b] = M[a,c] = true
#			end
#		end
#		# list of remaining edges (retrograde oriented)
#		L= sort([(i,j) for i in 1:N, j in 1:N if M[i,j] && ! M[j,i]], by=v->v[1])
#		next(i) = L[searchsorted(L, i, by=y->y[1])][1][2]
#		R = zeros(Int, length(L))
#		R[1:2] .= L[1] # initialize with first edge
#		for i in 3:length(R)
#			R[i] = next(R[i-1])
#		end
#		# returns in retrograde ordering (OpenSCAD convention):
#		points[R]
# end#»»
# 3d convex hull ««2

"""
    convex_hull(x::Geometry{3}...)

Returns the convex hull of the union of all the given solids, as a
pair `(points, faces)`. `faces` is a list of triangles.
"""
@inline convex_hull(x::Geometry{3}) =
	convex_hull(vcat([vertices(y) for y in x]...))

"""
    convex_hull(vector of 3d points)

Returns the convex hull of these points, as a pair `(points, faces)`.
All the faces are triangles.
"""
function convex_hull(p::AbstractVector{<:Point{3,T}}) where{T}
	M = hcat(Vector.(coordinates.(p))...)
	PH = Polyhedra
	poly = PH.polyhedron(PH.vrep(transpose(M)), polyhedra_lib(T))
	R = PH.removevredundancy!(poly)
	V = Point{3,T}.(collect(PH.points(poly)))

	triangles = Vec{3,Int}[]
	for i in PH.eachindex(PH.halfspaces(poly)) # index of halfspace
		h = PH.get(poly, i)
		pts = PH.incidentpointindices(poly, i) # vector of indices of points
		for t in triangulate_face(
				[Point(PH.get(poly, j)) for j in pts];
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

# 2d Minkowski sum««1
# Identify holes of PolygonXor
"""
    identify_polygons(s::PolygonXor)

Given a `PolygonXor` defined as the exclusive union of polygons and holes,
returns a vector mapping individual polygons of `s` to:
 - positive indices 1, 2, 3... for the exterior polygons of `s`;
 - negative indices for the holes, indicating in which exterior polygon is the corresponding hole.
"""
function identify_polygons(s::PolygonXor)
	n = length(paths(s))
	m = [ point_in_polygon(vertices(p)[1], q) for p in paths(s), q in paths(s) ]
	index = zeros(Int, n)
	hole = zeros(Int, n)
	c = 0 # number of connected components
	for i in 1:n
		# if this is a hole, then it will have m[i,j] > 0 for some i
		k = findfirst(m[i,:] .> 0)
		if k == nothing
			index[i] = (c+= 1)
		else
			hole[i] = k
		end
	end
	for i in 1:n
		if hole[i] > 0
			index[i] = -index[hole[i]]
		end
	end
	return index
end
# Convolution of polygons««2
# http://acg.cs.tau.ac.il/tau-members-area/general%20publications/m.sc.-theses/thesis-lienchapter.pdf
"""
    circularcmp(v1, v2, v3, [Val(:offset)])

Circular comparison predicate; returns true iff directions of vectors
`v1`, `v2`, `v3` are arranged in a trigonometric ordering along the unit
circle.

If `Val(:offset)` is passed then `v1`, `v3` are infinitesimally rotated
in the positive direction compared to `v2`.
"""
function circularcmp(v1, v2, v3)
	d1 = v2[1]*v3[2] ≥ v2[2]*v3[1]
	d2 = v3[1]*v1[2] ≥ v3[2]*v1[1]
	d3 = v1[1]*v2[2] ≥ v1[2]*v2[1]
	return (d1+d2+d3) ≥ 2
end
function circularcmp(v1, v2, v3, ::Val{:offset})
	d1 = v2[1]*v3[2] > v2[2]*v3[1]
	d2 = v3[1]*v1[2] ≥ v3[2]*v1[1]
	d3 = v1[1]*v2[2] ≥ v1[2]*v2[1]
	return (d1+d2+d3) ≥ 2
end

function convolution(p::AbstractVector{<:Vec{2}},
		q::AbstractVector{<:Vec{2}})
	(np, nq) = (length(p), length(q))
	ep = [p[cyclindex(i, p)]-p[i] for i in eachindex(p)] # edges of p
	eq = [q[cyclindex(i, q)]-q[i] for i in eachindex(q)]
	j0 = 0
	newpoly = similar(p, 0)
	for ip in eachindex(p)
		for iq in eachindex(q)
			iq0 = cyclindex(iq, q, -1)
			if circularcmp(eq[iq0], ep[ip], eq[iq], Val(:offset))
				push!(newpoly, p[ip]+q[iq])
				push!(newpoly, p[cyclindex(ip, p)]+q[iq])
			end
		end
	end
	newpoly
end
p⋆q = convolution(p, q)
# Minkowski sum of polygons and their unions ««2
function minkowski(p::AnyPath{2}, q::AnyPath{2}; fill=:nonzero)
	r = convolution(p, q)
	return simplify([r]; fill)
end
function minkowski(vp::AbstractVector{<:AbstractVector{<:Vec{2}}},
		vq::AbstractVector{<:AbstractVector{<:Vec{2}}}; fill=:nonzero)
	vr = vec([Point{2}.(convolution(p, q)) for p in vp, q in vq])
	return simplify(vr; fill)
end
function minkowski(p::Polygon, q::Polygon)
	return PolygonXor(minkowski(p, q; fill=:evenodd)...)
end
function minkowski(p::PolygonXor, q::PolygonXor)
	cp = [ coordinates.(vertices(x)) for x in paths(p) ]
	cq = [ coordinates.(vertices(y)) for y in paths(q) ]
	return PolygonXor(minkowski(cp, cq, fill=:evenodd)...)
end

# TODO: 3d Minkowski««2

# 2d meshing««1
# Primitive objects««2
function vertices(s::Square)
	# in trigonometric order:
	(u, v) = (minimum(s), maximum(s))
	return Point{2}.([
		SA[u[1],u[2]],
		SA[v[1],u[2]],
		SA[v[1],v[2]],
		SA[u[1],v[2]]])
end
@inline vertices(c::Circle, parameters) =
	[ c.center + p for p in unit_n_gon(c.radius, parameters) ]

mesh(s::PolygonXor, parameters) = s
mesh(s::Polygon, parameters) = PolygonXor(s)
mesh(s::Square, parameters) = PolygonXor(vertices(s))
mesh(s::Circle, parameters) = PolygonXor(vertices(s, parameters))

# Transforms««2
# Reduction of CSG operations««2
@inline clip(op, s::PolygonXor...) =
	reduce((p,q)->clip(op, vertices.(paths(p)), vertices.(paths(q)),
		fill=:evenodd), s)
mesh(s::CSGUnion{2}, parameters) =
	PolygonXor(clip(:union, [mesh(x, parameters) for x in children(s)]...)...)

mesh(s::CSGInter{2}, parameters) =
	PolygonXor(clip(:intersection,
		[mesh(x, parameters) for x in children(s)]...)...)

mesh(s::CSGDiff{2}, parameters) =
	PolygonXor(clip(:difference,
		mesh(s.children[1], parameters), mesh(s.children[2], parameters))...)

function mesh(s::CSGHull{2}, parameters)
	l = [mesh(x, parameters) for x in children(s)]
	return PolygonXor(convex_hull([vertices.(l)...;]))
end

function mesh(s::CSGMinkowski{2}, parameters)
	l = [mesh(x, parameters) for x in children(s)]
	global G = minkowski(l[1], l[2])
# 	return PolygonXor(reduce((p,q)->minkowski(p,q), l)...)
end

function mesh(s::AffineTransform{2}, parameters)
	g = mesh(s.child, parameters)
	f = s.data
	b = sign(f)
	@assert b ≠ 0 "Only invertible linear transforms are supported (for now)"
	return f(g) # reversal (if b < 0) is done polygon-by-polygon there
end
# Set-wise operations:
# # Minkowski sum:
# function (U::Type{<:PolyUnion})(s::ConstructedSolid{2,:minkowski},
# 	parameters)::U
# 	reduce((a,b)->U(minkowski(a.poly, b.poly)),
# 		_convert(U, s.children, parameters))
# end
# function _combine2(::Val{:minkowski}, a::PolyUnion{T}, b::PolyUnion{T}) where{T}
# 	# not implemented in Clipper.jl...
# end

# Offset and draw««2
function mesh(s::Offset, parameters)
	T = coordtype(s)
	m = mesh(s.child, parameters)
	ε = max(get_parameter(parameters,:accuracy), get_parameter(parameters,:precision) * s.data.r)
	return PolygonXor(offset(vertices.(paths(m)), s.data.r;
		join = s.data.join,
		ends = :fill,
		miter_limit = s.data.miter_limit,
		precision = ε)...)
end

function mesh(s::Draw, parameters)
	r = one_half(s.width)
	ε = max(get_parameter(parameters,:accuracy), get_parameter(parameters,:precision) * r)
	p = offset([s.path], r; join=s.join, ends=s.ends, miter_limit = s.miter_limit)
	return PolygonXor(p...)
end

# """
#     draw(path, width; kwargs...)
# 
#     ends=:round|:square|:butt|:closed
#     join=:round|:miter|:square
# """
# function draw(path::Path{2,T}, width::Real;
# 		ends::Symbol = :round, join::Symbol = :round,
# 		miter_limit::Float64 = 2.0, precision::Real = 0.2) where{T}
# 	CT = clipper_type(T)
# 	RT = clipper_rettype(T)
# 	c = ClipperOffset(miter_limit, clipper_float(CT, precision))
# 	println("join=$join, round=$round")
# 	Clipper.add_path!(c, clipper_path(path),
# 		JoinType(Val(join)), EndType(Val(ends)))
# 	println("$(clipper_type(T)) $(CT(1.)); prec=$(Float64(CT(precision)))")
# 	ret = clipper_unpath.(RT, Clipper.execute(c, clipper_float(CT, width)/2))
# 	return PolyUnion(ret)
# end
# 
# # Convex hull««2
# # """
# # 		convex_hull(x::Geometry{2}...)
# # 
# # Returns the convex hull of the union of all the given solids, as a
# # `PolyUnion` structure.
# # """
# @inline convex_hull(x::Geometry{2}...) =
# 	convex_hull(PolyUnion(union(x...)))
# 
# @inline convex_hull(u::PolyUnion) = convex_hull(Vec{2}.(vertices(u)))
# 

# 2d subsystem««1
# # PolyUnion««2
# # type and constructors from points««
# """
# 		PolyUnion
# 
# Represents a union of polygons. Each polygon is assumed to be simple and
# ordered in trigonometric ordering.
# """
# struct PolyUnion{T} <: Geometry{2,T}
# 	poly::Vector{Path{2,T}}
# 	@inline PolyUnion{T}(p::AbstractVector{<:AnyPath{2,T}}) where{T} =
# 		new{T}(Path{2,T}.(p))
# end
# @inline (U::Type{PolyUnion})(p::AbstractVector{<:AnyPath{2,T}}) where{T} =
# 		U{real_type(eltype.(eltype.(p))...)}(p)
# 
# @inline (U::Type{<:PolyUnion})(path::AnyPath{2}...) = U([path...])
# 
# @inline vertices(u::PolyUnion) = vcat(u.poly...)
# 
# # this is used to broadcast conversion for recursive conversion to PolyUnion:
# @inline _convert(U::Type{<:PolyUnion}, l, parameters) =
# 	[ U(s, parameters) for s in l ]
# 
# # »»
# # I/O««
# function scad(io::IO, u::PolyUnion{S}, spaces::AbstractString) where{S}
# 	print(io, spaces, "// union of $(length(u.poly)) polygon(s):\n")
# 	length(u.poly) != 1 && print(io, spaces, "union() {\n")
# 	for p in u.poly
# 		print(io, spaces, " polygon([")
# 		join(io, convert.(Vec{2,Float64}, p), ",")
# 		print(io, "]);\n")
# 	end
# 	length(u.poly) != 1 && print(io, spaces, "}\n")
# end
# #»»
# # Conversion from leaf 2d types««2
# @inline PolyUnion(l::Geometry{2}; kwargs...) =
# 	PolyUnion{real_type(eltype(l))}(l, merge(_DEFAULT_PARAMETERS, kwargs.data))
# 
# @inline (U::Type{<:PolyUnion})(x::Square, parameters) =
# 	U(Vec{2}.(vertices(x)))
# @inline (U::Type{<:PolyUnion})(x::Circle, parameters) =
# 	U(Vec{2}.(vertices(x, parameters)))
# @inline (U::Type{<:PolyUnion})(x::Polygon, parameters) =
# # FIXME: simplify and define orientation
# 	U(Vec{2}.(vertices(x)))
# 
# @inline function (U::Type{<:PolyUnion})(f::AffineTransform{2}, parameters)
# 	child = U(f.child, parameters)
# 	return U([ f.data.(path) for path in child.poly ])
# end
# @inline (U::Type{<:PolyUnion})(s::SetParameters{2}, parameters) =
# 	U(s.child, merge(parameters, s.data))
# # fall-back case (`color`, etc.):
# @inline (U::Type{<:PolyUnion})(s::Transform{S,2}, parameters) where{S} =
# 	U(s.child, parameters)
# 
# Reduction of CSG operations««2
# # Reduction of CSG operations««2
# @inline (clip(op, u::U...)::U) where{U<:PolyUnion} =
# 	reduce((a,b)->U(clip(op, a.poly, b.poly)), u)
# 
# # Set-wise operations:
# @inline (U::Type{<:PolyUnion})(s::ConstructedSolid{2,:union}, parameters) =
# 	clip(:union, _convert(U, s.children, parameters)...)
# 
# @inline (U::Type{<:PolyUnion})(s::ConstructedSolid{2,:intersection}, parameters) =
# 	clip(:intersection, _convert(U, s.children, parameters)...)
# 
# function ((U::Type{<: PolyUnion})(s::ConstructedSolid{2,:difference}, parameters)::U)
# 	length(s.children) == 1 && return U(s.children[1], parameters)
# 	L = _convert(U, s.children, parameters)
# 	r2= clip(:union, view(L,2:length(L))...)
# 	clip(:difference, L[1], r2)
# end
# 
# # Convex hull:
# function (U::Type{<:PolyUnion})(s::ConstructedSolid{2,:hull}, parameters)
# 	pts = points.(_convert(U, s.children, parameters))
# 	U(convex_hull([pts...;]))
# end
# 
# # Minkowski sum:
# function (U::Type{<:PolyUnion})(s::ConstructedSolid{2,:minkowski},
# 	parameters)::U
# 	reduce((a,b)->U(minkowski(a.poly, b.poly)),
# 		_convert(U, s.children, parameters))
# end
# # function _combine2(::Val{:minkowski}, a::PolyUnion{T}, b::PolyUnion{T}) where{T}
# # 	# not implemented in Clipper.jl...
# # end
# 
# 
# # Offset ««2
# """
# 		offset(P::Polygon, u::Real; options...)
# 
# Offsets polygon `P` by radius `u` (negative means inside the polygon,
# positive means outside). Options:
# 
#  - `join_type`: :round | :square | :miter
#  - `miter_limit` (default 2.0)
# """
# function offset(U::PolyUnion{T}, u::Real;
# 		join_type = :round,
# 		miter_limit::Float64 = 2.0,
# 		precision::Real = 0.2) where{T}
# 
# 	c = ClipperOffset(miter_limit, clipper_float(clipper_type(T), precision))
# 	add_paths!(c, U.poly, join_type, Clipper.EndTypeClosedPolygon)
# 	PolyUnion(execute(T, c, u))
# end
# @inline offset(x::Geometry{2}, args...; kwargs...) =
# 	offset(PolyUnion(x), args...; kwargs...)
# 
# # Draw ««2
# """
#     draw(path, width; kwargs...)
# 
#     ends=:round|:square|:butt|:closed
#     join=:round|:miter|:square
# """
# function draw(path::Path{2,T}, width::Real;
# 		ends::Symbol = :round, join::Symbol = :round,
# 		miter_limit::Float64 = 2.0, precision::Real = 0.2) where{T}
# 	CT = clipper_type(T)
# 	RT = clipper_rettype(T)
# 	c = ClipperOffset(miter_limit, clipper_float(CT, precision))
# 	println("join=$join, round=$round")
# 	Clipper.add_path!(c, clipper_path(path),
# 		JoinType(Val(join)), EndType(Val(ends)))
# 	println("$(clipper_type(T)) $(CT(1.)); prec=$(Float64(CT(precision)))")
# 	ret = clipper_unpath.(RT, Clipper.execute(c, clipper_float(CT, width)/2))
# 	return PolyUnion(ret)
# end
# 
# # Convex hull««2
# # """
# # 		convex_hull(x::Geometry{2}...)
# # 
# # Returns the convex hull of the union of all the given solids, as a
# # `PolyUnion` structure.
# # """
# @inline convex_hull(x::Geometry{2}...) =
# 	convex_hull(PolyUnion(union(x...)))
# 
# @inline convex_hull(u::PolyUnion) = convex_hull(Vec{2}.(vertices(u)))
# 
#————————————————————— Meshing (3d) —————————————————————————————— ««1
#=
mesh(geom, parameters = _DEFAULT_PARAMETERS) returns a Surface
mesh(x::Surface, parameters...) = x
=#
#»»1
# Basic 3d geometry««1
# Low-level functions««2
function face_normal(points)
	@assert length(points) >= 3
	return cross(points[2]-points[1], points[3]-points[1])
end
"""
    supporting_plane(t::Triangle)

Returns an equation (`a*x = b`) of the supporting plane, with `a`
pointing *outwards*.
"""
function supporting_plane(t::Triangle)
	(p1, p2, p3) = t.vertices
	c = cross(p2-p1, p3-p1)
	b = dot(c, p1.coords)
	return Polyhedra.HyperPlane(c, b)
end

"""
    circular_lt

Circular comparison of vectors, sorted according to their angle in
]-π,π]. Implemented with only integral arithmetic (no `atan2`, `√` or `/`).
"""
@inline function circular_lt(p,q)
	if p[2] < 0
		if q[2] ≥ 0 return true; end
	else # p[2] ≥ 0
		if q[2] < 0 return false; end
		if p[2] == 0 && q[2] == 0 return p[1] > q[1]; end
	end
	return det2(p, q) > 0
end
"""
    circular_sign(u,v)

Let `α` and `β` be the angles of `u`,`v` in [-π, π[.
This function returns a number <0 iff `α` < `β`, >0 iff `α` > `β`,
and `0` iff `α` == `β`.
"""
@inline function circular_sign(u, v)
# 16 cases to consider: u = (-1, -i, 1, i), same for v
	if u[2] > 0
		v[2] ≤ 0 && return -1 # (i,-1), (i,-i), (i,1)
	elseif u[2] < 0
		v[2] > 0 && return 1 #(-i,i)
	elseif u[2] == 0
		if v[2] == 0
			return sign(v[1]) - sign(u[1]) #(1,1) (1,-1) (-1,1) (-1,-1)
		elseif v[2] > 0
			return 1 #(-1,i) (1,i)
		else
			# the following line is not needed, but faster than computing det:
			return -sign(u[1]) #(-1,-i) (1,-i) 
		end
	end
	# determinant also works for the following cases:
	# (-1,-i), (-i, -1), (-i, 1), (1, i)
	return det2(u,v)
end



# Using a Box as a bounding box««2
BBox = AbstractMeshes.Box
@inline Base.min(b::BBox) = b.min; @inline Base.max(b::BBox) = b.max
@inline ∈(x::AbstractVector, b::BBox) = all(min(b) .≤ x .≤ max(b))
# @inline ∈(x::Point, b::BBox) = x.coords ∈ b
@inline isempty(b::BBox) = any(coordinates(min(b)) .> coordinates(max(b)))
@inline intersect(a::BBox, b::BBox) =
	BBox(Point(max.(coordinates(min(a)), coordinates(min(b)))),
	    Point(min.(coordinates(max(a)), coordinates(max(b)))))
@inline boundingbox(v::StaticVector{D,T}...) where{D,T} =
	BBox{D,T}(min.(v...), max.(v...))
@inline boundingbox(p::Point...) = boundingbox(coordinates.(p)...)
@inline boundingbox(g::Geometry) = boundingbox(AbstractMeshes.vertices(g)...)

# SpatialSorting interface:
@inline Base.min(p::Point...) = Point(min.(coordinates.(p)...))
@inline Base.max(p::Point...) = Point(max.(coordinates.(p)...))
@inline SpatialSorting.position(b::BBox) =
	coordinates(b.min) + coordinates(b.max)
@inline SpatialSorting.merge(b1::BBox, b2::BBox) =
	BBox(min(b1.min, b2.min), max(b1.max, b2.max))
# 3d -> 2d projections««2
const plus1mod3 = SA[2,3,1]
const plus2mod3 = SA[3,1,2]
@inline function project_2d(direction::AbstractVector, index::Val = Val(false))
	# we inline the `findmax` call since we know the length is 3:
	# (doing this about halves the running time of this function. Besides,
	# since the value `k` only takes three possible values, it enables the
	# compiler to propagate constants.)
	a1 = abs(direction[1]); a2=abs(direction[2]); a3=abs(direction[3])
	k = (a1 < a2) ? ((a2 < a3) ? 3 : 2) : ((a1 < a3) ? 3 : 1)
	v = direction[k]
	@inbounds p = (v > 0) ? SA[plus1mod3[k], plus2mod3[k]] :
		SA[plus2mod3[k], plus1mod3[k]]
	return _project_2d(index, p, k)
end
# @inline function project_2d1(direction::AnyVec{3}, index::Val = Val(false))
# 	# we inline the 'findmax' call since we know the length is 3:
# 	# (doing this about halves the running time of this function)
# 	a1 = abs(direction[1]); a2=abs(direction[2]); a3=abs(direction[3])
# 	# this part does not do any speed-up (constant propagation):
# 	if a1 < a2
# 		if a2 < a3 @goto max3; end
# 		if direction[2] > 0
# 			return _project_2d(index, SA[3,1], 2)
# 		else
# 			return _project_2d(index, SA[1,3], 2)
# 		end
# 	elseif a1 < a3
# 		@label max3
# 		if direction[3] > 0
# 			return _project_2d(index, SA[1,2], 3)
# 		else
# 			return _project_2d(index, SA[2,1], 3)
# 		end
# 	else
# 		if direction[1] > 0
# 			return _project_2d(index, SA[2,3], 1)
# 		else
# 			return _project_2d(index, SA[3,2], 1)
# 		end
# 	end
# end
@inline _project_2d(::Val{false}, p, _) = p
@inline _project_2d(::Val{true}, p, e) = (p, e)

"""
    project_2d(plane::Polyhedra.HyperPlane)

Returns a (named) tuple `(coordinates, linear, origin)` where
 - `coordinates` is the set of coordinates to keep for projection,
 - `linear`*x+`origin` is an affine section.
"""
function project_2d(plane::Polyhedra.HyperPlane)
	v = direction(plane)
  (coords, k) = project_2d(v, Val(true))
	f = inv(convert(real_type(eltype(v)), v[k]))
		# e=1: proj=(2,3) sect=[-b/a -c/a;1 0;0 1]
		# e=2: proj=(3,1) sect=[0 1;-a/b -c/b;1 0]
		# e=3: proj=(1,2) sect=[1 0;0 1;-a/c -b/c]
		# e=1: proj=(3,2) sect=[-c/a -b/a;0 1;1 0]
		# e=2: proj=(1,3) sect=[1 0;-a/b -c/b;0 1]
		# e=3: proj=(2,1) sect=[0 1;1 0;-b/c -a/c]
	m = SMatrix{3,2}((i==k) ? -f*v[coords[j]] : (i == coords[j])
			for i=1:3, j=1:2)
	c = SVector{3}((i==k) ? plane.β*f : 0 for i in 1:3)
	return (coordinates=coords, lift=Affine(m, c))
end

# Intersections (3d)««2
# function inter(s1::Segment{3}, s2::Segment{3}; thickness = 0)
# 	(a1, b1) = vertices(s1)
# 	(a2, b2) = vertices(s2)
# # 	d = det3(a1, b1, a2, b2)
# # 	@debug "inter: d=$d"
# # 	!iszero(d) && return nothing
# 	# compute supporting plane
# 	plane = supporting_plane(Triangle(a1, b1, a2))
# 	# check all points are coplanar
# 	abs(plane(b2)) > thickness && return nothing
# 	iszero(direction(plane)) && return a2
# 	(proj, lift) = project_2d(plane)
# 	int2 = inter(Segment(a1[proj],b1[proj]), Segment(a2[proj],b2[proj]);
# 		thickness)
# 	int2 == nothing && return nothing
# 	return lift(int2)
# end

"""
    inter(triangle1, triangle2; ε)

Returns a description of the intersection of two 3-dimensional triangles.
The description is given as a pair `(newpoints, newedges)`, where:
 * `newpoints[i] = (coordinates, position1, position2)`
  - `coordinates` is a `Point` object
  - `position1` is the position relative to triangle 1, encoded as:
    {1,2,3} for points on edges {(2,3),(3,1),(1,2)};
    4 for point inside the face; 0 for point outside the face;
  - `position2` is the same relative to triangle 2;
 * `newedges` is a list of added edges, represented as pairs (i,j),
   where -1,-2,-3 represent the vertices of t1;
   -4,-5,-6 the vertices of t2; and 1,2,3,... the new points.
"""
function inter(t1::Triangle{3}, t2::Triangle{3}; ε=_THICKNESS)#««
	# [Devillers, Guigue, _Faster triangle-triangle intersection tests_;
	#   https://hal.inria.fr/inria-00072100/document]
	# https://github.com/yusuketomoto/ofxCGAL/blob/master/libs/CGAL/include/CGAL/Triangle_3_Triangle_3_intersection.h
	# https://fossies.org/linux/CGAL/include/CGAL/Intersections_3/internal/Triangle_3_Triangle_3_do_intersect.h
	(p1, q1, r1) = vertices(t1)
	(p2, q2, r2) = vertices(t2)

	normal2 = cross(q2-p2, r2-p2)
	dp1 = dot(normal2, p1-p2)
	dq1 = dot(normal2, q1-p2)
	dr1 = dot(normal2, r1-p2)
	
	# test for coplanarity
	(dp1 == 0) && (dq1 == 0) && (dr1 == 0) && return inter_coplanar(t1, t2)

	# the permutation applied to points is stored in an integer
	# low part: {1,2,3} indicates which point is put first (3-cycle);
	# high part: +{0,3} indicates whether to then transpose two last points
	i1 = 1; i2 = 1
	# rotate triangle 1 as needed so that t2 separates p1 from (q1, r1)««
	if abs(dp1) ≤ ε
		if abs(dq1) ≤ ε
			if abs(dr1) ≤ ε return inter_coplanar(t1, t2) # 000
			elseif dr1 > 0  i2+=3     # 00+: change to 00-
			else            nothing   # 00-
			end
		elseif dq1 > 0
		  if dr1 < -ε     i1+=1     # 0+-: put q1 first
			else            i2+=3     # 0++ or 0+0: swap to 0-- or 0-0
			end
		else # dq1 < 0
			if dr1 > ε  i1+=2     # 0-+: put r1 first
			else            nothing   # 0-- or 0-0
			end
		end
	elseif dp1 > 0
		if abs(dq1) ≤ ε
			if dr1 > ε  i1+=1; i2+=3 # +0+: put q first and swap
			else            nothing   # +0- or +00: do nothing
			end
		elseif dq1 > 0
			if dr1 > ε  return ((),()) # +++, no intersection
			else        i1+= 2; i2+= 3 # ++- or ++0: put r1 first; swap q2, r2
			end
		else
			if dr1 > ε i1+=1; i2+= 3  # +-+: q1 first, swap q2, r2
			else                      # +-- or +-0: do nothing
			end
		end
	else
	  if abs(dq1) ≤ ε
			if dr1 > ε  i1+= 2 # -0+: put r first
			else        i1+= 1 # -0- or -00: put q first
			end
		elseif dq1 < 0
			if dr1 < -ε return ((),()) # ---, no intersection
			else        i1+= 2         # --+ or --0: r1 first
			end
		else
			if dr1 < 0 i1+= 1         # -+-: q1 first
			else       i2+= 3         # -++: swap p2, q2
			end
		end
	end #»»
	# likewise for second triangle  ««
	normal1 = cross(q1-p1, r1-p1)
	dp2 = dot(normal1, p2-p1)
	dq2 = dot(normal1, q2-p1)
	dr2 = dot(normal1, r2-p1)
	if dp2 > 0
		if dq2 > 0
			if dr2 > 0 return ((),()) # +++, no intersection
			else i2+= 2; i1+= 3       # ++-: put r2 first; swap q1, r1
			end
		else
			if dr2 > 0 i2+=1; i1+= 3  # +-+: q2 first, swap q1, r1
			else                      # +--: do nothing
			end
		end
	else
		if dq2 < 0
			if dr2 < 0 return ((),()) # ---, no intersection
			else       i2+= 2         # --+: r2 first
			end
		else
			if dr2 < 0 i2+= 1         # -+-: q2 first
			else       i1+= 3         # -++: swap p1, q1
			end
		end
	end #»»
	# apply both permutations ««
	perm_table = ((1,2,3),(2,3,1),(3,1,2),(1,3,2),(2,1,3),(3,2,1))
	@inbounds σ1 = perm_table[i1]
	@inbounds σ2 = perm_table[i2]

	@inbounds (a1, b1, c1) = (vertices(t1)[i] for i in σ1)
	@inbounds (a2, b2, c2) = (vertices(t2)[i] for i in σ2)
	# the permutations σ1 and σ2 will allow us to recover the real indices
	# of points (edges, etc):
	# e.g. assume a1 = q1, which means σ1(1) = 2
	# if the algorithm returns a value x=1 (meaning a1),
	# then indexing σ1[x] returns 2 (meaning this is really q1).

	# re-use already computed determinants as z-coordinates:
	@inbounds (za1, zb1, zc1) = ((dp1, dq1, dr1)[i] for i in σ1)
	@inbounds (za2, zb2, zc2) = ((dp2, dq2, dr2)[i] for i in σ2)
	
	# 1,2,3 represent even permutations, and 4,5,6 odd ones:
	(i1 > 3) && (normal1 = -normal1; za2 = -za2; zb2 = -zb2; zc2 = -zc2)
	(i2 > 3) && (normal2 = -normal2; za1 = -za1; zb1 = -zb1; zc1 = -zc1)

	# a laundry list of assertions to check that we are in a standard
	# configuration:
	@assert normal1 ≈ cross(b1-a1, c1-a1)
	@assert normal2 ≈ cross(b2-a2, c2-a2)
	@assert dot(normal2, a1-a2) ≈ za1
	@assert dot(normal2, b1-a2) ≈ zb1
	@assert dot(normal2, c1-a2) ≈ zc1
	@assert dot(normal1, a2-a1) ≈ za2
	@assert dot(normal1, b2-a1) ≈ zb2
	@assert dot(normal1, c2-a1) ≈ zc2
	@assert za1 ≥ 0
	@assert zb1 ≤ 0
	@assert zc1 ≤ 0
	@assert za2 ≥ 0
	@assert zb2 ≤ 0
	@assert zc2 ≤ 0
	# »»
	# coordinates of four intersection points bb1, cc1, bb2, cc2
	# (all four are aligned on the intersection of the two planes)
	@inline barycenter(p1, p2, λ) = p2 + λ*(p1-p2) # λp1+(1-λ)p2
	bb1 = barycenter(b1, a1, za1/(za1-zb1))
	cc1 = barycenter(c1, a1, za1/(za1-zc1))
	bb2 = barycenter(b2, a2, za2/(za2-zb2))
	cc2 = barycenter(c2, a2, za2/(za2-zc2))
	@assert abs(dot(normal2, bb1-a2)) ≤ 1e-10
	@assert abs(dot(normal2, cc1-a2)) ≤ 1e-10
	@assert abs(dot(normal1, bb2-a1)) ≤ 1e-10
	@assert abs(dot(normal1, cc2-a1)) ≤ 1e-10

	# relative position of these points on the line:
	# project on best coordinate, i.e. largest abs. value of direction
	# and flip sign if needed so that bb1 is before cc1:
	v = cc1 - bb1
	e = abs(v[1]) > abs(v[2]) ? (abs(v[1]) > abs(v[3]) ? 1 : 3) :
	    abs(v[2]) > abs(v[3]) ? 2 : 3
	# and map b1 to 0, c1 to 1:
	α = 1/(cc1[e]-bb1[e])
	xb2 = α*(bb2[e]-bb1[e]); xc2 = α*(cc2[e]-bb1[e])

	# we know _a priori_ that the ordering of (bb1, cc1) is opposed to that
	# of (bb2, cc2):
	@assert xc2 ≤ xb2
	@inline aligned(p,q,r) = norm(cross(q-p, r-p)) < ε
	@inline monotonic(p,q,r) = dot(p-q, r-q) < 0
	@inline bb1() = barycenter(cc2, bb2, -xb2/(xc2-xb2))
	@inline cc1() = barycenter(cc2, bb2, (1-xb2)/(xc2-xb2))

	if xc2 < 0-ε
		if xb2 < 0-ε return ((),())
		elseif xb2 ≤ 0+ε # bb2 ∈ edge [a1, b1]
			return (
				((cc2, 0, σ2[2]), (bb2, σ1[3], σ2[3])),
				((1,2),))
		elseif xb2 < 1-ε # [cc2, bb2] intersects [a1, b1]
			return (
				((cc2, 0, σ2[2]), (bb1(), σ1[3], 4), (bb2, 4, σ2[3])),
				((1,2), (2,3),))
		elseif xb2 ≤ 1+ε # bb2 ∈ edge [a1, c1]
			return (
				((cc2, 0, σ2[2]), (bb1(), σ1[3], 4), (bb2, 2, σ2[3])),
				((1,2), (2,3),))
		else               # [cc2, bb2] intersects [a1, b1] and [a1, c1]
			return (
				((bb1(), σ1[3], 4), (cc1(), σ1[2], 4)),
				((1,2),))
		end
	elseif xc2 ≤ 0+ε   # cc2 ∈ edge [a1, b1]...
		if xb2 < 1-ε     # ... and is the single intersection point
			return (
				((cc2, σ1[3], σ2[2]), (bb2, 4, σ2[3])),
				((1,2),))
		elseif xb2 ≤ 1+ε # ... and bb2 ∈ [a1, c1]
			return (
				((cc2, σ1[3], σ2[2]), (bb2, σ1[2], σ2[3])),
				((1,2),))
		else               # ... and [cc2, bb2] intersects [a1, c1]
			return (
				((cc2, σ1[3], σ2[2]), (cc1(), σ1[2], 4), (bb2, 0, σ2[3])),
				((1,2),(2,3),))
		end
	elseif xc2 < 1-ε
		if xb2 < 1-ε     # [cc2, bb2] is an inner segment of triangle t1
			return (
				((cc2, 4, σ2[2]), (bb2, 4, σ2[2])),
				((1,2),))
		elseif xb2 ≤ 1+ε # bb2 ∈ edge [a1, c1]
			return (
				((cc2, 4, σ2[2]), (bb2, σ1[2], σ2[3])),
				((1,2),))
		else               # [cc2, bb2] intersects [a1, c1]
			return (
				((cc2, 4, σ2[2]), (cc1(), σ1[2], 4), (bb2, 0, σ2[3])),
				((1,2),(2,3),))
		end
	elseif xc2 ≤ 1+ε   # cc2 ∈ edge [a1, c1]
		return (
			((cc2, σ1[2], σ2[2]), (bb2, 0, σ2[3])),
			((1,2),(2,3),))
	else # 1 < xc2 ≤ xb2, no intersection
		return ((),())
	end
end#»»


# Rays««2
"""
    AffineRay

An open ray, of the form `a + ]0,∞[ v`.
"""
struct AffineRay{D,T} <: Geometry{D,T}
	origin::Point{D,T}
	direction::Vec{D,T}
end

@inline direction(a::AffineRay) = a.direction
@inline origin(a::AffineRay) = a.origin
@inline AffineRay(origin::SVector{D}, direction::SVector{D}) where{D} =
	AffineRay{D,promote_type(eltype.((origin,direction))...)}(origin, direction)

"""
    intersects(a::AffineRay{3}, t::Triangle)

Returns 1 iff ray intersects triangle in given order,
-1 if it intersects in opposite order, otherwise 0.

Algorithm is inspired by Segura-Feito[1]:
after translating everything so that the ray starts at the origin,
we apply the linear transformation
x ↦ (⟨x,v2,v3⟩,⟨v1,x,v3⟩,⟨v1,v2,x⟩)  (⟨⟩ is the determinant)
This maps the three points of the triangle to
(δ,0,0), (0,δ,0), (0,0,δ), where δ=⟨v1,v2,v3⟩.
The sign of δ gives the position of the origin w.r.t the triangle: δ>0
iff the origin is *inside* the triangle.

Once the vertices of the triangle are aligned to the axes, the ray
intersects the triangle iff the 3 coordinates of its direction `u` are >0.
This means that ⟨u,v2,v3⟩, ⟨v2,u,v3⟩, ⟨v1,v2,u⟩ > 0.

[1] https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.2.2084&rep=rep1&type=pdf

"""
function intersects(a::AffineRay, t::Triangle, strict = false)
	u = a.direction
	(v1, v2, v3) = [x-a.origin for x ∈ vertices(t)]
	d = det3(v1, v2, v3)
	# FIXME: could be make slightly faster (but less readable)
	# by precomputing u∧v2
	if strict
		if (d > 0 && det3(u, v2, v3) > 0
			&& det3(v1, u, v3) > 0 && det3(v1, v2, u) > 0)
			return 1
		elseif (d < 0 && det3(u, v2, v3) < 0
			&& det3(v1, u, v3) < 0 && det3(v1, v2, u) < 0)
			return -1
		end
	else
		if (d ≥ 0 && det3(u, v2, v3) ≥ 0
			&& det3(v1, u, v3) ≥ 0 && det3(v1, v2, u) ≥ 0)
			return 1
		elseif (d ≤ 0 && det3(u, v2, v3) ≤ 0
			&& det3(v1, u, v3) ≤ 0 && det3(v1, v2, u) ≤ 0)
			return -1
		end
	end
	return 0
end


struct AffineRayInv{D,T} <: Geometry{D,T}
	origin::Vec{D,T}
	direction::Vec{D,T}
	inv_dir::Vec{D,T}
end
_inv0(x::Real) = iszero(x) ? zero(real_type(x)) : inv(to_real(x))
@inline inv(a::AffineRay{D,T}) where{D,T} =
	AffineRayInv{D,real_type(T)}(a.origin, a.direction, _inv0.(a.direction))

"""
    intersects(a::AffineRay{3}, b::BBox)

https://tavianator.com/2011/ray_box.html
"""
function intersects(a::AffineRayInv{3}, b::BBox; strict=false)
	tmin = zero(eltype(a.direction))
	tmax = typemax(eltype(a.direction))
	for i in 1:3
		# position = x + t*z, hence t = (position-x)*z1
		x = a.origin[i]; z = a.direction[i]; z1 = a.inv_dir[i]
		u = b.proj[i]
		if iszero(z)
			if x < u.low || x > u.high return false; end
		else
			# if z > 0 then p_min = x+tmin*z, p_max = x+tmax*z
			#   (thus tmin = (pmin-x)*z1, tmax = (pmax-x)*z1)
			# else tmin, tmax are swapped
			@inline f(v) = (v-x)*z1
			if z > 0
				tmax = min(tmax, f(u.high))
				tmin = min(tmin, f(u.low))
			else
				tmax = min(tmax, f(u.low))
				tmin = min(tmin, f(u.high))
			end
			if strict
				if tmin ≥ tmax return false; end
			else
				if tmin > tmax return false; end # early abort
			end
		end
	end
	return true
end

function intersections(a::AffineRay, s::AbstractSurface,
		minface = 0)
	return sum(intersects(a, Triangle(s,i), i>minface )
		for i in eachindex(faces(s)))
end
# Operations on triangulated surfaces««1
# Triangle iterator««2
struct TrianglesIterator{T,S} <:
		AbstractVector{Triangle{3,T,SVector{3,Point{3,T}}}}
	surface::S
end
function Base.getindex(tri::TrianglesIterator, i::Integer)
	s = tri.surface
	return Triangle(vertices(s)[faces(s)[i]])
end
"""
    triangles(s::AbstractSurface)

An iterator over all the `Triangle` geometric faces of `s`.
"""
triangles(s::AbstractSurface) =
	TrianglesIterator{coordtype(s),typeof(s)}(s)
Base.size(tri::TrianglesIterator) = (nfaces(tri.surface),)

# Face edges iterator««2
struct FaceEdgesIterator{T}
	vertices::T
end
function iterate(itr::FaceEdgesIterator, s::Int = 1)
	s > length(itr.vertices) && return nothing
	(e1, e2) = itr.vertices[[s, mod1(s+1, length(itr.vertices))]]
	return (e1 < e2 ? SA[e1, e2] : SA[e2, e1], s+1)
end
@inline length(itr::FaceEdgesIterator) = length(itr.vertices)
"""
    face_edges(f)

Returns a list of edges bordering this face, in standard form.
"""
@inline face_edges(f) = FaceEdgesIterator(f)

# Detecting opposite faces««2
# Returns true iff these two faces (vectors of three vertex indices) are
# opposite.
function opposite_faces(f::AbstractVector{<:Integer},
	g::AbstractVector{<:Integer})
	return ((g[1] == f[1] && g[2] == f[3] && g[3] == f[2])
		  ||(g[2] == f[2] && g[1] == f[3] && g[3] == f[1])
			||(g[3] == f[3] && g[1] == f[2] && g[2] == f[1]))
end

function remove_opposite_faces(flist)
	keep = trues(length(flist))
	for (i, f) in pairs(flist), j in 1:i-1
		g = flist[j]
		opposite_faces(f, g) && (keep[i] = keep[j] = false)
	end
	return flist[keep]
end
# Adjacency««2
"""
    adjacency_points(s::AbstractSurface)

Returns the adjacency matrix on points of s, indexed by entries of
`vertices(s)`.
"""
@inline adjacency_points(s::AbstractSurface) =
	adjacency_points(vertices(s), faces(s))
function adjacency_points(points, faces)
	n = length(points)
	m = spzeros(Bool,n, n)
	for f in faces
		for i in 1:length(f), j in 1:i-1
			m[f[i],f[j]] = m[f[j],f[i]] = true
		end
	end
	return m
end

"""
    edge_can

Puts an edge in canonical form. Returns `(sign, canonical)`,
where the boolean is `true` if the edge was reversed.
"""
@inline edge_can(e) =
	(e[1] < e[2]) ? (false, SA[e[1],e[2]]) : (true, SA[e[2],e[1]])
# this must return an array because we use array[array] shorthand

# Merging, simplifying and selecting««2
"""
    same_points(points; ε)

Returns all pairs of coincident points (up to ε in ‖⋅‖∞) in this set.
"""
function same_points(points::AbstractVector{<:Point{3}}; ε=_THICKNESS)
	boxes = [ BBox(p, p+SA[ε,ε,ε]) for p in points ]
	return SpatialSorting.intersections(boxes)
end
"""
    merge(s1::AbstractSurface, s2::AbstractSurface)

Combines both triangulations, renumbering points of `s2` as needed.
(Numbering in `s1` is preserved).
"""
@inline function merge(slist::AbstractSurface...; ε=_THICKNESS)
	return simplify(concatenate(slist...); ε)
end
function concatenate(slist::AbstractSurface...)
# 	renum = [ collect(offset[i]+1:offset[i]+nvertices(s))
# 		for (i,s) in pairs(slist) ]
	newfaces = sizehint!(copy(faces(first(slist))), sum(nfaces.(slist)))
	newpoints= sizehint!(copy(vertices(first(slist))), sum(nvertices.(slist)))
	offset = 0
	for i in 2:length(slist)
		offset+= nvertices(slist[i-1])
		push!(newfaces, [ f .+ offset for f in faces(slist[i]) ]...)
		push!(newpoints, vertices(slist[i])...)
	end
	return Surface(newpoints, newfaces)
end

"""
    simplify(s::AbstractSurface)

Removes duplicate points in `s`, renumbering as needed.
"""
function simplify(s::AbstractSurface; ε=_THICKNESS)
	(newpoints, reindex) = simplify_points(vertices(s); ε)
	return Surface(newpoints, [ reindex[f] for f in faces(s) ])
end

"""
    simplify_points(points; ε)

Removes duplicates from the set of points, returning (list of new points,
map from old index to new index).
"""
function simplify_points(points; ε=_THICKNESS)
	n = length(points)
	samepoints= extrema.(same_points(points; ε))
	merged = collect(1:n)
	# merged[j]: oldindex of point replacing this one
	for (i, j) in samepoints
		(merged[j] > merged[i]) && (merged[j] = merged[i])
	end
	newindex = similar(merged)
	newpoints = sizehint!(similar(points, 0), n)
	for (i, j) in pairs(merged)
		if i == j # this is a kept point
			push!(newpoints, points[i])
			newindex[i] = length(newpoints)
		else # this is a relabeled point
			newindex[i] = newindex[j]
		end
	end
	return (newpoints, newindex)
end


"""
    select_faces(f, s::AbstractSurface)
    select_faces(list, s::AbstractSurface)

Returns the subcomplex containing only the faces `i` for which `f(i)`
evaluates to a true value. Points are renamed.
"""
function select_faces(list::AbstractVector{<:Integer}, s::AbstractSurface)
	@debug "select_faces: keeping $(length(list)) out of $(nfaces(s)) faces««"
	keep = trues(length(list))
	for (i, fi) in pairs(list), j in 1:i-1
		fj = list[j]
# 		@debug "examine: ($fi, $fj) => $(faces(s)[fi]), $(faces(s)[fj])"
		if opposite_faces(faces(s)[fi], faces(s)[fj])
			@debug "  faces $fi=$(faces(s)[fi]) and $fj=$(faces(s)[fj]) are opposite"
			(keep[i] = keep[j] = false)
		end
	end
	list = list[keep]
	@debug "after removing opposites, $(length(list)) remain"
	renum = fill(0, eachindex(vertices(s)))
	newfaces = SVector{3,Int}[]
	newpoints = similar(vertices(s),0)
	for i in list
		if i > 0
			f = faces(s)[i]
		else
			f = reverse(faces(s)[-i])
		end
		for p in f
			if iszero(renum[p])
				push!(newpoints, vertices(s)[p])
				renum[p] = length(newpoints)
				@debug "point $p = $(vertices(s)[p]) renamed $(renum[p])"
			end
		end
		@debug "face $f renamed $(renum[f])"
		push!(newfaces, renum[f])
	end
	@debug "end select_faces»»"
	return (typeof(s))(
		newpoints,
		remove_opposite_faces(newfaces))
end
@inline select_faces(test::Function, s::AbstractSurface) =
	select_faces(filter(test,eachindex(faces(s))), s)


# Surfaces with incidence information««1
# Basic type««2
"""
    AbstractSurfaceIncidence

A triangulated surface with incidence data.

    `inc_pf(s)`  vertex -> faces
    `inc_ef(s)`  edge -> (oriented) faces
		`inc_ff(s)`  face -> (non-oriented) faces
"""
abstract type AbstractSurfaceIncidence{T} <: AbstractSurface{T} end
@inline inc_pf(s::AbstractSurfaceIncidence) = _inc_pf(incidence(s))
@inline inc_ef(s::AbstractSurfaceIncidence) = _inc_ef(incidence(s))
@inline inc_ff(s::AbstractSurfaceIncidence) = _inc_ff(incidence(s))

struct SurfaceIncidence{T, S<:AbstractSurface{T}} <: AbstractSurfaceIncidence{T}
	surface::S
	inc_pf:: Vector{Vector{Int}}
	inc_ef::DictOfLists{SVector{2,Int}, Int}
	inc_ff::Vector{Vector{Int}}
end

@inline surface(s::SurfaceIncidence) = s.surface
@inline incidence(s::SurfaceIncidence) = s
@inline _inc_pf(s::SurfaceIncidence) = s.inc_pf
@inline _inc_ef(s::SurfaceIncidence) = s.inc_ef
@inline _inc_ff(s::SurfaceIncidence) = s.inc_ff

# Constructor from generic triangulated surface««2
"""
    SurfaceIncidence(s::AbstractSurface)

Returns an incidence and ajacency structure for the simplicial complex s.
This returns a named tuple with fields:
 - `points`: adjacency for points;
 - `faces`: adjacency for faces;
 - `edge_faces`: incidence edge -> face;
 - `point_faces`: incidence point -> face;
"""
function SurfaceIncidence(s::AbstractSurface;
		vf=true, ef=true, ff=true)
  inc_pf = [Int[] for p in vertices(s)]
	# face adjacency needs edge-face:
	ef = ef||ff
	if vf
		for (i, f) in pairs(faces(s)), p in f
			push!(inc_pf[p], i)
		end
	end
	inc_ef = DictOfLists{SVector{2,Int}, Int}()
	if ef
		for (i, f) in pairs(faces(s)), u in 1:3
			e = (f[u], f[plus1mod3[u]])
			(b, c) = edge_can(e)
			listpush!(inc_ef, c => b ? -i : i)
		end
	end
	inc_ff = [Int[] for f in faces(s)]
	if ff
		for a in values(inc_ef), i in eachindex(a), j in 1:i-1
			(f, g) = abs(a[i]), abs(a[j])
			push!(inc_ff[f], g)
			push!(inc_ff[g], f)
		end
	end

	return SurfaceIncidence{coordtype(s), typeof(s)}(s, inc_pf, inc_ef, inc_ff)
end

# """
#     connected_components(s::AbstractSurface)
# 
# Returns a vector of objects (same type as `s`), each one of which is a
# (renumbered) connected component.
# """
# @inline connected_components(s::AbstractSurface) =
# 	[ typeof(s)(p,f)
# 		for (p,f) in connected_components(vertices(s), faces(s)) ]
# function connected_components(points, faces)
# 	# Build the incidence matrix from the list of faces
# 	N = length(points)
# 	G = LightGraphs.SimpleGraph(adjacency_points(points, faces))
# 	C = LightGraphs.connected_components(G)
# 	# C is a vector of vector of indices
# 	# newindex[oldindex] = [component, new index]
# 	component = zeros(Int, N)
# 	newindex = zeros(Int, N)
# 	for (i, c) in pairs(C)
# 		for (j, p) in pairs(c)
# 			component[p] = i
# 			newindex[p] = j
# 		end
# 	end
# 	return [ (typeof(s))(
# 		# all points in component i
# 		points[filter(p->component[p] == i, 1:N)],
# 		[ [newindex[p] for p in f] for f in faces if component[f[1]] == i ]
# 		) for i in eachindex(C) ]
# end

# Neighbours««2
# returns all vertices neighbours of vertex v
function neighbors(s::AbstractSurfaceIncidence, v)
	return union([filter(≠(v), faces(s)[f]) for f in inc_pf(s)[v]]...)
end

# Manifoldness test««2
"""
    ismanifold(s::AbstractSurface)

Returns `(value, text)`, where `value` is a Bool indicating whether this
is a manifold surface, and `text` explains, if this is not manifold,
where the problem lies.
"""
function ismanifold(s::AbstractSurfaceIncidence)
	# TODO: check that triangles do not intersect
# 	inc = incidence(s) # needed here: vf, ef, ff
	for (e, f) in pairs(inc_ef(s))
		if length(f) != 2
			# edge adjacent to wrong number of faces
			return (value=false, text=(:singular_edge, e, f))
		end
		# check orientability; we stored orientation in the sign bit, this
		# simplifies the following code **a lot**:
		if f[1]*f[2] > 0
			# incompatible orientations for both faces
			return (value=false, text=(:not_orientable, e, f))
		end
	end
	for (p, flist) in pairs(inc_pf(s))
		adj = [ flist ∩ inc_ff(s)[f] for f in flist ]
		for (i, a) in pairs(adj)
			if length(a) != 2
				# face adjacent to wrong number of faces around this vertex
				return (value=false, text=(:vertex_faces_adj, p, flist[i], a))
# 					text="vertex $p: face $(flist[i]) adjacent to $(length(a)) other faces ($a)")
			end
		end
		# we need to check that all these adjacent faces form a simple loop.
		# This is easy to do by checking connectedness of the graph:
		c = falses(length(flist)); c[1] = true; n = 1
		rev = Dict{Int,Int}()
		for (u,v) in pairs(flist) rev[v] = u; end
		while n > 0
			n = 0
			for i in vcat(adj[c]...)
				if !c[rev[i]] c[rev[i]] = true; n+= 1; end
			end
		end
		if(!all(c))
			# faces around this vertex do not form a simple loop
			return (value=false, text=(:singular_vertex, p,
				[flist[i]=>adj[i] for i in eachindex(flist)]))
# 			return (value=false, text="faces around vertex $p do not forc a connected graph: $([flist[i] => adj[i] for i in eachindex(flist)])")
		end
	end
	return (value=true, text="is manifold")
end

# Surfaces with patch information««1
# # Splitting into regular components««2
# """
#     edgewise_connected_components(s)
# 
# Returns a tuple `(components, label)`.
#  - `components[c]` is the list of face indexes in `c`-th connected comp.
#  - `label[i] = c` is the component to which faced `i` belongs.
# 
# Not used. (Working on the global structure allows us to completely
# dispense from ray tracing).
# """
# function edgewise_connected_components(s::AbstractSurfaceIncidence)
# 	label = [0 for _ in eachindex(faces(s))]
# 	components = Vector{Int}[]
# 	visit = Int[]
# 	@inline function mark_face(i, n)
# 		label[i] = n; push!(components[n], i)
# 		push!(visit, i)
# 	end
# 	for (i, f) in pairs(faces(s))
# 		if !iszero(label[i]) continue; end
# 		push!(components, Int[]); n = length(components)
# 		mark_face(i, n)
# 		while !isempty(visit)
# 			i = pop!(visit)
# 			for j in inc_ff(s)[i]
# 				if !iszero(label[j]) continue; end
# 				mark_face(j, n)
# 			end
# 		end
# 	end
# 	return (components=components, label=label)
# end
# Type definition««2
"""
    AbstractSurfacePatches

Contains information describing the partition of the surface
in manifold patches.

 - `components(s)`: vector of faces in this regular component.
 - `label(s)`: label assignment (as an index in `components`) for each face.
 - `adjacency(s)`: for each pair of adjacent components, one of the adjacent edges.
"""
abstract type AbstractSurfacePatches{T} <: AbstractSurfaceIncidence{T} end
struct SurfacePatches{T,S <: AbstractSurfaceIncidence{T}} <:
		AbstractSurfacePatches{T}
	incidence::S
	label::Vector{Int} # label[face] = patch
	components::Vector{Vector{Int}} # components[patch] = [face, face…]
	adjacency::Matrix{SVector{2,Int}} # adjacency[component, component] = edge
end

@inline incidence(s::SurfacePatches) = s.incidence
@inline surface(s::AbstractSurfacePatches) = surface(incidence(s))
# this defines faces() and vertices():
@inline label(s::SurfacePatches) = s.label
@inline components(s::SurfacePatches) = s.components
@inline adjacency(s::SurfacePatches) = s.adjacency

# Constructor from surface with incidence««2
function SurfacePatches(s::AbstractSurfaceIncidence)
	label = [0 for _ in eachindex(faces(s))]
	components = Vector{Int}[]
	visit = Int[]
	adjacency = zeros(SVector{2,Int},0,0)
	@inline function mark_face(i, n)
# 		println("   (marking face $i=$(faces(s)[i]) as $n)")
		label[i] = n; push!(components[n], i)
		push!(visit, i)
	end
	for (i₀, f₀) in pairs(faces(s))
		if !iszero(label[i₀]) continue; end
		push!(components, Int[]); n = length(components)
		adjacency = let new_adjacency = similar(adjacency, n, n)
			new_adjacency[1:n-1,1:n-1] .= adjacency
			fill!(view(new_adjacency, n, :), SA[0,0])
			fill!(view(new_adjacency, 1:n-1, n), SA[0,0])
			new_adjacency
		end
		mark_face(i₀, n)
		while !isempty(visit)
			i = pop!(visit); f = faces(s)[i]
# 			println(collect(face_edges(f)))
			for e in face_edges(f)
# 				println("  adjacent edge $e")
				adj = filter(!isequal(i), abs.(inc_ef(s)[e]))
# 				println("  faces = $adj")
				if length(adj) == 1
					# regular edge: 2 adjacent faces. One is f, mark the other.
					iszero(label[adj[1]]) && mark_face(adj[1], n)
				else # singular edge
				# populate adjacency matrix
					for g in adj
						l = label[g]
						iszero(l) || (adjacency[l,n] = adjacency[n,l] = e)
					end
				end
			end
		end
	end
	return SurfacePatches{coordtype(s),typeof(s)}(s,
		label, components, adjacency)
end

# # DirectedEdgesTriMesh««1
# # Basic types««2
# struct DirectedEdgesTriMesh{T}
# 	opposite::Vector{Int} # 3×n
# 	destination::Vector{Int}
# 	from::T # Vector or Dict
# 	@inline DirectedEdgesTriMesh(;opposite, destination, from) =
# 		new{typeof(from)}(opposite, destination, from)
# end
# 
# function DirectedEdgesTriMesh(
# 		faces::AbstractVector{<:AbstractVector{<:Integer}})
# 	@assert all(length.(faces) .== 3)
# 	# vf[p] = [faces containing point p]
# 	points = union(faces...)
# 	vf = Dict(p=>Int[] for p in points)
# 	from = Dict(p=>0 for p in points)
# # 	from = Vector{Int}(undef, length(points))
# 	# face i has 3i-2..3i
# 	for (i, f) in pairs(faces), (j, p) in pairs(f[1:3])
# 		push!(vf[p], i)
# 		from[p] = 3*i-3+j
# 	end
# 	function find_edge(p, q)#««
# 	# returns index of edge pq (in this direction)
# 		for i in vf[p]
# 			f = faces[i]
# 			if f[1]==p && f[2]==q return 3*i-2
# 			elseif f[2]==p && f[3]==q return 3*i-1
# 			elseif f[3]==p && f[1]==q return 3*i
# 			end
# 		end
# 		return 0 # opposite edge not found: we are on the boundary
# 		# (e.g. when dissecting a triangle)
# 	end#»»
# 	opposite = Vector{Int}(undef, 3*length(faces))
# 	destination=Vector{Int}(undef, 3*length(faces))
# 	for (i, f) in pairs(faces)
# 		destination[3*i-2] = f[2]
# 		destination[3*i-1] = f[3]
# 		destination[3*i  ] = f[1]
# 		opposite[3*i-2] = find_edge(f[2], f[1])
# 		opposite[3*i-1] = find_edge(f[3], f[2])
# 		opposite[3*i  ] = find_edge(f[1], f[3])
# 	end
# 	return DirectedEdgesTriMesh(;
# 		opposite=opposite,
# 		destination=destination,
# 		from=from)
# end
# 
# @inline next(::DirectedEdgesTriMesh, n) = n+1-3*(n%3==0)
# @inline prev(::DirectedEdgesTriMesh, n) = n-1+3*(n%3==1)
# @inline nfaces(m::DirectedEdgesTriMesh) = fld(length(m.opposite),3)
# # @inline vertices(m::DirectedEdgesTriMesh) = m.points
# @inline opposite(m::DirectedEdgesTriMesh, ab) = m.opposite[value(ab)]
# @inline opposite!(m::DirectedEdgesTriMesh, ab, x) =
# 	m.oposite[value(ab)] = x
# @inline destination(m::DirectedEdgesTriMesh, ab) = m.destination[value(ab)]
# @inline destination!(m::DirectedEdgesTriMesh, ab, x) =
# 	m.destination[value(ab)] = x
# @inline from(m::DirectedEdgesTriMesh, pt) = m.from[pt]
# @inline from!(m::DirectedEdgesTriMesh, pt, x) = m.from[pt] = x
# 
# # @inline function new_half_edges(m:::DirectedEdgesTriMesh, n::Integer)
# # 	l = length(m.opposite)
# # 	resize!(m.opposite, l+n)
# # 	resize!(m.destination, l+n)
# # 	return (l+1:l+n)
# # end
# 
# @inline destination(m::DirectedEdgesTriMesh, ab) = m.destination[value(ab)]
# 
# struct DirectedEdgesTriFaces <: AbstractVector{SVector{3,Int}}
# 	mesh::DirectedEdgesTriMesh
# end
# @inline Base.size(itr::DirectedEdgesTriFaces) = (nfaces(itr.mesh),)
# @inline Base.getindex(itr::DirectedEdgesTriFaces, n::Integer) =
# 	SVector{3,Int}(view(itr.mesh.destination, 3*n-2:3*n))
# @inline faces(m::DirectedEdgesTriMesh) = DirectedEdgesTriFaces(m)
# 
# # Splitting edges and faces««2
# """
#     split_edge!(m::DirectedEdgesTriMesh, ab, p)
# 
# Inserts points `p` in the middle of the half-edge `ab` and its opposite.
# """
# function split_edge!(m::DirectedEdgesTriMesh, ab, pt::Point{3})
# 	n = nfaces(m)
# 	# Grab all the edge and vertex info from structure
# 	bc = next(m, ab); ca = next(m, ab)
# 	cb = opposite(m, bc); ac = opposite(m, ca)
# 	ba = opposite(m, ab); ad = next(m, ba); db = next(m, ad)
# 	da = opposite(m, ad); bd = opposite(m, db)
# 	# use the inner half-edges for computing destination:
# 	# (outer half-edges might not be defined if we have a boundary...)
# 	b = destination(m, ab); a = destination(m, ba)
# 	c = destination(m, bc); d = destination(m, ad)
# 	# Define new values for point x and triangles xbc, xad
# 	push!(m.points, pt);
# 	x = length(m.points)
# 	resize!(m.halfedges, 3*n+6)
# 	xb = 3*n+1; new_bc = 3*n+2; cx = 3*n+3
# 	xa = 3*n+4; new_ad = 3*n+5; dx = 3*n+6
# 	push!(m.from, xb)
# 	# adjust structure to record all values
# 	@inline set!(edge, dest, opp) = (
# 		m.destination[edge] = dest; m.opposite[edge] = opp; )
# 	# triangle axc
# 	ax = ab; set!(ax, x, xa)
# 	xc = bc; set!(xc, c, cx) # here we overwrite `c` by `c`...
# 	# ca unchanged
# 	# triangle bxd
# 	bx = ba; set!(bx, x, xb)
# 	xd = ad; set!(xd, d, dx) # ditto
# 	# db unchanged
# 	# triangle xbc
# 	set!(xb, b, bx); set!(new_bc, c, cb); set!(cx, x, xc)
# 	# triangle xad
# 	set!(xa, a, ax); set!(new_ad, d, da); set!(dx, x, xd)
# 	return m # or x...
# end
# Triangulations««1
# Wrapping the Triangle package««2
# hide the `Triangle` module name to prevent namespace conflict:
module LibTriangle
	import Triangle: Triangle, basic_triangulation
	using ..ConstructiveGeometry: Point, coordinates, SVector

	function edges(loop) # works with 1:n or a vector
		n = length(loop)
		return [loop[mod1(i+j-1, n)] for i in 1:n, j in 1:2]
	end
	@inline edges(loops...) = [edges.(loops)...;]

	function constrained_triangulation(vertices::Matrix{Float64},
			vmap::Vector{Int}, edge_list::Matrix{Int})
		# XXX temporary: the libtriangle call tends to segfault whenever lines
		# cross, so we show what it is called with
		s = "constrained_triangulation($(size(vertices)), $(size(vmap)), $(size(edge_list)):\n$vertices\n$vmap\n$edge_list\n"
		for (i, v) in pairs(vmap)
			s*= "\n  $(vertices[i,1]) $(vertices[i,2]) $v"
		end
# 		println(s,"\n") # debug output is not always flushed in time before segfault
		s*= "\n\n"
		for e in eachrow(edge_list)
			v1 = view(vertices, findfirst(==(e[1]), vmap), :)
			v2 = view(vertices, findfirst(==(e[2]), vmap), :)
			s*=("\n  $(v1[1]) $(v1[2]) $(v2[1]-v1[1]) $(v2[2]-v1[2]) # $e")
		end
		@debug s
		isunique(array) = length(unique(array)) == length(array)
		@assert isunique(vmap) "points must be unique: $(vmap)"
		for i in 1:size(vertices,1), j in 1:i-1
			@assert vertices[i,:] != vertices[j,:] "points must be distinct: $i, $j"
		end
		return SVector{3,Int}.(Triangle.constrained_triangulation(vertices,
			vmap, edge_list))
	end
	constrained_triangulation(vertices::AbstractVector{<:Point{2}},
			vmap, edge_list) =
		constrained_triangulation(
			Matrix{Float64}([transpose.(coordinates.(vertices))...;]),
			vmap, edge_list)
	constrained_triangulation(vertices, vmap, loops...) =
		constrained_triangulation(vertices, collect(vmap), edges(loops...))
end; import .LibTriangle

"""
    triangulate(s::PolgyonXor)

Returns a triangulation of vertices of `s`, removing holes.
"""
function triangulate(s::PolygonXor)
	v = vertices(s)
	id = identify_polygons(s)
	peri = perimeters(s)
	is_hole = falses(length(v))
	for (i, p) in pairs(peri)
		if id[i] < 0 # is hole
			is_hole[p] .= true
		end
	end
			
	tri = LibTriangle.constrained_triangulation(v, 1:length(v), peri...)
	# remove triangles made entirely of hole vertices
	return tri[[!all(is_hole[t]) for t in tri]]
end

# 2d triangulation««2
"""
    triangulate_loop(path::Path{2})

Given a closed loop of points,
returns a Delaunay triangulation of the *inside* of the loop only.
(Constrained triangulation and all triangles lying outside the loop are
removed. Triangles are oriented in the same direction as the loop.)

The triangulation is returned as a vector of [i1, i2, i3] = integer indices
in the list of points.
"""
function triangulate_loop(points::AbstractVector{<:Point{2}})
	N = length(points)
	return LibTriangle.constrained_triangulation(points, 1:N, 1:N)
end
# @inline triangulate_loop(points::AbstractVector{<:Point{2}}) =
# 	triangulate_loop(Matrix{Float64}(vcat(transpose.(coordinates.(points))...)))
# 	N = length(points)
# 	m = Matrix{Float64}(vcat(transpose.(points)...))
# 
# 	ct = LibTriangle.constrained_triangulation(
# 		m,
# 		collect(1:N), # trivial map on points
# 		[mod1(i+j-1, N) for i in 1:N, j in 1:2])
# 	# returns triangles of the same orientation as the loop:
# 	# outer triangles are removed by Triangle library
# 	# each triangle is either inner or outer
# 	# this is determined e.g. by its barycenter
# # 	return filter(t->point_in_polygon(sum(points[t])/3, points) > 0, ct)
# 	return ct
# end

# 3d face triangulation««2
"""
    triangulate_face(points; direction, map, convex)

Returns a triangulation of the face (assumed convex; points in any order)
Optional keyword arguments:
 - `direction` is a normal vector (used for projecting to 2d).
 - `map` is a labeling of points (default is identity map).
 - `convex` is a Val(true) or Val(false).

The triangulation is returned as a vector of StaticVector{3,Int},
containing the labels of three points of each triangle.
"""
function triangulate_face(
		points::AbstractVector{<:Point{3}}
		;
		direction::AbstractVector = face_normal(points),
		map::AbstractVector{<:Integer} = [1:length(points)...],
		convex::Val = Val(false)
		)
	coords = project_2d(direction)
	N = length(points)
	# this common case deserves a shortcut:
	if N == 3 return [Vec{3}(map)]; end

	points2d = Matrix{Float64}(undef, N, 2)
	for (i, p) in pairs(points), (j, x) in pairs(coords)
		points2d[i, j] = p[x]
	end
	if convex isa Val{true}
		r= Vec{3,Int}.(LibTriangle.basic_triangulation(points2d, map))
		return r
	else
		edges = vcat(([map[i] map[mod1(i+1,N)]] for i in 1:N)...)
		r= Vec{3,Int}.(LibTriangle.constrained_triangulation(points2d,
			map, edges))
		return r
	end
end
# Triangulating surfaces««2
"""
    triangulate(points, faces)

Constructs a triangulated surface from the given lists of points and
(polygonal) faces.
This triangulates all faces and removes all degenerate (collinear) faces.
"""
function triangulate(points::AbstractVector{<:Point{3}},
		faces::AbstractVector{<:AbstractVector{<:Integer}})
	triangles = Vec{3,Int}[]
	for f in faces
		# kill all 2-faces
		length(f) <= 2 && continue
		# triangulate
		thisface = triangulate_face(points[f]; map=f)
		# kill all degenerate faces
		@inline is_degenerate(t) =
			iszero(face_normal(points[t]))
		push!(triangles, filter(!is_degenerate, thisface)...)
	end
	return Surface(points, triangles)
end
triangulate(points::AbstractVector{<:AbstractVector{<:Real}},
		faces::AbstractVector{<:AbstractVector{<:Integer}}) =
	triangulate(Point{3}.(points), faces)

triangulate(points::AbstractMatrix{<:Real}, faces) =
	triangulate(ViewRows(points), faces)
# 3d union and intersection««1
# After [Zhou, Grinspun, Zorin, Jacobson](https://dl.acm.org/doi/abs/10.1145/2897824.2925901)
# Self-intersection««2
# TODO: move all intersection computations to clean little functions
"""
    self_intersect(s::AbstractSurface)

Returns all self-intersections of `s`, as a `NamedTuple`:
 - `points`: all new points of intersection (as a vector of `Point{3}`,
   to be appended to the original geometric points of the structure).
 - `edge_points`: for all edges, a vector of indices of new points
   (sorted in the direction of the edge).
 - `face_points`: for all faces of `s`, the list of new points in this face,

Point indices are returned as indices in `vertices(s)` ∪ {new points}.

"""
function self_intersect(s::AbstractSurfaceIncidence)
# 	println("incidence...")
# 	inc = incidence(s; vf=false) # we only need edge_faces
# 	println("planes...")
	@debug "self-intersect ($(nvertices(s)) vertices, $(nfaces(s)) faces)««"
	@debug " Input surface:\n"*strscad(s)
	# we precompute all planes: we need them later for edge intersections
	planes = [ supporting_plane(tri) for tri in triangles(s) ]

	n = nvertices(s)
# 	println("self_intersect: $n points at beginning")
	new_points = similar(vertices(s), 0)
	T = eltype(eltype(vertices(s)))
	face_points = DictOfLists{Int,Int}()
	edge_points = Dict([ k=>Int[] for k in keys(inc_ef(s)) ])
	edge_coords = Dict([ k=>T[] for k in keys(inc_ef(s)) ])# used for sorting
	@inline function create_point!(p)
		@debug "inserting point $p"
		j = findfirst(isapprox(p;atol=_THICKNESS), new_points)
		if !isempty(new_points)
			m = findmin([distance²(p,q) for q in new_points])
		end
		if j isa Int; return j; end
		push!(new_points, p)
		return n + length(new_points)
	end
	@inline function add_point_edge!(e, k, p)#««
		@debug "adding point $k to edge $e ««"
		vec = vertices(s)[e[2]]-vertices(s)[e[1]]
		# fixme: unroll this loop to allow constant-propagation:
		i = argmax(abs.(vec))
		if length(edge_points[e]) == 0 # most common case
			push!(edge_points[e], k)
			push!(edge_coords[e], p[i])
			@debug "first point on this edge, insertion is trivial»"*"»"
			return
		end
# 		dprintln("  sorted by coordinate $i ($(vec[i]))")
# 		dprintln("  e=$e")
		rev = (vec[i] < 0)
		j = searchsorted(edge_coords[e], p[i]; rev=rev)
# 		dprintln("  inserting at position $j, $(first(j))")
		insert!(edge_points[e], first(j), k)
		insert!(edge_coords[e], first(j), p[i])
		@debug "  now edge_points[$e] = $(edge_points[e])\n»»"
	end#»»

# 	println("faces...")
	@debug "face-edge and face-vertex intersections: ««\n"
	# face-edge and face-vertex intersections
	for (i, f) in pairs(faces(s))
# 		println("  face $i: $f")
		# set up infrastructure for this face
		triangle = vertices(s)[f]
		bbox = boundingbox(triangle...)
		plane = planes[i]
		(proj, lift) = project_2d(plane)
		triangle2 = hrep([p[proj] for p in triangle]...) # hrep of projection

		# face-vertex intersections««
		for (j, p) in pairs(vertices(s))
			if j ∈ f || p ∉ bbox || !isapprox(plane(p), 0; atol=1e-10)
				continue
			end
			p2 = p[proj]
			if any(h(p2) <= 0 for h in triangle2)
				continue
			end
			# vertex is inside this face, mark it
			@debug "vertex $j is inside face $i=$f"
			listpush!(face_points, i => j)
		end#»»
		# face-edge intersections««
		for (e, flist) in pairs(inc_ef(s))
			segment = Segment(vertices(s)[e]...)
			if isempty(boundingbox(segment) ∩ bbox) || !isempty(e ∩ f)
				continue
			end
			# FIXME move this to a segment ∩ triangle function
			(z1, z2) = plane.(vertices(segment))
			# correct for some rounding errors for vertices in the plane
			(z1*z2 ≥ 0 || abs(z1) <= 1e-10 || abs(z2) <= 1e-10 ) && continue

			# the previous line ensures that z1 ≠ z2, so this never fails:
			(a1, a2) = coordinates.(vertices(segment))
			p2 = Point((z2*a1[proj] - z1*a2[proj])/(z2-z1))
			if any(h(p2) <= 0 for h in triangle2)
				continue
			end
			p3 = lift(p2)
			k = create_point!(p3)
			@debug "+p$k = $p3: edge $e ∩ face $i=$f"
			@debug "adding point $k to face $i=$f"
			@assert norm(cross(vertices(s)[e[1]]-vertices(s)[e[2]], p3-vertices(s)[e[2]])) <= 1e-3*norm(vertices(s)[e[1]]-vertices(s)[e[2]])
# 			dprintln("cross: $(cross(vertices(s)[e[1]]-vertices(s)[e[2]], p3-vertices(s)[e[2]]))")
			listpush!(face_points, i=>k)
			add_point_edge!(e, k, p3)
		end#»»
	end
	@debug "»»\n"
# 	println("edges...")
	@debug "edge-edge and edge-vertex intersections:««\n"
	# edge-edge and edge-vertex intersections
	for (e, flist) in pairs(inc_ef(s))
		# two equations define this edge:
		# first one is that of an adjacent face
		eq1 = planes[abs(flist[1])]
		# for the second equation, it happens (quite often) that edges
		# delimitate two parallel faces, so we cannot use a second face
		# equation. Instead we project on first plane and look for an
		# equation of the projection here.
		v = direction(eq1)
		(proj, kmax) = project_2d(v, Val(true))
		eq2 = line(vertices(s)[e[1]][proj] => vertices(s)[e[2]][proj])

		bbox = boundingbox(vertices(s)[e]...)
		# edge-vertex intersections:««
		for (j, p) in pairs(vertices(s))
			if p ∉ bbox || j ∈ e || eq2(p[proj]) ≠ 0 || eq1(p) ≠ 0
				continue
			end

			# edge (segment) is intersection of bbox and line,
			# therefore here we know that the point is on the edge:
			@debug "vertex $j is on edge $e"
			add_point_edge!(e, j, p)
		end#»»
		# edge-edge intersections:««
		for (e1, flist1) in pairs(inc_ef(s))

			# TODO: could this be made simpler by just checking if determinant
			# is zero?
			segment = Segment(vertices(s)[e]...)
			# this makes the iteration triangular:
			if e1 == e break; end
			seg1 = Segment(vertices(s)[e1]...)
			if isempty(boundingbox(seg1) ∩ bbox) continue; end
			if !isempty(e ∩ e1) continue; end
			p = inter(segment, seg1; thickness=_THICKNESS)
			if p isa Nothing continue; end
			# check that this is not a vertex of one of the edges
			(a1, b1, a2, b2) = vertices(s)[[e[1], e[2], e1[1], e1[2]]]
			isapprox(p, a1; atol = _THICKNESS) && continue
			isapprox(p, b1; atol = _THICKNESS) && continue
			isapprox(p, a2; atol = _THICKNESS) && continue
			isapprox(p, b2; atol = _THICKNESS) && continue
# 			@debug "eq2 => $(eq2(p[proj]))"
# 			if eq2(p[proj]) ≠ 0 continue; end
			# point p is a new point and on both edges e and e1
			@debug "edges $e and $e1 intersect"
			k = create_point!(p)
			@debug "+p$k = $p: edges $e ∩ $e1"
			add_point_edge!(e, k, p)
			add_point_edge!(e1, k, p)
		end#»»
	end
	@debug " end of edges»»\n"
	@debug face_points
# 	@debug join("computed intersections: ",
# 		["\n face $f=$(faces(s)[f]): $p" for(f,p) in pairs(face_points)])
	@debug join(["\n edge $e: $p" for(e,p) in pairs(edge_points)])
	@debug " end of self-intersect»»\n"
	return (points = new_points,
		edge_points = edge_points,
		face_points = face_points)
end
function self_int2(s::AbstractSurfaceIncidence; ε=_THICKNESS)
	boxes = [ boundingbox(t...) for t in triangles(s) ]
	for (i1, i2) in intersections(boxes)
		tri1 = triangle(s, i1)
		tri2 = triangle(s, i2)
	# we know that the bounding boxes of faces (i) and (j) intersect,
	# now determine the intersection type of those two triangles
	# possible types:
	# vv: ignore this (we will simplify points in the next step)
	# ve:
	# vf: 6 cases (vertices of i in face j, and conversely)
	# ee: 9 possibilities (3 × 3 edges)
	# ef: 6 cases (edges of i in face j)
	# ff: coincident supporting planes; counted as “vf” or “ve”
	end
end

# Sub-triangulation««2
# FIXME: thickness -> ε
function coplanar_faces(s::AbstractSurfaceIncidence,
		hyperplane::Polyhedra.HyperPlane,
		thickness = 0)
	# looking for connected-component faces is not enough (e.g.
	# intersection of two objects havin a common flat side).
	coplanar = Set{Int}()
	for (j, f) in pairs(faces(s))
		v = hyperplane.(vertices(s)[f])
		if all(<(thickness), abs.(v))
			n = direction(supporting_plane(Triangle(vertices(s)[f])))
			push!(coplanar, j*sign(dot(n, direction(hyperplane))))
		end
	end
	return coplanar
end

# FIXME: after [ZGZJ], this should be done in *clusters* of coplanar
# faces, so as to ensure compatible triangulation in exceptional cases.

"""
    subtriangulate(s::AbstractSurface)

Returns a refined triangulation of `s` with vertices at all
self-intersection points.
"""
function subtriangulate(s::AbstractSurfaceIncidence)
# 	println("self-intersect...")
	self_int = self_intersect(s)
# 	println("subtriangulate...")
# 	explain(s, "/tmp/before-subtriangulate.scad", scale=30)
	@debug "subtriangulate ($(nvertices(s)) vertices, $(nfaces(s)) faces)««"
	newpoints = [ vertices(s); self_int.points ]
	newfaces = SVector{3,Int}[]
	@inline edge_points(e1, e2) =
		e1 < e2 ? self_int.edge_points[SA[e1,e2]] :
		reverse(self_int.edge_points[SA[e2,e1]])
	
	# all faces to triangulate
	# they will be removed from this set as triangulation is complete
	faces_todo = Set(keys(self_int.face_points))
	for e in keys(self_int.edge_points)
		push!(faces_todo, abs.(inc_ef(s)[e])...)
	end
	@debug "$(length(faces_todo)) faces to do: $faces_todo" *
		join(["\n $i=$(faces(s)[i])" for i in faces_todo])
	while !isempty(faces_todo)
		i = first(faces_todo)
		plane = normalize(supporting_plane(triangles(s)[i]))
		@debug "from face $i=$(faces(s)[i]): $plane"
		cluster = coplanar_faces(s, plane, _THICKNESS)
		push!(cluster, i)
		str="cluster=$cluster"*
			join(["\n $i=$(faces(s)[abs(i)])" for i in cluster])
		@debug "««cluster=$cluster"*
			join(["\n $i=$(faces(s)[abs(i)])" for i in cluster])

		proj = project_2d(direction(plane))
		pset = Set{Int}()
		eset = Set{NTuple{2,Int}}()
		fp = Dict{Int,Set{Int}}()
		for j in cluster # i0 is a signed face
			i = abs(j)
			f = faces(s)[i]
			fp[j] = Set{Int}(f)
# 			push!(pset, f...)
			haskey(self_int.face_points, i) &&
				push!(fp[j], self_int.face_points[i]...)
			for k in 1:3
				k1 = plus1mod3[k]
				(_, e) = edge_can([f[k], f[k1]])
				line = [e[1]; self_int.edge_points[e]; e[2]]
				push!(fp[j], line...)
# 				@debug "face $i=$f: pushing edge $e => $line"
				for i in 1:length(line)-1
					@inbounds push!(eset, minmax(line[i], line[i+1]))
				end
			end
			@debug "face $i=$f contains points $(fp[j])"
			union!(pset, fp[j]...)
		end
# 		@debug "using points $pset and edges $eset"
		pvect = collect(pset)
		pcoord = [ newpoints[p][proj] for p in pvect ]
		emat = [ e[i] for e in eset, i in 1:2 ]
# 		@debug "triangulation: $pcoord, $pvect, $emat"
		tri = LibTriangle.constrained_triangulation(pcoord, pvect, emat)
		@debug "returned $tri"
		# triangles from each face are added separately;
		# this preserves multiplicity
		for (i, ps) in pairs(fp) # ps is a set of points
			@debug "adding triangles from face $i=$(faces(s)[abs(i)]): $ps««"
			for t in tri
# 				@debug "  examining $(Vector(t))"
				# FIXME: check orientation of face!
				if issubset(t, ps)
					t1 = i > 0 ? t : SA[t[1],t[3],t[2]]
					@debug "  (face $i) pushing $t1 to newtriangles"
					push!(newfaces, t1)
				end
			end
			@debug "»»"
		end
		@debug "cluster done»»"

		for j in cluster
			delete!(faces_todo, abs(j))
		end
	end
	@debug "newfaces=$(Vector.(newfaces))"

# 	for (i, f) in pairs(faces(s))#««
# 		extra = get(self_int.face_points, i, Int[])
# # 		println("face $i=$f: got extra points $extra")
# 		perimeter =
# 			[ f[1]; edge_points(f[1], f[2]);
# 			  f[2]; edge_points(f[2], f[3]);
# 				f[3]; edge_points(f[3], f[1]); ]
# 		# common case: nothing was added; in this case, skip this face:
# 		if length(extra) == 0 && length(perimeter) == 3
# # 			println("nothing to do, keeping face $f")
# 			push!(newfaces, f)
# 			continue
# 		end
# 		@debug "subtri f=$f: $extra, perim=$(edge_points(f[1], f[2])) + $(edge_points(f[2],f[3])) + $(edge_points(f[3],f[1]))"
# 		triangle = Triangle(vertices(s)[f]...)
# 		plane = supporting_plane(triangle)
# 		proj = project_2d(direction(plane))
# 		@debug "triangle = $triangle, plane = $plane, proj = $proj"
# 
# 		plist = [ perimeter; extra] # indices of points in face
# 		# as a matrix for `constrained_triangulation`:
# 		coords = [ newpoints[p][i] for p in plist, i in proj ]
# # 		println("perimeter = $perimeter")
# 		l = length(perimeter)
# 		cons = [perimeter[mod1(i+j,l)] for i in eachindex(perimeter), j in 0:1]
# # 		for (i, p) in pairs(plist)
# # 			println("$(coords[i,1]) $(coords[i,2]) $p")
# # 		end
# # 		println("($coords, $plist, $cons)")
# 		tri = LibTriangle.constrained_triangulation(coords, plist,
# 			LibTriangle.edges(perimeter))
# # 		println("triangulation = $tri")
# 		push!(newfaces, tri...)
# 		@debug "returned triangulation=$(Vector.(tri))"
# 	end#»»

	@debug "(end subtriangulate)»»"
# 	newfaces = remove_opposite_faces(newfaces)
	return Surface(newpoints, newfaces)
end

# Faces around an edge««2
"""
    faces_around_edge(s, edge, incidence, [vector = 0])

Returns a cyclically ordered list of all faces of `s` around edge `e`,
with sign indicating the orientation of the face. (The list starts at an arbitrary index).

If a `vector` is provided then this will return a (signed)
face matching this vector.
"""
function faces_around_edge(s::AbstractSurfaceIncidence,
		edge, vec3 = zero(Vec{3,coordtype(s)}))
	# we project the faces on the plane perpendicular to edge e;
	# the eye is at position e[2] looking towards e[1].
	dir3 = vertices(s)[edge[2]]-vertices(s)[edge[1]]
	(proj, k) = project_2d(dir3, Val(true))
	dir2 = dir3[proj]
	dir2scaled = dir2/norm²(dir3)
	flist = inc_ef(s)[edge]
	# for each adjacent face, compute a (3d) vector which, together with
	# the edge, generates the face (and pointing from the edge to the face):
	# 2d projection of face_vec3 (preserving orientation)
	face_vec2 = begin
		face_pt3 = [sum(faces(s)[abs(f)]) - sum(edge) for f in flist]
		face_vec3 = [ vertices(s)[p] - vertices(s)[edge[2]] for p in face_pt3 ]
		[ v[proj] - (v ⋅ dir3)*dir2scaled for v in face_vec3 ]
	end
	reorder = sort(eachindex(flist);
		lt=(i, j) -> let b = circular_sign(face_vec2[i], face_vec2[j])
			if !iszero(b) return (b > 0)
			# the use of **signed** face numbers guarantees consistent ordering
			# even for coincident faces
			# e.g. let 0 < f < g be two coincident faces
			# on an edge for which both are positive:
			# f < g, cell is f­(c)->g
			# on an edge for which both are negative:
			# -g < -f, cell is g<-(c)-f
			#
			# for opposite faces: 0 < f < g, opposite
			# on an edge with f positive, order is -g < 0 < f
			# cell is g<-(c)->f
			# on an edge with g positive, order is -f < 0 < g
			# cell is -f<-(c)->f
			else return flist[i] > flist[j]
			end end)
	if !iszero(vec3)
		vec2 = vec3[proj] - (vec3⋅dir3)*dir2scaled
		@assert !iszero(vec2) "edge $edge aligned with point $vec3"
		@debug "searching vec2=$vec2 in list: $(face_vec2[reorder])"
		k = searchsorted(face_vec2[reorder], vec2,
			lt = (u, v)->circular_sign(u, v) > 0)
		@assert k.start > k.stop "impossible to determine point location at this edge"
		if k.stop == 0 # k==(1:0): we insert before first edge
			return -flist[reorder[1]]
		else # k==(i+1:i): 
			return flist[reorder[k.stop]]
		end
	end
	str = "faces around edge $edge: $flist\n"
	for i in reorder
		f = abs(flist[i])
		p3 = sum(faces(s)[f]) - sum(edge)
		str*= "  face $(flist[i]) to vertex $p3: vec2=$(face_vec2[i])\n"
	end
	@debug str
	return flist[reorder]
end

# Level structure««2
"""
    LevelStructure

Records the graph of levels.

level[i] = (reference point, delta)
component[reference point] = [all i starting from this point]
"""
struct LevelStructure
	level::Vector{Tuple{Int,Int}}
	component::Vector{Set{Int}}
	@inline LevelStructure(n) = new(
		[(i, 0) for i in 1:n],
		[Set(i) for i in 1:n])
end
function connect!(l::LevelStructure, i1, i2, delta)
	@debug "connect($i1, $i2, $delta)"
	(c1, r1) = l.level[i1]
	(c2, r2) = l.level[i2]
	if c1 == c2
		@assert r1 + delta == r2 "connect($l, $i1, $i2, $delta): faces are already connected and multiplicity should be $(r1+delta) (it is $r2). Faces are probably crossed around one edge."
		return
	end
	# connect everything in component c = l.level[i2][1]
	for k in l.component[c2]
		# r is replaced by l.level[i1]+delta
		# r' is replaced by (r'-r) + (l.level[i1]+delta)
		l.level[k] = (c1, l.level[k][2] - r2 + r1 + delta)
	end
	union!(l.component[c1], l.component[c2])
	empty!(l.component[c2])
	return l
end
@inline connected(l::LevelStructure, i1, i2) =
	l.level[i1][1] == l.level[i2][1]

# Returns a list of connected components for this level structure.
# Each element of the list represents one connected component,
# arranged by multiplicity:
# [ [faces-with-multiplicity = 1], [faces-with-mult = 2], ... ]
function components(l::LevelStructure)
	ret = Vector{Vector{Int}}[]
	for c in l.component
		isempty(c) && continue
		(mmin, mmax) = extrema(i->l.level[i][2], c)
		r = [ Int[] for m in mmin:mmax ]
		for i in c
			m = l.level[i][2] # multiplicity of patch c
			push!(r[1+m-mmin], i)
		end
		push!(ret, r)
	end
	return ret
end

# Point location algorithm««2
# finds a good edge from point i, viewed from point k
function find_good_edge(s::AbstractSurfaceIncidence, i, vp)
	@debug "finding good edge from $i relative to $vp"
	l = neighbors(s, i)
	# it is possible that all edges lie in the same plane
	# (if this is a flat cell), so we pick, as a second point j,
	# the one which maximizes |y/x|, where
	# y = ‖pi∧pj‖, x=‖pi·ij‖/‖pi‖²
	# in other words, we maximize ‖pi∧pj‖²/‖pi·ij‖²
	# caution: it is possible that pi⋅ij = 0 (edge exactly orthogonal),
	# in this case we must return j immediately
# 	vp = vertices(s)[p];
	vi = vertices(s)[i]; vpi = vi - vp
	j = l[1]; vj = vertices(s)[j]; vij = vj - vi
	sp = vpi⋅vij
	iszero(sp) && return j
	xyj = (sp*sp, norm²(cross(vpi, vij)))
	@debug "initial vertex $vp to $i=$vi"
	@debug "    initial edge to $j = $vj: $xyj"
	for k in l[2:end]
		vk = vertices(s)[k]; vik = vj - vk
		sp = vpi⋅vik
		iszero(sp) && return k
		xyk = (sp*sp, norm²(cross(vpi, vik)))
		@debug "   trying edge $vp to $k=$vk: xy = $xyk"
		if xyk[2]*xyj[1] > xyk[1]*xyj[2]
			j = k; xyj = xyk
		end
	end
	return j
end
"""
    locate_point(s, c, point)

Determines if the point is inside (returns +1) or outside (returns -1) of connected component `c`.
"""
function locate_point(s::AbstractSurfaceIncidence, vlist, point)
	@debug "locate_point($vlist, $point)««"
	i = vlist[argmin([distance²(vertices(s)[t], point) for t in vlist])]
	@debug "  closest vertex is $i = $(vertices(s)[i])"
	j = find_good_edge(s, i, point)
	@debug "  best neighbor is $j = $(vertices(s)[j])"
	(_, edge) = edge_can([i, j])
	@debug "  best edge is $edge"
	f = faces_around_edge(s, edge, point - vertices(s)[i])
	@debug "  faces_around_edge returned $f, returning its sign"
	@debug "(end locate_point)»»"
	return sign(f)
end

# Computation of multiplicity levels««2
# `c` is a vector of regular component numbers;
# returns the list of all vertices belonging to one of the regular
# components in `c`
function vertices_in_components(s::AbstractSurfacePatches, c)
	flist = union(components(s)[[c...]]...)
	return union(faces(s)[flist]...)
end

"""
    multiplicity_levels(s)

Given a triangulated surface `s`,
returns the set of level surfaces enclosing each multiplicity component
of `s`.
"""
function multiplicity_levels(s::AbstractSurfacePatches)
	@debug "multiplicity_levels ($(nvertices(s)) vertices, $(nfaces(s)) faces)««"
# 	@debug " Input surface:\n"*strscad(s)
# 	explain(s, "/tmp/before-multiplicity.scad", scale=40)
# 	println("multiplicity...")
	@debug "regular components: $(components(s))"
	@debug "  labeling: $(label(s))"
	@debug "  adjacency: $(adjacency(s))"

	levels = LevelStructure(length(components(s)))
	for (i1, r1) in pairs(components(s)), i2 in 1:i1-1
		connected(levels, i1, i2) && continue
		edge = adjacency(s)[i1, i2]
		iszero(edge) && continue
		r2 = components(s)[i2]
		@debug "regular components $i1 and $i2 meet at edge $edge««"
		flist = faces_around_edge(s, edge)
		str="ordered faces at this edge (viewed from $(edge[2]) to $(edge[1])):\n"
		for f in flist
			f1 = abs(f)
			v = sum(faces(s)[f1]) - sum(edge)
			str*=string("  ", (f > 0) ? "↺" : "↻",
				"(f $f, c $(label(s)[f1]), v $v)\n")
		end
		@debug str
		for (j, f) in pairs(flist)
			# connect consecutive patches around an edge,
			# depending on their orientations
			j1 = mod1(j+1, length(flist)); f1 = flist[j1]
			dir = sign(f1) + sign(f)
			k = (sign(f) + sign(f1)) >> 1
			c = label(s)[abs(f)]; c1 = label(s)[abs(f1)]
			connect!(levels, c, c1, -k)
		end
		@debug "(end edge)»»"
	end
	@debug levels
	# connected components of level structure;
	# each component is arranged by multiplicity:
	cc = components(levels)
	@debug cc
	# now we rearrange these connected components by inclusion
	# vlist[i] = list of all vertices in connected component i:
	# vmax[i] = index of an extremal vertex in cc i
	# pmax[i] = Point corresponding to this vertex
	# nest[i] = number of components enclosing this one
	vlist = [ vertices_in_components(s, union(c...)) for c in cc ]
	vmax = [ v[argmax([coordinates(vertices(s)[i])[1] for i in v])]
		for v in vlist]
	pmax = vertices(s)[vmax]
	nest = fill(0, length(cc))
	
	@debug join(["cc $i: vmax=$(vmax[i]): $(pmax[i])" for i in 1:length(cc)], "\n")

	# iterate over distinct pairs of components to determine inclusion
	for i1 in 1:length(cc), i2 in 1:length(cc)
		i1 == i2 && continue
		k = locate_point(s, vlist[i2], pmax[i1])
		# if k > 0 then component i1 is inside i2
		nest[i1] += (k < 0)
	end
	@debug "cc = $cc"
	@debug "nest = $nest"
	# mpatches[i] = patches enclosing region with multiplicity i
	mpatches = [ union([get(c, i-n, Int[]) for (c, n) in zip(cc, nest)]...)
		for i in 1:maximum(nest)+maximum(length.(cc))]
	@debug "mpatches = $mpatches"
	@debug join(["region of multiplicity $i enclosed by patches $p" for (i,p) in pairs(mpatches)], "\n")

	mfaces = [ union(Int[], components(s)[p]...) for p in mpatches ]
	@debug join(["\n region $i: faces $f" for (i,f) in pairs(mfaces)],"")
	@debug "(end multiplicity_levels)»»"
	return mfaces
end
# Binary union and intersection««2
function select_multiplicity(m, s::AbstractSurface...)
	t = subtriangulate(SurfaceIncidence(merge(s...)))
	face_idx = multiplicity_levels(SurfacePatches(SurfaceIncidence(t)))
	@debug "select_multiplicity««"
	flist = get(face_idx, m, Int[])
	for f in flist
		@debug "keeping face $f=$(faces(t)[f])"
	end
	@debug "(end select_multiplicity)»»"
	return select_faces(flist, t)
end

# Extrusion ««1
# Linear extrusion««2
function mesh(s::LinearExtrude, parameters)
	g = mesh(s.child, parameters)
	@assert g isa PolygonXor
	pts2 = vertices(g)
	tri = triangulate(g)
	peri = perimeters(g)
	# perimeters are oriented ↺, holes ↻

	n = length(pts2)
	pts3 = vcat([[Point([coordinates(p); z]) for p in pts2]
		for z in [0, s.data.height]]...)
	# for a perimeter p=p1, p2, p3... outward: ↺
	# with top vertices q1, q2, q3...
	# faces = [p1,p2,q1], [p2,q2,q1], [p2,p3,q2],  [p2,q3,q2]...
	#  - bottom: identical to tri
	#  - top: reverse of tri + n
	#  - sides:
	faces = [ tri;
		[ reverse(f) .+ n for f in tri ];
		vcat([[SA[i,j,i+n] for (i,j) in consecutives(p) ] for p in peri]...);
		vcat([[SA[j,j+n,i+n] for (i,j) in consecutives(p) ] for p in peri]...);
	]
	return Surface(pts3, faces)
end
# Rotate extrusion««2
norm₁(s) = maximum(abs, [coordinates.(vertices(mesh(s)))...;])

"""
    _rotate_extrude(point, data, parameters)

Extrudes a single `Point{2}`, returning a vector of `Point{3}`.
(x,y) ↦ (x cosθ, x sinθ, y).
"""
function _rotate_extrude(p::Point{2}, data, parameters)
	@assert p[1] ≥ 0
	# special case: point is on the y-axis; returns a single point:
	p[1] == 0 && return [Point(p[1], p[1], p[2])]
	n = cld(sides(p[1], parameters) * data.angle, 360)

	T = real_type(coordtype(p))
	ω = Complex{T}(cosd(data.angle/n), sind(data.angle/n))
	z = Vector{Complex{T}}(undef, n+1)
	z[1] = one(T)
	for i in 2:n
		@inbounds z[i] = z[i-1]*ω; z[i]/= abs(z[i])
	end
	# close the loop:
	z[n+1] = Complex{T}(cosd(data.angle), sind(data.angle))
	return [Point{3,T}(p[1]*real(u), p[1]*imag(u), T(p[2])) for u in z]
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
function mesh(s::RotateExtrude, parameters)
	@debug "mesh(rotate_extrude)««"
	# take only right half of child:
	m1 = mesh(s.child); N₁ = norm₁(m1)
	g = mesh(intersect(m1,
		polygon([SA[0,-N₁],SA[N₁,-N₁],SA[N₁,N₁],SA[0,N₁]])))
	@assert g isa PolygonXor
	pts2 = vertices(g)
	tri = triangulate(g)
	peri= perimeters(g) # oriented ↺	
	n = length(pts2)
	@debug "perimeters: $peri"
	#
	pts3 = _rotate_extrude(pts2[1], s.data, parameters)
	firstindex = [1]
	arclength = [length(pts3)]
	@debug "newpoints[$(pts2[1])] = $(length(pts3))"
	for p in pts2[2:end]
		push!(firstindex, length(pts3)+1)
		newpoints = _rotate_extrude(p, s.data, parameters)
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
	@debug "end rotate_extrude»»"
	surf = Surface(pts3, triangles)
# 	explain(surf, "/tmp/mu.scad", scale=30)
	return simplify(surf)
end
# Path extrusion ««2
# triangulate_between: triangulate between two parallel paths««
"""
		triangulate_between(poly1, poly2, start1, start2)

Given two polygons `poly1` and `poly2`, both of them represented as a
vector of paths, and produced as offsets from a common path,
find a triangulation for the region between the two polygons.

This functions returns a pair `(triangulation, edge)`, where:

 - the triangulation is a vector of `SVector{3,Int}`,
where each point is represented by its index. Indices in `poly1` start at
value `start1`, and in `poly2` at `start2`.

 - the edge is a pair `(lastidx1, lastidx2)` corresponding to the last
	 points visited on each polygon. (this will be useful for closing the
	 extrusion).

"""
function triangulate_between(
		poly1::AbstractVector{<:Path},
		poly2::AbstractVector{<:Path},
		start1::Int = 1, start2::Int = 1)
	Big = typemax(coordtype(eltype(eltype(poly1))))
# 	Triangle = SVector{3,Int}
	triangles = SVector{3,Int}[]
	# head is the marker of current leading edge
	# headpoint[i] is the point marked to by head[i]
	# headidx is the new index for this marked point
	# status[i][j] is the number of last used point in j-th path of i-th poly
	head = [(1,1), (1,1)]
	headpoint = [poly1[1][1], poly2[1][1]]
	headidx = [start1, start2]
	status = zeros.(Int,length.((poly1, poly2)))
	# so far we used exactly one point on each side:
	status[1][1] = status[2][1] = 1

	# we need a way to convert (poly, path, index) to integer index««
	function first_indices(start::Int, l::Vector{Int})::Vector{Int}
		f = zeros.(Int, length(l))
		f[1] = start
		for i in 1:length(l)-1
			@inbounds f[i+1] = f[i] + l[i]
		end
		f
	end
	# firstindex[poly][path] is the first index for this path
	# firstindex[1][1] = start1
	# firstindex[1][2] = start1 + len(poly1[1]) etc.
	firstindex = (first_indices(start1, length.(poly1)),
								first_indices(start2, length.(poly2)))
	newindex(poly::Int, path::Int, index::Int)::Int =
		firstindex[poly][path] + index - 1
#»»
	# computing diagonal distances to find the smallest one:««
	distance(pt, path, i) =
		i > length(path) ? Big : distance²(pt, path[i])

	closest(pt, poly, status) =
		findmin([distance(pt, poly[i], status[i]+1) for i in eachindex(poly)])
#»»

	while true
		d1, i1 = closest(headpoint[2], poly1, status[1])
		d2, i2 = closest(headpoint[1], poly2, status[2])
		# if no more points are left, we return:
		(d1 == d2 == Big) && break

		if d1 < d2 # we append a point from poly1
			# add the triangle: head1, head2, newpoint
			s = status[1][i1] += 1
			newidx = newindex(1, i1, s)
			push!(triangles, SA[headidx[1], headidx[2], newidx])
			# update head1 to point to new point
			headidx[1] = newidx
			head[1] = (i1, s)
			headpoint[1] = poly1[i1][s]
		else
			# add the triangle: head1, head2, newpoint
			s = status[2][i2] += 1
			newidx = newindex(2, i2, s)
			push!(triangles, SA[headidx[1], headidx[2], newidx])
			# update head1 to point to new point
			headidx[2] = newidx
			head[2] = (i2, s)
			headpoint[2] = poly2[i2][s]
		end
	end
	(triangles, (headidx[1], headidx[2]))
end#»»
# path_extrude««
"""
		path_extrude(path, poly, options...)

Extrudes the given polygon (a path of points forming a simple loop)
along the given path. Both arguments are provided as a
`Vector{SVector{2}}`.

Returns a `Surface` (defined by points and a triangulation).
"""
function path_extrude(path::AbstractVector{Point{2,T}},
	poly::AbstractVector{<:Point{2}};
	join = :round,
	miter_limit::Float64 = 2.0,
	precision::Float64 = 0.25,
	closed::Bool = true
	) where{T}

	N = length(poly)
	# offset_path is a vector of vector of paths
	offset_path = offset([path], [pt[1] for pt in poly];
		join = join, ends = closed ? :fill : :butt,
		precision = precision)
	# new_points is a flat list of all 3d points produced
	new_points = [[
		[ Point([pt[1], pt[2], poly[i][2]]) for pt in [p...;] ]
		for (i, p) in pairs(offset_path)
	]...;]
# 	println("returning new_points:")

	# first index for each path
	first_face = cumsum([1; # initial
		map(p->sum(length.(p)), offset_path)])
# 	println("first_face=$first_face")

	triangles = map(1:N) do i
		i1 = mod1(i+1, N)
		triangulate_between(offset_path[i], offset_path[i1],
			first_face[i], first_face[i1])
		# XXX keep the last edge for closing the poly
	end
	# this completes the set of triangles for the tube:
	tube_triangles = vcat([ t[1] for t in triangles ]...)
	last_face = [ t[2][1] for t in triangles ]
# 	println("last_face=$last_face")
	# here we decide if it is closed or open
	# if open, triangulate the two facets
	# if closed, join them together
	if closed
		more_triangles = vcat(map(1:N) do i
			j = (i%N)+1
			[ SA[first_face[i], last_face[i], first_face[j]],
				SA[first_face[j], last_face[i], last_face[j]] ]
		end...)
# 		println("more_triangles=$more_triangles")
		tube_triangles = [ tube_triangles; more_triangles ]
	else
	# TODO: triangulate the surface
	# or, for now, close with two non-triangular facets...
		more_triangles = [ reverse(first_face), last_face ]
# 		println("more_triangles=$more_triangles")
	end
	return Surface(new_points, tube_triangles)
end#»»
# Converting 3d objects to Surfaces««1
# Primitive objects««2
mesh(s::Surface, parameters) = s

function vertices(s::Cube, parameters)
	(u,v) = (s.min, s.max)
	return Point{3}.([
		SA[u[1],u[2],u[3]],
		SA[u[1],u[2],v[3]],
		SA[u[1],v[2],u[3]],
		SA[u[1],v[2],v[3]],
		SA[v[1],u[2],u[3]],
		SA[v[1],u[2],v[3]],
		SA[v[1],v[2],u[3]],
		SA[v[1],v[2],v[3]],
	])
end
function vertices(c::Cylinder, parameters)
	p1 = unit_n_gon(c.r1, parameters)
	p2 = unit_n_gon(c.r2, parameters)
	return vcat([ c.origin + [ p; 0 ] for p in p1],
	            [ c.origin + [p; c.height ] for p in p2 ])
end
@inline vertices(s::Sphere, parameters) =
	[ s.center + p for p in fibonacci_sphere_points(s.radius, parameters) ]

# All of these are convex, so we use the lazy approach and just take
# convex hull of all the points.
function mesh(s::Cube, parameters)
	pts = vertices(s, parameters)
	return Surface(pts, [
	 SA[6, 5, 7], SA[7, 8, 6], SA[7, 3, 4], SA[4, 8, 7],
	 SA[4, 2, 6], SA[6, 8, 4], SA[5, 1, 3], SA[3, 7, 5],
	 SA[2, 1, 5], SA[5, 6, 2], SA[3, 1, 2], SA[2, 4, 3],
	])
end
function mesh(s::Union{Cylinder, Sphere}, parameters)
	p = vertices(s, parameters)
	(pts, faces) = convex_hull(p)
	return triangulate(pts, faces)
end
# CSG operations««2
function mesh(s::CSGUnion{3}, parameters)
	return select_multiplicity(1,
		[mesh(x, parameters) for x in children(s)]...)
end
function mesh(s::CSGInter{3}, parameters)
	return select_multiplicity(length(children(s)),
		[mesh(x, parameters) for x in children(s)]...)
end
function mesh(s::CSGComplement{3}, parameters)
	t = mesh(s.children[1], parameters)
	return (typeof(t))(vertices(t), reverse.(faces(t)))
end
@inline mesh(s::CSGDiff{3}, parameters) =
	mesh(intersect(s.children[1], complement(s.children[2])), parameters)
function mesh(s::CSGHull{3}, parameters)
	l = [mesh(x, parameters) for x in children(s)]
	(pts, faces) = convex_hull(vcat(vertices.(l)...))
	return Surface(pts, faces)
end
#————————————————————— Extra tools —————————————————————————————— ««1
#»»1
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
	global toplevel_objs = Geometry[]
	Base.include(x->expr_filter(obj_filter, x), Main, file)
	return union(toplevel_objs...)
end
# # TODO: somehow attach a comment indicating the origin of these objects
# # last_linenumber holds the last LineNumberNode value encountered before
# # printing this object; we use this to insert relevant comments in the
"""
    obj_filter(x)

Appends `x` to the global list of returned objects if `x` is a `Geometry`.
"""
@inline obj_filter(x) = x
@inline obj_filter(x::Geometry) =
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
# A Symbol might be an Geometry variable name, so we add it:
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
# 		find_anchor(x::Geometry, name)
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
# @inline find_anchor(x::Geometry, labels::NTuple{N,Symbol}) where{N} =
# 	find_anchor(x, sum([ labeled_anchor(x, l) for l in labels]))
# @inline function find_anchor(x::Geometry, label::Symbol)
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
export scad
# don't export include, of course
# »»1
function explain(s::AbstractSurface, io::IO = stdout; scale=1,
		offset=[0.,0.,0.], name=:m )
	println(io, """
module $name(pos=$offset, c="gray", s=$scale) {
translate($scale*pos) {
""")
	for (i, p) in pairs(ConstructiveGeometry.vertices(s))
		println(io, """
translate(s*$(Vector{Float64}(coordinates(p)))) {
	color("red") sphere(1);
	color("black", .8) linear_extrude(1) text("$i", size=5);
}
""")
	end
	println(io, "color(c, .7) polyhedron([")
	b = false
	for p in ConstructiveGeometry.vertices(s)
		print(io, b ? "," : ""); b = true
		print(io, " s*",Vector{Float64}(coordinates(p)))
	end
	println(io, "],[")
	b = false
	for f in ConstructiveGeometry.faces(s)
		print(io, b ? "," : ""); b = true
		println(io, " ", Vector{Int}(f) .- 1, " //", Vector{Int}(f))
	end
	println(io, "]); } }\n$name();")
end
@inline explain(s::AbstractSurface, f::AbstractString; kwargs...) = begin
	println("writing a surface with $(nvertices(s)) points to $f")
	open(f, "w") do io explain(s, io; kwargs...) end
end
function strscad(args...)
	b = IOBuffer()
	scad(b, args...)
	return String(take!(b))
end
end #««1 module
# »»1
# vim: fdm=marker fmr=««,»» noet:
