module Solids
# BUG: SA[1 2; 3 4]*HybridArray{Tuple{2,StaticArrays.Dynamic()}}[1 2; 3 4]
# is a (dynamic) Matrix
#

# using Printf
using LinearAlgebra, StaticArrays, FixedPointNumbers
using Polyhedra, GLPK, Triangle
# using MiniQhull
import Rotations
import Colors
import Clipper

import Base: show, print, length, getindex, size
import Base: union, intersect, setdiff, -
import Base: *, +, -

#————————————————————— Ideal objects —————————————————————————————— <<<1
#>>>1
# Types<<<1
# Numeric types <<<2

"""
    Solids._FIXED
The type used whenever a fixed-precision real number is needed (e.g.
when interfacing with `Clipper.jl`).
"""
const _FIXED = Fixed{Int64,16}
"""
    Solids._REAL

The default type used for computing coordinates (i.e. the type to which
integers are converted). Defaults to a `Fixed` type (fixed-precision
real), but this module should work just as well with floats.
"""
const _REAL = _FIXED
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

# Static vectors and paths<<<2
"""
    Vec{D,T}
The type used for representing `D`-dimensional vectors (or points) with
coordinates of type `T`. An alias for `SVector{D,T}`.
"""
const Vec{D,T} = SVector{D,T} # type alias
# this comes with the following constructors:
# Vec{2}(SA[1,2])		Vec(SA[1,2])
# Vec{2}((1,2))			Vec((1,2))
# Vec{2}(1,2)				Vec(1,2)
# Vec{2}([1,2])			- (needs explicit size) -
"""
    AnyVec{D, T}

A supertype of all plausible user input for `D`-dimensional vectors of
type `T`.
"""
const AnyVec{D,T<:Real} =
Union{StaticVector{D,<:T},AbstractVector{<:T},NTuple{D,<:T}}
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

const Path{D,T<:Real} = Vector{Vec{D,T}}
# the following functions are wrappers, which should allow us to
# rewrite `Path` objects as matrices instead (instead of `length` and
# `getindex`, return one matrix dimension and one column:)
# See also `apply` below.
@inline count(v::Path) = length(v)
@inline point(v::Path, i::Integer) = v[i]
@inline points(v::Path) = v

# This represents plausible user input for a `Path` object:
const AnyList{T} = Union{AbstractVector{<:T},NTuple{N,<:T} where{N}}
const AnyPath{D,T<:Real} = AnyList{AnyVec{D,T}}

# @inline Base.convert(::Type{Path{D,T}}, v::AnyPath) where{D,T} =
# 	[Vec{D,T}(x) for x in v]
# @inline (::Type{<:Path{D}})(v::AnyPath{D,T}, fill=zero(T)) where{D,T} =
# 	[Vec{D,T}(x) for x in v]

# Path{D,T} = HybridArray{Tuple{D,StaticArrays.Dynamic()}, T}
# @inline count(v::Path) = size(v, 2)
# @inline point(v::Path, i::Integer) = v[:,i]
# Path{D,T}(v::AnyPath{D,T}, fill = zero(T)) where{D,T} =
#		Path{D,T}([p[i] for i in 1:D, p in v])

# const Paths{T} = AbstractVector{<:Path{T}}

# Angle types<<<2
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

# Library parameters<<<2
@inline polyhedra_lib(T::Type{<:Real}) =
	Polyhedra.DefaultLibrary{T}(GLPK.Optimizer)

# AbstractSolid <<<1
# Some reasons not to use `GeometryBasics.jl`:
# - (as of 2021-01-14) not documented
# - matrix × rectangle multiplication is incorrectly defined (it returns
# a bounding box for the product)
# - inconsistency between types (e.g. Rectangle is always parallel to
# axes but Cylinder has an arbitrary orientation)
# the AbstractSolid type<<<2
"""
		AbstractSolid{D,T}

Abstract supertype for all solids.
`D` is the dimension and `T` the type used for coordinates.
A `AbstractSolid{D,T}` object may contain, as children,
objects with other coordinate types;
such coordinates will be converted to type `T` when rendering.
"""
abstract type AbstractSolid{D,T} end

"""
		dim(x::AbstractSolid)

Returns the *a priori* dimension of `x`, i.e. 2 if `x` is built purely
out of planar objects, and 3 otherwise.
"""
@inline dim(::AbstractSolid{D}) where{D} = D
@inline Base.eltype(::Type{<:AbstractSolid{D,T}}) where{D,T} = T

@inline children(::AbstractSolid) = NeutralSolid{:empty,Bool}[]
@inline parameters(::AbstractSolid) = NamedTuple()

# OpenSCAD conversion <<<2
@inline Base.show(io::IO, s::AbstractSolid) = scad(io, s, "")

"""
    scad(filename::AbstractString, s::AbstractSolid...
    scad(io::IO, s::AbstractSolid)

Prints an OpenSCAD-like representation of the given solid(s).
"""
@inline function scad(filename::AbstractString, L::AbstractSolid...)
	open(filename, "w") do f
		for s in L
			scad(f, s, "")
		end
	end
end

function scad(io::IO, s::AbstractSolid, spaces::AbstractString)
	print(io, spaces, scad_name(s))
	scad_parameters(io, s)
	L = children(s)
	if isempty(L)
		print(io, ";\n")
	else
		print(io, " {\n")
		for t in L scad(io, t, spaces*" ") end
		print(io, spaces, "}\n")
	end
end

@inline scad_parameters(io::IO, s::AbstractSolid) =
	scad_parameters(io, Val(scad_name(s)), parameters(s))

function scad_parameters(io::IO, name::Val, t::NamedTuple)
	first = true
	print(io, "(")
	for (key, val) in pairs(t)
		if first first = false; else print(io, ", "); end
		print(io, key, "=")
		# this allows a different treatment depending on (object name, key):
		scad_param_value(io, name, Val(key), val)
	end
	print(io, ")")
end
@inline scad_param_value(io::IO, n::Val, k::Val, val) = scad(io, val)
@inline scad(io::IO, val) = print(io, val)
@inline scad(io::IO, val::Fixed) = print(io, Float64(val))
function scad(io::IO, v::AbstractVector)
	first = true
	print(io, "[")
	for x in v
		if first first = false; else print(io, ", "); end
		scad(io, x)
	end
	print(io, "]")
end

# Mechanism for printing top-level objects <<<1

# FIXME: replace Main by caller module?
# FIXME: add some blob to represent function arguments
"""
		Solids.include(file::AbstractString, f::Function)

Reads given `file` and returns the union of all top-level `Solids` objects
(except the results of assignments) found in the file.

```
#### Example: contents of file `example.jl`
C=Solids.Cube(1)
S=Solids.Square(1)
Solids.Circle(3)
S

julia> Solids.include("example.jl")
union() {
 circle(radius=3.0);
 square(size=[1.0, 1.0], center=false);
}
```

"""
function include(file::AbstractString)
	global toplevel_objs = AbstractSolid[]
	Base.include(x->expr_filter(obj_filter, x), Main, file)
	return union(toplevel_objs...)
end
# # TODO: somehow attach a comment indicating the origin of these objects
# # last_linenumber holds the last LineNumberNode value encountered before
# # printing this object; we use this to insert relevant comments in the
# # output SCAD file:
# global last_linenumber = ""
# scad(io::IO = stdout) = function(x)
# 	if x isa AbstractSolid
# 		print(io, "//", last_linenumber, "\n")
# 		scad(io, x, "")
# 	end
# end
# @inline scad(x::Any) = x
# scad(x::AbstractSolid) = (scad(stdout, x, ""); x)
"""
    obj_filter(x)

Appends `x` to the global list of returned objects if `x` is a
`AbstractSolid`.
"""
@inline obj_filter(x) = x
@inline obj_filter(x::AbstractSolid) =
	(global toplevel_objs; push!(toplevel_objs, x); return x)

"""
		expr_filter(E)

Read top-level expression `E` and decides if a Solids object is returned.

The function `expr_filter` transforms top-level expressions by wrapping
each non-assignment expression inside a call to function `obj_filter`,
which can then dispatch on the type of the expression.
"""
# Numeric values, LineNumber expressions etc. will never be useful to us:
expr_filter(f::Function, e::Any) = e
# expr_filter(e::LineNumberNode) = (global last_linenumber = e)
# A Symbol might be an AbstractSolid variable name, so we add it:
expr_filter(f::Function, e::Symbol) = :($f($e))

# we add any top-level expressions, except assignments:
expr_filter(f::Function, e::Expr) = expr_filter(f, Val(e.head), e)
expr_filter(f::Function, ::Val, e::Expr) = :($f($e))
expr_filter(f::Function, ::Val{:(=)}, e::Expr) = e

# if several expressions are semicolon-chained, we pass this down
expr_filter(f::Function, ::Val{:toplevel}, x::Expr) =
	Expr(:toplevel, expr_filter.(f, x.args)...)

# Leaf solids<<<1
# PrimitiveSolid<<<2
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
"""
struct PrimitiveSolid{S,D,T,X} <: AbstractSolid{D,T}
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

Ortho{S,D,T} = PrimitiveSolid{S,D,T,
	@NamedTuple{size::Vec{D,T}, center::Bool}}
@inline _infer_type(::Type{<:Ortho};size, center) = real_type(size...)
@inline (H::Type{<:Ortho{S,D}})(v::AbstractVector;center=false) where{S,D}=
	H(size=Vec{D}(v), center=center)
@inline (H::Type{<:Ortho{S,D}})(v::Real;kwargs...) where{S,D} =
	H(SVector{D}(v for i in 1:D); kwargs...)

"""
    Square{T}
    Square(size, [center=false])
    Square([size1, size2], [center=false])

A square (or a rectangle) with coordinate type `T`, parallel to the axes.
If `center=true` is passed, then the origin is the center-of-mass;
otherwise it is the lower-left corner.
"""
Square = Ortho{:square,2}
"""
    Cube{T}
    Cube(size, [center=false])
    Cube([size1, size2, size3], [center=false])

A cube (or a rectangle parallelepiped) with coordinate type `T`, parallel
to the axes. If `center=true` is passed, then the origin is the center of
mass; otherwise it is the lower-front-left corner.
"""
Cube = Ortho{:cube,3}

Ball{S,D,T} = PrimitiveSolid{S,D,T, @NamedTuple{radius::T}}
@inline _infer_type(::Type{<:Ball};radius) = real_type(radius)
# groumpf, no dispatch on named arguments...
# @inline (T::Type{<:Ball})(;diameter) = T(one_half(diameter))
@inline (H::Type{<:Ball})(r) = H(radius=r)

"""
    Circle{T}
    Circle(r)

A circle with radius `r`, centered at the origin.
"""
Circle = Ball{:circle,2}
"""
    Sphere{T}
    Sphere(r)

A sphere with radius `r`, centered at the origin.
"""
Sphere = Ball{:sphere,3}

"""
    Cylinder(h, r1, r2 [, center=false])
    Cylinder(h, (r1, r2) [, center=false])
    Cylinder(h, r [, center=false])

**Warning:** `Cylinder(h,r)` is interpreted as `Cylinder(h,r,r)`,
not `(h,r,0)` as in OpenSCAD.
"""
Cylinder{T} = PrimitiveSolid{:cylinder,3,T,
	@NamedTuple{r1::T,r2::T,h::T,center::Bool}}
@inline _infer_type(::Type{<:Cylinder}; h,r1,r2,center) = real_type(h, r1, r2)
@inline (H::Type{<:Cylinder})(h, r1,r2; center=false) =
	H(r1=r1, r2=r2, h=h, center=center)
@inline (H::Type{<:Cylinder})(h, r; kwargs...) = H(h, r, r; kwargs...)
@inline (H::Type{<:Cylinder})(h, r::Union{AbstractVector, Tuple}; kwargs...) =
	H(h, r...; kwargs...)

"""
    Polygon{T}
    Polygon([point1, point2, ...])
    Polygon(point1, point2, ...)

A simple, closed polygon enclosed by the given vertices.
"""
Polygon{T} = PrimitiveSolid{:polygon,2,T,
	@NamedTuple{points::Path{2,T}}}
@inline _infer_type(::Type{<:Polygon}; points) = real_type(eltype.(points)...)
@inline (T::Type{<:Polygon})(points::AnyPath{2}) = T(points=points)
# @inline (T::Type{<:Polygon})(points::AnyVec{2,<:Real}...) = T([points...])

@inline points(p::Polygon) = p.points
"""
    Surface([points...], [faces...])
"""
Surface{T,V} = PrimitiveSolid{:polyhedron,3,T,
  @NamedTuple{points::Path{3,T}, faces::Vector{V}}}
@inline _infer_type(::Type{<:Surface}; points, faces) =
	real_type(eltype.(points)...)
# OpenSCAD is zero-indexed
@inline function scad_param_value(io::IO,
		::Val{:polyhedron}, ::Val{:faces}, faces)
	print(io, "[")
	join(io, [f .- 1 for f in faces], ",")
	print(io, "]")
end
@inline (::Type{Surface{T}})(;points, faces) where{T} =
	Surface{T,eltype(faces)}(points=points, faces=faces)
@inline (T::Type{<:Surface})(points, faces) = T(points=points, faces=faces)

abstract type LeafSolid{D,T} <: AbstractSolid{D,T} end
# Neutral solids<<<2
"""
    NeutralSolid{D,T}

A convenience type representing either an empty or full solid.
This exists mostly to provide a neutral element
for `union()` and `intersect()` operators, hence the name.
In those cases, it is impossible to know in advance
the dimension of the returned solid;
hence, as an exception to the general rule,
the `D` type parameter is either the symbol `:empty` or `:full`.
Since neutral objects are removed at compile-time
from corresponding CSG operations,
this should have no influence on the dimension of a top-level object.

The `T` type parameter is always `Bool`.
"""
struct NeutralSolid{D,T} <: LeafSolid{D,T} end
@inline scad_name(::NeutralSolid{D}) where{D} = D
# @inline parameters(::IO, ::NeutralSolid) = nothing

macro define_neutral(op, what, result)
	W = QuoteNode(what)
	F=:(Solids.$op)
	quote
	@inline $F(neutral::AbstractSolid, absorb::NeutralSolid{$W}) = $result
	@inline $F(absorb::NeutralSolid{$W}, neutral::AbstractSolid) = $result
	@inline $F(x::NeutralSolid{$W}, ::NeutralSolid{$W}) = x
	end
end
@inline union() = NeutralSolid{:empty,Bool}()
@inline intersect() = NeutralSolid{:full,Bool}()

# these are necessary for the following macros:
function minkowski end
function hull end
@define_neutral union empty neutral
@define_neutral union full  absorb
@define_neutral intersect empty absorb
@define_neutral intersect full  neutral
@define_neutral minkowski empty absorb
@define_neutral minkowski full  absorb
@define_neutral hull empty neutral
@define_neutral hull full  absorb

# Somewhat reduce type I/O clutter<<<1
Base.show(io::IO, ::Type{_FIXED}) = print(io, "_FIXED")
Base.show(io::IO, ::Type{Vec}) = print(io, "Vec")
Base.show(io::IO, ::Type{Vec{D}}) where{D} = print(io, Vec,"{$D}")
Base.show(io::IO, ::Type{Vec{D,T}}) where{D,T<:Number} = print(io, Vec,"{$D,$T}")
Base.show(io::IO, ::Type{Path}) = print(io, "Path")
Base.show(io::IO, ::Type{Path{D}}) where{D} = print(io, Path,"{$D}")
Base.show(io::IO, ::Type{Path{D,T}}) where{D,T<:Number} = print(io, Path,"{$D,$T}")
for type in (:Square, :Cube, :Circle, :Sphere, :Cylinder, :Polygon)
	@eval begin
		Base.show(io::IO, ::Type{$type}) = print(io, $(String(type)))
		Base.show(io::IO, ::Type{$type{T}}) where{T} =
			print(io, $(String(type)), "{$T}")
	end
end

#DerivedSolid<<<1
# https://www.usenix.org/legacy/event/usenix05/tech/freenix/full_papers/kirsch/kirsch.pdf
"""
		DerivedSolid{D,S}

A type representing CSG operations on solids. `D` is the dimension and
`S` is a symbol representing the operation (union, intersection etc.)
"""
struct DerivedSolid{D,S,T} <: AbstractSolid{D,T}
	children::Vector{<:AbstractSolid}
	DerivedSolid{D,S,T}(v::AbstractVector{<:AbstractSolid}) where{D,S,T} =
		new{D,S,T}(v)
	DerivedSolid{D,S}(v::AbstractVector{<:AbstractSolid}) where{D,S} =
		new{D,S,promote_type(eltype.(typeof.(v))...)}(v)
end
@inline DerivedSolid{D,S,T}(v::AbstractSolid...) where{D,S,T} =
	DerivedSolid{D,S,T}([v...])
# the previous method does *not* match the zero-length case:
@inline DerivedSolid{D,S,T}() where{D,S,T} =
	DerivedSolid{D,S,T}(AbstractSolid[])
@inline children(s::DerivedSolid) = s.children
@inline scad_name(::DerivedSolid{D, S}) where{D,S} = S
@inline scad_name(::DerivedSolid{D, :intersection}) where{D} = :intersection

# make operators associative; see definition of + in operators.jl
for op in (:union, :intersect, :minkowski, :hull)
	Q=QuoteNode(op)
	# union, intersection, minkowski are trivial on single objects:
	op != :hull &&  @eval ($op)(a::AbstractSolid) = a
	@eval begin
	# all of these are associative:
	# we leave out the binary case, which will be defined on a case-by-case
	# basis depending on the operators (see below).
#		($op)(a::AbstractSolid, b::AbstractSolid) =
#			DerivedSolid{$Q}([unroll(a, Val($Q)); unroll(b, Val($Q))])
	($op)(a::AbstractSolid, b::AbstractSolid, c::AbstractSolid, x...) =
		Base.afoldl($op, ($op)(($op)(a,b),c), x...)
	end
end
"""
    union(s::AbstractSolid...)
    s1 ∪ s2

Represents the union of given solids.
"""
@inline union(a1::AbstractSolid{D1}, a2::AbstractSolid{D2}) where{D1, D2} =
	DerivedSolid{max(D1,D2), :union}(unroll2(a1, a2, Val(:union)))
"""
    intersect(s::AbstractSolid...)
    s1 ∩ s2

Represents the intersection of given solids.
"""
@inline intersect(a1::AbstractSolid{D1}, a2::AbstractSolid{D2}) where{D1, D2} =
	DerivedSolid{min(D1,D2), :intersection}(unroll2(a1, a2, Val(:intersection)))
"""
    minkowski(s::AbstractSolid...)

Represents the Minkowski sum of given solids.
"""
@inline minkowski(a1::AbstractSolid{D1}, a2::AbstractSolid{D2}) where{D1, D2} =
	DerivedSolid{max(D1,D2), :minkowski}(unroll2(a1, a2, Val(:minkowski)))
"""
    hull(s::AbstractSolid...)

Represents the convex hull of given solids.
"""
@inline hull(s::AbstractSolid...) =
	DerivedSolid{maximum(dim.(s)), :hull}(
		[unroll(t, Val.((:hull, :union))...) for t in s])
# @inline hull(a1::AbstractSolid{D1}, a2::AbstractSolid{D2}) where{D1, D2} =
# 	DerivedSolid{max(D1,D2), :hull}(unroll2(a1, a2, Val.((:hull, :union))...))

"""
		unroll(x::AbstractSolid, Val(sym1), Val(sym2)...)

Returns either `[x]` or, if `x` is a `DerivedSolid` matching one of the
symbols `sym1`, `sym2`..., `children(x)`. (This helps reduce nesting).
"""
@inline unroll(s::AbstractSolid, ::Val, tail...) = unroll(s, tail...)
@inline unroll(s::AbstractSolid) = s
@inline unroll(s::DerivedSolid{D, S}, ::Val{S}, tail...) where{D, S} =
	children(s)
@inline unroll2(s::AbstractSolid, t::AbstractSolid, tail...) =
	[unroll(s, tail...); unroll(t, tail...)]

# minus operator is binary:
"""
    Solids.difference(s1, s2, ...)

Represents the difference `s1 ∖ (s2 ∪ ..)`.
"""
@inline difference(x::AbstractSolid{D1}, y::AbstractSolid...) where{D1} =
	DerivedSolid{D1, :difference}(
		[unroll(x, Val(:difference)); unroll.(y, Val(:union))...])
# added interface: difference([x...], [y...])
@inline difference(x::AbstractVector{<:AbstractSolid},
				y::AbstractVector{<:AbstractSolid}) =
	DerivedSolid{foldr(max,dim.(x);init=3),:difference}(union(x...), y...)

# General transforms<<<1
# Curry<<<2
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

# Transform type<<<2
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
struct Transform{S,D,T,X} <: AbstractSolid{D,T}
	child::AbstractSolid{D,T}
	data::X
	Transform{S}(data::X, child::AbstractSolid{D,T}) where{S,X,D,T} =
		new{S,D,T,X}(child, data)
end
# default values for I/O:
# (parameters in `data` are assumed to be stored in a NamedTuple).
@inline children(f::Transform) = [f.child]
@inline scad_name(f::Transform{S}) where{S} = S
@inline parameters(f::Transform) = f.data
@inline (T::Type{Transform{S}})(f, s::AbstractSolid...) where{S} =
	T(f, union(s...))
@inline (T::Type{Transform{S}})(f, s::Vector{<:AbstractSolid}) where{S} =
	T(f, s...)
@inline (T::Type{Transform{S}})(f) where{S} = Curry{S}((s...)->T(f, s...))
# We can extract the `f` value from the above in the following way:
"""
    extract(c::Curry)

Given a `Curry` object with function `s -> Transform{...}(f, s)`,
recovers the parameter `f`.
"""
function extract end
@inline (T::Type{<:Transform})(f, ::typeof(extract)) = f
@inline extract(c::Curry) = c.f(extract)

# SetParameters<<<2
SetParameters = Transform{:parameters}

@inline set_parameters(s...; parameters...) =
	SetParameters(parameters.data, s...)

function scad(io::IO, s::SetParameters, spaces::AbstractString = "")
	println(io, spaces, "{ // SetParameters")
	for (i,j) in pairs(s.data)
		println(io, spaces, "// ", i, "=", j, ";")
	end
	scad(io, s.child, spaces*" ")
	println(io, spaces, "} // SetParameters")
end
# Color<<<2
Paint = Transform{:color}

"""
    color(c::Colorant, s...)
    color(c::AbstractString, s...)
    color(c::AbstractString, α::Real, s...)
    color(c) * s...

Paints objects `s...` in the given color.
"""
@inline color(c::Colors.Colorant, s...) = Paint(c, s...)
@inline color(c::AbstractString, s...) =
	color(parse(Colors.Colorant, c), s...)
@inline color(c::AbstractString, a::Real, s...) =
	color(Colors.coloralpha(parse(Colors.Colorant, c), a), s...)

@inline scad_parameters(io::IO, c::Colors.Colorant) =
	print(io, round.(Float64.([red(c), green(c), blue(c), alpha(c)]), digits=3))

# Linear extrusion<<<2
LinearExtrude = Transform{:linear_extrude}
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

# Rotational extrusion<<<2
RotateExtrude = Transform{:rotate_extrude}
@inline rotate_extrude(angle::Real, s...; center=false) =
	RotateExtrude((angle=angle, center=center), s...)
@inline rotate_extrude(s...) = rotate_extrude(360, s...)
# Offset
Offset = Transform{:offset}
@inline offset(r::Real, s...; join=:round, miter_limit=2.) =
	Offset((r=r, join=join, miter_limit=miter_limit), s...)
@inline scad_parameters(io::IO, s::Offset) =
	scad_parameters(io, s, Val(parameters(s).join), parameters(s))
@inline scad_parameters(io::IO, ::Offset, ::Val{:round}, param) =
	scad_parameters(io, (r=param.r,))
@inline scad_parameters(io::IO, ::Offset, ::Val{:miter}, param) =
	scad_parameters(io, (delta=param.r, chamfer=false,))
@inline scad_parameters(io::IO, ::Offset, ::Val{:square}, param) =
	scad_parameters(io, (delta=param.r, chamfer=true,))

# Affine transforms<<<1
# Affine type<<<2
# type and constructors<<<3
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

@inline apply(f::Affine, v) = f.a * v + f.b
@inline apply(f::Affine, p::Path) = [apply(f, v) for v in p]

# neutral elements: <<<3
# this could in principle be defined for Val{T} where{T}, but we try
# to pirate a minimum number of functions in Base.
@inline Base.:*(::Val{true}, v) = v
@inline Base.:*(a, v::Val{false}) = v
@inline Base.:+(v, ::Val{false}) = v
@inline Base.:-(v, ::Val{false}) = v

# I/O: <<<3
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

# Reflections<<<2
"""
		Reflection

A type containing compressed information for an orthogonal reflection
matrix. Inspired by the `Diagonal` type.
"""
struct Reflection{D,T,V<:AbstractVector{T}} <: StaticMatrix{D,D,T}
	axis::V
	@inline Reflection{D,T,V}(v::V) where {D,T,V<:AnyVec{D,T}} =
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

# AffineTransform<<<2
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

# these two functions are now enough to pre-compose all affine transforms
# *before* applying them to objects:
@inline assoc(::Curry{:multmatrix}, ::Curry{:multmatrix}) = :left
@inline function compose(c1::Curry{:multmatrix}, c2::Curry{:multmatrix})
	(f1, f2) = (extract(c1), extract(c2))
	mult_matrix(f1.a*f2.a, f1.a*f2.b + f1.b)
end

# Translation, scaling, rotation, mirror<<<2
# FIXME change this '1' to a compile-time constant?
"""
    translate(v, s...)
    translate(v) * s

Translates solids `s...` by vector `v`.
"""
@inline translate(v::AnyVec, s...) = mult_matrix(1, v, s...)
"""
    scale(a, s...; center=0)
    scale(a; center=0) * s
Scales solids `s` by factor `a`. If `center` is given then this will be
the invariant point.

`a` may also be a vector, in which case coordinates will be multiplied by
the associated diagonal matrix.
"""
@inline scale(a::Real, s...; kwargs...) = mult_matrix(a, s...; kwargs...)
@inline scale(a::AnyVec, s...; kwargs...) =
	mult_matrix(Diagonal(a), s...; kwargs...)
"""
    mirror(v, s...; center=0)
    mirror(v; center=0) * s

Reflection with axis given by the hyperplane normal to `v`.
If `center` is given, then the affine hyperplane through this point will
be used.
"""
@inline mirror(v::AnyVec, s...; kwargs...) =
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

# Operators<<<2
@inline +(v::AnyVec, x::AbstractSolid) = translate(v,x)
@inline +(x::AbstractSolid, v::AnyVec) = translate(v,x)

# this purposely does not define a method for -(x::AbstractSolid).
@inline Base.:-(x::AbstractSolid, y::AbstractSolid, tail...) =
	difference(x, [y, tail...])
@inline Base.:-(x::AbstractSolid{D}) where{D} = difference(intersect(), x)
@inline Base.:-(x::AbstractVector{<:AbstractSolid},
                y::AbstractVector{<:AbstractSolid}) = difference(x, y)

# @inline *(f::AbstractAffineMap, x::AbstractSolid) = mult_matrix(f, x)
@inline *(s::Union{Real,AnyVec{3}}, x::AbstractSolid) = scale(s,x)

⋃ = Base.union
⋂ = Base.intersect

# Attachments<<<1
# Anchor system<<<2
"""
		find_anchor(x::AbstractSolid, name)

Returns the anchor (an affine rotation) found for the given `name` for
the solid `x`.

Some types of `name` that can be used include:
 - a symbol: either one of the six standard directions (`:left`, `:right`,
	 `:front`, `:back`, `:top`, `:bottom`, `:center`)
	 or a custom-defined label (TODO);
	 (*Note:* for 2-dimensional solids, `:bottom` is equivalent to `:front`
	 and `:top` is equivalent to `:back`)
 - a list (tuple) of symbols, which is interpreted as the sum of the
	 corresponding directions;
 - for standard convex bodies, a vector of the same dimension as `x` is
	 normalized to a point at the boundary of `x` (see below);
 - a way to designate a point at the boundary of `x` (see below).
"""
@inline find_anchor(x::AbstractSolid, labels::NTuple{N,Symbol}) where{N} =
	find_anchor(x, sum([ labeled_anchor(x, l) for l in labels]))
@inline function find_anchor(x::AbstractSolid, label::Symbol)
	y = labeled_anchor(x, label)
	y isa Missing && error("No anchor named '$label' found in $(scad_name(x))")
	return find_anchor(x, y)
end

default_positions = (
	left  = SA[-1,0,0],
	right = SA[+1,0,0],
	front = SA[0,-1,0],
	back  = SA[0,+1,0],
	bot	  = SA[0,0,-1],
	bottom= SA[0,0,-1],
	top	  = SA[0,0,+1],
	center= SA[0,0,0],
)
@inline labeled_anchor(x::AbstractSolid{3}, label::Symbol) =
	get(default_positions, label, missing)
@inline labeled_anchor(x::AbstractSolid{2}, label::Symbol) =
	_labeled_anchor_3to2(get(default_positions, label, missing))
@inline _labeled_anchor_3to2(::Missing) = missing
@inline _labeled_anchor_3to2(v::StaticVector{3}) =
# for 2d objects, allow (:left..:right, :bot..:top) as anchor names:
	(v[1] == 0 && v[2] == 0 && v[3] != 0) ? SA[v[1], v[3]] : SA[v[1], v[2]]
# Define named anchors <<<2
# first column is translation, second (if existing) rotation,
# spin is angle in 2d
struct AnchorData{D,T,R}
	origin::Vec{D,T}
	direction::R
	spin::T
	@inline AnchorData{3,T}(o::AnyVec{3}, r::AnyVec{3}, s::Real) where{T} =
		new{3,T,Vec{3,T}}(o, r, s)
	@inline AnchorData{2,T}(o::AnyVec{2}, s::Real) where{T} =
		new{2,T,Nothing}(o, nothing, s)
end

@inline AnchorData{3,T}(v::AnyVec{3}) where{T} =
	AnchorData{3,T}(v,SA[0,0,1], zero(T))
@inline AnchorData{3,T}(data::Tuple{<:AnyVec{3},<:AnyVec{3}}) where{T} =
	AnchorData{3,T}(data[1], data[2], zero(T))
@inline AnchorData{3,T}(data::Tuple{<:AnyVec{3},<:AnyVec{3},<:Real}) where{T} =
	AnchorData{3,T}(data[1], data[2], T(radians(data[3])))
@inline text(x::AnchorData{3}) =
	"origin=$(Float64.(x.origin)) direction=$(Float64.(x.direction)) spin=$(Float64(x.spin))"

@inline AnchorData{2,T}(v::AnyVec{2}) where{T} =
	AnchorData{2,T}(v, zero(T))
@inline AnchorData{2,T}(data::Tuple{<:AnyVec{2},<:Real}) where{T} =
	AnchorData{2,T}(data[1], T(radians(data[2])))
@inline text(x::AnchorData{2}) =
	"origin=$(Float64.(x.origin)) angle=$(Float64(x.spin))"

@inline affine(x::AnchorData) = AffineMap(linear(x), x.origin)
@inline linear(x::AnchorData{3}) =
	rotation_between(SA[0,0,1], x.direction) * RotZ(x.spin)
@inline rotation(x::AnchorData{2}) = Angle2d(x.spin)


"""
		struct NamedAnchors

Wraps an object, adding symbolic anchors to it.
"""
struct NamedAnchors{D,T} <: AbstractSolid{D,T}
	child::AbstractSolid{D,T}
	anchors::Dict{Symbol, AnchorData{D,T}}
end
"""
		named_anchors(x, label=>anchor...)

Adds symbolic anchors to an object. `anchor` may be either

 - a vector: the associated anchor is a translation.
 - a pair (origin, direction): the associated anchor is the affine
	 rotation to given direction.
 - a triple (origin, direction, spin).
 - (in 2d) a pair (origin, angle).
"""
@inline named_anchors(x::AbstractSolid{D,T},
	a::Pair{Symbol}...) where{D,T} =
	NamedAnchors{D,T}(x, Dict(k => AnchorData{D,T}(v) for (k,v) in a))

function scad(io::IO, x::NamedAnchors, spaces::AbstractString = "")
	println(io, spaces, "// Object with named anchors:")
	for (label, anchor) in x.anchors
		println(io, spaces, "// $label: $(text(anchor))")
	end
	scad(io, x.child, spaces)
end

function labeled_anchor(x::NamedAnchors, label::Symbol)
	y = get(x.anchors, label, missing)
	if y isa Missing
		return get(default_anchors, label, missing)
	end
end
# Anchors for convex bodies <<<2
"""
		find_anchor(x::Square, position::SVector{2})
		find_anchor(x::Circle, position::SVector{2})
		find_anchor(x::Cube, position::SVector{3})
		find_anchor(x::Cylinder, position::SVector{3})
		find_anchor(x::Sphere, position::SVector{3})

For a convex body, the anchor corresponding to `position` has its origin
at the surface fo the body and maps the unit vector `[0,1]` (resp.
`[0,0,1]` in dimension 3) to the normal vector at this position.

If `position` is zero, then the translation to the center of the body is
returned.

The rotation is computed using `Rotations.rotation_between`. This returns
a well-defined result even for the rotation mapping `[0,0,1]` to
`[0,0,-1]`.
"""
function find_anchor(x::Ortho{D}, pos::Vec{D,<:Real}) where{D}
	center = ~x.center*one_half(x.size) # the center of the square/cube
	if iszero(pos) return Translation(center) end
	maxc = findmax(abs.(pos))[2] # max abs value of a coordinate
	v = sum(pos[abs.(pos) .= maxc]), # sum of coordinates with max abs value

	# we don't need to normalize v since `rotation_between` does it for us:
	AffineMap(rotation_between(SA[0,0,1], v), center + v .* x.size)
end

function find_anchor(x::Circle, pos::Vec{2,<:Real})
	if iszero(pos) return Translation(SA[0,0]) end
	p1 = pos / sqrt(pos[1]^2+pos[2]^2)
	AffineMap(SA[p1[2] p1[1]; -p1[1] p1[2]], radius(x)*p1)
end

function find_anchor(x::Sphere, pos::Vec{3,<:Real})
	if iszero(pos) return Translation(SA[0,0,0]) end
	return AffineMap(rotation_between(SA[0,0,1], pos), pos*x.radius / norm(pos))
end

function find_anchor(x::Cylinder, pos::Vec{3,<:Real})
	center = ~x.center*one_half(x.h)
	if iszero(pos) return Translation(center) end
	r = sqrt(pos[1]*pos[1]+pos[2]*pos[2]) # convert to cylindrical coords
	if pos[3]*x.r2 > r # top face: normalize to pos[3]==1
		return Translation(SA[pos[1]/pos[3], pos[2]/pos[3], center+one_half(x.h)])
	elseif pos[3]*x.r1 < -r # bottom face: normalize to pos[3]==-1
		return AffineMap(rotation_between(SA[0,0,1], SA[0,0,-1]),
			SA[-pos[1]/pos[3], -pos[2]/pos[3], center-one_half(x.h)])
	end
	# the line equation is 2r = (r2-r1) z + (r1+r2)
	r3 = one_half(x.r1+x.r2+pos[3]*(x.r2-x.r1)) # radius at given z
	# in cyl coordinates, the contact point is (r=r3, z=pos[3]*h/2+center)
	p = SA[pos[1]*r3/r, pos[2]*r3/r, center + one_half(x.h)*pos[3]]
	# normal vector is 2 dr = (r2-r1) dz
	n = SA[2*pos[1]/r, 2*pos[2]/r, x.r2-x.r1]
	AffineMap(rotation_between(SA[0,0,1], n), p)
end

# Coordinates on circle & sphere <<<2
"""
		find_anchor(x::Circle, angle::Real)
Returns anchor at point `(-sin(angle), cos(angle))` (trig. orientation,
starting at top of circle; angle in **degrees**) with outwards normal vector
(the start at top guarantees that angle 0 preserves upwards-pointing
vector).
"""
function find_anchor(x::Circle{T}, angle::Real) where{T}
# 	a = T(radians(angle))
	(s, c) = T.(sincosd(a))
	AffineMap(SA[c -s;s c], SA[-radius(x)*s, radius(x)*c])
end
"""
		find_anchor(x::Sphere, (latitude, longitude))

Returns normal vector to sphere at this position (angles in **degrees**).
"""
function find_anchor(x::Sphere{T}, angle::AnyVec{2,<:Real}) where{T}
# 	(lat, lon) = T.(radians.(angle))
	(s1, c1) = T.(sincosd(lat))
	(s2, c2) = T.(sincosd(lon))
	r = x.radius
	AffineMap(SA[s1*c2 -s2 c1*c2;s1*s2 c2 c1*s2; -c1 0 s1],
						SA[r*c1*c2, r*c1*s2, r*s1])
end


# attach <<<2
"""
		attach(parent, {:label => child}...)

Moves (by rotations) all children so that their anchor matches the
anchors of `parent` defined by the given labels.
"""
function attach(parent::AbstractSolid, list::Pair{<:Any,<:AbstractSolid}...)
	union(parent, [ attach_at(parent, pos, child) for (pos,child) in list]...)
end
function attach_at(parent::AbstractSolid{D}, label, child::AbstractSolid) where{D}
	m = find_anchor(parent, label)
	mult_matrix(m, child)
end

# anchor <<<2
"""
		half_turn(q::UnitQuaternion)

Returns the half-turn rotation with same axis as `q`.
"""
@inline half_turn(q::Rotations.UnitQuaternion) =
	Rotations.UnitQuaternion(0, q.x, q.y, q.z) # the constructor normalizes this
"""
		half_turn_complement(q::Rotation{3})

Returns the unique rotation `r` such that `qr=rq` is a half-turn.
"""
@inline half_turn_complement(q::Rotations.Rotation{3}) =
	inv(q)*half_turn(q)

"""
		anchor(solid, label)

Translates the solid so that the anchor with name `label` is at origin
and the corresponding anchor vector points inward.
"""
function anchor(x::AbstractSolid, label)
# Ax + B maps [0,0,0] to anchor point p, and e_z=[0,0,1] to anchor vec v
# namely: A e_z = v, B = p
# we want to map p to [0,0,0] and -v to e_z
# i.e. A' v = - e_z and A' p + B' = 0, or B' = -A' p = -A' B
#

	m = find_anchor(x, label)
	a1= half_turn_complement(linear(m))
# ⚠ -a1*m is interpreted as (-a1)*m, and a1 is a quaternion ⇒ -a1≡a1 (as
# a rotation)
	mult_matrix(AffineMap(a1, a1*-translation(m)), x)
end

# # Reduction of mult_matrix operations <<<1
# @inline affine_reduce(m::AbstractAffineMap, x::LeafSolid) = mult_matrix(m, x)
# @inline affine_reduce(u::LinearMap{<:Diagonal}, x::Square) =
# 	Square([ u.linear.diag[i] * @inbounds x.size[i] for i in 1:2 ], x.center)
# @inline affine_reduce(u::LinearMap{<:Union{UniformScaling,Diagonal}},
# 	x::Union{Square,Cube}) = (typeof(x))(linear(u) * x.size, x.center)
# @inline affine_reduce(u::LinearMap{<:UniformScaling},
# 	x::Union{Circle,Sphere}) =
# 	(typeof(x))(linear(u) * x.r, x.frag)
# @inline affine_reduce(u::AbstractAffineMap, x::Polygon) =
# 	Polygon(apply(u, points(x)), x.path, x.convexity)
# @inline affine_reduce(u::AbstractAffineMap, x::Surface) =
# 	Surface(apply(u, points(x)), x.faces, x.convexity)
# @inline affine_reduce(m::AbstractAffineMap, x::T) where{T<:DerivedSolid} =
# 	T([ affine_reduce(m, y) for y in children(x) ])
# function affine_reduce(m1::AbstractAffineMap, m2::MultMatrix)
# 	m3 = compose(m1, m2.m)
# 	union((affine_reduce(m3, x) for x in children(m2))...)
# 	# TODO
# end
# affine_reduce(x::AbstractSolid) = affine_reduce(LinearMap(I), x)
#————————————————————— Meshing —————————————————————————————— <<<1
#>>>1
# Accuracy and precision <<<1
# Accuracy is the absolute deviation allowed.
# Default value is 2.0 (from OpenSCAD `$fs`), interpreted as 2mm.
#
# Precision is the relative deviation allowed.
# Default value is 0.02 (1-cos(180°/`$fa`)).

_DEFAULT_PARAMETERS = (accuracy = 2.0, precision = .02)

"""
    sides(radius, parameters, θ = 360°)

Returns the number of sides used to draw a circle (arc) of given angle.
The base value `n` is given by the minimum of:
 - accuracy: each side (length `θ r/n`) must not be smaller than accuracy,
   or n = θ r / accuracy;
 - precision: deviation=1/cos(θ/2n)≈ θ²/8n²,
	 or n ≈ θ/√(8 * precision)
"""
function sides(r::Real, parameters::NamedTuple, angle::AnyAngle = 360°)
	θ = radians(angle)
	acc = θ*r / (parameters.accuracy)
	pre = θ/sqrt(8*parameters.precision)
	base = min(acc, pre)
	# a circle always has at least 4 sides
	return round(Int, max(4, base), RoundUp)
end

"""
    unit_n_gon(T::Type, n::Int)
		unit_n_gon(r, n)

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
@inline unit_n_gon(r, parameters) =
	r*unit_n_gon(real_type(r), sides(r, parameters))

# Clipper.jl interface: clip, offset, simplify<<<1
# This is the only section in this file which contains code directly
# related to `Clipper.jl`. The entry points to this section are the
# functions `clip` and `offset` defined below.
# Types<<<2
# default number of bits for Clipper types
const _CLIPPER_BITS = FixedPointNumbers.nbitsfrac(_FIXED)
# this must be a 64-bit type, even if _FIXED is modified:
const _CLIPPER_FIXED = Fixed{Int64,_CLIPPER_BITS}
# constants
const _CLIPPER_ENUM = (#<<<
	clip=(
		union       =Clipper.ClipTypeUnion,
		intersection=Clipper.ClipTypeIntersection,
		difference  =Clipper.ClipTypeDifference,
		xor         =Clipper.ClipTypeXor,
	),
	ends=(
		closed=Clipper.EndTypeClosedPolygon,
		square=Clipper.EndTypeOpenSquare,
		round =Clipper.EndTypeOpenRound,
		butt  =Clipper.EndTypeOpenButt,
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
)#>>>

"""
    to_clipper(OriginalType, ...)
		from_clipper(OriginalType, ...)

Converts stuff (numbers, vectors, paths...) to and from `Clipper.jl` types.
"""
@inline to_clipper(::Type{<:Real}) = _CLIPPER_FIXED
@inline to_clipper(T::Type{<:FixedPoint{<:Int64}}) = T
@inline to_clipper(T::Type, x::Real) = reinterpret(convert(to_clipper(T), x))
@inline to_clipper(T::Type, v::AnyVec{2,<:Real}) =
	Clipper.IntPoint(to_clipper(T, v[1]), to_clipper(T, v[2]))
@inline to_clipper(T, p::AnyPath{2}) = [to_clipper(T, v) for v in p]
@inline to_clipper(T, p::Vector{<:AnyPath{2}}) =
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
	SA[from_clipper(T, p.X), from_clipper(T, p.Y)]
# paths...
@inline from_clipper(T::Type{<:Real}, p::Vector{Clipper.IntPoint}) =
	[ from_clipper(T, v) for v in p ]
@inline from_clipper(T::Type{<:Fixed{Int64}}, p::Vector{Clipper.IntPoint}) =
	reinterpret(Vec{2,T}, p)
# vectors of paths...
@inline from_clipper(T, polys::Vector{Vector{Clipper.IntPoint}}) =
	[ from_clipper(T, p) for p in polys ]
# Wrappers for Clipper calls<<<2
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

# Calls on Path values<<<2
function clip(op::Symbol,
		v1::AbstractVector{Path{2,T}},
		v2::AbstractVector{Path{2,T}};
		fill = :positive)::Vector{Path{2,T}} where {T}
	c = ClipperClip(T)
	add_paths!(c, v1, Clipper.PolyTypeSubject, true) # closed=true
	add_paths!(c, v2, Clipper.PolyTypeClip, true)

	f = _CLIPPER_ENUM.fill[fill]
	return execute(c, _CLIPPER_ENUM.clip[op], f, f)
end
function offset(v::AbstractVector{Path{2,T}}, r::Real;
		join = :round,
		ends = :closed,
		miter_limit = 2.,
		precision = 0.2
		)::Vector{Path{2,T}} where{T}
	c = ClipperOffset(T, miter_limit, precision)
	add_paths!(c, v, _CLIPPER_ENUM.join[join], _CLIPPER_ENUM.ends[ends])
	execute(c, r)
end
function offset(v::AbstractVector{Path{2,T}}, r::AbstractVector{<:Real};
		join = :round,
		ends = :closed,
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
#Convex hull<<<1
# 2d convex hull <<<2

"""
    convex_hull([2d points])

Returns the convex hull (as a vector of 2d points).
"""
function convex_hull(points::AbstractVector{<:Vec{2,T}}) where{T}
	M = hcat(Vector.(points)...)
# M is a matrix with the points as *columns*, hence the transpose
# below:
	PH = Polyhedra
	poly = PH.polyhedron(PH.vrep(transpose(M)), polyhedra_lib(T))
	PH.removevredundancy!(poly).points.points
end

# this is the version using MiniQhull: #<<<
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
# end#>>>
# Triangulate face <<<2
function triangulate_face_convex(face::Vector,
		pts::AbstractVector{<:Vec{3,T}}, map::AbstractVector{Int}) where{T}
	# project all points on two coordinates:
	# by default, (x,y), unless plane is vertical (face[3] == 0),
	proj = iszero(face[3]) ? ( iszero(face[2]) ? [2,3] : [1,3]) : [1,2]
	N = length(pts)
	pts1 = Matrix{T}(undef, N, 2)
	for (j, x) in pairs(proj)
		for (i, p) in pairs(pts)
			pts1[i, j] = p[x]
		end
	end
	Vec{3,Int}.(Triangle.basic_triangulation(pts1, map))
end
# 3d convex hull <<<2

"""
		convex_hull(x::AbstractSolid{3}...)

Returns the convex hull of the union of all the given solids, as a
`Surface` structure.
"""
@inline convex_hull(x::AbstractSolid{3}) =
	convex_hull(vcat([Vec{3}.(points(y)) for y in x]...))

function convex_hull(p::AbstractVector{<:Vec{3,T}}) where{T}
	M = hcat(Vector.(p)...)
	PH = Polyhedra
	poly = Polyhedra.polyhedron(vrep(transpose(M)), polyhedra_lib(T))
	R = removevredundancy!(poly)
	V = Vec{3,T}.(collect(PH.points(poly)))

	triangles = Vec{3,Int}[]
	for i in eachindex(halfspaces(poly)) # index of halfspace
		pts = incidentpointindices(poly, i) # vector of indices of points
		h = get(poly, i) # halfspace equation
		for t in triangulate_face_convex(h.a,
				[Vec{3,T}(get(poly, j)) for j in pts], [j.value for j in pts])
			(a,b,c) = (V[j] for j in t)
			k = det([b-a c-a h.a])
			push!(triangles, (k > 0) ? t : SA[t[1], t[3], t[2]])
		end
	end
	Surface(V, triangles)
end


# Minkowski sum<<<1
# Convolution of polygons<<<2
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

function convolution(p::AnyPath{2}, q::AnyPath{2})
	(np, nq) = (length(p), length(q))
	ep = [p[(i%np)+1]-p[i] for i in 1:np] # edges of p
	eq = [q[(i%nq)+1]-q[i] for i in 1:nq]
	j0 = 0
	newpoly = similar(p, 0)
	for ip in 1:np
		for iq in 1:nq
			iq0 = mod1(iq-1, nq)
			if circularcmp(eq[iq0], ep[ip], eq[iq], Val(:offset))
				push!(newpoly, p[ip]+q[iq])
				push!(newpoly, p[mod1(ip+1,np)]+q[iq])
			end
		end
	end
	newpoly
end
p⋆q = convolution(p, q)
# Minkowski sum of polygons and their unions <<<2
function minkowski(p::AnyPath{2}, q::AnyPath{2})
	r = convolution(p, q)
	return simplify([r]; fill=:nonzero)
end
function minkowski(vp::Vector{<:AnyPath{2}}, vq::Vector{<:AnyPath{2}})
	vr = vec([convolution(p, q) for p in vp, q in vq])
	global X = vr
	return simplify(vr; fill=:nonzero)
end

# 2d subsystem<<<1
# PolyUnion<<<2
# type and constructors from points<<<
"""
		PolyUnion

Represents a union of polygons. Each polygon is assumed to be simple and
ordered in trigonometric ordering.
"""
struct PolyUnion{T} <: LeafSolid{2,T}
	poly::Vector{Path{2,T}}
	@inline PolyUnion{T}(p::AbstractVector{<:AnyPath{2,T}}) where{T} =
		new{T}(Path{2,T}.(p))
end
@inline (U::Type{PolyUnion})(p::AbstractVector{<:AnyPath{2,T}}) where{T} =
		U{real_type(eltype.(eltype.(p))...)}(p)

@inline (U::Type{<:PolyUnion})(path::AnyPath{2}...) = U([path...])

@inline points(u::PolyUnion) = vcat(u.poly...)

# this is used to broadcast conversion for recursive conversion to PolyUnion:
@inline _convert(U::Type{<:PolyUnion}, l, parameters) =
	[ U(s, parameters) for s in l ]

# >>>
# I/O<<<
function scad(io::IO, u::PolyUnion{S}, spaces::AbstractString) where{S}
	print(io, spaces, "// union of $(length(u.poly)) polygon(s):\n")
	length(u.poly) != 1 && print(io, spaces, "union() {\n")
	for p in u.poly
		print(io, spaces, " polygon([")
		join(io, convert.(Vec{2,Float64}, p), ",")
		print(io, "]);\n")
	end
	length(u.poly) != 1 && print(io, spaces, "}\n")
end
#>>>
# Conversion from leaf 2d types<<<2
@inline PolyUnion(l::AbstractSolid{2}; kwargs...) =
	PolyUnion{real_type(eltype(l))}(l, merge(_DEFAULT_PARAMETERS, kwargs.data))

@inline (U::Type{<:PolyUnion})(x::Square, parameters) =
	U(Vec{2}.(points(x)))
@inline (U::Type{<:PolyUnion})(x::Circle, parameters) =
	U(Vec{2}.(points(x, parameters)))
@inline (U::Type{<:PolyUnion})(x::Polygon, parameters) =
# FIXME: simplify and define orientation
	U(Vec{2}.(points(x)))

@inline function (U::Type{<:PolyUnion})(f::AffineTransform{2}, parameters)
	child = U(f.child, parameters)
	return U([ apply(f.data, path) for path in child.poly ])
end
@inline (U::Type{<:PolyUnion})(s::SetParameters{2}, parameters) =
	U(s.child, merge(parameters, s.data))
# fall-back case (`color`, etc.):
@inline (U::Type{<:PolyUnion})(s::Transform{S,2}, parameters) where{S} =
	U(s.child, parameters)

function points(s::Square)
	# in trigonometric order:
	if s.center
		(a,b) = one_half(s.size)
		return [SA[-a,-b], SA[a,-b], SA[a,b], SA[-a,b]]
	else
		return [SA[0,0], SA[s.size[1],0], s.size, SA[0,s.size[2]]]
	end
end
@inline points(c::Circle, parameters) = unit_n_gon(c.radius, parameters)
# Reduction of CSG operations<<<2
@inline (clip(op, u::U...)::U) where{U<:PolyUnion} =
	reduce((a,b)->U(clip(op, a.poly, b.poly)), u)

# Set-wise operations:
@inline (U::Type{<:PolyUnion})(s::DerivedSolid{2,:union}, parameters) =
	clip(:union, _convert(U, s.children, parameters)...)

@inline (U::Type{<:PolyUnion})(s::DerivedSolid{2,:intersection}, parameters) =
	clip(:intersection, _convert(U, s.children, parameters)...)

function ((U::Type{<: PolyUnion})(s::DerivedSolid{2,:difference}, parameters)::U)
	length(s.children) == 1 && return U(s.children[1], parameters)
	L = _convert(U, s.children, parameters)
	r2= clip(:union, view(L,2:length(L))...)
	clip(:difference, L[1], r2)
end

# Convex hull:
function (U::Type{<:PolyUnion})(s::DerivedSolid{2,:hull}, parameters)
	pts = points.(_convert(U, s.children, parameters))
	U(convex_hull([pts...;]))
end

# Minkowski sum:
function (U::Type{<:PolyUnion})(s::DerivedSolid{2,:minkowski},
	parameters)::U
	reduce((a,b)->U(minkowski(a.poly, b.poly)),
		_convert(U, s.children, parameters))
end
# function _combine2(::Val{:minkowski}, a::PolyUnion{T}, b::PolyUnion{T}) where{T}
# 	# not implemented in Clipper.jl...
# end


# Offset and draw <<<2
"""
		offset(P::Polygon, u::Real; options...)

Offsets polygon `P` by radius `u` (negative means inside the polygon,
positive means outside). Options:

 - `join_type`: :round | :square | :miter
 - `miter_limit` (default 2.0)
"""
function offset(U::PolyUnion{T}, u::Real;
		join_type = :round,
		miter_limit::Float64 = 2.0,
		precision::Real = 0.2) where{T}

	global X=(clipper_type(T), precision)
	c = ClipperOffset(miter_limit, clipper_float(clipper_type(T), precision))
	add_paths!(c, U.poly, join_type, Clipper.EndTypeClosedPolygon)
	PolyUnion(execute(T, c, u))
end
@inline offset(x::AbstractSolid{2}, args...; kwargs...) =
	offset(PolyUnion(x), args...; kwargs...)

# Draw <<<2
"""
    draw(path, width; kwargs...)

    ends=:round|:square|:butt|:closed
    join=:round|:miter|:square
"""
function draw(path::Path{2,T}, width::Real;
		ends::Symbol = :round, join::Symbol = :round,
		miter_limit::Float64 = 2.0, precision::Real = 0.2) where{T}
	CT = clipper_type(T)
	RT = clipper_rettype(T)
	c = ClipperOffset(miter_limit, clipper_float(CT, precision))
	println("join=$join, round=$round")
	Clipper.add_path!(c, clipper_path(path),
		JoinType(Val(join)), EndType(Val(ends)))
	println("$(clipper_type(T)) $(CT(1.)); prec=$(Float64(CT(precision)))")
	ret = clipper_unpath.(RT, Clipper.execute(c, clipper_float(CT, width)/2))
	return PolyUnion(ret)
end

# Convex hull<<<2
# """
# 		convex_hull(x::AbstractSolid{2}...)
# 
# Returns the convex hull of the union of all the given solids, as a
# `PolyUnion` structure.
# """
@inline convex_hull(x::AbstractSolid{2}...; params=_DEFAULT_PARAMETERS) =
	convex_hull(PolyUnion(union(x...), params))

@inline convex_hull(u::PolyUnion) = convex_hull(Vec{2}.(points(u)))

# # 3d conversion<<<1
# # function points(s::Cube, parameters::NamedTuple)
# # 	if s.center
# # 		(a,b,c) = one_half(s.size)
# # 		return [SA[-a,-b,-c], SA[-a,b,-c], SA[a,b,-c], SA[a,-b,-c],
# # 			SA[-a,-b,c], SA[-a,b,c], SA[a,b,c], SA[a,-b,c]]
# # 	else
# # 		return [SA[0,0,0], SA[0,s.size[2],0],
# # 			SA[s.size[1],s.size[2],0], SA[s.size[1],0,0],
# # 			SA[0,0,s.size[3]], SA[0,s.size[2],s.size[3]],
# # 			s.size, SA[s.size[1],0,s.size[3]]]
# # 	end
# # end
# # function points(c::Cylinder, parameters)
# # 
# # 	p1 = unit_n_gon(c.r1, parameters)
# # 	p2 = unit_n_gon(c.r2, parameters)
# # 	h = one_half(c.h)
# # 	z = h*~c.center
# # 	return vcat([ [ p; h-z ] for p in p1], [[p; h+z ] for p in p2 ])
# # end
# 
# Extrusion <<<1
# Path extrusion <<<2
# poly_dist_type(::Type{Vec{2,T}}) where{T} = T
# poly_dist_type(::Type{Clipper.IntPoint}) = Int
distance2(x::Vec{2,T}, y::Vec{2,T}) where{T} = let z = x .- y
	z[1]*z[1] + z[2]*z[2]
end

# triangulate_between: triangulate between two parallel paths<<<
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
		poly1::AbstractVector{<:Path{D,T}},
		poly2::AbstractVector{<:Path{D,T}},
		start1::Int = 1, start2::Int = 1) where {D,T}
	println("T=$T")
	Big = typemax(T)
# 	Distance = poly_dist_type(T); Big = typemax(Distance)
	Triangle = SVector{3,Int}
	triangles = Triangle[]
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

	# we need a way to convert (poly, path, index) to integer index<<<
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
#>>>
	# computing diagonal distances to find the smallest one:<<<
	distance(pt, path, i) =
		i > length(path) ? Big : distance2(pt, path[i])

	closest(pt, poly, status) =
		findmin([distance(pt, poly[i], status[i]+1) for i in 1:length(poly)])
#>>>

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
end#>>>
# path_extrude<<<
"""
		path_extrude(path, poly, options...)

Extrudes the given polygon (a path of points forming a simple loop)
along the given path. Both arguments are provided as a
`Vector{SVector{2}}`.

Returns a `Surface` (defined by points and a triangulation).
"""
function path_extrude(path::Path{2,T},
	poly::Path{2};
	join = :round,
	miter_limit::Float64 = 2.0,
	precision::Float64 = 0.25,
	closed::Bool = true
	) where{T}

	N = length(poly)
	# offset_path is a vector of vector of paths
	offset_path = offset([path], [pt[1] for pt in poly],
		join = join, ends = closed ? :closed : :butt)
	global OP = offset_path
	# new_points is a flat list of all 3d points produced
	new_points = [[
		[ SA[pt[1], pt[2], poly[i][2]] for pt in [p...;] ]
		for (i, p) in pairs(offset_path)
	]...;]
	println("returning new_points:")

	# first index for each path
	first_face = cumsum([1; # initial
		map(p->sum(length.(p)), offset_path)])
	println("first_face=$first_face")

	triangles = map(1:N) do i
		i1 = mod1(i+1, N)
		triangulate_between(offset_path[i], offset_path[i1],
			first_face[i], first_face[i1])
		# XXX keep the last edge for closing the poly
	end
	# this completes the set of triangles for the tube:
	tube_triangles = vcat([ t[1] for t in triangles ]...)
	last_face = [ t[2][1] for t in triangles ]
	println("last_face=$last_face")
	# here we decide if it is closed or open
	# if open, triangulate the two facets
	# if closed, join them together
	if closed
		more_triangles = vcat(map(1:N) do i
			j = (i%N)+1
			[ SA[first_face[i], last_face[i], first_face[j]],
				SA[first_face[j], last_face[i], last_face[j]] ]
		end...)
		println("more_triangles=$more_triangles")
		tube_triangles = [ tube_triangles; more_triangles ]
	else
	# TODO: triangulate the surface
	# or, for now, close with two non-triangular facets...
		more_triangles = [ reverse(first_face), last_face ]
		println("more_triangles=$more_triangles")
	end
	Surface( new_points, tube_triangles )
end#>>>

# # Annotations <<<1
# abstract type AbstractAnnotation{D} end
# 
# struct DimensionArrow{D,T} <: AbstractAnnotation{D}
# 	center::Vec{D,T}
# 	vec::Vec{D,T}
# 	label::AbstractString
# 	fontsize::Float64
# 	offset::Float64
# end
# """
#     DiameterArrow{X}
# 
# Indicates that a diameter arrow should be drawn for the given object. The
# parameter `X` is a type indicating which type of arrow should be drawn.
# 
#  - `Circle`: parametrized by center (`Vec{2}`) and radius (scalar),
#  and possibly preferred orientation (vector if non-zero).
# 
#  - `Sphere`: parametrized by center (`Vec{3}`) and radius (scalar),
#   and possibly preferred orientation.
# 
#  - `Cylinder`: shows a circle in 3d space, parametrized by center (`Vec{3}`), normal vector (non-zero), radius (scalar), and preferred orientation (vector; should be in the circle plane).
# """
# struct DiameterArrow{X<:AbstractSolid,T,D,N} <: AbstractAnnotation{D}
# 	center::Vec{D,T}
# 	radius::T
# 	orientation::Vec{D,T}
# 	normal::N
# 	# inner constructors enforce the correct type for N
# 	DiameterArrow{Circle,T}(center, radius, orientation) where{T} =
# 		new{Circle, T, 2, Nothing}(center, radius, orientation, nothing)
# 	DiameterArrow{Sphere,T}(center, radius, orientation) where{T} =
# 		new{Sphere, T, 3, Nothing}(center, radius, orientation, nothing)
# 	DiameterArrow{Cylinder,T}(center, radius, orientation, normal) where{T} =
# 		new{Cylinder, T, 3, Vec{3,T}}(center, radius, orientation, normal)
# end
# 
# struct Annotate{D,T} <: AbstractSolid{D,T}
# 	annotations::Vector{<:AbstractAnnotation{D}}
# 	child::AbstractSolid{D,T}
# end
# # the offset is just a hint; we let the visualizer take care of using
# # this
# 
#
# Exports <<<1
export dim
export Square, Circle, Cube, Cylinder, Polygon
export PolyUnion
export difference
export ⋃, ⋂
export offset, hull, minkowski, convex_hull

end #<<<1 module
# >>>1

# macro use(m)#<<<
# 	N = filter(x->x != m, names(eval(m)))
# 	Expr(:block, :(using .$m),
# 		map(N) do x quote
# 			$(esc(x))(args...) = eval($(esc(m))).$x(args...)
# 		end
# 	end...)
# end#>>>

# vim: fdm=marker fmr=<<<,>>> noet ts=2:
