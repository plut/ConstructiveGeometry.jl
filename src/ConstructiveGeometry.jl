module ConstructiveGeometry
#————————————————————— Tools —————————————————————————————— ««1

#»»1
# Module imports««1
# using Printf
using LinearAlgebra
import LinearAlgebra: det
using StaticArrays
using FixedPointNumbers
using FastClosures
import Unitful: Unitful, °
import FileIO: FileIO, @format_str
# import MeshIO

import Rotations
import Colors: Colors, Colorant
import Makie: Makie, surface, plot, plot!, SceneLike, mesh

import Base: show, print
import Base: union, intersect, setdiff, symdiff, copy, isempty, merge
import Base: *, +, -, ∈, inv, sign, iszero

import AbstractTrees: AbstractTrees, children

# sub-packages:
include("ConvexHull.jl")
using .ConvexHull
include("Shapes.jl")
using .Shapes
include("TriangleMeshes.jl")
using .TriangleMeshes
include("Bezier.jl")
using .Bezier
include("Voronoi.jl")

# General tools««1
struct Consecutives{T,V} <: AbstractVector{T}
	parent::V
end
consecutives(v::AbstractVector{T}) where{T} =
	Consecutives{T,typeof(v)}(v)
Base.getindex(c::Consecutives, i::Integer) =
	(c.parent[i], c.parent[mod1(i+1, length(c.parent))])
Base.size(c::Consecutives) = size(c.parent)

# Geometry ««2
@inline one_half(x::Real) = x/2
@inline norm²(p) = sum(p .* p)
# round to a multiple of `m`:
@inline round(m::Integer, x::Integer, ::typeof(RoundUp)) = x + m - mod1(x, m)

@inline function unitnormal(a::StaticVector{3}, b::StaticVector{3})
	v = cross(a,b)
	return v/norm(v)
end

# Angle types««1
# All angles are internally stored as degrees: the user might (often)
# want exact degrees angles and probably never exact radians angles.
const _UNIT_DEGREE = typeof(°)
Angle{X,T} = Unitful.Quantity{T,Unitful.NoDims,X}
todegrees(x::Angle) = x
todegrees(x::Real) = x*°

# Affine transformations««1
struct AffineMap{A,B}
	a::A
	b::B
end

struct ZeroVector end # additive pendant of UniformScaling
+(a::ZeroVector, ::ZeroVector) = a
+(a, ::ZeroVector) = a
+(::ZeroVector, a) = a
-(z::ZeroVector) = z
-(a,::ZeroVector) = a
*(_, z::ZeroVector) = z
*(z::ZeroVector, _) = z
SVector{N}(::ZeroVector) where{N} = zero(SVector{N})

# identify size of matrices
@inline hassize(a::AbstractArray, s) = size(a) == s
@inline hassize(::UniformScaling, _) = true
@inline hassize(::ZeroVector, _) = true
@inline hassize(f::AffineMap, (m,n)) = hassize(f.a,(m,n)) && hassize(f.b,(m,))

@inline signdet(a::UniformScaling, d) = iszero(a.λ) ? 0 : iseven(d) ? +1 : a.λ
@inline signdet(a::AbstractMatrix, _) = det(a)

SAffineMap2{T} = AffineMap{SMatrix{2,2,T,4},SVector{2,T}}
SAffineMap3{T} = AffineMap{SMatrix{3,3,T,9},SVector{3,T}}
Base.convert(A::Type{SAffineMap3{T}}, u::UniformScaling) where{T} =
	A(SMatrix{3,3,T}(u), zero(SVector{3,T}))
SAffineMap3(f::AffineMap) = AffineMap(SMatrix{3,3}(f.a), SVector{3}(f.b))

@inline compose(f::AffineMap, g::AffineMap) = AffineMap(f.a*g.a, f.a*g.b+f.b)
@inline (f::AffineMap)(v::AbstractVector) = f.a*v + f.b
function (f::AffineMap)(s::Shapes.PolygonXor)
	d = signdet(f.a, 2)
	iszero(d)&& throw(SingularException(2))
	d > 0 && return Shapes.PolygonXor([ f.(p) for p in Shapes.paths(s) ])
	return Shapes.PolygonXor([reverse(f.(p)) for p in Shapes.paths(s) ])
end

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
# StaticArrays.similar_type(::Vector, T::Type) = Vector{T} # piracy!
@inline Reflection(v::AbstractVector) =
	Reflection{typeof(v[1]/1)}(v)
@inline Reflection{T}(v::AbstractVector) where{T} =
	Reflection{T, SVector{length(v),T}}(v, transpose(v*inv(dot(v,v))))

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

#————————————————————— Basic structure —————————————————————————————— ««1

# Abstract base type and interface««1
"""
    AbstractGeometry{D}

A `D`-dimensional geometric object.
Interface:
 - `children(x)`
 - `mesh(::MeshOptions, x)
 - `auxmeshes(::MeshOptions, x)` (optional)
 - `printnode(::IO, x)` (optional)
"""

abstract type AbstractGeometry{D} end
@inline children(t::AbstractGeometry) =
	error("children() not defined for type ", typeof(t))
abstract type AbstractGeometryLeaf{D} <: AbstractGeometry{D} end
@inline children(::AbstractGeometryLeaf) = AbstractGeometry[]
@inline AbstractTrees.printnode(io::IO, s::AbstractGeometry) =
	print(io, typeof(s).name.name)

include("scad.jl")
"""
    mesh(opt::MeshOptions, object, children_meshes)

This is the main method used by each concrete `AbstractGeometry` subtype
to compute the main mesh for this object type
from the (possibly empty) list of meshes of its children.

This method should not be called by the end user.
"""
function mesh end

"""
    AbstractMesh{D,T}

Base type for fully-explicit objects with dimension `D` and coordinate type `T`.
Derived types should implement the following methods:
 - `plot`
 - `points`

FIXME: this type is currently unused, should it be removed?
"""
abstract type AbstractMesh{D,T} <: AbstractGeometryLeaf{D} end
@inline Base.map!(f, m::AbstractMesh) = vertices(m) .= f.(vertices(m))
# Empty unions and intersects««2
"""
    EmptyUnion

A convenience type representing the union of nothing.
This is removed whenever it is `union()`-ed with anything else.

This is *not* a subtype of AbtractGeometry or anything,
since we overload `union()` (which is generally undefined) to this.
"""
struct EmptyUnion end
union() = EmptyUnion()
Base.show(io::IO, ::EmptyUnion) = print(io, "union()")
"""
    EmptyIntersect

A convenience type representing the intersection of nothing.
This is removed whenever it is `intersect()`-ed with anything else.

This is *not* a subtype of AbtractGeometry or anything,
since we overload `intersect()` (which is generally undefined) to this.
"""
struct EmptyIntersect end
intersect() = EmptyIntersect()
Base.show(io::IO, ::EmptyIntersect) = print(io, "intersect()")

# Meshing««1
# MeshOptions type««2
"""
    MeshOptions{T,C}

A structure storing the parameters for the (abstract object)->(mesh) conversion.
The coordinate type determines the Julia type of the returned object;
therefore, it is stored as a type parameter of `MeshOptions`
instead as data.

 - `T`: coordinate type of returned object
 - `C`: color type of returned object

NOTE: for now, only `Float64` is possible.
"""
struct MeshOptions{T<:Real,C,P<:NamedTuple}
	# at least `atol`, `rtol` and `symmetry`, but user data is allowed here:
	parameters::P
	color::C
	@inline (M::Type{<:MeshOptions{T}})(p::P,c::C) where{T,C,P} = new{T,C,P}(p,c)
end

@inline Base.merge(g::MeshOptions{T,C}; kwargs...) where{T,C} =
	MeshOptions{T,C}(merge(g.parameters, kwargs.data), g.color)
const MeshColor = Colors.RGBA{N0f8}
const _DEFAULT_COLOR = Colors.RGBA{N0f8}(.3,.4,.5) # bluish gray
const _DEFAULT_OPTIONS = MeshOptions{Float64}(NamedTuple(), _DEFAULT_COLOR)
const _DEFAULT_PARAMETERS = (
	atol = 0.1, rtol = .005, symmetry = 1, icosphere = 1000,
)

@inline MeshOptions(g::MeshOptions{T}, p::NamedTuple) where{T} =
	MeshOptions{T}(merge(g.parameters, p), g.color)
@inline Base.get(g::MeshOptions, name) =
	get(g.parameters, name, _DEFAULT_PARAMETERS[name])

# meshing functions ««2
# the mesh types (ShapeMesh, VolumeMesh) are defined below
# this is the infrastructure for recursive meshing:

struct FullMesh{T,U,A}
	main::T
	aux::Vector{Pair{A,U}}
end

# @inline FullMesh(main::T, aux::Vector{Pair{A,U}}) where{T,U,A} =
# 	FullMesh{T,U,A}(main, aux)
@inline auxtype(::Type{FullMesh{T,U,A}}) where{T,U,A} = Pair{A,U}

@inline fullmesh(s::AbstractGeometry; kwargs...) =
	mesh(merge(_DEFAULT_OPTIONS;kwargs...), s)

"""
    mesh(g::MeshOptions, object)

Recursively compute the main mesh of this object,
calling the `mesh` and `auxmeshes` functions as needed.
"""
function mesh(g::MeshOptions, s::AbstractGeometry)
	# we only ever overload Makie.mesh(::MeshOptions, ...), promise!
	l = [ mesh(g, x) for x in children(s) ]
	m = mesh(g, s, [ x.main for x in l ])
	a = auxmeshes(g, s, m, [ x.aux for x in l ])
	return MeshType(g,s)(m, a)
end
# special cases below: set_parameters

"""
    auxmeshes(opt::MeshOptions, object, children_mains, children_auxes)

Returns only auxiliary meshes for this object.
"""
@inline auxmeshes(g::MeshOptions, s::AbstractGeometry, m, l) =
	let T = auxtype(MeshType(g,s))
	eltype(eltype(l)) == T ? [l...; ] : T[]
end
# special cases below: highlight


# gridcells ««2
"""
    gridcells(g::MeshOptions, vertices, maxgrid)

Returns the number of cells subdividing `vertices`,
with a maximum of `maxgrid` (if it is nonzero).
"""

function gridcells(g::MeshOptions, vertices, maxgrid)
	bbox = [ extrema(v[i] for v in vertices) for i in SOneTo(3)]
	width= maximum(x[2]-x[1] for x in bbox)
	# theoretical number of points needed:
	N = ceil(Int, min(1/get(g, :rtol), width/get(g,:atol)))
	return maxgrid > 0 ? min(maxgrid, N) : N
end

#————————————————————— Primitive objects —————————————————————————————— ««1
# Explicit case: ShapeMesh««1
"""
    ShapeMesh{T}

The exclusive union of a number of simple, closed polygons.
"""
struct ShapeMesh{T} <: AbstractMesh{2,T}
	# encapsulates the following as an AbstractGeometry object:
	poly::Shapes.PolygonXor{T}
	# plus a 2d -> 3d transform:
	position::SAffineMap3{T}

	@inline ShapeMesh{T}(s::Shapes.PolygonXor, f=I) where{T}= new{T}(s,f)
	@inline ShapeMesh(s::Shapes.PolygonXor{T}, f=I) where{T}= new{T}(s,f)
end
@inline ShapeMesh(::MeshOptions{T}, paths, f=I) where{T} =
	ShapeMesh(Shapes.PolygonXor{T}(paths), f)
@inline poly(s::ShapeMesh) = s.poly
@inline paths(s::ShapeMesh) = poly(s).paths
@inline position(s::ShapeMesh) = s.position
@inline vertices(s::ShapeMesh) = Shapes.vertices(poly(s))

@inline AbstractTrees.printnode(io::IO, s::ShapeMesh) =
	print(io, "ShapeMesh # ", length(paths(s)), " polygon(s), ",
		sum(length.(paths(s))), " vertices")
@inline mesh(g::MeshOptions, m::ShapeMesh, _) = ShapeMesh(g, paths(m))

@inline MeshType(::MeshOptions{T,C},::AbstractGeometry{2}) where{T,C} =
	FullMesh{ShapeMesh{T},Shapes.PolygonXor{T},C}
@inline raw(m::ShapeMesh) = m.poly
@inline Base.empty(::ShapeMesh{T}) where{T} =
	ShapeMesh{T}(Shapes.PolygonXor{T}([]))

"""
    polygon(paths...; fill=:nonzero)

Filled polygon delimitated by the given vertices.

Crossing paths and reversed polygon are allowed; they will be simplified
upon polygon creation, using the fill method given as `fill`.
Possible methods are: `:nonzero`, `:evenodd`, `:positive`.
"""
function polygon(paths::AbstractVector...; fill=:nonzero)
	# we need a floating-point type here to be able to simplify crossings:
	T = float(promote_type((promote_type(Int,eltype.(p)...) for p in paths)...))
	ShapeMesh(Shapes.simplify(Shapes.PolygonXor{T}([paths...]); fill))
end
@inline polygon(s::AbstractGeometry{2}) = fullmesh(s).main
"""
    bezier(p0,p1,p2,p3,n)

Interpolates `n` points of the Bézier curve parametrized by points
`p0`..`p3`. The points are (parameter-wise) regularly spaced.
"""
@inline bezier(p0,p1,p2,p3,n::Integer) =
	Bezier.interpolate(Bezier.BezierCurve(p0,p1,p2,p3),n) # [points..., bezier()...]

# Square««1
struct Square{T} <: AbstractGeometryLeaf{2}
	size::SVector{2,T}
end
@inline scad_info(s::Square) = (:square, (size=s.size,))
@inline square_vertices(u, v) = [ SA[0,0], SA[u,0], SA[u,v], SA[0,v]]

@inline mesh(g::MeshOptions{T}, s::Square, _) where{T} =
	ShapeMesh(g, [square_vertices(T(s.size[1]), T(s.size[2]))])

# TODO: add syntactic sugar (origin, center, etc.) here:
"""
    square(size; origin, center=false)
    square(width, height; origin, center=false)

An axis-parallel square or rectangle  with given `size`
(scalar or vector of length 2).
"""
@inline square(a::Real, b::Real; center=nothing, anchor=nothing) =
	_square(a, b, center, anchor)
@inline square(a::Real; kwargs...) = square(a,a; kwargs...)
@inline square(a::AbstractVector; kwargs...) = square(a...; kwargs...)

@inline _square(a, b, ::Nothing, ::Nothing) = Square(SA[a, b])
@inline _square(a, b, center::Bool, anchor) =
	center ? _square(a,b,SA[0,0],anchor) : _square(a,b,nothing,anchor)
@inline _square(a, b, center::AbstractVector, anchor) =
	translate(center-SA[one_half(a),one_half(b)])*_square(a,b,nothing,anchor)

# Circle««1
struct Circle{T} <: AbstractGeometryLeaf{2}
	radius::T
	circumscribed::Bool
end

@inline scad_info(s::Circle) = (:circle, (r=s.radius,))

# Regular polygon««2
"""
    unit_n_gon(T::Type, n::Int)
		unit_n_gon(r, parameters::NamedTuple)

Returns the vertices of a regular n-gon inscribed in the unit circle
as points with coordinates of type `T`, while avoiding using too many
trigonometric computations.
"""
function unit_n_gon(r::T, n::Int) where{T<:Real}
	ω = cispi(2/n) # exp(2iπ/n)
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
@inline unit_n_gon(r::Rational{BigInt}, n::Int) =
	SVector{2,Rational{BigInt}}.(unit_n_gon(Float32(r), n))

# number of sides of a circle ««2

"""
    circle_nvertices(g::MeshOptions, radius)

Returns the number of sides used to draw a circle (arc) of given angle.
The base value `n` is given by the minimum of:
 - atol: each sagitta (s= r(1-cos 2π/n)) must not be smaller
 than `atol`, or n = π √(r/2 atol);
 - rtol: s/r = 1-cos(2π/n)  not smaller than rtol,
 or n = π /√(2*rtol).
"""
function circle_nvertices(g::MeshOptions, r)
	ε = max(get(g,:rtol), get(g,:atol)/r)
	base = ceil(Int, π/√(2ε))
	# a circle always has at least 4 sides
	return round(get(g,:symmetry), max(4, base), RoundUp)
end

function mesh(g::MeshOptions{T}, s::Circle, _) where{T}
	n = circle_nvertices(g, s.radius)
	r = T(s.circumscribed ? s.radius/cospi(1/n) : s.radius)
	return ShapeMesh(g, [unit_n_gon(r, n)])
end

# Interface««2
"""
    circle(r::Real, [circumscribed = false])

A circle with diameter `r`, centered at the origin.

The corresponding mesh is a regular polygon,
which is circumscribed to the ideal circle if `circumscribed == true`
and inscribed otherwise.
"""
@inline circle(a::Real; circumscribed=false) = Circle(a, circumscribed)

# Stroke ««1
struct Stroke{T} <: AbstractGeometryLeaf{2}
	points::Vector{SVector{2,T}}
	width::Float64
	ends::Symbol
	join::Symbol
	miter_limit::Float64
end
Stroke(points, width; ends=:round, join=:round, miter_limit=2.) =
	Stroke{Float64}(points, width, ends, join, miter_limit)
@inline AbstractTrees.printnode(io::IO, s::Stroke) =
	print(io, "Stroke # ", length(s.points), " points")

function mesh(g::MeshOptions{T}, s::Stroke, _) where{T}
	r = one_half(T(s.width))
	ε = max(get(g,:atol), get(g,:rtol) * r)
	return ShapeMesh(g, Shapes.offset([s.points], r;
		join=s.join, ends=s.ends, miter_limit = s.miter_limit))
end

"""
    stroke(points, width; kwargs)
    ends = :loop|:butt|:square|:round
    join = :round|:square|:miter
    miter_limit = 2.0

Draws a path of given width.
"""
stroke(points, width; kwargs...) = Stroke(points, width; kwargs...)


# VolumeMesh««1
"""
    VolumeMesh{T}

A surface delimited by the given faces.
"""
struct VolumeMesh{T,A} <: AbstractMesh{3,T}
	mesh::TriangleMesh{T,A}
	@inline VolumeMesh(m::TriangleMesh{T,A}) where{T,A} = new{T,A}(m)
end
@inline vertices(s::VolumeMesh) = TriangleMeshes.vertices(s.mesh)
@inline faces(s::VolumeMesh) = TriangleMeshes.faces(s.mesh)
@inline attributes(s::VolumeMesh) = TriangleMeshes.attributes(s.mesh)
@inline Base.empty(::VolumeMesh{T,A}) where{T,A} =
	VolumeMesh(TriangleMesh{T,A}([],[],[]))
@inline AbstractTrees.printnode(io::IO, s::VolumeMesh) =
	print(io, "VolumeMesh # ", length(vertices(s)), " vertices, ",
		length(faces(s)), " faces")

function apply(f::AffineMap, m::TriangleMesh{T,A}, d=signdet(f.a, 3)) where{T,A}
	iszero(d) && throw(SingularException(3))
	return TriangleMesh{T,A}(f.(m.vertices),
		d > 0 ? m.faces : reverse.(m.faces), m.attributes)
end
@inline apply(f::AffineMap, m::VolumeMesh, d=signdet(f.a, 3)) =
	VolumeMesh(apply(f, m.mesh, d))

@inline TriangleMesh(g::MeshOptions{T}, points, faces,
	attrs::AbstractVector{A} = fill(g.color, size(faces))) where{F,T,A} =
	TriangleMesh{T,A}(SVector{3,T}.(points), faces, attrs)
@inline VolumeMesh(g::MeshOptions{T,A}, points, faces) where{T,A} =
	VolumeMesh(TriangleMesh{T,A}(SVector{3,T}.(points), faces,
		fill(g.color, size(faces))))
@inline mesh(g::MeshOptions, s::VolumeMesh, _) =
	VolumeMesh(g, vertices(s), faces(s))

@inline MeshType(::MeshOptions{T,A},::AbstractGeometry{3}) where{T,A} =
	FullMesh{VolumeMesh{T,A},TriangleMesh{T,Nothing},A}
@inline raw(m::VolumeMesh) = TriangleMesh{Float64,Nothing}(
	m.mesh.vertices, m.mesh.faces, fill(nothing, size(m.mesh.faces)))

@inline scad_info(s::VolumeMesh) =
	(:surface, (points=s.points, faces = [ f .- 1 for f in s.faces ]))

"""
    surface(vertices, faces)

Produces a surface with the given vertices.
`faces` is a list of n-uples of indices into `vertices`.

Non-triangular faces are triangulated
(by being first projected on the least-square fit plane).
"""
function surface(vertices, faces)
	triangles = sizehint!(NTuple{3,Int}[], length(faces))
	for f in faces
		length(f) == 3 && (push!(triangles, (f[1], f[2], f[3])); continue)
		# to triangulate, we project onto the least-squares fit plane:
		a = float.(sum(vertices[i]*vertices[i]' for i in f))
		# The first eigenvector corresponds to the smallest eigenvalue:
		# (and has approximately unit norm)
		w = SVector{3}(eigvecs(a)[:,1])
		v = abs(w[1]) < abs(w[2]) ?
			(abs(w[1]) < abs(w[3]) ? SA[0,w[3],-w[2]] : SA[w[2],-w[1],0]) :
			(abs(w[2]) < abs(w[3]) ? SA[w[3],0,-w[1]] : SA[w[2],-w[1],0])
		u = cross(v, w)
		v = cross(w, u)
		# orientation-preserving triangulation:
		tri = Shapes.triangulate([
			SA[dot(u,vertices[i]), dot(v,vertices[i])] for i in f ])
		push!(triangles, tri...)
	end
	return VolumeMesh(_DEFAULT_OPTIONS, vertices, triangles)
end
@inline surface(s::AbstractGeometry{3}) = fullmesh(s).main

# Cube««1
struct Cube{T} <: AbstractGeometryLeaf{3}
	size::SVector{3,T}
end
@inline scad_info(s::Cube) = (:cube, (size=s.size,))

@inline cube_vertices(u, v, w) = [
		SA[0,0,0], SA[0,0,w], SA[0,v,0], SA[0,v,w],
		SA[u,0,0], SA[u,0,w], SA[u,v,0], SA[u,v,w]]

mesh(g::MeshOptions{T}, s::Cube, _) where{T} =
	VolumeMesh(TriangleMesh(g,
	cube_vertices(T(s.size[1]), T(s.size[2]), T(s.size[3])),
	[ # 12 triangular faces:
	 (6, 5, 7), (7, 8, 6), (7, 3, 4), (4, 8, 7),
	 (4, 2, 6), (6, 8, 4), (5, 1, 3), (3, 7, 5),
	 (2, 1, 5), (5, 6, 2), (3, 1, 2), (2, 4, 3),
	]))

"""
    cube(size; origin, center=false)
    cube(size_x, size_y, size_z; origin, center=false)

An axis-parallel cube (or sett) with given `size`
(scalar or vector of length 3).

The first vertex is at the origin and all vertices have positive coordinates.
If `center` is `true` then the cube is centered.
"""
@inline cube(a::Real, b::Real, c::Real; center=nothing, anchor=nothing) =
	_cube(a, b, c, center, anchor)

@inline _cube(a,b,c, ::Nothing, ::Nothing) = Cube(SA[a,b,c])
@inline _cube(a,b,c, center::Bool, anchor) =
	center ? _cube(a,b,c,SA[0,0,0],anchor) : _cube(a,b,c,nothing,anchor)
@inline _cube(a,b,c, center::AbstractVector, anchor) =
	translate(center-SA[one_half(a),one_half(b),one_half(c)])*
	_cube(a,b,c,nothing,anchor)

@inline cube(a::Real; kwargs...) = cube(a,a,a; kwargs...)
@inline cube(a::AbstractVector; kwargs...) = cube(a...; kwargs...)


# Sphere««1
struct Sphere{T} <: AbstractGeometryLeaf{3}
	radius::T
	circumscribed::Bool
end
@inline scad_info(s::Sphere) = (:sphere, (r=s.radius,))

# Fibonacci sphere««2
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
		ρ = T(√(1-z^2))
		(s,c) = sincos(θ)
		@inbounds v[i] = SVector{3,T}(r*c*ρ, r*s*ρ, r*z)
	end
	return v
end
@inline fibonacci_sphere_points(r::Rational{BigInt}, n::Int) =
	SVector{3,Rational{BigInt}}.(fibonacci_sphere_points(Float32(r), n))
# number of vertices of a sphere ««2
"""
    sphere_nvertices(g::MeshOptions, r::Real)

Returns the number `n` of points on a sphere according to these
parameters.

"""
function sphere_nvertices(g::MeshOptions, r)
	ε = max(get(g,:rtol), get(g,:atol)/r)
	base = 2 + ceil(Int, (π/√3)/ε)
	# a sphere always has at least 6 vertices
	return max(6, base)
end

"""
Returns (points, triangles) for an icosphere of radius `r`:
an icosahedron with each face subdivided into `n`^2 small triangles.
"""
function icosphere_points(r::T, n) where{T<:Real}#««
	# to each of the 20 faces we map the numbers for the three edges + vertices
	# (the edge label is negative when edge is reversed)
	elabel = (
	  (1,-5,21),(2,-1,22),(3,-2,23),(4,-3,24),(5,-4,25),
		(6,-21,-15),(7,-22,-11),(8,-23,-12),(9,-24,-13),(10,-25,-14),
		(11,-6,-26),(12,-7,-27),(13,-8,-28),(14,-9,-29),(15,-10,-30),
		(16,26,-20),(17,27,-16),(18,28,-17),(19,29,-18),(20,30,-19))
	# (lhs, rhs) vertices of each edge
	edges = (
		(11,1),(11,2),(11,3),(11,4),(11,5), (1,6),(2,7),(3,8),(4,9),(5,10),
		(1,7),(2,8),(3,9),(4,10),(5,6), (7,12),(8,12),(9,12),(10,12),(6,12),
		(1,5),(2,1),(3,2),(4,3),(5,4), (6,7),(7,8),(8,9),(9,10),(10,6))
	z, ρ = sqrt(.2), sqrt(.8)
	vertices = resize!([
		SA[ρ*T(cosd(0)), ρ*T(sind(0)), z],
		SA[ρ*T(cosd(72)), ρ*T(sind(72)), z],
		SA[ρ*T(cosd(144)), ρ*T(sind(144)), z],
		SA[ρ*T(cosd(216)), ρ*T(sind(216)), z],
		SA[ρ*T(cosd(288)), ρ*T(sind(288)), z],
		SA[ρ*T(cosd(324)), ρ*T(sind(324)), -z],
		SA[ρ*T(cosd(36)), ρ*T(sind(36)), -z],
		SA[ρ*T(cosd(108)), ρ*T(sind(108)), -z],
		SA[ρ*T(cosd(180)), ρ*T(sind(180)), -z],
		SA[ρ*T(cosd(252)), ρ*T(sind(252)), -z],
		SA[0,0,1.], SA[0,0,-1.]], 10*n^2+2);
	triangles = sizehint!(NTuple{3,Int}[], 20*n^2)
	# labeling of icosphere vertices:
	# 1:12 = (12) vertices of the icosahedron
	# 13:30n-18 = 30(n-1) on edges
	#    i.e. edge e is 13+(n-1)(e-1):12+(n-1)e
	# 30n-19:10n^2+2 = 20(n-1)(n-2)/2 inside faces  ((n-1)(n-2)/2 is integer)
	#    i.e. face f has F(f):F(f+1)-1, F(f) = 30*n-19+(n-1)(n-2)/2*(f-1)
	for e in 1:30, k in 1:n-1 # interpolate edge points
		v1, v2 = edges[e]
		vertices[12+(n-1)*(e-1)+k] = ((n-k)*vertices[v1] + k*vertices[v2])
	end
	for f in 1:20
		# we allocate an array for all the (n+1)(n+2)/2 points inside this face
		# (since we return a Θ(n^2) array eventually this does not add
		# complexity)
		idx = Vector{Int}(undef, (n+1)*(n+2)÷2)
		e1, e2, e3 = elabel[f]
		# \ side and top corner: vertices 1, k(k+1)÷2 ««
		if e1 > 0
			idx[1] = edges[e1][1]
			for k in 1:n-1
				idx[(k+1)*(k+2)÷2] = 12+(n-1)*(e1-1)+k
			end
		else
			e1 = -e1
			idx[1] = edges[e1][2]
			for k in 1:n-1
				idx[(k+1)*(k+2)÷2] = 12+(n-1)*(e1-1)+(n-k)
			end
		end#»»
		# / side and bottom-left corner: vertices (n²+n+2)/2, k^2-k+2 ««
		if e2 > 0
			idx[(n^2+n+2)÷2] = edges[e2][1]
			for k in 1:n-1
				idx[(n-k)*(n-k+1)÷2+1] = 12+(n-1)*(e2-1)+k
			end
		else
			e2 = -e2
			idx[(n^2+n+2)÷2] = edges[e2][2]
			for k in 1:n-1
				idx[(n-k)*(n-k+1)÷2+1] = 12+(n-1)*(e2-1)+(n-k)
			end
		end#»»
		# _ side and bottom-right corner: vertices (n+1)(n+2)÷2, (n²+n+2)÷2+k««
		if e3 > 0
			idx[(n+1)*(n+2)÷2] = edges[e3][1]
			for k in 1:n-1
				idx[(n+1)*(n+2)÷2 - k] = 12+(n-1)*(e3-1)+k
			end
		else
			e3 = -e3
			idx[(n+1)*(n+2)÷2] = edges[e3][2]
			for k in 1:n-1
				idx[(n+1)*(n+2)÷2 - k] = 12+(n-1)*(e3-1)+(n-k)
			end
		end#»»
		# central vertices««
		f0 = 30*n-18 + (n-1)*(n-2)*(f-1)÷2 # offset for this face
		v1, v2, v3 = vertices[idx[SA[1, (n^2+n+2)÷2, (n+1)*(n+2)÷2]]]
		for k in 1:n-2
			# build row (k²+3k+6)÷2 : (k²+5k+4)÷2
			idx[(k^2+3*k+6)÷2:(k^2+5*k+4)÷2] .= f0+(k^2-k+2)÷2 : f0+k*(k+1)÷2
			for i in 1:k
				vertices[f0+k*(k-1)÷2+i] = (n-k-1)*v1 + (k-i+1)*v2 + i*v3
			end
		end#»»
		# triangulate this face««
		v = 1
		for k in 1:n
			# write triangles between rows k and k+1
			# row k starts at v=(k^2-k+2)÷2
			# row k+1 starts at w=(k^2+k+2)÷2
			w = (k^2+k+2)÷2
			# △:  (w+i, w+i+1, v+i) i=0:k-1
			# ▽: (w+i+1, v+i+1, v+i) i=0:k-2
			for i in 0:k-2
				push!(triangles, (idx[w+i], idx[w+i+1], idx[v+i]),
					(idx[w+1+i], idx[v+1+i], idx[v+i]))
			end
			push!(triangles, (idx[w+k-1], idx[w+k], idx[v+k-1]))
			v=w
		end#»»
	end
	for i in 13:length(vertices)
		v = vertices[i]; vertices[i] = v * r/ √(v'*v)
	end
	return (vertices, triangles)
end#»»

function mesh(g::MeshOptions{T}, s::Sphere, _) where{T}
	n = sphere_nvertices(g, s.radius)
	r = s.circumscribed ? s.radius/sqrt(1-8π/√3/(n-2)) : s.radius
	if n < get(g, :icosphere) # Fibonacci case
		plist = fibonacci_sphere_points(T(r), n)
		(pts, faces) = convex_hull(plist)
		return VolumeMesh(g, pts, faces)
	else # Icosphere case
		(pts, faces) = icosphere_points(T(r), ceil(Int, sqrt(n ÷ 10)))
		return VolumeMesh(g, pts, faces)
	end
end

# Interface ««2
"""
    sphere(r::Real, [circumscribed=false])

A sphere with diameter `r`, centered at the origin.
"""
@inline sphere(a::Real; circumscribed=false) = Sphere(a, circumscribed)
#————————————————————— Transformations —————————————————————————————— ««1

# `AbstractTransform`:
# subtype of `AbstractGeometry` for objects with a single child
abstract type AbstractTransform{D} <: AbstractGeometry{D} end
@inline children(m::AbstractTransform) = (m.child,)

# Transform type««2
# An operator F may be called in either full or currified form:
# F(parameters, object) - applied
# F(parameters)*object - abstract transform (allows chaining)
"""
    Transform{S,F}

This type defines functions allowing to chain transforms; these are used by
`multmatrix`, `color` etc. operations (see below).
"""
struct Transform{S,F}
	f::F
	@inline Transform{S}(f) where{S} = new{S,typeof(f)}(f)
end
# end-case: operator applied on a single geometry object
@inline operator(f, x, s::AbstractGeometry; kw...) = f(x..., s; kw...)

# other cases reducing to this one:
@inline operator(f, x, s::Union{<:AbstractGeometry,EmptyUnion}...;kw...) =
	operator(f, x, union(s...);kw...)
@inline operator(f, x, s::AbstractVector; kw...) = operator(f, x, s...; kw...)
@inline operator(f, x; kw...) =
	Transform{f}(@closure s -> operator(f, x, s; kw...))

# Multiplicative notation:
@inline Base.:*(u::Transform, s) = u.f(s)
@inline (u::Transform)(s) = u*s
@inline (::Transform)(s::Union{EmptyUnion,EmptyIntersect}) = s
@inline Base.:*(u::Transform, v::Transform) = compose(u, v)
@inline compose(u::Transform, v::Transform) = Transform{typeof(∘)}(u.f∘v.f)

# Affine transforms««1
# AffineTransform type««2
struct AffineTransform{D,A,B} <: AbstractTransform{D}
	f::AffineMap{A,B}
	child::AbstractGeometry{D}
end
# AffineTransform(f, child) constructor is defined
@inline mesh(g::MeshOptions, s::AffineTransform{3}, (m,)) = apply(s.f, m)

function mesh(g::MeshOptions, s::AffineTransform{2}, (m,))
	if hassize(s.f, (2,2))
		return ShapeMesh(s.f(poly(m)), m.position)
	elseif hassize(s.f,(3,3))
		return ShapeMesh(poly(m), compose(SAffineMap3(s.f), position(m)))
	elseif hassize(s.f,(3,2)) # pad by one column:
		f = SAffineMap3(AffineMap([s.f.a SA[0;0;1]], s.f.b))
		return ShapeMesh(poly(m), compose(f, position(m)))
	end
	throw(DimensionMismatch("linear map * shape should have dimension (2,2) or (3,3)"))
end

function auxmeshes(g::MeshOptions, s::AffineTransform{3}, m, l)
	d = signdet(s.f.a, 3) # precompute the determinant
	return [ x => apply(s.f, y, d) for (x,y) in [l...;] ]
end

# affine_transform««2
affine_transform(f, s::AbstractGeometry) = AffineTransform(f, s)
# This composes matrices before applying them to the objects,
# thus saving time (and hierarchy complexity).
affine_transform(f, s::AffineTransform) =
	AffineTransform(compose(f, s.f), s.child)

# mult_matrix««2
"""
    mult_matrix(a, [center=c], solid...)
    mult_matrix(a, b, solid...)
    mult_matrix(a, b) * solid
    a * solid + b # preferred form

Represents the affine operation `x -> a*x + b`.

# Extended help
!!! note "Types of `mult_matrix` parameters"

    The precise type of parameters `a` and `b` is not specified.
    Usually, `a` will be a matrix and `b` a vector, but this is left open
    on purpose; for instance, `a` can be a scalar (for a scaling).
    Any types so that `a * Vector + b` is defined will be accepted.

    Conversion to a matrix will be done when meshing.

!!! note "Matrix multiplication"

    Chained affine transformations are composed before applying
    to the objects. This saves time: multiple
    (3 × n) matrix multiplications are replaced by
    (3 × 3) multiplications, followed by a single (3 × n).
"""
@inline mult_matrix(a, s...; center=ZeroVector()) =
	operator(affine_transform, (AffineMap(a, center-a*center),), s...)
# translate ««2
"""
    translate(v, s...)
    translate(v) * s
    v + s

Translates object(s) `s...` by vector `v`.
"""
@inline translate(v,s...)= operator(affine_transform,(AffineMap(I,v),), s...)
"""
    raise(z, s...)

For volumes: equivalent to `translate([0,0,z], s...)`.
For shapes: equivalent to `translate([0,z], s...)`.
"""
@inline raise(z, s...) = operator(_raise, (z,), s...)
@inline _raise(z, s::AbstractGeometry{2}) = translate(SA[0,z], s)
@inline _raise(z, s::AbstractGeometry{3}) = translate(SA[0,0,z], s)
"""
    lower(z, s...)

For volumes: equivalent to `translate([0,0,-z], s...)`.
For shapes: equivalent to `translate([0,z], s...)`.
"""
@inline lower(z, x...) = raise(-z, x...)

# scale««2
# `scale` is really the same as one of the forms of `mult_matrix`:
"""
    scale(a, s...; center=0)
    scale(a; center=0) * s
    a * s
Scales solids `s` by factor `a`. If `center` is given then this will be
the invariant point.

`a` may also be a vector, in which case coordinates will be multiplied by
the associated diagonal matrix.
"""
@inline scale(a,s...; kwargs...) = mult_matrix(scaling(a), s...; kwargs...)
@inline scaling(a::Real) = a*I
@inline scaling(a::AbstractVector) = Diagonal(a)

# reflect««2
"""
    reflect(v, s...; center=0)
    reflect(v; center=0) * s

Reflection with axis given by the hyperplane normal to `v`.
If `center` is given, then the affine hyperplane through this point will
be used.
"""
@inline reflect(v::AbstractVector, s...; kwargs...) =
	mult_matrix(Reflection(v), s...; kwargs...)
# rotate««2
"""
    rotate(θ, [center=center], [solid...])

Rotation around the Z-axis (in trigonometric direction, i.e.
counter-clockwise).

    rotate(θ, axis=axis, [center=center], [solid...])

Rotation around any axis, given by a vector directing this axis.
"""
@inline rotate(θ::Number, s...; axis=nothing, kwargs...) =
	operator(_rotate1, (todegrees(θ), axis,), s...; kwargs...)
# we dispatch on the dimension of `s` to make (Angle2d) act on z-axis:
@inline _rotate1(θ, ::Nothing, s::AbstractGeometry{2}; kwargs...) =
	mult_matrix(Rotations.Angle2d(float(θ)), s; kwargs...)
@inline _rotate1(θ, axis::AbstractVector, s::AbstractGeometry{3}; kwargs...) =
	mult_matrix(Rotations.AngleAxis(θ, axis...), s; kwargs...)
@inline _rotate1(θ, ::Nothing, s::AbstractGeometry{3}; kwargs...) =
	mult_matrix(Rotations.RotZ(float(θ)), s; kwargs...)
"""
    rotate((θ,φ,ψ), [center=center], [solid...])

Rotation given by Euler angles (ZYX; same ordering as OpenSCAD).
"""
@inline rotate(angles, s...; kwargs...) =
	mult_matrix(Rotations.RotZYX(todegrees.(angles)...), s...; kwargs...)


# Projection««1
struct Project <: AbstractTransform{2}
	child::AbstractGeometry{3}
end

function mesh(g::MeshOptions, ::Project, (m,))
	triangles = TriangleMeshes.project(m.mesh)
	# FIXME: strangely, this makes some tiny “holes”
	# (projection should be exact?!). Could be Clipper's fault?
	c = (Shapes.clip(:union,
		[Shapes.PolygonXor([t]) for t in triangles]...))
	filter!(p->abs(Shapes.area(p)) > 1e-9, c.paths)
	return ShapeMesh(c)
end

"""
    project(s...)

Computes the (3d to 2d) projection of a shape on the horizontal plane.
"""
@inline project(s...) = operator(Project, (), s...)

# Slicing ««1
struct Slice{T} <: AbstractTransform{2}
	z::T
	child::AbstractGeometry{3}
end

function mesh(g::MeshOptions, s::Slice, (m,))
	(pts, seg) = TriangleMeshes.plane_slice(s.z, m.mesh)
	return ShapeMesh(Shapes.glue_segments(pts, seg))
end

"""
    slice(z, s...)
    slice(s...)

Computes the (3d to 2d) intersection of a shape and the given horizontal plane
(at `z=0` if not precised).
"""
@inline slice(z::Real, s...) = operator(Slice, (z,), s...)
@inline slice(s...) = operator(Slice, (0,), s...)

# TODO: slice at different z value


# Half-plane and half-space intersection ««1
"""
    Half{D,T}

Intersection with half-space `direction`*x == `origin`.
"""
struct Half{D,T}<: AbstractTransform{D}
	direction::SVector{D,T}
	origin::T
	child::AbstractGeometry{D}
end

symbolic_dir = ( (), # 1 is ignored
	(right = SA[1,0], left = SA[-1,0], top = SA[0,1], bottom=SA[0,-1],
	 front = SA[0,-1], back=SA[0,1],),
	(right = SA[1,0,0], left=SA[-1,0,0], top=SA[0,0,1], bottom=SA[0,0,-1],
	 front = SA[0,-1,0], back=SA[0,1,0],),
)
@inline Half(dir::Symbol, origin, s::AbstractGeometry{D}) where{D} =
	Half(symbolic_dir[D][dir], origin, s)
@inline Half(dir::AbstractVector{T}, ::ZeroVector, s) where{T} =
	Half(dir, zero(T), s)
@inline Half(dir::AbstractVector{T}, origin::AbstractVector, s) where{T} =
	Half(dir, dot(dir, origin), s)

@inline origin_vec(dir, off) = dir*(off / norm²(dir))

mesh(g::MeshOptions, s::Half{3}, (m,)) =
	VolumeMesh(TriangleMeshes.halfspace(
	-s.direction, origin_vec(s.direction, s.origin), m.mesh,g.color))

mesh(g::MeshOptions, s::Half{2}, (m,)) =
	ShapeMesh(intersect(Shapes.HalfPlane(s.direction,-s.origin), m.poly))

"""
    half(direction, s...; origin = 0)
    half(direction; origin = 0) * s

Keeps only the part of objects `s` lying in the halfspace/halfplane
with given `direction` and `origin`.

`direction` may be either a vector, or one of the six symbols  `:top`, `:bottom`, `:left`, `:right`, `:front`, `:back`.

`origin` may be either a point (i.e. one point on the hyperplane)
or a scalar (b in the equation a*x=b of the hyperplane).
"""
@inline half(dir, s...; origin=ZeroVector()) =
	operator(Half, (dir, origin,), s...)

# Prisms and cones ««1
# Prism meshes ««2
function tube_mesh(loops, slices)#««
	# this computes the combinatorial mesh of an open tube,
	# given the loops of the base and the number of slices (≥ 1)
	layer = sum(length.(loops)) # total number of points in a layer
	faces = sizehint!(NTuple{3,Int}[], layer*2*slices)
	p = 0
	for s in 1:slices
		q = p + layer
		for l in loops, (i,j) in consecutives(l)
			push!(faces, (i+p, j+p, i+q), (j+p,j+q,i+q))
		end
		p = q
	end
	return faces
end#»»
function prism_mesh(g::MeshOptions, h0, h1, twist, scale, m::ShapeMesh)#««
	pts2 = Shapes.vertices(m.poly)
	tri = Shapes.triangulate(m.poly)
	peri = Shapes.loops(m.poly)
	# perimeters are oriented ↺, holes ↻
	n = length(pts2)

	rmax = sqrt(maximum(norm².(pts2)))
	nslices = max(Int(cld(circle_nvertices(g, rmax)*twist, 360°)), 1)
	pts3 = sizehint!([SA[p..., 0] for p in pts2], nslices*n)
	dz = (h1-h0) / nslices
	dq = cis(twist/nslices)
	dr = (scale .- 1)/nslices
	z = h0; q = 1; r = SA[1.,1.]
	N = nslices - iszero(scale)
	for _ in 1:N
		z+= dz; q*= dq; r+= dr
		push!(pts3, (SA[r[1]*(p[1]*real(q)-p[2]*imag(q)),
			r[2]*(p[1]*imag(q)+p[2]*real(q)), z] for p in pts2)...)
	end
	faces = [reverse.(tri); tube_mesh(peri, N); ]
	if iszero(scale)
		# end by a cone
		push!(pts3, [0,0,h1])
		p = n*N
		faces = [faces;
			[(i+p,j+p,n+1+p) for l in peri for (i,j) in consecutives(l)]]
	else
		# close top face
		faces = [faces; [f .+ n*N for f in tri ]; ]
	end
	# reposition according to 3d position of child mesh
	a = let a=m.position.a, n=unitnormal(a[:,1], a[:,2])
		SA[a[1,1] a[1,2] n[1]; a[2,1] a[2,2] n[2]; a[3,1] a[3,2] n[3]]
	end
	return apply(AffineMap(a, m.position.b), VolumeMesh(g, pts3, faces))
end#»»
# Prism ««2
struct Prism{T} <: AbstractTransform{3}
	height::T
	twist::Angle
	scale::SVector{2,T}
	child::AbstractGeometry{2}
	@inline Prism(height::T, twist, scale::AbstractVector{S}, child) where{T,S} =
		new{promote_type(T,S)}(height, todegrees(twist), scale, child)
end

@inline mesh(g::MeshOptions, s::Prism, (m,)) =
	prism_mesh(g, 0, s.height, s.twist, s.scale, m)

"""
    prism(h, s...; twist, scale) # preferred form
    prism(h; twist, scale) * s
    linear_extrude(h, s...; twist, scale)
    linear_extrude(h) * s...

Build a prism (linear extrusion) of height `h` on the given base shape.
"""
@inline prism(height, s...; center=false, twist=0, scale=1) =
	operator(_prism, (height, center, twist, scale), s...)
@inline _prism(h,c,t,s::Real,x) = _prism(h,c,t,SA[s,s],x)
@inline function _prism(height, center, twist, scale::AbstractVector,
		s::AbstractGeometry{2})
	m = Prism(height, twist, scale, s)
	center ? translate([0,0,-one_half(height)], m) : m
end
linear_extrude = prism

# Cylinder ««2
#     cylinder(h, (r1, r2) [, center=false])
"""
    cylinder(h, r , [center=false], [circumscribed=false])
    cylinder(h, r1, r2 , [center=false], [circumscribed=false])

A cylinder (or cone frustum)
with basis centered at the origin, lower radius `r1`, upper radius `r2`,
and height `h`.

The mesh is a regular prism,
circumscribed to the cylinder if `circumscribed == true`
and inscribed otherwise.

!!!    note "OpenSCAD compatibility"

    `cylinder(h,r)` is interpreted as `cylinder(h,r,r)`,
    not `(h,r,0)` as in OpenSCAD. To input a cone, use `cone(h,r)` instead.

"""
@inline cylinder(h::Real, r::Real;center::Bool=false, circumscribed=false) =
	let c = prism(h)*circle(r;circumscribed)
	center ? lower(one_half(h))*c : c
	end

function cylinder(h::Real, r1::Real, r2::Real; center::Bool = false, circumscribed=false)
	a1, a2 = minmax(r1, r2)
	m = prism(h, scale=a1/a2)*circle(a2;circumscribed)
	b = 2(r1<r2)
	return SA[0,0,one_half(h)*(b-center)] + SA[1 0 0;0 1 0;0 0 1-b]*m
end

# Cone ««2
struct Cone{T} <: AbstractTransform{3}
	apex::SVector{3,T}
	child::AbstractGeometry{2}
	@inline Cone(a::AbstractVector{T}, c::AbstractGeometry{2}) where{T} =
		new{T}(a,c)
end

function mesh(g::MeshOptions{T}, s::Cone, (m,)) where{T}
	m1 = prism_mesh(g, 0, 1, 0°, SA[0,0], m)
	vertices(m1)[end] = s.apex
	return m1
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
    cone(h, r; circumscribed=false)

Circular right cone with basis centered at the origin,
radius `r`, and height `h`. Equivalent to `cone([0,0,h])*circle(r)`.

    cone(apex, r; circumscribed=false)

Circular, possibly oblique, cone with given apex point
and radius `r` around the origin.
"""
@inline cone(v::AbstractVector, r::Real; circumscribed=false) =
	Cone(v, circle(r; circumscribed))
@inline cone(h::Real, r::Real; circumscribed=false) =
	cone(SA[0,0,h], r; circumscribed)

# Rotate extrusion««1
struct Revolution{T,X} <: AbstractTransform{3}
	angle::Angle{_UNIT_DEGREE,T}
	slide::X
	child::AbstractGeometry{2}
end

# ladder triangles ««2
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

# _revolution««2
"""
    _revolution(point, data, parameters)

Extrudes a single point, returning a vector of 3d points
(x,y) ↦ (x cosθ, x sinθ, y).
"""
function _revolution(g::MeshOptions, p::StaticVector{2,T}, angle, slide) where{T}
	@assert p[1] ≥ 0
	# special case: point is on the y-axis; returns a single point:
	iszero(p[1]) && iszero(slide) && return [SA[p[1], p[1], p[2]]]
	n = Int(cld(circle_nvertices(g, p[1])*angle, 360°))

	ω = Complex{T}(cosd(angle/n), sind(angle/n))
	z = Vector{Complex{T}}(undef, n+1)
	z[1] = one(T)
	for i in 2:n
		@inbounds z[i] = z[i-1]*ω; z[i]/= abs(z[i])
	end
	# close the loop:
	z[n+1] = Complex{T}(cosd(angle), sind(angle))
	s = slide / n
	return [SA[p[1]*real(u), p[1]*imag(u), p[2]+(i-1)*s] for (i, u) in pairs(z)]
end

# Meshing ««2
function mesh(g::MeshOptions{T}, s::Revolution, (m,)) where{T}
	# right half of child:
	m1 = intersect(Shapes.HalfPlane(SA[1,0],0), m.poly)
	pts2 = Shapes.vertices(m1)
	tri = Shapes.triangulate(m1)
	peri = Shapes.loops(m1) # oriented ↺
	n = length(pts2)
	
	pts3 = _revolution(g, pts2[1], s.angle, s.slide)
	firstindex = [1]
	arclength = Int[length(pts3)]
	@debug "newpoints[$(pts2[1])] = $(length(pts3))"
	for p in pts2[2:end]
		push!(firstindex, length(pts3)+1)
		newpoints = _revolution(g, p, s.angle, s.slide)
		@debug "newpoints[$p] = $(length(newpoints))"
		push!(arclength, length(newpoints))
		pts3 = vcat(pts3, newpoints)
	end
	@debug "pts3: $pts3"
	@debug "firstindex: $firstindex"
	@debug "arclength: $arclength"
	# point i ∈ polygonxor: firstindex[i]:(firstindex[i]-1+arclength[i])
	get3 = @closure (a,t) -> (a[t[1]], a[t[2]], a[t[3]])
	triangles = vcat(
		[ get3(firstindex,t) for t in tri ],
		[ get3(firstindex,t) .+ get3(arclength,t) .- 1 for t in reverse.(tri) ]
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
			nt = ladder_triangles(
				arclength[p2], arclength[p1],
				firstindex[p2], firstindex[p1],
				)
			push!(triangles, nt...)
			@debug "new triangles: $(nt)"
		end
		@debug "(end perimeter)»»"
	end
	@debug "triangles = $triangles"
	return VolumeMesh(g, pts3, triangles)
end

# Interface ««2
"""
    revolution([angle = 360°], shape...; slide=0)
    revolution([angle = 360°]; [slide=0]) * shape
    rotate_extrude([angle = 360°], shape...; slide=0)
    rotate_extrude([angle = 360°]; [slide=0]) * shape

Similar to OpenSCAD's `rotatre_extrude` primitive.

The `slide` parameter is a displacement along the `z` direction.
"""
@inline revolution(angle::Number, s...; slide=0) =
	operator(Revolution, (todegrees(angle), slide,), s...)
@inline revolution(s...; kwargs...) = revolution(360°, s...; kwargs...)
rotate_extrude = revolution


# Sweep (surface and volume)««1
# ShapeSweep ««2
struct ShapeSweep <: AbstractGeometry{3}
	trajectory::AbstractGeometry{2}
	profile::AbstractGeometry{2}
# 	closed::Bool
	join::Symbol
	miter_limit::Float64
	@inline ShapeSweep(path, child; join=:round, miter_limit=2.0) =
		new(path, child, join, miter_limit)
end
@inline children(s::ShapeSweep) = (s.trajectory, s.profile,)

function mesh(g::MeshOptions, s::ShapeSweep, (mt,mp,))
	tlist = paths(mt) # trajectory paths
	plist = paths(mp) # profile paths
	m = nothing
	for t in tlist # each loop in the trajectory
		rmax = sqrt(maximum(norm².(t)))
		ε = max(get(g, :atol), get(g, :rtol)*rmax)
		v = [ surface(Shapes.path_extrude(t, p;
			closed=true, join=s.join, miter_limit=s.miter_limit, precision=ε)...)
			for p in plist ]
		r = reduce(symdiff, v)
		m = (m == nothing) ? r : csgunion(m, r)
	end
	return m
end

# VolumeSweep ««2
struct VolumeSweep <: AbstractTransform{3}
	transform::Function
	nsteps::Int
	maxgrid::Int
	isolevel::Float64
	child::AbstractGeometry{3}
end

function mesh(g::MeshOptions, s::VolumeSweep, (m,))
	return VolumeMesh(TriangleMeshes.swept_volume(m.mesh, s.transform,
		min(ceil(Int,1/get(g,:rtol)),s.nsteps),
		gridcells(g, vertices(m), s.maxgrid), s.isolevel))
end
# Interface ««2
"""
    sweep(path, shape...)

Extrudes the given `shape` by
1) rotating perpendicular to the path (rotating the unit *y*-vector to the direction *z*), and 
2) sweeping it along the `path`, with the origin on the path.

---

    sweep(transform, volume)

Sweeps the given `volume` by applying the `transform`:

  V' = ⋃ { f(t) ⋅ V | t ∈ [0,1] }

`f` is a function mapping a real number in [0,1] to a pair (matrix, vector)
defining an affine transform.

 - nsteps:  upper bound on the number of steps for subdividing the [0,1] interval
 - gridsize: subdivision for marching cubes
 - isolevel: optional distance to add/subtract from swept volume

# Extended help
!!! warning "Swept surfaces: limitations"

    The `sweep` feature currently has two main limitations:
    1. the trajectory may only contain closed loops (no open paths);
    2. the extrusion of each profile vertex must preserve this topology.

    Both of these restrictions are due to the `clipper` library, which
    [does not support single-side offset](http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Types/EndType.htm).
    For an (**experimental**) remedy, see `path_extrude`.

!!! warning "Complexity of swept volume"
    Volume sweep is a costly operation,
    implemented using a marching cubes algorithm on a grid defined
    by `gridsize`.
		Thus, its complexity is **cubic** in the parameter `gridsize`.

FIXME: unify `gridsize` and offset's `maxgrid` parameters.
"""
function sweep end

@inline sweep(path, s...; kwargs...) =
	operator(_sweep, (path,), s...; kwargs...)

@inline _sweep(path, s::AbstractGeometry{2}; kwargs...) =
	ShapeSweep(path, s; kwargs...)
@inline _sweep(transform, s::AbstractGeometry{3};
	nsteps=32, maxgrid=32, isolevel=0) =
	VolumeSweep(transform, nsteps, maxgrid, isolevel, s)


# Path extrusion ««2
struct PathExtrude <: AbstractGeometry{3}
	trajectory::Vector{Vector{SVector{2,Float64}}}
	profile::AbstractGeometry{2}
end
@inline children(s::PathExtrude) = (s.profile,)
function mesh(g::MeshOptions, s::PathExtrude, (m,))
	s = surface(Voronoi.extrude(s.trajectory, paths(m), get(g, :atol))...)
	TriangleMeshes.close_loops!(s.mesh, _DEFAULT_COLOR)
	return csgunion(s, empty(s))
end

"""    path_extrude(trajectory)*profile

Extrusion of the given `profile` along the trajectory.

!!! note "Warning: path extrusion"

    This feature is still experimental. Success is not guaranteed.
    Also, API should not be considered stable.

"""
@inline path_extrude(trajectory, s...) =
	operator(PathExtrude, (trajectory,), s...)
# Generic coordinate transform ««1
struct Deformation{D} <: AbstractTransform{D}
	transform::Function
	isvalid::Function
	distance2::Function
	child::AbstractGeometry{D}
end

function mesh(g::MeshOptions, s::Deformation{3}, (m,))
	mref = TriangleMeshes.splitedges(m.mesh, get(g,:atol)^2;
		distance2=s.distance2)
	vlist = [ (@assert s.isvalid(v); s.transform(v)) for v in mref.vertices ]
	return self_union(VolumeMesh(TriangleMesh(g, vlist, mref.faces, mref.attributes)))
end

"""
    deform(f, s...; isvalid, distance2)

Image of the volume `s` by the vertex-wise transformation `f`
(as a function on `SVector{3}` points).

After applying the transformation, overlapping parts of the volume are
cleaned by a self-union operation. **It is however the user's
responsibility to ensure that the image still forms a valid,
positively-oriented mesh.**

Since `f` is (in principle) not a linear function,
all edges longer than the the current `atol` meshing parameter
will be split, using the `refine` transformation,
before applying `f` to the vertices.

Optional parameters:
 - `isvalid`: a predicate which will be asserted on all points of the solid
   before applying the transformation;
 - `distance2`: a function for evaluating which edges should be split.
   The default is to use the Euclidean distance.

"""
@inline deform(transform, s...; isvalid=(@closure x->true),
	distance2=TriangleMeshes.distance2) =
	operator(Deformation, (transform,isvalid, distance2), s...)

# function wrapd2(r,p0,p1)
# 	v = p1-p0
# 	l0= sqrt(p0[1]^2*v[2]^2/r^2+v[1]^2+v[3]^2)
# 	l1= sqrt(p1[1]^2*v[2]^2/r^2+v[1]^2+v[3]^2)
# 	(iszero(v[1]) || iszero(v[2])) && return l0
# 	return (1/(2*v[1])*(l1*p1[1]-l0*p0[1])
# 	+ (v[1]^2*v[2]^2+v[3]^2*r^2)/(2*v[1]*v[2]*r)
# 	  *log((v[2]*p1[1]+r*l1)/(v[2]*p0[1]+r*l0)))
# end
function wrapsagitta(r,p,q)
# c = chord, s = sagitta, R = radius
# R^2 = c^2/4 + (R-s)^2
# c^2/4 = R^2-R^2-s^2+2Rs = s(2R-s), hence c^2 ≈ 8Rs or s ≈ c^2/8R
# we want d2(p,q) == sagitta^2, hence d2(p,q) = c^4/64R^2
# ≈ ‖pq‖^4/64r^2
	chord =(p[1]-q[1])^2+(p[2]-q[2])^2
	chord*=(1+max(p[1],q[1])/r)^2
	radius = 8*(r+min(p[1],q[1]))
	return (chord/radius)^2
end

"""
    wrap(r, s...)

Wraps the solid `s` around a cylinder with radius `r`
by applying the coordinate transformation
``(x \\cos(y/r), x \\sin(y/r), z)``.
Long edges will be split so that their image resembles the correct spirals.
"""
@inline wrap(r::Real, s...) = deform(
	(@closure p->SA[(r+p[1])*cos(p[2]/r),(r+p[1])*sin(p[2]/r), p[3]]), s...;
	isvalid=(@closure p->p[1]+r ≥ 0),
	distance2=(@closure (p,q)->wrapsagitta(r,p,q)))


# Offset««1

struct SurfaceOffset <: AbstractTransform{2}
	radius::Float64
	ends::Symbol
	join::Symbol
	miter_limit::Float64
	child::AbstractGeometry{2}
	@inline SurfaceOffset(radius::Real, ends::Symbol, join::Symbol,
		miter_limit::Real, child::AbstractGeometry{2}) =
		new(radius, ends, join, miter_limit, child)
end

function mesh(g::MeshOptions, s::SurfaceOffset, (m,))
	ε = max(get(g,:atol), get(g,:rtol) * s.radius)
	return ShapeMesh(Shapes.offset(poly(m), s.radius;
	join = s.join, ends = s.ends, miter_limit = s.miter_limit, precision = ε))
end

struct VolumeOffset <: AbstractTransform{3}
	radius::Float64
	maxgrid::Int
	child::AbstractGeometry{3}
	@inline VolumeOffset(radius::Real, maxgrid::Real, s::AbstractGeometry{3}) =
		new(radius, maxgrid, s)
end

mesh(g::MeshOptions, s::VolumeOffset, (m,)) =
	VolumeMesh(TriangleMeshes.offset(m.mesh, s.radius,
		gridcells(g, vertices(m), s.maxgrid)))

# Front-end««2
"""
    offset(r, solid...; kwargs...)
    offset(r; kwargs...) * solid

Offsets by given radius.
Positive radius is outside the shape, negative radius is inside.

Parameters for 2d shapes:

    ends=:round|:square|:butt|:loop
    join=:round|:miter|:square
    miter_limit=2.0

Parameter for 3d solids:

    maxgrid = 32 # upper bound on the number of cubes used in one direction

!!! warning "Complexity"
    Offset of a volume is a costly operation;
    it is realized using a marching cubes algorithm on a grid defined
    by `maxgrid`.
		Thus, its complexity is **cubic** in the parameter `maxgrid`.
"""
@inline offset(r::Real, s...; kw...) = operator(_offset,(r,), s...; kw...)
@inline _offset(r, s::AbstractGeometry{2};
	ends=:fill, join=:round, miter_limit=2.) =
	SurfaceOffset(r, ends, join, miter_limit, s)
@inline _offset(r, s::AbstractGeometry{3}; maxgrid=32) =
	VolumeOffset(r, maxgrid, s)

"""
    opening(r, shape...; kwargs...)

[Morphological opening](https://en.wikipedia.org/wiki/Opening_(morphology)):
offset(-r) followed by offset(r).
Removes small appendages and rounds convex corners.
"""
@inline opening(r::Real, s...;kwargs...) =
	offset(r;kwargs...)*offset(-r, s...; kwargs...)
"""
    closing(r, shape...; kwargs...)

[Morphological closing](https://en.wikipedia.org/wiki/Closing_(morphology)):
offset(r) followed by offset(-r).
Removes small holes and rounds concave corners.
"""
@inline closing(r::Real, s...; kwargs...) = opening(-r, s...; kwargs...)


# Decimate/subdivide««1
# decimate««2
struct Decimate <: AbstractTransform{3}
	max_faces::Int
	child::AbstractGeometry{3}
end

@inline mesh(g::MeshOptions, s::Decimate, (m,)) =
	VolumeMesh(TriangleMeshes.decimate(m.mesh, s.max_faces))
"""
    decimate(n, surface...)

Decimates a 3d surface to at most `n` triangular faces.
"""
@inline decimate(n::Real, s...) = operator(Decimate, (Int(n),), s...)

# subdivide««2
struct LoopSubdivide <: AbstractTransform{3}
	count::Int
	child::AbstractGeometry{3}
end

@inline mesh(g::MeshOptions, s::LoopSubdivide, (m,)) =
	VolumeMesh(TriangleMeshes.loop(m.mesh, s.count))
"""
    loop_subdivide(n, shape...)

Applies `n` iterations of loop subdivision to the solid.
This does not preserve shape;
instead, it tends to “round out” the solid.
"""
@inline loop_subdivide(n::Integer, s...) = operator(LoopSubdivide, (n,), s...)

# refine««2
struct Refine{D} <: AbstractTransform{D}
	maxlen::Float64
	child::AbstractGeometry{D}
end
@inline mesh(g::MeshOptions, s::Refine{3}, (m,)) =
	VolumeMesh(TriangleMeshes.splitedges(m.mesh, s.maxlen^2))

"""
    refine(maxlen, volume...)

Splits all edges of `volume` repeatedly (preserving global geometry),
until no edge is longer than `maxlen`.
"""
@inline refine(maxlen::Real, s...) = operator(Refine, (maxlen,), s...)

# SetParameters etc.««1
# SetParameters ««2
struct SetParameters{D} <: AbstractTransform{D}
	parameters
	child::AbstractGeometry{D}
	@inline SetParameters(parameters, child::AbstractGeometry{D}) where{D} =
		new{D}(parameters, child)
end

# this object introduces a special case for recursive meshing:
mesh(g::MeshOptions, s::SetParameters) =
	mesh(MeshOptions(g, s.parameters), s.child)

"""
    set_parameters(;atol, rtol, symmetry, icosphere) * solid...

A transformation which passes down the specified parameter values to its
child. Roughly similar to setting `\$fs` and `\$fa` in OpenSCAD.

See `meshing.md` for documentation about specific parameters.
"""
@inline set_parameters(s...; parameters...) =
	operator(SetParameters, (parameters.data,), s...)

# Color ««2
struct Color{D,C<:Colorant} <: AbstractTransform{D}
	color::C
	child::AbstractGeometry{D}
end

# FIXME: `color` is ignored for shapes
@inline mesh(g::MeshOptions, c::Color{2}, (m,)) = m
@inline mesh(g::MeshOptions, c::Color{3}, (m,)) =
	VolumeMesh(TriangleMesh(vertices(m), faces(m), fill(c.color, size(faces(m)))))

"""
    color(c::Colorant, s...)
    color(c::AbstractString, s...)
    color(c::AbstractString, α::Real, s...)
    color(c) * s...
    colorant"color" * s...

Colors objects `s...` in the given color.
"""
@inline color(c::Colors.RGBA, s...) = operator(Color, (c,), s...)
@inline color(c::Colors.RGB, s...) = color(Colors.RGBA(c,1.), s...)
@inline color(c::Union{Symbol,AbstractString}, s...) =
	color(parse(Colorant, c), s...)
@inline color(c::Union{Symbol,AbstractString}, a::Real, s...) =
	color(Colors.coloralpha(parse(Colorant, c), a), s...)

struct RandomColor{D} <: AbstractTransform{D}
	child::AbstractGeometry{D}
end
@inline mesh(g::MeshOptions, c::RandomColor{3}, (m,)) =
	VolumeMesh(TriangleMesh(vertices(m), faces(m),
		Colors.RGBA{N0f8}.(rand(Colors.RGB{N0f8}, size(faces(m))))))
"""
    randomcolor(s...)

Paints each triangle of `s` in a random color. Intended for debugging purposes.
"""
@inline randomcolor(s...) = operator(RandomColor, (), s...)

# Highlight ««2
struct Highlight{D,C<:Colorant} <: AbstractTransform{D}
	color::C
	child::AbstractGeometry{D}
end

@inline mesh(g::MeshOptions, c::Highlight, (m,)) = m
@inline auxmeshes(g::MeshOptions, s::Highlight, m, l) =
	[l...; s.color => raw(m) ]

const _HIGHLIGHT_ALPHA=.15
"""
    highlight(c::Colorant, s)
    highlight(c::AbstractString, s)
    (c::Colorant) % s

Marks an object as highlighted.
This means that the base object will be displayed (in the specified color)
at the same time as all results of operations built from this object.
"""
@inline highlight(c::Colors.RGBA, s...) = operator(Highlight, (c,), s...)
@inline highlight(c::Union{Symbol,AbstractString}, s...) =
	highlight(Colors.parse(Colorant, c), s...)
@inline highlight(c::Colors.RGB, s...) =
	highlight(Colors.RGBA(c, _HIGHLIGHT_ALPHA), s...)

#————————————————————— CSG operations —————————————————————————————— ««1
# https://www.usenix.org/legacy/event/usenix05/tech/freenix/full_papers/kirsch/kirsch.pdf

abstract type AbstractConstructed{S,D} <: AbstractGeometry{D} end

# default interface:
@inline children(s::AbstractConstructed) = s.children
@inline scad_info(::AbstractConstructed{S}) where{S} = (S, ())
@inline AbstractTrees.printnode(io::IO, ::AbstractConstructed{S}) where{S} =
	print(io, S)

"""
		ConstructedSolid{S,V,D}

A type representing CSG operations on solids. `D` is the dimension and
`S` is a symbol representing the operation (union, intersection etc.)
"""
struct ConstructedSolid{S,V,D} <: AbstractConstructed{S,D}
	children::V # Vector{<:AbstractGeometry}, or tuple etc.
end

constructed_solid_type(S::Symbol, T=@closure A->Vector{A}) =
	ConstructedSolid{S,T(AbstractGeometry{D}),D} where{D}

# Generic code for associativity etc.««2
# Macros for defining neutral element etc.
macro mkneutral(op, what, result); quote
	@inline $(esc(op))(neutral, absorb::$what) = $result
	@inline $(esc(op))(absorb::$what, neutral) = $result
	@inline $(esc(op))(x::$what, ::$what) = x
end end
macro mkunit(op,T=AbstractGeometry); quote $(esc(op))(a::$(esc(T))) = a end end
# make operators associative; see definition of + in operators.jl
macro mkassoc(op,T=AbstractGeometry)
	op=esc(op); T=esc(T); quote
	($op)(a::$T, b::$T, c::$T, x...) = Base.afoldl($op, ($op)(($op)(a,b),c),x...)
end end

# Unrolling: prevent nesting of similar associative constructions
"""
		unroll(x::AbstractGeometry, Val(sym1), Val(sym2)...)

Returns either `[x]` or, if `x` is a `ConstructedSolid` matching one of the
symbols `sym1`, `sym2`..., `children(x)`.
"""
@inline unroll(s::AbstractGeometry, ::Val, tail...) = unroll(s, tail...)
@inline unroll(s::AbstractGeometry) = s
@inline unroll(s::ConstructedSolid{S}, ::Val{S}, tail...) where{S} =
	children(s)
@inline unroll2(s::AbstractGeometry, t::AbstractGeometry, tail...) =
	[unroll(s, tail...); unroll(t, tail...)]

# Booleans ««1
# Complement««2
CSGComplement{D} = ConstructedSolid{:complement,Tuple{<:AbstractGeometry{D}},D}
# This is symbolic and replaced at construction time
@inline mesh(g::MeshOptions, ::CSGComplement, (m,)) =
	error("Complements are not supposed to be meshed...")

"""
    complement(x::AbstractGeometry)
    ~x

Returns the complement of `x`, i.e. an object X such that y ∩ X = y ∖ x.

!!! note "Warning: complement"

    Complements are symbolic and only supported as shortcuts
    for some Boolean operations: ∪, ∩, ∖.
    (They would not make any sense in most other constructions anyway).

"""
@inline complement(x::AbstractGeometry{D}) where{D} = CSGComplement{D}((x,))
@inline complement(::EmptyUnion) = EmptyIntersect()
@inline complement(::EmptyIntersect) = EmptyUnion()

# Union««2
CSGUnion = constructed_solid_type(:union)

@inline mesh(g::MeshOptions, s::CSGUnion, mlist) = reduce(csgunion, mlist)
@inline csgunion(s1::S, s2::S) where {S<:VolumeMesh} =
	VolumeMesh(TriangleMeshes.boolean(0, s1.mesh, s2.mesh))
@inline csgunion(s1::ShapeMesh, s2::ShapeMesh) =
	ShapeMesh(Shapes.clip(:union, s1.poly, s2.poly), I)

"""
    union(a::AbstractGeometry...)

Computes the union of several solids. The dimensions must match.
"""
@inline union(a::AbstractGeometry, b::AbstractGeometry) =
	throw(DimensionMismatch("union of 2d and 3d objects not allowed"))
@inline union(a::AbstractGeometry, b::CSGComplement) =
	complement(setdiff(b.children[1],a))
@inline union(a::CSGComplement, b::AbstractGeometry) = union(b,a)
@inline union(a::CSGComplement, b::CSGComplement) =
	complement(intersect(a.children[1], b.children[1]))
@inline union(a::AbstractGeometry{D}, b::AbstractGeometry{D}) where{D} =
	CSGUnion{D}(unroll2(a, b, Val(:union)))

# TODO: make self_union a mesh cleaning function?
@inline self_union(m::ShapeMesh) = ShapeMesh(Shapes.simplify(poly(m)))
@inline self_union(m::VolumeMesh) = csgunion(m, m)

@mkunit union Union{AbstractGeometry,EmptyUnion}
@mkassoc union Union{AbstractGeometry,EmptyUnion}
@mkneutral union EmptyUnion neutral
@mkneutral union EmptyIntersect  absorb

# Intersection««2
CSGInter = constructed_solid_type(:intersection)

@inline mesh(g::MeshOptions, s::CSGInter, mlist) = reduce(csginter, mlist)
@inline csginter(s1::S, s2::S) where {S<:VolumeMesh} =
	VolumeMesh(TriangleMeshes.boolean(1, s1.mesh, s2.mesh))
@inline csginter(s1::ShapeMesh, s2::ShapeMesh) =
	ShapeMesh(Shapes.clip(:intersection, s1.poly, s2.poly), I)

"""
    intersect(a::AbstractGeometry...)

Computes the intersection of several solids.
Mismatched dimensions are allowed; 3d solids will be intersected
with the horizontal plane (as if by the `slice()` operation)
and a 2d intersection will be returned.
"""
@inline intersect(a::AbstractGeometry{D}, b::AbstractGeometry{D}) where{D} =
	CSGInter{D}(unroll2(a, b, Val(:intersection)))
@inline intersect(a::AbstractGeometry{3}, b::AbstractGeometry{2}) =
	intersect(slice(a), b)
@inline intersect(a::AbstractGeometry{2}, b::AbstractGeometry{3}) =
	intersect(a, slice(b))
@inline intersect(a::AbstractGeometry, b::CSGComplement) = setdiff(a, b.children[1])
@inline intersect(a::CSGComplement, b::AbstractGeometry) = intersect(b,a)
@inline intersect(a::CSGComplement, b::CSGComplement) =
	complement(union(a.children[1], b.children[1]))

@mkunit intersect Union{AbstractGeometry,EmptyIntersect}
@mkassoc intersect Union{AbstractGeometry,EmptyIntersect}
@mkneutral intersect EmptyUnion absorb
@mkneutral intersect EmptyIntersect  neutral

# Difference««2
# this is a binary operator:
CSGDiff = constructed_solid_type(:difference, A->SVector{2,<:A})

@inline mesh(g::MeshOptions, s::CSGDiff, (m1,m2,)) = csgdiff(m1, m2)
@inline csgdiff(s1::S, s2::S) where {S<:VolumeMesh} =
	VolumeMesh(TriangleMeshes.boolean(2, s1.mesh, s2.mesh))
@inline csgdiff(s1::ShapeMesh, s2::ShapeMesh) =
	ShapeMesh(Shapes.clip(:difference, s1.poly, s2.poly), I)

"""
    setdiff(a::AbstractGeometry, b::AbstractGeometry)

Computes the difference of two solids.
The following dimensions are allowed: (2,2), (3,3), and (2,3).
In the latter case, the 3d object will be intersected with the horizontal
plane via the `slice()` operation.
"""
@inline setdiff(a::AbstractGeometry{D}, b::AbstractGeometry{D}) where{D} =
	CSGDiff{D}(SA[a,b])
@inline setdiff(a::AbstractGeometry, b::AbstractGeometry,
	c::AbstractGeometry...) = setdiff(setdiff(a,b), c...)
@inline setdiff(a::AbstractGeometry{2}, b::AbstractGeometry{3}) =
	setdiff(a, slice(b))
@inline setdiff(a::AbstractGeometry{3}, b::AbstractGeometry{2}) =
	throw(DimensionMismatch("difference (3d) - (2d) not allowed"))
@inline setdiff(a::AbstractGeometry, b::CSGComplement) =
	intersect(a, b.children[1])
@inline setdiff(a::CSGComplement, b::AbstractGeometry) = setdiff(b,a)
@inline setdiff(a::CSGComplement, b::CSGComplement) =
	complement(setdiff(b.children[1], a.children[1]))

"""
    setdiff([a...], [b...])

Shorthand for `setdiff(union(a...), union(b...))`.
"""
@inline setdiff(x::AbstractVector{<:Union{AbstractGeometry,EmptyUnion}},
                y::AbstractVector{<:Union{AbstractGeometry,EmptyUnion}}) =
	setdiff(union(x...), union(y...))

# Xo(used by extrusion) ««2
@inline symdiff(s1::S, s2::S) where {S<:VolumeMesh} =
	VolumeMesh(TriangleMeshes.boolean(3, s1.mesh, s2.mesh))
# Convex hull ««1
struct CSGHull{D} <: AbstractConstructed{:hull,D}
	children::Vector{<:AbstractGeometry}
end

@inline mesh(g::MeshOptions, ::CSGHull{2}, mlist) =
	ShapeMesh(g, [convex_hull(reduce(vcat, vertices.(mlist)))])

@inline vertices3(m::VolumeMesh) = vertices(m)
@inline vertices3(s::ShapeMesh) = (position(s)([v;0]) for v in vertices(s))

function mesh(g::MeshOptions{T}, ::CSGHull{3}, mlist) where{T}
	v = SVector{3,T}[]
	for m in mlist; push!(v, vertices3(m)...); end
	(p, f) = convex_hull(v)
	return VolumeMesh(g, p, f)
end

"""
    hull(s::AbstractGeometry...)
    hull(s::AbstractGeometry | StaticVector...)

Represents the convex hull of given solids (and, possibly, individual
points).
Mixing dimensions (and points) is allowed.
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

# no @mkunit here!
@mkassoc hull Union{AbstractGeometry,EmptyUnion}
@mkneutral hull EmptyUnion neutral
@mkneutral hull EmptyIntersect  absorb

# Minkowski sum ««1
# Minkowski sum in unequal dimensions is allowed;
# we place the higher-dimensional summand first,
# and its dimension is equal to the dimension of the sum:
struct MinkowskiSum{D,D2} <: AbstractGeometry{D}
	first::AbstractGeometry{D}
	second::AbstractGeometry{D2}
end

@inline children(m::MinkowskiSum) = (m.first, m.second,)

@inline mesh(g::MeshOptions, ::MinkowskiSum{2,2}, (m1,m2)) =
	ShapeMesh(Shapes.minkowski_sum(m1.poly, m2.poly))
@inline mesh(g::MeshOptions, ::MinkowskiSum{3,3}, (m1,m2)) =
	VolumeMesh(TriangleMeshes.minkowski_sum(m1.mesh, m2.mesh))
function mesh(g::MeshOptions{T}, ::MinkowskiSum{3,2}, (m1,m2)) where{T}
	tri = Shapes.triangulate(poly(m2))
	v2 = collect(vertices3(m2))
	return VolumeMesh(TriangleMeshes.minkowski_sum(m1.mesh,
		TriangleMesh{T,MeshColor}(v2, tri,
			fill(first(m1.mesh.attributes), size(tri)))))
end

"""
    minkowski(s1::AbstractGeometry, s2::AbstractGeometry)

Represents the Minkowski sum of given solids.
Mixing dimensions is allowed (and returns a three-dimensional object).
"""
@inline minkowski(a::AbstractGeometry{D}, b::AbstractGeometry{D}) where{D} =
	MinkowskiSum{D,D}(a, b)
@inline minkowski(a::AbstractGeometry{3}, b::AbstractGeometry{2}) =
	minkowski(a,prism(0)*b)
# 	MinkowskiSum{3,2}(a,b)
@inline minkowski(a::AbstractGeometry{2}, b::AbstractGeometry{3})=minkowski(b,a)

@mkunit minkowski
@mkassoc minkowski Union{AbstractGeometry,EmptyUnion}
@mkneutral minkowski EmptyUnion absorb
@mkneutral minkowski EmptyIntersect  absorb

# Overloading Julia operators««1
# backslash replaces U+2216 ∖ SET MINUS, which is not an allowed Julia operator
@inline Base.:\(x::AbstractGeometry, y::AbstractGeometry) = setdiff(x, y)
@inline Base.:~(x::AbstractGeometry) = complement(x)
@inline Base.:-(x::AbstractGeometry, y::AbstractGeometry,
	tail::AbstractGeometry...) = setdiff(x, union(y, tail...))
# this purposely does not define a method for -(x::AbstractGeometry),
# which could also be interpreted as multiplication
@inline Base.:-(x::AbstractGeometry{D}) where{D} = setdiff(intersect(), x)
@inline Base.:-(x::AbstractVector{<:AbstractGeometry},
                y::AbstractVector{<:AbstractGeometry}) =
	setdiff(union(x), union(y))
⋃ = Base.union
⋂ = Base.intersect

@inline Base.:+(v::AbstractVector, x::AbstractGeometry) = translate(v, x)
@inline Base.:+(x::AbstractGeometry, v::AbstractVector) = translate(v, x)
@inline Base.:+(x::AbstractGeometry, y::AbstractGeometry) = minkowski(x,y)
@inline Base.:-(x::AbstractGeometry, v::AbstractVector) = translate(-v, x)

@inline Base.:*(c::Real, x::AbstractGeometry) = scale(c, x)
@inline Base.:*(c::AbstractVector, x::AbstractGeometry) = scale(c, x)
@inline Base.:*(a::AbstractMatrix, x::AbstractGeometry) = mult_matrix(a, x)
@inline _to_matrix(z::Complex) = [real(z) -imag(z); imag(z) real(z)]
@inline Base.:*(z::Complex, x::AbstractGeometry) = _to_matrix(z)*x
@inline Base.:*(a::Transform,b::AbstractMatrix) = a*mult_matrix(b)
@inline Base.:*(a::AbstractMatrix,b::Transform) = mult_matrix(a)*b
@inline Base.:*(a::Transform,b::Real) = a*scale(b)
@inline Base.:*(a::Transform,z::Complex) = a*mult_matrix(_to_matrix(z))

@inline ×(a::Real, x::AbstractGeometry) = prism(a, x)
×(a::AbstractVector{<:Real}, x::AbstractGeometry) =
	if length(a) == 1
		prism(a[1], x)
	elseif length(a) == 2
		raise(a[1], prism(a[2]-a[1], x))
	else
		throw("vector × AbstractGeometry only defined for lengths 1 and 2")
	end

@inline +(a::Angle, x::AbstractGeometry) = rotate(a, x)
@inline ×(a::Angle, x::AbstractGeometry) = revolution(a, x)
@inline ×(a::AbstractVector{<:Angle}, x::AbstractGeometry) =
	if length(a) == 1
		revolution(a[1], x)
	elseif length(a) == 2
		rotate(a[1], revolution(a[2]-a[1], x))
	else
		throw("vector{Angle} × AbstractGeometry only defined for lengths 1 and 2")
	end

@inline Base.:*(c::Symbol, x::AbstractGeometry) = color(String(c), x)
@inline Base.:*(c::Colorant, x::AbstractGeometry) = color(c, x)
@inline Base.:%(c::Symbol, x::AbstractGeometry) = highlight(String(c), x)
@inline Base.:%(c::Colorant, x::AbstractGeometry) = highlight(c, x)

#————————————————————— I/O —————————————————————————————— ««1
function __init__loadsave()
	# FileIO has a central registry: this is bad!
	# The PLY, STL etc. slots are already taken by MeshIO: this is worse!
	# We don't use MeshIO since it needs its own data type (bad!)
	# instead of an interface, and is badly documented.
	# All those bad APIs force us to play dirty:
	for s in (:STL_ASCII, :STL_BINARY, :PLY_ASCII, :PLY_BINARY)
		FileIO.sym2loader[s] = [ConstructiveGeometry]
		FileIO.sym2saver[s] = [ConstructiveGeometry]
	end
end

"""
    load(file)

Reads an object in a data file. Currently supported formats are:
 - `.stl` (ascii and binary)
 - `.ply` (ascii and binary)
"""
@inline load(f::FileIO.File{F}; kw...) where{F<:FileIO.DataFormat} =
	# also, open(::FileIO.File) somehow does not return an IO!? sad!
	open(f) do x; load(F(), x.io; kw...); end

function FileIO.save(f::AbstractString, m::AbstractGeometry; kw...)
	# this is a method of `FileIO.save` since we know the type of `m`
	endswith(f, r".stl"i) && return open(f, "w") do io
		save(format"STL_ASCII"(), io, m; kw...) end
	endswith(f, r".ply"i) && return open(f, "w") do io
		save(format"PLY_ASCII"(), io, m; kw...) end
	endswith(f, r".svg"i) && return open(f, "w") do io
		save(format"SVG"(), io, m; kw...) end
	return Makie.save(f, plot(m; kw...)) # using CairoMakie
end
"""
    FileIO.save(file, object; kw...)

Saves `set_parameters(kw...)*object` in the given file.
Currently supported file formats are:
 - `.stl` (ascii)
 - `.ply` (ascii)
 - `.svg` for shapes.

Other image formats (e.g. `.png`) are supported through
the Makie `plot` function.
"""
@inline save(f::FileIO.DataFormat, io::IO, m::AbstractGeometry; kwargs...) =
	save(f, io, fullmesh(m; kwargs...).main)
@inline save(f::FileIO.DataFormat, io::IO, m::VolumeMesh; kwargs...) =
	error("recursion detected in save($f)")

# .jl file reading (experimental)««1

# FIXME: replace Main by caller module?
# FIXME: add some blob to represent function arguments
"""
		ConstructiveGeometry.load_jl(file::AbstractString, f::Function)

Reads given `file` and returns the union of all top-level `Geometry`
objects (except the results of assignments) found in the file.

```
#### Example: contents of file `example.jl`
C=ConstructiveGeometry.cube(1)
S=ConstructiveGeometry.square(1)
ConstructiveGeometry.circle(3)
S

julia> ConstructiveGeometry.load_jl("example.jl")
union() {
 circle(radius=3.0);
 square(size=[1.0, 1.0], center=false);
}
```

"""
function load_jl(file::AbstractString)
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
function save(::format"STL_ASCII", io::IO, m::VolumeMesh)#««
	print(io, "solid Julia_ConstructiveGeometry_jl_model\n")
	points = TriangleMeshes.vertices(m.mesh)
	faces = TriangleMeshes.faces(m.mesh)
	for f in faces
		tri = (points[f[1]], points[f[2]], points[f[3]])
		n = cross(tri[2]-tri[1], tri[3]-tri[1])
		print(io, "facet normal ", n[1], " ", n[2], " ", n[3], "\nouter loop\n")
	  for p in tri
			print(io, "    vertex ", p[1], " ", p[2], " ", p[3], "\n")
		end
		print(io, "\nendloop\nendfacet\n")
	end
	print(io, "endsolid\n")
end#»»
function load(::format"STL_BINARY", io) #««
	skip(io, 80) # header
	n = read(io, UInt32)
	f = sizehint!(NTuple{3,Int}[], n)
	v = sizehint!(SVector{3,Float64}[], 3n)

	i = 0
	while !eof(io)
		push!(io, (3i+1, 3i+2, 3i+3))
		skip(io, 12) # normal
		for _ in 1:3
			push!(v, Float64.(read(io, NTuple{3,Float32})))
		end
		skip(io, 2)
	end
	return surface(v, f)
end#»»
function load(::format"STL_ASCII", io) #««
	n = 0
	f = NTuple{3,Int}[]
	v = SVector{3,Float64}[]
	while !eof(io)
		l = split(lowercase(readline(io)))
		if !isempty(l) && first(l) == "facet"
			readline(io) # outer loop
			push!(f, length(v) .+ (1,2,3))
			for _ in 1:3
				push!(v, parse.(Float64, split(lowercase(readline(io)))[2:4]))
			end
			readline(io); readline(io) # end loop, end facet
		end
	end
	return surface(v, f)
end#»»
# PLY ««1
function load(::Union{format"PLY_ASCII",format"PLY_BINARY"}, io)#««
	nv = nf = 0
	readline(io) # "ply"
	format = 0
	while !eof(io)
		l = split(lowercase(readline(io)))
		if l[1] == "format"
			if l[2] == "ascii"
				format = 1
			elseif l[2] == "binary_little_endian"
				format = 2
			elseif l[2] == "binary_big_endian"
				format = 3
			else
				error("bad format: $l[2]")
			end
		elseif l[1] == "element"
			if l[2] == "vertex"
				nv = parse(Int, l[3])
			elseif l[2] == "face"
				nf = parse(Int, l[3])
			end
		elseif l[1] == "end_header"
			break
		end
	end
	if format == 1
		v = [ SVector{3}(parse.(Float64, split(lowercase(readline(io)))))
			for _ in 1:nv ]
		f = [ 1 .+ parse.(Int, split(lowercase(readline(io)))[2:end]) for _ in 1:nf]
	elseif format == 2 # little-endian
		v = [ bswap.(ntoh.(read(io, SVector{3,Float32}))) for _ in 1:nv ]
		f = [ 1 .+ Int64.(bswap.(ntoh.(read(io, NTuple{3,Int32})))) for _ in 1:nf ]
	else # format = 3: big-endian
		v = [ ntoh.(read(io, SVector{3,Float32})) for _ in 1:nv ]
		f = [ 1 .+ Int64.(ntoh.(read(io, NTuple{3,Int32}))) for _ in 1:nf ]
	end
	return surface(v,f)
end#»»
function save(::format"PLY_ASCII", io::IO, m::VolumeMesh)#««
	v,f = TriangleMeshes.vertices(m.mesh), TriangleMeshes.faces(m.mesh)
	print(io, "ply\nformat ascii 1.0\nelement vertex ",
		length(v),
	"\nproperty float32 x\nproperty float32 y\nproperty float32 z\nelement face ",
		length(f),
		"\nproperty list uint8 int32 vertex_indices\nend_header\n"
	)
	for x in v
		print(io, x[1], " ", x[2], " ", x[3], "\n")
	end
	for x in f 
		print(io, "3 ", x[1]-1, " ", x[2]-1, " ", x[3]-1, "\n")
	end
end#»»
# SVG ««1
function save(::format"SVG", io::IO, m::ShapeMesh)
	(x,y) = paths(m)[1][1]; rect = MVector(x, x, y, y)
	for (x,y) in Shapes.vertices(m.poly)
		x < rect[1] && (rect[1] = x)
		x > rect[2] && (rect[2] = x)
		y < rect[3] && (rect[3] = y)
		y > rect[4] && (rect[4] = y)
	end
	# a point (x,y) is displayed as (x,-y)
	dx = rect[2]-rect[1]; dy = rect[4]-rect[3]
	λ = .05
	viewbox = (rect[1]-λ*dx, -(rect[4]+λ*dy), (1+2λ)*dx, (1+2λ)*dy)
	print(io, """
<svg xmlns="http://www.w3.org/2000/svg"
  viewBox="$(viewbox[1]) $(viewbox[2]) $(viewbox[3]) $(viewbox[4])">
<!-- A shape with $(length(paths(m))) paths and $(length.(paths(m))) vertices -->
<path fill-rule="evenodd" fill="#999" stroke-width="0"
  d="\n""")
	for p in paths(m)
		print(io, "M ", p[1][1], ",", -p[1][2], " ")
		for q in p[2:end]
			print(io, "L ", q[1], ",", -q[2], " ")
		end
		print(io, "Z\n")
	end
	print(io, """ " /> </svg>\n""")
end
Base.show(io::IO, ::MIME"image/svg+xml", s::AbstractGeometry{2}) =
	save(format"SVG"(), io, s)
# Viewing««1
@inline Base.display(m::AbstractGeometry) = AbstractTrees.print_tree(m)

# TODO: figure out how to rewrite the following using native Makie recipes
@inline plot(g::AbstractGeometry; kwargs...) =
	plot!(Makie.Scene(), g; kwargs...)

@inline plot!(scene::SceneLike, s::AbstractGeometry; kwargs...)=
	plot!(scene, fullmesh(s); kwargs...)

function plot!(scene::SceneLike, m::FullMesh; kwargs...)
	rawmain = raw(m.main)
	for (c,s) in m.aux
		plot!(scene, auxdiff(s, rawmain); color=c, kwargs...)
	end
	plot!(scene, m.main; kwargs...)
end

@inline auxdiff(aux::TriangleMesh, main::TriangleMesh) =
	TriangleMeshes.boolean(2, aux, main)
@inline auxdiff(aux::Shapes.PolygonXor, main::Shapes.PolygonXor) =
	Shapes.clip(:difference, aux, main)

@inline plot!(scene::SceneLike, m::VolumeMesh; kwargs...) =
	plot!(scene, m.mesh; kwargs...)

function plot!(scene::SceneLike, m::TriangleMesh; kwargs...)
	v = TriangleMeshes.vertices(m); f = TriangleMeshes.faces(m);
	a = TriangleMeshes.attributes(m)
	vmat = similar(first(v), 3*length(f), 3)
	for (i, f) in pairs(f), j in 1:3, k in 1:3
		vmat[3*i+j-3, k] = v[f[j]][k]
	end
	fmat = collect(1:size(vmat, 1))
	attr = [ a[fld1(i,3)] for i in 1:size(vmat, 1)]
	Makie.mesh!(scene, vmat, fmat, color=attr; # opaque!
		lightposition=Makie.Vec3f0(5e3,1e3, 10e3), kwargs... )
	return scene
end

function plot!(scene::SceneLike, m::TriangleMesh{<:Real,Nothing};
		kwargs...)
	# this is a special case for aux meshes: we plot them as transparent,
	# and all vertices are the same color
	v = TriangleMeshes.vertices(m); f = TriangleMeshes.faces(m);
	vmat = [ p[j] for p in v, j in 1:3 ]
	fmat = [ q[j] for q in f, j in 1:3 ]
	Makie.mesh!(scene, vmat, fmat; transparency=true,
		lightposition=Makie.Vec3f0(5e3,1e3,10e3), kwargs...)
end

@inline plot!(scene::SceneLike, s::ShapeMesh; kwargs...) =
	plot!(scene, s.poly; kwargs...)

# function Makie.convert_arguments(T::Type{<:Makie.Mesh},p::Shapes.PolygonXor)
# 	# FIXME: what to do if `p.paths` is empty
# 	v = Shapes.vertices(p)
# 	tri = Shapes.triangulate(p)
# 	f = [ t[j] for t in tri, j in 1:3 ]
# 	return Makie.convert_arguments(T, v, f)
# end
# Makie.plottype(::Shapes.PolygonXor) = Makie.Mesh
function plot!(scene::SceneLike, p::Shapes.PolygonXor;
		color=_DEFAULT_COLOR, kwargs...)
	isempty(p.paths) && return scene
	v = Shapes.vertices(p)
	tri = Shapes.triangulate(p)
	m = [ t[j] for t in tri, j in 1:3 ]
	Makie.mesh!(scene, v, m,
		specular=Makie.Vec3f0(0,0,0), diffuse=Makie.Vec3f0(0,0,0),
		color=color; kwargs...)
	return scene
end

# OpenSCAD output««1
# Annotations ««1
# Abstract type ««2
struct Annotate{D,P,A}<:AbstractTransform{D}
	annotation::A
	points::Vector{P}
	child::AbstractGeometry{D}
@inline Annotate(ann::A, p::AbstractVector{P}, g::AbstractGeometry{D}
	) where{D,P,A} = new{D,P,A}(ann, p, g)
end

struct AnnotatedMesh{A,D,T,P} <: AbstractMesh{D,T}
	annotation::A
	points::Vector{P}
	child::AbstractMesh{D,T}
end

@inline Base.map!(f, m::AnnotatedMesh) =
	(map!(f, m.child); map!(f, m.points, m.points))
@inline mesh(g::MeshOptions, a::Annotate) =
	AnnotatedMesh(a.annotation, a.points, mesh(g, a.child))

function plot!(scene::SceneLike, m::AnnotatedMesh)
	plot!(scene, m.child)
	annotate!(scene, m.annotation, m.points)
	return scene
end

# Text annotation ««2
@inline annotate(s::AbstractString, p::AbstractVector{<:Real}, x...) =
	operator(Annotate,(s,[p],), x...)

function annotate!(scene::SceneLike, s::AbstractString, points)
	text!(scene, [s]; position=[(points[1]...,)], align=(:center,:center))
	return scene
end
# # Arrow annotation ««2
# struct ArrowAnnotation{S}
# 	text::S
# end
# @inline arrows(s::AbstractString, pos1::AbstractVector{<:Real},
# 	pos2::AbstractVector{<:Real}, x...) =
# 	operator(Annotate,(ArrowAnnotation(s), [pos1,pos2]), x...)
# 
# function annotate!(scene::SceneLike, a::ArrowAnnotation, points)
# 	p1, p2 = points
# 	c = (p1+p2)/2
# 	v = (p2-p1)/2
# 	text!(scene, [a.text]; position= [(c...,)], align=(:left, :center))
# 	arrows!(scene, [c[1],c[1]], [c[2],c[2]], [c[3],c[3]],
# 		[v[1],-v[1]],[v[2],-v[2]],[v[3],-v[3]])
# 	return scene
# end
# 
# Arrow annotation (TODO) ««2
# arrows!(s,[0],[0],[0],[1],[1],[1],color=:truc)
# ! arrowtail does not seem to work...
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
# Convenience functions ««1
# function rounded_square(dims::AbstractVector, r)
# 	iszero(r) && return square(r)
# 	return offset(r, [r,r]+square(dims - [2r,2r]))
# end
# @inline rounded_square(s::Real, r) = rounded_square(SA[s,s], r)
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
# 	print(io, spaces, "// Object with named anchors:\n")
# 	for (label, anchor) in x.anchors
# 		print(io, spaces, "// $label: $(text(anchor))\n")
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
# 		attach(parent, [:label => child]...)
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
# Exports ««1
export square, circle, stroke, polygon
export cube, sphere, cylinder, cone, surface
export offset, opening, closing, hull, minkowski
export mult_matrix, translate, scale, rotate, reflect, raise, lower
export project, slice, half
export decimate, loop_subdivide, refine
export linear_extrude, prism, revolution, rotate_extrude
export sweep, deform, wrap, path_extrude
export color, randomcolor, highlight, set_parameters
export mesh
export ×

# »»1
function __init__()
	__init__loadsave()
end
end #««1 module
# »»1
# vim: fdm=marker fmr=««,»» noet:
