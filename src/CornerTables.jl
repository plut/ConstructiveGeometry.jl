# TODO:
#  - simplify (retriangulate) faces
#  - maybe precompute all face normals
#  - clarify API
#  - make level computation output-sensitive
"""
    CornerTables

This module implements basic functions for operating with triangular meshes.

The data structure used in this module is an extension of corner tables
(as defined in [Rossignac 2001](https://www.cc.gatech.edu/~jarek/papers/CornerTableSMI.pdf)), allowing the use of non-manifold meshes
(inspired by [Shin et al 2004](https://www.researchgate.net/profile/Hayong_Shin/publication/4070748_Efficient_topology_construction_from_triangle_soup/links/55efd5b408ae199d47c02cd2.pdf)).

This module exports the `CornerTable` data type and defines the following functions (not exported):
 - `faces`, `points`: exports a mesh to a list of triangles.
 - `concatenate`: disjoint union of meshes.
 - `combine`: compute the result of a Boolean operation (intersection, union)
   on meshes.
 - `Base.reverse`: computes the complement of a mesh; together with `combine`,
   this allows computation of a Boolean difference of meshes.
 - `Base.union`, `Base.intersect`, `Base.setdiff`: shorthand for usual combinations.

Additional bibliography:
 - [Rhodes 2013](https://ubm-twvideo01.s3.amazonaws.com/o1/vault/gdc2013/slides/822512 Rhodes_Graham_Math_for_Game \\(2\\).pdf): a general introduction to half-edge data structures;

"""
module CornerTables
using StaticArrays
using LinearAlgebra
using FastClosures
using DataStructures
using StructArrays
module LibTriangle
	using Triangle
end
include("Projectors.jl")
using .Projectors
include("TriangleIntersections.jl")
include("SpatialSorting.jl")
include("EquivalenceStructures.jl")
using .EquivalenceStructures
include("SegmentGraphs.jl")

const _INDEX_TYPE = Int
const _DEFAULT_EPSILON=1/65536

# For a mesh having `n` vertices, we expect the following counts:
#  - 6n + O(1) half-edges;
#  - 3n + O(1) edges;
#  - 2n + O(1) faces.
# The self-intersection of a mesh is generally a curve; hence we expect
# this self-intersection to involve Θ(√n) faces (and the same order of
# magnitude for edges and vertices).
#
# Most algorithms here are output-aware

# tools ««1
# uniquesort!««2
@inline uniquesort!(a) = Base._groupedunique!(sort!(a))
# Box««2
@inline boundingbox(v::AbstractVector...) =
	SpatialSorting.Box{eltype(v)}(min.(v...), max.(v...))

@inline equivalent_points(points, ε=0) =
	equivalence_structure(length(points), SpatialSorting.duplicates(points, ε))
"""
    simplify_points(points, ε)

Removes duplicates from the set of points, returning the set of
substitutions (as an (oldindex => newindex) iterator).
"""
function simplify_points(points, ε=0)
	eqv = equivalent_points(points, ε)
	return (i => j for (i, j) in pairs(representatives(eqv)) if i ≠ j)
end

# initialize one entry in a dictionary (e.g. of vectors) #««2
function create_entry!(dict::SortedDict, key, init)
	idx = findkey(dict, key)
	status((dict, idx)) == 1 && return dict[idx]
	push!(dict, key => init)
	return dict[key]
end
@inline push_entry!(dict::SortedDict, key, value...) =
	push!(create_entry!(dict, key, []), value...)
# LazyMap ««2
struct LazyMap{Y,X,F} <: AbstractVector{Y}
	f::F
	source::X
end
@inline Base.size(l::LazyMap) = size(l.source)
@inline Base.getindex(l::LazyMap, i::Integer) = l.f(l.source[i])

LazyMap{Y}(f, v) where{Y} = LazyMap{Y,typeof(v),typeof(f)}(f, v)
lazymap(f,v) = LazyMap{Base.return_types(f,Tuple{eltype(v)})[1]}(f, v)

# Geometry««2
@inline norm²(v) = dot(v,v)
# NamedIndex««2
struct NamedIndex{S,T<:Signed} i::T; end
@inline Int(i::NamedIndex) = i.i
@inline Base.show(io::IO, i::NamedIndex{S}) where{S} = print(io, S, Int(i))
@inline NamedIndex{S}(x) where{S} = NamedIndex{S,typeof(x)}(x)
@inline Base.isless(x::NamedIndex{S}, y::NamedIndex{S}) where{S} = Int(x)<Int(y)
macro NamedIndex(name,S)
	quote
	const $(esc(name)) = NamedIndex{$(QuoteNode(S))}
	Base.show(io::IO, ::Type{<:$(esc(name))}) = print(io, $(string(name)))
end end
@NamedIndex Vertex v
@NamedIndex Corner c
@NamedIndex Face f
@NamedIndex Fan fan
@NamedIndex Side s
@inline Base.iterate(::Type{Side}, s...) =
	Base.iterate((Side(1),Side(2),Side(3)), s...)

# corner table mesh ««1
# data structure ««2
"""
    CornerTable{I,P,A}

 - `I` is an (integer) index type (usually `Int`; `Int32` is also plausible);
 - `P` is the point type (e.g. `SVector{3,Float64}`);
 - `A` is the per-face attribute type (e.g. a `Colorant` type).

This is an extension of the corner table type of [Rossignac 2001],
allowing to store arbitrary meshes (including: non-manifold, with boundary,
or non-orientable), using ideas from [Shin et al 2004].

In addition to the usual data (geometric point, opposite corner,
and one corner per poitn), the `CornerTable` type also defines *fans*,
which are set of connected corners around a single vertex.
A fan may be either open (i.e. with definite first and last corners)
or closed (i.e. a simple loop of corners around a vertex).
A regular point has a single closed fan.
For simplification, this fan (“implicit fan”) is usually not stored
in the data structure.

                    fan ↺ 
                     ↑↓
  point ←— vertex —→ corner

The fields are as follow:
  - `vertex[v]`:
  
    - `point` = geometric point,
    - `corner` = one corner at this vertex.
	 
  - `corner[c]`: `fan` is fan id, `opposite` is either

    - `+opposite_corner` if the base is a simple edge,
    - `0` if the base is a boundary edge,
    - `-next_corner` if the base is a part of a multiple edge.
  
  - `fan[k]`:

    - `next` is the next fan in cyclic chained list at this vertex,
    - `vertex` is the vertex for this fan,

    - `corner` is either

      - `+first_corner` if the fan is closed,
      - `-first_corner` if the fan is open.

 - `attribute[face]` is the attribute (e.g. color) for this face.

Total size (in number of references) of topological info:
2*ncorners+nvertices+2*nfans = 13*nv + 2*nfans + O(1),
or 13*nv + O(1) for a manifold mesh (only implicit fans).
(+3*nv for geometry, +≈2*nv for attributes)

With explicit fans:
vertex -> firstfan (+nv)
fan->nextfan, firstcorner, vertex (+3nv)
corner -> fan, opposite (+2*6nv)
total = 16*nv+O(1); +23%
"""
struct CornerTable{I<:Signed,P,A,VS<:StructArray,CS<:StructArray,KS<:StructArray}
	vertex::VS
	corner::CS
	fan::KS
	attribute::Vector{A}
	@inline (T::Type{<:CornerTable{I,P,A}})(vertex::VS, corner::CS, fan::KS,
			attribute::AbstractVector{A}) where{I,P,A,VS<:StructArray,
			CS<:StructArray, KS<:StructArray} =
		new{I,P,A,VS,CS,KS}(vertex, corner, fan, attribute)
	@inline (T::Type{<:CornerTable{I,P,A}})(points::AbstractVector) where{I,P,A} =
		T(
			StructArray((point=points, fan=zeros(I,length(points)))),
			StructArray((opposite=I[], fan=I[])),
			StructArray((next=I[], start=I[], vertex=I[])),
			A[])
	@inline CornerTable{I,P,A}() where{I,P,A} = CornerTable{I,P,A}(P[])
	@inline CornerTable{I}(points::AbstractVector) where{I} =
		CornerTable{I,eltype(points),Nothing}(points)
	@inline CornerTable(points::AbstractVector) = CornerTable{Int}(points)
end
# struct CornerTable{I<:Signed,P,A} # I is index type (integer), P is point type
# 	points::Vector{P}
# 	corner::Vector{I}
# 	vertex::Vector{I}
# 	opposite::Vector{I}
# 	fan_next::Vector{I}
# 	fan_first::Vector{I}
# 	attribute::Vector{A}
# 	@inline CornerTable{I,P,A}(points, opp, dest, ef, cn=[], cs=[],a=[]) where{I,P,A} = # TEMPORARY
# 		new{I,P,A}(points, opp, dest, ef, cn, cs,a)
# 	@inline CornerTable(points, opp, dest, ef, cn, cs, a) = 
# 		CornerTable{Int,eltype(points),eltype(a)}(points, opp, dest, ef, cn, cs, a)
# 	@inline CornerTable{I,P,A}(points::AbstractVector) where{I,P,A} =
# 		new{I,P,A}(points, zeros(I, length(points)), [], [], [], [], [])
# 	@inline CornerTable{I,P,A}() where{I,P,A} = CornerTable{I,P,A}(P[])
# 	@inline CornerTable{I}(points::AbstractVector) where{I} =
# 		CornerTable{I,eltype(points),Nothing}(points)
# 	@inline CornerTable(points::AbstractVector) = CornerTable{Int}(points)
# 	# Possible extensions:
# 	# - fan[corner]: number of the (unique) fan containing this corner
# 	#   (we don't keep it: for a manifold mesh, it's a big table of zeros)
# 	# - fan_prev[fan]: to make the list doubly-chained
# end

@inline index_type(::Type{<:CornerTable{I}}) where{I} = I
@inline index_type(m::CornerTable{I}) where{I} = I
# several names for the same type, for clarity:
@inline NamedIndex{S}(m::CornerTable) where{S} = NamedIndex{S,index_type(m)}
point_type(::Type{<:CornerTable{I,P}}) where{I,P} = P
point_type(m::CornerTable) = point_type(typeof(m))

# size and resizing functions««2
# vertices««3
@inline nvertices(m::CornerTable) = length(m.vertex)
@inline nvertices!(m::CornerTable, n::Integer) = (resize!(m.vertex, n); m)
@inline verticeshint!(m::CornerTable, n::Integer) = (sizehint!(m.vertex, n); m)
# @inline nvertices(m::CornerTable) = length(m.points)
# @inline lastvertex(m::CornerTable{I}) where{I} = Vertex{I}(nvertices(m))
# @inline function nvertices!(m::CornerTable, nv::Integer)
# 	resize!(m.points, nv)
# 	resize!(m.corner, nv)
# 	return m
# end
# @inline function verticeshint!(m::CornerTable, nv::Integer)
# 	sizehint!(m.points, nv)
# 	sizehint!(m.corner, nv)
# 	return m
# end
@inline allvertices(m::CornerTable) = (Vertex(v) for v in 1:nvertices(m))
# faces««3
@inline ncorners(m::CornerTable) = length(m.corner)
@inline nfaces(m::CornerTable) = length(m.attribute)
@inline nfaces!(m::CornerTable, n) =
	(resize!(m.corner, 3n); resize!(m.attribute, n); m)
@inline faceshint!(m::CornerTable, n) =
	(sizehint!(m.corner, 3n); sizehint!(m.attribute, 3n); m)


# @inline ncorners(m::CornerTable) = length(m.opposite)
# @inline nfaces(m::CornerTable) = length(m.opposite) ÷ 3
@inline lastface(m::CornerTable{I}) where{I} = Face{I}(nfaces(m))
# @inline function ncorners!(m::CornerTable, nc::Integer)
# 	resize!(m.opposite, nc)
# 	resize!(m.vertex, nc)
# 	resize!(m.attribute, nc ÷ 3)
# 	return m
# end
# @inline nfaces!(m::CornerTable, nf::Integer) = ncorners!(m, 3nf)
# @inline function faceshint!(m::CornerTable, nf::Integer)
# 	sizehint!(m.opposite, 3nf)
# 	sizehint!(m.vertex, 3nf)
# 	sizehint!(m.attribute, nf)
# 	return m
# end
@inline attribute(m::CornerTable, f::Face) = m.attribute[Int(f)]
@inline attribute!(m::CornerTable, f::Face, a) = m.attribute[Int(f)] = a
@inline attribute!(m::CornerTable, a) =
	for f in allfaces(m); attribute!(m, f, a); end

@inline allcorners(m::CornerTable) = (Corner(c) for c in 1:ncorners(m))
@inline allfaces(m::CornerTable) = (Face(f) for f in 1:nfaces(m))
# fans««3
@inline nfans(m::CornerTable) = length(m.fan)
@inline nfans!(m::CornerTable, n::Integer) = (resize!(m.fan, n); m)
@inline fanshint!(m::CornerTable, n::Integer) = (sizehint!(m.fan, n); m)
# @inline nfans(m::CornerTable) = length(m.fan_first)
@inline lastfan(m::CornerTable{I}) where{I} = Fan{I}(nfans(m))
# @inline function nfans!(m::CornerTable, nk::Integer)
# 	resize!(m.fan_first, nk)
# 	resize!(m.fan_next, nk)
# 	return m
# end
# @inline function fanshint!(m::CornerTable, nk::Integer)
# 	sizehint!(m.fan_first, nk)
# 	sizehint!(m.fan_next, nk)
# 	return m
# end
@inline allfans(m::CornerTable) = (Fan(k) for k in 1:nfans(m))

# simple accessors ««2
# corners ««3
@inline next(c::Corner) = Corner(Int(c)+1 - 3*(Int(c)%3==0))
@inline prev(c::Corner) = Corner(Int(c)-1 + 3*(Int(c)%3==1))
@inline getc(m::CornerTable, c::Corner) = LazyRow(m.corner, Int(c))
@inline opposite(m::CornerTable, c::Corner) = Corner(getc(m, c).opposite)
@inline opposite!(m::CornerTable, c::Corner, x::Corner) =
	getc(m, c).opposite = Int(x)
@inline fan(m::CornerTable, c::Corner) = Fan(getc(m, c).fan)
@inline fan!(m::CornerTable, c::Corner, k::Fan) = getc(m, c).fan = Int(k)

# @inline vertex(m::CornerTable, c::Corner) = Vertex(m.vertex[Int(c)])
# @inline vertex!(m::CornerTable, c::Corner, v::Vertex) =
# 	m.vertex[Int(c)] = Int(v)
# @inline getv(m::CornerTable, v::Vertex) = LazyRow(m.vertex, Int(v))
@inline apex(m::CornerTable, c::Corner) = vertex(m, fan(m, c))
@inline right(m::CornerTable, c::Corner) = apex(m, next(c))
@inline left(m::CornerTable, c::Corner) = apex(m, prev(c))
@inline base(m::CornerTable, c::Corner) = (right(m, c), left(m, c))
@inline after(m::CornerTable, c::Corner) = next(opposite(m, next(c)))
@inline after!(m::CornerTable, c::Corner, x::Corner) =
	opposite!(m, next(c), prev(x))
@inline before(m::CornerTable, c::Corner) = prev(opposite(m, prev(c)))

# properties of opposite(c)
@inline issimple(c::Corner) = Int(c) > 0
@inline isboundary(c::Corner) = iszero(Int(c))
@inline ismultiple(c::Corner) = Int(c) < 0

# properties of anyfan(v)
@inline isisolated(k::Fan) = iszero(Int(k))
# @inline isisolated(c::Corner) = iszero(Int(c))


# faces ««3
@inline face(c::Corner) = Face(fld1(Int(c), 3))
@inline side(c::Corner) = Side(mod1(Int(c), 3))
@inline corner(f::Face, s::Side) = Corner(3*Int(f)-3+Int(s))
@inline corners(f::Face) =
	(corner(f,Side(1)), corner(f,Side(2)), corner(f,Side(3)))

@inline vertex(m::CornerTable, f::Face, s::Side) = apex(m, corner(f,s))
@inline point(m::CornerTable, f::Face, s::Side) = point(m, vertex(m, f, s))
@inline vertices(m::CornerTable, f::Face) =
	(vertex(m,f,Side(1)), vertex(m,f,Side(2)), vertex(m,f,Side(3)))
@inline edge(m::CornerTable, f::Face, s::Side) = base(m, corner(f, s))
@inline faces(m::CornerTable) = (vertices(m, f) for f in allfaces(m))
# @inline facevertices(m::CornerTable) =
# 	(Face(f) => vertices(m,Face(f)) for f in 1:nfaces(m))
# @inline function vertices!(m::CornerTable, f::Face, tri::NTuple{3,<:Vertex})
# 	vertex!(m, corner(f,Side(1)), tri[1])
# 	vertex!(m, corner(f,Side(2)), tri[2])
# 	vertex!(m, corner(f,Side(3)), tri[3])
# 	return m
# end

@inline adjacent(m::CornerTable, f::Face, s::Side) =
	face(opposite(m, corner(f, s)))
@inline adjacent(m::CornerTable, f::Face) =
	(adjacent(m, f, Side(1)), adjacent(m, f, Side(2)), adjacent(m, f, Side(3)))
@inline isadjacent(m::CornerTable, f1::Face, f2::Face) =
	any(==(f2), adjacent(m, f1))

# vertices ««3
@inline getv(m::CornerTable, v::Vertex) = LazyRow(m.vertex, Int(v))
@inline anyfan(m::CornerTable, v::Vertex) = Fan(getv(m, v).fan)
@inline anyfan!(m::CornerTable, v::Vertex, k::Fan) = getv(m,v).fan= Int(k)
@inline point(m::CornerTable, v::Vertex) = getv(m, v).point
@inline point!(m::CornerTable, v::Vertex, p) = getv(m, v).point = p
@inline points(m::CornerTable) = m.vertex.point
@inline anycorner(m::CornerTable, v::Vertex) = corner(m, anyfan(m, v))
# @inline anycorner(m::CornerTable, v::Vertex) = getv(m, v).corner
# @inline anycorner!(m::CornerTable, v::Vertex, c::Corner) = getv(m, v).corner = c

@inline Base.map!(f, m::CornerTable) = points(m) .= f.(points(m))

@inline triangle(m::CornerTable, f::Face) =
	(point(m, f, Side(1)), point(m,f,Side(2)), point(m,f,Side(3)))
@inline triangles(m::CornerTable) = (triangle(m, Face(f)) for f in 1:nfaces(m))
@inline function normal(m::CornerTable, f::Face)
	t = triangle(m, f)
	return cross(t[2]-t[1], t[3]-t[2])
end
@inline Projectors.main_axis(m::CornerTable, f::Face) = main_axis(normal(m, f))
function normalized_plane(m::CornerTable, f::Face; absolute=false)
	t = triangle(m, f)
	d = cross(t[2]-t[1], t[3]-t[2])
	a = main_axis(d); absolute && (a = abs(a))
	c = @inbounds inv(d[abs(a)])
	return SA[oftype(c, a), project2d(a, d) .* c..., dot(d, t[1])*c]
end
# manifold without boundary:
# @inline ismanifold(m::CornerTable) = all(m.corner .> 0) && all(m.opposite .> 0)
@inline volume(m::CornerTable) =
	nfaces(m) > 0 ? sum(dot(u, cross(v, w)) for (u,v,w) in triangles(m))/6 :
	zero(eltype(point_type(m)))

# fans ««3
@inline getk(m::CornerTable, k::Fan) = LazyRow(m.fan, Int(k))
@inline nextfan(m::CornerTable, k::Fan) = Fan(getk(m, k).next)
@inline nextfan!(m::CornerTable, k::Fan, x::Fan) = getk(m, k).next = Int(x)
@inline fanstart(m::CornerTable, k::Fan) = getk(m, k).start
@inline fanstart!(m::CornerTable, k::Fan, x::Integer) = getk(m, k).start = x
@inline vertex(m::CornerTable, k::Fan) = Vertex(getk(m, k).vertex)
@inline vertex!(m::CornerTable, k::Fan, v::Vertex) = getk(m, k).vertex = Int(v)

@inline isclosedfan(c::Int) = Int(c) > 0
@inline isopenfan(c::Int) = Int(c) < 0
@inline fan_open(c::Corner)   = -Int(c)
@inline fan_closed(c::Corner) = +Int(c)

const implicitfan = Fan(0)
# @inline fan_first(m::CornerTable, k::Fan) = m.fan_first[Int(k)]
# @inline fan_first!(m::CornerTable, k::Fan, x::Integer) = m.fan_first[Int(k)]= x
# @inline fan_firstcorner(m::CornerTable, k::Fan) = Corner(abs(fan_first(m,k)))
# @inline fanvertex(m::CornerTable, k::Fan) = vertex(m, fan_firstcorner(m, k))

# multiple edges««3
@inline Multi(c::Corner) = Corner(-Int(c))
# create multiple edge - from (at least two) half-edges:
@inline function multi!(m::CornerTable, clist)
	c1 = clist[1]
	for c in clist[2:end]
		opposite!(m, c, Multi(c1))
		c1 = c
	end
	opposite!(m, clist[1], Multi(c1))
end
@inline function add_multi!(m::CornerTable, c::Corner, x::Corner)
	opposite!(m, x, opposite(m, c))
	opposite!(m, c, Multi(x))
end

# iterators««2
struct MeshIterator{S,I,T}
	mesh::CornerTable{I}
	start::T
	@inline MeshIterator{S}(m::CornerTable{I}, s::T) where{S,I,T} =
		new{S,I,T}(m, s)
end
# global COUNT=0
@inline Base.iterate(it::MeshIterator) = (global COUNT=0; (it.start, it.start))
@inline function Base.iterate(it::MeshIterator, s)
# 	global COUNT+=1; @assert COUNT <= 1000
	s = next(it, it.mesh, s)
	stop(it, it.mesh, s) && return nothing
	return (s, s)
end
@inline Base.eltype(::MeshIterator{S,I,T}) where{S,I,T} = T
@inline Base.IteratorSize(::MeshIterator) = Base.SizeUnknown()


@inline star(m::CornerTable, k::Fan) =
	MeshIterator{star}(m, Corner(abs(fanstart(m, k))))
# @inline star(m::CornerTable, c::Corner) = MeshIterator{star}(m, c)
# @inline star(m::CornerTable, v::Vertex) = MeshIterator{star}(m, corner(m, v))
@inline next(it::MeshIterator{star}, m, c) = after(m, c)
@inline stop(it::MeshIterator{star}, m, c) = (c==it.start)||!issimple(c)

@inline fans(m::CornerTable, v::Vertex) = MeshIterator{fans}(m, anyfan(m, v))
@inline Base.iterate(it::MeshIterator{fans}) =
	(iszero(Int(it.start)) ? nothing : (it.start, it.start))
@inline next(it::MeshIterator{fans}, m, k) = nextfan(m, k)
@inline stop(it::MeshIterator{fans}, m, k) = k == it.start

function prevfan(m::CornerTable, k::Fan)
	k2 = k
	while true
		u = nextfan(m, k2)
		u == k && return k2
		k2 = u
	end
end

@inline corners(m::CornerTable, v::Vertex) =
	(c for k in fans(m, v) for c in star(m, k))

# @inline fans(m::CornerTable, k::Fan) = MeshIterator{fans}(m, k)
# @inline next(it::MeshIterator{fans}, m, k) = fan_next(m, k)
# @inline stop(it::MeshIterator{fans}, m, k) = k == it.start
# 
# @inline fancorners(m::CornerTable, k::Fan) =
# 	MeshIterator{fancorners}(m, Corner(abs(fan_first(m, k))))
# @inline next(it::MeshIterator{fancorners}, m, c) = after(m, c)
# @inline stop(it::MeshIterator{fancorners}, m, c) = (c==it.start)||!issimple(c)

@inline radial(m::CornerTable, c::Corner) = MeshIterator{radial}(m, c)
@inline next(it::MeshIterator{radial}, m, c) = Multi(opposite(m, c))
@inline stop(it::MeshIterator{radial}, m, c) = c == it.start
@inline function genradial(m::CornerTable, c::Corner)
	op = opposite(m, c)
	if issimple(op)
		return (c, op)
	else
		return radial(m, c)
	end
end

# "iterator on (fans, corners) around a given vertex"
# @inline function vertexcorners(m::CornerTable, v::Vertex)
# 	c0 = corner(m, v)
# 	isregular(c0) && return ((implicitfan, c) for c in star(m, c0))
# 	isisolated(c0) && return ((implicitfan, c0) for _ in 1:0)
# 	return ((k,c) for k in fans(m, Fan(c0)) for c in fancorners(m, k))
# end
# 
# "finds in which fan a corner lies"
# function fan(m::CornerTable, c0::Corner)
# 	v = vertex(m, c0)
# 	c1 = corner(m, v)
# 	issingular(c1) || return implicitfan
# 	for k in fans(m, Fan(c1)), c in fancorners(m, k)
# 		c == c0 && return k
# 	end
# 	return implicitfan
# end
# 
# function fan_lastcorner(m::CornerTable, k::Fan)
# 	c0 = fan_first(m, k)
# 	for c in fancorners(m, k);  c0 = c; end
# 	return c0
# end

# mesh modification««1
# move_xxx(m, x, y): moves object from index x to index y, overwriting
#                    previous value at y
# delete_xxx(m, x):  deletes object at index x (by moving last onto it)
# vertices««2
function append_points!(m::CornerTable, plist)#««
	Base.IteratorSize(plist) ≠ Base.SizeUnknown() &&
		verticeshint!(m, nvertices(m) + length(plist))
	push!(m.vertex, ((point=p, fan=0) for p in plist)...)
end#»»
@inline function rename_vertex!(m::CornerTable, v::Vertex, x::Vertex, l=nothing)
	for c in corners(m, fan(m, v))
		vertex!(m, c, x)
	end
	l ≠ nothing && replace!(l, v=>x)
	fan!(m, x, fan(m, v))
	point!(m, x, point(m,v))
end
"deletes O(1) vertices at once, by moving the last ones to their place"
function delete_vertex!(m::CornerTable, vlist::Vertex...)
	vlist = MVector(vlist)
	n = nvertices(m)
	for (i, v) in pairs(vlist)
		Int(v)≠n && rename_point!(m, Vertex(n), v, view(vlist, i+1:length(vlist)))
		n-= 1
	end
	nvertices!(m, n)
end
function merge_point!(m::CornerTable, v::Vertex, x::Vertex)
	# deletes vertex v, replacing all references by x
	# we need to precompute the list of incident corners
	# collect() won't work, because Generators has a stupid haslength() type.
	clist = Corner{index_type(m)}[]
	for (_,c) in vertexcorners(m,v); push!(clist, c); end
	for c in clist
		f = face(c); vlist = MVector{3}(vertices(m, f))
		cut_face!(m, f)
		vlist[Int(side(c))] = x
		set_face!(m, f, vlist, attribute(m, f))
	end
	delete_vertex!(m, v)
end

struct Renaming{T}
	newname::SortedDict{T,T,Base.Order.ForwardOrdering}
	oldnames::SortedDict{T,Vector{T},Base.Order.ForwardOrdering}
	@inline Renaming{T}() where{T} =
		new{T}(SortedDict{T,T}(), SortedDict{T,Vector{T}}())
end
@inline newname(r::Renaming, a) = get(r.newname, a, a)
@inline oldnames(r::Renaming, a) = get(r.oldnames, a, [a])
@inline function link!(r::Renaming, (l, a))
	for i in l; push!(r.newname, i => a); end
	push!(r.oldnames, a => l)
end

"removes ε-duplicates from point of m, and returns list of substitutions
(as `oldindex=>newindex` iterator)"
function simplify_points!(m::CornerTable{I}, ε = 0) where{I}
	repl = simplify_points(points(m), ε)
	r = Renaming{Vertex{I}}()
	for (a, b) in repl
		a = Vertex(a); b = Vertex(b)
		a1 = newname(r, a)
		b1 = newname(r, b)
		l = lastvertex(m) # call this before merge_point!
		a1 == b1 && continue
		(u,v) = minmax(a1, b1)
		link!(r, [oldnames(r, a1); oldnames(r, b1)] => u)

		v ≠ l && link!(r, oldnames(r, l) => v)
		merge_point!(m, v, u)
	end
	return r.newname
end
# fans««2
function create_fan!(m::CornerTable, v::Vertex, c::Corner)
# 	println("\e[34;1mcreate_fan!($v, $c)\e[m")
	nfans!(m, nfans(m)+1)
	k = lastfan(m)
	nextfan!(m, k, k)
	fanstart!(m, k, fan_open(c))
	anyfan!(m, v, k)
	vertex!(m, k, v)
	fan!(m, c, k)
	return k
end
# "creates the first fan at vertex `v`"
# function create_fan!(m::CornerTable, v::Vertex, x::Integer)
# # 	println("\e[34;1mcreate_fan!($v)\e[m")
# 	@assert !issingular(corner(m, v)) # this must be the *first* fan
# 	push!(m.fan_first, x)
# 	push!(m.fan_next, nfans(m))
# 	return lastfan(m)
# end
# "appends a new fan to same vertex as `k`"
# function append_fan!(m::CornerTable, k::Fan, x::Integer)
# # 	println("\e[34;1mappend_fan!($k, $x)\e[m")
# 	nfans!(m, nfans(m)+1)
# 	l = lastfan(m)
# 	nextfan!(m, l, nextfan(m, k))
# 	nextfan!(m, k, l)
# 	vertex!(m, l, vertex(m, k))
# 	fanstart!(m, k, x)
# # 	push!(m.fan_first, x)
# # 	push!(m.fan_next, m.fan_next[Int(k)])
# # 	m.fan_next[Int(k)] = nfans(m)
# 	return l
# end
"creates a new single-corner fan at this vertex; returns fan number"
function new_fan!(m::CornerTable, v::Vertex, c::Corner)
# 	println("\e[34;1mnew_fan!($v, $c)\e[m")
	k0 = anyfan(m, v)
	k = create_fan!(m, v, c)
	if !isisolated(k0)
		nextfan!(m, k, nextfan(m, k0))
		nextfan!(m, k0, k)
	end
	return k
# 	elseif isregular(cv)
# 		k = create_fan!(m, v, fan_closed(cv))
# 		corner!(m, v, Corner(k))
# 		return append_fan!(m, k, fan_open(c))
# 	else # issingular(cv)
# 		return append_fan!(m, Fan(cv), fan_open(c))
end
"moves a fan from index k to index x"
function move_fan!(m::CornerTable, k::Fan, x::Fan, t = nothing)
# 	println("\e[34;1mmove_fan!($k, $x)\e[m")
	nextfan!(m, prevfan(m, k), x)
	nextfan!(m, x, nextfan(m, k))
	fanstart!(m, x, fanstart(m, k))
	vertex!(m, x, vertex(m, k))
	anyfan!(m, vertex(m, k), x)
	for c in star(m, k)
		fan(m, c) ≠ k && break
		fan!(m, c, x)
	end
	t ≠ nothing && replace!(t, k => x)
end
function delete_fan!(m::CornerTable, k::Fan, t = nothing)
# 	println("\e[34;1mdelete_fan!($k)\e[m")
	@assert k ≠ implicitfan
	n = nextfan(m, k)
	v = vertex(m, k)
	if n == k # this is the single fan at this point
		anyfan!(m, v, Fan(0))
	else
		# remove from cyclic list
		anyfan!(m, v, n)
		nextfan!(m, prevfan(m, k), n)
	end
	# replace all references to last fan with references to k
	l = lastfan(m)
	l ≠ k && move_fan!(m, l, k, t)
	nfans!(m, nfans(m)-1)
end
"connects two open fans together (k2 after k1)."
function glue_fans!(m::CornerTable, k1, k2, t = nothing)
# 	println("\e[34;1mglue_fans!($k1, $k2) $t\e[m")
	@assert vertex(m, k1) == vertex(m, k2)
	@assert isopenfan(fanstart(m, k1))
	@assert isopenfan(fanstart(m, k2))
	if k1 == k2 # glue an open fan with itself by making it closed
			fanstart!(m, k2, -fanstart(m, k2))
	else # remove fan k2 from linked list
		for c in star(m, k2)
			fan!(m, c, k1)
		end
# 		t ≠ nothing && replace!(t, k2 => k1)
		delete_fan!(m, k2, t)
	end
end

"splits fan `k` at a simple edge, making `c` its new first corner. Returns the pair of fans (before, after)"
@inline function split_fan!(m::CornerTable, k::Fan, c::Corner, t = nothing)
	@assert k ≠ implicitfan
	if isclosedfan(fan_first(m, k))
		fan_first!(m, k, fan_open(c))
	else
		append_fan!(m, k, fan_open(c))
	end
end

# edges««2
function edge_flip!(m::CornerTable, ab)
	#   c        c
	#  ╱ ╲      ╱│╲
	# a———b => a │ b
	#  ╲ ╱      ╲│╱
	#   d        d
	ba = opposite(m, ab); @assert issimple(ba)
	bc = next(ab); ca = next(bc)
	ad = next(ba); db = next(ad)
	a = apex(m, bc)
	b = apex(m, ca)
	c = apex(m, ab)
	d = apex(m, ba)
	cb = opposite(m, bc); @assert issimple(cb)
	da = opposite(m, ad); @assert issimple(da)
	# (ca) and (db) are kept in place (but change opposite vertices)
# 	vertex!(m, ca, d)
# 	vertex!(m, db, c)
	# prev(ca) = bc becomes dc
	# next(ca) = ab becomes ad
	# prev(db) = ad becomes cd
	# next(db) = ba becomes bc
	opposite!(m, bc, ad); opposite!(m, ad, bc)
	opposite!(m, ab, da); opposite!(m, da, ab)
	opposite!(m, ba, cb); opposite!(m, cb, ba)
	# vertices don't change:
	@assert apex(m, bc) == a
	@assert apex(m, ad) == b
	@assert apex(m, ab) == c
	@assert apex(m, ba) == d
	return m
end
function edge_collapse!(m::CornerTable, ab)
	rad = collect(genradial(m, ab))
	a = right(m, ab); b = left(m, ab)
	merge_point!(m, b, a)
	for c in rad
		cut_face!(m, face(c))
	end
	delete_face!(m, face.(rad)...)
end
function face_collapse!(m::CornerTable, f)
	ab = corner(f, Side(3))
	ca = corner(f, Side(2))
	a = apex(m, corner(f, Side(1)))
	b = apex(m, ca)
	c = apex(m, ab)
	rad_ab = collect(genradial(m, ab))
	rad_ac = collect(genradial(m, ca))
	merge_point!(m, c, a)
	merge_point!(m, b, a)
	for x in rad_ab
		cut_face!(m, face(x))
	end
	for x in rad_ac
		cut_face!(m, face(x))
	end
	delete_face!(m, face.(rad_ab)..., face.(rad_ac)...)
end
"replaces vertex b by vertex a"
# function edge_collapse!(m::CornerTable, ab)#««
# 	#   c         c
# 	#  ╱ ╲       ╱
# 	# a———b  => a
# 	#  ╲ ╱       ╲
# 	#   d         d
# 	bc = next(ab); ca = next(bc)
# 	# replace vertex b by a:
# 	a = vertex(m, bc); b = vertex(m, ca)
# 	for x in star(m, ca); vertex!(m, x, a); end
# 	ba = opposite(m, ab); @assert issimple(ba)
# 	ad = next(ba); db = next(ad)
# 	ac = opposite(m, ca); @assert issimple(ac)
# 	cb = opposite(m, bc); @assert issimple(cb)
# 	bd = opposite(m, db); @assert issimple(bd)
# 	da = opposite(m, ad); @assert issimple(da)
# 	# cb becomes ca
# 	opposite!(m, cb, ac); opposite!(m, ac, cb);
# 	# bd becomes ad
# 	opposite!(m, bd, da); opposite!(m, da, bd)
# 	# delete 2 faces and 1 vertex
# 	corner!(m, a, next(ac))
# 	delete_face!(m, face(ab), face(ba))
# 	delete_vertex!(m, b)
# 	return m
# end#»»
"given a corner abc, collapses b, then c, onto a"
# function face_collapse!(m::CornerTable, n)#««
# 	#      e             e
# 	#     ╱ ╲           ╱
# 	#    a———c    =>   a
# 	#   ╱ ╲ ╱ ╲       ╱ ╲
# 	#  f———b———d     f   d
# 	bc = corner(n, Side(1))
# 	ca = corner(n, Side(2))
# 	ab = corner(n, Side(3))
# 	a = vertex(m, bc); b = vertex(m, ca); c = vertex(m, ab)
# 	# this guarantees that all edges from these three vertices are simple:
# 	@assert isregular(corner(m, a))
# 	@assert isregular(corner(m, b))
# 	@assert isregular(corner(m, c))
# 	for x in star(m, ab); vertex!(m, x, a); end
# 	for x in star(m, ca); vertex!(m, x, a); end
# 	cb = opposite(m, bc)
# 	ac = opposite(m, ca)
# 	ba = opposite(m, ab)
# 	bd = next(cb); dc = next(bd)
# 	ce = next(ac); ea = next(ce)
# 	af = next(ba); fb = next(af)
# 	ae = opposite(m, ea)
# 	ec = opposite(m, ce)
# 	cd = opposite(m, dc)
# 	db = opposite(m, bd)
# 	bf = opposite(m, fb)
# 	fa = opposite(m, af)
# 	d = vertex(m, cb)
# 	e = vertex(m, ac)
# 	f = vertex(m, ba)
# 	# (ae, ec) ⇒ (ae, ea)
# 	opposite!(m, ae, ec); opposite!(m, ec, ae)
# 	pec = prev(ec); vertex!(m, pec, a); corner!(m, a, pec)
# 	corner!(m, e, next(ec))
# 	# (cd, db) ⇒ (ad, da)
# 	opposite!(m, cd, db); opposite!(m, db, cd)
# 	vertex!(m, next(cd), a)
# 	vertex!(m, prev(db), a)
# 	corner!(m, d, next(db))
# 	# (bf, fa) ⇒ (af, fa)
# 	opposite!(m, bf, fa); opposite!(m, fa, bf)
# 	vertex!(m, next(bf), a)
# 	corner!(m, f, next(fa))
# 
# 	# delete 4 faces and 2 vertices
# 	delete_face!(m, n, face(cb), face(ba), face(ac))
# 	delete_vertex!(m, b, c)
# end#»»
# mesh queries««2
function findedge(m::CornerTable, u, v)
# 	println("\e[31;1m findedge($u, $v)\e[m")
# 	verbose(m, u)
	for c in corners(m, u)
		right(m, c) == v && return prev(c)
	end
	return Corner(0)
end

"returns (corner, fan at u, fan at v)"
# function findedge(m::CornerTable, u, v, w, clist, klist, i)#««««
# 	# Corner annotation:
# 	# i  : sees (u, v)
# 	# i+3: sees (u, w)
# 	# Fan labeling:
# 	# i  : last out of u before v
# 	# i+6: first out of u after w
# 
# 	# rotates around vertex u
# 	# annotates in clist the following corners:
# 	# i  : a corner which sees (u, v)
# 	# i+3: a corner which sees (u, w)
# 	# i+6: same as next(i) (corner at u after v)
# 	# i+9: same as next(i+3) (corner at u after w)
# 	# and also the following fans:
# 	# i: fan before v
# 	# i+6: fan after w
# 	#
# # 	klist[i] = klist[i+6] = implicitfan
# # 	clist[i] = clist[i+3] = clist[i+6] = clist[i+9] = Corner(0)
# 	klist[i] = klist[i+6] = implicitfan
# 	clist[i] = clist[i+3] = Corner(0)
# 	c0 = corner(m, u)
# 	isisolated(c0) && return
# 	# we can safely assume that vertex `u` has an explicit fan:
# 	@assert issingular(c0) "vertex $u has corner $c0"
# # 	println("rotating around $u/$(nvertices(m)), looking for ($v, $w) (i=$i)")
# 	for k in fans(m, Fan(c0)), c in fancorners(m, k)
# 		n = next(c); r = apex(m, n)
# 		p = prev(c); l = apex(m, p)
# # 		println("  at ($k, $c): $(vertex(m,c)) $(base(m,c)) r=$r, l=$l")
# # 		(r == v) && println("   right is v=$v, storing clist[$i] ← $p")
# # 		(l == v) && println("   left is v=$v, storing klist[$i] ← $k")
# # 		(r == w) && println("   right is w=$w, storing clist[$(i+3)] ← p, klist[$(i+6)] ← $k")
# 		(r == v) && (clist[i] = p)
# 		(l == v) && (klist[i] = k)
# 		(r == w) && (clist[i+3] = p; klist[i+6] = k)
# 	end
# end#»»»»
# corners««2
"moves a corner from index c to index x (used by move_face!)"
function move_corner!(m::CornerTable, c::Corner, x::Corner)
# 	println("move_corner!($c, $x)")
	# replace all references to corner c by x
	v = apex(m, c)
	cv = corner(m, v)
	if cv == c
		corner!(m, v, x)
	elseif issingular(cv)
		for k in fans(m, Fan(cv))
			fan_first(m, k) == fan_closed(c) && fan_first!(m, k, fan_closed(x))
			fan_first(m, k) == fan_open(c) && fan_first!(m, k, fan_closed(x))
		end
	end

	co = opposite(m, c)
	if issimple(co)
		opposite!(m, co, x)
	elseif ismultiple(co)
		for r in radial(m, Multi(co))
			println("in radial($co): $r")
			opposite(m, r) == Multi(c) || continue
			opposite!(m, r, Multi(x))
			break
		end
	end # if isboundary: nothing to do
	opposite!(m, x, opposite(m, c))
	vertex!(m, x, apex(m, c))
end
# edges««2

function match_edge!(m::CornerTable, c, cin, cout)
	println("\e[34;7m match_edge($c $(base(m,c))) cout=$cout cin=$cin\e[m")
# 	verbose(m)
	if !isboundary(cin) # there is already an inner corner
		op_in = opposite(m, cin)
		if isboundary(op_in) # cin is the last edge in its fan
# 			println("\e[31;7mmake new multiple 2-edge $([c,cin])\e[m")
			multi!(m, (c, cin))
		elseif issimple(op_in) # we break a simple edge
# 		println("\e[1mmake new multiple 3-edge $([c,cin,cout]) $(base(m,c))\e[m")
			@assert op_in == cout
			multi!(m, (c, cin, cout))
			split_fan!(m, fan(m, prev(cout)), prev(cout))
			split_fan!(m, fan(m, prev(cin)),  prev(cin))
		else # already a multiple edge
# 			println("\e[1mappend to multiple edge $c\e[m")
			add_multi!(m, Multi(op_in), c)
		end
	elseif !isboundary(cout) # there is an outer corner, but no inner
		op_out = opposite(m, cout)
		@assert !issimple(op_out)
		if ismultiple(op_out)
# 			println("\e[1mappending $c to multiple edge ($cout, $cin)\e[m")
			add_multi!(m, Multi(op_out), c)
		else # op_out is a boundary: we are creating a simple edge
# 			println("\e[1mcreate simple edge $c -- $cout ($cin)\e[m")
			opposite!(m, c, cout); opposite!(m, cout, c)
			glue_fans!(m, fan(m, prev(cout)), fan(m, next(c)))
			glue_fans!(m, fan(m, prev(c)), fan(m, next(cout)))
		end
# 	else
# 		println("\e[1mcin=$cin, cout=$cout, nothing to do\e[m")
	end
end
# function match_edge!(m::CornerTable, klist, clist, i1)#««
# 	i2 = @inbounds (2,3,1)[i1]; i3=@inbounds (3,1,2)[i1]
# 	# corner id:
# 	# i2  :  previous inner corner
# 	# i3+3:  previous outer corner
# 	# fan id:
# 	# i1:    fan before corner
# 	# i1+3:  (new) fan inside corner
# 	# i1+6:  fan after corner
# 	c = clist[i1+6]; cout = clist[i3+3]; cin = clist[i2]
# # 	println("\e[34;7m glue_edge($c $(base(m,c))) cout=$cout cin=$cin\e[m");verbose(m); global ME=deepcopy(m)
# 	if !isboundary(cin) # there is already an inner corner
# 		op_in = opposite(m, cin)
# 		if isboundary(op_in) # op_in is the last edge in its fan
# # 			println("\e[31;7mmake new multiple 2-edge $([c,cin])\e[m")
# 			multi!(m, (c, cin))
# 		elseif issimple(op_in) # we break a simple edge
# # 		println("\e[1mmake new multiple 3-edge $([c,cin,cout]) $(base(m,c))\e[m")
# 			@assert op_in == cout
# 			multi!(m, (c, cin, cout))
# 			split_fan!(m, klist[i2], next(cin))
# 			split_fan!(m, klist[i3+6], next(cout))
# 		else # already a multiple edge
# # 			println("\e[1mappend to multiple edge $c\e[m")
# 			add_multi!(m, Multi(op_in), c)
# 		end
# 	elseif !isboundary(cout) # there is an outer corner, but no inner
# 		op_out = opposite(m, cout)
# 		@assert !issimple(op_out)
# 		if ismultiple(op_out)
# # 			println("\e[1mappending $c to multiple edge ($cout, $cin)\e[m")
# 			add_multi!(m, Multi(op_out), c)
# 		else # op_out is a boundary: we are creating a simple edge
# # 			println("\e[1mcreate simple edge $c -- $cout ($cin)\e[m")
# 			opposite!(m, c, cout); opposite!(m, cout, c)
# 			glue_fans!(m, klist[i2], klist[i2+3], klist)
# 			glue_fans!(m, klist[i3+3], klist[i3+6], klist)
# 		end
# # 	else
# # 		println("\e[1mcin=$cin, cout=$cout, nothing to do\e[m")
# 	end
# end#»»
function cut_edge!(m::CornerTable, c::Corner, t = nothing)#««
	println("\e[35;7m  cut_edge $c = $(base(m,c))\e[m"); # verbose(m)
	# 1. cut fans if needed
	p = prev(c)
	n = next(c)
	@assert isregular(c)
	op = opposite(m, c)
	if isboundary(op) # nothing to do
	elseif issimple(op)
		# simple edge: we need to split fans at both vertices v2 and v3
# 		println("\e[1m  cutting *simple edge* $c--$op = $(base(m,c))\e[m")
		k2 = fan(m, n)
		k3 = fan(m, p)
		split_fan!(m, k2, n)
		split_fan!(m, k3, next(op))
# 		println("\e[33;7mafter split_fans:\e[m"); verbose(m)
		opposite!(m, op, Corner(0))
		opposite!(m, c, Corner(0))
		# delete a simple edge: replace by two boundaries
	else # multiple edge; removing one edge may change it to boundary or simple
		# to simple: only if there are three edges (c, c1, c2) in the loop...
		c1 = Multi(op)
		c2 = Multi(opposite(m, c1))
		println("\e[1mremoving multiple edge $c ($c1, $c2) $(base(m,c))\e[m")
		if c2 == c # create boundary edge
			@assert left(m, c1) == left(m, c)
			println("\e[1m  => create boundary edge $c1 => nothing\e[m")
			opposite!(m, c1, Corner(0))
		else # does this create a simple edge?
			c3 = Multi(opposite(m, c2))
			println("  c3=$c3")
			# ... and both remaining edges (c1 and c2) are opposed
			if c3 == c && left(m, c1) ≠ left(m, c2)
				println("\e[1m  => create simple edge $c1 $c2\e[m")
				println((opposite(m, c1), opposite(m, c2)))
				verbose(m, apex(m, c1))
				verbose(m, apex(m, c2))
				println("*** glue fans: $(fan(m, prev(c2))) $(fan(m, next(c1)))")
				verbose(m, fan(m, prev(c2)))
				verbose(m, fan(m, next(c1)))
				glue_fans!(m, fan(m,prev(c2)), fan(m, next(c1)))
				println("*** glue fans: $(fan(m, prev(c1))) $(fan(m, next(c2)))")
				verbose(m, fan(m, prev(c1)))
				verbose(m, fan(m, next(c2)))
				glue_fans!(m, fan(m,prev(c1)), fan(m, next(c2)))
				opposite!(m, c1, c2)
				opposite!(m, c2, c1)
			else # otherwise, this edge remains multiple, we unlink c
				for c2 in radial(m, c1)
					opposite(m, c2) == Multi(c) || continue
					opposite!(m, c2, op)
				end
			end
		end
	end
	opposite!(m, c, Corner(0))
end#»»
# faces««2
"attaches all three corners and edges of a face"
function set_face!(m::CornerTable{I}, f, vlist, a = nothing) where{I}
	attribute!(m, f, a)
	n = 3*Int(f)-3
	println("\e[32;7mset_face!($f, $vlist)\e[m") #; verbose(m)
	v1 = vlist[1]; v2=vlist[2]; v3=vlist[3]
	# this table holds 12 previous corners + the 3 new we create
	# around vertex 1, old corners are 1,4,7,10, new is 13, etc.
	klist = MVector{9,Fan{I}}(undef)
# 	clist = MVector{9,Corner{I}}(undef)
	cin = (findedge(m, v2, v3), findedge(m, v3, v1), findedge(m, v1, v2))
	cout= (findedge(m, v3, v2), findedge(m, v1, v3), findedge(m, v2, v1))

	for (i, v) in pairs(vlist)
		c = Corner(n+i)
		k = new_fan!(m, v, c)
		opposite!(m, c, Corner(0))
# 		fan!(m, c, k)
	end
	for i in 1:3
		match_edge!(m, Corner(n+i), cin[i], cout[i])
	end
# 	@assert nfaces(m) < 4
# 	verbose(m)
# 	println("""
# $v1 fans: $(klist[[1,4,7]]) in=$(clist[2]) out=$(clist[6])
# $v2 fans: $(klist[[2,5,8]]) in=$(clist[3]) out=$(clist[4])
# $v3 fans: $(klist[[3,6,9]]) in=$(clist[1]) out=$(clist[5])
# """)
# 	println("clist=$clist")
# 	for c in clist; Int(c) > 0 && println("$c: $(vertex(m,c)) $(base(m,c))"); end
# 	println("klist=$klist")
# 	match_edge!(m, klist, clist, 1)
# 	match_edge!(m, klist, clist, 2)
# 	match_edge!(m, klist, clist, 3)
	return m
end

function append_face!(m::CornerTable{I}, vlist, a = nothing) where{I}#««
	nfaces!(m, nfaces(m)+1)
	return set_face!(m, lastface(m), vlist, a)
end#»»
"disconnects all edges and corners of a face"
function cut_face!(m::CornerTable, f::Face)
# 	println("\e[31;7mcut_face!($f $(vertices(m,f)))\e[m"); verbose(m)
	c1 = corner(f, Side(1)); v1 = apex(m, c1)
	c2 = corner(f, Side(2)); v2 = apex(m, c2)
	c3 = corner(f, Side(3)); v3 = apex(m, c3)
	cut_edge!(m, c1);
	cut_edge!(m, c2)
	cut_edge!(m, c3)
	@assert opposite(m, c1) == Corner(0)
	@assert opposite(m, c2) == Corner(0)
	@assert opposite(m, c3) == Corner(0)
	delete_fan!(m, fan(m, c1))
	delete_fan!(m, fan(m, c2))
	delete_fan!(m, fan(m, c3))
# 	valfans(m, f)
# 	println("\e[31;7mafter cut_face!($f $(vertices(m,f))):\e[m")
# 	verbose(m)
# 	@assert f ∈ (Face(1),Face(5),Face(2),Face(6),)
end
"moves a face from index f to index x"
@inline function move_face!(m::CornerTable, f::Face, x::Face, t=nothing)
	for s in Side; move_corner!(m, corner(f, s), corner(x, s)); end
	attribute!(m, x, attribute(m, f))
	t ≠ nothing && replace!(t, f => x)
end
"deletes O(1) faces at once, by moving last faces in their place"
function delete_face!(m::CornerTable, flist::Face...)#««
	flist = MVector(flist)
	n = nfaces(m)
	for (i, f) in pairs(flist)
		Int(f) ≠ n && move_face!(m, Face(n), f, view(flist, i+1:length(flist)))
		n-= 1
	end
	nfaces!(m, n)
end#»»
function replace_face!(m::CornerTable, f::Face, vlist, a = nothing)
	# don't move last face since we are immediately replacing:
	cut_face!(m, f)
	set_face!(m, f, vlist, a)
end
# Constructor««2
function CornerTable{I,P,A}(points, faces,
	attributes = Iterators.repeated(nothing)) where{I,P,A}#««
	m = CornerTable{I,P,A}(points)
	println("faces = $faces")
	for ((v1,v2,v3), a) in zip(faces, attributes)
# 		verbose(m)
		append_face!(m, (Vertex(v1), Vertex(v2), Vertex(v3)), a)
	end
	return m
end#»»
@inline CornerTable{I,P}(points, faces, attributes) where{I,P} =
	CornerTable{I,P,eltype(attributes)}(points, faces, attributes)

@inline CornerTable{I}(points, faces, attributes) where{I} =
	CornerTable{I,eltype(points)}(points, faces, attributes)

@inline CornerTable(points, faces, attributes) =
	CornerTable{Int}(points, faces, attributes)
# coplanar and opposite faces««2
# returns the equivalence relation, as pairs of indices in `flist`:
function coplanar_faces(m::CornerTable, flist, ε = 0)
	@inline box(l,ε) = SpatialSorting.Box(l, l .+ ε)
	boxes = [ box(normalized_plane(m, f, absolute=true), ε) for f in flist ]
	return SpatialSorting.intersections(boxes)
end

function opposite_faces(m::CornerTable)
	r = NTuple{2,Face{index_type(m)}}[]
	for f in allfaces(m)
		c1 = corner(f, Side(1))
		ismultiple(opposite(m, c1)) || continue
		v2 = apex(m, f, Side(2))
		v1 = apex(m, f, Side(1))
		for c in radial(m, c1)
			c <= c1 && continue
			# similarly-oriented faces (including this one) don't count:
			right(m, c) == v2 && continue
			apex(m, c) == v1 && push!(r, (face(c1), face(c)))
		end
	end
	return r
end

function coplanar_clusters(m::CornerTable, ε = 0)
	rel = coplanar_faces(m, allfaces(m), ε)
	rel2= [(Face(x),Face(y)) for (x,y) in rel if isadjacent(m, Face(x),Face(y))]
	return collect(classes(equivalence_structure(rel2)))
end
# reverse, concatenate««2
function Base.reverse(m::CornerTable)#««
	newc = @closure c -> c-1+(c%3) # (1,3,2) permutation
	newopp = @closure c -> sign(c)*newc(abs(c))
	newcorner = @closure c -> (c > 0) ? newc(c) : c

	r = (typeof(m))(points(m))
	ncorners!(r, ncorners(m))
	r.corner .= newcorner.(m.corner)
	r.attribute .= m.attribute
	for n in 3:3:ncorners(r)
		r.opposite[n-2] = newopp(m.opposite[n-2]); r.vertex[n-2] = m.vertex[n-2]
		r.opposite[n-1] = newopp(m.opposite[n-0]); r.vertex[n-1] = m.vertex[n-0]
		r.opposite[n  ] = newopp(m.opposite[n-1]); r.vertex[n  ] = m.vertex[n-1]
	end
	nfans!(r, nfans(m))
	r.fan_next .= m.fan_next
	r.fan_first .= newopp.(m.fan_first)
	return r
end#»»
function concatenate(mlist::CornerTable...)#««
	r = (typeof(first(mlist)))(vcat(points.(mlist)...))
	println("\e[34;7mconcatenate: $(nvertices.(mlist)) $(nfaces.(mlist))\e[m")
	nfaces!(r, sum(nfaces(m) for m in mlist))
	nfans!(r, sum(nfans(m) for m in mlist))
	voffset = 0
	coffset = 0
	koffset = 0
	@inline shift(x,y) = sign(x)*(abs(x)+y)
	for m in mlist
		(nv, nc, nk) = (nvertices(m), ncorners(m), nfans(m))
		for i in 1:nv
			r.corner[i+voffset] = shift(m.corner[i], coffset)
		end
		for i in 1:nc
			r.vertex[i+coffset] = m.vertex[i] + voffset
			r.opposite[i+coffset] = shift(m.opposite[i], coffset)
		end
		for i in 1:nk
			r.fan_next[i+koffset] = m.fan_next[i] + koffset
			r.fan_first[i+koffset] = shift(m.fan_first[i], coffset)
		end
		voffset+= nv
		coffset+= nc
		koffset+= nk
	end
	regularize!(r)
	return r
end#»»
# select faces!««2
"""
    select_faces(m, fkept)

Select faces in `m` from vector `fkept`.
"""
function select_faces(m::CornerTable{I}, fkept) where{I}
	# fkept is a list of kept faces
	# first step: compute vertex/face renaming
	vkept = sizehint!(Vertex{I}[], 3*length(fkept))
	for f in fkept
		push!(vkept, vertices(m, f)...)
	end

	uniquesort!(vkept)

	# FIXME do something about fans!
	r = (typeof(m))(m.points[Int.(vkept)])
	for f in fkept
		(v1, v2, v3) = vertices(m, f)
		u1 = Vertex(searchsortedfirst(vkept, v1))
		u2 = Vertex(searchsortedfirst(vkept, v2))
		u3 = Vertex(searchsortedfirst(vkept, v3))
		append_face!(r, (u1, u2, u3), attribute(m, f))
	end
	return r
end

# self-intersection««1
"""
    insert_point(m, si)

Creates (if needed) a new point in `allpoints` for p,
and returns (in all cases) the index of the new point.

Info stored in `in_face` and `in_edge` is updated according
to intersection type `t`.

If this point already has an index (`idx` ≠ `nothing`)
then nothing new is created
(but `in_face` and `in_edge` are still updated).
"""
function insert_point(m, si, p, idx, f, t, ε)#««
	TI = TriangleIntersections

	# easy case: this is a pre-existing vertex
	TI.isvertex(t) && return something(idx,
		apex(m, f, Side(TI.index(t, TI.isvertex))))

	# other cases: we create a new point
	if TI.isedge(t)
		idx == nothing &&
			(push!(si.points, p); idx = Vertex(length(si.points) + nvertices(m)))
		k = Side(TI.index(t, TI.isedge))
		push!(si.in_edge, (edge(m, f, k)..., idx))
		return idx
	end
	# point is interior
	if idx == nothing
		for j in vertices(m, f)
			norm(point(m, j) - p, Inf) ≤ ε && (idx = j; @goto has_idx)
		end
		push!(si.points, p); idx = Vertex(length(si.points) + nvertices(m))
	end
	@label has_idx
	in_face_f = create_entry!(si.in_face, f, [])
	push!(in_face_f, idx)
	return idx
end#»»
"""
    self_intersect(s)

Returns `(points, in_face, in_edge, faces)` describing the
self-intersection graph of `s`.
"""
function self_intersect(m::CornerTable{I}, ε=0) where{I}#««
	regularize!(m)
	println("REGULARIZED $(nvertices(m)), $(nfaces(m))")
	boxes = [ boundingbox(t...) for t in triangles(m) ]
	si = (points=similar(points(m), 0),
		in_face = SortedDict{Face{I},Vector{Vertex{I}}}(),
		in_edge=NTuple{3,Vertex{I}}[],
		faces = Face{I}[],
		edges = NTuple{2,Vertex{I}}[],
		)

	for (f1, f2) in SpatialSorting.intersections(boxes)
		f1 = Face(f1); f2 = Face(f2)
		isadjacent(m, f1, f2) && continue
		println("inter: ($f1 = $(vertices(m, f1)) $(triangle(m,f1))\n      ($f2=$(vertices(m,f2)) $(triangle(m,f2)))")
		it = TriangleIntersections.inter(triangle(m, f1), triangle(m, f2), ε)
		isempty(it) && continue
		
		# create points as needed, and store their index:
		push!(si.faces, f1, f2)
		idx = MVector{6,Vertex{I}}(undef)
		for (i, (p, (t1, t2))) in pairs(it)
			idx[i] = insert_point(m, si, p, nothing, f1, t1, ε)
			idx[i] = insert_point(m, si, p, idx[i],  f2, t2, ε)
			Int(idx[i]) in (56,410,416,409,417) && println("\e[34;3m created $(idx[i]) from $f1=$(vertices(m,f1)) ∩ $f2=$(vertices(m,f2)):\e[m\n$it")
			if Int(idx[i]) == 410
				for i in (56, 88, 75, 84)
					println("\e[33m$i\e[m $(points(m)[i])")
				end
			end
		end
		length(it) == 2 && push!(si.edges, (idx[1], idx[2]))
	end
	uniquesort!(si.faces)
	return si
end#»»
# subtriangulation««1
# project_and_triangulate ««2
function project_and_triangulate(m::CornerTable, proj, vlist,elist, ε = 0)
	plist = [ proj(point(m, v)) for v in vlist ]
	elist2 = map(e->map(v->Int(searchsortedfirst(vlist, v)), e), elist)

	SegmentGraphs.simplify!(plist, elist2, ε)
	newpoints = inv(proj).(plist[length(vlist)+1:end])
	vlist = [vlist; Vertex.(nvertices(m)+1:nvertices(m)+length(plist)-length(vlist))]
	append_points!(m, newpoints)

	# convert to format used by constrained_triangulation
	vmat = Float64[ p[i] for p in plist, i in 1:2 ]
	emat = [ e[i] for e in elist2, i in 1:2 ]
# 	io = open("/tmp/a", "w")
# 	for (i, v) in pairs(vlist)
# 		(x,y) = vmat[i,:]
# 		println(io,"$x\t$y\t$v")
# 	end
# 	println(io,"\n\n")
# 	for (i1, i2) in eachrow(emat)
# 		(x1,y1) = vmat[i1,:]; (x2,y2) = vmat[i2,:]
# 		println(io,"$x1\t$y1\t$(x2-x1)\t$(y2-y1) # $(vlist[i1])--$(vlist[i2])")
# 	end
# 	close(io)
	tri = LibTriangle.constrained_triangulation(vmat, [1:length(plist)...], emat)
	return [ (vlist[a], vlist[b], vlist[c]) for (a,b,c) in tri ]
end

# subtriangulate««2
function edge_inserts(m, in_edge)
	sort!(in_edge) # this has the effect of grouping points by edge
	start = 1
	V = eltype(eltype(in_edge))
	inserts = (k = NTuple{2,V}[], v = Vector{V}[])
	while start <= length(in_edge)
		# group by in_edge
		stop = start
		v1 = in_edge[start][1]
		v2 = in_edge[start][2]
		while stop≤length(in_edge) && in_edge[stop][1]==v1 && in_edge[stop][2]==v2
			stop+= 1
		end
		# find best coordinate and sort
		vec = point(m,v2) - point(m,v1)
		proj = main_axis(vec)
		ins = [ project1d(proj, point(m,in_edge[i][3])) for i in start:stop-1 ]
		perm = sortperm(ins)

		# push unique vertices in list
		i = in_edge[start-1+perm[1]][3]
		vlist = [i]
		for p in perm[2:end]
			j = in_edge[start-1+p][3]
			j ≠ i && push!(vlist, j)
			i = j
		end

		push!(inserts.k, (v1,v2)); push!(inserts.v, vlist)
		start = stop
	end
	return inserts
end
function subtriangulate!(m::CornerTable{I}, ε=0) where{I}
	si = self_intersect(m, ε)
	# first renumber points, removing duplicates, including in self-intersect:
	println("\e[32;7m  new points: $(length(si.points))\e[m")
	append_points!(m, si.points)
	vmap = simplify_points!(m, ε)
	# antécédents: 56=>56; (410, 416)=>150; (409,417)=>161
	for (k, v) in pairs(vmap)
		Int(v) in (56,150,161) && println("  \e[1m$k => $v\e[m")
	end
	println("\e[33;7m  after simplify_points!: $(nvertices(m))\e[m")
	explain(m, "/tmp/x.scad", scale=3)

	for i in eachindex(si.in_edge)
		(u, v, w) = map(x->get(vmap, x, x), si.in_edge[i])
		(u, v) = minmax(u, v)
		si.in_edge[i] = (u, v, w)
	end
	uniquesort!(si.in_edge)
	for i in eachindex(si.in_face)
		si.in_face[i] = map(x->get(vmap, x, x), si.in_face[i])
	end
	for i in eachindex(si.edges)
		si.edges[i] = extrema(map(x->get(vmap, x, x), si.edges[i]))
	end
	# insert points in edges
	in_edge = edge_inserts(m, si.in_edge)
	uniquesort!(si.edges)

	# determine clusters of coplanar (broken) faces:
	coplanar_rel = coplanar_faces(m, si.faces, ε)
	clusters = equivalence_structure(length(si.faces), coplanar_rel)
	
	# compute points by face
	in_face_v = [ get(si.in_face, f, Vertex{I}[]) for f in si.faces ]
	in_face_e = [ NTuple{2,Vertex{I}}[] for f in si.faces ]
	for (f, vlist, elist) in zip(si.faces, in_face_v, in_face_e)
		ff = vertices(m, f)
		push!(vlist, ff...)
		v = 3
		for u in (1,2,3)
			e = minmax(ff[u], ff[v])
			idx = searchsorted(in_edge.k, e)
			if isempty(idx)
				push!(elist, e)
			else
				in_e = in_edge.v[first(idx)]
				push!(vlist, in_e...)
				j = in_e[1]
				push!(elist, minmax(j, e[1]))
				for i in in_e[2:end]
					push!(elist, minmax(i, j))
					j = i
				end
				push!(elist, minmax(j, e[2]))
			end
			v = u
		end
		uniquesort!(vlist)
		for v in vlist
			r = searchsorted(si.edges, (v,v); by=first)
			for e in si.edges[r]
				isempty(searchsorted(vlist, e[2])) && continue
				push!(elist, e)
			end
		end
	end
	# iterate over all clusters of broken faces
	for icluster in classes(clusters)
		f0 = si.faces[first(icluster)]
		nf0 = normal(m, f0)
		proj = Projector(nf0, point(m, apex(m, corner(f0, Side(1)))))

		# type-stable versions of vcat:
		allvertices = uniquesort!(reduce(vcat, view(in_face_v, icluster)))
		alledges = uniquesort!(reduce(vcat, view(in_face_e, icluster)))

# 		for i in icluster; f = si.faces[i]; println(" $f: $(vertices(m,f)) $(main_axis(m,f)>0) $(in_face_v[i]) $([(Int(x),Int(y)) for (x,y) in in_face_e[i]])"); end

		alltriangles = project_and_triangulate(m, proj, allvertices, alledges, ε)
# 			println("\e[1mcluster $icluster\e[m $proj\e[m")
# 			for i in icluster; f = si.faces[i]; println("  $f: $(vertices(m,f)) $(in_face_v[i])"); end
# 			println("triangles: $alltriangles")
		# determine which points are inside which face:
		ftri = [ proj.(triangle(m,f)) for f in view(si.faces, icluster) ]
		boxes = Vector{SpatialSorting.Box{eltype(eltype(ftri))}}(undef,
			length(allvertices)+length(icluster))
		for (i,v) in pairs(allvertices)
			p = proj(point(m, v))
			boxes[i] = SpatialSorting.Box(p, p .+ ε)
		end
		for (i, t) in pairs(ftri)
			boxes[i+length(allvertices)] = boundingbox(t...)
		end
		vlist = [ Vertex{I}[] for _ in icluster ]
		for (i, j) in extrema.(SpatialSorting.intersections(boxes))
			j -= length(allvertices)
			(i ≤ length(allvertices) && j > 0) || continue
			# vertex i is inside bounding box of face i
			v = allvertices[i]
			p = proj(point(m, v))
			(a, b, c) = (ftri[j][1] - p, ftri[j][2] - p, ftri[j][3] - p)
			f = si.faces[icluster[j]]
			if main_axis(m,f) == main_axis(proj)
				a[1]*b[2]-a[2]*b[1] ≥ -ε || continue
				b[1]*c[2]-b[2]*c[1] ≥ -ε || continue
				c[1]*a[2]-c[2]*a[1] ≥ -ε || continue
			else
				a[1]*b[2]-a[2]*b[1] ≤ -ε || continue
				b[1]*c[2]-b[2]*c[1] ≤ -ε || continue
				c[1]*a[2]-c[2]*a[1] ≤ -ε || continue
			end
			push!(vlist[j], v)
		end
		# apply face refinement:
		for (k, i) in pairs(icluster)
			f = si.faces[i]
			a = attribute(m, f)
			orient = main_axis(m,f) == main_axis(proj)
			isfirst = true
			for tri in alltriangles
				issubset(tri, vlist[k]) || continue
				orient || (tri = reverse(tri))
# 				println("  triangle $tri: (ε=$ε)")
				(p,q,r) = (point(m,tri[1]), point(m,tri[2]), point(m,tri[3]))
# 				println("     $p, $q ,$r\n   $(norm(r-q,Inf)) $(norm(p-r,Inf)) $(norm(p-q,Inf))")
				if norm(cross(q-p, r-p), Inf) ≤ ε
					println("\e[31;1m $tri\e[m")
					println("\e[35;1m $(norm(cross(q-p,r-p),Inf))\e[m")
# 					explain(m, "/tmp/x.scad", scale=3)
# 					@assert norm(cross(q-p, r-p), Inf) > ε
				end
				if isfirst
					cut_face!(m, f); isfirst = false
				else
					nfaces!(m, nfaces(m)+1); f = lastface(m)
					attribute!(m, f, a)
				end
				set_face!(m, f, tri, a)
			end
		end
	end
	regularize!(m)
	return m
end
# intersect with plane ««2
function plane_cut(m::CornerTable, ε = 0)
	# returns a pair (points, segments)
	# with points de-duplicated
	r = maximum(max(abs(p[1]), abs(p[2])) for p in points(m))
	# slice ⊂ [-r,r]×[-r,r] ⊊ triangle [3r,0], [0,±2r]
	bigtriangle = (SA[4r,0,0], SA[-2r,4r,0], SA[-2r,-4r,0])

	pts = SVector{2,eltype(point_type(m))}[]
	seg = Tuple{Int,Int}[]

	for t in triangles(m)
		it = TriangleIntersections.inter(t, bigtriangle)
		n = length(it)
		n ≠ 2 && continue
		push!(pts, SA[it[1][1][1], it[1][1][2]], SA[it[2][1][1], it[2][1][2]])
		push!(seg, (length(pts)-1, length(pts)))
	end
	eqv = equivalent_points(pts, ε)
	for i in 1:length(seg)
		seg[i] = (EquivalenceStructures.class_idx(eqv, seg[i][1]),
		          EquivalenceStructures.class_idx(eqv, seg[i][2]))
	end
	pts2 = [ pts[first(c)] for c in classes(eqv) ]
	return (pts2, uniquesort!(seg))
end

# multiplicity««1
# regular_patches««2
"""
    regular_patches(mesh)

Returns a pair `(label, adjacency)`, where `label` is a labeling of the
faces of `s` in regular (manifold) patches, and `adjacency` is an adjacency
matrix between these patches (containing a half-edge where both patches meet).
"""
function regular_patches(m::CornerTable{I}) where{I}
	label = zeros(Int, nfaces(m))
	todo = Face{I}[]
	adjacency = Matrix{Corner{I}}(undef,0,0)
	n = 0
	for start_face in allfaces(m)
		!iszero(label[Int(start_face)]) && continue
		label[Int(start_face)] = n+= 1
		push!(todo, start_face)
		adjacency = let new_adjacency = similar(adjacency, n, n)
			new_adjacency[1:n-1, 1:n-1] .= adjacency
			fill!(view(new_adjacency, n, :), Corner(0))
			fill!(view(new_adjacency, 1:n-1, n), Corner(0))
			new_adjacency
		end
		while !isempty(todo)
			current_face = pop!(todo)
			for s in Side
				c0 = corner(current_face, s)
				op = opposite(m, c0)
				if issimple(op)
					next_face = face(op)
					if iszero(label[Int(next_face)])
						label[Int(next_face)] = n
						push!(todo, next_face)
					end
				elseif ismultiple(op) # singular edge
					for c1 in radial(m, c0)
						println("radial($c0): $c1")
						l = label[Int(face(c1))]
						!iszero(l) && (adjacency[l,n] = adjacency[n,l] = c0)
					end
				end
			end
		end
	end
	return (label=label, adjacency=adjacency)
end
# circular_sign ««2
"""
    circular_sign(u,v)

Let `α` and `β` be the angles of `u`,`v` in [-π, π[.
This function returns a number <0 iff `α` < `β`, >0 iff `α` > `β`,
and `0` iff `α` == `β`.
"""
@inline function circular_sign(u, v)
# 16 cases to consider: u = (-1, -i, 1, i), same for v
	if u[2] > 0
		v[2] ≤ 0 && return -one(u[2]) # (i,-1), (i,-i), (i,1)
	elseif u[2] < 0
		v[2] > 0 && return one(u[2]) #(-i,i)
	elseif u[2] == 0
		if v[2] == 0
			return sign(v[1]) - sign(u[1]) #(1,1) (1,-1) (-1,1) (-1,-1)
		elseif v[2] > 0
			return one(u[2]) #(-1,i) (1,i)
		else
			# the following line is not needed, but faster than computing det:
			return -sign(u[1]) #(-1,-i) (1,-i) 
		end
	end
	# determinant also works for the following cases:
	# (-1,-i), (-i, -1), (-i, 1), (1, i)
	return u[1]*v[2]-u[2]*v[1]
end
# sort_radial_loop««2
function sort_radial_loop(m::CornerTable, c0, pt3 = nothing)#««
	(v1, v2) = base(m, c0)
	dir3 = point(m, v2) - point(m, v1)
	axis = main_axis(dir3)
	dir2 = project2d(axis, dir3)
	dir2scaled = dir2 ./dot(dir3, dir3)
	# for each adjacent face, compute a (3d) vector which, together with
	# the edge, generates the face (and pointing from the edge to the face):
	# 2d projection of face_vec3 (preserving orientation)
	clist = collect(genradial(m, c0)) # corners around this multiple edge
	# we could use this to determine edge orientation:
	p1 = point(m, v1)
	fv = [ point(m, apex(m, c)) - p1 for c in clist ]
	# face vector, projected in 2d (as viewed from v2 to v1):
	fv2= [ project2d(axis, v) .- dot(v, dir3) .* dir2scaled for v in fv ]
# 	println("\e[7msort_radial_loop $c0=$(base(m,c0)): $dir3 axis $axis\e[m")
# 	for (c, y, z) in zip(clist, fv, fv2)
# 		println("# $c = $(vertex(m,c)) $(base(m, c)): fv3=$y, fv2=$z")
# 	end
	face_cmp = @closure (i1, i2) -> let b = circular_sign(fv2[i1], fv2[i2])
		!iszero(b) && return (b > 0)
		# we want a consistent ordering between coincident faces.
		# the rule is: 1. all arrows point inward the thin wedge (this is
		# useful for later determining multiplicity): positive-oriented faces
		# come *before* negative-oriented faces
		# 2. sort by face number (i.e. decreasing face number for rule 1).
		# we have as input an alternated list starting with a positive edge:
		return (-1)^i1*Int(clist[i1]) < (-1)^i2*Int(clist[i2])
		# half-edge number is proportional to face number
	end
	reorder = sort(1:length(clist), lt=face_cmp)
# 	for i in reorder
# 		c=clist[i]; u=fv2[i]
# 		println("  $c=$(vertex(m,c)) $(base(m,c)): $u")
# 	end
	pt3 == nothing && return clist[reorder]

	# find where `pt3 - v1` inserts in radial loop:
	vec3 = pt3 - p1
	vec2 = project2d(axis, vec3) .- dot(vec3, dir3) .*dir2scaled
	@assert !all(iszero,vec2) "half-edge $h aligned with point $vec3"
	k = searchsorted(fv2[reorder], vec2,
		lt = (u,v)->circular_sign(u,v) > 0)
	@assert k.start > k.stop "impossible to determine location at this edge"
	# possibilities are: (i+1:i) for i in 0..n
	k.stop == 0 && return clist[reorder[end]]
	return clist[reorder[k.stop]]
end#»»
# find_good_edge««2
function find_good_edge(m::CornerTable, i, p)
	# it is possible that all edges lie in the same plane
	# (if this is a flat cell), so we pick, as a second point j,
	# the one which maximizes |y/x|, where
	# y = ‖pi∧pj‖, x=‖pi·ij‖/‖pi‖²
	# in other words, we maximize ‖pi∧pj‖²/‖pi·ij‖²
	# caution: it is possible that pi⋅ij = 0 (edge exactly orthogonal),
	# in this case we must return j immediately
	vpi = point(m, i) - p
	fc = vertexcorners(m, i)
	# in triangle uvw, c is vertex u (edge vw), we want edge uv (corner w)
	((_, c), state) = iterate(fc); c = prev(c); cj = left(m, c)
	vpj = point(m, cj) - p
	xj = dot(vpi, vpj)
	iszero(xj) && return c
	xyj = (xj*xj, norm²(cross(vpi, vpj)))
	best = c
	while true
		u = iterate(fc, state)
		u == nothing && return best
		((_, c), state) = u; c = prev(c); ck = left(m, c)
		vpk = point(m, ck) - p
		xk = dot(vpi, vpk)
		iszero(xk) && return c
		xyk = (xk*xk, norm²(cross(vpi, vpk)))
		if xyk[2]*xyj[1] > xyk[1]*xyj[2]
			best = c; xyj = xyk
		end
	end
	return best
end
# locate_point««2
"""
    locate_point(s, labels, comp, p)

Returns `(face, flag)`, where `flag` is zero if p lies outside this face, and one if p lies inside this face.
"""
function locate_point(m::CornerTable, cc_label, c, p)
	# find closest vertex to p
	closest = 0; z = zero(p[1])
	for (i, q) in pairs(points(m))
		cc_label[i] == c || continue
		z1 = norm²(q-p)
		(closest > 0 && z1 ≥ z) && continue
		closest = i; z = z1
	end
	# find a good edge from closest vertex
	c = find_good_edge(m, Vertex(closest), p)
	y = sort_radial_loop(m, c, p) # y is a half-edge in the radial loop of h
	return (face(y), right(m, y) ≠ right(m, c))
end
# multiplicity ««2
function multiplicity(m::CornerTable{I}) where{I}#««
	rp = regular_patches(m)
	ncomp = size(rp.adjacency, 1)
	# `cells`: allows identification of cells bounded by regular patches
	# cells[2i-1] is inside patch i, cells[2i] is outside it
	# `levels`: level[i] is the multiplicity of cells[2i-1]
	levels = LevelStructure(ncomp)
	# cells is kept implicit for now (this *might* be used later, for
	# reconnecting edges when selecting faces)
	for i1 in 2:ncomp, i2 in 1:i1-1
		c0 = rp.adjacency[i1, i2]
		isboundary(c0) && continue
		# regular components i and j meet at edge eindex
		clist = sort_radial_loop(m, c0)
		n = length(clist)
		# plist is the sorted list of regular patches at this edge
		# dlist is the list of edge orientations
		plist = rp.label[Int.(face.(clist))]
		v2 = right(m, c0)
		olist = [right(m, c) == v2 for c in clist]
		p1 = plist[1]; o1 = olist[1]
# 		println("\e[35;1mcomp ($i1, $i2), corner $c0\e[m sorted radial loop is $clist")
# 		for i in 1:n
# 			println("  $(clist[i]) $(olist[i]):   face $(face(clist[i]))=$(vertices(m,face(clist[i])))  p$(plist[i])")
# 		end
		for i in 2:n
			p2 = plist[i]; o2 = olist[i]
			# patch p2 is positively oriented iff d2==v2, etc.:
			# if p1 is positively oriented (d1==v2) then cell between p1 and p2
			# is 2p1 (otherwise 2p1-1);  if p2 is positively oriented (d2==v2),
			# this same cell is 2p2-1 (else 2p2)
			k = 1-o1-o2
			connect!(levels, p1, p2, k)
			p1 = p2; o1 = o2
		end
	end
	# at this stage, `levels` contains full relative multiplicity info for
	# each connected component of the surface. To complete the data, it
	# remains to organize the components themselves. This needs to be done
	# for each pair of connected components (not as a graph, and both
	# directions of a pair), e.g.:
	# - if c1, c2 ⊂ c3, we don't learn the relation between c1 and c2;
	# - it is possible that both c1 ⊂ c2 and c2 ⊂ c1 (e.g. reversed faces).

	# first identify the connected components:
	# face_cc[i] is (the index of a regular component identifying)
	#               the connected component of face i
	# cc_list is the list of all unique labels of connected components
	# vertex_cc[i] is the label for connected component of vertex i
	# vmax[i] is the vertex with maximal x in connected comp. cc_list[i]
	# cc_nest[i] is the nesting level of connected comp. cc_list[i]
	root!(levels)
	face_cc = levels.parent[rp.label]
	cc_list = unique!(sort(levels.parent))
	vertex_cc = similar(cc_list, nvertices(m))
	for (j, f) in pairs(faces(m)), v in f
		vertex_cc[v] = levels.parent[rp.label[j]]
	end
	vmax = similar(cc_list, Vertex{I})
	for (i, c) in pairs(cc_list)
		b = false; z = point(m,Vertex(1))[1]
		for (v, l) in pairs(vertex_cc); v = Vertex(v)
			l ≠ c && continue
			z1 = point(m, v)[1]
			(b && z1 ≤ z) && continue
			vmax[i] = v; z = z1; b = true
		end
	end
	# compute the nesting level for all connected components
	cc_nest = zeros(Int, length(cc_list))
	for i1 in 1:length(cc_list), i2 in 1:length(cc_list)
		i1 == i2 && continue
		(f, b) = locate_point(m, vertex_cc, i2, point(m, vmax[i1]))
		k = b + levels.level[rp.label[Int(f)]]
		# if k > 0 then i1 inside i2
		cc_nest[i1]+= k
	end

	face_level = Vector{Int}(undef, nfaces(m))
	for f in 1:nfaces(m)
		c = searchsortedfirst(cc_list, face_cc[f])
		face_level[f] = 1 + levels.level[rp.label[f]] + cc_nest[c]
	end
	for (f1, f2) in opposite_faces(m)
		face_level[Int(f1)]-= 1
		face_level[Int(f2)]-= 1
	end
	return face_level
end#»»
# simplification (retriangulate faces) ««1
function regularize!(m::CornerTable, ε = _DEFAULT_EPSILON)
	ε² = ε*ε
	ε3 = ε^(1/3)
	println("regularizing...")
	TI = TriangleIntersections
	# remove all triangular faces with area less than ε
	for f in allfaces(m)
		d = TI.degeneracy(triangle(m, f), ε3)
		d == TI.Constants.invalid && continue
		println("  regularizing degenerate face $f=$(vertices(m,f)): $d")
		if TI.isedge(d)
			edge_collapse!(m, corner(f, Side(TI.index(d, TI.isedge))))
		elseif TI.isvertex(d)
			edge_flip!(m, corner(f, Side(TI.index(d, TI.isvertex))))
		else # interior: all three points confounded
			println("   3 confounded points, merging")
			verbose(m, f)
			verbose(m, apex(m, corner(f, Side(1))))
			verbose(m, apex(m, corner(f, Side(2))))
			verbose(m, apex(m, corner(f, Side(3))))
			face_collapse!(m, f)
		end
	end
end
function simplify!(m::CornerTable, ε = _DEFAULT_EPSILON)
	for cluster in coplanar_clusters(m, ε)
		direction = main_axis(m, first(cluster))
		vlist = Vertex{index_type(m)}[]
		# TODO: detect cluster border and triangulate *that* only
		# (this would work with non-connected cluster as well)
		for f in cluster; push!(vlist, vertices(m, f)...); end
		tri = collect(project_and_triangulate(m, abs(direction),  vlist))
# 		length(tri) == length(cluster) && continue
# 		println("$cluster: $(length(cluster)) => $(length(tri))\n  $cluster\n  $([vertices(m,f) for f in cluster])\n  $tri")
	end
end
# operations««1
"""
    combine(meshes, μ, ε = 0)

Combines all given `meshes` and return a new mesh containing only
those faces with multiplicity `μ`.
The parameter `ε` is the precision used for intersection computation.
"""
function combine(meshes, μ, ε = _DEFAULT_EPSILON)
	newmesh = subtriangulate!(concatenate(meshes...), ε)
	levels = multiplicity(newmesh)
	fkept = [ Face(i) for (i,l) in pairs(levels) if l == μ ]
	return select_faces(newmesh,fkept)
end
@inline Base.union(meshes::CornerTable...; ε = _DEFAULT_EPSILON) =
	combine(meshes, 1, ε)
@inline Base.intersect(meshes::CornerTable...; ε = _DEFAULT_EPSILON) =
	combine(meshes, length(meshes), ε)
@inline Base.setdiff(m1::CornerTable, m2::CornerTable; ε = _DEFAULT_EPSILON)=
	combine([m1, reverse(m2)], 2, ε)
#»»1
function valfans(m::CornerTable, f = Face(0))
	global VF=m
	# check for duplicates
	klist = [ implicitfan for _ in allcorners(m) ]
	for v in allvertices(m), (k, c) in vertexcorners(m, v)
		k1 = klist[Int(c)]
		if k1 ≠ implicitfan
			verbose(m)
			println("\e[31;7;1mCorner $c present in fans $k and $k1\e[m")
			@assert false
		end
		klist[Int(c)] = k
	end
	for c in allcorners(m)
		face(c) == f && continue
		v = apex(m, c)
		k = klist[Int(c)]
		if isregular(corner(m, v))
			k == implicitfan && continue
			println("regular corner $c ($v) has explicit fan $k")
			@assert false
		else
			k == implicitfan || continue
			println("singular corner $c ($v) has implicit fan $k")
			@assert false
		end
	end
end
function validate(s::CornerTable)
	for (i, j) in pairs(s.opposite)
		j ∈ keys(s.opposite) || println("edge h$i: opposite = $j, invalid")
		h = s.opposite[j]
		halfedge(s, j) == reverse(halfedge(s, i)) ||
		println("edge h$i =$(halfedge(s,i)): opposite h$j = $(halfedge(s,j))")
	end
	for (v,h) in pairs(s.edgefrom)
		j = s.destination[prev(s, h)]
		j == v || println("vertex $v: edgefrom = $h comes from $j")
	end
end
function explain(s::CornerTable, h)
	opp=opposite(s, h)
	j=fld1(h, 3)
	print("""
half-edge $h: ->$(destination(s,h)), opp=$opp -> $(destination(s,opp))
  in triangle $j with $(map(i->destination(s,i),(3j-2,3j-1,3j)))
""")
end
function explain(s::CornerTable, io::IO = stdout;
	scale=1, offset=[0.,0.,0.], name="m")
	println(io, """
module $name(pos=$offset, c="gray", s=$scale) {
translate(s*pos) {
""")
	for (i, p) in pairs(points(s))
		println(io, """
translate(s*$(Vector{Float64}(p))) {
	color("red") sphere(1);
	color("black", .8) linear_extrude(1) text("$i", size=5);
}
""")
	end
	println(io, "color(c, .7) polyhedron([")
	join(io, [ " s*$(Vector{Float64}(p))" for p in points(s) ], ",\n")
	println(io, "],[")
	join(io, [ " $(Vector{Int}([f...]) .- 1)" for f in faces(s) ], ",\n")
	println(io, "]); } }\n$name();")
end
@inline explain(s::CornerTable, f::AbstractString; kwargs...) =
	open(f, "w") do io explain(s, io; kwargs...) end
function verbose(m::CornerTable, k::Fan)
	print("\e[35;1m$k\e[m: ")
	print("\e[34m$(vertex(m, k))\e[m")
	print(isopenfan(fanstart(m, k)) ? "open " : "closed ")
	print(join(("$c ($(base(m,c)))" for c in star(m, k)), ","))
	print("  next=\e[35m$(nextfan(m,k))\e[m")
	println()
end
function verbose(m::CornerTable, v::Vertex)
	c = anyfan(m, v)
	print("\e[34;1m", v, "\e[m: ")
	if Int(c) == 0
		println("isolated")
	else
		for k in fans(m, v)
			print("\e[35m $k\e[m(")
			for c in star(m, k)
				print(c, " ")
			end
			println(")")
		end
	end
		
# 	elseif isregular(c)
# 		print("regular first corner $c, ")
# 		print("star = (", join(("$u" for u in star(m, c)),","), ")")
# 		println()
# 	else
# 		print("\e[34;1m$v\e[m: singular [$c] (")
# 		println(join(("$k (first=$(fan_first(m,k)))" for k in fans(m,Fan(c))), " "),")")
# 	end
end
function verbose(m::CornerTable, c::Corner)
	o = Int(opposite(m, c))
	println("  \e[33;1m$c\e[m: \e[32m", fan(m,c), " \e[34m", apex(m, c), "\e[m",
# 				" \e[35m", fan(m, c), "\e[m",
		" edge(", apex(m, next(c)), ",", apex(m, prev(c)), ") ",
		" opposite=", (o > 0 ? "\e[32m" : o < 0 ? "\e[31m" : "\e[38;5;8m"),
		opposite(m, c), "\e[m")
end
function verbose(m::CornerTable, f::Face)
	println("\e[32;1m$f\e[m: ", join(vertices(m, f), ","))
	for c in corners(f)
		verbose(m, c)
	end
end
function verbose(m::CornerTable)
	global V=deepcopy(m)
	println("\e[7mtable with $(ncorners(m)) corners = $(nfaces(m)) faces, $(nvertices(m)) vertices, $(nfans(m)) fans:\e[m")
	for v in allvertices(m)
		verbose(m, v)
	end
	for f in allfaces(m)
		verbose(m, f)
	end
	done = falses(ncorners(m))
	for c0 in allcorners(m)
		done[Int(c0)] && continue
		ismultiple(opposite(m,c0)) || continue
		print("\e[31;1mmultiple edge\e[m:")
		for c in radial(m, c0)
			done[Int(c)] = true
			print(" $c ($(base(m,c))) ")
		end
		println()
	end
	for k in allfans(m)
		verbose(m, k)
	end
end

export CornerTable

end # module
