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
@inline Base.convert(T::Type{<:Integer}, i::NamedIndex) = T(Int(i))
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
@inline nextside(s::Side) = Side(mod1(Int(s)+1, 3))
@inline Base.iterate(::Type{Side}, s...) =
	Base.iterate((Side(1),Side(2),Side(3)), s...)

# data structure and accessors ««1
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

@inline indextype(::Type{<:CornerTable{I}}) where{I} = I
@inline indextype(m::CornerTable{I}) where{I} = I
# several names for the same type, for clarity:
@inline NamedIndex{S}(m::CornerTable) where{S} = NamedIndex{S,indextype(m)}
point_type(::Type{<:CornerTable{I,P}}) where{I,P} = P
point_type(m::CornerTable) = point_type(typeof(m))

# size and resizing functions««2
# vertices««3
@inline nvertices(m::CornerTable) = length(m.vertex)
@inline nvertices!(m::CornerTable, n::Integer) = (resize!(m.vertex, n); m)
@inline function verticeshint!(m::CornerTable, n::Integer)
	sizehint!(m.vertex.point, n)
	sizehint!(m.corner.fan, n)
	return m
end
@inline lastvertex(m::CornerTable{I}) where{I} = Vertex{I}(nvertices(m))
@inline allvertices(m::CornerTable) = (Vertex(v) for v in 1:nvertices(m))
# faces««3
@inline ncorners(m::CornerTable) = length(m.corner)
@inline nfaces(m::CornerTable) = length(m.attribute)
@inline nfaces!(m::CornerTable, n) =
	(resize!(m.corner, 3n); resize!(m.attribute, n); m)
@inline function faceshint!(m::CornerTable, n)
	sizehint!(m.corner.opposite, 3n)
	sizehint!(m.conrer.fan, 3n)
	sizehint!(m.attribute, n)
	return m
end
@inline lastface(m::CornerTable{I}) where{I} = Face{I}(nfaces(m))
struct AllFaces{T<:CornerTable} m::T; end
@inline Base.IteratorSize(a::AllFaces) = Base.SizeUnknown()
@inline Base.iterate(a::AllFaces, i=1) =
	i ≤ nfaces(a.m) ? (Face(i), i+1) : nothing
@inline allfaces(m::CornerTable) = AllFaces(m)

@inline attribute(m::CornerTable, f::Face) = m.attribute[Int(f)]
@inline attribute!(m::CornerTable, f::Face, a) = m.attribute[Int(f)] = a
@inline attribute!(m::CornerTable, a) =
	for f in allfaces(m); attribute!(m, f, a); end

@inline allcorners(m::CornerTable{I}) where{I}=
	(Corner{I}(c) for c in 1:ncorners(m))
# fans««3
@inline nfans(m::CornerTable) = length(m.fan)
@inline nfans!(m::CornerTable, n::Integer) = (resize!(m.fan, n); m)
@inline lastfan(m::CornerTable{I}) where{I} = Fan{I}(nfans(m))
@inline allfans(m::CornerTable{I}) where{I} = (Fan{I}(k) for k in 1:nfans(m))

# simple accessors ««2
# corners ««3
@inline next(c::Corner{I}) where{I} = Corner{I}(Int(c)+1 - 3*(Int(c)%3==0))
@inline prev(c::Corner{I}) where{I} = Corner{I}(Int(c)-1 + 3*(Int(c)%3==1))
@inline next2(c::Corner) = (next(c), prev(c)) # FIXME
@inline getc(m::CornerTable, c::Corner) = LazyRow(m.corner, Int(c))
@inline opposite(m::CornerTable{I}, c::Corner) where{I} =
	Corner{I}(getc(m, c).opposite)
@inline opposite!(m::CornerTable, c::Corner, x::Corner) =
	getc(m, c).opposite = Int(x)
@inline fan(m::CornerTable{I}, c::Corner) where{I} = Fan{I}(getc(m, c).fan)
@inline fan!(m::CornerTable, c::Corner, k::Fan) = begin
  getc(m, c).fan = Int(k)
# 	Int(c) == 1049 && check1049(m)
end
function checkn(m::CornerTable, n)
	ncorners(m) < n && return true
	@assert Int(fan(m, Corner(n))) ≤ nfans(m) "fan($n)=$(fan(m,Corner(n)))/$(nfans(m))"
end
	
# function checkfans(m::CornerTable)
# 	star1 = [ Int[] for _ in allfans(m) ]
# 	for c in allcorners(m)
# 		k = Int(fan(m,c))
# 		k > 0 && push!(star1[k], Int(c))
# 	end
# 	for k in allfans(m)
# 		stark = Int.(collect(star(m,k)))
# 		issubset(star1[Int(k)], stark) && continue
# 		println("\e[31;7m error for fan $k:\n  star=$stark\n  corners $(star1[Int(k)])\e[m")
# 		@assert false
# 	end
# end
function check1049(m::CornerTable)
	@assert isvalid(m)
	checkn(m, 1049)
	checkn(m, 630)
end

# @inline apex(m::CornerTable, c::Corner) = vertex(m, fan(m, c))
@inline apex(m::CornerTable, c::Corner) =
begin
# 	check1049(m)
	k = fan(m, c)
	println("      apex($c): fan=$k/$(nfans(m))")
	vertex(m, k)
end
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

# faces ««3
@inline face(c::Corner) = Face(fld1(Int(c), 3))
@inline side(c::Corner) = Side(mod1(Int(c), 3))
@inline corner(f::Face{I}, s::Side) where{I} = Corner{I}(3*Int(f)-3+Int(s))
@inline corners(f::Face) =
	(corner(f,Side(1)), corner(f,Side(2)), corner(f,Side(3)))

@inline vertex(m::CornerTable, f::Face, s::Side) = apex(m, corner(f,s))
@inline point(m::CornerTable, f::Face, s::Side) = point(m, vertex(m, f, s))
@inline vertices(m::CornerTable, f::Face) =
	(vertex(m,f,Side(1)), vertex(m,f,Side(2)), vertex(m,f,Side(3)))
@inline edge(m::CornerTable, f::Face, s::Side) = base(m, corner(f, s))
struct Faces{I,T<:CornerTable} <: AbstractVector{NTuple{3,I}}
	mesh::T
end
@inline Base.size(f::Faces) = (nfaces(f.mesh),)
@inline Base.getindex(f::Faces, i::Integer) = Int.(vertices(f.mesh, Face(i)))
@inline faces(m::CornerTable) = Faces{indextype(m),typeof(m)}(m)

@inline adjacent(m::CornerTable, f::Face, s::Side) =
	face(opposite(m, corner(f, s)))
@inline adjacent(m::CornerTable, f::Face) =
	(adjacent(m, f, Side(1)), adjacent(m, f, Side(2)), adjacent(m, f, Side(3)))
@inline isadjacent(m::CornerTable, f1::Face, f2::Face) =
	any(==(f2), adjacent(m, f1))

# vertices ««3
@inline getv(m::CornerTable, v::Vertex) = LazyRow(m.vertex, Int(v))
@inline anyfan(m::CornerTable{I}, v::Vertex) where{I} = Fan{I}(getv(m, v).fan)
@inline anyfan!(m::CornerTable, v::Vertex, k::Fan) = getv(m,v).fan= Int(k)
@inline point(m::CornerTable, v::Vertex) = getv(m, v).point
@inline point!(m::CornerTable, v::Vertex, p) = getv(m, v).point = p
@inline points(m::CornerTable) = m.vertex.point

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
@inline nextfan(m::CornerTable{I}, k::Fan) where{I} = Fan{I}(getk(m, k).next)
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

# multiple edges««3
@inline Multi(c::Corner{I}) where{I} = Corner{I}(-Int(c))
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
struct TypedIterator{T,X}
	it::X
end
@inline Base.iterate(it::TypedIterator, s...) = iterate(it.it, s...)
@inline Base.eltype(::TypedIterator{T}) where{T} = T
@inline TypedIterator{T}(it) where{T} = TypedIterator{T,typeof(it)}(it)

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
@inline Base.eltype(::Type{MeshIterator{S,I,T}}) where{S,I,T} = T
@inline Base.eltype(it::MeshIterator) = eltype(typeof(it))
@inline Base.IteratorSize(::MeshIterator) = Base.SizeUnknown()
# @inline Base.eltype(x::Base.Generator{T}) where{T<:MeshIterator} =
# 	first(return_types(x.f,(eltype(T),)))


@inline star(m::CornerTable{I}, k::Fan) where{I} =
	MeshIterator{star}(m, Corner{I}(abs(fanstart(m, k))))
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
	Iterators.flatten(star(m, k) for k in fans(m, v))
@eval Base 

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


# mesh modification (low-level)««1
# move_xxx(m, x, y): moves object from index x to index y, overwriting
#                    previous value at y
# delete_xxx(m, x):  deletes object at index x (by moving last onto it)
# vertices««2
function append_points!(m::CornerTable, plist)#««
	length(plist) == 0 && return
	@debug "appending $(length(plist)) points: $(nvertices(m)+1)..$(nvertices(m)+length(plist))"
	Base.IteratorSize(plist) ≠ Base.SizeUnknown() &&
		verticeshint!(m, nvertices(m) + length(plist))
	push!(m.vertex, ((point=p, fan=0) for p in plist)...)
end#»»
@inline function rename_vertex!(m::CornerTable, v::Vertex, x::Vertex, l=nothing)#««
	for c in corners(m, v)
		k = fan(m, c)
		Int(k) ≠ 0 && vertex!(m, k, x)
	end
	l ≠ nothing && replace!(l, v=>x)
	anyfan!(m, x, anyfan(m, v))
	point!(m, x, point(m,v))
end#»»
"deletes a vertex by swap-and-pop; modifies `l` to account for this"
function delete_vertex!(m::CornerTable, v::Vertex, l=nothing)#««
	@debug "deleting vertex $v/$(nvertices(m))"
	# FIXME: this may leak fans
	w = lastvertex(m)
	if v ≠ w
		rename_vertex!(m, w, v, l)
	end
	nvertices!(m, nvertices(m)-1)
end#»»

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
end
"moves a fan from index k to index x"
function move_fan!(m::CornerTable, k::Fan, x::Fan, t = nothing)
	println("\e[34;1mmove_fan!($k, $x)\e[m")
	verbose(m, k)
	verbose(m, x)
	nextfan!(m, prevfan(m, k), x)
	nextfan!(m, x, nextfan(m, k))
	fanstart!(m, x, fanstart(m, k))
	vertex!(m, x, vertex(m, k))
	anyfan!(m, vertex(m, k), x)
# 	Int(k) == 231 && println("star($k) = $(collect(star(m,k)))")
# 	Int(k) == 231 && println("fan(630) = $(fan(m,Corner(630)))")
	for c in star(m, k)
		println("    change fan for corner $c (fan=$(fan(m,c))) from $k to $x")
		fan(m, c) ≠ k && break
		println("    calling fan!($c, $x)")
		fan!(m, c, x)
		println("    now fan($c) = $(fan(m, c))")
	end
	t ≠ nothing && replace!(t, k => x)
	println("\e[34;3m move_fan($k) done\e[m")
end
function delete_fan!(m::CornerTable, k::Fan, t = nothing)
	println("\e[34;1mdelete_fan!($k / $(nfans(m)))\e[m")
	@assert k ≠ implicitfan
	verbose(m, k)
	n = nextfan(m, k)
	v = vertex(m, k)
	verbose(m, v)
	println("vertex $v, next fan is $n")
	if n == k # this is the single fan at this point
		anyfan!(m, v, Fan(0))
	else
		# remove from cyclic list
		anyfan!(m, v, n)
		nextfan!(m, prevfan(m, k), n)
		verbose(m, v)
	end
	# replace all references to last fan with references to k
	l = lastfan(m)
# 	Int(k) == 133 && nfans(m) == 231 &&
# 		println("\e[1;7m  fan(630)=$(fan(m,Corner(630)))\e[m")
	l ≠ k && move_fan!(m, l, k, t)
# 	Int(k) == 133 && nfans(m) == 231 &&
# 		println("\e[1;7m  fan(630)=$(fan(m,Corner(630)))\e[m")
	nfans!(m, nfans(m)-1)
	println("\e[34;3mfan $k deleted\e[m")
end
@inline Base.show(io::IO, ::Type{<:CornerTable}) = print(io, "CTable")
"connects two open fans together (k2 after k1).
Assumes the corresponding edges are already connected (for `star`)."
function glue_fans!(m::CornerTable, k1, k2, t = nothing)
	println("\e[34;1mglue_fans!($k1, $k2) $t\e[m")
	println("   common vertex is $(vertex(m, k1))")
	println("   star($k1) is $(collect(star(m,k1)))")
	println("   star($k2) is $(collect(star(m,k2)))")
	@assert vertex(m, k1) == vertex(m, k2)
	@assert isopenfan(fanstart(m, k1))
	@assert isopenfan(fanstart(m, k2))
	if k1 == k2 # glue an open fan with itself by making it closed
			fanstart!(m, k2, -fanstart(m, k2))
	else # remove fan k2 from linked list
		println(" modifying corners in star($k2)=$(collect(star(m,k2)))")
		for c in star(m, k2)
			println("    fan!($c) <- $k1")
			fan!(m, c, k1)
		end
		println("  all corners are modified")
		verbose(m, vertex(m, k1))
# 		t ≠ nothing && replace!(t, k2 => k1)
		delete_fan!(m, k2, t)
	# Note: k1 may now point to a non-existent fan (if k1 was == last,
	# delete_fan! has moved it).
	end
	println("\e[34;3m glue_fans($k1,$k2) done\e[m")
end

"splits fan `k` at a simple edge, making `c` its new first corner."
@inline function split_fan!(m::CornerTable, k::Fan, c::Corner)
	@assert k ≠ implicitfan
# 	println("cutting fan $k (around $(vertex(m,k)), start = $(fanstart(m, k)) at $c")
	if isclosedfan(fanstart(m, k))
# 		println("  => making it open at $c")
		fanstart!(m, k, fan_open(c))
# 		verbose(m, k)
	else
# 		println("  => making new fan: ")
		k1 = new_fan!(m, vertex(m, k), c)
		for c1 in star(m, k1)
			fan!(m, c1, k1)
		end
	end
end
"removes the corner `c` from its fan and adjusts (splits/removes) it."
function unfan_corner!(m::CornerTable, c::Corner)
	TODO
	# TODO
end

# mesh queries««2
function findedge(m::CornerTable, u::Vertex, v::Vertex)
	println("\e[33;1;7m findedge($u, $v)/$(nvertices(m))\e[m")
	ncorners(m) ≥ 1049 && println("  fan(1049)=$(fan(m,Corner(1049)))")
	for k in fans(m, u); println("$k: $(collect(star(m, k)))"); end
# 	verbose(m, u)
# 	for k in fans(m, u)
# 		verbose(m, k)
# 	end
	for k in fans(m, u), c in star(m, k)
		println("** $k -> $c/$(ncorners(m))")
		println("   after $c = $(after(m, c))")
	ncorners(m) ≥ 1049 && println("    fan(1049)=$(fan(m,Corner(1049)))")

# 		println("  corner $c = $(apex(m,c))->$(base(m,c)) has fan $(fan(m, c))")
		right(m, c) == v && return prev(c)
	end
	return Corner(0)
end

# corners««2
"moves a corner from index c to index x, including all references to it"
function move_corner!(m::CornerTable, c::Corner, x::Corner)
	println("  move_corner!($c =>  $x)/$(ncorners(m))")
	# replace all references to corner c by x
	v = apex(m, c)
	for k in fans(m, v)
		fanstart(m, k) == fan_closed(c) && fanstart!(m, k, fan_closed(x))
		fanstart(m, k) == fan_open(c) && fanstart!(m, k, fan_closed(x))
	end

	co = opposite(m, c)
	if issimple(co)
		opposite!(m, co, x)
	elseif ismultiple(co)
		for r in radial(m, Multi(co))
# 			println("in radial($co): $r")
			opposite(m, r) == Multi(c) || continue
			opposite!(m, r, Multi(x))
			break
		end
	end # if isboundary: nothing to do
	opposite!(m, x, opposite(m, c))
	fan!(m, x, fan(m, c))
end
# mesh modification (high level)««1
# edges««2
function match_edge!(m::CornerTable, c, cin, cout)#««
	println("\e[32;7m match_edge($c $(base(m,c))) cout=$cout cin=$cin\e[m")
	@assert isvalid(m)
# 	Int(cout) > 0 && println(" cout: $(apex(m, cout)) $(base(m, cout))")
# 	Int(cin) > 0 && println(" cin: $(apex(m, cin)) $(base(m, cin))")
# 	verbose(m)

	if !isboundary(cin) # there is already an inner corner
		op_in = opposite(m, cin)
		if isboundary(op_in) # cin is the last edge in its fan
			@assert c ≠ cin
# 			println("\e[31;7mmake new multiple 2-edge $([c,cin])\e[m")
			multi!(m, (c, cin))
		elseif issimple(op_in) # we break a simple edge
# 		println("\e[1mmake new multiple 3-edge $([c,cin,cout]) op_in = $op_in $(base(m,c))\e[m")
# 		println(base(m, c))
# 		println(base(m, cin))
# 		println(base(m, op_in))
# 			if op_in ≠ cout
# 				verbose(m)
# 			end
			@assert op_in == cout
			multi!(m, (c, cin, cout))
			ni = next(cin); no = next(cout)
			split_fan!(m, fan(m, no), no)
			split_fan!(m, fan(m, ni), ni)
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
	@assert isvalid(m)
	println("\e[32;3m match_edge($c) done\e[m")
	end
end#»»
function cut_edge!(m::CornerTable, c::Corner, t = nothing)#««
	println("\e[35;7m  cut_edge $c = $(base(m,c))\e[m"); # verbose(m)
	@assert isvalid(m)
	# 1. cut fans if needed
	(n, p) = next2(c)
	@assert Int(c) > 0
	op = opposite(m, c)
	if isboundary(op) # nothing to do
	elseif issimple(op)
		# simple edge: we need to split fans at both vertices v2 and v3
		println("\e[1m  cutting *simple edge* $c--$op = $(base(m,c))\e[m")
		k2 = fan(m, n)
		k3 = fan(m, p)
		split_fan!(m, k2, n)
		split_fan!(m, k3, next(op))
		println("\e[33;7mafter split_fans:\e[m")
		opposite!(m, op, Corner(0))
		opposite!(m, c, Corner(0))
		# delete a simple edge: replace by two boundaries
	else # multiple edge; removing one edge may change it to boundary or simple
		# to simple: only if there are three edges (c, c1, c2) in the loop...
		c1 = Multi(op)
		c2 = Multi(opposite(m, c1))
		println("\e[1mremoving multiple edge $c ($c1, $c2) $(base(m,c))\e[m")
		println(collect(genradial(m, c)))
		if c2 == c # create boundary edge
			@assert left(m, c1) == left(m, c)
			println("\e[1m  => create boundary edge $c1 => nothing\e[m")
			opposite!(m, c1, Corner(0))
		else # does this create a simple edge?
			c3 = Multi(opposite(m, c2))
			println("  c3=$c3")
			# ... and both remaining edges (c1 and c2) are opposed
			if c3 == c && left(m, c1) ≠ left(m, c2)
				println("\e[1m  => create simple edge opposite=$c1 sameside=$c2\e[m")
				println("  base($c1) = $(base(m,c1)) base($c2) = $(base(m,c2))")
				verbose(m, fan(m, c1))
				verbose(m, fan(m, next(c1)))
				verbose(m, fan(m, prev(c1)))
				verbose(m, fan(m, c2))
				verbose(m, fan(m, next(c2)))
				verbose(m, fan(m, prev(c2)))
				println("  nfans=$(nfans(m))")
				Int(c) == 452 && println("  fan(843)=$(fan(m,Corner(843)))")
				println("making $c1 and $c2 opposite...")
				opposite!(m, c1, c2)
				opposite!(m, c2, c1)
				glue_fans!(m, fan(m,prev(c2)), fan(m, next(c1)))
				Int(c) == 452 && println("  fan(843)=$(fan(m,Corner(843)))")
				glue_fans!(m, fan(m,prev(c1)), fan(m, next(c2)))
				Int(c) == 452 && println("  fan(843)=$(fan(m,Corner(843)))")
				println("now fan($(prev(c1))) = $(fan(m,prev(c1)))")
				println("now fan($(next(c2))) = $(fan(m,next(c2)))")
				println("now fan($(prev(c2))) = $(fan(m,prev(c2)))")
				println("now fan($(next(c1))) = $(fan(m,next(c1)))")
				Int(c) == 452 && println("  fan(843)=$(fan(m,Corner(843)))")
				println("  nfans=$(nfans(m))")
				verbose(m, c1)
				verbose(m, c2)
			else # otherwise, this edge remains multiple, we unlink c
			println("\e[1m => remains multiple...\e[m")
				for c2 in radial(m, c1)
					opposite(m, c2) == Multi(c) || continue
					opposite!(m, c2, op)
				end
			end
		end
	end
	opposite!(m, c, Corner(0))
	@assert isvalid(m)
	println("\e[35;3m cut_edge($c) done\e[m")
end#»»
# faces««2
"creates a new, detached face"
function new_face!(m::CornerTable)
	nfaces!(m, nfaces(m)+1)
	for c in corners(lastface(m))
		fan!(m, c, Fan(0))
		opposite!(m, c, Corner(0))
	end
	return lastface(m)
end
"creates a new face and connects it"
@inline append_face!(m::CornerTable, vlist, a=nothing) =
	attach_face!(m, new_face!(m), vlist, a)
"attaches all three corners and edges of a face"
function attach_face!(m::CornerTable{I}, f, vlist, a = nothing) where{I}
	attribute!(m, f, a)
	n = 3*Int(f)-3
# 	println("\e[32;7mattach_face!($f, $vlist)\e[m") #; verbose(m)
	v1 = vlist[1]; v2=vlist[2]; v3=vlist[3]
	c   = ( corner(f, Side(1)),  corner(f, Side(2)),  corner(f, Side(3)))
	cin = (findedge(m, v2, v3), findedge(m, v3, v1), findedge(m, v1, v2))
	cout= (findedge(m, v3, v2), findedge(m, v1, v3), findedge(m, v2, v1))
	new_fan!(m, vlist[1], Corner(n+1))
	new_fan!(m, vlist[2], Corner(n+2))
	new_fan!(m, vlist[3], Corner(n+3))
	for i in 1:3
		match_edge!(m, Corner(n+i), cin[i], cout[i])
	end
	return m
end
@inline isattachedface(m::CornerTable, f::Face) =
	!iszero(Int(fan(m, corner(f, Side(1)))))

"disconnects all three edges of a face"
function detach_face!(m::CornerTable, f::Face)
	println("\e[31;7mdetach_face!($f $(vertices(m,f)))\e[m")
	check1049(m)
# 	@assert isvalid(m)
# 	verbose(m)
	c1 = corner(f, Side(1)); v1 = apex(m, c1)
	c2 = corner(f, Side(2)); v2 = apex(m, c2)
	c3 = corner(f, Side(3)); v3 = apex(m, c3)
	println("fan(m, $c1)=$(fan(m,c1))")
	println("fan(m, $c2)=$(fan(m,c2))")
	println("fan(m, $c3)=$(fan(m,c3))")
	verbose(m, fan(m, c1))
	verbose(m, fan(m, c2))
	verbose(m, fan(m, c3))
	check1049(m)
	cut_edge!(m, c1)
	check1049(m)
	cut_edge!(m, c2)
	check1049(m)
	cut_edge!(m, c3)
	check1049(m)
# 	println("\e[31;1m after cut_edge:\e[m"); verbose(m)
	@assert opposite(m, c1) == Corner(0)
	@assert opposite(m, c2) == Corner(0)
	@assert opposite(m, c3) == Corner(0)
	# FIXME: don't delete fans, but remove one corner from the fan
	# (possibly deleting/splitting/moving the beginning).
	verbose(m, fan(m, c1))
	verbose(m, fan(m, c2))
	verbose(m, fan(m, c3))
	check1049(m)
	k = fan(m, c1); fan!(m, c1, Fan(0)); delete_fan!(m, k)
	k = fan(m, c2); fan!(m, c2, Fan(0)); delete_fan!(m, k)
	k = fan(m, c3); fan!(m, c3, Fan(0)); delete_fan!(m, k)
# 	delete_fan!(m, fan(m, c1)); fan!(m, c1, Fan(0))
# 	delete_fan!(m, fan(m, c2)); fan!(m, c2, Fan(0))
# 	delete_fan!(m, fan(m, c3)); fan!(m, c3, Fan(0))
	println("\e[1m face $f detached\e[m")
end
"deletes O(1) (detached) faces at once, by moving last faces in their place"
function delete_faces!(m::CornerTable, flist)#««
	println("\e[31;1;7mdelete_faces($flist)\e[m")
	check1049(m)
	isempty(flist) && return
	l = lastface(m)
	for (i, f) in pairs(flist)
		if isattachedface(m, l) # move last face `l` to position `f`
			a = attribute(m, l)
			vlist = vertices(m, l)
			detach_face!(m, l)
			attach_face!(m, f, vlist, a)
		end
		replace!(view(flist, i+1:length(flist)), l => f)
		l = Face(Int(l)-1)
	end
	nfaces!(m, Int(l))
end#»»
function replace_face!(m::CornerTable, f::Face, vlist, a = nothing)
	# don't move last face since we are immediately replacing:
	detach_face!(m, f)
	attach_face!(m, f, vlist, a)
end
# vertices««2
"merge_point!(m, v, x): replaces all references to `v` by `x`; deletes `v`"
function merge_point!(m::CornerTable, v::Vertex{I}, x::Vertex, l=nothing) where{I}#««
	check1049(m)
# 	println("\e[7mmerge_point($v => $x)\e[m")
# 	verbose(m)
	# we need to precompute the list of incident corners
	# collect() won't work, because Generators has a stupid haslength() type.
# 	clist = Corner{indextype(m)}[]
	clist = collect(c for k in fans(m, v) for c in star(m, k))
# 	println("  list of corners around $v=$clist")
	delfaces = Face{I}[]
# 	for (_,c) in vertexcorners(m,v); push!(clist, c); end
	for c in clist
		b = false
		f = face(c)
		println("  examining corner $c in face $f")
		check1049(m)
		if side(c) == Side(1)
			u1 = vertex(m, f, Side(2))
			u2 = vertex(m, f, Side(3))
			vlist = (x, u1, u2)
		elseif side(c) == Side(2)
			u1 = vertex(m, f, Side(3))
			u2 = vertex(m, f, Side(1))
			vlist = (u2, x, u1)
		else # side(c) == Side(3)
			u1 = vertex(m, f, Side(1))
			u2 = vertex(m, f, Side(2))
			vlist = (u1, u2, x)
		end
		checkface(m,f)
		check1049(m)
# 		println("  vertices = $(vertices(m,f))")
		detach_face!(m, f)
		check1049(m)
		if u1 == x || u2 == x
# 				println("    mark $f for deletion")
			push!(delfaces, f)
# 			delete_face!(m, f)
		else
# 			println("    set face $f to $vlist")
			check1049(m)
			attach_face!(m, f, vlist, attribute(m, f))
			check1049(m)
			checkface(m,f)
		end
	end
	delete_faces!(m, delfaces)
	delete_vertex!(m, v, l)
end#»»
# Constructor««2
function CornerTable{I,P,A}(points, faces, attributes) where{I,P,A}#««
	m = CornerTable{I,P,A}(points)
# 	println("faces = $faces")
	for ((v1,v2,v3), a) in zip(faces, attributes)
# 		verbose(m)
		append_face!(m, (Vertex(v1), Vertex(v2), Vertex(v3)), a)
	end
	return m
end#»»
@inline CornerTable{I,P}(points, faces,
		attributes = Iterators.repeated(nothing)) where{I,P} =
	CornerTable{I,P,eltype(attributes)}(points, faces, attributes)

@inline CornerTable{I}(points, faces, attributes...) where{I} =
	CornerTable{I,eltype(points)}(points, faces, attributes...)

@inline CornerTable(points, faces, attributes...) =
	CornerTable{Int}(points, faces, attributes...)
# coplanar and opposite faces««2
# returns the equivalence relation, as pairs of indices in `flist`:
function coplanar_faces(m::CornerTable, flist, ε = 0)
	@inline box(l,ε) = SpatialSorting.Box(l, l .+ ε)
	boxes = [ box(normalized_plane(m, f, absolute=true), ε) for f in flist ]
	return SpatialSorting.intersections(boxes)
end

function opposite_faces(m::CornerTable)
	r = NTuple{2,Face{indextype(m)}}[]
	for f in allfaces(m)
		c1 = corner(f, Side(1))
		ismultiple(opposite(m, c1)) || continue
		v2 = vertex(m, f, Side(2))
		v1 = vertex(m, f, Side(1))
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
# simplification (retriangulate faces) ««2
function checkface(m::CornerTable, f::Face)
	v =vertices(m,f)
	@assert v[1] ≠ v[2] && v[2] ≠ v[3] && v[3] ≠ v[1]
end
@inline function median(x1, x2, x3)
	if x1 ≤ x2
		x2 ≤ x3 && return x2
		x3 ≤ x1 && return x1
		return x3
	else # x2 < x1
		x1 ≤ x3 && return x1
		x3 ≤ x2 && return x2
		return x3
	end
end

function collapse_face!(m::CornerTable, f::Face)
	(v1, v2, v3) = vertices(m, f)
	p1 = point(m, v1); p2 = point(m, v2); p3 = point(m, v3)
	point!(m, v1, SVector{3}(median(p1[i], p2[i], p3[i]) for i in 1:3))
	refv3 = MVector(v3)
	merge_point!(m, v2, v1, refv3)
	verbose(m)
	merge_point!(m, refv3[1], v1)
end
"inserts vertex `s` in opposite edge"
function flatten_face!(m::CornerTable, c0::Corner)
	clist = collect(genradial(m, c0))
	v0 = apex(m, c0)
	f0 = face(c0)
	detach_face!(m, f0)
	for c in clist
		(c == c0) && continue # don't touch the original face (it is detached)
		(u1, u2, u3) = vertices(m, face(c))
		a = attribute(m, face(c))
		if side(c) == Side(1)
			(x1, x2, x3) = (u1, u2, v0)
			(y1, y2, y3) = (u1, v0, u3)
		elseif side(c) == Side(2)
			(x1, x2, x3) = (u2, u3, v0)
			(y1, y2, y3) = (u2, v0, u1)
		else # side(c) == Side(3)
			(x1, x2, x3) = (u3, u1, v0)
			(y1, y2, y3) = (u3, v0, u2)
		end
		detach_face!(m, face(c))
		attach_face!(m, face(c), (x1, x2, x3), a)
		append_face!(m, (y1, y2, y3), a)
	end
	delete_faces!(m, SA[f0])
end
function regularize!(m::CornerTable, ε = _DEFAULT_EPSILON)
	ε² = ε*ε
	ε3 = ε^(1/3)
	println("regularizing...")
	TI = TriangleIntersections
	# remove all triangular faces with area less than ε
	for f in allfaces(m)
# 		checkface(m, f)
		d = TI.degeneracy(triangle(m, f), ε3)
		d == TI.Constants.invalid && continue
		println("  regularizing degenerate face $f=$(vertices(m,f)): $d")
		if TI.isedge(d)
			println("    two merged points in position $(TI.index(d,TI.isedge))")
			c = corner(f, Side(TI.index(d, TI.isedge)))
			println("  merging vertices $(apex(m, prev(c))) and $(apex(m, next(c)))")
			merge_point!(m, apex(m, prev(c)), apex(m, next(c)))
# 			edge_collapse!(m, corner(f, Side(TI.index(d, TI.isedge))))
		elseif TI.isvertex(d)
			println("    vertex on edge $(TI.index(d,TI.isvertex))")
			s = Side(TI.index(d, TI.isvertex))
			flatten_face!(m, corner(f, s))
		else # interior: all three points confounded
			println("   3 confounded points, merging")
			verbose(m, f)
			(v1, v2, v3) = vertices(m, f)
			refv3 = MVector(v3)
			merge_point!(m, v2, v1, refv3)
			merge_point!(m, refv3[1], v1)
			# TODO: move to barycenter? or median coordinates?
# 			verbose(m, apex(m, corner(f, Side(1))))
# 			verbose(m, apex(m, corner(f, Side(2))))
# 			verbose(m, apex(m, corner(f, Side(3))))
# 			face_collapse!(m, f)
		end
	end
end
function simplify!(m::CornerTable, ε = _DEFAULT_EPSILON)
	for cluster in coplanar_clusters(m, ε)
		direction = main_axis(m, first(cluster))
		vlist = Vertex{indextype(m)}[]
		# TODO: detect cluster border and triangulate *that* only
		# (this would work with non-connected cluster as well)
		for f in cluster; push!(vlist, vertices(m, f)...); end
		tri = collect(project_and_triangulate(m, abs(direction),  vlist))
# 		length(tri) == length(cluster) && continue
# 		println("$cluster: $(length(cluster)) => $(length(tri))\n  $cluster\n  $([vertices(m,f) for f in cluster])\n  $tri")
	end
end
# create new meshes««1
"""
    reverse(m::CornerTable)

Returns a new corner table in which all faces are reversed (inside-out).
"""
function Base.reverse(m::CornerTable)#««
	newc = @closure c -> c-1+(c%3) # (1,3,2) permutation
	newopp = @closure c -> Corner(sign(Int(c))*newc(abs(Int(c))))
	newcorner = @closure c -> (c > 0) ? newc(c) : c

	r = (typeof(m))(points(m))
	r.vertex.fan .= m.vertex.fan
	nfaces!(r, nfaces(m))
	r.attribute .= m.attribute
	for n in 3:3:ncorners(r)
		opposite!(r, Corner(n-2), newopp(opposite(m, Corner(n-2))))
		opposite!(r, Corner(n-1), newopp(opposite(m, Corner(n-0))))
		opposite!(r, Corner(n-0), newopp(opposite(m, Corner(n-1))))
		fan!(r, Corner(n-2), fan(m, Corner(n-2)))
		fan!(r, Corner(n-1), fan(m, Corner(n-0)))
		fan!(r, Corner(n-0), fan(m, Corner(n-1)))
	end
	nfans!(r, nfans(m))
	r.fan.next .= m.fan.next
	r.fan.start .= newcorner.(m.fan.start)
	r.fan.vertex .= m.fan.vertex
	return r
end#»»
"""
    concatenate(m::CornerTable...)

Creates a new corner table from all given objects (renaming vertices).
"""
function concatenate(mlist::CornerTable...)#««
	r = (typeof(first(mlist)))(vcat(points.(mlist)...))
# 	println("\e[34;7mconcatenate: $(nvertices.(mlist)) $(nfaces.(mlist))\e[m")
	nfaces!(r, sum(nfaces(m) for m in mlist))
	nfans!(r, sum(nfans(m) for m in mlist))
	voffset = 0
	coffset = 0
	koffset = 0
	@inline shift(x,y) = sign(x)*(abs(x)+y)
	for m in mlist
		(nv, nc, nk) = (nvertices(m), ncorners(m), nfans(m))
		for i in 1:nv
			v0 = getv(m, Vertex(i))
			v1 = getv(r, Vertex(i+voffset))
			v1.point = v0.point
			v1.fan = v0.fan + koffset
		end
		for i in 1:nc
			c0 = getc(m, Corner(i))
			c1 = getc(r, Corner(i+coffset))
			c1.opposite = shift(c0.opposite, coffset)
			c1.fan = c0.fan + koffset
		end
		for i in 1:nk
			k0 = getk(m, Fan(i))
			k1 = getk(r, Fan(i+koffset))
			k1.next = k0.next + koffset
			k1.start = shift(k0.start, coffset)
			k1.vertex = k0.vertex + voffset
		end
		voffset+= nv
		coffset+= nc
		koffset+= nk
	end
	return r
end#»»
# self-intersection««1
function identify_point(m::CornerTable, f, t, p, idx, ε)
	!iszero(Int(idx)) && return idx
	TI = TriangleIntersections
	TI.isvertex(t) && return vertex(m, f, Side(TI.index(t, TI.isvertex)))
	for v in vertices(m, f)
		norm(point(m, v) - p, Inf) ≤ ε && return v
	end
	return idx
end
"""
    self_intersect(s)

Returns `(points, in_face, in_edge, faces)` describing the
self-intersection graph of `s`.
"""
function self_intersect(m::CornerTable{I}, ε=0) where{I}#««
	TI = TriangleIntersections
	regularize!(m)
	boxes = [ boundingbox(t...) for t in triangles(m) ]
	si = (points=similar(points(m), 0),
		in_face = SortedDict{Face{I},Vector{Vertex{I}}}(),
		faces = Face{I}[],
		edges = NTuple{2,Vertex{I}}[],
		)

	for (i1, i2) in SpatialSorting.intersections(boxes)
		f1 = Face(i1); f2 = Face(i2)
		isadjacent(m, f1, f2) && continue
# 		println("inter: ($f1 = $(vertices(m, f1)) $(triangle(m,f1))\n      ($f2=$(vertices(m,f2)) $(triangle(m,f2)))")
		it = TI.inter(triangle(m, f1), triangle(m, f2), ε)
		isempty(it) && continue
		
		# create points as needed, and store their index:
		idx = MVector{6,Vertex{I}}(undef)
		for (i, (p, (t1, t2))) in pairs(it)
			idx[i] = Vertex(0)
			idx[i] = identify_point(m, f1, t1, p, idx[i], ε)
			idx[i] = identify_point(m, f2, t2, p, idx[i], ε)
			if iszero(Int(idx[i]))
				push!(si.points, p)
				@debug "««creating point $p at $f1 ∩ $f2"
				for s in Side, f in (f1,f2)
					@debug "    distance to $(vertex(m,f,s)) = $(norm(p-point(m,vertex(m,f,s)),1))"
					(norm(p-point(m,vertex(m,f,s)),1)) < 1e-3 && @debug "$p, $(point(m,vertex(m,f,s)))"

				end
				@debug "»»"
				idx[i] = Vertex(length(si.points) + nvertices(m))
			end
			push_entry!(si.in_face, f1, idx[i])
			push_entry!(si.in_face, f2, idx[i])
		end
		# if length is >2 then the triangles are coplanar and this will be
		# detecting when regularizing the segment graph:
		length(it) == 2 && push!(si.edges, (idx[1], idx[2]))
		push!(si.faces, f1, f2)
	end
	uniquesort!(si.faces)
	return si
end#»»
# select_faces««2
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
	@debug  "selecting $(length(fkept)) faces, $(length(vkept)) vertices"

	# FIXME do something about fans!
	r = (typeof(m))(points(m)[Int.(vkept)])
	for f in fkept
		(v1, v2, v3) = vertices(m, f)
		u1 = Vertex(searchsortedfirst(vkept, v1))
		u2 = Vertex(searchsortedfirst(vkept, v2))
		u3 = Vertex(searchsortedfirst(vkept, v3))
		append_face!(r, (u1, u2, u3), attribute(m, f))
	end
	return r
end

# reverse, concatenate««2
function Base.reverse(m::CornerTable)#««
	newc = @closure c -> c-1+(c%3) # (1,3,2) permutation
	newopp = @closure c -> Corner(sign(Int(c))*newc(abs(Int(c))))
	newcorner = @closure c -> (c > 0) ? newc(c) : c

	r = (typeof(m))(points(m))
	r.vertex.fan .= m.vertex.fan
	nfaces!(r, nfaces(m))
	r.attribute .= m.attribute
	for n in 3:3:ncorners(r)
		opposite!(r, Corner(n-2), newopp(opposite(m, Corner(n-2))))
		opposite!(r, Corner(n-1), newopp(opposite(m, Corner(n-0))))
		opposite!(r, Corner(n-0), newopp(opposite(m, Corner(n-1))))
		fan!(r, Corner(n-2), fan(m, Corner(n-2)))
		fan!(r, Corner(n-1), fan(m, Corner(n-0)))
		fan!(r, Corner(n-0), fan(m, Corner(n-1)))
	end
	nfans!(r, nfans(m))
	r.fan.next .= m.fan.next
	r.fan.start .= newcorner.(m.fan.start)
	r.fan.vertex .= m.fan.vertex
	return r
end#»»
function concatenate(mlist::CornerTable...)#««
	r = (typeof(first(mlist)))(vcat(points.(mlist)...))
# 	println("\e[34;7mconcatenate: $(nvertices.(mlist)) $(nfaces.(mlist))\e[m")
	nfaces!(r, sum(nfaces(m) for m in mlist))
	nfans!(r, sum(nfans(m) for m in mlist))
	voffset = 0
	coffset = 0
	koffset = 0
	@inline shift(x,y) = sign(x)*(abs(x)+y)
	for m in mlist
		(nv, nc, nk) = (nvertices(m), ncorners(m), nfans(m))
		for i in 1:nv
			v0 = getv(m, Vertex(i))
			v1 = getv(r, Vertex(i+voffset))
			v1.point = v0.point
			v1.fan = v0.fan + koffset
		end
		for i in 1:nc
			c0 = getc(m, Corner(i))
			c1 = getc(r, Corner(i+coffset))
			c1.opposite = shift(c0.opposite, coffset)
			c1.fan = c0.fan + koffset
		end
		for i in 1:nk
			k0 = getk(m, Fan(i))
			k1 = getk(r, Fan(i+koffset))
			k1.next = k0.next + koffset
			k1.start = shift(k0.start, coffset)
			k1.vertex = k0.vertex + voffset
		end
		voffset+= nv
		coffset+= nc
		koffset+= nk
	end
	return r
end#»»
# subtriangulation««1
# project_and_triangulate ««2
function project_and_triangulate(m::CornerTable, proj, vlist,elist, ε = 0)
	plist = [ proj(point(m, v)) for v in vlist ]
	elist2 = map(e->map(v->Int(searchsortedfirst(vlist, v)), e), elist)

# 	global K=(plist, elist2, ε); @assert false
	println("plist = $plist\nelist=$elist2\n")
	SegmentGraphs.simplify!(plist, elist2, ε)
	newpoints = inv(proj).(plist[length(vlist)+1:end])
	vlist = [vlist; Vertex.(nvertices(m)+1:nvertices(m)+length(plist)-length(vlist))]
	append_points!(m, newpoints)

	# convert to format used by constrained_triangulation
	vmat = Float64[ p[i] for p in plist, i in 1:2 ]
	emat = [ e[i] for e in elist2, i in 1:2 ]
	io = open("/tmp/a", "w")
	for (i, v) in pairs(vlist)
		(x,y) = vmat[i,:]
		println(io,"$x\t$y\t$v\t# ($i)")
	end
	println(io,"\n\n")
	for (i1, i2) in eachrow(emat)
		(x1,y1) = vmat[i1,:]; (x2,y2) = vmat[i2,:]
		println(io,"$x1\t$y1\t$(x2-x1)\t$(y2-y1) # $(vlist[i1])--$(vlist[i2])")
	end
	println(io, "# using Triangle; constrained_triangulation($vmat, [1:$(length(plist))...], $emat)")
	close(io)
	tri = LibTriangle.constrained_triangulation(vmat, [1:length(plist)...], emat)
	return [ (vlist[a], vlist[b], vlist[c]) for (a,b,c) in tri ]
end

# subtriangulate««2
function subtriangulate!(m::CornerTable{I}, ε=0) where{I}
	si = self_intersect(m, ε)
	# first renumber points, removing duplicates, including in self-intersect:
	append_points!(m, si.points)
	vmap = simplify_points!(m, ε)
# 	explain(m, "/tmp/x.scad", scale=30)

	for i in eachindex(si.in_face)
		si.in_face[i] = map(x->get(vmap, x, x), si.in_face[i])
	end
	for i in eachindex(si.edges)
		si.edges[i] = extrema(map(x->get(vmap, x, x), si.edges[i]))
	end
	uniquesort!(si.edges)

	# determine clusters of coplanar (broken) faces:
	coplanar_rel = coplanar_faces(m, si.faces, ε)
	clusters = equivalence_structure(length(si.faces), coplanar_rel)
	
	# compute points by face
	in_face_v = [ get(si.in_face, f, Vertex{I}[]) for f in si.faces ]
	in_face_e = [ NTuple{2,Vertex{I}}[] for f in si.faces ]
	for (f, vlist, elist) in zip(si.faces, in_face_v, in_face_e)
		ff = vertices(m, f)
# 		println("examining $f=$ff; in_face=$vlist, in_face_e=$elist")
		push!(vlist, ff...)
		for s in Side
			push!(elist, base(m, corner(f, s)))
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
# 			println("\e[1mcluster $icluster\e[m\n  $allvertices\n  $alledges")
# 			for i in icluster; f = si.faces[i]; println("  $f=$(vertices(m,f)) ∋ $(in_face_v[i])"); end
# 			println("  tri= $alltriangles")
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
				a[1]*b[2]-a[2]*b[1] ≤ ε || continue
				b[1]*c[2]-b[2]*c[1] ≤ ε || continue
				c[1]*a[2]-c[2]*a[1] ≤ ε || continue
			end
			push!(vlist[j], v)
		end
		# apply face refinement:
		for (k, i) in pairs(icluster)
			f = si.faces[i]
			a = attribute(m, f)
			orient = main_axis(m,f) == main_axis(proj)
			isfirst = true
# 			println("  \e[1mface $f = $(vertices(m, f)):\e[m $(vlist[k])")
			for tri in alltriangles
# 				println("    checking $tri")
				issubset(tri, vlist[k]) || continue
# 				println("    is subset!")
# 				println("    \e[31;1m $tri\e[m")
				orient || (tri = reverse(tri))
				(p,q,r) = (point(m,tri[1]), point(m,tri[2]), point(m,tri[3]))
				if norm(cross(q-p, r-p), Inf) ≤ ε
# 					println("\e[35;1m $(norm(cross(q-p,r-p),Inf))\e[m")
				end
				if isfirst
					detach_face!(m, f); isfirst = false
				else
					f = new_face!(m)
				end
				attach_face!(m, f, tri, a)
			end
		end
	end
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
# 						println("radial($c0): $c1")
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
# 		println("# $c = $(apex(m,c)) $(base(m, c)): fv3=$y, fv2=$z")
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
# 		println("  $c=$(apex(m,c)) $(base(m,c)): $u")
# 	end
	pt3 == nothing && return clist[reorder]

	# find where `pt3 - v1` inserts in radial loop:
	vec3 = pt3 - p1
	vec2 = project2d(axis, vec3) .- dot(vec3, dir3) .*dir2scaled
	println("searching for pt3=$pt3; vec3=$vec3, vec2=$vec2")
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
	fc = corners(m, i)
	# in triangle uvw, c is vertex u (edge vw), we want edge uv (corner w)
	(c, state) = iterate(fc); c = prev(c); cj = left(m, c)
	vpj = point(m, cj) - p
	xj = dot(vpi, vpj)
	iszero(xj) && return c
	xyj = (xj*xj, norm²(cross(vpi, vpj)))
	best = c
	while true
		u = iterate(fc, state)
		u == nothing && return best
		(c, state) = u; c = prev(c); ck = left(m, c)
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
		println("locate_point($(vmax[i1]))")
		println("  point=$(point(m, vmax[i1]))")
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
function isvalid(m::CornerTable)
	global J=m
	for v in allvertices(m)
		k = anyfan(m, v)
		if isisolated(k)
			for k in allfans(m)
				if vertex(m, k) == v
					println("vertex($k) = $v, which is isolated")
					return false
				end
			end
		else
			v2 = vertex(m, k)
			if v2 ≠ v
				println("vertex($k) = $v2, should be $v")
				return false
			end
		end
	end
	for c in allcorners(m)
		k = Int(fan(m,c))
		if (k < 0 || k > nfans(m))
			println(" fan($c) is invalid: $k (should be 0:$(nfans(m)))")
			return false
		end
		cop = opposite(m, c)
		if issimple(cop)
			cop2 = opposite(m, cop)
			# if c is simple then fan(n) = fan(po) and fan(p) = fan(no)
			if cop2 ≠ c
				println("opposite²($c) = $cop2, should be $c")
				return false
			end
			if fan(m, next(c)) ≠ fan(m, prev(cop))
				println("edge ($c ↔ $cop) is simple but fan(next($c)) = $(fan(m, next(c))) ≠ $(fan(m, prev(cop))) = fan(m, prev($cop))")
				return false
			end
		else
			# multiple edge
		end
	end
	allstars = [ Int[] for _ in allfans(m) ]
	for c in allcorners(m)
		k = Int(fan(m,c))
		k > 0 && push!(allstars[k], Int(c))
	end
	for k in allfans(m)
		star_k = Int.(collect(star(m,k)))
		issubset(allstars[Int(k)], star_k) && continue
		println("\e[31;7m error for fan $k:\n  star=$star_k\n  but fan found in corners $(allstars[Int(k)])\e[m")
		return false
	end
	for k in allfans(m)
		v1 = vertex(m, k)
		v2 = vertex(m, nextfan(m, k))
		if v1 ≠ v2
			println("Fans $k and next=$(nextfan(m,k)) have different vertices $v1 ≠ $v2")
			return false
		end
		for c in star(m, k)
			kc = fan(m, c)
			if kc ≠ k
				println("fan($c) = $kc, should be $k")
				return false
			end
		end
	end
	return true
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
	for c in star(m, k)
		print(c, " (")
		kr = fan(m, next(c))
		print(Int(kr) == 0 ? Vertex(0) : vertex(m, kr), ",")
		kl = fan(m, prev(c))
		print(Int(kl) == 0 ? Vertex(0) : vertex(m, kl), ") ")
	end
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
			print(")")
		end
	end
	println()
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
	print("  \e[33;1m$c\e[m: \e[32m", fan(m,c), "\e[m")
	global M=m
	if isisolated(fan(m, c))
		print("  (isolated corner)")
	else
		print(" \e[34;7m", apex(m, c), "\e[m")
	end
	o = Int(opposite(m, c))
	print(" opposite=", (o > 0 ? "\e[32m" : o < 0 ? "\e[31m" : "\e[38;5;8m"),
		opposite(m, c), "\e[m")
	println()
end
function verbose(m::CornerTable, f::Face)
	print("\e[32;1m$f\e[m: ")
	vl = [isisolated(fan(m,c)) ? Vertex(0) : apex(m,c) for c in corners(f)]
	print(join(vl, ","))
	println()
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
# 	@assert isvalid(m)
end

export CornerTable

end # module
