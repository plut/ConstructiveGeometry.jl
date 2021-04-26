# TODO:
#  - locate_point
#  - reverse all faces of a mesh
#  - simplify (retriangulate) faces
#  - add per-face attributes
#  - maybe precompute all face normals
#  - clarify API
#  - make level computation output-sensitive
"""
    CornerTables

This module contains the basic function for operating with half-edge meshes
with triangular faces. For a general introduction to half-edge data structures,
see [Rhodes 2013](https://ubm-twvideo01.s3.amazonaws.com/o1/vault/gdc2013/slides/822512 Rhodes_Graham_Math_for_Game \\(2\\).pdf).
The structure used here is inspired by corner tables [Rossignac 2001].

This module exports the `CornerTable` data type and defines the following functions (not exported):
 - `faces`, `points`: exports a mesh to a list of triangles.
 - `concatenate`: disjoint union of meshes.
 - `combine`: compute the result of a Boolean operation (intersection, union)
   on meshes.
 - `Base.reverse`: computes the complement of a mesh; together with `combine`,
   this allows computation of a Boolean difference of meshes.
 - `Base.union`, `Base.intersect`, `Base.setdiff`: shorthand for usual combinations.
"""
module CornerTables
using StaticArrays
using LinearAlgebra
using FastClosures
using DataStructures
module LibTriangle
	using Triangle
end
include("TriangleIntersections.jl")
include("SpatialSorting.jl")
include("EquivalenceStructures.jl")
using .EquivalenceStructures

const _INDEX_TYPE = Int

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
@inline replace0!(::Nothing, args...) = nothing
@inline replace0!(t, args...) = replace!(t, args...)
# BBox««2
struct BBox{P}
	min::P
	max::P
	@inline BBox{P}(min, max) where{P} = new{P}(min, max)
	@inline BBox(min, max) = BBox{typeof(min)}(min, max)
end

@inline Base.eltype(::Type{BBox{P}}) where{P} = P
@inline Base.eltype(b::BBox) = eltype(typeof(b))

@inline Base.min(b::BBox) = b.min; @inline Base.max(b::BBox) = b.max
@inline Base.:∈(x::AbstractVector, b::BBox) = all(min(b) .≤ x .≤ max(b))
@inline Base.isempty(b::BBox) = any(min(b) .> max(b))
@inline Base.intersect(a::BBox, b::BBox) =
	BBox((max.(min(a), min(b))), (min.(max(a), max(b))))
@inline boundingbox(v::AbstractVector...) =
	BBox{eltype(v)}(min.(v...), max.(v...))
# `SpatialSorting` interface:
SpatialSorting.position(b::BBox) = b.min + b.max
SpatialSorting.merge(b1::BBox, b2::BBox) =
	BBox{eltype(b1)}(min.(b1.min, b2.min), max.(b1.max, b2.max))

"""
    simplify_points(points, ε)

Removes duplicates from the set of points, returning the set of
substitutions (as an (oldindex => newindex) iterator).
"""
function simplify_points(points, ε=0)
	n = length(points)
	boxes = [ BBox(p, p .+ ε) for p in points ]
	# `extrema` guarantees that all pairs (i,j) are sorted i < j
	samepoints = extrema.(SpatialSorting.intersections(boxes))
	eqv = equivalence_structure(n, samepoints)
	return (Vertex(i) => Vertex(j)
		for (i, j) in pairs(representatives(eqv)) if i ≠ j)
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
function main_axis(direction)#««
	@inbounds (a1, a2, a3) = abs.(direction)
	if a2 < a1
		a1 < a3 && @goto max3
		return direction[1] > 0 ? 1 : -1
	elseif a2 > a3
		return direction[2] > 0 ? 2 : -2
	end; @label max3
		return direction[3] > 0 ? 3 : -3
end#»»
function project2d(axis, vec)#««
	# orientation-preserving projection:
	axis == 1 && return (vec[2], vec[3])
	axis == 2 && return (vec[3], vec[1])
	axis == 3 && return (vec[1], vec[2])
	axis ==-1 && return (vec[3], vec[2])
	axis ==-2 && return (vec[1], vec[3])
	@assert axis == -3
	             return (vec[2], vec[1])
end#»»
function project1d(axis, vec)#««
	axis == 1 && return vec[1]
	axis == 2 && return vec[2]
	axis == 3 && return vec[3]
	axis ==-1 && return -vec[1]
	axis ==-2 && return -vec[2]
	@assert axis == -3
	             return -vec[3]
end#»»
# half-edge mesh««1
# half-edge data structure ««2
struct CornerTable{I<:Signed,P} # I is index type (integer), P is point type
#= Contents of the tables:
  points[vertex] = geometry for this point
	corner[vertex] = +(corner index)    if regular point (one implicit closed fan)
	                 0                  if isolated point (no triangle),
	                 -(first fan index) if singular point
	opposite[corner] = +(opposite corner)            if simple edge
	                   0                             if boundary edge
	                   -(next corner in radial loop) if multiple edge
	vertex[corner] = vertex index for this corner
	fan_next[fan]  = next fan around this vertex
	fan_first[fan] = +(first corner) if closed fan
	                 -(first corner) if open fan
=#
	points::Vector{P}
	corner::Vector{I}
	vertex::Vector{I}
	opposite::Vector{I}
	fan_next::Vector{I}
	fan_first::Vector{I}
	@inline CornerTable(points, opp, dest, ef, cn=[], cs=[]) = # TEMPORARY
		new{Int,eltype(points)}(points, opp, dest, ef, cn, cs)
	@inline CornerTable{I,P}(points::AbstractVector) where{I,P} =
		new{I,P}(points, zeros(I, length(points)), [], [], [], [])
	@inline CornerTable{I,P}() where{I,P} = CornerTable{I,P}(P[])
	@inline CornerTable{I}(points::AbstractVector) where{I} =
		CornerTable{I,eltype(points)}(points)
	@inline CornerTable(points::AbstractVector) = CornerTable{Int}(points)
end

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

@inline index_type(::Type{<:CornerTable{I}}) where{I} = I
@inline index_type(m::CornerTable{I}) where{I} = I
# several names for the same type, for clarity:
@inline NamedIndex{S}(m::CornerTable) where{S} = NamedIndex{S,index_type(m)}
point_type(::Type{<:CornerTable{I,P}}) where{I,P} = P
point_type(m::CornerTable) = point_type(typeof(m))

# size and resizing functions««2
# vertices««3
@inline nvertices(m::CornerTable) = length(m.points)
@inline lastvertex(m::CornerTable{I}) where{I} = Vertex{I}(nvertices(m))
@inline function nvertices!(m::CornerTable, nv::Integer)
	resize!(m.points, nv)
	resize!(m.corner, nv)
	return m
end
@inline function verticeshint!(m::CornerTable, nv::Integer)
	sizehint!(m.points, nv)
	sizehint!(m.corner, nv)
	return m
end
@inline allvertices(m::CornerTable) = (Vertex(v) for v in 1:nvertices(m))
# faces««3
@inline nfaces(m::CornerTable) = length(m.opposite) ÷ 3
@inline lastface(m::CornerTable{I}) where{I} = Face{I}(nfaces(m))
@inline function nfaces!(m::CornerTable, nf::Integer)
	resize!(m.opposite, 3nf)
	resize!(m.vertex, 3nf)
	return m
end
@inline function faceshint!(m::CornerTable, nf::Integer)
	sizehint!(m.opposite, 3nf)
	sizehint!(m.vertex, 3nf)
	return m
end

@inline ncorners(m::CornerTable) = length(m.opposite)

@inline allcorners(m::CornerTable) = (Corner(c) for c in 1:length(m.opposite))
@inline allfaces(m::CornerTable) = (Face(f) for f in 1:nfaces(m))
# fans««3
@inline nfans(m::CornerTable) = length(m.fan_first)
@inline lastfan(m::CornerTable{I}) where{I} = Fan{I}(nfans(m))
@inline function nfans!(m::CornerTable, nk::Integer)
	resize!(m.fan_first, nk)
	resize!(m.fan_next, nk)
	return m
end
@inline function fanshint!(m::CornerTable, nk::Integer)
	sizehint!(m.fan_first, nk)
	sizehint!(m.fan_next, nk)
	return m
end
@inline allfans(m::CornerTable) = (Fan(k) for k in 1:nfans(m))

# simple accessors ««2
# corners ««3
@inline next(c::Corner) = Corner(Int(c)+1 - 3*(Int(c)%3==0))
@inline prev(c::Corner) = Corner(Int(c)-1 + 3*(Int(c)%3==1))
@inline opposite(m::CornerTable, c::Corner) = Corner(m.opposite[Int(c)])
@inline opposite!(m::CornerTable, c::Corner, x::Corner) =
	m.opposite[Int(c)] = Int(x)
@inline vertex(m::CornerTable, c::Corner) = Vertex(m.vertex[Int(c)])
@inline vertex!(m::CornerTable, c::Corner, v::Vertex) =
	m.vertex[Int(c)] = Int(v)
@inline right(m::CornerTable, c::Corner) = vertex(m, next(c))
@inline left(m::CornerTable, c::Corner) = vertex(m, prev(c))
@inline base(m::CornerTable, c::Corner) = (right(m, c), left(m, c))
@inline after(m::CornerTable, c::Corner) = next(opposite(m, next(c)))
@inline after!(m::CornerTable, c::Corner, x::Corner) =
	opposite!(m, next(c), prev(x))
@inline before(m::CornerTable, c::Corner) = prev(opposite(m, prev(c)))

# properties of opposite(c)
@inline issimple(c::Corner) = Int(c) > 0
@inline isboundary(c::Corner) = Int(c) == 0
@inline ismultiple(c::Corner) = Int(c) < 0

# faces ««3
@inline face(c::Corner) = Face(fld1(Int(c), 3))
@inline corner(f::Face, s::Side) = Corner(3*Int(f)-3+Int(s))
@inline corners(f::Face) =
	(corner(f,Side(1)), corner(f,Side(2)), corner(f,Side(3)))
@inline vertex(m::CornerTable, f::Face, s::Side) = vertex(m, corner(f,s))
@inline point(m::CornerTable, f::Face, s::Side) = point(m, vertex(m, f, s))
@inline vertices(m::CornerTable, f::Face) =
	(vertex(m,f,Side(1)), vertex(m,f,Side(2)), vertex(m,f,Side(3)))
@inline edge(m::CornerTable, f::Face, s::Side) = base(m, corner(f, s))
@inline faces(m::CornerTable{I}) where{I} = reinterpret(NTuple{3,I}, m.vertex)
@inline facevertices(m::CornerTable) =
	(Face(f) => vertices(m,Face(f)) for f in 1:nfaces(m))
@inline function vertices!(m::CornerTable, f::Face, tri::NTuple{3,<:Vertex})
	vertex!(m, corner(f,Side(1)), tri[1])
	vertex!(m, corner(f,Side(2)), tri[2])
	vertex!(m, corner(f,Side(3)), tri[3])
	return m
end

@inline adjacent(m::CornerTable, f::Face, s::Side) =
	face(opposite(m, corner(f, s)))
@inline adjacent(m::CornerTable, f::Face) =
	(adjacent(m, f, Side(1)), adjacent(m, f, Side(2)), adjacent(m, f, Side(3)))
@inline isadjacent(m::CornerTable, f1::Face, f2::Face) =
	any(==(f2), adjacent(m, f1))

# vertices ««3
@inline corner(m::CornerTable, v::Vertex) = Corner(m.corner[Int(v)])
@inline corner!(m::CornerTable, v::Vertex, c::Corner)= m.corner[Int(v)]=Int(c)
@inline points(m::CornerTable) = m.points
@inline point(m::CornerTable, v::Vertex) = m.points[Int(v)]
@inline point!(m::CornerTable, v::Vertex, p) = m.points[Int(v)] = p

@inline triangle(m::CornerTable, f::Face) =
	(point(m, f, Side(1)), point(m,f,Side(2)), point(m,f,Side(3)))
@inline triangles(m::CornerTable) = (triangle(m, Face(f)) for f in 1:nfaces(m))
@inline function normal(m::CornerTable, f::Face)
	t = triangle(m, f)
	return cross(t[2]-t[1], t[3]-t[2])
end
@inline main_axis(m::CornerTable, f::Face) = main_axis(normal(m, f))
function normalized_plane(m::CornerTable, f::Face)
	t = triangle(m, f)
	d = cross(t[2]-t[1], t[3]-t[2])
	a = main_axis(d)
	c = @inbounds inv(d[abs(a)])
	return SA[oftype(c, a), project2d(a, d) .* c..., dot(d, t[1])*c]
end
@inline volume(m::CornerTable) =
	sum(dot(u, cross(v, w)) for (u,v,w) in triangles(s))/6

# fans ««3
# properties of `corner(::Vertex)`:
@inline issingular(c::Corner) = Int(c) < 0
@inline isregular(c::Corner) = Int(c) > 0
@inline isisolated(c::Corner) = Int(c) == 0

@inline fan_next(m::CornerTable, k::Fan) = Fan(m.fan_next[Int(k)])
@inline fan_next!(m::CornerTable, k::Fan, x::Fan) = m.fan_next[Int(k)] = Int(x)

@inline Fan(c::Corner) = Fan(-Int(c))
@inline Corner(k::Fan) = Corner(-Int(k))
@inline isclosedfan(c::Int) = Int(c) > 0
@inline isopenfan(c::Int) = Int(c) < 0
@inline fan_open(c::Corner)   = -Int(c)
@inline fan_closed(c::Corner) = +Int(c)

const implicitfan = Fan(0)
@inline fan_first(m::CornerTable, k::Fan) = m.fan_first[Int(k)]
@inline fan_first!(m::CornerTable, k::Fan, x::Integer) = m.fan_first[Int(k)]= x
@inline fan_firstcorner(m::CornerTable, k::Fan) = Corner(abs(fan_first(m,k)))
@inline fanvertex(m::CornerTable, k::Fan) = vertex(m, fan_firstcorner(m, k))

function fan_prev(m::CornerTable, k1::Fan)
	k2 = k1
	while true
		u = fan_next(m, k2)
		u == k1 && return k2
		k2 = u
	end
end
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
global COUNT=0
@inline Base.iterate(it::MeshIterator) = (global COUNT=0; (it.start, it.start))
@inline function Base.iterate(it::MeshIterator, s)
	global M = it.mesh
	global COUNT+=1; @assert COUNT <= 10
	s = next(it, it.mesh, s)
	stop(it, it.mesh, s) && return nothing
	return (s, s)
end
@inline Base.eltype(::MeshIterator{S,I,T}) where{S,I,T} = T
@inline Base.IteratorSize(::MeshIterator) = Base.SizeUnknown()

@inline star(m::CornerTable, c::Corner) = MeshIterator{star}(m, c)
@inline star(m::CornerTable, v::Vertex) = MeshIterator{star}(m, corner(m, v))
@inline next(it::MeshIterator{star}, m, c) = after(m, c)
@inline stop(it::MeshIterator{star}, m, c) = c == it.start

@inline fans(m::CornerTable, k::Fan) = MeshIterator{fans}(m, k)
@inline next(it::MeshIterator{fans}, m, k) = fan_next(m, k)
@inline stop(it::MeshIterator{fans}, m, k) = k == it.start

@inline fancorners(m::CornerTable, k::Fan) =
	MeshIterator{fancorners}(m, Corner(abs(fan_first(m, k))))
@inline next(it::MeshIterator{fancorners}, m, c) = after(m, c)
@inline stop(it::MeshIterator{fancorners}, m, c) = (c==it.start)||!issimple(c)

# FIXME: give a way to iterate over a simple edge
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
	

# @inline closedfan(m::CornerTable, c::Corner) = MeshIterator{closedfan}(m, g)
# @inline next(it::MeshIterator{closedfan},m, c) = after(m, c)
# @inline stop(it::MeshIterator{closedfan},m, c) = c == it.start
# 
# @inline openfan(m::CornerTable, c::Corner) = MeshIterator{openfan}(m, g)
# @inline next(it::MeshIterator{openfan}, m, c) = after(m, c)
# @inline stop(it::MeshIterator{openfan}, m, c) = isboundary(c)

"finds in which fan a corner lies"
function fan(m::CornerTable, c0::Corner)
	for (k, c) in vertexcorners(m, vertex(m, c0))
		c == c0 && return k
	end
	return implicitfan
end
"iterator on (fans, corners) around a given vertex"
function vertexcorners(m::CornerTable, u::Vertex)
	c0 = corner(m, u)
	isregular(c0) && return ((implicitfan, c) for c in star(m, c0))
	isisolated(c0) && return ((implicitfan, c0) for _ in 1:0)
	return ((k,c) for k in fans(m, Fan(c0)) for c in fancorners(m, k))
end
# function vertexcorners(m::CornerTable, u::Vertex)
# 	c0 = corner(m, u)
# 	isregular(c0) && return star(m, c0)
# 	isisolated(c0) && return (c0 for _ in 1:0) # empty with Corner eltype
# 	# c0 is a fan
# 	return (c for k in fans(m, Fan(c0)) for c in fancorners(m, k))
# end


# mesh modification««1
# move_xxx(m, x, y): moves object from index x to index y, overwriting
#                    previous value at y
# delete_xxx(m, x):  deletes object at index x (by moving last onto it)
# vertices««2
function append_points!(m::CornerTable, plist)#««
	if Base.IteratorSize(plist) ≠ Base.SizeUnknown()
		verticeshint!(m, nvertices(m) + length(plist))
	end
	for p in plist
		push!(m.points, p)
		push!(m.corner, 0)
	end
	return m
end#»»
@inline move_vertex!(m::CornerTable, v::Vertex, x::Vertex) =
	for (_,c) in vertexcorners(m, v); vertex!(m, c, x); end
@inline move_point!(m::CornerTable, v::Vertex, x::Vertex) =
	(move_vertex!(m, v, x); corner!(m, x, corner(m, v)); point!(m, x, point(m,v)))
function merge_point!(m::CornerTable, v::Vertex, x::Vertex)
	# deletes vertex v, replacing all references by x
	move_vertex!(m, v, x)
	l = lastvertex(m)
	v ≠ l && move_point!(m, l, v)
	nvertices!(m, nvertices(m)-1)
end

"removes ε-duplicates from point of m, and returns list of substitutions
(as `oldindex=>newindex` iterator)"
function simplify_points!(m::CornerTable{I}, ε = 0) where{I}
	repl = simplify_points(points(m), ε)
	newname = SortedDict{Vertex{I},Vertex{I}}()
	for (a, b) in repl
		a1 = get(newname, a, a)
		b1 = get(newname, b, b)
		l = lastvertex(m) # call this before merge_point!
		merge_point!(m, a1, b1)

		(a == a1) ? push!(newname, a => b1) : (newname[a] = b1)
		l ≠ a1 && push!(newname, l => a1)
	end
	return newname
end
# fans««2
"creates the first fan at vertex `v`"
function create_fan!(m::CornerTable, v::Vertex, x::Integer)
# 	println("\e[34;1mcreate_fan!($v)\e[m")
	global M=deepcopy(m)
	@assert !issingular(corner(m, v)) # this must be the *first* fan
	push!(m.fan_first, x)
	push!(m.fan_next, nfans(m))
	return lastfan(m)
end
"appends a new fan to same vertex as `k`"
function append_fan!(m::CornerTable, k::Fan, x::Integer)
# 	println("\e[34;1mappend_fan!($k, $x)\e[m")
	push!(m.fan_first, x)
	push!(m.fan_next, m.fan_next[Int(k)])
	m.fan_next[Int(k)] = nfans(m)
	return lastfan(m)
end
"creates a new single-corner fan at this vertex; returns fan number"
function new_fan!(m::CornerTable, v::Vertex, c::Corner)
# 	println("\e[34;1mnew_fan!($v, $c)\e[m")
	cv = corner(m, v)
	if isisolated(cv)
		k = create_fan!(m, v, fan_open(c))
		corner!(m, v, Corner(k))
		return k
	elseif isregular(cv)
		k = create_fan!(m, v, fan_closed(cv))
		corner!(m, v, Corner(k))
		return append_fan!(m, k, fan_open(c))
	else # issingular(cv)
		return append_fan!(m, Fan(cv), fan_open(c))
	end
end
"if v has an implicit fan, make it explicit"
function explicit_fan!(m::CornerTable, v::Vertex)
	c = corner(m, v)
	isregular(c) || return
	k = create_fan!(m, v, fan_closed(c))
	corner!(m, v, Corner(k))
end
"if v has a single closed fan, make it implicit"
function implicit_fan!(m::CornerTable, v::Vertex)
	k = Fan(corner(m, v))
	@assert Int(k) > 0
	c = fan_first(m, k)
	isclosedfan(c) || return
	fan_next(m, k) == k || return
	corner!(m, v, Corner(c))
	delete_fan!(m, k)
end

"moves a fan from index k to index x"
function move_fan!(m::CornerTable, k::Fan, x::Fan, t = nothing)
# 	println("\e[34;1mmove_fan!($k, $x)\e[m")
	corner!(m, fanvertex(m, k), Corner(x))
	fan_next!(m, fan_prev(m, k), x)
	fan_next!(m, x, fan_next(m, k))
	fan_first!(m, x, fan_first(m, k))
	t ≠ nothing && replace!(t, k => x)
end
function delete_fan!(m::CornerTable, k::Fan, t = nothing)
# 	println("\e[34;1mdelete_fan!($k)\e[m")
	@assert k ≠ implicitfan
	n = fan_next(m, k)
	v = fanvertex(m, k)
	if n == k
		c = fan_first(m, k)
		if isclosedfan(c)
			# if this is the only fan for the point, make it a regular point:
			corner!(m, v, Corner(c))
		else
			# make it an isolated point:
			corner!(m, c, Corner(0))
		end
	else
		# remove from cyclic list
		corner!(m, v, Corner(n))
		fan_next!(m, fan_prev(m, k), n)
	end
	# replace all references to last fan with references to k
	l = lastfan(m)
	l ≠ k && move_fan!(m, l, k, t)
	nfans!(m, nfans(m)-1)
end
"connects two open fans together (k2 after k1)."
function glue_fans!(m::CornerTable, k1, k2, t = nothing)
# 	println("\e[34;1mglue_fans!($k1, $k2) $t\e[m")
	@assert fanvertex(m, k1) == fanvertex(m, k2)
	@assert isopenfan(fan_first(m, k1))
	@assert isopenfan(fan_first(m, k2))
	if k1 == k2 # glue an open fan with itself by making it closed
# 		if fan_next(m, k2) == k2 # a single closed fan is made implicit
# 			replace0!(t, k2 => implicitfan)
# 			delete_fan!(m, k2, t)
# 		else
			fan_first!(m, k2, -fan_first(m, k2))
# 		end
	else # remove fan k2 from linked list
		replace0!(t, k2 => k1)
		delete_fan!(m, k2, t)
	end
end

"splits fan `k` at a simple edge, making `c` its new first corner. Returns the pair of fans (before, after)"
@inline function split_fan!(m::CornerTable, k::Fan, c::Corner, t = nothing)
# 	println("  \e[34msplit_fan!($k, $c) $(vertex(m,c))\e[m")
	if k == implicitfan
		@assert false
# 		v = vertex(m, c)
# 		println("  splitting implicit fan at vertex $(vertex(m,c)); current corner($v)=$(corner(m,v))")
# 		k = create_fan!(m, v, fan_open(c))
# 		corner!(m, v, Corner(k))
# 		return (k, k)
	elseif isclosedfan(fan_first(m, k))
		fan_first!(m, k, fan_open(c))
		return (k, k)
	else
		return (k, append_fan!(m, k, fan_open(c)))
		# create a new fan
# 		println("splitting (open) fan $k at corner $c")
# 		verbose(m, k)
# 		println("vertex is $(vertex(m,c))")
# 		create_fan!(m, vertex(m, c), fan_open(c))
	end
end

# edges««2
# mesh queries««2
"returns (corner, fan at u, fan at v)"
function findedge(m::CornerTable, u, v, w, clist, klist, i)#««
	# rotates around vertex u
	# annotates in clist the following corners:
	# i  : a corner which sees (u, v)
	# i+3: a corner which sees (u, w)
	# i+6: same as next(i) (corner at u after v)
	# i+9: same as next(i+3) (corner at u after w)
	# and also the following fans:
	# i: fan before v
	# i+6: fan after w
	#
	klist[i] = klist[i+6] = implicitfan
	clist[i] = clist[i+3] = clist[i+6] = clist[i+9] = Corner(0)
	for (k, c) in vertexcorners(m, u)
		n = next(c); r = vertex(m, n)
		p = prev(c); l = vertex(m, p)
		(r == v) && (clist[i] = p; clist[i+6] = c)
		(l == v) && (klist[i] = k)
		(r == w) && (clist[i+3] = p; clist[i+9] = c; klist[i+6] = k)
	end
end#»»
# corners««2
"moves a corner from index c to index x (used by move_face!)"
function move_corner!(m::CornerTable, c::Corner, x::Corner)
# 	println("move_corner!($c, $x)")
	# replace all references to corner c by x
	v = vertex(m, c)
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
		for r in radial(m, co)
			opposite(m, r) == Multi(c) || continue
			opposite!(m, r, Multi(x))
			break
		end
	end # if isboundary: nothing to do
	opposite!(m, x, opposite(m, c))
	vertex!(m, x, vertex(m, c))
end
# edges««2
function edge!(m::CornerTable, klist, clist, i1)#««
	i2 = @inbounds (2,3,1)[i1]; i3=@inbounds (3,1,2)[i1]
	# corner id:
	# i2:    previous inner corner for our edge
	# i3+3:  previous outer corner
	# i1+6:  first inner corner if regular fan
	# i1+9:  first outer corner if regular fan
	# fan id:
	# i1:    fan before corner
	# i1+3:  fan inside corner
	# i1+6:  fan after corner
	c = clist[i1+12]; cout = clist[i3+3]; cin = clist[i2]
# 	println("\e[34;7mglue_edge($c $(base(m,c))) klist[5]=$(klist[5])\e[m")
# 	verbose(m)
# 	println("""
# \e[1m$c (i1=$i1) apex=$(vertex(m,c)) cin=$cin, out=$cout\e[m
#   right=$(right(m,c)) fans out=$(klist[i2])→in=$(klist[i2+3])
#   left=$(left(m,c)) fans in=$(klist[i3+3])→out=$(klist[i3+6])""")
	if !isboundary(cin) # there is already an inner corner
		op_in = opposite(m, cin)
		if isboundary(op_in) # op_in is the last edge in its fan
# 			println("\e[1mmake new multiple 2-edge $([c,cin])\e[m")
			multi!(m, [c, cin])
		elseif issimple(op_in) # we break a simple edge
# 		println("\e[1mmake new multiple 3-edge $([c,cin,cout]) $(base(m,c))\e[m")
# 		println("around $(left(m,c)): $(klist[[i3,i3+3,i3+6]]) around $(right(m,c)): $(klist[[i2,i2+3,i2+6]])")
			@assert op_in == cout
			multi!(m, [c, cin, cout])
# 			println("edge ($cin, $cout) is simple")
# 			verbose(m, klist[i2])
# 			verbose(m, klist[i3+3])
			(u1, u2) = split_fan!(m, klist[i2+3], clist[i2+6])
# 			println("result of split_fan!: ($u1, $u2)")
			(u3, u4) = split_fan!(m, klist[i3+6], clist[i3+9])
# 			println("\e[34;1mresults of split_fan!: $u1 $u2, $u3 $u4\e[m")
		else # already a multiple edge
# 			println("\e[1mappend $c to (inner) multiple edge $op_in\e[m")
			add_multi!(m, Multi(op_in), c)
		end
	elseif !isboundary(cout) # there is an outer corner, but no inner
		op_out = opposite(m, cout)
		@assert !issimple(op_out)
		if ismultiple(op_out)
# 			println("\e[1mappend $c to (outer) multiple edge $op_out\e[m")
			add_multi!(m, Multi(op_out), c)
		else # op_out is a boundary: we are creating a simple edge
# 			println("\e[1mcreate simple edge $c -- $cout ($cin)\e[m")
			opposite!(m, c, cout); opposite!(m, cout, c)
			glue_fans!(m, klist[i2], klist[i2+3], klist)
			glue_fans!(m, klist[i3+3], klist[i3+6], klist)
		end
	end
end#»»
function cut_edge!(m::CornerTable, c::Corner, t = nothing)#««
# 	println("\e[35;7mcut_edge $c = $(base(m,c))\e[m")
# 	verbose(m)
# 	println("  first corner for $(vertex(m,c)) = $(corner(m,vertex(m,c)))")
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
		opposite!(m, op, Corner(0))
		opposite!(m, c, Corner(0))
# 		verbose(m)
		# delete a simple edge: replace by two boundaries
	else # multiple edge; removing one edge may change it to boundary or simple
		# to simple: only if there are three edges (c, c1, c2) in the loop...
		c1 = Multi(op)
		c2 = Multi(opposite(m, c1))
# 		println("removing multiple edge $c ($c1, $c2) $(base(m,c))")
		if c2 == c # create boundary edge
			@assert left(m, c1) == left(m, c)
# 			println("  => create boundary edge $c1 => nothing ")
			opposite!(m, c1, Corner(0))
		else # does this create a simple edge?
			c3 = Multi(opposite(m, c2))
# 			println("  c3=$c3")
			# ... and both remaining edges (c1 and c2) are opposed
			if c3 == c && left(m, c1) ≠ left(m, c2)
# 				println("  => create simple edge $c1 $c2")
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
function set_face!(m::CornerTable{I}, nf, vlist) where{I}
	n = 3*Int(nf)-3
# 	println("\e[32;7mset_face!($nf, $vlist)\e[m")
# 	verbose(m)
	v1 = vlist[1]; v2=vlist[2]; v3=vlist[3]
	explicit_fan!(m, v1)
	explicit_fan!(m, v2)
	explicit_fan!(m, v3)
	# this table holds 12 previous corners + the 3 new we create
	# around vertex 1, old corners are 1,4,7,10, new is 13, etc.
	klist = MVector{9,Fan{I}}(undef)
	clist = MVector{15,Corner{I}}(undef)
	findedge(m, v1, v2, v3, clist, klist, 1)
	findedge(m, v2, v3, v1, clist, klist, 2)
	findedge(m, v3, v1, v2, clist, klist, 3)
	for i in (1,2,3)
		clist[i+12] = Corner(n+i)
		opposite!(m, clist[i+12], Corner(0))
		vertex!(m, clist[i+12], vlist[i])
		klist[i+3] = new_fan!(m, vlist[i], clist[i+12])
	end
	edge!(m, klist, clist, 1)
	edge!(m, klist, clist, 2)
	edge!(m, klist, clist, 3)
	implicit_fan!(m, v1)
	implicit_fan!(m, v2)
	implicit_fan!(m, v3)
	return m
end

function append_face!(m::CornerTable{I}, vlist) where{I}#««
	nfaces!(m, nfaces(m)+1)
	return set_face!(m, lastface(m), vlist)
end#»»
"disconnects all edges and corners of a face"
function cut_face!(m::CornerTable, f::Face)
# 	println("\e[31;7mcut_face!($f $(vertices(m,f)))\e[m")
# 	verbose(m)
	c1 = corner(f, Side(1)); v1 = vertex(m, c1); explicit_fan!(m, v1)
	c2 = corner(f, Side(2)); v2 = vertex(m, c2); explicit_fan!(m, v2)
	c3 = corner(f, Side(3)); v3 = vertex(m, c3); explicit_fan!(m, v3)
# 	klist = MVector(fan(m, c1), fan(m, c2), fan(m, c3))
	cut_edge!(m, c1); cut_edge!(m, c2); cut_edge!(m, c3)
	@assert opposite(m, c1) == Corner(0)
	@assert opposite(m, c2) == Corner(0)
	@assert opposite(m, c3) == Corner(0)
	delete_fan!(m, fan(m, c1))
	delete_fan!(m, fan(m, c2))
	delete_fan!(m, fan(m, c3))
	implicit_fan!(m, v1)
	implicit_fan!(m, v2)
	implicit_fan!(m, v3)
# 	println("\e[31;7mafter cut_face!($f $(vertices(m,f))):\e[m")
# 	verbose(m)
	global M=deepcopy(m)
# 	@assert f ∈ (Face(1),Face(5),Face(2),Face(6),)
end
"moves a face from index f to index x"
@inline move_face!(m::CornerTable, f::Face, x::Face) =
	for s in Side; move_corner!(m, corner(f, s), corner(x, s)); end
"deletes a face (disconnect + move last onto it)"
function delete_face!(m::CornerTable, f::Face)#««
	cut_face!(m, f)
	l = lastface(m)
	f ≠ l && move_face!(m, l, f)
	nfaces!(m, Int(l)-1)
end#»»
function replace_face!(m::CornerTable, f::Face, vlist)
	# don't move last face since we are immediately replacing:
	cut_face!(m, f)
	set_face!(m, f, vlist)
end
# reversing
# function radial_loops!(s::CornerTable, edge_loc)#««
# 	for (e, hlist1) in edge_loc
# 		e[1] > e[2] && continue # we do two this only once per edge
# 		isempty(hlist1) && continue
# 		hlist2 = edge_loc[reverse(e)]
# 		unique!(sort!(hlist1))
# 		unique!(sort!(hlist2))
# # 		println("e$e: +h$hlist1 -h$hlist2")
# # 		for h in hlist1; println(" +h$h=$(halfedge(s,h))"); end
# # 		for h in hlist2; println(" -h$h=$(halfedge(s,h))"); end
# 		@assert length(hlist1) == length(hlist2)
# 		for i in 1:length(hlist1)-1
# 			opposite!(s, hlist1[i], hlist2[i])
# 			opposite!(s, hlist2[i], hlist1[i+1])
# 		end
# 		opposite!(s, last(hlist1), last(hlist2))
# 		opposite!(s, last(hlist2), first(hlist1))
# 	end
# end#»»
# Constructor««2
function CornerTable(points, faces)#««
	m = CornerTable(points)
	for (v1,v2,v3) in faces
# 		verbose(m)
		append_face!(m, (Vertex(v1), Vertex(v2), Vertex(v3)))
	end
	return m
end#»»
# function CornerTable(points, faces)#««
# 	nv = length(points)
# 	@assert all(length.(faces) .== 3)
# 	s = CornerTable(undef, points, length(faces))
# 	edge_loc = edge_locator(s)
# 
# 	for (i, f) in pairs(faces)
# 		# (v3v1), (v1v2), (v3v2)
# 		mark_face!(s, edge_loc, i, f[1], f[2], f[3])
# 	end
# 	radial_loops!(s, edge_loc)
# 	return s
# end#»»
# find coplanar faces««2
# returns the equivalence relation, as pairs of indices in `flist`:
function coplanar_faces(s::CornerTable, flist, ε = 0)
	@inline box(l,ε) = BBox(l, l .+ ε)
	boxes = [ box(normalized_plane(s, f), ε) for f in flist ]
	return SpatialSorting.intersections(boxes)
end

# opposite faces««2
function opposite_faces(m::CornerTable)
	r = NTuple{2,Vertex{index_type(m)}}[]
	for f in allfaces(m)
		c1 = corner(f, Side(1))
		ismultiple(opposite(m, c1)) || continue
		v2 = vertex(m, f, Side(2))
		v1 = vertex(m, f, Side(1))
		global M=m
		for c in radial(m, c1)
			c <= c1 && continue
			# similarly-oriented faces (including this one) don't count:
			right(m, c) == v2 && continue
			vertex(m, c) == v1 && push!(r, (face(c1), face(c)))
		end
	end
	return r
end
# reverse mesh««2
function Base.reverse(s::CornerTable)
	r = (typeof(s))(undef, points(s), nfaces(s))
	# reverse all half-edges:
	# (31)(12)(23) to (21)(13)(32)
	#
	# edgefrom: 3f-2 <-> 3f, and 3f-1 remains in place
	# destination: 3f-2 in place, 3f-1 ↔ 3f
	# opposite: 3f-2  in place, 3f-1 ↔ 3f
	@inline newedgeno(h) = let m = mod(h,3)
		(m == 1) ? h + 1 : (m == 2) ? h - 1 : h end
	for f in 1:nfaces(s)
		r.destination[3*f  ] = s.destination[3*f-1]
		r.destination[3*f-1] = s.destination[3*f  ]
		r.destination[3*f-2] = s.destination[3*f-2]
		r.opposite[3*f-2] = newedgeno(s.opposite[3*f-1])
		r.opposite[3*f-1] = newedgeno(s.opposite[3*f-2])
		r.opposite[3*f  ] = newedgeno(s.opposite[3*f  ])
	end
	for v in 1:nvertices(s)
		h = s.edgefrom[v]; m = mod(h,3)
		r.edgefrom[v] = (m == 1) ? h+2 : (m == 2) ? h : h-2
	end
	return r
end

# select faces««2
"""
    select_faces!(s, fkept)

Select faces in `s` from vector `fkept`.
"""
function select_faces!(m::CornerTable{I}, fkept) where{I}
	# fkept is a list of kept faces
	global M = m
	global K = fkept
	# first step: compute vertex/face renaming
	vkept = sizehint!(Vertex{I}[], 3*length(fkept))
	for f in fkept
		push!(vkept, vertices(m, f)...)
	end

	unique!(sort!(vkept))
	println("kept vertices: $vkept")

	# FIXME do something about fans!
	r = CornerTable{I}(m.points[Int.(vkept)])
	for f in fkept
		(v1, v2, v3) = vertices(m, f)
		u1 = Vertex(searchsortedfirst(vkept, v1))
		u2 = Vertex(searchsortedfirst(vkept, v2))
		u3 = Vertex(searchsortedfirst(vkept, v3))
		println("append_face!($u1, $u2, $u3)")
		append_face!(r, (u1, u2, u3))
	end
	return r

	# first step: reconnect edge loops as needed
	# FIXME: this will fail if any edge loop has (≠2) selected faces
	# (which should not happen in real cases) (?).
	ekept = similar(fkept, 3*length(fkept))
	for (i, k) in pairs(fkept)
		ekept[3i-2] = ekept[3i-1] = ekept[3i] = k
	end
	
	for e1 in 1:3*nfaces(m)
		ekept[e1] || continue # skip skipped faces
		e2 = opposite(m, e1)
		# case 1: (e1, e2) -> regular edge, do nothing
		# case 2: fkept[e2] -> do nothing
		# case 3: iterate (2-step) on e2 untli
		#   a) e1 is found -> stop
		#   b) fkept[e2] -> set opposite[e1] = e2
		ekept[e2] && continue
		e3 = opposite(m, e2)
		e3 == e1 && continue
		while true
			e2 = opposite(m, e3)
			if ekept[e2]
				opposite!(m, e1, e2)
				break
			end
			e3 = opposite(m, e2)
			e3 == e1 && break
		end
	end
	# second step: compute face/vertex renaming
	# `emap[i]` is the new name (if `ekept[i]`) of half-edge `i`
	emap = cumsum(ekept)
	vkept = falses(nvertices(m))
	for i in 1:length(emap)
		vkept[m.destination[i]] |= ekept[i]
	end
	vmap = cumsum(vkept)
	# third step: proceed in renaming everything
	for i in 1:nvertices(m)
		vkept[i] || continue
		m.edgefrom[vmap[i]] = emap[m.edgefrom[i]]
		m.points[vmap[i]] = m.points[i]
	end
	for i in 1:length(emap)
		ekept[i] || continue
		m.destination[emap[i]] = vmap[m.destination[i]]
		m.opposite[emap[i]] = emap[m.opposite[i]]
	end
	resize_points!(m, last(vmap))
	resize_faces!(m, last(emap) ÷ 3)
	return m
end

# split_faces!««2
"""
    split_faces

Replaces faces in `s`, as indicated by iterator `fsplit`
(as (face number) => (replacement triangles)).
"""
function split_faces!(s::CornerTable, fsplit)
	edge_loc = edge_locator(s)
	# n is number of face being written
	n = nfaces(s)
	for (f, tri) in fsplit
		# we don't skip trivial triangulations: they might have singular edges
		length(tri) ≤ 0 && continue # skipped face, or trivial retriangulation

		# before overwriting this face, we save the location of opposed half-edges,
		(v1,v2,v3) = face(s, f)
		h = opposite(s, 3*f-2)
		halfedge(s, h) == (v1,v3) && mark_edge!(s, edge_loc, h, v1, v3)
		h = opposite(s, 3*f-1)
		halfedge(s, h) == (v2,v1) && mark_edge!(s, edge_loc, h, v2, v1)
		h = opposite(s, 3*f-0)
		halfedge(s, h) == (v3,v2) && mark_edge!(s, edge_loc, h, v3, v2)
		unmark_edge!(s, edge_loc, 3*f-2)
		unmark_edge!(s, edge_loc, 3*f-1)
		unmark_edge!(s, edge_loc, 3*f  )

		resize_faces!(s, nfaces(s) + length(tri) - 1)

# 		println("face $f = $((v1,v2,v3)) => $tri")
		# now iterate over new triangles:
		((a1,a2,a3), state) = iterate(tri)
		# mark the new face
		mark_face!(s, edge_loc, f, a1,a2,a3)
		# append new faces
		while true
			u = iterate(tri, state)
			u == nothing && break
			((a1,a2,a3), state) = u
			n+= 1
			mark_face!(s, edge_loc, n, a1,a2,a3)
		end
	end
	radial_loops!(s, edge_loc)
end
# concatenate««2
function concatenate(mlist::CornerTable...)#««
	r = (typeof(first(mlist)))(vcat(points.(mlist)...))
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
	return r
end#»»
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
		vertex(m, f, Side(TI.index(t, TI.isvertex))))

	# other cases: we create a new point
	if TI.isedge(t)
		idx == nothing &&
			(push!(si.points, p); idx = Vertex(length(si.points) + nvertices(m)))
		k = Side(TI.index(t, TI.isedge))
		# CornerTable has edges in order: edge31, edge12, edge23
		# index(edge31<<2) = index(edge23) = 1 fixes this:
# 		k = Side(TI.index(t<<2, TI.isedge))
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
	T = eltype(point_type(m))

	boxes = [ boundingbox(t...) for t in triangles(m) ]
	si = (points=similar(points(m), 0),
		in_face = SortedDict{Face{I},Vector{Vertex{I}}}(),
		in_edge=NTuple{3,Vertex{I}}[],
		faces = Face{I}[])

	for (f1, f2) in SpatialSorting.intersections(boxes)
		f1 = Face(f1); f2 = Face(f2)
		isadjacent(m, f1, f2) && continue
		it = TriangleIntersections.inter(triangle(m, f1), triangle(m, f2), ε)
		isempty(it) && continue
		
		# create points as needed, and store their index:
		push!(si.faces, f1, f2)
# 		vindex = MVector{6,Vertex{I}}(undef)
		for (i, (p, (t1, t2))) in pairs(it)
			idx = nothing
			idx = insert_point(m, si, p, idx, f1, t1, ε)
			idx = insert_point(m, si, p, idx, f2, t2, ε)

# 			vindex[i] = idx
		end
	end
	unique!(sort!(si.faces))
	return si
end#»»
# subtriangulation««1
# project_and_triangulate ««2
function project_and_triangulate(m::CornerTable, direction, vlist, elist)
	# build matrix of coordinates according to `vlist`:
	    if direction == 1 d = [2,3]
	elseif direction ==-1 d = [3,2]
	elseif direction == 2 d = [3,1]
	elseif direction ==-2 d = [1,3]
	elseif direction == 3 d = [1,2]
	else @assert direction ==-3; d = [2,1]
	end
	vmat = Matrix{Float64}(undef, length(vlist), 2)
	vmat[:,1] .= (point(m, v)[d[1]] for v in vlist)
	vmat[:,2] .= (point(m, v)[d[2]] for v in vlist)

	vref = Int.(vlist)

	emat = Matrix{Int}(undef, length(elist), 2)
	emat[:,1] .= collect(Int(e[1]) for e in elist)
	emat[:,2] .= collect(Int(e[2]) for e in elist)
# 	println("triangulate: $vmat $vlist $emat")

#=
plot '/tmp/a' index 0 u 1:2 w p pt 5, '' index 1 u 1:2:3:4:0 w vectors lc palette lw 3, '' index 0 u 1:2:3 w labels font "bold,14"
=#
# 	for (i, v) in pairs(vlist)
# 		println("$(vmat[i,1])\t$(vmat[i,2])\t$v")
# 	end
# 	println("\n\n")
# 	for i in 1:size(emat,1)
# 		v1 = findfirst(==(emat[i,1]), vlist)
# 		v2 = findfirst(==(emat[i,2]), vlist)
# 		println("$(vmat[v1,1])\t$(vmat[v1,2])\t$(vmat[v2,1]-vmat[v1,1])\t$(vmat[v2,2]-vmat[v1,2])")
# 	end
# 	println("\n\n")
	return ((Vertex(t[1]), Vertex(t[2]), Vertex(t[3]))
		for t in LibTriangle.constrained_triangulation(vmat, vref, emat))
end

# subtriangulate««2
function edge_inserts(s, in_edge)
	sort!(in_edge) # this has the effect of grouping points by edge
	start = 1
	V = eltype(eltype(in_edge))
	inserts = (k = NTuple{2,V}[], v = Vector{V}[])
	while start <= length(in_edge)
		# group by in_edge
		stop = start
		v1 = in_edge[start][1]
		v2 = in_edge[start][2]
		while stop <= length(in_edge) &&
			in_edge[stop][1] == v1 && in_edge[stop][2] == v2
			stop+= 1
		end
		# find best coordinate and sort
		vec = point(s,v2) - point(s,v1)
		proj = main_axis(vec)
		ins = [ project1d(proj, point(s,in_edge[i][3])) for i in start:stop-1 ]
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
	append_points!(m, si.points)
	vmap = simplify_points!(m, ε)
	println("vmap=$vmap")
	explain(m, "/tmp/x.scad", scale=30)

	for i in eachindex(si.in_edge)
		si.in_edge[i] = map(x->get(vmap, x, x) , si.in_edge[i])
	end
	unique!(sort!(si.in_edge))
	for i in eachindex(si.in_face)
		si.in_face[i] = map(x->get(vmap, x, x), si.in_face[i])
	end
	# insert points in edges
	in_edge = edge_inserts(m, si.in_edge)
	println("in_edge: $(si.in_edge) => $in_edge")

	# determine clusters of coplanar (broken) faces:
	coplanar_rel = coplanar_faces(m, si.faces, ε)
	for (i, j) in coplanar_rel
# 		println("
# \e[36;1mcoplanar\e[m: f$i=$(face(s,i)): p$(normalized_plane(s,i)), $(triangle(s,i))
#         ≡ f$j=$(face(s,j)): p$(normalized_plane(s,j)), $(triangle(s,j))")
	end
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
		unique!(sort!(vlist))
		println("$f $ff: $vlist, $elist")
	end
	faces_tri = [ NTuple{3,Vertex{I}}[] for _ in si.faces ]
		
	# iterate over all clusters of broken faces
	for icluster in classes(clusters)
		direction = main_axis(m, si.faces[first(icluster)])[1]

		# for type-stability, we need empty vectors here:
		allvertices = vcat(similar(in_face_v, 0), view(in_face_v, icluster)...)
		alledges = vcat(similar(in_face_e, 0), view(in_face_e, icluster)...)
		unique!(sort!(allvertices))
		unique!(sort!(alledges))

		alltriangles = project_and_triangulate(m, abs(direction),
			allvertices, alledges)
		println("\e[1mcluster $icluster\e[m $(in_face_v[icluster]) => $(collect(alltriangles))")

		for i in icluster
			orient = main_axis(m,si.faces[i])[1] > 0
			isfirst = true
			for tri in alltriangles
				issubset(tri, in_face_v[i]) || continue
				if isfirst
					f = si.faces[i]; cut_face!(m, f); isfirst = false
				else
					nfaces!(m, nfaces(m)+1); f = lastface(m)
				end
				set_face!(m, f, orient ? tri : reverse(tri))
			end
# 			tri = [ t for t in alltriangles if issubset(t, in_face_v[i]) ]
# 			newfaces = (main_axis(m,si.faces[i])[1] > 0) ? tri : reverse.(tri)
# 			println("$i $(si.faces[i]) $(in_face_v[i]): $newfaces)
# 			faces_tri[i] = (main_axis(m,si.faces[i])[1] > 0) ? tri : reverse.(tri)
# 			println("$i $(si.faces[i]) $(in_face_v[i]): $(faces_tri[i])")
		end
	end
	# apply refinement computed above
# 	println(collect(zip(si.faces, faces_tri)))
# 	split_faces!(m, zip(si.faces, faces_tri))
	return m
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
# used for locate_point:
# vec3, proj, dir3, dir2scaled, order
	# prepare geometry information
	(v1, v2) = base(m, c0)
	dir3 = point(m, v2) - point(m, v1)
	axis = main_axis(dir3)
	dir2 = project2d(axis, dir3)
	dir2scaled = dir2 ./dot(dir3, dir3)
	# collect half-edges and corresponding opposed vertices
	# for each adjacent face, compute a (3d) vector which, together with
	# the edge, generates the face (and pointing from the edge to the face):
	# 2d projection of face_vec3 (preserving orientation)
	clist = collect(genradial(m, c0)) # corners around this multiple edge
	# we could use this to determine edge orientation:
# 	dv = [ destination(m, h) == v2 for h in hlist ]
	p1 = point(m, v1)
	fv = [ point(m, vertex(m, c)) - p1 for c in clist ]
# 	fv = [ point(m, v) - p1 for v in vlist ] # face vector
	# face vector, projected in 2d:
	fv2= [ project2d(axis, v) .- dot(v, dir3) .* dir2scaled for v in fv ]
	println("\e[7msort_radial_edge $c0=$(base(m,c0)): $dir3 axis $axis\e[m")
	for (c, y, z) in zip(clist, fv, fv2)
		println("# $c = $(vertex(m,c)) $(base(m, c)): fv3=$y, fv2=$z")
	end
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
	pt3 == nothing && return clist[reorder]

	# find where `pt3 - v1` inserts in radial loop:
	vec3 = pt3 - p1
	vec2 = project2d(axis, vec3) .- dot(vec3, dir3) .*dir2scaled
	println("\e[34mpt3=$pt3, vec3=$vec3, vec2=$vec2\e[m")
	@assert !all(iszero,vec2) "half-edge $h aligned with point $vec3"
	k = searchsorted(fv2[reorder], vec2,
		lt = (u,v)->circular_sign(u,v) > 0)
	@assert k.start > k.stop "impossible to determine location at this edge"
	# possibilities are: (i+1:i) for i in 0..n
	k.stop == 0 && return clist[reorder[end]]
	return clist[reorder[k.stop]]
end#»»
# find_good_halfedge««2
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
	closest = 0; b = false; z = zero(p[1])
	for (i, q) in pairs(points(m))
		cc_label[i] == c || continue
		z1 = norm²(q-p)
		(b && z1 ≥ z) && continue
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
# 	cells = UnionFind(2*ncomp)
	for i1 in 2:ncomp, i2 in 1:i1-1
		c0 = rp.adjacency[i1, i2]
		isboundary(c0) && continue
# 		eindex = rp.adjacency[i1, i2]
# 		iszero(eindex) && continue
		# regular components i and j meet at edge eindex
		clist = sort_radial_loop(m, c0)
# 		hlist = sort_radial_loop(m, eindex)
		n = length(clist)
		# plist is the sorted list of regular patches at this edge
		# dlist is the list of edge orientations
		plist = rp.label[Int.(face.(clist))]
		v2 = right(m, c0)
		olist = [right(m, c) == v2 for c in clist]
		p1 = plist[1]; o1 = olist[1]
		println("sorted radial loop is $clist")
# 		for i in 1:n
# 			println("  $(clist[i]) $(dlist[i]):   face $(face(clist[i]))=$(vertices(m,face(clist[i])))  p$(plist[i])")
# 		end
		for i in 2:n
			p2 = plist[i]; o2 = olist[i]
			# patch p2 is positively oriented iff d2==v2, etc.:
			# if p1 is positively oriented (d1==v2) then cell between p1 and p2
			# is 2p1 (otherwise 2p1-1);  if p2 is positively oriented (d2==v2),
			# this same cell is 2p2-1 (else 2p2)
# 			union!(cells, 2*p1-(d1≠v2), 2*p2-(d2==v2))
			k = 1-o1-o2
			connect!(levels, p1, p2, k)
			p1 = p2; o1 = o2
		end
		# close the loop by identifying both sides of the last cell
		# (the level is already computed by the level structure):
# 		p2 = plist[1]; d2 = dlist[1]
# 		union!(cells, 2*p1-(d1≠v2), 2*p2-(d2==v2))
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
		(f, b) = locate_point(m, face_cc, i2, point(m, vmax[i1]))
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
		face_level[f1]-= 1
		face_level[f2]-= 1
	end
	return face_level
end#»»
# operations««1
const _DEFAULT_EPSILON=1/65536
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
	println("fkept=$fkept")
	select_faces!(newmesh,fkept)
end
@inline Base.union(meshes::CornerTable...; ε = _DEFAULT_EPSILON) =
	combine(meshes, 1, ε)
@inline Base.intersect(meshes::CornerTable...; ε = _DEFAULT_EPSILON) =
	combine(meshes, length(meshes), ε)
@inline Base.setdiff(m1::CornerTable, m2::CornerTable; ε = _DEFAULT_EPSILON)=
	combine([m1, reverse(m2)], 2, ε)
#»»1
function validate(s::CornerTable)
	for (i, j) in pairs(s.opposite)
		j ∈ keys(s.opposite) || println("edge h$i: opposite = $j, invalid")
		h = s.opposite[j]
		halfedge(s, j) == reverse(halfedge(s, i)) ||
		println("edge h$i =$(halfedge(s,i)): opposite h$j = $(halfedge(s,j))")
# 		h == i || println("edge $i: opposite² = $h")
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
	for (i, p) in pairs(s.points)
		println(io, """
translate(s*$(Vector{Float64}(p))) {
	color("red") sphere(1);
	color("black", .8) linear_extrude(1) text("$i", size=5);
}
""")
	end
	println(io, "color(c, .7) polyhedron([")
	join(io, [ " s*$(Vector{Float64}(p))" for p in s.points ], ",\n")
	println(io, "],[")
	join(io, [ " $(Vector{Int}([f...]) .- 1)" for f in faces(s) ], ",\n")
	println(io, "]); } }\n$name();")
end
@inline explain(s::CornerTable, f::AbstractString; kwargs...) =
	open(f, "w") do io explain(s, io; kwargs...) end
function verbose(m::CornerTable, k::Fan)
	print("\e[35;1m$k\e[m: ")
	print(isopenfan(fan_first(m, k)) ? "open " : "closed ")
	print(join(("$c (\e[34m$(vertex(m,c))\e[m $(base(m,c)))" for c in fancorners(m, k)), ","))
	print("  next=\e[35m$(fan_next(m,k))\e[m")
	println()
end
function verbose(m::CornerTable)
	println("\e[7mtable with $(ncorners(m)) corners, $(nvertices(m)) vertices, $(nfans(m)) fans:\e[m")
	for v in allvertices(m)
		c = corner(m, v)
		if isisolated(c)
			println("\e[34;1m$v\e[m: isolated")
		elseif isregular(c)
			print("\e[34;1m$v\e[m: regular first corner $c, ")
			print("star = (", join(("$u" for u in star(m, c)),","), ")")
			println()
		else
			println("\e[34;1m$v\e[m: singular (",
				join(("$k (first=$(fan_first(m,k)))" for k in fans(m,Fan(c))), " "),")")
# 			for k in fans(m, Fan(c))
# 				verbose(m, k)
# 			end
		end
	end
	for f in allfaces(m)
		println("\e[32;1m$f\e[m: ", join(vertices(m, f), ","))
		for c in corners(f)
			o = Int(opposite(m,c))
			println("  \e[33;1m$c\e[m: \e[34m", vertex(m, c), "\e[m",
# 				" \e[35m", fan(m, c), "\e[m",
				" edge(", vertex(m, next(c)), ",", vertex(m, prev(c)), ") ",
				" opposite=", (o > 0 ? "\e[32m" : o < 0 ? "\e[31m" : "\e[38;5;8m"),
				opposite(m, c), "\e[m")
		end
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
