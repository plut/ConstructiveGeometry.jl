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

const _INDEX_TYPE = Int

# For a mesh having `n` vertices, we expect the following counts:
#  - 6n + O(1) half-edges;
#  - 3n + O(1) edges;
#  - 2n + O(1) faces.
# The self-intersection of a mesh is generally a curve; hence we expect
# this self-intersection to involve Θ(√n) faces (and the same order of
# magnitude for edges and vertices).

# tools ««1
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
	merged = lowest_representatives(samepoints, 1:n)
	# thus this gives a descending substitution:
	return (i => j for (i, j) in pairs(merged) if i ≠ j)
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
# several names for the same type, for clarity:
@inline NamedIndex{S}(m::CornerTable) where{S} = NamedIndex{S,index_type(m)}
point_type(::Type{<:CornerTable{I,P}}) where{I,P} = P
point_type(m::CornerTable) = point_type(typeof(m))

# size and resizing functions««2
@inline nvertices(m::CornerTable) = length(m.points)
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
@inline nfaces(m::CornerTable) = length(m.opposite) ÷ 3
@inline ncorners(m::CornerTable) = length(m.opposite)
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
@inline nfans(m::CornerTable) = length(m.fan_first)
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
@inline allvertices(m::CornerTable) = (Vertex(v) for v in 1:nvertices(m))
@inline allcorners(m::CornerTable) = (Corner(c) for c in 1:length(m.opposite))
@inline allfaces(m::CornerTable) = (Face(f) for f in 1:nfaces(m))
@inline allfans(m::CornerTable) = (Fan(k) for k in 1:nfans(m))


# simple accessors ««2
# corners
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
@inline base(m::CornerTable, c::Corner) = (vertex(m,next(c)), vertex(m,prev(c)))
@inline after(m::CornerTable, c::Corner) = next(opposite(m, next(c)))
@inline after!(m::CornerTable, c::Corner, x::Corner) =
	opposite!(m, next(c), prev(x))
@inline before(m::CornerTable, c::Corner) = prev(opposite(m, prev(c)))

# properties of opposite(c)
@inline issimple(c::Corner) = Int(c) > 0
@inline isboundary(c::Corner) = Int(c) == 0
@inline ismultiple(c::Corner) = Int(c) < 0

# faces
@inline face(c::Corner) = Face(fld1(c.i, 3))
@inline corner(f::Face, s::Side) = Corner(3*Int(f)-3+Int(s))
@inline corners(CornerTable, f::Face) =
	(corner(f,Side(1)), corner(f,Side(2)), corner(f,Side(3)))
@inline vertex(m::CornerTable, f::Face, s::Side) = vertex(m, corner(f,s))
@inline vertices(m::CornerTable, f::Face) =
	(vertex(m,f,Side(1)), vertex(m,f,Side(2)), vertex(m,f,Side(3)))
@inline faces(m::CornerTable) =
	(Face(f) => vertices(m,Face(f)) for f in 1:nfaces(m))
@inline function vertices!(m::CornerTable, f::Face, tri::NTuple{3,<:Vertex})
	vertex!(m, corner(f,Side(1)), tri[1])
	vertex!(m, corner(f,Side(2)), tri[2])
	vertex!(m, corner(f,Side(3)), tri[3])
	return m
end

@inline adjacent(m::CornerTable, f::Face, s::Side) =
	face(m, opposite(m, corner(f, s)))
@inline adjacent(m::CornerTable, f::Face) =
	(adjacent(m, f, Side(1)), adjacent(m, f, Side(2)), adjacent(m, f, Side(3)))
@inline isadjacent(m::CornerTable, f1::Face, f2::Face) =
	any(==(f2), adjacent(m, f1))

# vertices
@inline corner(m::CornerTable, v::Vertex) = Corner(m.corner[v.i])
@inline corner!(m::CornerTable, v::Vertex, c::Corner)= m.corner[v.i]=c.i
@inline points(m::CornerTable) = m.points
@inline point(m::CornerTable, v::Vertex) = m.points[v.i]

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

# fans
# properties of `corner(::Vertex)`:
@inline issingular(c::Corner) = Int(c) < 0
@inline isregular(c::Corner) = Int(c) > 0
@inline isisolated(c::Corner) = Int(c) == 0

@inline fan_next(m::CornerTable, k::Fan) = Fan(m.fan_next[Int(k)])
@inline fan_next!(m::CornerTable, k::Fan, x::Fan) = m.fan_next[Int(k)] = Int(x)

@inline Fan(c::Corner) = Fan(-Int(c))
@inline Corner(k::Fan) = Corner(-Int(k))
@inline isclosedfan(c::Corner) = Int(c) > 0
@inline isopenfan(c::Corner) = Int(c) < 0
@inline fan_first(m::CornerTable, k::Fan) = Corner(m.fan_first[Int(k)])
@inline fan_first!(m::CornerTable, k::Fan, x::Integer) =
	m.fan_first[Int(k)] = x
@inline fan_firstcorner(m::CornerTable, k::Fan) =
	Corner(abs(m.fan_first[Int(k)]))
@inline fanvertex(m::CornerTable, k::Fan) = vertex(m, fan_firstcorner(m, k))
@inline vertex_fan!(m::CornerTable, v::Vertex,k::Fan)= m.corner[Int(v)]=-Int(k)

@inline fan_open(c::Corner)   = -Int(c)
@inline fan_closed(c::Corner) = +Int(c)

function fan_prev(m::CornerTable, k1::Fan)
	k2 = k1
	while true
		u = fan_next(m, k2)
		u == k1 && return k2
		k2 = u
	end
end
@inline function create_fan!(m::CornerTable, v::Vertex, x::Integer)
	# if no fan at v
	c = corner(m, v)
	println("create_fan!($v, $x): $c")
	if isisolated(c)
		push!(m.fan_first, x)
		push!(m.fan_next, nfans(m))
		m.corner[Int(v)] = -nfans(m)
		println("nfans = $(nfans(m))")
		println("  starts with $x")
		println("  fan_firstcorner is $(fan_firstcorner(m, Fan(nfans(m))))")
		println("  next corner is $(after(m, Corner(abs(x))))")
	elseif isregular(c)
		# in this case, create two fans, one for the (ex-)regular corners
		push!(m.fan_first, x)
		push!(m.fan_first, Fan(c))
		push!(m.fan_next, nfans(m)+2)
		push!(m.fan_next, nfans(m))
		m.corner[Int(v)] = -nfans(m)
	else # this is a singular vertex: add one fan
		push!(m.fan_first, x)
		push!(m.fan_next, m.fan_next[-Int(c)])
		m.fan_next[-Int(c)] = nfans(m)
	end
	return Fan(nfans(m))
end
function delete_fan!(m::CornerTable, k::Fan)
	# overwrite fan k with last fan
	nk = fan_next(m, k)
	vk = fanvertex(m, k)
	if nk == k
		# if this is the only fan for the point, make it a regular point:
		corner!(m, vk, fan_firstcorner(m, k))
	else
		# remove from cyclic list
		corner!(m, vk, Corner(nk))
		fan_next!(m, fan_prev(m, k), nk)
	end
	# replace all references to last fan with references to k
	l = Fan(nfans(m))
	if l ≠ k
		vertex_fan!(m, fanvertex(m, l), k)
		fan_next!(m, fan_prev(m, l), k)
		fan_next!(m, k, fan_next(m, l))
		fan_first!(m, k, Int(fan_first(m, l)))
	end
	nfans!(m, nfans(m)-1)
end
@inline rename_fans(m::CornerTable, t::AbstractVector{<:Fan}, x::Fan) =
	for i in 1:length(t);
	(Int(t[i]) > nfans(m)) && (t[i] = x); end
	
"glues two open fans together (k1 after k2). This deletes fan k2.
This also assumes that all corners have already been glued."
@inline function glue_fans!(m::CornerTable, ft::AbstractVector{<:Fan}, i1, i2)
	@assert fanvertex(m, ft[i1]) == fanvertex(m, ft[i2])
	@assert isopenfan(fan_first(m, ft[i1])) # we can only glue open fans
	@assert isopenfan(fan_first(m, ft[i2]))
	# glue edges
	if ft[i1] == ft[i2] # we are gluing an open fan with itself: make it closed
		if fan_next(m, ft[i1]) == ft[i1] # any regular closed fan is made implicit
			delete_fan!(m, ft[i1])
		else
			m.fan_first[Int(ft[i1])]*= -1
		end
	else # remove fan ft[i2] from the list
		delete_fan!(m, ft[i2])
	end
	for i in 1:length(ft)
		Int(ft[i]) > nfans(m) && (ft[i] = ft[i2])
	end
	ft[i2] = ft[i1]
end
"splits fan `k` at a simple edge, making `c` its new first corner"
@inline function split_fan!(m::CornerTable, k::Fan, c::Corner)
	if isclosedfan(fan_first(m, k))
		fan_first!(m, k, fan_open(c))
	else
		# create a new fan
		println("splitting (open) fan $k at corner $c")
		verbose(m, k)
		println("vertex is $(vertex(m,c))")
		create_fan!(m, vertex(m, c), fan_open(c))
	end
end

# multiple edges
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
@inline Base.iterate(it::MeshIterator) = (it.start, it.start)
@inline function Base.iterate(it::MeshIterator, s)
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
	MeshIterator{fancorners}(m, Corner(abs(Int(fan_first(m, k)))))
@inline next(it::MeshIterator{fancorners}, m, c) = after(m, c)
@inline stop(it::MeshIterator{fancorners}, m, c) = (c==it.start)||!issimple(c)

@inline radial(m::CornerTable, c::Corner) = MeshIterator{radial}(m, c)
@inline next(it::MeshIterator{radial}, m, c) = Multi(opposite(m, c))
@inline stop(it::MeshIterator{radial}, m, c) = c == it.start

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
	return Fan(0)
end
"iterator on (fan, corner) pairs around a given vertex"
function vertexcorners(m::CornerTable, u::Vertex)
	c0 = corner(m, u)
	if isregular(c0)
		return ((Fan(0), c) for c in star(m, c0))
	elseif issingular(c0)
		k0 = Fan(c0)
		return ((k, c) for k in fans(m, k0) for c in fancorners(m, k))
	else
		return ((Fan(c0), c0) for i in 1:0) # empty with defined eltype
	end
end


# mesh modification««1
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
	klist[i] = klist[i+6] = Fan(0)
	clist[i] = clist[i+3] = clist[i+6] = clist[i+9] = Corner(0)
	for (k, c) in vertexcorners(m, u)
		n = next(c); r = vertex(m, n)
		p = prev(c); l = vertex(m, p)
		(r == v) && (clist[i] = p; clist[i+6] = c)
		(l == v) && (klist[i] = k)
		(r == w) && (clist[i+3] = p; clist[i+9] = c; klist[i+6] = k)
	end
end#»»
# creating vertices, edges, and faces««2
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
	println("""\e[1mglue_edge\e[m
view from $c ($(vertex(m,c))): in=$cin, out=$cout
	fan before=$(klist[i1]), fan after=$(klist[i1+6])""")
	if !isboundary(cin) # there is already an inner corner
		op_in = opposite(m, cin)
		if isboundary(op_in) # op_in is the last edge in its fan
			multi!(m, [c, cin])
		elseif issimple(op_in) # we break a simple edge
			@assert op_in == cout
			multi!(m, [c, cin, cout])
			println("edge ($cin, $cout) is simple")
			verbose(m, klist[i2])
			verbose(m, klist[i3+3])
			split_fan!(m, klist[i2], clist[i2+6])
			split_fan!(m, klist[i3+3], clist[i3+9])
		else # already a multiple edge
			add_multi!(m, Multi(op_in), c)
		end
	elseif !isboundary(cout) # there is an outer corner, but no inner
		op_out = opposite(m, cout)
		@assert !issimple(op_out)
		if ismultiple(op_out)
			add_multi!(m, Multi(op_out), c)
		else # op_out is a boundary: we are creating a simple edge
			println("\e[1mcreate simple edge $c -- $cout\e[m")
			opposite!(m, c, cout); opposite!(m, cout, c)
			glue_fans!(m, klist, i2  , i2+3)
			glue_fans!(m, klist, i3+3, i3+6)
# 			verbose(m, klist[i2])
# 			verbose(m, klist[i3+3])
		end
	end
end#»»
function points!(m::CornerTable, plist)#««
	if Base.IteratorSize(plist) ≠ Base.SizeUnknown()
		verticeshint!(m, nvertices(m) + length(plist))
	end
	for p in plist
		push!(m.points, p)
		push!(m.corner, 0)
	end
	return m
end#»»
function face!(m::CornerTable{I}, vlist) where{I}#««
	n = ncorners(m)
	v1 = Vertex(vlist[1]); v2=Vertex(vlist[2]); v3=Vertex(vlist[3])
	# this table holds 12 previous corners + the 3 new we create
	# around vertex 1, old corners are 1,4,7,10, new is 13, etc.
	klist = MVector{9,Fan{I}}(undef)
	clist = MVector{15,Corner{I}}(undef)
	findedge(m, v1, v2, v3, clist, klist, 1)
	findedge(m, v2, v3, v1, clist, klist, 2)
	findedge(m, v3, v1, v2, clist, klist, 3)

	nfaces!(m, nfaces(m)+1)
	for i in (1,2,3)
		clist[i+12] = Corner(n+i)
		opposite!(m, clist[i+12], Corner(0))
		vertex!(m, clist[i+12], Vertex(vlist[i]))
		klist[i+3] = create_fan!(m, Vertex(vlist[i]), fan_open(clist[i+12]))
	end
	println("klist = $klist\nclist=$(collect(filter(((x,y),)->Int(y)≠0,pairs(clist))))")
	edge!(m, klist, clist, 1)
	edge!(m, klist, clist, 2)
	edge!(m, klist, clist, 3)
	return m
end#»»
function CornerTable(points, faces)#««
	m = CornerTable(points)
	for f in faces
		verbose(m)
		face!(m, f)
	end
	return m
end#»»
# deleting vertices, edges, and faces««2
function kill_edge!(m::CornerTable, c::Corner)#««
	# 1. cut fans if needed
	p = prev(c)
	n = next(c)
	@assert isregular(c)
	op = opposite(m, c)
	if isboundary(op) # nothing to do
	elseif issimple(op)
		# simple edge: we need to split fans at both vertices v2 and v3
		k2 = fan(m, next(op)); split_fan!(m, k2, n)
		k3 = fan(m, prev(op)); split_fan!(m, k3, next(op))
		# delete a simple edge: replace by two boundaries
		opposite!(m, op, Corner(0))
	else # multiple edge; removing one edge may change it to boundary or simple
		# this can work only if there are three edges (c, c1, c2) in the loop...
		c1 = Multi(op)
		c2 = Multi(opposite(m, c1))
		println("removing multiple edge $c ($c1, $c2) $(base(m,c))")
		if c2 == c
			@assert left(m, c1) == left(m, c)
# 			println("  => create boundary edge $c1 => nothing ")
			opposite!(m, c1, Corner(0))
		else
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
function kill_face!(m::CornerTable, f::Face)#««
	c1 = corner(f, Side(1))
	c2 = corner(f, Side(2))
	c3 = corner(f, Side(3))
	kill_edge!(m, c1); kill_edge!(m, c2); kill_edge!(m, c3)
	@assert opposite(m, c1) == Corner(0)
	@assert opposite(m, c2) == Corner(0)
	@assert opposite(m, c3) == Corner(0)
	delete_fan!(m, fan(m, c1))
	delete_fan!(m, fan(m, c2))
	delete_fan!(m, fan(m, c3))
	l = Face(nfaces(m))
	if f ≠ l # replace all references to last face with references to f
		# corners are pointed to as
		#  -   corner(m, vertex)
		#  - ± opposite(m, corner)
		#  - ± fan_first(m, fan)
		for s in Side
			cl = corner(l, s); cf = corner(f, s)
			vl = vertex(m, cl)

			if corner(m, vl) == cl
				corner!(m, vl, ck)
			else
				for k in fans(m, vl)
					fan_first(m, k) == vl && fan_first!(m, k, cf)
					fan_first(m, k) == Fan(vl) && fan_first!(m, k, Fan(cf))
				end
			end

			clop = opposite(m, cl)
			if issimple(clop) # only one edge to change
				opposite!(m, clop, cf)
			elseif ismultiple(clop) # go along radial loop
				for c in radial(m, clop)
					opposite(m, c) == Multi(cl) || continue
					opposite!(m, c, Multi(cf)); break
				end
			end # if isboundary, nothing to do
			
			opposite!(m, cf, opposite(m, cl))
			vertex!(m, cf, vertex(m, cl))
		end
	end
	nfaces!(m, Int(l)-1)
	return m
end#»»
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
function opposite_faces(s::CornerTable)
	r = NTuple{2,face_type(s)}[]
	for f in 1:nfaces(s)
		e1 = 3*f
		v1 = destination(s, e1)
		v3 = destination(s, next(s, e1))
		for e2 in radial_loop(s, e1)
			# similarly-oriented faces (including this one) don't count:
			e2 < e1 && continue
			destination(s, e2) == v1 && continue
			destination(s, next(s, e2)) == v3 && push!(r, fld1.((e1,e2), 3))
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

# change points««2
# function resize_faces!(s::CornerTable, nf)
# 	resize!(s.opposite, 3nf)
# 	resize!(s.destination, 3nf)
# 	s
# end
# function resize_points!(s::CornerTable, nv)
# 	resize!(s.points, nv)
# 	resize!(s.edgefrom, nv)
# 	s
# end
function points!(s::CornerTable, newpoints, ε = 0)
	n = length(newpoints)
	vmap = Vector{Int}(undef, n)
	replace = simplify_points(newpoints, ε)
	start = 1; offset = 0
	for (x, y) in replace
		current = x-offset
		for i in start:current-1; vmap[i+offset] = i; end
		vmap[x] = vmap[y]
		start = current; offset+= 1
	end
	for i in start:n-offset; vmap[i+offset] = i; end
	# offset is length of replace
	
	s.destination .= view(vmap, s.destination)
	resize_points!(s, n)
	for (i, j) in pairs(vmap)
		s.points[j] = newpoints[i]
		s.edgefrom[i] = s.edgefrom[j]
	end
	# we do this *after*, in case s.points == newpoints
	resize_points!(s, n-offset)
	return vmap
end
# select faces««2
"""
    select_faces!(s, fkept)

Select faces in `s` from vector `fkept`.
"""
function select_faces!(s::CornerTable, fkept)
	# first step: reconnect edge loops as needed
	# FIXME: this will fail if any edge loop has (≠2) selected faces
	# (which should not happen in real cases) (?).
	ekept = similar(fkept, 3*length(fkept))
	for (i, k) in pairs(fkept)
		ekept[3i-2] = ekept[3i-1] = ekept[3i] = k
	end
	
	for e1 in 1:3*nfaces(s)
		ekept[e1] || continue # skip skipped faces
		e2 = opposite(s, e1)
		# case 1: (e1, e2) -> regular edge, do nothing
		# case 2: fkept[e2] -> do nothing
		# case 3: iterate (2-step) on e2 untli
		#   a) e1 is found -> stop
		#   b) fkept[e2] -> set opposite[e1] = e2
		ekept[e2] && continue
		e3 = opposite(s, e2)
		e3 == e1 && continue
		while true
			e2 = opposite(s, e3)
			if ekept[e2]
				opposite!(s, e1, e2)
				break
			end
			e3 = opposite(s, e2)
			e3 == e1 && break
		end
	end
	# second step: compute face/vertex renaming
	# `emap[i]` is the new name (if `ekept[i]`) of half-edge `i`
	emap = cumsum(ekept)
	vkept = falses(nvertices(s))
	for i in 1:length(emap)
		vkept[s.destination[i]] |= ekept[i]
	end
	vmap = cumsum(vkept)
	# third step: proceed in renaming everything
	for i in 1:nvertices(s)
		vkept[i] || continue
		s.edgefrom[vmap[i]] = emap[s.edgefrom[i]]
		s.points[vmap[i]] = s.points[i]
	end
	for i in 1:length(emap)
		ekept[i] || continue
		s.destination[emap[i]] = vmap[s.destination[i]]
		s.opposite[emap[i]] = emap[s.opposite[i]]
	end
	resize_points!(s, last(vmap))
	resize_faces!(s, last(emap) ÷ 3)
	return s
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
    insert_point!(s, etc.)

Creates (if needed) a new point in `allpoints` for p,
and returns (in all cases) the index of the new point.

Info stored in `in_face` and `in_edge` is updated according
to intersection type `t`.

If this point already has an index (`idx` ≠ `nothing`)
then nothing new is created
(but `in_face` and `in_edge` are still updated).
"""
function insert_point!(s, si, p, idx, f, t, ε)
	TI = TriangleIntersections

	# easy case: this is a pre-existing vertex
	TI.isvertex(t) && return something(idx, face(s, f)[TI.index(t, TI.isvertex)])

	# other cases: we create a new point
	if TI.isedge(t)
		# CornerTable has edges in order: edge31, edge12, edge23
		# index(edge31<<2) = index(edge23) = 1 fixes this:
		idx == nothing && (push!(si.points, p); idx = length(si.points))
		k = TI.index(t<<2, TI.isedge)
		push!(si.in_edge, (halfedge(s, 3*f-3+k)..., idx))
		return idx
	end
	# point is interior
	if idx == nothing
		for j in f
			norm(si.points[j] - p, Inf) ≤ ε && (idx = j; @goto has_idx)
		end
		push!(si.points, p); idx = length(si.points)
	end
	@label has_idx
	in_face_f = create_entry!(si.in_face, f, [])
	push!(in_face_f, idx)
	return idx
end
"""
    self_intersect(s)

Returns `(points, in_face, in_edge, faces)` describing the
self-intersection graph of `s`.
"""
function self_intersect(s::CornerTable, ε=0)#««
	T = eltype(point_type(s))

	boxes = [ boundingbox(t...) for t in triangles(s) ]
	si = (points=copy(points(s)),
		in_face = SortedDict{face_type(s),Vector{vertex_type(s)}}(),
		in_edge=NTuple{3,vertex_type(s)}[],
		faces = face_type(s)[])

	for (f1, f2) in SpatialSorting.intersections(boxes)
		adjacent_faces(s, f1, f2) && continue
		it = TriangleIntersections.inter(triangle(s, f1), triangle(s, f2), ε)
		isempty(it) && continue
		
		# create points as needed, and store their index:
		push!(si.faces, f1, f2)
		vindex = MVector{6,vertex_type(s)}(undef)
		for (i, (p, (t1, t2))) in pairs(it)
			idx = nothing
			idx = insert_point!(s, si, p, idx, f1, t1, ε)
			idx = insert_point!(s, si, p, idx, f2, t2, ε)

			vindex[i] = idx
		end
	end
	unique!(sort!(si.faces))
	return si
end#»»
# subtriangulation««1
# project_and_triangulate ««2
function project_and_triangulate(points, direction, vlist, elist)
	# build matrix of coordinates according to `vlist`:
	    if direction == 1 d = [2,3]
	elseif direction ==-1 d = [3,2]
	elseif direction == 2 d = [3,1]
	elseif direction ==-2 d = [1,3]
	elseif direction == 3 d = [1,2]
	else @assert direction ==-3; d = [2,1]
	end
	vmat = Matrix{Float64}(undef, length(vlist), 2)
	vmat[:,1] .= (points[v][d[1]] for v in vlist)
	vmat[:,2] .= (points[v][d[2]] for v in vlist)

	emat = Matrix{Int}(undef, length(elist), 2)
	emat[:,1] .= collect(e[1] for e in elist)
	emat[:,2] .= collect(e[2] for e in elist)
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
	return LibTriangle.constrained_triangulation(vmat, vlist, emat)
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
function subtriangulate!(s::CornerTable, ε=0)
	si = self_intersect(s, ε)
	# first renumber points, removing duplicates, including in self-intersect:
	vmap = points!(s, si.points, ε)
	explain(s, "/tmp/x.scad", scale=30)

	for i in eachindex(si.in_edge)
		si.in_edge[i] = map(x->vmap[x], si.in_edge[i])
	end
	unique!(sort!(si.in_edge))
	for i in eachindex(si.in_face)
		si.in_face[i] = map(x->vmap[x], si.in_face[i])
	end
	# insert points in edges
	in_edge = edge_inserts(s, si.in_edge)

	# determine clusters of coplanar (broken) faces:
	coplanar_rel = coplanar_faces(s, si.faces, ε)
	for (i, j) in coplanar_rel
# 		println("
# \e[36;1mcoplanar\e[m: f$i=$(face(s,i)): p$(normalized_plane(s,i)), $(triangle(s,i))
#         ≡ f$j=$(face(s,j)): p$(normalized_plane(s,j)), $(triangle(s,j))")
	end
	clusters = equivalence_structure(length(si.faces), coplanar_rel)
	
	# compute points by face
	in_face_v = [ get(si.in_face, f, vertex_type(s)[]) for f in si.faces ]
	in_face_e = [ NTuple{2,vertex_type(s)}[] for f in si.faces ]
	for (f, vlist, elist) in zip(si.faces, in_face_v, in_face_e)
		ff = face(s, f)
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
		sort!(vlist)
	end
	faces_tri = [ NTuple{3,vertex_type(s)}[] for _ in si.faces ]
		
	# iterate over all clusters of broken faces
	for icluster in classes(clusters)
		direction = main_axis(s, si.faces[first(icluster)])[1]

		# for type-stability, we need empty vectors here:
		allvertices = vcat(vertex_type(s)[], view(in_face_v, icluster)...)
		alledges = vcat(NTuple{2,vertex_type(s)}[], view(in_face_e, icluster)...)
		unique!(sort!(allvertices))
		unique!(sort!(alledges))

		alltriangles = project_and_triangulate(points(s), abs(direction),
			allvertices, alledges)

		for i in icluster
			tri = [ (t[1],t[2],t[3]) for t in alltriangles
				if issubset(t, in_face_v[i]) ]
			faces_tri[i] = (main_axis(s,si.faces[i])[1] > 0) ? tri : reverse.(tri)
		end
	end
	# apply refinement computed above
	split_faces!(s, zip(si.faces, faces_tri))
	return s
end

# multiplicity««1
# regular_patches««2
"""
    regular_patches(mesh)

Returns a pair `(label, adjacency)`, where `label` is a labeling of the
faces of `s` in regular (manifold) patches, and `adjacency` is an adjacency
matrix between these patches (containing a half-edge where both patches meet).
"""
function regular_patches(s::CornerTable)
	label = zeros(Int, nfaces(s))
	todo = face_type(s)[]
	adjacency = zeros(halfedge_type(s), 0, 0)
	n = 0
	for start_face in 1:nfaces(s)
		!iszero(label[start_face]) && continue
		label[start_face] = n+= 1
		push!(todo, start_face)
		adjacency = let new_adjacency = similar(adjacency, n, n)
			new_adjacency[1:n-1, 1:n-1] .= adjacency
			fill!(view(new_adjacency, n, :), 0)
			fill!(view(new_adjacency, 1:n-1, n), 0)
			new_adjacency
		end
		while !isempty(todo)
			current_face = pop!(todo)
			for k in -2:0
				h = 3*current_face+k
				if regular_edge(s, h)
					next_face = fld1(opposite(s, h), 3)
					if iszero(label[next_face])
						label[next_face] = n
						push!(todo, next_face)
					end
				else # singular edge
					for e1 in radial_loop(s, h)
						l = label[fld1(e1,3)]
						!iszero(l) && (adjacency[l,n] = adjacency[n,l] = h)
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
function sort_radial_loop(s::CornerTable, h, pt3 = nothing)
# used for locate_point:
# vec3, proj, dir3, dir2scaled, order
	# prepare geometry information
	(v1, v2) = halfedge(s, h)
	dir3 = point(s, v2) - point(s, v1)
	axis = main_axis(dir3)
	dir2 = project2d(axis, dir3)
	dir2scaled = dir2 ./dot(dir3, dir3)
	# collect half-edges and corresponding opposed vertices
	# for each adjacent face, compute a (3d) vector which, together with
	# the edge, generates the face (and pointing from the edge to the face):
	# 2d projection of face_vec3 (preserving orientation)
	hlist = collect(radial_loop(s, h)) # half-edges
	ov = [ opposite_vertex(s, x) for x in hlist ] #opposite vertices
	# we could use this to determine edge orientation:
# 	dv = [ destination(s, h) == v2 for h in hlist ]
	p1 = point(s, v1)
	fv = [ point(s, v) - p1 for v in ov ] # face vector
	# face vector, projected in 2d:
	fv2= [ project2d(axis, v) .- dot(v, dir3) .* dir2scaled for v in fv ]
# 	println("edge h$h = (v$v1,v$v2): $dir3 axis $axis")
# 	for (e1, x, y, z) in zip(hlist, ov, fv, fv2)
# 		println("# $e1 = $(halfedge(s, e1)): opp.v=$x, fv3=$y, fv2=$z")
# 	end
	face_cmp = @closure (i1, i2) -> let b = circular_sign(fv2[i1], fv2[i2])
		!iszero(b) && return (b > 0)
		# we want a consistent ordering between coincident faces.
		# the rule is: 1. all arrows point inward the thin wedge (this is
		# useful for later determining multiplicity): positive-oriented faces
		# come *before* negative-oriented faces
		# 2. sort by face number (i.e. decreasing face number for rule 1).
		# we have as input an alternated list starting with a positive edge:
		return (-1)^i1*hlist[i1] < (-1)^i2*hlist[i2]
		# half-edge number is proportional to face number
	end
	reorder = sort(1:length(hlist), lt=face_cmp)
	pt3 == nothing && return hlist[reorder]

	# find where `pt3 - v1` inserts in radial loop:
	vec3 = pt3 - p1
	vec2 = project2d(axis, vec3) .- dot(vec3, dir3) .*dir2scaled
	@assert !all(iszero,vec2) "half-edge $h aligned with point $vec3"
	k = searchsorted(fv2[reorder], vec2,
		lt = (u,v)->circular_sign(u,v) > 0)
	@assert k.start > k.stop "impossible to determine location at this edge"
	# possibilities are: (i+1:i) for i in 0..n
	k.stop == 0 && return hlist[reorder[end]]
	return hlist[reorder[k.stop]]
end
# find_good_halfedge««2
function find_good_halfedge(s::CornerTable, i, p)
	# it is possible that all edges lie in the same plane
	# (if this is a flat cell), so we pick, as a second point j,
	# the one which maximizes |y/x|, where
	# y = ‖pi∧pj‖, x=‖pi·ij‖/‖pi‖²
	# in other words, we maximize ‖pi∧pj‖²/‖pi·ij‖²
	# caution: it is possible that pi⋅ij = 0 (edge exactly orthogonal),
	# in this case we must return j immediately
	vpi = point(s, i) - p
	nv = neighbours(s, i)
	((h, j), state) = iterate(nv)
	vpj = point(s, j) - p
	xj = dot(vpi, vpj)
	iszero(xj) && return h
	xyj = (xj*xj, norm²(cross(vpi, vpj)))
	best = h
	while true
		u = iterate(nv, state)
		u == nothing && return best
		((h, k), state) = u
		vpk = point(s, k) - p
		xk = dot(vpi, vpk)
		iszero(xk) && return h
		xyk = (xk*xk, norm²(cross(vpi, vpk)))
		if xyk[2]*xyj[1] > xyk[1]*xyj[2]
			best = h; xyj = xyk
		end
	end
	return best
end
# locate_point««2
"""
    locate_point(s, labels, comp, p)

Returns `(face, flag)`, where `flag` is zero if p lies outside this face, and one if p lies inside this face.
"""
function locate_point(s::CornerTable, cc_label, c, p)
	# find closest vertex to p
	closest = 0; b = false; z = zero(p[1])
	for (i, q) in pairs(points(s))
		cc_label[i] == c || continue
		z1 = norm²(q-p)
		(b && z1 ≥ z) && continue
		closest = i; z = z1
	end
	# find a good edge from closest vertex
	h = find_good_halfedge(s, closest, p)
	y = sort_radial_loop(s, h, p) # y is a half-edge in the radial loop of h
	return (fld1(y, 3), destination(s, y) ≠ destination(s, h))
end
# multiplicity ««2
function multiplicity(s::CornerTable)#««
	rp = regular_patches(s)
	ncomp = size(rp.adjacency, 1)
	# `cells`: allows identification of cells bounded by regular patches
	# cells[2i-1] is inside patch i, cells[2i] is outside it
	# `levels`: level[i] is the multiplicity of cells[2i-1]
	levels = LevelStructure(ncomp)
	# cells is kept implicit for now (this *might* be used later, for
	# reconnecting edges when selecting faces)
# 	cells = UnionFind(2*ncomp)
	for i1 in 2:ncomp, i2 in 1:i1-1
		eindex = rp.adjacency[i1, i2]
		iszero(eindex) && continue
		# regular components i and j meet at edge eindex
		hlist = sort_radial_loop(s, eindex)
		n = length(hlist)
		# plist is the sorted list of regular patches at this edge
		# dlist is the list of edge orientations
		plist = rp.label[fld1.(abs.(hlist), 3)]
		dlist = [destination(s, h) for h in hlist]
		v2 = destination(s, eindex)
		h1 = hlist[1]; p1 = plist[1]; d1 = dlist[1]
# 		println("sorted radial loop is $hlist")
# 		for i in 1:n
# 			println("  h$(hlist[i]) to v$(dlist[i]):   f$(face(s, fld1(hlist[i],3))) p$(plist[i])")
# 		end
		for i in 2:n
			h2 = hlist[i]; p2 = plist[i]; d2 = dlist[i]
			# patch p2 is positively oriented iff d2==v2, etc.:
			# if p1 is positively oriented (d1==v2) then cell between p1 and p2
			# is 2p1 (otherwise 2p1-1);  if p2 is positively oriented (d2==v2),
			# this same cell is 2p2-1 (else 2p2)
# 			union!(cells, 2*p1-(d1≠v2), 2*p2-(d2==v2))
			k = 1-(d1==v2)-(d2==v2)
			connect!(levels, p1, p2, k)
			h1 = h2; p1 = p2; d1 = d2
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
	vertex_cc = similar(cc_list, nvertices(s))
	for (j, f) in pairs(faces(s)), v in f
		vertex_cc[v] = levels.parent[rp.label[j]]
	end
	vmax = similar(cc_list)
	for (i, c) in pairs(cc_list)
		b = false; z = zero(point(s,1)[1])
		for (v, l) in pairs(vertex_cc)
			l ≠ c && continue
			z1 = point(s, v)[1]
			(b && z1 ≤ z) && continue
			vmax[i] = v; z = z1; b = true
		end
	end
	# compute the nesting level for all connected components
	cc_nest = zeros(Int, length(cc_list))
	for i1 in 1:length(cc_list), i2 in 1:length(cc_list)
		i1 == i2 && continue
		(f, b) = locate_point(s, face_cc, i2, point(s, vmax[i1]))
		k = b + levels.level[rp.label[f]]
		# if k > 0 then i1 inside i2
		cc_nest[i1]+= k
	end

	face_level = Vector{Int}(undef, nfaces(s))
	for f in 1:nfaces(s)
		c = searchsortedfirst(cc_list, face_cc[f])
		face_level[f] = 1 + levels.level[rp.label[f]] + cc_nest[c]
	end
	for (f1, f2) in opposite_faces(s)
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
	select_faces!(newmesh, levels .== μ)
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
	print("  \e[35;1m$k\e[m: ")
	print(isopenfan(fan_first(m, k)) ? "open " : "closed ")
	print(join(("$c" for c in fancorners(m, k)), ","))
	println()
end
function verbose(m::CornerTable)
	println("\e[7mtable with $(ncorners(m)) corners, $(nvertices(m)) vertices, $(nfans(m)) fans:\e[m")
	for v in allvertices(m)
		c = corner(m, v)
		if isisolated(c)
			println("\e[34;1m$v\e[m: isolated")
		elseif isregular(c)
			println("\e[34;1m$v\e[m: regular vertex, first corner $c, star = (",
				join(("$u" for u in star(m, c)),","), ")")
		else
			println("\e[34;1m$v\e[m: singular vertex, fans($c) = ",
				join(("$k (first=$(fan_first(m,k)))" for k in fans(m, Fan(c))), " "))
			for k in fans(m, Fan(c))
				verbose(m, k)
			end
		end
	end
	for f in allfaces(m)
		println("\e[32;1m$f\e[m: ", join(vertices(m, f), ","))
		for c in corners(m,f)
			o = Int(opposite(m,c))
			println("  \e[33;1m$c\e[m: ", vertex(m, c),
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
end

export CornerTable
end # module
