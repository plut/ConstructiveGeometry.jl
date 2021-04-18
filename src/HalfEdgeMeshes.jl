"""
    HalfEdgeMeshes

This module contains the basic function for operating with half-edge meshes with triangular faces.
For a general introduction to half-edge data structures,
see [Rhodes 2013](https://ubm-twvideo01.s3.amazonaws.com/o1/vault/gdc2013/slides/822512 Rhodes_Graham_Math_for_Game \\(2\\).pdf).

This module exports the `HalfEdgeMesh` data type and defines the following functions (not exported):
 - `faces`, `points`: exports a mesh to a list of triangles.
 - `concatenate`: disjoint union of meshes.
 - `combine`: compute the result of a Boolean operation (intersection, union)
   on meshes.
 - `Base.reverse`: computes the complement of a mesh; together with `combine`,
   this allows computation of a Boolean difference of meshes.
 - `Base.union`, `Base.intersect`, `Base.setdiff`: shorthand for usual combinations.
"""
module HalfEdgeMeshes
using StaticArrays
using LinearAlgebra
using FastClosures
using DataStructures
module LibTriangle
	using Triangle
end
include("TriangleIntersections.jl")
include("SpatialSorting.jl")

const _HALFEDGE_INDEX_TYPE = Int

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
# equivalence relations««1
# uniquenames ««2
"""
    uniquenames(list)

Returns `(names, idx)`, where `names = unique(sort(list))`,
and `list == names[idx]`.

"""
function uniquenames(list)
	# used by: equivalence_structure
	perm = sortperm(list)
	x = list[perm[1]]; names = [x] # smallest name
	idx = similar(perm)
	idx[perm[1]] = j = 1
	for i in 2:length(perm)
		k = perm[i]
		y = list[k]
		if y ≠ x
			x = y
			push!(names, x)
			j+= 1
		end
		idx[k] = j
	end
	return (names, idx)
end

# union-find structure««2
"""
    AbstractUnionFind

Supertype for union-find-like algorithms.
Implements: parent tree, encoded as a vector `u.parent`.
"""
abstract type AbstractUnionFind{T} end

function root(u::AbstractUnionFind, i)
	while true
		j = u.parent[i]
		i == j && return i
		i = j
	end
end
@inline representatives(u::AbstractUnionFind) =
	[ root(u, i) for i in 1:length(u.parent) ]
# @inline equivalent(u::AbstractUnionFind, i, j) = root(u, i) == root(u, j)

struct UnionFind{T} <: AbstractUnionFind{T}
	parent::Vector{T}
	treesize::Vector{T}
	@inline UnionFind(n::Integer) =
		new{typeof(n)}(collect(1:n), ones(typeof(n), n))
end
function Base.union!(u::UnionFind, i, j)
	x = root(u, i)
	y = root(u, j)
	x == y && return
	if u.treesize[x] < u.treesize[y]
		u.parent[x] = u.parent[y]; u.treesize[y]+= u.treesize[x]
	else
		u.parent[y] = u.parent[x]; u.treesize[x]+= u.treesize[y]
	end
	return u
end

struct LevelStructure{T} <: AbstractUnionFind{T}
	parent::Vector{T}
	treesize::Vector{T}
	level::Vector{T} # level[i] is difference between i and parent[i]
	@inline LevelStructure(n::Integer) =
		new{typeof(n)}(collect(1:n), ones(typeof(n), n), zeros(typeof(n), n))
end
@inline Base.eltype(u::LevelStructure) = eltype(u.parent)
function root_and_level(u::LevelStructure, i)
	l = 0
	while true
		j = u.parent[i]
		l+= u.level[i]
		i == j && return (i, l)
		i = j
	end
end
"""
    root!(u::LevelStructure)

Returns a standard form for `u`, where each entry has, as its parent, the
first, lowest-level equivalent entry.
"""
function root!(u::LevelStructure)
	for i in 1:length(u.parent)
		(r, l) = root_and_level(u, i)
		u.parent[i] = r
		u.level[i] = l
	end
	minlevel = SortedDict{eltype(u),NTuple{2,eltype(u)}}()
	for i in 1:length(u.parent)
		r = u.parent[i]; l = u.level[i]
		idx = findkey(minlevel, r)
		if status((minlevel, idx)) == 1 # real token
			m = minlevel[idx]
			m[2] > l && (minlevel[idx] = (i, l))
		else # past-end
			push!(minlevel, r => (i, l))
		end
	end
	for i in 1:length(u.parent)
		r = u.parent[i]; l = u.level[i]; m = minlevel[r]
		u.parent[i] = m[1]
		u.level[i] = l - m[2]
	end
	return u
end

# sets level[i2] = level[i1] + δ
function connect!(u::LevelStructure, i1, i2, δ)
	(r1, d1) = root_and_level(u, i1) # l(i1) = l(r1) + d1
	(r2, d2) = root_and_level(u, i2) # l(i2) = l(r2) + d2
	# we now want l(i2) = l(i1) + δ, hence
	# l(r2) + d2 = l(r1) + d1 + δ
	# l(r2) = l(r1) + d1 - d2 + δ
	if r1 == r2
		if d2 ≠ d1 + δ
			str="connect($i1, $i2, $δ): faces are already connected and multiplicity should be $(d1+δ) (it is $d2). Faces are likely crossed around an edge."
			l = 0; i =i1
			str*="\n path for $i1:"
			while true
				str*="\n   ($i, level $l) parent[$i]=$(u.parent[i]), level[$i]=$(u.level[i])"
				j = u.parent[i]; l+= u.level[i]
				i == j && break
				i = j
			end
			l = 0; i =i2
			str*="\n path for $i2:"
			while true
				str*="\n   ($i, level $l) parent[$i]=$(u.parent[i]), level[$i]=$(u.level[i])"
				j = u.parent[i]; l+= u.level[i]
				i == j && break
				i = j
			end
			@assert d2==d1+δ str
		end
		return
	end
	if u.treesize[r1] < u.treesize[r2]
		u.parent[r1] = u.parent[r2]; u.treesize[r2]+= u.treesize[r1]
		u.level[r1] = d2 - d1 -δ
	else
		u.parent[r2] = u.parent[r1]; u.treesize[r1]+= u.treesize[r2]
		u.level[r2] = d1 - d2 + δ
	end
	return u
end

# find representative given equivalence relation ««2
# here we assume that we are already provided with the transitive cloture
# of the relation (as returned by SpatialSorting.intersections);
# hence at some point we *will* encounter the lowest representative.
function lowest_representatives(relation, elements)
	rep = collect(elements) # works for 1:n or a vector
	for (x, y) in relation
		i = searchsorted(elements, x)
		j = searchsorted(elements, y)
		if rep[i] < rep[j]
			rep[j] = rep[i]
		else
			rep[i] = rep[j]
		end
	end
	return rep
end

function lowest_representatives(relation)
	elements = unique!(sort!([r[i] for r in relation for i in 1:2]))
	return (elements, lowest_representatives(relation, elements))
end

# Equivalence structure for arbitrary (sortable) objects««2
"""
    EquivalenceStructure

Equivalence structure for arbitrary (sortable) objects

Useful methods:
 - `elements(s)`: returns (unique+sorted) list of all elements
 - `equivalent(s, x,y)`: true iff x, y are equivalent
 - `class(s, x)`: returns all elements equivlaent to x
 - `representative(s,x)`: first(class(s,x))
"""
# from an element, find an equivalence class identifier (integer)
# from an equivalence class id, find (sorted) list of equiv. elements
struct EquivalenceStructure{E}
	elements::E # element list: either a sorted Vector, or OneTo (for Int)
	class_idx::Vector{Int}
	class::Vector{Vector{Int}}
end
@inline elements(s::EquivalenceStructure) = s.elements
@inline classes(s::EquivalenceStructure) = (s.elements[c] for c in s.class)
@inline function index(s::EquivalenceStructure, x)
	idx = searchsortedfirst(s.elements, x)
	@boundscheck (idx > length(s.elements) || s.elements[idx] ≠ x) &&
		return nothing
# 		throw(KeyError(x))
	return idx
end
@inline class_idx(s::EquivalenceStructure, x) = s.class_idx[index(s, x)]
@inline equivalent(s::EquivalenceStructure, x, y) = index(s, x) == index(s, y)
@inline class(s::EquivalenceStructure, x) = s.elements[s.class[class_idx(s, x)]]
@inline representative(s::EquivalenceStructure, x) =
	s.elements[first(s.class[class_idx(s,x)])]

function equivalence_structure_flat(n::Integer, rel_flat)#««
	# rel_flat is a flat list: r[1] ∼ r[2], r[3] ∼ r[4] etc.
	r = length(rel_flat)
	# compute transitive closure and representatives of each class:
	uf = UnionFind(n)
	for i in 1:2:length(rel_flat)
		union!(uf, rel_flat[i], rel_flat[i+1])
	end
	rep = representatives(uf)

	# now index the representatives
	(rep_names, rep_idx) = uniquenames(rep)

	# index the classes
	class = [ Int[] for _ in rep_names ]
	@inbounds for (i, r) in pairs(rep_idx)
		push!(class[r], i)
	end
	map(sort!, class)
	return EquivalenceStructure(Base.OneTo(n), rep_idx, class)
end#»»
function equivalence_structure(n::Integer, relation)#««
	isempty(relation) &&
		return EquivalenceStructure(Base.OneTo(0), Int[], Vector{Int}[])
	eq1 = equivalence_structure_flat(n, reinterpret(Int, relation))
	return EquivalenceStructure(Base.OneTo(n), eq1.class_idx, eq1.class)
end#»»
function equivalence_structure(relation)#««
	isempty(relation) &&
		return EquivalenceStructure(eltype(relation)[], Int[], Vector{Int}[])
	# flatten relation:
	rel_flat = reinterpret(eltype(eltype(relation)), relation)
# 	# first compute names
# 	rel_flat = [ r[i] for r in relation for i in 1:2 ]

	# `names` is the unique list of names
	# `rel_idx` is the renamed equivalence relation (as small integers):
	(names, indexed_relation) = uniquenames(rel_flat)
	eq1 = equivalence_structure_flat(length(names), indexed_relation)
	return EquivalenceStructure(names, eq1.class_idx, eq1.class)
end#»»
# half-edge mesh««1
# half-edge data structure ««2
struct HalfEdgeMesh{I<:Signed,P} # I is index type (integer), P is point type
	                       # of type:
	points::Vector{P}      # v->P
	opposite::Vector{I}    # h->h (or -h if singular edge)
	destination::Vector{I} # h->v
	edgefrom::Vector{I}    # v->h (or -c if singular vertex)
# 	cone_next::Vector{I}   # c->c
# 	cone_edge::Vector{I}   # c->h
end

# struct MeshException <: Exception # cone
# 	what
# end

@inline index_type(::Type{<:HalfEdgeMesh{I}}) where{I} = I
# several names for the same type, for clarity:
@inline vertex_type(m::HalfEdgeMesh) = index_type(typeof(m))
@inline halfedge_type(m::HalfEdgeMesh) = index_type(typeof(m))
@inline face_type(m::HalfEdgeMesh) = index_type(typeof(m))
# @inline cone_type(m::HalfEdgeMesh) = index_type(typeof(m))
point_type(::Type{<:HalfEdgeMesh{I,P}}) where{I,P} = P
point_type(m::HalfEdgeMesh) = point_type(typeof(m))

@inline HalfEdgeMesh{I,P}(::UndefInitializer, points, nf::Integer
	) where{I,P} = HalfEdgeMesh{I,P}(
	points,
	Vector{I}(undef, 3nf),
	Vector{I}(undef, 3nf),
	Vector{I}(undef, length(points)),
# 	I[], # cone
# 	I[],
	)

@inline HalfEdgeMesh(::UndefInitializer, points, args...) =
	HalfEdgeMesh{_HALFEDGE_INDEX_TYPE,eltype(points)}(undef, points, args...)

# simple accessors ««2
@inline nfaces(m::HalfEdgeMesh) = length(m.opposite) ÷ 3
@inline nvertices(m::HalfEdgeMesh) = length(m.edgefrom)
# @inline ncones(m::HalfEdgeMesh) = length(m.cone_next)
# operations on half-edges
@inline next(::HalfEdgeMesh, h) = h+1 - 3*(h%3==0)
@inline prev(::HalfEdgeMesh, h) = h-1 + 3*(h%3==1)
@inline opposite(m::HalfEdgeMesh, h) = m.opposite[h]
@inline opposite!(m::HalfEdgeMesh, h, x) = m.opposite[h] = x
@inline destination(m::HalfEdgeMesh, h) = m.destination[h]
@inline destination!(m::HalfEdgeMesh, h, v) = m.destination[h] = v
@inline hisregular(m::HalfEdgeMesh, h) = opposite(m, h) > 0
# operations on vertices
@inline edgefrom(m::HalfEdgeMesh, v) = m.edgefrom[v]
@inline edgefrom!(m::HalfEdgeMesh, v, h) = m.edgefrom[v] = h
@inline visregular(m::HalfEdgeMesh, v) = edgefrom(m, v) > 0
@inline points(m::HalfEdgeMesh) = m.points
@inline point(m::HalfEdgeMesh, v) = m.points[v]
# operatoins on cones
# @inline firstedge(m::HalfEdgeMesh, c) = m.cone_edge[c]

# O(1) operators««2
# faces
@inline adjacent_faces(s::HalfEdgeMesh, i, j) =
	fld1(opposite(s, 3i-2), 3) == j ||
	fld1(opposite(s, 3i-1), 3) == j ||
	fld1(opposite(s, 3i  ), 3) == j
@inline face_edge(s::HalfEdgeMesh, f, j) = halfedge(s, 3*f-3+j)
@inline face_vertex(s::HalfEdgeMesh, f, j) = destination(s, 3*f-3+j)
@inline adjacent_face(s::HalfEdgeMesh, f, j) = fld1(opposite(s, 3*f-3+j), 3)
@inline adjacent_faces(s::HalfEdgeMesh, f) =
	map(j->adjacent_face(s, f, j), (1,2,3))

# half-edges
@inline origin(m::HalfEdgeMesh, h) = destination(m, prev(m, h))
@inline halfedge(m::HalfEdgeMesh, h) = (origin(m, h), destination(m, h))
@inline opposite_vertex(m::HalfEdgeMesh, h) =
	destination(m, next(m, h))
# 1 = regular, 0 = border, -1 = singular
@inline edge_type(m::HalfEdgeMesh, h) = sign(opposite(m,h))
@inline regular_edge(m::HalfEdgeMesh, h) = opposite(m, opposite(m, h)) == h

# cones
# @inline cone_vertex(m::HalfEdgeMesh, c) = origin(m, firstedge(m, c))

# face/triangle iterators««2
struct HalfEdgeFaceIterator{I} <: AbstractVector{NTuple{3,I}}
	mesh::HalfEdgeMesh{I}
end

@inline Base.size(h::HalfEdgeFaceIterator) = (nfaces(h.mesh),)
@inline Base.getindex(h::HalfEdgeFaceIterator, i::Integer) =
	(destination(h.mesh,3i-2), destination(h.mesh,3i-1), destination(h.mesh,3i))
@inline faces(s::HalfEdgeMesh) = HalfEdgeFaceIterator(s)
@inline face(s::HalfEdgeMesh, f) = faces(s)[f]

struct HalfEdgeTriangleIterator{I,P} <: AbstractVector{NTuple{3,P}}
	mesh::HalfEdgeMesh{I,P}
end

@inline Base.size(it::HalfEdgeTriangleIterator) = (nfaces(it.mesh),)
@inline Base.getindex(it::HalfEdgeTriangleIterator, i::Integer) = (
	point(it.mesh,destination(it.mesh, 3i-2)),
	point(it.mesh,destination(it.mesh, 3i-1)),
	point(it.mesh,destination(it.mesh, 3i  )),
	)
@inline triangles(s::HalfEdgeMesh) =
	HalfEdgeTriangleIterator{halfedge_type(s),point_type(s)}(s)
@inline triangle(s::HalfEdgeMesh, f) = triangles(s)[f]
@inline function normal(s::HalfEdgeMesh, f)
	t = triangle(s, f)
	return cross(t[2]-t[1], t[3]-t[2])
end
@inline main_axis(s::HalfEdgeMesh, f) = main_axis(normal(s, f))
function normalized_plane(s::HalfEdgeMesh, f)
	t = triangle(s, f)
	d = cross(t[2]-t[1], t[3]-t[2])
	a = main_axis(d)
	c = @inbounds inv(d[abs(a)])
	return SA[oftype(c, a), project2d(a, d) .* c..., dot(d, t[1])*c]
end
@inline volume(s::HalfEdgeMesh) =
	sum(dot(u, cross(v, w)) for (u,v,w) in triangles(s))/6

# star iterator««2
struct HalfEdgeVertexEdges{I}
	mesh::HalfEdgeMesh{I}
	start::I # half-edge
end
@inline function neighbours(m::HalfEdgeMesh, v, h = edgefrom(m, v))
	# the parameter `h` allows to pass the first half-edge of a singular cone
# 	h > 0 || throw(MeshException(("singular vertex", v)))
	return HalfEdgeVertexEdges(m, h)
end
@inline Base.eltype(::Type{HalfEdgeVertexEdges{I}}) where{I} = Pair{I,I}
@inline Base.iterate(it::HalfEdgeVertexEdges) =
	(it.start => destination(it.mesh, it.start), it.start)
@inline function Base.iterate(it::HalfEdgeVertexEdges, s)
	o = opposite(it.mesh, s); iszero(o) && return nothing
	h = next(it.mesh, o); h == it.start && return nothing
	return (h => destination(it.mesh, h), h)
end
Base.IteratorSize(::HalfEdgeVertexEdges) = Base.SizeUnknown()
function edge(m::HalfEdgeMesh, (v1, v2))
	for (h, dest) in neighbours(m, v1)
		dest == v2 && return h
	end
	throw(KeyError((v1,v2))); return edgefrom(v1) # type-stability
end

# # cone iterator««2
# struct ConeIterator{I}
# 	mesh::HalfEdgeMesh{I}
# 	start::I
# end
# @inline Base.eltype(it::ConeIterator) = cone_type(it.mesh)
# @inline Base.IteratorSize(::ConeIterator) = Base.SizeUnknown()
# 
# @inline function cones(m::HalfEdgeMesh, v)
# 	c = edgefrom(m, v)
# 	c < 0 || throw(MeshException(("regular vertex", v)))
# 	return ConeIterator(m, -c)
# end
# function Base.iterate(it::ConeIterator) = (it.start, it.start)
# function Base.iterate(it::ConeIterator, s)
# 	t = it.mesh.cone_next(s)
# 	t == it.start && return nothing
# 	return (t, t)
# end
# radial iterator««2
struct HalfEdgeRadialIterator{I}
	mesh::HalfEdgeMesh{I}
	start::I
end
@inline radial_loop(m::HalfEdgeMesh, h) = HalfEdgeRadialIterator(m, h)
@inline Base.iterate(it::HalfEdgeRadialIterator) = (it.start, it.start)
function Base.iterate(it::HalfEdgeRadialIterator, s)
	h = opposite(it.mesh, s)
# 	h > 0 || throw(MeshException(("non-singular edge", s, h))) # cone
	h == it.start && return nothing
	return (h, h)
end
@inline Base.IteratorSize(::HalfEdgeRadialIterator) = Base.SizeUnknown()
@inline Base.eltype(::HalfEdgeRadialIterator{I}) where{I} = I

# constructor from points and faces««2
@inline edge_locator(s::HalfEdgeMesh) =
	SortedDict{NTuple{2,vertex_type(s)},Vector{halfedge_type(s)}}()
function mark_edge!(s::HalfEdgeMesh, edge_loc, h, v1, v2)#««
	destination!(s, h, v2)
	edgefrom!(s, v1, h)
	push_entry!(edge_loc, (v1,v2), h)
end#»»
function mark_face!(s::HalfEdgeMesh, edge_loc, i, v1, v2, v3)#««
	mark_edge!(s, edge_loc, 3i-2, v3, v1)
	mark_edge!(s, edge_loc, 3i-1, v1, v2)
	mark_edge!(s, edge_loc, 3i  , v2, v3)
end#»»
function unmark_edge!(s::HalfEdgeMesh, edge_loc, h)
	l = get(edge_loc, halfedge(s,h), halfedge_type(s)[])
	filter!(≠(h), l)
end
function radial_loops!(s::HalfEdgeMesh, edge_loc)#««
	for (e, hlist1) in edge_loc
		e[1] > e[2] && continue # we do two this only once per edge
		isempty(hlist1) && continue
		hlist2 = edge_loc[reverse(e)]
		unique!(sort!(hlist1))
		unique!(sort!(hlist2))
# 		println("e$e: +h$hlist1 -h$hlist2")
# 		for h in hlist1; println(" +h$h=$(halfedge(s,h))"); end
# 		for h in hlist2; println(" -h$h=$(halfedge(s,h))"); end
		@assert length(hlist1) == length(hlist2)
		for i in 1:length(hlist1)-1
			opposite!(s, hlist1[i], hlist2[i])
			opposite!(s, hlist2[i], hlist1[i+1])
		end
		opposite!(s, last(hlist1), last(hlist2))
		opposite!(s, last(hlist2), first(hlist1))
	end
end#»»
function HalfEdgeMesh(points, faces)#««
	nv = length(points)
	@assert all(length.(faces) .== 3)
	s = HalfEdgeMesh(undef, points, length(faces))
	edge_loc = edge_locator(s)

	for (i, f) in pairs(faces)
		# (v3v1), (v1v2), (v3v2)
		mark_face!(s, edge_loc, i, f[1], f[2], f[3])
	end
	radial_loops!(s, edge_loc)
	return s
end#»»

# find coplanar faces««2
# returns the equivalence relation, as pairs of indices in `flist`:
function coplanar_faces(s::HalfEdgeMesh, flist, ε = 0)
	@inline box(l,ε) = BBox(l, l .+ ε)
	boxes = [ box(normalized_plane(s, f), ε) for f in flist ]
	return SpatialSorting.intersections(boxes)
end

# opposite faces««2
function opposite_faces(s::HalfEdgeMesh)
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
function Base.reverse(s::HalfEdgeMesh)
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
function resize_faces!(s::HalfEdgeMesh, nf)
	resize!(s.opposite, 3nf)
	resize!(s.destination, 3nf)
	s
end
function resize_points!(s::HalfEdgeMesh, nv)
	resize!(s.points, nv)
	resize!(s.edgefrom, nv)
	s
end
function points!(s::HalfEdgeMesh, newpoints, ε = 0)
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
function select_faces!(s::HalfEdgeMesh, fkept)
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
function split_faces!(s::HalfEdgeMesh, fsplit)
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
function concatenate(slist::HalfEdgeMesh...)
	r = (typeof(first(slist)))(undef,
		vcat(points.(slist)...),
		sum(nfaces.(slist)))
	eoffset = 0
	voffset = 0
	for s in slist
		(nv, ne) = length(s.edgefrom), length(s.opposite)
		r.opposite[eoffset+1:eoffset+ne] .= s.opposite .+ eoffset
		r.destination[eoffset+1:eoffset+ne] .= s.destination .+ voffset
		r.edgefrom[voffset+1:voffset+nv] .= s.edgefrom .+ eoffset
		eoffset += ne
		voffset += nv
	end
	return r
end
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
		# HalfEdgeMesh has edges in order: edge31, edge12, edge23
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
function self_intersect(s::HalfEdgeMesh, ε=0)#««
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
function subtriangulate!(s::HalfEdgeMesh, ε=0)
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
function regular_patches(s::HalfEdgeMesh)
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
function sort_radial_loop(s::HalfEdgeMesh, h, pt3 = nothing)
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
function find_good_halfedge(s::HalfEdgeMesh, i, p)
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
function locate_point(s::HalfEdgeMesh, cc_label, c, p)
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
function multiplicity(s::HalfEdgeMesh)#««
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
@inline Base.union(meshes::HalfEdgeMesh...; ε = _DEFAULT_EPSILON) =
	combine(meshes, 1, ε)
@inline Base.intersect(meshes::HalfEdgeMesh...; ε = _DEFAULT_EPSILON) =
	combine(meshes, length(meshes), ε)
@inline Base.setdiff(m1::HalfEdgeMesh, m2::HalfEdgeMesh; ε = _DEFAULT_EPSILON)=
	combine([m1, reverse(m2)], 2, ε)
#»»1
function validate(s::HalfEdgeMesh)
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
function explain(s::HalfEdgeMesh, h)
	opp=opposite(s, h)
	j=fld1(h, 3)
	print("""
half-edge $h: ->$(destination(s,h)), opp=$opp -> $(destination(s,opp))
  in triangle $j with $(map(i->destination(s,i),(3j-2,3j-1,3j)))
""")
end
function explain(s::HalfEdgeMesh, io::IO = stdout;
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
@inline explain(s::HalfEdgeMesh, f::AbstractString; kwargs...) =
	open(f, "w") do io explain(s, io; kwargs...) end

export HalfEdgeMesh
end # module
