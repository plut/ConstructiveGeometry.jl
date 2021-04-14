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
"""
module HalfEdgeMeshes
using StaticArrays
using LinearAlgebra
using FastClosures
using Dictionaries
module LibTriangle
	using Triangle
end
include("TriangleIntersections.jl")
include("SpatialSorting.jl")

# using DataStructures
const HalfEdgeId = Int
const VertexId = Int

# For a mesh having `n` vertices, we expect the following counts:
#  - 6n + O(1) half-edges;
#  - 3n + O(1) edges;
#  - 2n + O(1) faces.
# The self-intersection of a mesh is generally a curve; hence we expect
# this self-intersection to involve Θ(√n) faces (and the same order of
# magnitude for edges and vertices).
# This has implication on the data structures used for storing this
# intersection. We first define two types for sparse tables:
#  - `SparseTable0` for n×n tables with o(n) entries,
#  - `SparseTable1` for n×n tables with Θ(n) entries.
# More generally, assuming the matrix has `n` entries,
# `SparseTable0` has size O(n+q) and acces O(1+q/n),
# while `SparseTable1` has size O(q) and access O(q/n).

# tools ««1
@inline norm²(v) = dot(v,v)
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
# SortedVector««2
struct SortedVector{K} <: AbstractVector{K}
	elements::Vector{K}
end
Base.getindex(s::SortedVector, i) = getindex(s.elements, i)
Base.size(s::SortedVector) = size(s.elements)
Base.convert(T::Type{<:SortedVector}, v::AbstractVector) = T(v)
@inline Base.in(x, s::SortedVector) = !isempty(searchsorted(s, x))
	
@inline function find(s::SortedVector, k)
	i = searchsorted(s, k)
	!isempty(i) && return first(i)
	return nothing
end
@inline function find(s::SortedVector, k, ::typeof(insert!))
	i = searchsorted(s, k)
	f = first(i)
	!isempty(i) && return (f, false)
	insert!(s.elements, f, k)
	return (f, true)
end
@inline function find(s::SortedVector, k, ::typeof(push!))
	(i, c) = find(s, k, insert!)
	c || throw(KeyError(k))
	return i
end

@inline find(s::AbstractVector, k) = findfirst(isequal(k), s)
@inline function find(s::Vector, k, ::typeof(push!))
	push!(s, k)
	return length(s)
end
@inline function find(s::Vector, k, ::typeof(insert!))
	i = findfirst(isequal(k), s)
	i ≠ nothing && return (i, false)
	push!(s, k)
	return (length(s), true)
end

@inline function find(s::Base.OneTo, k)
	@boundscheck (k ∈ s) || return nothing
	return k
end

# Associative tables ««2

# three types of associative maps are used:
# - generic Assoc: (vector of keys, vector of values)
#   append O(1), lookup O(n)
# - sorted Assoc: (sorted vector of keys, vector of values)
#   append O(n), lookup O(log n)
# - (for dense int keys): (vector of values)
#   append O(1), lookup O(1)


"""
    AbstractAssoc

 - `setindex!(assoc, val, key)`: update existing key
 - `push!(assoc, key => val)`: create new key
 - `insert!(assoc, key => val)`: update or create
 - `getindex(assoc, key)`: throws error if inexistent
 - `get(assoc, key, default)`
 - `index(keys(assoc), key, getindex)`: returns real index (or `nothing`)
 - `index(keys(assoc), key, push!)`: creates new index
 - `index(keys(assoc), key, insert!)`: returns existing or create new
 - `set!(assoc, idx, val)`: set at known index (O(1))
"""
abstract type AbstractAssoc{K,V} <: AbstractDict{K,V} end

struct Assoc{K,V} <: AbstractAssoc{K,V}
	keys::Vector{K}
	values::Vector{V}
end
# TODO: write a TwinAssoc subtype of AbstractAssoc
# replacing SparseTable

@inline (T::Type{<:AbstractAssoc})(entries::Pair...) =
	T([first.(entries)...], [last.(entries)...])

@inline Base.keys(a::AbstractAssoc) = a.keys
@inline Base.values(a::AbstractAssoc) = a.values
@inline Base.haskey(a::AbstractAssoc, u) = u ∈ keys(a)
@inline Base.length(a::AbstractAssoc) = length(keys(a))
@inline Base.convert(T::Type{<:AbstractAssoc}, ::Tuple{}) = T()

@inline Base.haskey(a::AbstractAssoc, key) = find(keys(a), key) ≠ nothing
@inline function Base.iterate(a::AbstractAssoc, s = 1)
	s > length(a) && return nothing
	return (keys(a)[s] => values(a)[s], s+1)
end
@inline function set!(a::AbstractAssoc, idx, value)
	@boundscheck idx == nothing && throw(KeyError(nothing))
	values(a)[i] = value
	return a
end

@inline function Base.getindex(a::AbstractAssoc, key)
	i = find(keys(a), key)
	i == nothing && throw(KeyError(key))
	return values(a)[i]
end
@inline function Base.get(a::AbstractAssoc, key, default)
	i = find(keys(a), key)
	i == nothing && return default
	return values(a)[i]
end
@inline function Base.setindex!(a::AbstractAssoc, value, key)
	i = find(keys(a), key)
	i == nothing && throw(KeyError(key))
	values(a)[i] = value
	return a
end
@inline function Base.push!(a::AbstractAssoc, (key, value)::Pair)
	i = find(keys(a), key, push!)
	insert!(values(a), i, value)
	return a
end
@inline function Base.insert!(a::AbstractAssoc, (key, value)::Pair)
	(i, create) = find(keys(a), key, insert!)
	if create
		insert!(values(a), i, value)
	else
		values(a)[i] = value
	end
end

"""
    Assoc{K,V}

Trivial implementation of an associative list.
Use only for *very small* lists.
"""
struct Assoc{K,V} <: AbstractAssoc{K,V}
	keys::Vector{K}
	values::Vector{V}
end
@inline Assoc{K,V}() where{K,V} = Assoc{K,V}(K[], V[])
@inline Assoc{K,V}(a::Assoc{K,V}) where{K,V} = Assoc{K,V}(a.keys, a.values)
@inline Assoc(keys::AbstractVector, values::AbstractVector) =
	Assoc{eltype(keys),eltype(values)}(keys, values)

"""
    SortedAssoc{K,V}

Implementation of an associative list with a sorted vector of keys.
Use only if keys are computed all at once (and pre-sorted).
"""
struct SortedAssoc{K,V} <: AbstractAssoc{K,V}
	keys::SortedVector{K}
	values::Vector{V}
end
@inline SortedAssoc{K,V}() where{K,V} = SortedAssoc{K,V}(K[], V[])
@inline SortedAssoc(keys::AbstractVector, values::AbstractVector) =
	SortedAssoc{eltype(keys),eltype(values)}(keys, values)

"""
    VectorAssoc{K,V}

Special case of an associative list where keys are a dense set of integers.
"""
struct VectorAssoc{K,V} <: AbstractAssoc{K,V}
	values::Vector{V}
end
@inline Base.keys(a::VectorAssoc) = Base.OneTo(length(values(a)))
VectorAssoc{K,V}(n::Integer) where{K,V<:AbstractAssoc} =
	VectorAssoc{K,V}([ V() for _ in 1:n])
@inline function Base.resize!(a::VectorAssoc, n)
	resize!(values(a), n)
	return a
end
@inline function Base.push!(a::VectorAssoc, (key, value)::Pair)
	resize!(values(a), key)
	a[key] = value
	return a
end

# SparseTable««2
struct TwinAssoc{K,V,X} <: AbstractDict{NTuple{2,K},V}
	data::X
end
make_twin_assoc(A,B) = TwinAssoc{K,V,A{K,B{K,V}}} where{K,V}
SparseTable0 = make_twin_assoc(SortedAssoc,Assoc)
SparseTable1 = make_twin_assoc(VectorAssoc,Assoc)

@inline (T::Type{<:SparseTable1})(n::Integer) =
	let D = fieldtype(T, :data)
	T(D([ valtype(D)() for _ in 1:n ]))
end
@inline (T::Type{<:SparseTable0})() = T(fieldtype(T,:data)())
@inline Base.convert(T::Type{<:SparseTable0},::Tuple{}) = T()

@inline Base.haskey(a::TwinAssoc, (i,j)) =
	haskey(a.data, i) && haskey(a.data[i], j)
@inline Base.getindex(a::TwinAssoc, i, j) =
	getindex(getindex(a.data, i), j)
@inline Base.setindex!(a::TwinAssoc, v, i, j) =
	setindex!(getindex(a.data, i), v, j)
@inline function Base.get(a::TwinAssoc, (i, j), default)
	haskey(a.data, i) || return default
	return get(a.data[i], j, default)
end
@inline function Base.push!(a::TwinAssoc, ((i, j), v)::Pair)
	if haskey(a.data, i) # FIXME: we search twice here
		push!(a.data[i], j => v)
	else
		push!(a.data, i => valtype(a.data)(j => v))
	end
	return a
end
@inline function Base.insert!(a::TwinAssoc, ((i, j), v)::Pair)
	if haskey(a.data, i)
		insert!(a.data[i], j => v)
	else
		push!(a.data, i => valtype(a.data)(j => v))
	end
	return a
end
@inline Base.length(m::TwinAssoc) = sum(length(r) for r in m.data)
@inline function Base.iterate(m::TwinAssoc) #««
	u = iterate(m.data)
	while true
		u == nothing && return u
		((i, a), s) = u # a is the alist
		v = iterate(a)
		if v ≠ nothing
			((j, y), t) = v
			return ((i,j) => y, (i, s, a, t))
		end
		u = iterate(m.data, s)
	end
end
@inline function Base.iterate(m::TwinAssoc, (i, s, a, t))
	while true
		v = iterate(a, t)
		if v ≠ nothing
			((j, y), t) = v
			return ((i,j) => y, (i, s, a, t))
		end
		while true
			u = iterate(m.data, s)
			u == nothing && return u
			((i, a), s) = u
			v = iterate(a)
			if v ≠ nothing
				((j, y), t) = v
				return ((i,j) => y, (i, s, a, t))
			end
		end
	end
end #»»
# function modify!(f, a::SortedAList, k, init)
# 	i = searchsorted(a.keys, k)
# 	if isempty(i)
# 		insert!(a.keys, first(i), k)
# 		insert!(a.values, first(i), init)
# 	else
# 		f(a.values, first(i))
# 	end
# 	return a
# end

# # SparseTable1««2
# # Associative 2-dimensional tables
# # this is different from the format used in `SparseArrays`: it is worse
# # for linear algebra (not needed here) and better for insertion.
# # this type is for when total number of elements is Θ(n), and indices are int
# # abstract type AbstractSparseTable{V,T} <: AbstractDict{NTuple{2,V},T} end««
# # # implements:
# # # rowkeys(m) - iterator for row keys
# # # data(m) - iterator for data in i-th row
# # 
# # @inline Base.length(m::AbstractSparseTable) = sum(length.(data(m)))
# # @inline Base.keys(m::AbstractSparseTable) =
# # 	[(x,y) for (x,v) in zip(rowkeys(m), data(m)) for y in keys(v)]
# # @inline function Base.iterate(m::AbstractSparseTable)
# # 	z = zip(rowkeys(m), data(m))
# # 	u = iterate(z)
# # 	while true
# # 		u == nothing && return u
# # 		((i, a), s) = u # a is the alist
# # 		v = iterate(a)
# # 		if v ≠ nothing
# # 			((j, y), t) = v
# # 			return ((i,j) => y, (i, s, a, t))
# # 		end
# # 		u = iterate(z, s)
# # 	end
# # end
# # @inline function Base.iterate(m::AbstractSparseTable, (i, s, a, t))
# # 	z = zip(rowkeys(m), data(m))
# # 	while true
# # 		v = iterate(a, t)
# # 		if v ≠ nothing
# # 			((j, y), t) = v
# # 			return ((i,j) => y, (i, s, a, t))
# # 		end
# # 		while true
# # 			u = iterate(z, s)
# # 			u == nothing && return u
# # 			((i, a), s) = u
# # 			v = iterate(a)
# # 			if v ≠ nothing
# # 				((j, y), t) = v
# # 				return ((i,j) => y, (i, s, a, t))
# # 			end
# # 		end
# # 	end
# # end»»
# """
#     SparseTable1{V,T}
# 
# Sparse 2-dimensional array with index type `V<:Integer` and value type `T`.
# This is appropriate when the total number of entries is Θ(matrix size).
# """
# struct SparseTable1{V<:Integer,T} <: AbstractDict{NTuple{2,V},T}
# 	data::Vector{Assoc{V,T}}
# 	SparseTable1{V,T}(data::AbstractVector) where{V,T} = new{V,T}(data)
# end
# 
# # empty constructor:
# @inline SparseTable1{V,T}(n::Integer) where{V,T} = 
# 	SparseTable1{V,T}([ Assoc{V,T}() for _ in 1:n ])
# @inline Base.convert(T::Type{<:SparseTable1}, n::Integer) = T(n)
# @inline Base.length(l::SparseTable1) = sum(length.(l.data))
# @inline Base.keys(l::SparseTable1) =
# 	[(i,j) for i in 1:length(l) for j in l.data[i].keys]
# @inline function Base.iterate(l::SparseTable1, s=(1,1))
# 	s[1] > length(l.data) && return nothing
# 	while s[2] > length(l.data[s[1]])
# 		s = (s[1]+1, 1)
# 		s[1] > length(l.data) && return nothing
# 	end
# 	a = l.data[s[1]]
# 	return ((s[1], a.keys[s[2]]) => a.values[s[2]], (s[1], s[2]+1))
# end
# @inline Base.haskey(l::SparseTable1, (v1, v2)) =
# 	v1 ∈ 1:length(l.data) && v2 ∈ keys(l.data[v1])
# 
# @inline Base.get(l::SparseTable1, (v1, v2), d) = get(l.data[v1], v2, d)
# @inline Base.getindex(l::SparseTable1, i) = get(l, i, nothing)
# @inline Base.setindex!(l::SparseTable1, e, (v1, v2)) =
# 	setindex!(l.data[v1], e, v2)
# @inline function Base.resize!(l::SparseTable1, n)
# 	m = length(l.data)
# 	resize!(l, n)
# 	for i in m+1:n l[i] = []; end
# 	return l
# end
# # SparseTable0««2
# """
#     SparseTable0{V,T}
# 
# Sparse 2-dimensional table with keys of type `(V,V)` and values of type `T`.
# This is appropriate when the number of entries is o(matrix size).
# """
# struct SparseTable0{V,T} <: AbstractDict{NTuple{2,V},T}
# 	# rows are stored in a SortedAssoc; this could just as well be a
# 	# SortedDict
# 	data::SortedAssoc{V,Assoc{V,T}}
# end
# @inline SparseTable0{V,T}() where{V,T} =
# 	SparseTable0{V,T}(SortedAssoc{V,Assoc{V,T}}())
# @inline SparseTable0{V,T}(::Nothing) where{V,T} = SparseTable0{V,T}()
# @inline Base.convert(T::Type{<:SparseTable0}, ::Tuple{}) = T()
# @inline Base.length(m::SparseTable0) = sum(length(x) for x in m.data)
# @inline Base.getindex(m::SparseTable0, (i,j)) =
# 	getindex(getindex(m.data, i), j)
# @inline function Base.get(m::SparseTable0, (i,j), d)
# 	row = get(m.data, i, nothing)
# 	row == nothing && return d
# 	return get(row, j, d)
# end
# @inline Base.setindex!(m::SparseTable0, y, (i,j)) =
# 	modify!(m.data, i, Assoc(j=>y)) do v, k; v[k][j] = y end
# @inline function Base.iterate(m::SparseTable0, s=(1,1))
# 	s[1] > length(m.data) && return nothing
# 	while s[2] > length(m.data.values[s[1]])
# 		s = (s[1]+1, 1)
# 		s[1] > length(m.data) && return nothing
# 	end
# 	a = m.data.values[s[1]]
# 	return ((m.data.keys[s[1]], a.keys[s[2]]) => a.values[s[2]], (s[1], s[2]+1))
# end
# 
# LazyMap ««2
struct LazyMap{Y,X,F} <: AbstractVector{Y}
	f::F
	source::X
end
@inline Base.size(l::LazyMap) = size(l.source)
@inline Base.getindex(l::LazyMap, i::Integer) = l.f(l.source[i])

LazyMap{Y}(f, v) where{Y} = LazyMap{Y,typeof(v),typeof(f)}(f, v)
lazymap(f,v) = LazyMap{Base.return_types(f,Tuple{eltype(v)})[1]}(f, v)

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
	minlevel = SortedAssoc{eltype(u),NTuple{2,eltype(u)}}()
	for i in 1:length(u.parent)
		r = u.parent[i]; l = u.level[i]
		m = get(minlevel, r, nothing)
		if m == nothing
			push!(minlevel, r => (i, l))
		elseif m[2] > l
			minlevel[r] = (i, l)
		end
# 		if (m == nothing) || (m[2] > l)
# 			push!(minlevel, r => (i, l))
# 		end
	end
	for i in 1:length(u.parent)
		r = u.parent[i]; l = u.level[i]; m = minlevel[r]
		u.parent[i] = m[1]
		u.level[i] = l - m[2]
	end
	return u
end

# level[i2] = level[i1] + δ
function connect!(u::LevelStructure, i1, i2, δ)
	(r1, d1) = root_and_level(u, i1) # l(i1) = l(r1) + d1
	(r2, d2) = root_and_level(u, i2) # l(i2) = l(r2) + d2
	# we now want l(i2) = l(i1) + δ, hence
	# l(r2) + d2 = l(r1) + d1 + δ
	# l(r2) = l(r1) + d1 - d2 + δ
	if r1 == r2
		@assert d2 == d1 + δ "connect($i1, $i2, $δ): faces are already connected and multiplicity should be $(d1+δ) (it is $d2). Faces are likely crossed around an edge."
		return
	end
	if u.treesize[r1] < u.treesize[r2]
		u.parent[r1] = u.parent[r2]; u.treesize[r2]+= u.treesize[r1]
		u.level[r1] = d2 - d1 -δ
	else
		u.parent[r2] = u.parent[r1]; u.treesize[r1]+= u.treesize[r2]
		u.level[r2] = d1 + d2 + δ
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

# function lowest_representatives(n, relation)
# 	rep = collect(1:n)
# 	for (i, j) in relation
# 		if rep[i] < rep[j]
# 			rep[j] = rep[i]
# 		else
# 			rep[i] = rep[j]
# 		end
# 	end
# 	return rep
# end

function lowest_representatives(relation)
	elements = unique!(sort!([r[i] for r in relation for i in 1:2]))
	return (elements, lowest_representatives(relation, elements))
end

# enumerate equivalence classes««2
function equivalence_classes(n, relation;
	representatives = lowest_representatives(relation, 1:n))
	g = SparseTable1{Int,Nothing}(n)
	for (i,j) in pairs(representatives)
		push!(g, (j,i) => nothing)
	end
	return [ keys(u) for u in g.data.values if !isempty(u) ]
end

# Equivalence structure for arbitrary (sortable) objects««2
"""
    EquivalenceStructure

Equivalence structure for arbitrary (sortable) objects

Useful methods:
 - `elements(s)`: returns uniquesorted list of all elements
 - `equivalent(s, x,y)`: true iff x, y are equivalent
 - `class(s, x)`: returns all elements equivlaent to x
 - `representative(s,x)`: first(class(s,x))
"""
# from an element, find an equivalence class identifier (integer)
# from an equivalence class id, find (sorted) list of equiv. elements
struct EquivalenceStructure{T,A}
	representative::A # SortedAssoc{T,Int}, etc.
	class::Vector{Vector{Int}}
	@inline EquivalenceStructure(rep::AbstractAssoc, class) =
		new{keytype(rep),typeof(rep)}(rep, class)
end
@inline function index(s::EquivalenceStructure, x)
	i = searchsortedfirst(keys(s.representative), x)
	i > length(keys(s.representative)) && throw(KeyError(x))
	return i
end
@inline class(s::EquivalenceStructure, x) =
	keys(s.representative)[classidx(s, x)]
@inline classidx(s::EquivalenceStructure, x) =
	s.class[s.representative[x]]
@inline equivalent(s::EquivalenceStructure, x, y) =
	s.representative[x] == s.representative[y]
@inline classes(s::EquivalenceStructure) =
	(s.representative.keys[i] for i in s.class)
@inline elements(s::EquivalenceStructure) = keys(s.representative)
@inline representative(s::EquivalenceStructure, x) = first(class(s,x))
@inline function Base.rem(x, s::EquivalenceStructure)
	i = find(keys(s.representative), x)
	i == nothing && return x
	return keys(s.representative)[first(s.class[values(s.representative)[i]])]
end

function equivalence_structure(relation)
	# first compute names
	rel_flat = [ r[i] for r in relation for i in 1:2 ]
	isempty(relation) && return EquivalenceStructure(
		SortedAssoc{eltype(eltype(relation)), Int}(),
		Vector{Int}[])

	(names, rel_idx) = uniquenames(rel_flat)
	# `names` is the unique list of names
	# `rel_idx` is the renamed equivalence relation (as small integers)
	n = length(names)
	r = length(rel_idx)
	# compute transitive closure and representatives of each class:
	uf = UnionFind(n)
	for i in 1:2:length(rel_idx)
		union!(uf, rel_idx[i], rel_idx[i+1])
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
	return EquivalenceStructure(SortedAssoc(names, rep_idx),class)
end
# half-edge mesh««1
# half-edge data structure ««2
struct HalfEdgeMesh{H,V,P}
	opposite::Vector{H}
	destination::Vector{V}
	edgefrom::Vector{H} # indexed by VertexId
	points::Vector{P}
	plane::Vector{Tuple{Int8,P}} # normalized plane equations
	# FIXME: convert to structure-of-array?
end

halfedge_type(::Type{<:HalfEdgeMesh{H}}) where{H,V,P} = H
halfedge_type(s::HalfEdgeMesh) = halfedge_type(typeof(s))
@inline face_type(s) = halfedge_type(s)
vertex_type(::Type{<:HalfEdgeMesh{H,V}}) where{H,V,P} = V
vertex_type(s::HalfEdgeMesh) = vertex_type(typeof(s))
point_type(::Type{<:HalfEdgeMesh{H,V,P}}) where{H,V,P} = P
point_type(s::HalfEdgeMesh) = point_type(typeof(s))
dir_type(::HalfEdgeMesh) = Int8

@inline HalfEdgeMesh{H,V,P}(::UndefInitializer, points, nf::Integer
	) where{H,V,P} = HalfEdgeMesh{H,V,P}(
	Vector{H}(undef, 3nf),
	Vector{V}(undef, 3nf),
	Vector{H}(undef, length(points)),
	points,
	Vector{Tuple{Int8,P}}(undef, nf),
	)

@inline HalfEdgeMesh(::UndefInitializer, points, args...) =
	HalfEdgeMesh{HalfEdgeId,VertexId,eltype(points)}(undef, points, args...)

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
function project2d(axis, vec)
	# orientation-preserving projection:
	axis == 1 && return (vec[2], vec[3])
	axis == 2 && return (vec[3], vec[1])
	axis == 3 && return (vec[1], vec[2])
	axis ==-1 && return (vec[3], vec[2])
	axis ==-2 && return (vec[1], vec[3])
	@assert axis == -3
	             return (vec[2], vec[1])
end
function normalized_plane_eq(p1,p2,p3)#««
	# normalized plane equation through three points
	direction = cross(p2-p1, p3-p1)
	@inbounds u1 = direction[1]
	@inbounds u2 = direction[2]
	@inbounds u3 = direction[3]
  @inbounds a1 = abs(u1)
	@inbounds a2 = abs(u2)
	@inbounds a3 = abs(u3)
	P = typeof(p1)
	if a2 < a1
		a1 < a3 && @goto max3
		v = (u2/u1, u3/u1)
		u1 > 0 && return ( 1, P(v[1], v[2], p1[1]+v[1]*p1[2]+v[2]*p1[3]))
		          return (-1, P(v[1], v[2], p1[1]+v[1]*p1[2]+v[2]*p1[3]))
	elseif a2 > a3
		v = (u1/u2, u3/u2)
		u2 > 0 && return ( 2, P(v[1], v[2], v[1]*p1[1]+p1[2]+v[2]*p1[3]))
		          return (-2, P(v[1], v[2], v[1]*p1[1]+p1[2]+v[2]*p1[3]))
	end; @label max3
		v = (u1/u3, u2/u3)
		u3 > 0 && return ( 3, P(v[1], v[2], v[1]*p1[1]+v[2]*p1[2]+p1[3]))
		          return (-3, P(v[1], v[2], v[1]*p1[1]+v[2]*p1[2]+p1[3]))
end#»»
function mark_edge!(s::HalfEdgeMesh, edge_loc, k, v1, v2)#««
	println("   mark edge ($v1, $v2) at position $k:", get(edge_loc, (2,4), :NOTFOUND))
	s.destination[k] = v2
	s.edgefrom[v1] = k
# 	haskey(edge_loc, (v1,v2)) && error("haskey: $((v1,v2))")
	insert!(edge_loc, (v1,v2) => k)
	j = get(edge_loc, (v2,v1), nothing)
	if j ≠ nothing
		s.opposite[j] = k
		s.opposite[k] = j
	end
end#»»
function mark_face!(s::HalfEdgeMesh, edge_loc, i, v1, v2, v3)#««
	println("  mark_face($v1,$v2,$v3)")
	mark_edge!(s, edge_loc, 3i-2, v3, v1)
	mark_edge!(s, edge_loc, 3i-1, v1, v2)
	mark_edge!(s, edge_loc, 3i  , v2, v3)
	s.plane[i] = normalized_plane_eq(point(s,v1), point(s,v2), point(s,v3))
end#»»

function HalfEdgeMesh(points, faces)
	nv = length(points)
	@assert all(length.(faces) .== 3)
	s = HalfEdgeMesh(undef, points, length(faces))
	edge_loc = SparseTable1{halfedge_type(s),vertex_type(s)}(nv)

	for (i, f) in pairs(faces)
		(v1,v2,v3) = f
		# (v3v1), (v1v2), (v3v2)
		mark_face!(s, edge_loc, i, v1, v2, v3)
	end
	return s
end

# simple accessors ««2
@inline nfaces(s::HalfEdgeMesh) = length(s.opposite) ÷ 3
@inline nvertices(s::HalfEdgeMesh) = length(s.edgefrom)
@inline next(::HalfEdgeMesh, n) = n+1 - 3*(n%3==0)
@inline prev(::HalfEdgeMesh, n) = n-1 + 3*(n%3==1)
@inline opposite(s::HalfEdgeMesh, e) = s.opposite[e]
@inline opposite!(s::HalfEdgeMesh, e, x) = s.opposite[e] = x
@inline destination(s::HalfEdgeMesh, e) = s.destination[e]
@inline destination!(s::HalfEdgeMesh, e, v) = s.destination[e] = v
@inline edgefrom(s::HalfEdgeMesh, v) = s.edgefrom[v]
@inline edgefrom!(s::HalfEdgeMesh, v, e) = s.edgefrom[v] = e
@inline points(s::HalfEdgeMesh) = s.points
@inline point(s::HalfEdgeMesh, i) = s.points[i]
@inline planes(s::HalfEdgeMesh) = s.plane
@inline plane(s::HalfEdgeMesh, i) = s.plane[i]

@inline adjacent_faces(s::HalfEdgeMesh, i, j) =
	fld1(opposite(s, 3i-2), 3) == j ||
	fld1(opposite(s, 3i-1), 3) == j ||
	fld1(opposite(s, 3i  ), 3) == j

@inline halfedge(s::HalfEdgeMesh, e) =
	(destination(s, prev(s, e)), destination(s, e))
@inline opposite_vertex(s::HalfEdgeMesh, e) =
	destination(s, next(s, e))
@inline face_edge(s::HalfEdgeMesh, f, j) = halfedge(s, 3*f-3+j)
@inline face_vertex(s::HalfEdgeMesh, f, j) = destination(s, 3*f-3+j)
@inline regular_edge(s::HalfEdgeMesh, e) = opposite(s, opposite(s, e)) == e

@inline adjacent_face(s::HalfEdgeMesh, f, j) =
	fld1(opposite(s, 3*f-3+j), 3)
@inline adjacent_faces(s::HalfEdgeMesh, f) =
	map(j->adjacent_face(s, f, j), (1,2,3))

function Base.resize!(s::HalfEdgeMesh, nv, nf)
	resize_points!(s, nv)
	resize_faces!(s, nf)
	s
end

function Base.reverse(s::HalfEdgeMesh)
	r = HalfEdgeMesh{halfedge_type(s), vertex_type(s), point_type(s)}(
		undef, points(s), nfaces(s))
	# reverse all half-edges:
	# (31)(12)(23) to (21)(13)(32)
	#
	# edgefrom: 3f-2 <-> 3f, and 3f-1 remains in place
	# destination: 3f-2 in place, 3f-1 ↔ 3f
	# opposite: 3f-2  in place, 3f-1 ↔ 3f
	@inline newedgeno(e) = let m = mod(e,3)
		(m == 1) ? e + 1 : (m == 2) ? e - 1 : e end
	for f in 1:nfaces(s)
		r.destination[3*f  ] = s.destination[3*f-1]
		r.destination[3*f-1] = s.destination[3*f  ]
		r.destination[3*f-2] = s.destination[3*f-2]
		r.opposite[3*f-2] = newedgeno(s.opposite[3*f-1])
		r.opposite[3*f-1] = newedgeno(s.opposite[3*f-2])
		r.opposite[3*f  ] = newedgeno(s.opposite[3*f  ])
		(d, v) = plane(s, f); r.plane[f] = (-d, v)
	end
	for v in 1:nvertices(s)
		e = s.edgefrom[v]; m = mod(e,3)
		r.edgefrom[v] = (m == 1) ? e+2 : (m == 2) ? e : e-2
	end
	return r
end

# faces iterator««2
struct HalfEdgeFaceIterator{H,V} <: AbstractVector{NTuple{3,V}}
	mesh::HalfEdgeMesh{H,V}
end

@inline Base.size(h::HalfEdgeFaceIterator) = (nfaces(h.mesh),)
@inline Base.getindex(h::HalfEdgeFaceIterator, i::Integer) =
	(destination(h.mesh,3i-2), destination(h.mesh,3i-1), destination(h.mesh,3i))
@inline faces(s::HalfEdgeMesh) =
	HalfEdgeFaceIterator{halfedge_type(s),vertex_type(s)}(s)
@inline face(s::HalfEdgeMesh, f) = faces(s)[f]

struct HalfEdgeTriangleIterator{H,V,P} <: AbstractVector{NTuple{3,P}}
	mesh::HalfEdgeMesh{H,V,P}
end

@inline Base.size(h::HalfEdgeTriangleIterator) = (nfaces(h.mesh),)
@inline Base.getindex(h::HalfEdgeTriangleIterator, i::Integer) = (
	point(h.mesh,destination(h.mesh, 3i-2)),
	point(h.mesh,destination(h.mesh, 3i-1)),
	point(h.mesh,destination(h.mesh, 3i  )),
	)
@inline triangles(s::HalfEdgeMesh) =
	HalfEdgeTriangleIterator{halfedge_type(s),vertex_type(s),point_type(s)}(s)
@inline triangle(s::HalfEdgeMesh, f) = triangles(s)[f]

# vertex iterator««2
struct HalfEdgeVertexEdges{H,V}
	mesh::HalfEdgeMesh{H,V}
	start::H
end
@inline neighbours(s::HalfEdgeMesh, v) =
	HalfEdgeVertexEdges{halfedge_type(s),vertex_type(s)}(s, edgefrom(s,v))
@inline Base.eltype(::Type{HalfEdgeVertexEdges{H,V}}) where{H,V} = Pair{H,V}
@inline Base.iterate(it::HalfEdgeVertexEdges) =
	(it.start => destination(it.mesh, it.start), it.start)
@inline Base.iterate(it::HalfEdgeVertexEdges, s) =
	let e = next(it.mesh, opposite(it.mesh, s))
	e == it.start && return nothing
	return (e => destination(it.mesh, e), e)
end
Base.IteratorSize(::HalfEdgeVertexEdges) = Base.SizeUnknown()
function edge(s::HalfEdgeMesh, (v1, v2))
	for (he, dest) in neighbours(s, v1)
		dest == v2 && return he
	end
	throw(KeyError((v1,v2))); return edgefrom(v1) # type-stability
end

struct HalfEdgeRadialIterator{H,V}
	mesh::HalfEdgeMesh{H,V}
	start::H
end
@inline radial_loop(s::HalfEdgeMesh, e) =
	HalfEdgeRadialIterator{halfedge_type(s),vertex_type(s)}(s, e)
@inline Base.iterate(it::HalfEdgeRadialIterator) = (it.start, it.start)
function Base.iterate(it::HalfEdgeRadialIterator, s)
	e = opposite(it.mesh, s)
	e == it.start && return nothing
	return (e, e)
end
@inline Base.IteratorSize(::HalfEdgeRadialIterator) = Base.SizeUnknown()
@inline Base.eltype(::HalfEdgeRadialIterator{H,V}) where{H,V} = H

# find coplanar faces««2
# returns the equivalence relation
function coplanar_faces(s::HalfEdgeMesh; ε = 0)
	@inline loc((dir, v)) = SA[abs(dir), v...]
	boxes = [ BBox(p, p .+ ε) for p in loc.(planes(s)) ]
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
# structure extension««2
function concatenate(slist::HalfEdgeMesh...)
	r = HalfEdgeMesh{halfedge_type(first(slist)),
		vertex_type(first(slist)),point_type(first(slist))}(undef,
		vcat(points.(slist)...),
		sum(nfaces.(slist)))
	eoffset = 0
	voffset = 0
	foffset = 0
	for s in slist
		(nv, ne, nf) = length(s.edgefrom), length(s.opposite), length(s.plane)
		r.opposite[eoffset+1:eoffset+ne] .= s.opposite .+ eoffset
		r.destination[eoffset+1:eoffset+ne] .= s.destination .+ voffset
		r.edgefrom[voffset+1:voffset+nv] .= s.edgefrom .+ eoffset
		r.plane[foffset+1:foffset+nf] .= s.plane
		eoffset += ne
		voffset += nv
		foffset += nf
	end
	return r
end
function extend_points!(s::HalfEdgeMesh, points)
	nv1 = nvertices(s)
	nv2 = length(points)
	resize!(s.edgefrom, nv2)
	resize!(s.points, nv2)
	s.points[nv1+1:nv2] .= points[nv1+1:nv2]
	s
end
function resize_faces!(s::HalfEdgeMesh, nf)
	resize!(s.opposite, 3nf)
	resize!(s.destination, 3nf)
	resize!(s.plane, nf)
	s
end
# remove points««2
"""
    shift_points!(s, vmap, start, stop, offset)

Deletes points in `s`, replacing interval `start:stop`
by `start+offset:stop+offset`.
`vmap` is a table holding the renaming of points (for `destination`).
"""
function shift_points!(s::HalfEdgeMesh, vmap, start, stop, offset)
# 		println("points $start:$stop shifted by $offset")
	for i in start:stop
		s.points[i] = s.points[i+offset]
		s.edgefrom[i] = s.edgefrom[i+offset]
		vmap[i+offset] = i
	end
end


function deduplicate_vertices!(s::HalfEdgeMesh, same_points)
	nv = nvertices(s)
	vmap = Vector{Int}(undef, nv)

	replace = keys(same_points.representative)[first.(
		same_points.class[values(same_points.representative)])]
	@assert replace == [x % same_points for x in elements(same_points)]
	start = 1; offset = 0
	for (x,y) in zip(elements(same_points), replace)
		x == y && continue
		current = x-offset
		shift_points!(s, vmap, start, current-1, offset)
		start = current; offset+= 1
	end
	shift_points!(s, vmap, start, nv-offset, offset)
	# replace points by correct values:
	vmap[elements(same_points)] .= vmap[replace]

	s.destination .= view(vmap, s.destination)
	resize!(s.points, nv-offset)
	resize!(s.edgefrom, nv-offset)
	return s
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
	for i in 1:length(fkept)
		fkept[i] || continue
		s.plane[i] = s.plane[fld(emap[3*i], 3)]
	end
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
	resize!(s.opposite, last(emap))
	resize!(s.destination, last(emap))
	resize!(s.edgefrom, last(vmap))
	resize!(s.points, last(vmap))
	resize!(s.plane, last(emap) ÷ 3)
	return s
end

# split faces!««2
"""
    split_faces

Replaces faces in `s`, as indicated by iterator `fsplit`
(as (face number) => (replacement triangles)).
"""
function split_faces(s::HalfEdgeMesh, points, fsplit)
	s1 = deepcopy(s)
	# localize new edges only — this is O(√n), hence SparseTable0:
	edge_loc = SparseTable0{halfedge_type(s1),vertex_type(s1)}()
	extend_points!(s1, points)
	# n is number of face being written
	n = nfaces(s1)
	println("\e[1msplit_faces\e[m")
	for (f, tri) in fsplit
		length(tri) ≤ 1 && continue # skipped face, or trivial retriangulation

		println("\e[34m$f $(face(s,f)) => $tri\e[m")
		resize_faces!(s1, nfaces(s1) + length(tri) - 1)
		((a1,b1,c1), state) = iterate(tri)
		(v1,v2,v3) = face(s1, f)
		# before overwriting this face, we save the location of opposed half-edges:
		# and save location of half-edges from v1, v2, v3:
		println("h-e $(3*f-2) is $(halfedge(s1, 3*f-2))")
		println("h-e $(3*f-1) is $(halfedge(s1, 3*f-1))")
		println("h-e $(3*f-0) is $(halfedge(s1, 3*f-0))")
		e = opposite(s1, 3*f-2); insert!(edge_loc, (v1,v3) => e); edgefrom!(s1,v1,e)
		e = opposite(s1, 3*f-1); insert!(edge_loc, (v2,v1) => e); edgefrom!(s1,v2,e)
		e = opposite(s1, 3*f  ); insert!(edge_loc, (v3,v2) => e); edgefrom!(s1,v3,e)
		mark_face!(s1, edge_loc, f, a1, b1, c1)
		# append new faces
		while true
			u = iterate(tri, state)
			u == nothing && break
			((a,b,c), state) = u
			n+= 1
			mark_face!(s1, edge_loc, n, a, b, c)
		end
	end
	return s1
end
#
# self-intersection««1
# self_intersect««2
# indexed by pairs of vertices
# other possibility would be to index by half-edge number
struct EdgeInserts{H,V,P,Z}
	mesh::HalfEdgeMesh{H,V,P}
	vertices::SparseTable0{V,Vector{V}}
	# i is the index of the useful coordinate for sorting
	sort::SparseTable0{V,@NamedTuple{dir::Int8,z::Vector{Z}}}
	@inline EdgeInserts{H,V,P,Z}(mesh, n::Integer) where{H,V,P,Z} =
		new{H,V,P,Z}(mesh, (), ())
	@inline EdgeInserts(mesh::HalfEdgeMesh) =
		EdgeInserts{halfedge_type(mesh),vertex_type(mesh),point_type(mesh),
			eltype(point_type(mesh))}(mesh, nvertices(mesh))
end
function Base.insert!(s::EdgeInserts, i, allpoints, p, ε = 0)
# i: half-edge index
# returns id of vertex
	e = extrema(halfedge(s.mesh, i))
	if !haskey(s.vertices, e)
		# no previous point on this edge: insertion is trivial
		v = point(s.mesh, e[2]) - point(s.mesh, e[1])
		k = findmax(abs.(v))[2]
		push!(s.sort, e => (dir = dir_type(s.mesh)(v[k] > 0 ? k : -k), z = [p[k]]))
		push!(allpoints, p)
		newvertex = vertex_type(s.mesh)(length(allpoints))
		push!(s.vertices, e => [newvertex])
		return newvertex
	else
		t = s.sort[e]
		z = p[abs(t.dir)]*sign(t.dir)
		ins = searchsortedfirst(t.z, z-ε)
		if (ins ≤ length(t.z) && t.z[ins] ≤ z+ε)
			return s.vertices[e][ins]
		else
			push!(allpoints, p)
			newvertex = vertex_type(s.mesh)(length(allpoints))
			insert!(s.vertices[e], ins, newvertex)
			insert!(t.z, ins, z)
			return newvertex
		end
	end
end

function register_point!(s::HalfEdgeMesh, allpoints, in_face, in_edge,#««
		i, t, p, ε = 0)
	# i is face number, t is point type, p is geometric point
	# returns the point number
	TI = TriangleIntersections

	TI.isvertex(t) && return face(s, i)[TI.index(t, TI.isvertex)]

# 	push!(allpoints, p)
# 	vindex = vertex_type(s)(nvertices(s) + length(allpoints))

	if TI.isedge(t)
		# HalfEdgeMesh has edges in order: edge31, edge12, edge23
		# index(edge31<<2) = index(edge23) = 1 fixes this:
		k = TI.index(t<<2, TI.isedge)

		return insert!(in_edge, 3*i-3+k, allpoints, p, ε)
	end
	# this is an iterior point:
	if !haskey(in_face, i)
		push!(allpoints, p)
		newvertex = vertex_type(s)(length(allpoints))
		push!(in_face, i => [newvertex])
		return newvertex
	end
	# look for a previous, matched point:
	for j in in_face[i]
		norm(allpoints[j] - p, Inf) ≤ ε && return j
	end
	push!(allpoints, p)
	newvertex = vertex_type(s)(length(allpoints))
	push!(in_face[i], p)
	return newvertex
end#»»

"""
    self_intersect(s)

Returns `(points, in_face, in_edge, same_edges)` describing the
self-intersection graph of `s`.
These values are *not* filtered for repetition.
"""
function self_intersect(s::HalfEdgeMesh, ε=0)#««
	T = eltype(point_type(s))
# 	edge_points = [ Int[] for _ in edges(s) ]
# 	edge_coords = [ T[] for _ in edges(s) ]

	boxes = [ boundingbox(t...) for t in triangles(s) ]
	# plist: list of new points
	# elist: list of new edges
	# elist[i] always connects points (2i-1) and (2i)
	# these will be filtered later to be made unique
	pts = copy(points(s))
	same_points = NTuple{2,vertex_type(s)}[]
	in_face = Assoc{face_type(s),Vector{vertex_type(s)}}()
	in_edge = EdgeInserts(s)
	same_edges = NTuple{2,NTuple{2,vertex_type(s)}}[]

	for (i1, i2) in SpatialSorting.intersections(boxes)
		adjacent_faces(s, i1, i2) && continue
		# this can fail only for a flat triangle:
		it = TriangleIntersections.inter(triangle(s, i1), triangle(s, i2), ε)
		isempty(it) && continue
		
		# create points as needed, and store their index:
		vindex = MMatrix{6,2,vertex_type(s)}(undef)
		for (i, (p, (t1, t2))) in pairs(it)
			idx1 = register_point!(s, pts, in_face, in_edge, i1, t1, p, ε)
			idx2 = register_point!(s, pts, in_face, in_edge, i2, t2, p, ε)
			idx1 ≠ idx2 && push!(same_points, (idx1, idx2))
			vindex[i,1] = idx1; vindex[i,2] = idx2
		end
		v1 = vindex[1,:]
		for i in 2:length(it)
			v2 = vindex[i,:]
			for j in 1:2
				if v1[j] < v2[j]
					push!(same_edges, ((v1[j], v2[j]), (v1[3-j], v2[3-j])))
				else
					push!(same_edges, ((v2[j], v1[j]), (v2[3-j], v1[3-j])))
				end
			end
		end
	end
	return (
		points=pts, in_face=in_face,
		in_edge=in_edge.vertices,
		same_points=equivalence_structure(same_points),
		same_edges=equivalence_structure(same_edges),
		)
end#»»

# subtriangulation««1
# faces_to_triangulate««2
"""
    faces_to_triangulate(s)

Returns `(faces_todo, faces_todo_vertices)`, such that
`faces_todo` is a list of all faces to retriangulate,
and `faces_todo_vertices[i]` is the set of all vertices
included inside `faces_todo[i]`.
"""
function faces_to_triangulate(s, in_face, in_edge)#««3
	# detect which faces we actually need to process
	# faces_todo: list of all faces to subdivide
	# faces_todo_vertices[i]: list of all vertices included in faces_todo[i]
	faces_todo = face_type(s)[]
	faces_todo_vertices = Vector{vertex_type(s)}[]
	#  - those bordering an (old) edge in `inserts`:
	for (e, vlist) in in_edge
		he = edge(s, e)
		@assert he isa halfedge_type(s)
		f1 = fld1(he, 3)
		push!(faces_todo, f1)
		push!(faces_todo_vertices, [vlist...; face(s, f1)...])
		f2 = fld1(opposite(s, he), 3)
		push!(faces_todo, f2)
		push!(faces_todo_vertices, [vlist...; face(s, f2)...])
	end
	for (f, vlist) in in_face
		push!(faces_todo, f)
		push!(faces_todo_vertices, copy(vlist))
	end
	isempty(faces_todo) && return Assoc(faces_todo,
		[SortedVector(sort!(x)) for x in faces_todo_vertices])
	# sort and unique faces_todo, and combine corresponding entries
	# in faces_todo_vertices
	p = sortperm(faces_todo)
	k = faces_todo[p[1]]
	newfaces_todo = [k]
	newfaces_todo_vertices = [faces_todo_vertices[p[1]]]
	for i in p[2:end]
		if k == faces_todo[i]
			push!(last(newfaces_todo_vertices), faces_todo_vertices[i]...)
		else
			k = faces_todo[i]
			push!(newfaces_todo, k)
			push!(newfaces_todo_vertices, copy(faces_todo_vertices[i]))
		end
	end
	for l in newfaces_todo_vertices
		unique!(sort!(l))
	end
	return Assoc(newfaces_todo, newfaces_todo_vertices)
end
# project_and_triangulate ««2
function project_and_triangulate(points, direction, vlist, elist, same_points)
	# this has duplicate points (for matched edges),
	# which libtriangle does not like => deduplicate points:
	uniquevertices = [ v % same_points for v in vlist ]
	unique!(sort!(uniquevertices))
	uniqueedges = [ minmax(e1 % same_points, e2 % same_points)
		for (e1, e2) in elist ]
	unique!(sort!(uniqueedges))
	# build matrix of coordinates according to `vlist`:
	vmat = Matrix{Float64}(undef, length(uniquevertices), 2)
	if direction == 1
		vmat[:,1] .= (points[v][2] for v in uniquevertices)
		vmat[:,2] .= (points[v][3] for v in uniquevertices)
	elseif direction ==-1
		vmat[:,1] .= (points[v][3] for v in uniquevertices)
		vmat[:,2] .= (points[v][2] for v in uniquevertices)
	elseif direction == 2
		vmat[:,1] .= (points[v][3] for v in uniquevertices)
		vmat[:,2] .= (points[v][1] for v in uniquevertices)
	elseif direction ==-2
		vmat[:,1] .= (points[v][1] for v in uniquevertices)
		vmat[:,2] .= (points[v][3] for v in uniquevertices)
	elseif direction == 3
		vmat[:,1] .= (points[v][1] for v in uniquevertices)
		vmat[:,2] .= (points[v][2] for v in uniquevertices)
	else @assert direction ==-3
		vmat[:,1] .= (points[v][2] for v in uniquevertices)
		vmat[:,2] .= (points[v][1] for v in uniquevertices)
	end
	emat = Matrix{Int}(undef, length(uniqueedges), 2)
	emat[:,1] .= collect(e[1] for e in uniqueedges)
	emat[:,2] .= collect(e[2] for e in uniqueedges)
# 	println("triangulate: $vmat $uniquevertices $emat")
	return LibTriangle.constrained_triangulation(vmat, uniquevertices, emat)

end

"""
    triangles_in

Returns those triangles in `triangulation` (as triple of vertices)
which have all their vertices in (sorted list) `elist`,
modulo equivalence relation `same_points`.
"""
function triangles_in(triangulation, vlist, same_points)
	# build associative list in reverse order, to identify points in
	# `vlist` by their lowest representative:
	vlist2 = [ x % same_points for x in vlist ]
	perm = sortperm(vlist2)
	assoc = SortedAssoc(vlist2[perm], vlist[perm])
	
	ret = NTuple{3,eltype(eltype(triangulation))}[]
	vlist2s = SortedVector(sort(vlist2))
	for tri in triangulation
		(v1, v2, v3) = tri
		v1 ∉ vlist2s && continue
		v2 ∉ vlist2s && continue
		v3 ∉ vlist2s && continue
		push!(ret, (assoc[v1], assoc[v2], assoc[v3]))
	end
	return ret
end

# subtriangulate««2
function subtriangulate(s::HalfEdgeMesh, ε=0)
	println("subtriangulate(v=$(nvertices(s)), f=$(nfaces(s)))")
	si = self_intersect(s, ε)
	# determine clusters of coplanar faces:
	# (using compact type for equivalence structure since O(n) faces...)
	coplanar_rel = coplanar_faces(s; ε)
	coplanar_rep = lowest_representatives(coplanar_rel, 1:nfaces(s))
	clusters = equivalence_classes(nfaces(s), coplanar_rel;
		representatives = coplanar_rep)

	faces_todo = faces_to_triangulate(s, si.in_face, si.in_edge)
	faces_done = falses(length(faces_todo))
	faces_tri = [ NTuple{3,vertex_type(s)}[] for _ in faces_todo ]

	# iterate over all faces todo
	for (f, done) in zip(keys(faces_todo), faces_done)
		done && continue
		lowrep = coplanar_rep[f]
		cluster_index =
			searchsortedfirst(LazyMap{face_type(s)}(first, clusters), lowrep)
		current_cluster = clusters[cluster_index]

		# mark these faces as done
		start = 1; stop = length(faces_todo)
		for f in current_cluster
			j = searchsortedfirst(keys(faces_todo), f, start,stop,Base.Order.Forward)
			j > stop && break
			faces_done[j] = true
		end

		# extract all vertices in this cluster:
		# FIXME this could also be accessed using values(faces_todo)...
		clustervertices = sizehint!(vertex_type(s)[], 3*length(current_cluster))
		clusteredges = sizehint!(NTuple{2,vertex_type(s)}[],
				2*length(current_cluster))

		# extract all kept edges in this cluster:
		for f in current_cluster
			push!(clustervertices, face(s, f)...)
			haskey(si.in_face, f) && push!(clustervertices, si.in_face[f]...)
			for i in 1:3
				g = adjacent_face(s, f, i)
				g < f && !isempty(searchsorted(current_cluster, g)) && continue
				e = extrema(face_edge(s, f, i))
				if haskey(si.in_edge, e)
					vlist = si.in_edge[e]
					push!(clustervertices, vlist...)
					push!(clusteredges, minmax(e[1], first(vlist)))
					push!(clusteredges, minmax(e[2], last(vlist)))
					for i in 2:length(vlist)
						push!(clusteredges, minmax(vlist[i-1], vlist[i]))
					end
				else
					push!(clusteredges, e)
				end
			end
		end
		unique!(sort!(clustervertices))
		length(clustervertices) == 3 && continue

		unique!(sort!(clusteredges))

		# project and triangulate
		newtriangles = project_and_triangulate(si.points,
			abs(plane(s,f)[1]), clustervertices, clusteredges, si.same_points)
# 		println("newtriangles for $current_cluster=$newtriangles")

		# extract triangulation for each face
		for f in current_cluster
			j = searchsortedfirst(keys(faces_todo), f)
			j > length(faces_todo) && continue
			keys(faces_todo)[j] ≠ f && continue
			vlist = values(faces_todo)[j]
			tri = triangles_in(newtriangles, vlist, si.same_points)
			faces_tri[j] = (plane(s,f)[1] > 0) ? tri : reverse.(tri)
# 			println(" replacing face $j=\e[34m$(face(s,j)) ($vlist)\e[m => $(faces_tri[j])")
		end
	end
	# apply refinement computed above
	s1 = split_faces(s, si.points, zip(keys(faces_todo), faces_tri))

	# build radial cycles for singular half-edges:
	global S = s1
	global SI = si
	he1 = [ edge(s1, e) for e in elements(si.same_edges) ]
	he2 = [ opposite(s1, e) for e in he1 ]
	for c in si.same_edges.class
		c1 = he1[c]; c2 = he2[c]
		# interleave the two lists:
		for i in 1:length(c)-1
			opposite!(s1, c1[i], c2[i])
			opposite!(s1, c2[i], c1[i+1])
		end
		opposite!(s1, last(c1), last(c2))
		opposite!(s1, last(c2), first(c1))
	end
	deduplicate_vertices!(s1, si.same_points)
	return s1
	# TODO: merge unused points here? this will eventually need to be done
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
				e = 3*current_face+k
				if regular_edge(s, e)
					next_face = fld1(opposite(s, e), 3)
					if iszero(label[next_face])
						label[next_face] = n
						push!(todo, next_face)
					end
				else # singular edge
					for e1 in radial_loop(s, e)
						l = label[fld1(e1,3)]
						!iszero(l) && (adjacency[l,n] = adjacency[n,l] = e)
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
function sort_radial_loop(s::HalfEdgeMesh, e, pt3 = nothing)
# used for locate_point:
# vec3, proj, dir3, dir2scaled, order
	# prepare geometry information
	(v1, v2) = halfedge(s, e)
	dir3 = point(s, v2) - point(s, v1)
	proj = main_axis(dir3)
	dir2 = project2d(proj, dir3)
	dir2scaled = dir2 ./dot(dir3, dir3)
	# collect half-edges and corresponding opposed vertices
	# for each adjacent face, compute a (3d) vector which, together with
	# the edge, generates the face (and pointing from the edge to the face):
	# 2d projection of face_vec3 (preserving orientation)
	he = collect(radial_loop(s, e)) # half-edges
	ov = [ opposite_vertex(s, e) for e in he ] #opposite vertices
	# we could use this to determine edge orientation:
# 	dv = [ destination(s, e) == v2 for e in he ]
	p1 = point(s, v1)
	fv = [ point(s, v) - p1 for v in ov ] # face vector
	# face vector, projected in 2d:
	fv2= [ project2d(proj, v) .- dot(v, dir3) .* dir2scaled for v in fv ]
# 	println("edge $e = ($v1,$v2): $dir3 proj=$proj")
# 	for (e1, x, y, z) in zip(he, ov, fv, fv2)
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
		return (-1)^i1*he[i1] < (-1)^i2*he[i2]
		# half-edge number is proportional to face number
	end
	reorder = sort(1:length(he), lt=face_cmp)
	pt3 == nothing && return he[reorder]

	# find where `pt3 - v1` inserts in radial loop:
	vec3 = pt3 - p1
	vec2 = project2d(proj, vec3) .- dot(vec3, dir3) .*dir2scaled
	@assert !all(iszero,vec2) "half-edge $e aligned with point $vec3"
	k = searchsorted(fv2[reorder], vec2,
		lt = (u,v)->circular_sign(u,v) > 0)
	@assert k.start > k.stop "impossible to determine location at this edge"
	# possibilities are: (i+1:i) for i in 0..n
	k.stop == 0 && return he[reorder[end]]
	return he[reorder[k.stop]]
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
	((e, j), state) = iterate(nv)
	vpj = point(s, j) - p
	xj = dot(vpi, vpj)
	iszero(xj) && return e
	xyj = (xj*xj, norm²(cross(vpi, vpj)))
	best = e
	while true
		u = iterate(nv, state)
		u == nothing && return best
		((e, k), state) = u
		vpk = point(s, k) - p
		xk = dot(vpi, vpk)
		iszero(xk) && return e
		xyk = (xk*xk, norm²(cross(vpi, vpk)))
		if xyk[2]*xyj[1] > xyk[1]*xyj[2]
			best = e; xyj = xyk
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
	e = find_good_halfedge(s, closest, p)
	f = sort_radial_loop(s, e, p)
	# f is a half-edge in the radial loop of e
	return (fld1(f, 3), destination(s, f) ≠ destination(s, e))
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
		elist = sort_radial_loop(s, eindex)
		n = length(elist)
		# plist is the sorted list of regular patches at this edge
		# dlist is the list of edge orientations
		plist = rp.label[fld1.(abs.(elist), 3)]
		dlist = [destination(s, e) for e in elist]
		v2 = destination(s, eindex)
		e1 = elist[1]; p1 = plist[1]; d1 = dlist[1]
# 		println("sorted radial loop is $elist")
		for i in 2:n
			e2 = elist[i]; p2 = plist[i]; d2 = dlist[i]
			# patch p2 is positively oriented iff d2==v2, etc.:
			# if p1 is positively oriented (d1==v2) then cell between p1 and p2
			# is 2p1 (otherwise 2p1-1);  if p2 is positively oriented (d2==v2),
			# this same cell is 2p2-1 (else 2p2)
# 			union!(cells, 2*p1-(d1≠v2), 2*p2-(d2==v2))
			k = 1-(d1==v2)-(d2==v2)
			connect!(levels, p1, p2, k)
			e1 = e2; p1 = p2; d1 = d2
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
	# face_cc[i] is the connected component of face i
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
	cc_nest_level = SortedAssoc(cc_list, cc_nest)

	face_level = Vector{Int}(undef, nfaces(s))
	for f in 1:nfaces(s)
		face_level[f] = 1 + levels.level[rp.label[f]] + cc_nest_level[face_cc[f]]
	end
	for (f1, f2) in opposite_faces(s)
		face_level[f1]-= 1
		face_level[f2]-= 1
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
function combine(meshes, μ, ε = 0)
	newmesh = subtriangulate(concatenate(meshes...), ε)
	levels = multiplicity(newmesh)
	println("newmesh: faces=$(collect(pairs(faces(newmesh)))) levels=$(collect(pairs(levels)))")
	select_faces!(newmesh, levels .== μ)
end
@inline Base.union(meshes::HalfEdgeMesh...; ε = 0) =
	combine(meshes, 1, ε)
@inline Base.intersect(meshes::HalfEdgeMesh...; ε = 0) =
	combine(meshes, length(meshes), ε)
@inline Base.setdiff(m1::HalfEdgeMesh, m2::HalfEdgeMesh; ε = 0) =
	combine([m1, reverse(m2)], 2, ε)
#»»1
function validate(s::HalfEdgeMesh)
	for (i, j) in pairs(s.opposite)
		j ∈ keys(s.opposite) || println("edge e$i: opposite = $j, invalid")
		e = s.opposite[j]
		halfedge(s, j) == reverse(halfedge(s, i)) ||
		println("edge e$i =$(halfedge(s,i)): opposite e$j = $(halfedge(s,j))")
# 		e == i || println("edge $i: opposite² = $e")
	end
	for (i,e) in pairs(s.edgefrom)
		j = s.destination[prev(s, e)]
		j == i || println("vertex $i: edgefrom = $v comes from $j")
	end
end
function explain(s::HalfEdgeMesh, e)
	opp=opposite(s, e)
	j=fld1(e, 3)
	print("""
half-edge $e: ->$(destination(s,e)), opp=$opp -> $(destination(s,opp))
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
