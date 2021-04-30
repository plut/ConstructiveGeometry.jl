module EquivalenceStructures
using DataStructures
# uniquenames ««1
"""
    uniquenames(list)

Returns `(names, idx)`, where `names = unique(sort(list))`,
and `list == names[idx]`.
"""
function uniquenames(list)
	# used by: equivalence_structure
	isempty(list) && return (similar(list,0), Int[])
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

# union-find structure««1
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

# find representative given equivalence relation ««1
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

# Equivalence structure for arbitrary (sortable) objects««1
"""
    EquivalenceStructure

Equivalence structure for arbitrary (sortable) objects

Useful methods:
 - `elements(s)`: returns (unique+sorted) list of all elements
 - `equivalent(s, x,y)`: true iff x, y are equivalent
 - `class(s, x)`: returns all elements equivalent to `x`
 - `classes(s)`: returns an iterator of all classes of `s`
 - `representative(s,x)` is same as `first(class(s,x))`
"""
# from an element, find an equivalence class identifier (integer)
# from an equivalence class id, find (sorted) list of equiv. elements
struct EquivalenceStructure{E}
	elements::E # element list: either a sorted Vector, or OneTo (for Int)
	class_idx::Vector{Int}
	class::Vector{Vector{Int}}
end
@inline elements(s::EquivalenceStructure) = s.elements
@inline classes(s::EquivalenceStructure) = [s.elements[c] for c in s.class]
@inline function index(s::EquivalenceStructure, x)
	idx = searchsortedfirst(s.elements, x)
	@boundscheck (idx > length(s.elements) || s.elements[idx] ≠ x) &&
		return nothing
	return idx
end
@inline class_idx(s::EquivalenceStructure, x) = s.class_idx[index(s, x)]
@inline equivalent(s::EquivalenceStructure, x, y) = index(s, x) == index(s, y)
@inline class(s::EquivalenceStructure, x) = s.elements[s.class[class_idx(s, x)]]
@inline representative(s::EquivalenceStructure, x) =
	s.elements[first(s.class[class_idx(s,x)])]
struct Representatives{E}
	eqv::EquivalenceStructure{E}
end
@inline representatives(s::EquivalenceStructure) = Representatives(s)
@inline Base.keys(r::Representatives) = elements(r.eqv)
@inline Base.values(r::Representatives) =
	(r.eqv.elements[first(r.eqv.class[i])] for i in r.eqv.class_idx)
@inline Base.length(r::Representatives) = length(elements(r.eqv))

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
	T = eltype(eltype(relation))
	isempty(relation) && return EquivalenceStructure(T[], Int[], Vector{Int}[])
	# flatten relation:
	rel_flat = reinterpret(T, relation)

	# `names` is the unique list of names
	# `rel_idx` is the renamed equivalence relation (as small integers):
	(names, indexed_relation) = uniquenames(rel_flat)
	eq1 = equivalence_structure_flat(length(names), indexed_relation)
	return EquivalenceStructure(names, eq1.class_idx, eq1.class)
end#»»
#««1
export EquivalenceStructure, equivalence_structure
export elements, class, classes, equivalent, representative, representatives

export LevelStructure, connect!, root!
end # module
