# using AbstractTrees
# include("printing.jl")
using StaticArrays
using StructArrays

# import AbstractTrees: children

# Utilities««1
@inline median(v::AbstractVector; kwargs...) = partialsort(v, (length(v)+1)>>1)

# Bounding box««1
struct BoundingBox{N,T}
	min::SVector{N,T}
	max::SVector{N,T}
end

one_half(x)=x/2

BoundingBox{N}(min::AbstractVector, max::AbstractVector) where{N} =
	BoundingBox{N,promote_type(eltype(min), eltype(max))}(min, max)
BoundingBox(min::StaticVector{N}, max::StaticVector{N}) where{N} =
	BoundingBox{N}(min, max)
function BoundingBox(min::AbstractVector, max::AbstractVector)
	length(min) == length(max) || throw(DimensionMismatch(
	"min and max must have same length ($(length(min)),$(length(max)))"))
	return BoundingBox{length(min)}(min, max)
end

@inline ndims(::Type{BoundingBox{N,T}}) where{N,T} = N
@inline ndims(b::BoundingBox) = ndims(typeof(b))
@inline coordtype(::Type{BoundingBox{N,T}}) where{N,T} = T
@inline coordtype(b::BoundingBox) = coordtype(typeof(b))

@inline Base.in(point::SVector{N}, box::BoundingBox{N}) where{N} =
	all(point .≥ box.min) && all(point .≤ box.max)
@inline ⊂(box1::BoundingBox{N}, box2::BoundingBox{N}) where{N} =
	all(box1.min .≥ box2.min) && all(box1.max .≤ box2.max)
@inline ∨(box1::BoundingBox{N}, box2::BoundingBox{N}) where{N} =
	BoundingBox{N}(min.(box1.min, box2.min), max.(box1.max, box2.max))
@inline Base.:∩(box1::BoundingBox{N}, box2::BoundingBox{N}) where{N} =
	BoundingBox{N}(max.(box1.min, box2.min), min.(box1.max, box2.max))
@inline Base.isempty(box::BoundingBox) = any(box.min .> box.max)
@inline volume(box::BoundingBox) = prod(box.max .- box.min)
# perimeter = (N-1)-volume of boundary
@inline perimeter(box::BoundingBox) = perimeter(box.max .- box.min)
@inline _perimeter(vec::SVector{N}) where{N} =
	sum(prod(vec[j] for j in 1:N if j ≠ i) for i in 1:N)
@inline twicecenter(box::BoundingBox) = box.min + box.max

@inline printbox(io::IO, b::BoundingBox) =
	print(io, b.min, "⋯", b.max)

# AABBTree ««1
# Abstract type ««2
abstract type AbstractAABBTree{N,T} end
@inline ndims(::Type{<:AbstractAABBTree{N}}) where{N} = N
@inline ndims(tree::AbstractAABBTree) = ndims(typeof(tree))
@inline coordtype(::Type{<:AbstractAABBTree{N,T}}) where{N,T} = T
@inline coordtype(tree::AbstractAABBTree) = coordtype(typeof(tree))

# Indexed type ««2

struct SizedUndef
	n::Int
end
@inline Base.convert(T::Type{<:AbstractVector}, s::SizedUndef) = T(undef, s.n)
@inline nothings(n) = Vector{Nothing}(undef, n)

struct IndexedAABBTree{N,T} <: AbstractAABBTree{N,T}
	# inner nodes have {left, right} child pointers + a bounding box
	# negative index indicates leaf, positive index indicates inner node
	# leaf nodes have data::D + a bounding box
	inner::StructVector{@NamedTuple{left::Int, right::Int, box::BoundingBox{N,T}}}
	leaf::Vector{BoundingBox{N,T}}


	@inline (X::Type{IndexedAABBTree{N,T}})(::UndefInitializer,
			n::Int) where{N,T} =
		# this is a binary tree; it has (n-1) inner nodes for (n) leaves
		new{N,T}(SizedUndef(n-1), SizedUndef(n))
	@inline (X::Type{<:IndexedAABBTree})(::UndefInitializer,
		boxes::AbstractVector{BoundingBox{N,T}}) where {N,T} =
		new{N,T}(SizedUndef(length(boxes)-1), boxes)
end

# AbstractTrees.treekind(::IndexedAABBTree) = AbstractTrees.IndexedTree()
@inline Base.getindex(tree::IndexedAABBTree, node::Int) =
	node > 0 ? tree.inner[node] : tree.leaf[-node]
# @inline AbstractTrees.rootstate(::IndexedAABBTree) = 1
@inline childindices(tree::IndexedAABBTree, node::Int) =
	_childindices(tree[node])
@inline _childindices(t::NamedTuple)::NTuple{2,Int} = (t.left, t.right)
@inline _childindices(t::BoundingBox) = ()

@inline spread(v) = ( e = extrema(v); e[2]-e[1])

function largest_spread_coord(points)
	N = length(first(points))
	s = collect(spread(p[i] for p in points) for i in 1:N)
	return findmax(s)[2]
end

function build_tree(boxes::AbstractVector{<:BoundingBox})
	tree = IndexedAABBTree(undef, boxes)
	length(boxes) ≤ 1 && return tree
	for i in 1:length(boxes)-1
		LazyRow(tree.inner, i).box =
			BoundingBox(zero(boxes[1].min), zero(boxes[1].max))
	end
	points = twicecenter.(boxes)
	# leaves are managed by swapping their indices around in this table:
	redirect= collect(1:length(boxes))
	# list of tree parts to sort
	todo = [(1, 1, length(boxes))]
	while !isempty(todo)
		(a, b, n) = pop!(todo)
		if n == 2
			LazyRow(tree.inner, b).left = -redirect[a]
			LazyRow(tree.inner, b).right = -redirect[a+1]
			continue
		elseif n == 1
			continue # terminal case
		end
		idx = view(redirect, a:a+n-1)
		j = largest_spread_coord(view(points, idx))
		k = (n+1)>>1 # fld(n, 2): number of guys to put on the left side
		# the leaves on left side will be renumbered (a:a+k-1),
		# those on the right side (a+k:a+n-1);
		# inner root is b, left side is (b+1:b+k-1), right (b+k:b+n-2)
		u = partialsort!(idx, k; by=i->points[i][j])

		# left-hand side always has at least 2 nodes in it
		LazyRow(tree.inner, b).left = (b+1)
		push!(todo, (a, b+1, k))
		# right-hand side might be a singleton, thus a leaf:
		if n == 3 # in this case the split is (a, a+1 || a+2):
			LazyRow(tree.inner, b).right = -redirect[a+2]
		else
			LazyRow(tree.inner, b).right = (b+k)
			push!(todo, (a+k, b+k, n-k))
		end
	end
	# this tree is monotonic, i.e. the children of an inner node always
	# have larger labels than the inner node
	# This means that we can go right-to-left in computing the bounding
	# boxes:
	for b in length(tree.inner):-1:1
		l = LazyRow(tree.inner, b)
		box1 = (l.left < 0) ? tree.leaf[-l.left] : tree.inner[l.left].box
		box2 = (l.right < 0) ? tree.leaf[-l.right] : tree.inner[l.right].box
		l.box = box1 ∨ box2
	end

	return tree
end

function leaves_inter(tree::IndexedAABBTree, box::BoundingBox)
	todo = [1]
	leaves = Int[]
	while !isempty(todo)
		b = pop!(todo)
		s = tree.inner[b]
		x = s.left
		if x < 0
			!isempty(tree.leaf[-x] ∩ box) && push!(leaves, -x)
		else
			!isempty(tree.inner[x].box ∩ box) && push!(todo, x)
		end
		x = s.right
		if x < 0
			!isempty(tree.leaf[-x] ∩ box) && push!(leaves, -x)
		else
			!isempty(tree.inner[x].box ∩ box) && push!(todo, x)
		end
	end
	return leaves
end



#»»1
Base.show(b::BoundingBox) = print(io, "BoundingBox(", b.min, ",", b.max, ")")
b = BoundingBox([7,0,0],[10,10,10])
using Random; Random.seed!(1789)
boxes = [ BoundingBox(SVector(rand(1:6,3)...),SVector(rand(4:10,3)...))
	for _ in 1:10 ]
t=build_tree(boxes);
function _print_tree(io::IO, t::IndexedAABBTree, i=1, depth="")
	print(io, depth)
	if i > 0
		print(io, "inner($i)")
		printbox(io, t.inner[i].box)
		println(io)
		_print_tree(io, t, t.inner[i].left, depth*"  ")
		_print_tree(io, t, t.inner[i].right, depth*"  ")
	else
		print(io, "leaf($(-i))")
		printbox(io, t.leaf[-i])
		println(io)
	end
end

print_tree(io::IO, t::IndexedAABBTree) = _print_tree(io, t, 1, "")
print_tree(t::IndexedAABBTree) = print_tree(stdout, t)

# l1=AABBTree(SA[1,1], "leaf 1")
# l2=AABBTree(BoundingBox([0,2],[5,3]), "leaf 2")
# root=AABBTree(l1, l2)
