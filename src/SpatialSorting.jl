"""
    SpatialSorting

Contains two useful functions:

 * `tree = tree(boxes)`: returns a bounded volume hierarchy made from the given bounding boxes.
 * `search(tree, box)`: returns the vector of indices of those `boxes` which intersect the `box` parameter.

The `boxes` may be of any type, as long as the following operations exist:
 - `SpatialSorting.position(box)`: converts the box to point (`AbstractVector`)
 representing the position of the box (e.g. the box center);
 - `Base.merge(box1, box2)`: (any superset of) the union of the two boxes;
 - `SpatialSorting.intersects(box1, box2)`: returns true iff intersection is not empty
 (default is `!isempty(box1 ∩ box2)`).
"""
module SpatialSorting
using StructArrays
using StaticArrays

@inline intersects(box1, box2) = !isempty(box1 ∩ box2)
@inline position(box) = error("`position` not implemented for $(typeof(box))")

@inline spread(v) = (e = extrema(v); e[2]-e[1])

function largest_spread_coord(points)
	N = length(first(points))
	@inbounds s = [spread(p[i] for p in points) for i in 1:N]
	return findmax(s)[2]
end

# Complete-tree type««1

"""
    CompleteAABBTree{B}

A tree describing a bounded volume hierarchy, using boxes of type `B`.

This hierarchy is implemented using a complete binary tree:
inner node `k` has children `2k` and `2k+1`.
In addition, a `leaf` table indicates the rearrangement of boxes
along the tree.
"""
struct CompleteAABBTree{B}
# this is a complete binary tree, i.e. inner node k has children 2k and
# 2k+1. The `leaf` vector indicates the renumbering of leaf nodes.
# lbox[i] and rbox[i] are the bounding boxes for the two subtrees.
	lbox::Vector{B}
	rbox::Vector{B}
	leaf::Vector{Int}
	firstleaf::Int # precomputation
	@inline CompleteAABBTree{B}(::UndefInitializer, n::Int) where{B} =
		new{B}(Vector{B}(undef, n-1), Vector{B}(undef, n-1), collect(1:n),
		prevpow(2, n-1)<<1)
end

@inline nleaves(t::CompleteAABBTree) = length(t.leaf)

function to_leaf(t::CompleteAABBTree, i::Integer)
	d = t.firstleaf
	i ≥ d && return t.leaf[i-d+1]
	return t.leaf[i-d+1+nleaves(t)]
end

"""
    SpatialSorting.tree(boxes; position)

Returns a tree representing the spatial disposition of the boxes.

 - `position`: a function mapping a box to a point. This should respect
 the spatial disposition of the boxes.
"""
function tree(boxes; position=position)
	n = length(boxes)
	@inline leftpart(r::Integer) =
		let m = 1 << ((sizeof(r)<<3)-leading_zeros(r-1)-1)
			min(m, r-m>>1)
		end
	@inline leftpart(r::UnitRange) =
		first(r):first(r) + leftpart(length(r))-1
	@inline rightpart(r::UnitRange) =
		first(r)+leftpart(length(r)) : last(r)

	firstleaf = prevpow(2, n-1)
	@inline leftleaf(k::Integer) =
		(k ≥ firstleaf) ? 2(k-firstleaf)+1 : n-2(firstleaf-k)+1

	# initialize the tree structure
	tree = CompleteAABBTree{eltype(boxes)}(undef, n)
	# `interval[k]` is the range of leaves accessed by inner node `k`.
	interval = Vector{UnitRange{Int}}(undef, n-1)
	interval[1] = 1:n
	@inbounds for i in 2:n-1
		interval[i] =
			iseven(i) ? leftpart(interval[i>>1]) : rightpart(interval[i>>1])
	end
	points = position.(boxes)

	# spatial sort of indices
	for i in 1:n-1
		subtree = view(tree.leaf, @inbounds interval[i])
		j = largest_spread_coord(view(points, subtree))
		k = leftpart(length(subtree))
		partialsort!(subtree, k; by=i->points[i][j])
	end
	@inline box(x) = x ≥ n ? boxes[to_leaf(tree, x)] :
			merge(tree.lbox[x], tree.rbox[x])
	@inbounds for i in n-1:-1:1
		tree.lbox[i] = box(2i)
		tree.rbox[i] = box(2i+1)
	end
	return tree
end

function search(t::CompleteAABBTree, box)
	todo = [1]
	leaves = Int[]
	n = nleaves(t)
	while !isempty(todo)
		i = pop!(todo)
		if intersects(t.lbox[i], box)
			if 2i ≥ n
				push!(leaves, to_leaf(t, 2i))
			else
				push!(todo, 2i)
			end
		end
		if intersects(t.rbox[i], box)
			if 2i+1 ≥ n
				push!(leaves, to_leaf(t, 2i+1))
			else
				push!(todo, 2i+1)
			end
		end
	end
	return leaves
end

function print_tree(io::IO, t::CompleteAABBTree, i=1)
	print(io, " "^ndigits(i, base=2))
	if i ≤ length(t.lbox)
		println(io, "node(",i,"):", t.lbox[i], t.rbox[i])
		print_tree(io, t, 2i)
		print_tree(io, t, 2i+1)
	else
		println("$i -> leaf $(to_leaf(t, i))")
	end
end

Base.show(io::IO, t::CompleteAABBTree) = print_tree(io, t)
print_tree(t::CompleteAABBTree) = print_tree(stdout, t)

#»»1
end # module
