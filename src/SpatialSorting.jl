"""
    SpatialSorting

This module provides a way to compute intersections among a set of boxes.
The two following functions are provided (not exported):

 * `tree = SpatialSorting.tree(boxes)`: returns a bounded volume hierarchy
   computed from the given bounding boxes.
 * `search(tree, box)`: returns the vector of indices of those `boxes`
   which intersect the `box` parameter.
 * `intersections(boxes)`: returns the set of all intersecting pairs of boxes in   the list.

The `tree` function has an average-time quasi-linear complexity
(w.r.to the number of boxes), while `search(tree, box)` is quasi-constant.

The `boxes` may be of any type, as long as the following operations exist:
 - `SpatialSorting.position(box)`: converts the box to a point
   (`AbstractVector`) representing the position of the box
   (e.g. the box center, or any affine transformation thereof);
 - `Base.merge(box1, box2)`: (any superset of) the union of the two boxes;
 - `SpatialSorting.intersects(box1, box2)`: returns true iff intersection is not empty (default is `!isempty(box1 ∩ box2)`).
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

"""
    CompleteBinTree{B}

A tree describing a bounded volume hierarchy, using boxes of type `B`.

This hierarchy is represented using a complete binary tree:
inner node `k` has children `2k` and `2k+1`.
Left-side box of node `i` is `box[2i-1]` and right-side is `box[2i]`.

In addition, a `leaf` table indicates the permutation of the leaves of the tree.
"""
struct CompleteBinTree{B}
# replace lbox[i] by box[2i-1] and rbox[i] by box[2i]
	box::Vector{B}
	leaf::Vector{Int} # permutation of leaves
	firstleaf::Int # precomputation
	@inline CompleteBinTree{B}(::UndefInitializer, n::Int) where{B} =
		new{B}(Vector{B}(undef, 2n-2), 
# 		Vector{B}(undef, n-1),
		collect(1:n),
		prevpow(2, n-1)<<1)
end

@inline nleaves(t::CompleteBinTree) = length(t.leaf)

function to_leaf(t::CompleteBinTree, i::Integer)
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
	tree = CompleteBinTree{eltype(boxes)}(undef, n)
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
		merge(tree.box[2x-1], tree.box[2x])
	@inbounds for i in 2n-2:-1:1
		tree.box[i] = box(i+1)
	end
	return tree
end

"""
    SpatialSorting.search(tree, box)

Returns the list of indices of given `tree` which
non-trivially intersect the `box`.
"""
function search(t::CompleteBinTree, box)
	todo = [1]
	leaves = Int[]
	n = nleaves(t)
	while !isempty(todo)
		i = pop!(todo)
		for j in (2i, 2i+1) if intersects(t.box[j-1], box)
			if j ≥ n push!(leaves, to_leaf(t, j)); else push!(todo, j); end
		end end
	end
	return leaves
end

function print_tree(io::IO, t::CompleteBinTree, i=1)
	print(io, " "^ndigits(i, base=2))
	if i < length(t.leaf)
		println(io, "node(",i,"):", t.box[2i-1], t.box[2i])
		print_tree(io, t, 2i)
		print_tree(io, t, 2i+1)
	else
		println("$i -> leaf $(to_leaf(t, i))")
	end
end

"""
    SpatialSorting.intersections(boxes)

Returns the set of all non-empty intersections of given `boxes`,
as a `Vector` of pairs of indices.
Complexity is quasi-linear in the number of boxes
(but also depends on the number of intersecting pairs).
"""
function intersections(boxes)
	t = tree(boxes)
	n = length(boxes)
	r = NTuple{2,Int}[]
	todo = [(1,1)]
	while !isempty(todo)
		(a,b) = pop!(todo)
		# (a, b) are either branches (<n) or leaves (≥n)
		if a < n
			if b < n
			# two branches
				intersects(t.box[2a-1], t.box[2b-1]) && push!(todo, (2a, 2b))
				intersects(t.box[2a], t.box[2b]) && push!(todo, (2a+1, 2b+1))
				intersects(t.box[2a-1], t.box[2b]) && push!(todo, (2a, 2b+1))
				# special case: remove duplicate pairs (i,j) and (j,i)
				a≠b && intersects(t.box[2a], t.box[2b-1]) && push!(todo, (2a+1, 2b))
			else # b is a node
				intersects(t.box[2a-1], t.box[b-1]) && push!(todo, (2a, b))
				intersects(t.box[2a], t.box[b-1]) && push!(todo, (2a+1, b))
			end
		else
			if b < n
				intersects(t.box[a-1], t.box[2b-1]) && push!(todo, (a, 2b))
				intersects(t.box[a-1], t.box[2b]) && push!(todo, (a, 2b+1))
			else
				# we computed in the previous step that these two boxes intersect:
				a ≠ b && push!(r, (to_leaf(t, a), to_leaf(t, b)))
			end
		end
	end
	return r
end

Base.show(io::IO, t::CompleteBinTree) = print_tree(io, t)
print_tree(t::CompleteBinTree) = print_tree(stdout, t)

end # module
