using AbstractTrees
import AbstractTrees: children
using StaticArrays
struct BoundingBox{N,T}
	min::SVector{N,T}
	max::SVector{N,T}
end

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

@inline printbox(io::IO, b::BoundingBox) =
	print(io, b.min, "⋯", b.max)

# AABBTree ««1
# Type and properties««2
abstract type AbstractAABBTree{N,T} end
# 2 possibilities:
# branch: box, branch= { left, right } where left, right = branch|leaf
# leaf: { data }
struct AABBTree{N,T,D} <: AbstractAABBTree{N,T}
	box::BoundingBox{N,T}
	# contains either 2 child nodes, or leaf data:
	data::Union{NTuple{2,AABBTree{N,T}}, D}

	# leaf
	AABBTree{N,T,D}(box::BoundingBox{N,T}, data::D) where{N,T,D} =
		new{N,T,D}(box, data)
	# non-leaf
	AABBTree{N,T}(box::BoundingBox{N,T},
		left::AABBTree{N,T,D}, right::AABBTree{N,T,D}) where{N,T,D} =
		new{N,T,D}(box, (left, right))
end

@inline ndims(::Type{<:AbstractAABBTree{N}}) where{N} = N
@inline ndims(tree::AbstractAABBTree) = ndims(typeof(tree))
@inline coordtype(::Type{<:AbstractAABBTree{N,T}}) where{N,T} = T
@inline coordtype(tree::AbstractAABBTree) = coordtype(typeof(tree))
@inline datatype(::Type{AABBTree{N,T,D}}) where{N,T,D} = D
@inline datatype(tree::AbstractAABBTree) = datatype(typeof(tree))

# type inference on parameters:
@inline AABBTree(box::BoundingBox, args...) =
	AABBTree{ndims(box)}(box, args...)
@inline AABBTree{N}(box::BoundingBox{N}, args...) where{N} =
	AABBTree{N,coordtype(box)}(box, args...)
@inline AABBTree{N,T}(box::BoundingBox{N,T}, data) where{N,T} =
	AABBTree{N,T,typeof(data)}(box, data)

# special case: trivial box
@inline (X::Type{<:AABBTree})(point::StaticVector, args...) =
	X(BoundingBox(point, point), args...)

# special case: mismatched dimension
@inline AABBTree{N}(box::BoundingBox{N1}, data...) where{N,N1} = throw(
DimensionMismatch("tree of dimension $N initialized with box of dimension $N1"))
@inline (T::Type{<:AABBTree})(left::AABBTree{N}, right::AABBTree{N}) where{N}=
	AABBTree(box(left) ∨ box(right), left, right)


@inline isleaf(t::AABBTree) = !isa(t.data, NTuple{2,<:AABBTree})
@inline data(t::AABBTree) = isleaf(t) ? t.data : nothing
@inline children(t::AABBTree) = isleaf(t) ? () : t.data
box(t::AABBTree) = t.box
printnode(io::IO, t::AABBTree) =
	(printbox(io, box(t)); isleaf(t) && print(io, ": ", t.data))

l1=AABBTree(SA[1,1], "leaf 1")
l2=AABBTree(BoundingBox([0,2],[5,3]), "leaf 2")
root=AABBTree(l1, l2)
