module Bezier


using StaticArrays
using LinearAlgebra


struct BezierCurve{N,D,T}
	p::SVector{N,SVector{D,T}}
	@inline BezierCurve(p::Vararg{Any,N}) where{N} = BezierCurve{N}(p...)
	@inline BezierCurve{N}(p::StaticVector{D}...) where{N,D} =
		new{N,D,promote_type(eltype.(p)...)}(p)
	@inline BezierCurve{N}(p::AbstractVector...) where{N} =
		BezierCurve{N}(SVector{length(first(p))}.(p)...)
end

@inline degree(::Type{<:BezierCurve{N}}) where{N} = N
@inline dim(::Type{<:BezierCurve{N,D}}) where{N,D} = D
@inline pttype(::Type{<:BezierCurve{N,D,T}}) where{N,D,T} = SVector{D,T}
@inline degree(b::B) where{B<:BezierCurve} = degree(B)
@inline dim(b::B) where{B<:BezierCurve} = dim(B)
@inline pttype(b::B) where{B<:BezierCurve} = pttype(B)
@inline controlpoint(b::BezierCurve, i::Integer) = b.p[i]

"""
    interpolate(b::BezierCurve, t::Real)

Computes the position of the point at parameter `t` on the curve `b`.
"""
@inline interpolate(b::BezierCurve, t::Real) = _interpolate(b.p, t)

# recursive definition for de Casteljau's algorithm:
@inline _interpolate(p::SVector{1}, _) = p[1]
@inline _interpolate(p::SVector{N,P}, t) where{N,P} =
	_interpolate(SVector{N-1,P}(t*p[i+1]+(1-t)*p[i] for i in 1:N-1), t)

"""
    interpolate(b::BezierCurve, v::AbstractVector{<:Real})

Interpolates `b` at all values t ∈ `v`.
"""
@inline interpolate(b::BezierCurve, v::AbstractVector{<:Real}) =
	(interpolate(b, t) for t in v)

"""
    interpolate(b::BezierCurve, n::Integer)

Interpolates as `n` regularly-spaced points (in parameter space).
"""
@inline interpolate(b::BezierCurve, n::Integer) = interpolate(b, (0:n)/n)

end
