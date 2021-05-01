"""
    Projectors

This module contains tool to define 3d → 2d projections:
 - `main_axis`, `project2d`, `project1d`
 - `Projector` class (encapsulating projection and lift).
"""
module Projectors
using LinearAlgebra
using StaticArrays
# main_axis««1
"""
    main_axis(direction)

Returns that of the six basic axes (±x,y,z, labelled as integers ±1,2,3)
which is closest to the given `direction`.
"""
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
"""
    project2d(axis, vector)

Given an axis ±1,2,3 computed by `main_axis`, computes
a faithful, orientation-preserving projection of the `vector`
as a normal to the axis, given by selecting two coordinates.
"""
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
"""
    project1d(axis, vector)

Given an axis ±1,2,3 computed by `main_axis`,
computes the projection of the `vector` along the axis.
"""
function project1d(axis, vec)#««
	axis == 1 && return vec[1]
	axis == 2 && return vec[2]
	axis == 3 && return vec[3]
	axis ==-1 && return -vec[1]
	axis ==-2 && return -vec[2]
	@assert axis == -3
	             return -vec[3]
end#»»

# Projector ««1
"""
    Projector

A structure containing information about a 3d -> 2d projection
defined via selecting coordinates (so not orthogonal),
and optionally about the lift back to a certain plaine in 3d.

Can be applied to any type by defining methods for `lift` and `project`.
"""
struct Projector{T}
	dir::Int8
	lift::NTuple{3,T}
end
struct FormalInv{X}
	x::X
end
"""
    Projector(direction)

Computes a linear projector,
i.e. just a pair of coordinates to extract for
orientation-preserving, faithful 3d→2d projection.
"""
@inline function Projector(direction::AbstractVector) # no lift
	@inbounds u1 = direction[1]
	@inbounds u2 = direction[2]
	@inbounds u3 = direction[3]
  @inbounds a1 = abs(u1)
	@inbounds a2 = abs(u2)
	@inbounds a3 = abs(u3)
	P = Projector{eltype(direction)}
	if a1 < a2
		if a2 < a3 @goto max3; end
		u2 > 0 && return P( 2, (0,0,0))
		          return P(-2, (0,0,0))
	elseif a1 < a3
		@label max3
		# (x,y) lifts to (y,x,z) st ay+bx+cz=d i.e. z=(d-bx-ay)/c
		u3 > 0 && return P( 3, (0,0,0))
		          return P(-3, (0,0,0))
	else # max1
		# (x,y) lifts to (z,x,y) st az+bx+cy=d i.e. z=(d-bx-cy)/a
		u1 > 0 && return P( 1, (0,0,0))
		          return P(-1, (0,0,0))
	end
end

"""
    Projector(direction, point)

Computes an affine projector, i.e. defines the affine lift
as well as the projection.
"""
@inline function Projector(direction::AbstractVector, point)
	@inbounds u1 = direction[1]
	@inbounds u2 = direction[2]
	@inbounds u3 = direction[3]
  @inbounds a1 = abs(u1)
	@inbounds a2 = abs(u2)
	@inbounds a3 = abs(u3)
	d = dot(direction, point)
	P = Projector{eltype(direction)}
	if a1 < a2
		if a2 < a3 @goto max3; end
		# (x,y) lifts to (y,z,x) st ay+bz+cx=d, ie z=(d-cx-ay)/b
		u2 > 0 && return P( 2, (-u3/u2, -u1/u2, d/u2))
		# (x,y) lifts to (x,z,y) st ax+bz+cy=d, ie z=(d-ax-cy)/b
		          return P(-2, (-u1/u2, -u3/u2, d/u2))
	elseif a1 < a3
		@label max3
		# (x,y) lifts to (y,x,z) st ay+bx+cz=d i.e. z=(d-bx-ay)/c
		u3 > 0 && return P( 3, (-u1/u3, -u2/u3, d/u3))
		          return P(-3, (-u2/u3, -u1/u3, d/u3))
	else # max1
		# (x,y) lifts to (z,x,y) st az+bx+cy=d i.e. z=(d-bx-cy)/a
		u1 > 0 && return P( 1, (-u2/u1, -u3/u1, d/u1))
		          return P(-1, (-u3/u1, -u2/u1, d/u1))
	end
end
"""
    project(p::Projector, u)

Computes the projection for a given vector.
This has, by default, methods for `Tuple`, `AbstractVector`,
and `StaticVector{3}`. To make this module work with other types,
these methods should be extended (cf. the `StaticVector{3}` method definition).
"""
@inline function project(p::Projector, u)
	p.dir == 1 && return (u[2], u[3])
	p.dir == 2 && return (u[3], u[1])
	p.dir == 3 && return (u[1], u[2])
	p.dir ==-1 && return (u[3], u[2])
	p.dir ==-2 && return (u[1], u[3])
	@assert p.dir ==-3
	              return (u[2], u[1])
	@assert false
end
"""
    lift(p::Projector, u)

Computes the affine lift for a vector.
This has, by default, methods for `Tuple`, `AbstractVector`,
and `StaticVector{2}`. To make this module work with other types,
these methods should be extended (cf. the `StaticVector{2}` method definition).
"""
@inline function lift(p::Projector, u)
	z = p.lift[1]*u[1]+p.lift[2]*u[2]+p.lift[3]
	p.dir == 1 && return (z, u[1], u[2])
	p.dir == 2 && return (u[2], z, u[1])
	p.dir == 3 && return (u[1], u[2], z)
	p.dir ==-1 && return (z, u[2], u[1])
	p.dir ==-2 && return (u[1], z, u[2])
	# for type-stability, no `if` here:
	@assert p.dir ==-3 
	              return (u[2], u[1], z)
end
# syntactic sugar: p(u) and inv(p)(u)
@inline (p::Projector)(u) = project(p,u)
Base.inv(p::Projector) = FormalInv{typeof(p)}(p)
@inline (l::FormalInv{<:Projector})(u) = lift(l.x, u)

@inline lift(p::Projector, u::AbstractVector) = [lift(p,(u...,))...]
@inline project(p::Projector, u::AbstractVector) = [project(p,(u...,))...]
@inline lift(p::Projector, u::StaticVector{2}) =
	SVector{3,eltype(u)}(lift(p,u.data))
@inline project(p::Projector, u::StaticVector{3}) =
	SVector{2,eltype(u)}(project(p,(u.data)))

#»»1
export main_axis, project2d, project1d
export Projector, project, lift

end # module

