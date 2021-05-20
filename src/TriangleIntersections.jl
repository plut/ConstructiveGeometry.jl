"""
    TriangleIntersections

Routines for computing the intersection of two 2d or 3d triangles.
The triangle vertices may be of any type which supports:
indexation, `cross`, `dot`, and interaction with `Projector` type.
"""
module TriangleIntersections
using LinearAlgebra
using StaticArrays
using FastClosures
include("Projectors.jl")
using .Projectors

"""
    TriangleIntersections.Constants

Constants describing position of a point in a triangle.
All computations in this module are branch-less
(implemented using binary logic).
"""
module Constants#««
# const Position = Int8
# 	interior=0
# 	vertex1=1
# 	vertex2=2
# 	vertex3=3
# 	edge23=4
# 	edge31=5
# 	edge12=6
# @inline isinterior(a::Position) = a == interior
# @inline isvertex(a::Position) = a ≤ 3
# @inline isedge(a::Position) = a ≥ 4
# @inline same_edge(a::Position, b::Position) =
# @inline next(a::Position) = typeof(a)(Integer(a)<<1|Integer(a)>>2)
# @inline prev(a::Position) = typeof(a)(Integer(a)<<2|Integer(a)>>1)

@enum Position::Int8 begin
	interior=0
	vertex1=2|4
	vertex2=4|1
	vertex3=1|2
	edge23=1
	edge31=2
	edge12=4
	invalid=-1
end
@inline Base.show(io::IO, a::Position) =
	print(io, (a == interior) ? "interior" :
	(a == vertex1) ? "v1" :
	(a == vertex2) ? "v2" :
	(a == vertex3) ? "v3" :
	(a == edge23) ? "e23" :
	(a == edge31) ? "e31" :
	(a == edge12) ? "e12" :
	"invalid")
@inline Base.iszero(a::Position) = a == interior
@inline isvertex(a::Position) = count_ones(Integer(a)) == 2
@inline isedge(a::Position) = count_ones(Integer(a)) == 1
@inline Base.:&(a::Position, b::Position) = Position(Integer(a)&Integer(b))
@inline Base.:|(a::Position, b::Position) = Position(Integer(a)|Integer(b))
@inline Base.:^(a::Position, b::Position) = Position(Integer(a)^Integer(b))
# rotating a triangle (1->2->3) implemented without branches
@inline function Base.:<<(a::Position, i::Integer)
	i == 0 && return a
	i == 1 && return (typeof(a))((Integer(a)<<1|Integer(a)>>2)&7)
	i == 2 && return (typeof(a))((Integer(a)<<2|Integer(a)>>1)&7)
	return a<<(i%3)
end
@inline index(a::Position) = Integer(a) ≤ 3 ? Integer(a) : 7-Integer(a)
# branch-less computation of index for vertices and edges:
@inline index(a::Position, ::typeof(isvertex)) = 6-Integer(a)+Integer(a)>>2
@inline index(a::Position, ::typeof(isedge)) = Integer(a)-Integer(a)>>2
@inline same_edge(a::Position, b::Position) = !iszero(a & b)
end#»»
import .Constants: isvertex, isedge, index, same_edge
# Preamble««1
# Generic stuff««2
# predicates for assertions:
@inline barycenter(p1, p2, λ) = p2 + λ*(p1-p2) # λp1+(1-λ)p2
@inline collinear3(a,b,c) = norm(cross(a-b,c-b),1) ≤ 1e-8
@inline collinear2(a,b,c) = abs(det2(a,b,c)) ≤ 1e-8
@inline monotonic(a,b,c) = dot(a-b,c-b) ≤ 0
@inline samepoint(a,b) = norm(a-b,1) ≤ 1e-8

# @inline plus1mod3(i,j=0)=(i==1) ? 2+j : (i==2) ? 3+j : 1+j
# @inline plus2mod3(i,j=0)=(i==1) ? 3+j : (i==2) ? 1+j : 2+j

# struct IntersectionData««2
struct SizedUndef{N} end
@inline sizedundef(N::Integer) = SizedUndef{N}()
@inline Base.convert(::Type{SizedVector{N,P}}, ::SizedUndef{N}) where{N,P} =
	SizedVector{N}(Vector{P}(undef, N))
"""
    struct IntersectionData

Describes an intersection. Each point of the returned polygon
is marked with its type w.r.t the two input simplexes,
using enums defined in `TriangleIntersections.Constants`
"""
struct IntersectionData{N,P} <:
		AbstractVector{Tuple{P,NTuple{2,Constants.Position}}}
# maximum of N points, constant-size
	npoints::Int
	pttype::SizedVector{N,NTuple{2,Constants.Position}}
	points::SizedVector{N,P}
	@inline IntersectionData{N,P}(::UndefInitializer, npoints::Integer) where{N,P} =
		new{N,P}(npoints, sizedundef(N), sizedundef(N))
# 			SizedVector{N,NTuple{2,Constants.Position}}(undef),
# 			MVector{N,P}(undef))
	@inline IntersectionData{N,P}(npoints::Int, pttype::MVector{N,P},
		points::MVector{N,P}) where{N,P} = new{N,P}(npoints, pttype, points)
	@inline function IntersectionData{N,P}(npoints::Int, pttype, points
		) where{N,P}
		it = IntersectionData{N,P}(undef, npoints)
		# instead of broadcasting, we write the loop: this allows us to
		# over-specify the set of points
		for i in 1:npoints
			it.pttype[i] = pttype[i]; it.points[i] = points[i]
		end
		return it
	end
	@inline IntersectionData{N,P}() where{N,P} = IntersectionData{N,P}(undef, 0)
	@inline IntersectionData{N,P}(pt::Pair...) where{N,P} =
		IntersectionData{N,P}(length(pt), last.(pt), first.(pt))
end

@inline allocsize(::Type{<:IntersectionData{N}}) where{N} = N
@inline allocsize(t::IntersectionData) = allocsize(typeof(t))
@inline Base.size(it::IntersectionData) = (it.npoints,)
@inline Base.getindex(it::IntersectionData, i) = (it.points[i], it.pttype[i])
@inline function Base.show(io::IO, it::IntersectionData)
	if isempty(it) print(io, "ø"); return; end
	print(io, length(it), " points:",
		join(["\n  $a => $b" for (a,b) in it]))
end
function swap!(it::IntersectionData)
	@inbounds for i in 1:length(it)
		it.pttype[i] = it.pttype[i][[2,1]]
	end
	return it
end
@inline function rename1!(it::IntersectionData, i::Integer, old_new::Pair...)
	for (x,y) in old_new
		if it.pttype[i][1] == x
			it.pttype[i] = (y, it.pttype[i][2])
			break
		end
	end
end
@inline function rename1!(it::IntersectionData, old_new::Pair...)
	for i in 1:length(it)
		rename1!(it, i, old_new...)
	end
	return it
end
@inline function shift1!(it::IntersectionData, turn::Integer)
	for i in 1:length(it)
		it.pttype[i] = (it.pttype[i][1]<<turn, it.pttype[i][2])
	end
end
@inline function rename!(it::IntersectionData, side::Int, table)
	if side == 1
		for i in 1:length(it)
			it.pttype[i] = (table[1+it.pttype[i][1]], it.pttype[i][2])
		end
	else
		for i in 1:length(it)
			it.pttype[i] = (it.pttype[i][1], table[1+it.pttype[i][2]])
		end
	end
end


macro tree27(vars,args...)#««2
	@assert Meta.isexpr(vars, :tuple)
	@assert length(vars.args) == 3
	vars = esc.(vars.args); ε = esc(:ε)
	expr = Dict{String, Any}()
	for arg in args
		@assert arg.args[1] == :(=>)
		@assert arg.args[2] isa String
		s = arg.args[2]
		@assert length(s) == 3
		e = arg.args[3]
		expr[s] = Meta.isexpr(e, :block) ? esc.(e.args) : [esc(e)]
	end
	neg=Dict('+'=>'-','-'=>'+','0'=>'0')
	function exprfor(s)
		(a,b,c) = s
		(d,e,f) = (neg[i] for i in s)
		ij=esc(:(turn, flip))
		haskey(expr, "$a$b$c") && return (:($ij=(0,false)),expr["$a$b$c"]...)
		haskey(expr, "$b$c$a") && return (:($ij=(1,false)),expr["$b$c$a"]...)
		haskey(expr, "$c$a$b") && return (:($ij=(2,false)),expr["$c$a$b"]...)
		haskey(expr, "$d$e$f") && return (:($ij=(0,true)),expr["$d$e$f"]...)
		haskey(expr, "$e$f$d") && return (:($ij=(1,true)),expr["$e$f$d"]...)
		haskey(expr, "$f$d$e") && return (:($ij=(2,true)),expr["$f$d$e"]...)
		(:($(esc(:position)) = $s),) # :(error("case not implemented: "*$s)),)
	end
	build(i,f,p) = quote
		if     $(vars[i]) > $ε; $(f(p*'+')...)
		elseif $(vars[i]) <-$ε; $(f(p*'-')...)
		else                  ; $(f(p*'0')...)
		end
	end
	expr1(p) = (build(3,exprfor,p),)
	expr2(p) = (build(2,expr1,p),)
	build(1,expr2,"")
end
macro cycle3!(i,vars)#««2
	i = esc(i)
	@assert Meta.isexpr(vars, :tuple)
	n = length(vars.args)
	@assert n % 3 == 0
	v0 = esc(vars)
	v1 = Expr(:tuple, [esc(vars.args[i+1-3*(i%3==0)]) for i in 1:n]...)
	v2 = Expr(:tuple, [esc(vars.args[i+2-3*(i%3!=1)]) for i in 1:n]...)
	quote
		@assert $i ∈ 0:2
		if $i == 1; $v0 = $v1
		elseif $i == 2; $v0 = $v2
		end
	end
end
macro permute3!(i,vars)#««2
	i = esc(i)
	@assert Meta.isexpr(vars, :tuple)
	n = length(vars.args)
	@assert n % 3 == 0
	v0 = esc(vars)
	@inline permute(f) = Expr(:tuple, [esc(vars.args[f(i)]) for i in 1:n]...)
	v1 = permute(i->i+1-3*(i%3==0))
	v2 = permute(i->i+2-3*(i%3!=1))
	v3 = permute(i->i+(i%3==2)-(i%3==0))
	v4 = permute(i->i+(i%3==1)-(i%3==2))
	v5 = permute(i->i-2(i%3==0)+2(i%3==1))
	quote
		@assert $i ∈ 0:5
		    if $i == 1; $v0 = $v1
		elseif $i == 2; $v0 = $v2
		elseif $i == 3; $v0 = $v3
		elseif $i == 4; $v0 = $v4
		elseif $i == 5; $v0 = $v5
		end
	end
end
# degeneracy test««2
"""
    degeneracy(triangle, ε²)

Returns an encoding of the degeneracy type of `triangle`, as a constant:
 - `vertex1` if vertex 1 is less than ε from opposite edge,
 - `edge1` if edge 1 has length less than ε,
 - `interior` if at least two edges are shorter than ε,
 - `invalid` if the triangle is nondegenerate.
"""
function degeneracy(triangle, ε² = 0)
	(a, b, c) = triangle
	bc=c-b; bc²=dot(bc,bc)
	ca=a-c; ca²=dot(ca,ca)
	ab=b-a; ab²=dot(ab,ab)
	if ab² ≤ ε²
		(bc² ≤ ε² || ca² ≤ ε²) && return Constants.interior # a ≈ b ≈ c
		return Constants.edge3 # a ≈ b
	elseif bc² ≤ ε²
		ca² ≤ ε² && return Constants.interior
		return Constants.edge1
	elseif ca² ≤ ε²
		return Constants.edge2
	end
	abc = norm(cross(ab, bc)); abc² = abc*abc

	abc² ≤ ε²*bc² && return Constants.vertex1 # a ∈ ]bc[
	abc² ≤ ε²*ca² && return Constants.vertex2
	abc² ≤ ε²*ab² && return Constants.vertex3
	return Constants.invalid # not degenerate
end
# 2d intersections««1
@inline det2(u,v)=u[1]*v[2]-u[2]*v[1]
@inline det2(p,q,r)=det2(q-p,r-p)
# inter_segment2_halfplane ««2
"""
    inter_segment2_halfplane(s::IntersectionData, (p,q), t)

Intersects `s` with the closed half-plane defined by `det(p,q,x) ≥ 0`,
marking intersection as if the line `(p,q)` is the `Constant` `t`.
"""
function inter_segment2_halfplane(it::IntersectionData, halfplane, t; ε=0)
	(p,q) = halfplane
	ID=typeof(it)
	isempty(it) && return it
	(u1, (a1,b1)) = it[1]
	d1 = det2(p,q,u1)
	if length(it) == 1
		d1 <-ε && return ID()
		d1 ≤ ε && (Int(b1)|Int(t) == 7) && println("b1=$b1, t=$t")
		d1 ≤ ε && return ID(u1 => (a1, b1|t))
		return ID(u1 => (a1, b1))
	end
	# length(it) == 2
	(u2, (a2,b2)) = it[2]
	d2 = det2(p,q,u2)
	# (uv) ∩ (pq) = u+(v-u)*[pqu]/[pqu-pqv]
	inter_point = @closure () -> u1 + (u2-u1) * d1/(d1-d2)
	if d1 < -ε
		d2 <-ε && return ID()
		d2 ≤ ε && return ID(u2 => (a2, b2|t))
		# opposite signs: compute intersection
		# FIXME this is wrong, should use (a
		return ID(inter_point() => (a1&a2, (b1&b2)|t), u2 => (a2, b2))
	elseif d1 ≤ ε # u1 ∈ (pq)
		d2 <-ε && return ID(u1 => (a1, b1|t))
		d2 ≤ ε && return ID(u1 => (a1, b1|t), u2 => (a2, b2|t))
		return ID(u1 => (a1, b1|t), u2 => (a2, b2))
	end
	# d1 > ε
	d2 <-ε && return ID(u1 => (a1, b1), inter_point() => (a1&a2, (b1&b2)|t))
	d2 ≤ ε && return ID(u1 => (a1, b1), u2 => (a2, b2|t))
	return ID(u1 => (a1, b1), u2 => (a2, b2))
end
# inter_segment2_triangle2««2
function inter_segment2_triangle2(seg, tri; turn=0, ε=0)
	(u1,v1) = seg
	(p2,q2,r2) = tri
	it0 = IntersectionData{2,typeof(u1)}(
		u1 => (Constants.vertex1, Constants.interior),
		v1 => (Constants.vertex2, Constants.interior))
	it1 = inter_segment2_halfplane(it0, (q2, r2), Constants.edge23; ε)
	it2 = inter_segment2_halfplane(it1, (r2, p2), Constants.edge31; ε)
	it3 = inter_segment2_halfplane(it2, (p2, q2), Constants.edge12; ε)
	return it3
end
		

# inter_triangle2 ««2
function connect2(n, pttype, points, itpq, itqr)
	isempty(itpq) && return n
	m = length(itpq) - (!isempty(itqr) && last(itpq)[2] == first(itqr)[2])
	for i in 1:m
		points[n+i] = itpq.points[i]
		pttype[n+i] = itpq.pttype[i]
	end
	return n+m
end
"""
    inter_triangle2

Computes the intersection of two 2d triangles.
"""
function inter_triangle2((p1,q1,r1),(p2,q2,r2); ε=0)
	itpq = inter_segment2_triangle2((p1,q1),(p2,q2,r2); ε)
	itqr = inter_segment2_triangle2((q1,r1),(p2,q2,r2); ε)
	itrp = inter_segment2_triangle2((r1,p1),(p2,q2,r2); ε)
	# connect those three together:
	shift1!(itqr, 1)
	shift1!(itrp, 2)
	points = similar(itqr.points, 6)
	pttype = similar(itqr.pttype, 6)
	n = connect2(0, pttype, points, itpq, itqr)
	n = connect2(n, pttype, points, itqr, itrp)
	n = connect2(n, pttype, points, itrp, itpq)
	return IntersectionData{6,eltype(points)}(n, pttype, points)
end

# 3d intersections««1
# Projectors interaction««2


# lifting intersection data
@inline function Projectors.lift(p::Projector, it::IntersectionData)
	newpoints = inv(p).(first.(it))
	return IntersectionData{allocsize(it),eltype(newpoints)}(length(it),
		last.(it), newpoints)
end

# inter(triangle1, triangle2) ««2
"""
    inter(triangle1, triangle2, ε)

Returns a description of the intersection of two 3-dimensional triangles.
This is returned as a `IntersectionData` structure.

Both arguments may be of any types, as long as that type supports enumeration
to three vertices, and those are compatible with basic geometry operations.
"""
@inline inter(tri1::SVector{3}, tri2::SVector{3}, args...) =
	inter(tri1.data, tri2.data, args...)
function inter(tri1::NTuple{3,SVector{3,T}}, tri2::NTuple{3,SVector{3,T}},
	ε=0) where{T}
	# loosely inspired by
	# [Devillers, Guigue, _Faster triangle-triangle intersection tests_;
	#   https://hal.inria.fr/inria-00072100/document]

	(p1,q1,r1) = tri1
	(p2,q2,r2) = tri2
	ID = IntersectionData{6,typeof(p1)}

	normal2 = cross(q2-p2, r2-p2)
	@assert norm(normal2, Inf) > ε "degenerate triangle2"

	dp1 = dot(normal2, p1-p2)
	dq1 = dot(normal2, q1-p2)
	dr1 = dot(normal2, r1-p2)

	# permute both triangles as needed so that t2 separates p1 from q1, r1
	# this guarantees that line (bb2 cc2) intersects segments (a1b1) and (a1c1).
# 	@debug ("signs for p1,q1,r1", Int.(sign.((dp1,dq1,dr1))))
	@tree27((dp1,dq1,dr1),
		"+++" => (return ID()),
		"0--" => begin
			@permute3! turn (p1,q1,r1)
			return inter_touch((p1,q1,r1), turn, (p2,q2,r2), normal2, ε)
		end,
		"0+-" => begin
			@permute3! turn (p1,q1,r1,dp1,dq1,dr1)
			return inter_arrow((p1,q1,r1), turn, (p2,q2,r2), flip,
				(dp1,dq1,dr1), normal2, ε)
		end,
		"+00" => begin
# 			@debug "config inter_border 1 (+00), turn=$turn, flip=$flip"
			@permute3! turn (p1,q1,r1)
		return inter_border((p1,q1,r1), turn, (p2,q2,r2), normal2, ε)
		end,
		"+--" => (i1=turn; i2=flip ? 3 : 0),
		"000" => begin 
			return inter_coplanar((p1,q1,r1), (p2,q2,r2), normal2, ε)
		end,
	)
# 	@debug "now (i1,i2)=($i1,$i2)"

	# likewise for second triangle
	normal1 = cross(q1-p1, r1-p1)
	@assert norm(normal1, Inf) > ε "degenerate triangle1"
	dp2 = dot(normal1, p2-p1)
	dq2 = dot(normal1, q2-p1)
	dr2 = dot(normal1, r2-p1)
# 	@debug "signs for dp2 dq2 dr2: $(Int.(sign.((dp2,dq2,dr2))))"
	@tree27((dp2,dq2,dr2),
		"+++" => (return ID()),
		"0--" => begin
			i2+= turn; (i2 ≥ 3) && (i2 -= 3);
			@permute3! i2 (p2,q2,r2)
			return swap!(inter_touch((p2,q2,r2), i2, (p1,q1,r1), normal1, ε))
		end,
		"0+-" => begin
			@permute3! turn (p2,q2,r2,dp2,dq2,dr2)
			return swap!(inter_arrow((p2,q2,r2), turn, (p1,q1,r1), flip,
				(dp2,dq2,dr2), normal1, ε))
		end,
		"+00" => begin
# 			@debug "config inter_border 2, turn=$turn"
			@permute3! turn (p2,q2,r2)
		return swap!(inter_border((p2,q2,r2), turn, (p1,q1,r1), normal1, ε))
		end,
		"+--" => (i2+=turn),
		"-++" => (i2+=turn; i1+=3),
	)
# 	@debug "now i1,i2 = $((i1,i2))"

	# apply both permutations ««
	(lp1, lq1, lr1) = (lp2,lq2,lr2) =
		(Constants.vertex1, Constants.vertex2, Constants.vertex3)
	@permute3! i1 (p1,q1,r1, dp1,dq1,dr1, lp1,lq1,lr1)
	@permute3! i2 (p2,q2,r2, dp2,dq2,dr2, lp2,lq2,lr2)
	# 1,2,3 represent even permutations, and 4,5,6 odd ones:
	(i1 ≥ 3) && (normal1 = -normal1; dp2 = -dp2; dq2 = -dq2; dr2 = -dr2)
	(i2 ≥ 3) && (normal2 = -normal2; dp1 = -dp1; dq1 = -dq1; dr1 = -dr1)

	# a laundry list of assertions to check that we are in a standard
	# configuration:
	@assert normal1 ≈ cross(q1-p1, r1-p1)
	@assert normal2 ≈ cross(q2-p2, r2-p2)
	@assert dot(normal2, p1-p2) ≈ dp1
	@assert dot(normal2, q1-p2) ≈ dq1
	@assert dot(normal2, r1-p2) ≈ dr1
	@assert dot(normal1, p2-p1) ≈ dp2
	@assert dot(normal1, q2-p1) ≈ dq2
	@assert dot(normal1, r2-p1) ≈ dr2
	@assert dp1 ≥ 0; @assert dq1 ≤ 0; @assert dr1 ≤ 0
	@assert dp2 ≥ 0; @assert dq2 ≤ 0; @assert dr2 ≤ 0
	# »»
	# two cases where we know that the intersection is empty:
	# dq1q2 is a determinant comparing the positions of q1 and q2
	p1p2q1 = cross(p2-p1, q1-p1)
	dq1q2 = dot(p1p2q1, q2-p1)
# 	@debug "dq1q2=$dq1q2"
	dq1q2 < -ε && return ID()

	p1p2r1 = cross(p2-p1, r1-p1)
	dr1r2 = dot(p1p2r1, r2-p1)
# 	@debug "dr1r2=$dr1r2"
	dr1r2 > ε && return ID()

	P1Q1=lp1&lq1; P1R1=lp1&lr1;
	P2Q2=lp2&lq2; P2R2=lp2&lr2;
# # 	P1Q1=lr1; P1R1=lq1; P2Q2=lr2; P2R2=lq2
	Q1 = @closure x-> barycenter(q1, p1, dp1/(dp1-dq1))=>(P1Q1,x)
	R1 = @closure x-> barycenter(r1, p1, dp1/(dp1-dr1))=>(P1R1,x)
	Q2 = @closure x-> barycenter(q2, p2, dp2/(dp2-dq2))=>(x,P2Q2)
	R2 = @closure x-> barycenter(r2, p2, dp2/(dp2-dr2))=>(x,P2R2)

	dq1r2 = dot(p1p2q1, r2-p1)
# 	@debug "dq1r2 = $dq1r2"
	if dq1r2 < -ε
		dq1q2 ≤ ε && return ID(Q2(P1Q1))
		dr1q2 = dot(p1p2r1, q2-p1)
		dr1q2 <-ε && return ID(Q1(Constants.interior),Q2(Constants.interior))
		dr1q2 ≤ ε && return ID(Q1(Constants.interior),Q2(P1R1))
		             return ID(Q1(Constants.interior),R1(Constants.interior))
	elseif dq1r2 ≤ ε # cc2 ∈ edge [a1,q1]
		@assert collinear3(p1,R2(Constants.interior)[1],q1)
		@assert monotonic(p1,R2(Constants.interior)[1],q1)
		dr1q2 = dot(p1p2r1, q2-p1)
		dr1q2 <-ε && return ID(R2(P1Q1),Q2(Constants.interior))
		dr1q2 ≤ ε && return ID(R2(P1Q1),Q2(P1R1))
		             return ID(R2(P1Q1),R1(Constants.interior))
	elseif dr1r2 <-ε
		dr1q2 = dot(p1p2r1, q2-p1)
		dr1q2 <-ε && return ID(R2(Constants.interior),Q2(Constants.interior))
		dr1q2 ≤ ε && return ID(R2(Constants.interior),Q2(P1R1))
		             return ID(R2(Constants.interior),R1(Constants.interior))
	end
	@assert dr1r2 ≤ ε # empty case was already detected above
	return ID(R2(P1R1)) # cc2 ∈ edge [a1,r1]
end
# inter_touch««2
"""
    inter_touch

Returns intersection in the case where triangle `t1` touches plane of `t2`
at the single point `p1`
# (oriented so that q1, r1 are *below* the plane).
"""
function inter_touch((p1,q1,r1), i1, (p2,q2,r2), normal2, ε)
	# Possible cases here:
	# - nothing
	# - vertex-vertex (return empty)
	# - vertex1 on edge2
	# - vertex1 inside face2
	# To decide: project everything on the plane of face2
	# and compute three 2d determinants
	proj = Projector(normal2)
# 	@debug "entering inter_touch $i1"
	ID=IntersectionData{6,typeof(p1)}

	a1 = proj(p1)
	(a2,b2,c2) = proj.((p2,q2,r2))
	dabc = abs(normal2[abs(proj.dir)])
	@assert det2(a2,b2,c2) ≈ dabc

	dab = cross(b2-a2, a1-a2)
	dab < -ε && return ID()
	dbc = cross(c2-b2, a1-b2)
	dbc < -ε && return ID()
	dca = dabc - dab - dbc
	@assert abs(dca-cross(a2-c2,a1-c2)) < 1e-9
	dca < -ε && return ID()

# 	@debug "point $a1 in ($a2,$b2,$c2): $dab, $dbc, $dca"
	if dab ≤ ε # it cannot be <-ε, so a1 ∈ [a2,b2]
		if dbc ≤ ε
			@assert samepoint(a1,b2)
			return ID(p1 => (Constants.vertex1<<i1, Constants.vertex2))
		elseif dca ≤ ε
			@assert samepoint(a1,a2)
			return ID(p1 => (Constants.vertex1<<i1, Constants.vertex1))
		end
		# a1 ∈ open segment ]a2, b2[
		return ID(p1 => (Constants.vertex1<<i1, Constants.edge12))
	elseif dbc ≤ ε # a1 ∈ ]b2,c2] (b2 was already treated above)
		if dca ≤ ε # a1 == c2
			@assert samepoint(a1,c2)
			return ID(p1 => (Constants.vertex1<<i1, Constants.vertex3))
		end # open segment ]b2,c2[
		return ID(p1 => (Constants.vertex1<<i1, Constants.edge23))
	elseif dca ≤ ε
		return ID(p1 => (Constants.vertex1<<i1, Constants.edge31))
	end
	# a1 ∈ interior of other face
	return ID(p1 => (Constants.vertex1<<i1, Constants.interior))
end
# inter_arrow ««2
"""
    inter_arrow

Computes intersection, assuming that the plane of one triangle cuts the
other through vertex `p1`.
"""
function inter_arrow((p1,q1,r1), turn, (p2,q2,r2), flip,
	(zp1, zq1, zr1), normal2, ε)
	proj = Projector(normal2, p2)
	ID = IntersectionData{6,typeof(p1)}
	@assert abs(zp1) ≤ ε

	(a1,b1,c1,a2,b2,c2) = proj.((p1,q1,r1,p2,q2,r2))
	dpqr = abs(normal2[abs(proj.dir)])
	# we know that interior of segment (q1,r1) intersects plane 2:
	v = barycenter(b1, c1, zr1/(zr1-zq1))
	it = inter_segment2_triangle2((a1, v), (a2,b2,c2); ε)
	length(it) == 0 && return ID()

	(pt1, (u1, v1)) = it[1]; z1 = inv(proj)(pt1)
	w1 = (u1 == Constants.edge12) ? Constants.interior :
		(u1 == Constants.vertex1) ? Constants.vertex1<<turn :
		Constants.edge23<<turn
	length(it) == 1 && return ID(z1 => (w1, v1))

	(pt2, (u2,v2)) = it[2]; z2 = inv(proj)(pt2)
	w2 = (u2 == Constants.edge12) ? Constants.interior :
		(u2 == Constants.vertex1) ? Constants.vertex1<<turn :
		Constants.edge23<<turn
	return ID(z1 => (w1, v1), z2 => (w2, v2))
end
# inter_border ««2
"""
    inter_border

Computes intersection, assuming that (exactly) vertices q1, r1
are on the plane of triangle 2.
"""
function inter_border((p1,q1,r1), i1, (p2,q2,r2), normal2, ε)
	ID = IntersectionData{6,typeof(p1)}
	proj = Projector(normal2, p2)
	(u1,v1,a2,b2,c2) = proj.((q1,r1,p2,q2,r2))
	dpqr = abs(normal2[abs(proj.dir)])
	it = inter_segment2_triangle2((u1,v1), (a2,b2,c2); ε)
	it3 = ID(undef, length(it))
	for i in 1:length(it)
		it3.pttype[i] = it.pttype[i]
		it3.points[i] = inv(proj)(it.points[i])
	end
	rename1!(it3, Constants.edge12 => Constants.edge23<<i1,
		Constants.vertex1 => Constants.vertex2<<i1,
		Constants.vertex2 => Constants.vertex3<<i1)
# 	rename!(it, 1, (i1, plus1mod3(i1,3), plus2mod3(i1,3)))
	return it3
end
# inter_coplanar ««2

@inline function inter_coplanar((p1,q1,r1), (p2,q2,r2), normal2, ε)
	proj = Projector(normal2, p2)
	(a1,b1,c1,a2,b2,c2) = proj.((p1,q1,r1,p2,q2,r2))
	dpqr = abs(normal2[abs(proj.dir)])
	it = inter_triangle2((a1,b1,c1),(a2,b2,c2); ε)
	return inv(proj)(it)
end
#»»1
end # module
