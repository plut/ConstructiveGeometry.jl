"""
    TriangleIntersections

Routines for computing the intersection of two 2d or 3d triangles.
The triangle vertices may be of any type which supports:
indexation, `cross`, `dot`, and interaction with `Projector` type.
"""
module TriangleIntersections
using LinearAlgebra
using StaticArrays

"""
    TriangleIntersections.Position

Encoding of position in a triangle.
"""
module Position#««
# const Encoding = Int8
# 	interior = 0
# 	interior=0
# 	vertex1=1
# 	vertex2=2
# 	vertex3=3
# 	edge23=4
# 	edge31=5
# 	edge12=6
# @inline isinterior(a::Encoding) = a == interior
# @inline isvertex(a::Encoding) = a ≤ 3
# @inline isedge(a::Encoding) = a ≥ 4
# @inline same_edge(a::Encoding, b::Encoding) =
# @inline next(a::Encoding) = typeof(a)(Integer(a)<<1|Integer(a)>>2)
# @inline prev(a::Encoding) = typeof(a)(Integer(a)<<2|Integer(a)>>1)

@enum Encoding::Int8 begin
	interior=0
	vertex1=2|4
	vertex2=4|1
	vertex3=1|2
	edge23=1
	edge31=2
	edge12=4
	invalid=-1
end
@inline Base.iszero(a::Encoding) = a == interior
@inline isvertex(a::Encoding) = count_ones(Integer(a)) == 1
@inline isedge(a::Encoding) = count_ones(Integer(a)) == 2
# @inline same_edge(a::Encoding, b::Encoding) = (Integer(a)&Integer(b)) ≠ 0
@inline Base.:&(a::Encoding, b::Encoding) = Encoding(Integer(a)&Integer(b))
@inline Base.:|(a::Encoding, b::Encoding) = Encoding(Integer(a)|Integer(b))
@inline Base.:^(a::Encoding, b::Encoding) = Encoding(Integer(a)^Integer(b))
@inline function Base.:<<(a::Encoding, i::Integer)
	i == 0 && return a
	i == 1 && return (typeof(a))((Integer(a)<<1|Integer(a)>>2)&7)
	i == 2 && return (typeof(a))((Integer(a)<<2|Integer(a)>>1)&7)
	return a<<(i%3)
end

end#»»
# Preamble««1
# Generic stuff««2
# predicates for assertions:
@inline barycenter(p1, p2, λ) = p2 + λ*(p1-p2) # λp1+(1-λ)p2
@inline collinear3(a,b,c) = norm(cross(a-b,c-b),1) ≤ 1e-8
@inline collinear2(a,b,c) = abs(det2(a,b,c)) ≤ 1e-8
@inline monotonic(a,b,c) = dot(a-b,c-b) ≤ 0
@inline samepoint(a,b) = norm(a-b,1) ≤ 1e-8

_THICKNESS=1e-8
# @inline plus1mod3(i,j=0)=(i==1) ? 2+j : (i==2) ? 3+j : 1+j
# @inline plus2mod3(i,j=0)=(i==1) ? 3+j : (i==2) ? 1+j : 2+j

# struct IntersectionData««2
"""
    struct IntersectionData

Describes an intersection. Each point of the returned polygon
is marked with its type w.r.t the two input simplexes.
"""
struct IntersectionData{N,P} <:
		AbstractVector{Tuple{P,NTuple{2,Position.Encoding}}}
# maximum of N points, constant-size
	npoints::Int
	pttype::MVector{N,NTuple{2,Position.Encoding}}
	points::MVector{N,P}
	@inline IntersectionData{N,P}(npoints::Int, ::UndefInitializer) where{N,P} =
		new{N,P}(npoints, MVector{N,NTuple{2,Position.Encoding}}(undef),
			MVector{N,P}(undef))
	@inline IntersectionData{N,P}(npoints::Int, pttype::MVector{N,P},
		points::MVector{N,P}) where{N,P} =
		new{N,P}(npoints, pttype, points)
	@inline function IntersectionData{N,P}(npoints::Int, pttype, points
		) where{N,P}
		it = IntersectionData{N,P}(npoints, undef)
		# instead of broadcasting, we write the loop: this allows us to
		# over-specify the set of points
		for i in 1:npoints
			it.pttype[i] = pttype[i]; it.points[i] = points[i]
		end
		return it
	end
	@inline IntersectionData{N,P}() where{N,P} = IntersectionData{N,P}(0, undef)
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
		ij=esc(:(turn, direct))
		haskey(expr, "$a$b$c") && return (:($ij=(0,true)),expr["$a$b$c"]...)
		haskey(expr, "$b$c$a") && return (:($ij=(1,true)),expr["$b$c$a"]...)
		haskey(expr, "$c$a$b") && return (:($ij=(2,true)),expr["$c$a$b"]...)
		haskey(expr, "$d$e$f") && return (:($ij=(0,false)),expr["$d$e$f"]...)
		haskey(expr, "$e$f$d") && return (:($ij=(1,false)),expr["$e$f$d"]...)
		haskey(expr, "$f$d$e") && return (:($ij=(2,false)),expr["$f$d$e"]...)
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
# 2d intersections««1
@inline det2(u,v)=u[1]*v[2]-u[2]*v[1]
@inline det2(p,q,r)=det2(q-p,r-p)
# inter2_segment_triangle««2
"""
    Computes the intersection between 2d segment and triangle

Returns an `IntersectionData` structure, encoded in the following way:
 - if there are 2 points, they are guaranteed to be in the same order as the segments;
 - on segment, (1,2) are vertices (u,v), 0 is interior
 - on triangle, (i2,i2+1,i2+2)(mod3) are edges (qr, rp, pq),
   and likewise (+3) for vertices (p,q,r); 0 is interior
"""
function inter2_segment_triangle((u1,v1),(p2,q2,r2),i2=0;
	dpqr = missing, ε=_THICKNESS)
	(dpqr isa Missing) && (dpqr = det2(p2,q2,r2)) # FIXME, make something cleaner
	@debug "segment-triangle intersection\n($u1,$v1)\n($p2,$q2,$r2)\nd=$dpqr"
	ID=IntersectionData{2,typeof(u1)}
	# compute position of u:
	dpq_u = det2(p2,q2,u1)
	dqr_u = det2(q2,r2,u1)
	drp_u = dpqr - dpq_u - dqr_u
	@assert abs(drp_u-det2(r2,p2,u1)) ≤ 1e-9
	# likewise for v:
	dpq_v = det2(p2,q2,v1)
	dqr_v = det2(q2,r2,v1)
	drp_v = dpqr - dpq_v - dqr_v
	# easy reject cases: segment lies fully outside one side of triangle
	dpq_u < -ε && dpq_v < -ε && return ID()
	dqr_u < -ε && dqr_v < -ε && return ID()
	drp_u < -ε && drp_v < -ε && return ID()

	# symbolic names for returned points: ««3
	@inline P2(x)=p2=>(x,Position.vertex1<<i2)
	@inline Q2(x)=q2=>(x,Position.vertex2<<i2)
	@inline R2(x)=r2=>(x,Position.vertex3<<i2)
	@inline U1(x)=u1=>(Position.vertex1,x)
	@inline V1(x)=v1=>(Position.vertex2,x)
	@inline QR() = Position.edge23<<i2
	@inline RP() = Position.edge31<<i2
	@inline PQ() = Position.edge12<<i2
	@inline UVQR()= barycenter(v1, u1, dqr_u/(dqr_u-dqr_v))=>
		(Position.interior,QR())
	@inline UVRP()= barycenter(v1, u1, drp_u/(drp_u-drp_v))=>
		(Position.interior,RP())
	@inline UVPQ()= barycenter(v1, u1, dpq_u/(dpq_u-dpq_v))=>
		(Position.interior,PQ())

	# rotate triangle to standard configuration: ««3
	@debug "position of u: $(Int.(sign.((dqr_u,drp_u,dpq_u))))"
	@tree27((dqr_u, drp_u, dpq_u),
		"+++" => (@goto inside),
		"-++" => begin
			@cycle3! turn (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
			i2+= turn; (i2≥3) && (i2-=3)
			@goto outside_p
		end,
		"-+0" => begin
			@cycle3! turn (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
			i2+= turn; (i2≥3) && (i2-=3)
			@goto outside_p
		end,
		"0++" => begin
			@cycle3! turn (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
			i2+= turn; (i2≥3) && (i2-=3)
			@goto edge_qr
		end,
		"+00" => begin
			@cycle3! turn (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
			i2+= turn; (i2≥3) && (i2-=3)
			@goto vertex_p
		end,
		# leave those here to prevent @tree27 from converting - to +:
		"000" => (error("impossible case")),
		"00-" => (error("impossible case")),
		"0--" => (error("impossible case")),
		"---" => (error("impossible case")),
	)
	@assert false "not reachable"

	@label inside #««3
	@tree27((dqr_v, drp_v, dpq_v),
		"+++" => (return ID(U1(Position.interior), V1(Position.interior))),
		"000" => (error("impossible case")),
		"00-" => (error("impossible case")),
		"0--" => (error("impossible case")),
		"---" => (error("impossible case")),
		# v is a triangle vertex:
		"00+" => (return ID(U1(Position.interior),
			turn == 0 ? P2(Position.vertex2) :
			turn == 1 ? Q2(Position.vertex2) :
			         R2(Position.vertex2))),
		# v is inside a triangle edge:
		"0++" => (return ID(U1(Position.interior), V1(Position.edge23 << turn))),
		"-0+" => (return ID(U1(Position.interior),
				turn==0 ? UVQR()  : turn == 1 ? UVRP() : UVPQ())),
		"-++" => (return ID(U1(Position.interior),
				turn==0 ? UVQR()  : turn == 1 ? UVRP() : UVPQ())),
		"--+" => begin
		@cycle3! turn (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
		i2+= turn; (i2≥3) && (i2-=3)
		# which side is crossed depends on ruv
		druv = det2(r2,u1,v1)
		druv > ε && return ID(U1(Position.interior), UVQR())
		druv <-ε && return ID(U1(Position.interior), UVRP())
		return ID(U1(Position.interior), R2(Position.interior))
		end
	)
	@assert false "not reachable"

	@label vertex_p #««3
	@assert samepoint(u1, p2)
	dpq_v < -ε && return ID(P2(Position.vertex1))
	drp_v < -ε && return ID(P2(Position.vertex1))
	if dpq_v ≤ ε
		dqr_v > ε && return ID(P2(Position.vertex1), V1(PQ()))
		dqr_v <-ε && return ID(P2(Position.vertex1), Q2(Position.interior))
		@assert samepoint(v1, q2)
		return ID(P2(Position.vertex1), Q2(Position.vertex2))
	end
	if drp_v ≤ ε
		dqr_v > ε && return ID(P2(Position.vertex1), V1(RP()))
		dqr_v <-ε && return ID(P2(Position.vertex1), R2(Position.interior))
		@assert samepoint(v1, r2)
		return ID(P2(Position.vertex1), R2(Position.vertex2))
	end
	dqr_v > ε && return ID(P2(Position.vertex1), V1(Position.interior))
	dqr_v <-ε && return ID(P2(Position.vertex1), UVQR())
	@assert collinear2(q,v,r)
	@assert monotonic(q,v,r)
	return ID(P2(Position.vertex1), V1(QR()))

	@label edge_qr #««3
	@debug "edge configuration; rotating(i=$i)..."
	@debug "position of u=$(Int.(sign.((dqr_u,drp_u,dpq_u))))\nposition of v=$(Int.(sign.((dqr_v,drp_v,dpq_v))))"
	@assert dqr_v ≥ -ε
	if dqr_v ≤ ε # u, v on same edge qr
		dpq_v <-ε && return ID(U1(QR()), Q2(Position.interior))
		drp_v <-ε && return ID(U1(QR()), R2(Position.interior))
		dpq_v ≤ ε && return ID(U1(QR()), Q2(Position.vertex2))
		drp_v ≤ ε && return ID(U1(QR()), R2(Position.vertex2))
		return ID(U1(QR()), V1(QR()))
	end
	@assert dqr_v > ε
	dpq_v > ε && drp_v > ε && return ID(U1(QR()), V1(Position.interior))
	dpuv = det2(u1,v1,p2)
	dpuv > ε && return ID(U1(QR()), UVRP())
	dpuv <-ε && return ID(U1(QR()), UVPQ())
	@assert collinear2(u1,v1,p2)
	dpq_v <-ε && return ID(U1(QR()), P2(Position.interior))
	@assert samepoint(v1,p2)
	return ID(U1(QR()), P2(Position.vertex2))

	@label outside_p #««3
	@assert dqr_u < -ε # u1 is far from p
	@assert dpq_u > ε
	@assert dqr_v ≥ -ε
	
	if dqr_v ≤ ε # edge case: v1 on (qr) line««
		if dpq_v > ε
			drp_v > ε && return ID(V1(QR())) # ]qr[ segment
			drp_v <-ε && return ID() # outside on r-side
			@assert samepoint(v1, r2)
			return ID(R2(Position.vertex2))
		end
		dpq_v <-ε && return ID() # outside on q-side
		@assert samepoint(v1, q2)
		return ID(Q2(Position.vertex2))
	end#»»
	# generic case: v1 is on close side of (qr)
	# we might cross (in this order) two of (qr), (rp), (pq)
	@debug "position of v: $dqr_v $drp_v $dpq_v"
	@assert dqr_v > ε
	# start by examining what happens on the u side:
	# given that ]uv[ crosses (qr), it crosses [qr[ iff ruv < 0
	druv = det2(r2,u1,v1)
	if druv > ε # don't cross [qr[
		# easy case: if u is on the -++ side and uv misses r:
		drp_u ≥ -ε && return ID()
		@assert drp_v ≥ -ε # because of easy reject
		# since drp_u <-ε and drp_v ≥-ε, uv crosses rp. Where does this happen?
		if drp_v > ε # generic case: v is (stricly) inside (rp).
			dpuv = drp_u - drp_v + druv
			dpuv > ε && return ID() # missed p and the triangle
			dpuv ≥-ε && return ID(P2(Position.interior)) # touched p

			@assert dpuv <-ε # we hit the ]rp[ side, and thus the triangle
			dpq_v > ε && return ID(UVRP(), V1(Position.interior)) # v is an inner point
			dpq_v <-ε && return ID(UVRP(), UVPQ()) # v is outside (pq)
			@assert collinear2(p2,v1,q2)
			@assert monotonic(p2,v1,q2)
			return ID(UVRP(), V1(PQ()))
		end
		@assert collinear2(r2,p2,v1)
		dpq_v > ε && return ID(V1(PQ()))
		dpq_v < ε && return ID()
		@assert samepoint(v1,p2)
		return ID(P2(Position.vertex2))
	elseif druv ≥ -ε # (uv) touches r
		@assert collinear2(u1,r2,v1)
		@assert monotonic(u1,r2,v1)
		drp_v <-ε && return ID(R2(Position.interior))
		if drp_v ≤ ε # (uv) == (rp)
			@assert collinear2(u1,v1,p2)
			dpq_v > ε && return ID(R2(Position.interior),V1(RP()))
			dpq_v <-ε && return ID(R2(Position.interior),P2(Position.interior))
			@assert samepoint(v1,p2)
			return ID(R2(Position.interior),P2(Position.vertex2))
		end
		@assert drp_v > ε # (uv) might yet cross (pq):
		dpq_v > ε && return ID(R2(Position.interior),V1(Position.interior))
		dpq_v < ε && return ID(R2(Position.interior),UVPQ())
		@assert collinear2(p2,q2,v1)
		return ID(R2(Position.interior),V1(PQ()))
	end
	# now check on the q-side:
	dquv = dqr_v - dqr_u + druv
	dquv < -ε && return ID() # we missed on the q-side
	if dquv ≤ ε # (uv) touches q
		@assert collinear2(u1,q2,v1)
		@assert monotonic(u1,q2,v1)
		@assert dpq_v ≤ ε
		dpq_v <-ε && return ID(Q2(Position.interior)) # far
		@assert samepoint(v1,q2)
		return ID(Q2(Position.vertex2))
	end
	@assert dquv > ε # uv crosses the ]q,r[ segment...
	if dpq_v > ε
		drp_v > ε && return ID(UVQR(), V1(Position.interior)) # inner point
		drp_v <-ε && return ID(UVQR(), UVRP())
		@assert collinear2(r2,v1,p2)
		@assert monotonic(r2,v1,p2)
		return ID(UVQR(), V1(RP()))
	elseif dpq_v ≥-ε
		drp_v > ε && return ID(UVQR(), V1(PQ())) # pq edge
		drp_v <-ε && return ID(UVQR(), UVRP()) # outside p on pq edge
		@assert samepoint(v1,p2)
		return ID(UVQR(), P2(Position.vertex2)) # 
	end
	@assert dpq_v <-ε
	# easy reject case where we don't cross rp:
	drp_v ≥-ε && return ID(UVQR(), UVPQ())
	@assert drp_v <-ε
	# worst case: we cross diagonally from (-++) to (+--) region;
	# decide which way around p we go:
	dpuv = drp_u - drp_v + druv
	dpuv > ε && return ID(UVQR(), UVRP())
	dpuv <-ε && return ID(UVQR(), UVPQ())
	@assert collinear2(u1,p2,v1)
	@assert monotonic(u1,p2,v1)
	return ID(UVQR(), P2(Position.interior))
end

# inter2 ««2
"""
    inter2

Computes the intersection of two 2d triangles.
dpqr is the determinant of second triangle (assumed > 0).
"""
function inter2((p1,q1,r1),(p2,q2,r2), dpqr, ε=_THICKNESS)
	@assert dpqr ≈ det2(p2,q2,r2)
	@assert dpqr > ε

	itpq = inter2_segment_triangle((p1,q1),(p2,q2,r2); dpqr, ε)
	itqr = inter2_segment_triangle((q1,r1),(p2,q2,r2); dpqr, ε)
	itrp = inter2_segment_triangle((r1,p1),(p2,q2,r2); dpqr, ε)
	# glue those three together:
	rename1!(itpq, Position.interior=>Position.edge12,
		Position.vertex1=>Position.vertex1,
		Position.vertex2=>Position.vertex2)
	rename1!(itqr, Position.interior=>Position.edge23,
		Position.vertex1=>Position.vertex2,
		Position.vertex2=>Position.vertex3)
	rename1!(itrp, Position.interior=>Position.edge31,
		Position.vertex1=>Position.vertex3,
		Position.vertex2=>Position.vertex1)
# 	rename!(itpq, 1, (3,4,5))
# 	rename!(itqr, 1, (1,5,6))
# 	rename!(itrp, 1, (2,6,4))
	points = similar(itqr.points, 6)
	pttype = similar(itqr.pttype, 6)
	n = 0

	@inline function push!(x, j::Int)
		@inbounds for i in 1:length(x)-j
			pttype[n+i] = x.pttype[i]
			points[n+i] = x.points[i]
		end
		n += (length(x)-j)
	end
	push!(itpq, last(itpq)[2] == first(itqr)[2] ? 1 : 0)
	push!(itqr, last(itqr)[2] == first(itrp)[2] ? 1 : 0)
	push!(itrp, last(itrp)[2] == first(itpq)[2] ? 1 : 0)
	return IntersectionData{6,eltype(points)}(n, pttype, points)
end
# 3d intersections««1
# basic geometry stuff««2


# Projector ««2
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
@inline lift(p::Projector, u::StaticVector{2}) = SVector(lift(p,(u...,)))
@inline project(p::Projector, u::StaticVector{3}) = SVector(project(p,(u...,)))

# lifting intersection data
@inline function lift(p::Projector, it::IntersectionData)
	newpoints = inv(p).(first.(it))
	return IntersectionData{allocsize(it),eltype(newpoints)}(length(it),
		last.(it), newpoints)
end

# inter(triangle1, triangle2) ««2
"""
    inter(triangle1, triangle2, ε)

Returns a description of the intersection of two 3-dimensional triangles.
This is returned as a `IntersectionData` structure,
encoded in the following way:
 - 0: point is interior to one triangle
 - 1,2,3: point is on an edge qr, rp, pq
 - 4,5,6: point is on a vertex p, q, r

Both arguments may be of any types, as long as that type supports enumeration
to three vertices, and those are compatible with basic geometry operations.
"""
function inter((p1,q1,r1),(p2,q2,r2), ε=_THICKNESS)
	# loosely inspired by
	# [Devillers, Guigue, _Faster triangle-triangle intersection tests_;
	#   https://hal.inria.fr/inria-00072100/document]
	@debug "computing triangle-triangle intersection\n($p1,$q1,$r1)\n($p2,$q2,$r2)"

	# return type:
	ID = IntersectionData{6,typeof(p1)}

	normal2= cross(q2-p2, r2-p2)

	dp1 = dot(normal2, p1-p2)
	dq1 = dot(normal2, q1-p2)
	dr1 = dot(normal2, r1-p2)

	# permute both triangles as needed so that t2 separates p1 from q1, r1
	# this guarantees that line (bb2 cc2) intersects segments (a1b1) and (a1c1).
	@debug "signs for p1,q1,r1: $(Int.(sign.((dp1,dq1,dr1))))"
	@tree27((dp1,dq1,dr1),
		"+++" => (return ID()),
		"0--" => begin
			@permute3! turn (p1,q1,r1)
			return inter_touch((p1,q1,r1), turn, (p2,q2,r2), normal2, ε)
		end,
		"0+-" => begin
			@debug "config inter_arrow 1 with ($p1,$q1,$r1,$dp1,$dq1,$dr1) turn=$turn"
			@permute3! turn (p1,q1,r1,dp1,dq1,dr1)
# 			(direct > 3) && (normal2 = -normal2; dp1 = -dp1; dq1 = -dq1; dr1 = -dr1)
# 			@permute3! direct (p2,q2,r2)
			@debug " -> ($p1,$q1,$r1,$dp1,$dq1,$dr1)"
		return inter_arrow((p1,q1,r1), turn, (p2,q2,r2), direct,
			(dp1,dq1,dr1), normal2, ε)
		end,
		"+00" => begin
			@debug "config inter_border 1, turn=$turn"
			@permute3! turn (p1,q1,r1)
		return inter_border((p1,q1,r1), turn, (p2,q2,r2), normal2, ε)
		end,
		"+--" => (i1=turn; i2=direct ? 0 : 3),
		"000" => begin 
			return inter_coplanar((p1,q1,r1), (p2,q2,r2), normal2, ε)
		end,
	)
	@debug "now (i1,i2)=($i1,$i2)"

	# likewise for second triangle
	normal1 = cross(q1-p1, r1-p1)
	dp2 = dot(normal1, p2-p1)
	dq2 = dot(normal1, q2-p1)
	dr2 = dot(normal1, r2-p1)
	@debug "signs for dp2 dq2 dr2: $(Int.(sign.((dp2,dq2,dr2))))"
	@tree27((dp2,dq2,dr2),
		"+++" => (return ID()),
		"0--" => begin
			i2+= turn; (i2 ≥ 3) && (i2 -= 3);
			println("now i2=$i2")
			@permute3! i2 (p2,q2,r2)
			return swap!(inter_touch((p2,q2,r2), i2, (p1,q1,r1), normal1, ε))
		end,
		"0+-" => begin
			@permute3! turn (p2,q2,r2,dp2,dq2,dr2)
			return swap!(inter_arrow((p2,q2,r2), turn, (p1,q1,r1), direct,
				(dp2,dq2,dr2), normal1, ε))
		end,
		"+00" => begin
			@assert false
		end,
		"+--" => (i2+=turn),
		"-++" => (i2+=turn; i1+=3),
	)
	@debug "now i1,i2 = $((i1,i2))"

	# apply both permutations ««
	(lp1, lq1, lr1) = (lp2,lq2,lr2) =
		(Position.vertex1, Position.vertex2, Position.vertex3)
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
	@debug "dq1q2=$dq1q2"
	dq1q2 < -ε && return ID()

	p1p2r1 = cross(p2-p1, r1-p1)
	dr1r2 = dot(p1p2r1, r2-p1)
	@debug "dr1r2=$dr1r2"
	dr1r2 > ε && return ID()

	P1Q1=lp1&lq1; P1R1=lp1&lr1;
	P2Q2=lp2&lq2; P2R2=lp2&lr2;
# 	P1Q1=lr1; P1R1=lq1; P2Q2=lr2; P2R2=lq2
	@inline Q1(x)= barycenter(q1, p1, dp1/(dp1-dq1))=>(P1Q1,x)
	@inline R1(x)= barycenter(r1, p1, dp1/(dp1-dr1))=>(P1R1,x)
	@inline Q2(x)= barycenter(q2, p2, dp2/(dp2-dq2))=>(x,P2Q2)
	@inline R2(x)= barycenter(r2, p2, dp2/(dp2-dr2))=>(x,P2R2)

	dq1r2 = dot(p1p2q1, r2-p1)
	@debug "dq1r2 = $dq1r2"
	if dq1r2 < -ε
		dq1q2 ≤ ε && return ID(Q2(P1Q1))
		dr1q2 = dot(p1p2r1, q2-p1)
		dr1q2 <-ε && return ID(Q1(Position.interior),Q2(Position.interior))
		dr1q2 ≤ ε && return ID(Q1(Position.interior),Q2(P1R1))
		             return ID(Q1(Position.interior),R1(Position.interior))
	elseif dq1r2 ≤ ε # cc2 ∈ edge [a1,q1]
		@assert collinear3(p1,R2(Position.interior)[1],q1)
		@assert monotonic(p1,R2(Position.interior)[1],q1)
		dr1q2 = dot(p1p2r1, q2-p1)
		dr1q2 <-ε && return ID(R2(P1Q1),Q2(Position.interior))
		dr1q2 ≤ ε && return ID(R2(P1Q1),Q2(P1R1))
		             return ID(R2(P1Q1),R1(Position.interior))
	elseif dr1r2 <-ε
		dr1q2 = dot(p1p2r1, q2-p1)
		dr1q2 <-ε && return ID(R2(Position.interior),Q2(Position.interior))
		dr1q2 ≤ ε && return ID(R2(Position.interior),Q2(P1R1))
		             return ID(R2(Position.interior),R1(Position.interior))
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
	@debug "entering inter_touch $i1"
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

	@debug "point $a1 in ($a2,$b2,$c2): $dab, $dbc, $dca"
	if dab ≤ ε # it cannot be <-ε, so a1 ∈ [a2,b2]
		if dbc ≤ ε
			@assert samepoint(a1,b2)
			return ID(p1 => (Position.vertex1<<i1, Position.vertex2))
		elseif dca ≤ ε
			@assert samepoint(a1,a2)
			return ID(p1 => (Position.vertex1<<i1, Position.vertex1))
		end
		# a1 ∈ open segment ]a2, b2[
		return ID(p1 => (Position.vertex1<<i1, Position.edge12))
	elseif dbc ≤ ε # a1 ∈ ]b2,c2] (b2 was already treated above)
		if dca ≤ ε # a1 == c2
			@assert samepoint(a1,c2)
			return ID(p1 => (Position.vertex1<<i1, Position.vertex3))
		end # open segment ]b2,c2[
		return ID(p1 => (Position.vertex1<<i1, Position.edge23))
	elseif dca ≤ ε
		return ID(p1 => (Position.vertex1<<i1, Position.edge31))
	end
	# a1 ∈ interior of other face
	return ID(p1 => (Position.vertex1<<i1, Position.interior))
end
# inter_arrow ««2
"""
    inter_arrow

Computes intersection, assuming that the plane of one triangle cuts the
other through vertex `p1`.
"""
function inter_arrow((p1,q1,r1), i1, (p2,q2,r2), j,
	(zp1, zq1, zr1), normal2, ε)
	proj = Projector(normal2, p2)
	ID = IntersectionData{6,typeof(p1)}
# 	@assert j ∈ (1,4)
	if j; @assert zq1 > ε; @assert zr1 <-ε
	else; @assert zq1 <-ε; @assert zr1 > ε
	end
	@assert abs(zp1) ≤ ε
# 	@assert zq1 > ε
# 	@assert zr1 <-ε
# 	@debug "arrow($i1,$i2): $p2,$q2,$r2; $normal2 -> $(proj.dir)"

	(a1,b1,c1,a2,b2,c2) = proj.((p1,q1,r1,p2,q2,r2))
	@debug "$a2,$b2,$c2"
	dpqr = abs(normal2[abs(proj.dir)])
	# we know that interior of segment (q1,r1) intersects plane 2:
	v = barycenter(b1, c1, zr1/(zr1-zq1))
	it = inter2_segment_triangle((a1, v), (a2,b2,c2); dpqr, ε)
	length(it) == 0 && return ID()

	(pt1, (u1, v1)) = it[1]; z1 = inv(proj)(pt1)
	w1 = (u1 == Position.interior) ? u1 :
		(u1 == Position.vertex1) ? Position.vertex1<<i1 :
		Position.edge23<<i1
# 	w1 = (u1 == 0) ? 0 : (u1 == 1) ? (i1+3) : i1
	length(it) == 1 && return ID(z1 => (w1, v1))

	(pt2, (u2,v2)) = it[2]; z2 = inv(proj)(pt2)
	w2 = (u2 == Position.interior) ? u2 :
		(u2 == Position.vertex1) ? Position.vertex1<<i1 :
		Position.edge23<<i1
# 	w2 = (u2 == 0) ? 0 : (u2 == 1) ? (i1+3) : i1
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
	@debug "in inter_border\n($u1,$v1)\n($p2,$q2,$r2)\nproj=$proj"
	dpqr = abs(normal2[abs(proj.dir)])
	it = inter2_segment_triangle((u1,v1), (a2,b2,c2), 0; dpqr, ε)
	rename1!(it, Position.interior => Position.edge23<<i1,
		Position.vertex1 => Position.vertex2<<i1,
		Position.vertex2 => Position.vertex3<<i1)
# 	rename!(it, 1, (i1, plus1mod3(i1,3), plus2mod3(i1,3)))
	return inv(proj)(it)
end
# inter_coplanar ««2

@inline function inter_coplanar((p1,q1,r1), (p2,q2,r2), normal2, ε)
	ID = IntersectionData{6,typeof(p1)}
	proj = Projector(normal2, p2)
	(a1,b1,c1,a2,b2,c2) = proj.((p1,q1,r1,p2,q2,r2))
	dpqr = abs(normal2[abs(proj.dir)])
	it = inter2((a1,b1,c1),(a2,b2,c2), dpqr, ε)
	return inv(proj)(it)
end
#»»1
end # module
