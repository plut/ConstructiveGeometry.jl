module TriangleIntersections
# import Meshes: Point, Triangle, coordtype, vertices, coordinates
using LinearAlgebra
using StaticArrays

# Preamble««1
# Generic stuff««2
# predicates for assertions:
@inline barycenter(p1, p2, λ) = p2 + λ*(p1-p2) # λp1+(1-λ)p2
@inline collinear3(a,b,c) = norm(cross(a-b,c-b),1) ≤ 1e-8
@inline collinear2(a,b,c) = abs(det2(a,b,c)) ≤ 1e-8
@inline monotonic(a,b,c) = dot(a-b,c-b) ≤ 0
@inline samepoint(a,b) = norm(a-b,1) ≤ 1e-8

@inline real_type(x...)=Float64
_THICKNESS=1e-8
@inline plus1mod3(i,j=0)=(i==1) ? 2+j : (i==2) ? 3+j : 1+j
@inline plus2mod3(i,j=0)=(i==1) ? 3+j : (i==2) ? 1+j : 2+j

# struct IntersectionType««2
"""
    struct IntersectionType

Describes an intersection. Each point of the returned polygon
is marked with its type w.r.t the two input simplexes.
"""
struct IntersectionType{N,P} <: AbstractVector{Tuple{P,NTuple{2,Int8}}}
# maximum of N points, constant-size
	npoints::Int
	pttype::MVector{N,NTuple{2,Int8}}
	points::MVector{N,P}
	@inline function IntersectionType{N,P}(pt::Pair...) where{N,P}
		it = new{N,P}(length(pt), MVector{N,NTuple{2,Int8}}(undef),
			MVector{N,P}(undef))
		it.pttype[1:length(pt)] .= last.(pt)
		it.points[1:length(pt)] .= first.(pt)
		return it
	end
end
@inline Base.size(it::IntersectionType) = (it.npoints,)
@inline Base.getindex(it::IntersectionType, i) = (it.points[i], it.pttype[i])
@inline function Base.show(io::IO, it::IntersectionType)
	if isempty(it) print(io, "ø"); return; end
	print(io, length(it), " points:",
		join(["\n  $a => $b" for (a,b) in it]))
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
		ij=:($(esc(:i)),$(esc(:j)))
		haskey(expr, "$a$b$c") && return (:($ij=(1,1)),expr["$a$b$c"]...)
		haskey(expr, "$b$c$a") && return (:($ij=(2,1)),expr["$b$c$a"]...)
		haskey(expr, "$c$a$b") && return (:($ij=(3,1)),expr["$c$a$b"]...)
		haskey(expr, "$d$e$f") && return (:($ij=(1,4)),expr["$d$e$f"]...)
		haskey(expr, "$e$f$d") && return (:($ij=(2,4)),expr["$e$f$d"]...)
		haskey(expr, "$f$d$e") && return (:($ij=(3,4)),expr["$f$d$e"]...)
		(:(error("case not implemented: "*$s)),)
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
macro cycle3(i,vars)#««2
	i = esc(i)
	@assert Meta.isexpr(vars, :tuple)
	n = length(vars.args)
	@assert n % 3 == 0
	v1 = esc(vars)
	v2 = Expr(:tuple,
		[esc(vars.args[i+1-3*(i%3==0)]) for i in 1:length(vars.args)]...)
	v3 = Expr(:tuple,
		[esc(vars.args[i+2-3*(i%3!=1)]) for i in 1:length(vars.args)]...)
	quote
		@assert $i ∈ (1,2,3)
		if $i == 2; $v1 = $v2
		elseif $i == 3; $v1 = $v3
		end
	end
end
# 2d intersections««1
@inline det2(u,v)=u[1]*v[2]-u[2]*v[1]
@inline det2(p,q,r)=det2(q-p,r-p)
# inter2_segment_triangle««2
"""
    Computes the intersection between 2d segment and triangle

Returns an `IntersectionType` structure, encoded in the following way:
 - if there are 2 points, they are guaranteed to be in the same order as the segments;
 - on segment, (1,2) are vertices (u,v), 0 is interior
 - on triangle, (i2,i2+1,i2+2)(mod3) are edges (qr, rp, pq),
   and likewise (+3) for vertices (p,q,r); 0 is interior
"""
function inter2_segment_triangle((u1,v1),(p2,q2,r2),i2=1;
	dpqr = missing, ε=_THICKNESS)
	(dpqr isa Missing) && (dpqr = det2(p2,q2,r2)) # FIXME, this causes 2 compilations...
	IT=IntersectionType{2,typeof(u1)}
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
	dpq_u < -ε && dpq_v < -ε && return IT()
	dqr_u < -ε && dqr_v < -ε && return IT()
	drp_u < -ε && drp_v < -ε && return IT()

	# symbolic names for returned points: ««3
	@inline P2(x)=p2=>(x,3+i2)
	@inline Q2(x)=q2=>(x,3+plus1mod3(i2))
	@inline R2(x)=r2=>(x,3+plus2mod3(i2))
	@inline U1(x)=u1=>(1,x)
	@inline V1(x)=v1=>(2,x)
	@inline QR() = i2
	@inline RP() = plus1mod3(i2)
	@inline PQ() = plus2mod3(i2)
	@inline UVQR()= barycenter(v1, u1, dqr_u/(dqr_u-dqr_v))=>(0,QR())
	@inline UVRP()= barycenter(v1, u1, drp_u/(drp_u-drp_v))=>(0,RP())
	@inline UVPQ()= barycenter(v1, u1, dpq_u/(dpq_u-dpq_v))=>(0,PQ())

	# rotate triangle to standard configuration: ««3
	@debug "position of u: $dqr_u $drp_u $dpq_u"
	@tree27((dqr_u, drp_u, dpq_u),
		"+++" => (@goto inside),
		"-++" => (@goto outside_p),
		"--+" => (@goto outside_p),
		"0++" => (@goto edge_qr),
		"00+" => (@goto vertex_p),
		"-+0" => (@goto outside_p),
		"-0+" => (@goto outside_p),
		"000" => (error("impossible case")),
		"00-" => (error("impossible case")),
		"0--" => (error("impossible case")),
		"---" => (error("impossible case")),
	)
	@assert false "not reachable"

	@label inside #««3
	@tree27((dqr_v, drp_v, dpq_v),
		"+++" => (return IT(U1(0), V1(0))),
		"000" => (error("impossible case")),
		"00-" => (error("impossible case")),
		"0--" => (error("impossible case")),
		"---" => (error("impossible case")),
		# v is a triangle vertex:
		"00+" => (return IT(U1(0), i == 1 ? P2(2) : i == 2 ? Q2(2) : R2(2))),
		# v is inside a triangle edge:
		"0++" => (return IT(U1(0), V1(i == 1 ? QR() : i == 2 ? RP() : PQ()))),
		"-0+" => (return IT(U1(0), i==1 ? UVQR()  : i == 2 ? UVRP() : UVPQ())),
		"-++" => (return IT(U1(0), i==1 ? UVQR()  : i == 2 ? UVRP() : UVPQ())),
		"--+" => begin
		@cycle3 i (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
		i2 = mod1(i2+i-1,3)
		# which side is crossed depends on ruv
		druv = det2(r2,u1,v1)
		druv > ε && return IT(U1(0), UVQR())
		druv <-ε && return IT(U1(0), UVRP())
		return IT(U1(0), R2(0))
		end
	)
	@assert false "not reachable"

	@label vertex_p #««3
	@cycle3 i (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
	@assert samepoint(u1, p2)
	i2 = mod1(i2+i-1,3)
	dpq_v < -ε && return IT(P2(1))
	drp_v < -ε && return IT(P2(1))
	if dpq_v ≤ ε
		dqr_v > ε && return IT(P2(1), V1(PQ()))
		dqr_v <-ε && return IT(P2(1), Q2(0))
		@assert samepoint(v1, q2)
		return IT(P2(1), Q2(2))
	end
	if drp_v ≤ ε
		dqr_v > ε && return IT(P2(1), V1(RP()))
		dqr_v <-ε && return IT(P2(1), R2(0))
		@assert samepoint(v1, r2)
		return IT(P2(1), R2(2))
	end
	dqr_v > ε && return IT(P2(1), V1(0))
	dqr_v <-ε && return IT(P2(1), UVQR())
	@assert collinear2(q,v,r)
	@assert monotonic(q,v,r)
	return IT(P2(1), V1(QR()))

	@label edge_qr #««3
	@cycle3 i (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
	@assert dqr_v ≥ -ε
	if dqr_v ≤ ε # u, v on same edge qr
		dpq_v <-ε && return IT(U1(PQ()), Q2(0))
		drp_v <-ε && return IT(U1(PQ()), R2(0))
		dpq_v ≤ ε && return IT(U1(PQ()), Q2(2))
		drp_v ≤ ε && return IT(U1(PQ()), R2(2))
		return IT(U1(PQ()), V1(PQ()))
	end
	@assert dqr_v > ε
	dpq_v > ε && drp_v > ε && return IT(U1(PQ()), V1(0))
	dpuv = det2(u1,v1,p2)
	dpuv > ε && return IT(U1(PQ()), UVRP())
	dpuv <-ε && return IT(U1(PQ()), UVPQ())
	@assert collinear2(u1,v1,p2)
	dpq_v <-ε && return IT(U1(PQ()), P2(0))
	@assert samepoint(v1,p2)
	return IT(U1(PQ()), P2(2))

	@label outside_p #««3
	@cycle3 i (p2,q2,r2, dqr_u,drp_u,dpq_u, dqr_v,drp_v,dpq_v)
	i2 = mod1(i2+i-1,3)

	@assert dqr_u < -ε # u1 is far from p
	@assert dpq_u > ε
	@assert dqr_v ≥ -ε
	
	if dqr_v ≤ ε # edge case: v1 on (qr) line««
		if dpq_v > ε
			drp_v > ε && return IT(V1(QR())) # ]qr[ segment
			drp_v <-ε && return IT() # outside on r-side
			@assert samepoint(v1, r2)
			return IT(R2(2))
		end
		dpq_v <-ε && return IT() # outside on q-side
		@assert samepoint(v1, q2)
		return IT(Q2(2))
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
		drp_u ≥ -ε && return IT()
		@assert drp_v ≥ -ε # because of easy reject
		# since drp_u <-ε and drp_v ≥-ε, uv crosses rp. Where does this happen?
		if drp_v > ε # generic case: v is (stricly) inside (rp).
			dpuv = drp_u - drp_v + druv
			dpuv > ε && return IT() # missed p and the triangle
			dpuv ≥-ε && return IT(P2(0)) # touched p

			@assert dpuv <-ε # we hit the ]rp[ side, and thus the triangle
			dpq_v > ε && return IT(UVRP(), V1(0)) # v is an inner point
			dpq_v <-ε && return IT(UVRP(), UVPQ()) # v is outside (pq)
			@assert collinear2(p2,v1,q2)
			@assert monotonic(p2,v1,q2)
			return IT(UVRP(), V1(PQ()))
		end
		@assert collinear2(r2,p2,v1)
		dpq_v > ε && return IT(V1(PQ()))
		dpq_v < ε && return IT()
		@assert samepoint(v1,p2)
		return IT(P2(2))
	elseif druv ≥ -ε # (uv) touches r
		@assert collinear2(u1,r2,v1)
		@assert monotonic(u1,r2,v1)
		drp_v <-ε && return IT(R2(0))
		if drp_v ≤ ε # (uv) == (rp)
			@assert collinear2(u1,v1,p2)
			dpq_v > ε && return IT(R2(0),V1(RP()))
			dpq_v <-ε && return IT(R2(0),P2(0))
			@assert samepoint(v1,p2)
			return IT(R2(0),P2(2))
		end
		@assert drp_v > ε # (uv) might yet cross (pq):
		dpq_v > ε && return IT(R2(0),V1(0))
		dpq_v < ε && return IT(R2(0),UVPQ())
		@assert collinear2(p2,q2,v1)
		return IT(R2(0),V1(PQ()))
	end
	# now check on the q-side:
	dquv = dqr_v - dqr_u + druv
	dquv < -ε && return IT() # we missed on the q-side
	if dquv ≤ ε # (uv) touches q
		@assert collinear2(u1,q2,v1)
		@assert monotonic(u1,q2,v1)
		@assert dpq_v ≤ ε
		dpq_v <-ε && return IT(Q2(0)) # far
		@assert samepoint(v1,q2)
		return IT(Q2(2))
	end
	@assert dquv > ε # uv crosses the ]q,r[ segment...
	if dpq_v > ε
		drp_v > ε && return IT(UVQR(), V1(0)) # inner point
		drp_v <-ε && return IT(UVQR(), UVRP())
		@assert collinear2(r2,v1,p2)
		@assert monotonic(r2,v1,p2)
		return IT(UVQR(), V1(RP()))
	elseif dpq_v ≥-ε
		drp_v > ε && return IT(UVQR(), V1(PQ())) # pq edge
		drp_v <-ε && return IT(UVQR(), UVRP()) # outside p on pq edge
		@assert samepoint(v1,p2)
		return IT(UVQR(), P2(2)) # 
	end
	@assert dpq_v <-ε
	# easy reject case where we don't cross rp:
	drp_v ≥-ε && return IT(UVQR(), UVPQ())
	@assert drp_v <-ε
	# worst case: we cross diagonally from (-++) to (+--) region;
	# decide which way around p we go:
	dpuv = drp_u - drp_v + druv
	dpuv > ε && return IT(UVQR(), UVRP())
	dpuv <-ε && return IT(UVQR(), UVPQ())
	@assert collinear2(u1,p2,v1)
	@assert monotonic(u1,p2,v1)
	return IT(UVQR(), P2(0))
end

# function inter(t1::Triangle{2}, t2::Triangle{2};
# 	ε=_THICKNESS, d2=det(t2), common_edge = 0)#««
# 	(p1, q1, r1) = vertices(t1)
# 	(p2, q2, r2) = vertices(t2)
# 	u2 = q2-p2; v2=r2-p2
# 
# 	# FIXME: this is already computed in inter(::Triangle{3}),
# 	# as one of the components of the `normal2` vector...
# 	# and, since we projected on the “right” coordinates, we know that area2>0
# 	@assert d2 == det(t2)
# 	@assert d2 == det(u2, v2)
# 	@assert d2 > ε
# end#»»

end # module
# 3d intersections««1
# basic geometry stuff««2


# """
#     supporting_plane(t::Triangle)
# 
# Returns an equation (`a*x = b`) of the supporting plane, with `a`
# pointing *outwards*.
# """
# function supporting_plane(t::Triangle)
# 	(p1, p2, p3) = t.vertices
# 	c = cross(p2-p1, p3-p1)
# 	b = dot(c, p1.coords)
# 	return HyperPlane(c, b)
# end

# struct HyperPlane{X,Y}
# 	a::X
# 	b::Y
# 	HyperPlane(a,b) = new{typeof(a),typeof(b)}(a,b)
# end
# direction(h::HyperPlane) = h.a

@inline function project_2d(direction::AbstractVector)
  # we inline the `findmax` call since we know the length is 3:
  # (doing this about halves the running time of this function. Besides,
  # since the value `k` only takes three possible values, it enables the
  # compiler to propagate constants.)
  @inbounds a1 = abs(direction[1]);
	@inbounds a2 = abs(direction[2]);
	@inbounds a3 = abs(direction[3])
  k = (a1 < a2) ? ((a2 < a3) ? 3 : 2) : ((a1 < a3) ? 3 : 1)
  v = direction[k]
  @inbounds p = (v > 0) ? SA[plus1mod3(k), plus2mod3(k)] :
    SA[plus2mod3(k), plus1mod3(k)]
	return (p, k)
end
# @inline Base.getindex(p::Point, i...) = getindex(coordinates(p), i...)

# # Permuting a triangle
# @inline permute3(i, (p,q,r)) =
# 	i==1 ? (p,q,r) : i==2 ? (q,r,p) : i==3 ? (r,p,q) :
# 	i==4 ? (p,r,q) : i==5 ? (q,p,r) :        (r,q,p)

# TriangleIntersection: structure containing intersection info ««2
"""
    struct TriangleIntersection

Describes the intersection of two 3d triangles, in a way which is
independent from permutations of both triangles.
Possibilities:
 - empty: indicated by `vertex[1] < 0`, tested by `isempty()`;
 - vertex-adjacent and edge-adjacent: represented as empty;
 - vertex-in-edge: `vertex` = (v, e1, e2);
 - vertex-in-face: `vertex` = (v, face, 0)
 (this is the single case depending on permutation, and only of faces);
 - any new points are represented with the info about the edge(s)
 to which they belong: (e1, e2) or (0,0) for interior points;
 - if a new edge (i.e. not a segment of an existing edge) is created,
 it is given as (e1, e2), indices of two points.
"""
struct TriangleIntersection{T} # T is a point type
	face_insert::NTuple{2,Int} # at most 1 point in 1 face
	edge_insert::Vector{NTuple{3,Int}}
	new_edge::NTuple{2,Int} # at most 1 new edge (not part of previous edge)
	newpoints::Vector{T} # added new points with coordinates
	@inline TriangleIntersection{T}(f::NTuple{2,Int},
		e::AbstractVector, ne::NTuple{2,Int}, np::AbstractVector) where{T} =
		new{T}(f, e, ne, np)
end

@inline (T::Type{<:TriangleIntersection})(::typeof(empty)) =
	T((-1,-1), [], (0,0), [])
@inline Base.isempty(t::TriangleIntersection) = t.face_insert[1] < 0
# vf case
@inline (T::Type{<:TriangleIntersection})(v::Int, f::Int) =
	T((v,f), [], (0,0), [])
# ve case
@inline (T::Type{<:TriangleIntersection})(v::Int, e::NTuple{2,Int}) =
	T((0,0), [(v,e[1],e[2])], (0,0), [])
# polygon constructor:
# start = number of first point
# edge = new edge added (if any)
# points: of the form (point-object => ((edge1), (edge2))),
# with edge2 == (0,0) if point is on a single edge
@inline _insert_edges(i, e1, e2) =
	iszero(e2[1]) ? [(i, e1...)] : [(i, e1...), (i, e2...)]
@inline function (T::Type{<:TriangleIntersection})(
		start::Int, edge::NTuple{2,Int}, points::Pair...)
	inserts = vcat((_insert_edges(i, pt[2]...)
		for (i,pt) in zip(start:start+5, points))...)
	return T((0,0), inserts, edge, [ first(p) for p in points ])
end
# which also includes the edge-edge “cross” case, and planar cases:
# @inline (T::Type{<:TriangleIntersection})(e::NTuple{2,Int},
# 	start, points::Pair...) =
# 	T(false, (0,0),
# 	# new_edges:
# 	[(i, last(p)[1]..., last(p)[2]...) for (i,p) in zip(start:start+5, point)],
# 	e, collect(first.(point)))

# SegmentTriangleInter««2
# for uv ∩ pqr = [ij] (in same order):
#  each point is either: a segment extremity 12 / interior0
#  a triangle vertex 123/edge 456/interior 0 (not exterior)
#  so 21^2 possibilities (representable on 4 octets)
#  in segment: inter = ∅, {u}, [u?], [uv], {?}, [??], [?v], {v} (8 cases)
#  i is either a vertex (123), an edge, or inner
#  j: same
#
#  if i, j == same vertex: equal
#  if i, j == same edge: 
# inter(triangle1, triangle2) ««2
"""
    inter(triangle1, triangle2; ε, indices, count)

Returns a description of the intersection of two 3-dimensional triangles.
This is returned as a IntersectionType structure, encoded in the following way:
 - 0: point is interior to one triangle
 - 1,2,3: point is on an edge qr, rp, pq
 - 4,5,6: point is on a vertex p, q, r
"""
function inter((p1,q1,r1),(p2,q2,r2), ε=_THICKNESS)
	# [Devillers, Guigue, _Faster triangle-triangle intersection tests_;
	#   https://hal.inria.fr/inria-00072100/document]
	# https://github.com/yusuketomoto/ofxCGAL/blob/master/libs/CGAL/include/CGAL/Triangle_3_Triangle_3_intersection.h
	# https://fossies.org/linux/CGAL/include/CGAL/Intersections_3/internal/Triangle_3_Triangle_3_do_intersect.h
	(p1, q1, r1) = vertices(t1)
	(p2, q2, r2) = vertices(t2)
	@debug "computing triangle-triangle intersection"

	# return types:
	IT = IntersectionType{6,typeof(p1)}

	normal2= normalize(cross(q2-p2, r2-p2))

	dp1 = dot(normal2, p1-p2)
	dq1 = dot(normal2, q1-p2)
	dr1 = dot(normal2, r1-p2)
# 	println("dp1 dq1 dr1 = $dp1 $dq1 $dr1")

	# filter for adjacent faces: these have opposed edges (detected with
	# vertex labels) and are not collinear
	# (they might also be collinear, but this is a bit harder to detect).
# 	    if idx1[1] == idx2[1] && idx1[2] == idx2[3]
# 		abs(dr1) > ε && return IT()
# 	elseif idx1[2] == idx2[2] && idx1[3] == idx2[1]
# 		abs(dp1) > ε && return IT()
# 	elseif idx1[3] == idx2[3] && idx1[1] == idx2[2]
# 		abs(dq1) > ε && return IT()
# 	end

	@inline _touch1(i1,i2) =
		inter_touch1(permute3(i1, (p1,q1,r1)), i1,
		             permute3(i2, (p2,q2,r2)), i2, 2,
								 (i2 ≤ 3) ? normal2 : -normal2, ε)
	@inline _arrow1(i1,i2) =
		inter_arrow1(permute3(i1, (p1,q1,r1)), i1,
		             permute3(i2, (p2,q2,r2)), i2,
								 permute3(i1, (dp1, dq1, dr1)),
								 normal2, count, ε)
	@inline _border1(i1,i2) =
		inter_border1(permute3(i1, (p1,q1,r1)), i1,
		              permute3(i2, (p2,q2,r2)), i2,
									count, ε)
	
	# permute both triangles as needed so that t2 separates p1 from q1, r1
	# this guarantees that line (bb2 cc2) intersects segments (a1b1) and (a1c1).
	@debug "signs for p1,q1,r1: $(sign(dp1)) $(sign(dq1)) $(sign(dr1))"
	@tree27((dp1,dq1,dr1),
		"+++" => (return IT()),
		"0--" => (return _touch1(i,j)),
		"0+-" => (return _arrow1(i,j)),
		"+00" => (return _border1(i,j)),
		"+--" => (i1=i; i2=j),
		"000" => begin 
			(coord, lift) = project_2d(plane2)
			t1proj = Triangle(p1[coord], q1[coord], r1[coord])
			t2proj = Triangle(p2[coord], q2[coord], r2[coord])
			intproj = inter(t1proj, t2proj; ε, common_edge)
			(newpoints2, newedges) = intproj
			newpoints3 = map((point,p1,p2)->(lift(point),p1,p2), newpoints2)
			# FIXME: lift polygon
			return IT()
		end,
	)
	@debug "now (i1,i2)=($i1,$i2)"

	@inline _touch2(i1,i2) =
		inter_touch1(permute3(i1, (p1,q1,r1)), i1,
		             permute3(i2, (p2,q2,r2)), i2,
		             1, (i2 ≤ 3) ? normal2 : -normal2, ε)
	# likewise for second triangle
	normal1 = cross(q1-p1, r1-p1)
	dp2 = dot(normal1, p2-p1)
	dq2 = dot(normal1, q2-p1)
	dr2 = dot(normal1, r2-p1)
	@debug "dp2 dq2 dr2: $(sign(dp2)) $(sign(dq2)) $(sign(dr2))"
	@tree27((dp2,dq2,dr2),
		"+++" => (return IT()),
		"0--" => (return _touch2(i,1)),
# 		"0+-" => (return _arrow2(i,1)),
# 		"+00" => (return _border2(i,1)),
		"+--" => (i2+=i-1),
		"-++" => (i2+=i-1; i1+=3),
	)
	@debug "now i1,i2 = $((i1,i2))"

	# apply both permutations ««
	@inbounds (a1, b1, c1) = permute3(i1, (p1,q1,r1))
	@inbounds (a2, b2, c2) = permute3(i2, (p2,q2,r2))
	# give symbolic names to what we return...
	(la1,lb1,lc1)= i1==1 ? (1,2,3) : i1==2 ? (2,3,1) : (3,1,2)
	(la2,lb2,lc2)= i2==1 ? (1,2,3) : i2==2 ? (2,3,1) : (3,1,2)
	# re-use already computed determinants as z-coordinates:
	@inbounds (za1, zb1, zc1) = permute3(i1, (dp1, dq1, dr1))
	@inbounds (za2, zb2, zc2) = permute3(i2, (dp2, dq2, dr2))
	# 1,2,3 represent even permutations, and 4,5,6 odd ones:
	(i1 > 3) && (normal1 = -normal1; za2 = -za2; zb2 = -zb2; zc2 = -zc2)
	(i2 > 3) && (normal2 = -normal2; za1 = -za1; zb1 = -zb1; zc1 = -zc1)

	# a laundry list of assertions to check that we are in a standard
	# configuration:
	@assert normal1 ≈ cross(b1-a1, c1-a1)
	@assert normal2 ≈ cross(b2-a2, c2-a2)
	@assert dot(normal2, a1-a2) ≈ za1
	@assert dot(normal2, b1-a2) ≈ zb1
	@assert dot(normal2, c1-a2) ≈ zc1
	@assert dot(normal1, a2-a1) ≈ za2
	@assert dot(normal1, b2-a1) ≈ zb2
	@assert dot(normal1, c2-a1) ≈ zc2
	@assert za1 ≥ 0; @assert zb1 ≤ 0; @assert zc1 ≤ 0
	@assert za2 ≥ 0; @assert zb2 ≤ 0; @assert zc2 ≤ 0
	# »»
	# two cases where we know that the intersection is empty:
	# db1b2 is a determinant comparing the positions of b1 and b2
	a1a2b1 = cross(a2-a1, b1-a1)
	db1b2 = dot(a1a2b1, b2-a1)
	@debug "db1b2=$db1b2"
	db1b2 < -ε && return IT()

	a1a2c1 = cross(a2-a1, c1-a1)
	dc1c2 = dot(a1a2c1, c2-a1)
	@debug "dc1c2=$dc1c2"
	dc1c2 > ε && return IT()

	# coordinates of four intersection points bb1, cc1, bb2, cc2
	# (all four are aligned on the intersection of the two planes)
	@inline bb1() = barycenter(b1, a1, za1/(za1-zb1))
	@inline cc1() = barycenter(c1, a1, za1/(za1-zc1))
	@inline bb2() = barycenter(b2, a2, za2/(za2-zb2))
	@inline cc2() = barycenter(c2, a2, za2/(za2-zc2))

	@inline segment(p1,p2) = TI(count+1, (count+1, count+2), p1, p2)
	@inline B1(t1,t2) = bb1() => ((la1, lb1), (t1,t2))
	@inline B2(t1,t2) = bb2() => ((la2, lb2), (t1,t2))
	@inline C1(t1,t2) = cc1() => ((la1, lc1), (t1,t2))
	@inline C2(t1,t2) = cc2() => ((la2, lc2), (t1,t2))

	# TODO: +00 cases (3a or 3b in https://www.polibits.gelbukh.com/2013_48/Triangle-Triangle%20Intersection%20Determination%20and%20Classification%20to%20Support%20Qualitative%20Spatial%20Reasoning.pdf)
	db1c2 = dot(a1a2b1, c2-a1)
	@debug "db1c2 = $db1c2"
	if db1c2 < -ε
		db1b2 ≤ ε && return single(B2(la1,lb1))
		dc1b2 = dot(a1a2c1, b2-a1)
		dc1b2 <-ε && return segment(B1(0,0), B2(0,0))
		dc1b2 ≤ ε && return segment(B1(0,0), B2(la1,lc1))
		             return segment(B1(0,0), C1(0,0))
	elseif db1c2 ≤ ε # cc2 ∈ edge [a1,b1]
		dc1b2 = dot(a1a2c1, b2-a1)
		dc1b2 <-ε && return segment(C2(la1,lb1), B2(0,0))
		dc1b2 ≤ ε && return segment(C2(la1,lb1), B2(la1,lc1))
		             return segment(C2(la1,lb1), C1(0,0))
	elseif dc1c2 <-ε
		dc1b2 = dot(a1a2c1, b2-a1)
		dc1b2 <-ε && return segment(C2(0,0), B2(0,0))
		dc1b2 ≤ ε && return segment(C2(0,0), B2(la1,lc1))
		             return segment(C2(0,0), C1(0,0))
	end
	@assert dc1c2 ≤ ε # empty case was already detected above
	return single(C2(a1c1)) # cc2 ∈ edge [a1,c1]
end
# inter_touch1««2
"""
    inter_touch1

Returns intersection in the case where triangle `t1` touches plane of `t2`
at the single point `p1` (oriented so that q1, r1 are *below* the plane).
"""
function inter_touch1((p1,q1,r1), i1, (p2,q2,r2), i2, face, normal2, ε)
	# Possible cases here:
	# - nothing
	# - vertex-vertex (return empty)
	# - vertex1 on edge2
	# - vertex1 inside face2
	# To decide: project everything on the plane of face2
	# and compute three 2d determinants
	(coord, e) = project_2d(normal2, Val(true))
	@debug "entering inter_touch1 $i1, $i2"
# 	println("point $p1 touches other plane:")
# 	println("tri2 = ($p2,$q2,$r2) $idx2")
# 	println("tri1 = ($p1,$q1,$r1) $idx1")
	IT=IntersectionType{6,typeof(p1)}

	a1 = p1[coord]
	a2 = p2[coord]
	b2 = q2[coord]
	c2 = r2[coord]
	dabc = abs(normal2[e])
	@assert cross(b2-a2,c2-a2) ≈ dabc
	dab = cross(b2-a2, a1-a2)
	dab < -ε && return IT()
	dbc = cross(c2-b2, a1-b2)
	dbc < -ε && return IT()
	dca = dabc - dab - dbc
# 	println("proj = ($coord, $e) from $normal2")
# 	println("projections: $p1 => $a1 vs. ($a2, $b2, $c2)")
# 	println("determinants: $dabc =  $dab, $dbc, $dca")
	dca < -ε && return IT()
	@assert abs(dca-cross(a2-c2,a1-c2)) < 1e-9
	@debug "point $a1 in ($a2,$b2,$c2): $dab, $dbc, $dca"
	if dab ≤ ε # it cannot be <-ε, so a1 ∈ [a2,b2]
		if dbc ≤ ε # a1 == b2
			@assert norm(a1-b2) ≤ ε
			return IT(a1 => (i1+3, plus1mod3(i2, 3)))
		elseif dca ≤ ε
			@assert norm(a1-a2) ≤ ε
			return IT(a1 => (i1+3, i2+3))
		# a1 ∈ open segment ]a2, b2[
		return IT(a1 => (i1+3, plus2mod3(i2)))
	elseif dbc ≤ ε # a1 ∈ ]b2,c2] (b2 was already treated above)
		if ca ≤ ε # a1 == c2
			@assert norm(a1-c2) ≤ ε
			return IT(a1 => (i1+3, plus2mod3(i2, 3)))
		end # open segment ]b2,c2[
		return IT(a1 => (i1+3, i2))
	elseif dca ≤ ε
		return IT(a1 => (i1+3, plus1mod3(i2)))
	end
	# a1 ∈ interior of other face
	return IT(a1 => (i1+3, 0))
end
# inter_arrow1 ««2
"""
    inter_arrow1

Computes intersection, assuming that the plane of one triangle cuts the
other through vertex `p1` (with `q1` on the + side and `r1` on the - side).

The difference with the generic case is that this has a symmetry:
(q1↔r1, q2↔r2) leaves '0+-' unchanged.
"""
function inter_arrow1(((p1,q1,r1),idx1), ((p2,q2,r2),idx2),
	(zp1, zq1, zr1), normal2, count, ε)
	(coord, e) = project_2d(normal2, Val(true))
	println("point $p1 crosses other plane (normal=$normal2, $coord, $e):")
	println("tri2 = ($p2,$q2,$r2) $idx2")
	println("tri1 = ($p1,$q1,$r1) $idx1")
	TI = TriangleIntersection{typeof(p1)}

	@assert abs(zp1) ≤ ε
	@assert zq1 > ε
	@assert zr1 <-ε

	p1b = p1[coord]; q1b = q1[coord]; r1b = r1[coord]
	p2b = p2[coord]; q2b = q2[coord]; r2b = r2[coord]
	# we know that interior of segment (q1,r1) intersects plane 2:
	u = barycenter(q1b, r1b, zr1/(zr1-zq1))
	println("projection on plane 2: ($p1b, $q1b, $r1b), ($p2b, $q2b, $r2b)")
	println("other segment is $u")
	dpqr = abs(normal2[e])
	@assert cross(q2b-p2b,r2b-p2b) ≈ dpqr
	# the position of u w.r.t the triangle is given by these three determinants:
	dpq_p = cross(q2b-p2b, p1b-p2b)
	dqr_p = cross(r2b-q2b, p1b-q2b)
	drp_p = dpqr - dpq_p - dqr_p
	@assert abs(drp_p-cross(p2b-r2b,p1b-r2b)) < 1e-9
	# likewise for p1b:
	dpq_u = cross(q2b-p2b, u-p2b)
	dqr_u = cross(r2b-q2b, u-q2b)
	drp_u = dpqr - dpq_u - dqr_u
	@assert abs(drp_u-cross(p2b-r2b,u-r2b)) < 1e-9
	# Sutherland-Cohen to intersect segment [a,u] with triangle a2b2c2:
	# trivial reject cases
	dpq_p < -ε && dpq_u < -ε && return IT()
	dqr_p < -ε && dqr_u < -ε && return IT()
	drp_p < -ε && drp_u < -ε && return IT()
	
	
end
