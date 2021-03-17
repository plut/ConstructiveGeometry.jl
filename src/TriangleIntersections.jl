module TriangleIntersections
import Meshes: Point, Triangle, coordtype, vertices, coordinates
using LinearAlgebra
using StaticArrays

# predicates for assertions:
@inline barycenter(p1, p2, λ) = p2 + λ*(p1-p2) # λp1+(1-λ)p2
@inline collinear(a,b,c) = norm(cross(a-b,c-b),1) ≤ ε
@inline monotonic(a,b,c) = dot(a-b,c-b) ≤ 0
@inline samepoint(a,b) = norm(a-b,1) ≤ ε

@inline real_type(x...)=Float64
@inline det(t::Triangle{2}) =
	cross(vertices(t)[2]-vertices(t)[1], vertices(t)[3]-vertices(t)[1])

macro tree27(vars,args...)#««
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
	println(build(1,expr2,""))
	build(1,expr2,"")
# 	quote#««
# 		if $(vars[1]) > $ε
# 			if $(vars[2]) > $ε
# 				if     $(vars[3]) > $ε; $(exprfor("+++")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("++-")...)
# 				else                  ; $(exprfor("++0")...)
# 				end
# 			elseif $(vars[2]) <-$ε
# 				if     $(vars[3]) > $ε; $(exprfor("+-+")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("+--")...)
# 				else                  ; $(exprfor("+-0")...)
# 				end
# 			else
# 				if     $(vars[3]) > $ε; $(exprfor("+0+")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("+0-")...)
# 				else                  ; $(exprfor("+00")...)
# 				end
# 			end
# 		elseif $(vars[1]) <-$ε
# 			if $(vars[2]) > $ε
# 				if     $(vars[3]) > $ε; $(exprfor("-++")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("-+-")...)
# 				else                  ; $(exprfor("-+0")...)
# 				end
# 			elseif $(vars[2]) <-$ε
# 				if     $(vars[3]) > $ε; $(exprfor("--+")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("---")...)
# 				else                  ; $(exprfor("--0")...)
# 				end
# 			else
# 				if     $(vars[3]) > $ε; $(exprfor("-0+")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("-0-")...)
# 				else                  ; $(exprfor("-00")...)
# 				end
# 			end
# 		else
# 			if $(vars[2]) > $ε
# 				if     $(vars[3]) > $ε; $(exprfor("0++")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("0+-")...)
# 				else                  ; $(exprfor("0+0")...)
# 				end
# 			elseif $(vars[2]) <-$ε
# 				if     $(vars[3]) > $ε; $(exprfor("0-+")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("0--")...)
# 				else                  ; $(exprfor("0-0")...)
# 				end
# 			else
# 				if     $(vars[3]) > $ε; $(exprfor("00+")...)
# 				elseif $(vars[3]) <-$ε; $(exprfor("00-")...)
# 				else                  ; $(exprfor("000")...)
# 				end
# 			end
# 		end
# 	end#»»
end#»»

_THICKNESS=1e-8

"""
    supporting_plane(t::Triangle)

Returns an equation (`a*x = b`) of the supporting plane, with `a`
pointing *outwards*.
"""

function supporting_plane(t::Triangle)
	(p1, p2, p3) = t.vertices
	c = cross(p2-p1, p3-p1)
	b = dot(c, p1.coords)
	return HyperPlane(c, b)
end

struct HyperPlane{X,Y}
	a::X
	b::Y
	HyperPlane(a,b) = new{typeof(a),typeof(b)}(a,b)
end
direction(h::HyperPlane) = h.a

const plus1mod3 = SA[2,3,1]
const plus2mod3 = SA[3,1,2]

@inline function project_2d(direction::AbstractVector, index::Val = Val(false))
  # we inline the `findmax` call since we know the length is 3:
  # (doing this about halves the running time of this function. Besides,
  # since the value `k` only takes three possible values, it enables the
  # compiler to propagate constants.)
  @inbounds a1 = abs(direction[1]);
	@inbounds a2 = abs(direction[2]);
	@inbounds a3 = abs(direction[3])
  k = (a1 < a2) ? ((a2 < a3) ? 3 : 2) : ((a1 < a3) ? 3 : 1)
  v = direction[k]
  @inbounds p = (v > 0) ? SA[plus1mod3[k], plus2mod3[k]] :
    SA[plus2mod3[k], plus1mod3[k]]
	return _project_2d(index, p, k)
end
@inline _project_2d(::Val{false}, p, _) = p
@inline _project_2d(::Val{true}, p, e) = (p, e)
@inline Base.getindex(p::Point, i...) = getindex(coordinates(p), i...)


#=
function _build_triangle3_decision_tree(
	changes::Pair{<:AbstractString,<:AbstractString}...)
	@inline to_idx(s::AbstractString) =
		[ c == '+' ? 1 : c == '-' ? 2 : 0 for c in s ]
	@inline function apply_perm3(i, j, s)
		tmp=s[[i+1, mod1(i+2,3), mod1(i+3, 3)]]
		(j == 0) ? tmp : map(x->mod(3-x,3), tmp)
	end
	@inline function find_perm3(s1, s2)
		for i in 0:2, j in (0,3)
			apply_perm3(i, j, s1) == s2 && return (i, j)
		end
		error("No permutation found transforming \"$s1\" to \"$s2\"")
	end
	tree = fill((0,0), 27)
	for (s1, s2) in changes
		(u1, u2) = to_idx(s1), to_idx(s2)
		@assert s2 ∈ ("+--", "+0-", "+00", "0--", "---", "+++", "000") "bad output: $s2"
		i = mod1( sum((9,3,1) .* u1), 27)
		(j1, j2) = find_perm3(u1, u2)
		@assert (s2[1] == '0') == (i ∈ (27, 4,8,10,12,20,24)) "$s2: $i"
		println("[$i] $s1 ⇒ $s2 by ($j1, $j2)")
		tree[i] = (j1, j2)
	end
	return tree
end
const _triangle3_decision_tree = _build_triangle3_decision_tree(
	"000" => "000", "00+" => "+00", "00-" => "+00", 
	"0+0" => "+00", "0++" => "0--", "0+-" => "+0-", 
	"0-0" => "+00", "0-+" => "+0-", "0--" => "0--", 
	"+00" => "+00", "+0+" => "0--", "+0-" => "+0-", 
	"++0" => "0--", "+++" => "+++", "++-" => "+--", 
	"+-0" => "+0-", "+-+" => "+--", "+--" => "+--", 
	"-00" => "+00", "-0+" => "+0-", "-0-" => "0--", 
	"-+0" => "+0-", "-++" => "+--", "-+-" => "+--", 
	"--0" => "0--", "--+" => "+--", "---" => "---", 
)
==#

# Permuting a triangle
@inline permute3(i, (p,q,r)) =
	i==1 ? (p,q,r) : i==2 ? (q,r,p) : i==3 ? (r,p,q) :
	i==4 ? (p,r,q) : i==5 ? (q,p,r) :        (r,q,p)

# TriangleIntersection: structure containing intersection info ««
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
# end
# # 	vertex::NTuple{3,Int}
# 	# up to 6 new points:
# 	# (new point, edge1 or 0, edge2 or 0)
# 	newpoints::Vector{Tuple{T,Int,Int,Int,Int}}
# 	# cross case: one new edge max
# 	# coplanar case: new edges are all implicit (they are segments of old
# 	# edges); no need to fill this structure (hence (0,0))
# 	edge::NTuple{2,Int}
# 	@inline TriangleIntersection(v::NTuple{3,Int},
# 		p,e) = new{coordtype(eltype(newpoints).parameters[1])}(v,p,e)
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
#»»

# inter(triangle1, triangle2) ««
"""
    inter(triangle1, triangle2; ε, indices, count)

Returns a description of the intersection of two 3-dimensional triangles.
`indices` is the list of names of vertices of the two triangles
(default `(1,2,3,4,5,6)`) and `count` the total number of points
(default `6`).
 * `face2`: same.

FIXME: missing
 - edge-inside-edge (Fig 3a: edge_edgeCollinear)
"""
function inter(t1::Triangle{3}, t2::Triangle{3};
	idx1=(1,2,3), idx2=(4,5,6), count=6,
	ε=_THICKNESS)
	# [Devillers, Guigue, _Faster triangle-triangle intersection tests_;
	#   https://hal.inria.fr/inria-00072100/document]
	# https://github.com/yusuketomoto/ofxCGAL/blob/master/libs/CGAL/include/CGAL/Triangle_3_Triangle_3_intersection.h
	# https://fossies.org/linux/CGAL/include/CGAL/Intersections_3/internal/Triangle_3_Triangle_3_do_intersect.h
	(p1, q1, r1) = vertices(t1)
	(p2, q2, r2) = vertices(t2)

	# return types:
	TI = TriangleIntersection{typeof(p1)}

	plane2 = (supporting_plane(t2))
	normal2= direction(plane2)

	dp1 = dot(normal2, p1-p2)
	dq1 = dot(normal2, q1-p2)
	dr1 = dot(normal2, r1-p2)
# 	println("dp1 dq1 dr1 = $dp1 $dq1 $dr1")

	# filter for adjacent faces: these have opposed edges (detected with
	# vertex labels) and are not collinear
	# (they might also be collinear, but this is a bit harder to detect).
	    if idx1[1] == idx2[1] && idx1[2] == idx2[3]
		abs(dr1) > ε && return TI(empty)
	elseif idx1[2] == idx2[2] && idx1[3] == idx2[1]
		abs(dp1) > ε && return TI(empty)
	elseif idx1[3] == idx2[3] && idx1[1] == idx2[2]
		abs(dq1) > ε && return TI(empty)
	end

	@inline _touch1(i1,i2) =
		inter_touch1(permute3.(i1, ((p1,q1,r1), idx1)),
		             permute3.(i2, ((p2,q2,r2), idx2)),
		             2, (i2 ≤ 3) ? normal2 : -normal2, ε)
	@inline _arrow1(i1,i2) =
		inter_arrow1(permute3.(i1, ((p1,q1,r1), idx1)),
		             permute3.(i2, ((p2,q2,r2), idx2)),
								 permute3(i1, (dp1, dq1, dr1)),
								 normal2, count, ε)
	@inline _border1(i1,i2) =
		inter_border1(permute3.(i1, ((p1,q1,r1), idx1)),
		              permute3.(i2, ((p2,q2,r2), idx2)),
									count, ε)
	
	# permute both triangles as needed so that t2 separates p1 from q1, r1
	# this guarantees that line (bb2 cc2) intersects segments (a1b1) and (a1c1).
	@tree27((dp1,dq1,dr1),
		"+++" => (return TI(empty)),
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
			return TI(empty)
		end,
	)
# 	println("dp1 dq1 dr1: $(sign(dp1)) $(sign(dq1)) $(sign(dr1)) $((i1,i2))")

	@inline _touch2(i1,i2) =
		inter_touch1(permute3.(i1, ((p1,q1,r1), idx1)),
		             permute3.(i2, ((p2,q2,r2), idx2)),
		             1, (i2 ≤ 3) ? normal2 : -normal2, ε)
	# likewise for second triangle
	normal1 = cross(q1-p1, r1-p1)
	dp2 = dot(normal1, p2-p1)
	dq2 = dot(normal1, q2-p1)
	dr2 = dot(normal1, r2-p1)
	@tree27((dp2,dq2,dr2),
		"+++" => (return TI(empty)),
		"0--" => (return _touch2(i,1)),
# 		"0+-" => (return _arrow2(i,1)),
# 		"+00" => (return _border2(i,1)),
		"+--" => (i2+=i-1),
		"-++" => (i2+=i-1; i1+=3),
	)
# 	println("dp2 dq2 dr2: $(sign(dp2)) $(sign(dq2)) $(sign(dr2)) $((i2,i2))")

	# apply both permutations ««
	@inbounds (a1, b1, c1) = permute3(i1, (p1,q1,r1))
	@inbounds (a2, b2, c2) = permute3(i2, (p2,q2,r2))
	# give symbolic names to what we return...
	@inbounds (la1,lb1,lc1)= permute3(i1, idx1)
	@inbounds (la2,lb2,lc2)= permute3(i2, idx2)
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
	println("  db1b2=$db1b2")
	db1b2 < -ε && return TI(empty)

	a1a2c1 = cross(a2-a1, c1-a1)
	dc1c2 = dot(a1a2c1, c2-a1)
	println("  dc1c2=$dc1c2")
	dc1c2 > ε && return TI(empty)

	# coordinates of four intersection points bb1, cc1, bb2, cc2
	# (all four are aligned on the intersection of the two planes)
	@inline bb1() = barycenter(b1, a1, za1/(za1-zb1))
	@inline cc1() = barycenter(c1, a1, za1/(za1-zc1))
	@inline bb2() = barycenter(b2, a2, za2/(za2-zb2))
	@inline cc2() = barycenter(c2, a2, za2/(za2-zc2))

	@inline single(p) = TI(count+1, (0,0), p)
	@inline segment(p1,p2) = TI(count+1, (count+1, count+2), p1, p2)
	@inline B1(t1,t2) = bb1() => ((la1, lb1), (t1,t2))
	@inline B2(t1,t2) = bb2() => ((la2, lb2), (t1,t2))
	@inline C1(t1,t2) = cc1() => ((la1, lc1), (t1,t2))
	@inline C2(t1,t2) = cc2() => ((la2, lc2), (t1,t2))

	# TODO: +00 cases (3a or 3b in https://www.polibits.gelbukh.com/2013_48/Triangle-Triangle%20Intersection%20Determination%20and%20Classification%20to%20Support%20Qualitative%20Spatial%20Reasoning.pdf)
	db1c2 = dot(a1a2b1, c2-a1)
	println("db1c2 = $db1c2")
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
end#»»

# inter_touch1««
"""
    inter_touch1

Returns intersection in the case where triangle `t1` touches plane of `t2`
at the single point `p1` (oriented so that q1, r1 are *below* the plane).
"""
function inter_touch1(((p1,q1,r1),idx1), ((p2,q2,r2),idx2), face, normal2, ε)
	# Possible cases here:
	# - nothing
	# - vertex-vertex (return empty)
	# - vertex1 on edge2
	# - vertex1 inside face2
	# To decide: project everything on the plane of face2
	# and compute three 2d determinants
	(coord, e) = project_2d(normal2, Val(true))
# 	println("point $p1 touches other plane:")
# 	println("tri2 = ($p2,$q2,$r2) $idx2")
# 	println("tri1 = ($p1,$q1,$r1) $idx1")
	TI = TriangleIntersection{typeof(p1)}

	a1 = p1[coord]
	a2 = p2[coord]
	b2 = q2[coord]
	c2 = r2[coord]
	dabc = abs(normal2[e])
	@assert cross(b2-a2,c2-a2) ≈ dabc
	dab = cross(b2-a2, a1-a2)
	dab < -ε && return TI(empty)
	dbc = cross(c2-b2, a1-b2)
	dbc < -ε && return TI(empty)
	dca = dabc - dab - dbc
# 	println("proj = ($coord, $e) from $normal2")
# 	println("projections: $p1 => $a1 vs. ($a2, $b2, $c2)")
# 	println("determinants: $dabc =  $dab, $dbc, $dca")
	dca < -ε && return TI(empty)
	@assert abs(dca-cross(a2-c2,a1-c2)) < 1e-9
# 	println("point $a1 in ($a2,$b2,$c2): $dab, $dbc, $dca")
	if dab ≤ ε
		if dbc ≤ ε
			@assert norm(a1-b2) ≤ ε
			return TI(empty)
		elseif dca ≤ ε
			@assert norm(a1-a2) ≤ ε
			return TI(empty)
		end
		# a1 ∈ open segment ]a2, b2[
		return (idx1[1], idx2[1], idx2[2])
	elseif dbc ≤ ε
		if ca ≤ ε
			@assert norm(a1-c2) ≤ ε
			return TI(empty)
		end
		return TI(idx1[1], idx2[2], idx2[3])
	elseif dca ≤ ε
		return TI(idx1[1], idx2[3], idx2[1])
	end
	# a1 ∈ interior of other face
	return TI(idx1[1], 2, 0)
end#»»
# inter_arrow1 ««
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
	dpq_p < -ε && dpq_u < -ε && return TI(empty)
	dqr_p < -ε && dqr_u < -ε && return TI(empty)
	drp_p < -ε && drp_u < -ε && return TI(empty)
	# Rotate p2q2r2 so that p1 is either inside, or outside (pq) but inside
	# (qr), i.e. change to signs -+? if possible:
# 	if dpq_p ≥ ε
# 		if dqr_p ≥ ε
# 			if     drp_p ≥ ε i=1 # inside the triangle
# 			elseif drp_p ≤-ε i=2 
	
	
end #»»

"""
    Computes the intersection between 2d segment and triangle

Returns a `TriangleIntersection` structure.
"""
function inter2_segment_triangle((u1,v1),idx1,(p2,q2,r2),idx2;
	dpqr = cross(q2-p2,r2-p2), ε)
	# compute position of u:
	dpq_u = cross(q2-p2, u-p2)
	dqr_u = cross(r2-q2, u-q2)
	drp_u = dpqr - dpq_u - dqr_u
	@assert abs(drp_u-cross(p2-r2,u-r2)) < 1e-9
	# likewise for v:
	dpq_v = cross(q2-p2, v-p2)
	dqr_v = cross(r2-q2, v-q2)
	drp_v = dpqr - dpq_v - dqr_v
	# easy reject cases: segment lies fully outside one side of triangle
	dpq_u < -ε && dpq_v < -ε && return TI(empty)
	dqr_u < -ε && dqr_v < -ε && return TI(empty)
	drp_u < -ε && drp_v < -ε && return TI(empty)
	# rotate triangle to standard configuration:
	# -+?, i.e. out
	@inline _out1(i) = inter2_segment_triangle_out1( # -++
		(u1,v1),idx1, permute3.(i,(p2,q2,r2),idx2)...; dpqr, ε)
	@inline _out2(i) = inter2_segment_triangle_out2( # --+
		(u1,v1),idx1, permute3.(i,(p2,q2,r2),idx2)...; dpqr, ε)
	@inline _edge(i) = inter2_segment_triangle_edge( # 0++
		(u1,v1),idx1, permute3.(i,(p2,q2,r2),idx2)...; dpqr, ε)
	@inline _vertex(i) = inter2_segment_triangle_vertex(
		(u1,v1),idx1, permute3.(i,(p2,q2,r2),idx2)...; dpqr, ε)
# 	if dqr_u > ε
# 		if drp_u > ε
# 			if     dpq_u > ε return _inside()
# 			elseif dpq_u <-ε return _out1(3)  # ++- => -++
# 			else             return _edge(3)  # ++0 => 0++
# 			end
# 		elseif drp_u < -ε # +-
# 			if     dpq_u > ε return _out1(2)  # +-+ => -++
# 			elseif dpq_u <-ε return _out2(2)  # +-- => --+
# 			else             @assert false    # +-0 ?
# 			end
# 		else # +0
# 			if     dpq_u > ε return _edge(2)  # +0+ => 0++
# 			elseif dpq_u <-ε @assert false # +0- ?
# 			else             return _vertex(1)# +00
# 			end
# 		end
# 	elseif  dqr_u <-ε
# 		if drp_u > ε
# 			if     dpq_u > ε return _out1(2)  # +-+ => -++
# 			elseif dpq_u <-ε return _out2(2)  # +-- => --+
# 			else             @assert false    # -+0 ?
# 			end
# 		elseif drp_u <-ε
# 			if 
end

function inter(t1::Triangle{2}, t2::Triangle{2};
	ε=_THICKNESS, d2=det(t2), common_edge = 0)#««
	(p1, q1, r1) = vertices(t1)
	(p2, q2, r2) = vertices(t2)
	u2 = q2-p2; v2=r2-p2

	# FIXME: this is already computed in inter(::Triangle{3}),
	# as one of the components of the `normal2` vector...
	# and, since we projected on the “right” coordinates, we know that area2>0
	@assert d2 == det(t2)
	@assert d2 == det(u2, v2)
	@assert d2 > ε
end#»»

end # module
