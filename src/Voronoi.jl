# With some ideas taken from VRONI as described by [Huber 2008]:
# https://www.sciencedirect.com/science/article/pii/S0925772101000037
# https://www.sthu.org/research/publications/files/mscthesis.pdf
#
# Offset algorithm: [Kim 1998]
# https://www.sciencedirect.com/science/article/abs/pii/S0010448598000633
#
# FIXME:
#
# TODO:
# - remove segdistance²
#
"""    Voronoi

Computation of Voronoi diagrams for planar polygons and polygonal paths,
and of offset paths using these diagrams."""
module Voronoi
using StaticArrays
using FastClosures
using LinearAlgebra
using Random
using Printf
using HypergeometricFunctions
module LibTriangle
	using Triangle
end

include("CornerTables.jl")
using .CornerTables
import .CornerTables: triangulation

const DEFAULT_ATOL=1e-1

# CONVENTIONS FOR EDGES ««1
# as a triangulation: next is (head->left), prev is (left->tail)
#         left
# ╲         ╱╲        ╱      e: an edge in the graph
#  ╲      ↗╱╱↖╲      ╱       n: next(e); p: prev(e)
#   ╲    a╱p  n╲    ╱        a: after(e) = opposite(p)
#    ╲  ╱╱↙ ●  ╲╲  ╱         b: before(e) = next(o)
#     ╲ ╱  —e→   ╲╱          o: opposite(e)
#    tail ——————head         ●: node(e)
#       ╲  ←o—  ╱
#        ╲b    ╱
#
#  as a Voronoi diagram:
#       ╲n  head    ╱ 
#  left  ╲____o____╱  right
#        ╱●   e    ╲
#      p╱a  tail   b╲

# Geometry ««1
# Elementary geometry ««2

abstract type GeometryException <: Exception end
struct CrossingSegments <: GeometryException end
struct PointInSegment <: GeometryException end
struct NotImplemented <: GeometryException end
struct ConcurrentLines <: GeometryException end

@inline det2(u,v) = u[1]*v[2]-u[2]*v[1]
@inline det2(u,v,w) = det2(v-u, w-u)
@inline norm²(v) = v[1]^2+v[2]^2
@inline unit(v) = v/√(norm²(v))
@inline distance²(a,b) = norm²(a-b)
@inline quarterturn(v) = SA[-v[2], v[1]]
@inline sqrtp(x) = √(max(x,zero(x)))

const Segment{T} = NTuple{2,<:AbstractVector{T}}

function lineinter(a,b,c,d)
	D = det2(a-b, c-d)
	t = det2(a-c, c-d)
	z = a+(t/D)*(b-a)
	return z
end

function segdistance²(a,b,c)
	ab,ac = b-a,c-a
	d = dot(ab,ac); ab2 = norm²(ab)
	d < 0 && return norm²(ac)
	d > ab2 && return norm²(c-b)
	return det2(ab,ac)^2/norm²(ab)
end

"    iscloser(a,b,c): is d(a,b) ≤ d(a,c) ?"
@inline iscloser(a,b,c) = dot(2a-b-c, b-c) ≥ 0
@inline isleft(a,b,c) = det2(a,b,c) > 0
@inline isleft(u,v) = det2(u,v) > 0

"""
    isincircle(a,b,c,x)

Returns `true` iff point x is in circumcircle of oriented triangle `(a,b,c)`.
"""
function isincircle(a,b,c,x)
	a,b,c = a-x,b-x,c-x
	@assert !isleft(a,c,b) "incircle: triangle ($a,$b,$c) has wrong orientation"
	m = SA[a[1] a[2] norm²(a); b[1] b[2] norm²(b); c[1] c[2] norm²(c)]
	return det(m) > 0
end

# Approximation of parabolic arc««2
H(x)=x*_₂F₁(1/4,1/2,3/2,-x^2)
H′(x)=_₂F₁(1/4,1/2,3/2,-x^2)-1/6*x^2*_₂F₁(5/4,3/2,5/2,-x^2)

function Hinv(y)
	iszero(y) && return y
	# use a rough first approximation to initialize Newton's method
	x = abs(y) < 3 ? y*(1.050302+.046546*y^2) : sign(y)*(y+1.1981402347355918)^2/4
	for _ in 1:4
		x = x - (H(x)-y)/H′(x)
	end
	return x
end

"""    approxparabola(a, x1, x2, δ)
Approximates the parabola y = a/2 + x²/2a by a polygonal chain
with maximal Hausdorff distance δ on the interval [x1, x2];
returns a vector of abscissas."""
function approxparabola(a::Real,x1::Real,x2::Real, δ)
	s1, s2 = H(x1/a), H(x2/a)
	n = ceil(Int,abs(s2-s1)*√(a/8δ))
	v = sizehint!([float(x1)],n+1)
	for i in 1:n-1
		push!(v, a*Hinv(s1+(i/n)*(s2-s1)))
	end
	push!(v, x2)
	return v
end
# x ↔ r reparametrization:
# a ↔ rmin
# x = norm(b.tangent) * √(r-rmin)
# r = rmin + (x/norm(b.tangent))^2

# Triangulation via libtriangle««2
function triangulate_loop(points, idx)
	n = length(idx)
	for i in idx
		println("$(points[i][1])\t$(points[i][2])\t$i")
	end
	vmat = Float64[ points[i][j] for i in idx, j in 1:2 ]
	elist = Int[ idx[mod1(i+j,n)] for j in 1:n, i in 0:1]
	println("vmat=",vmat)
	println("idx=",idx)
	println("elist=",elist)
	return LibTriangle.constrained_triangulation(vmat, idx, elist)
end

# split_loop ««2
""" split_loop(points)

Splits a loop in the (x,y) plane at the (x=0) axis;
returns matched lists (list of new loops), (list of sides)
(with sign encoded as `true` for left side and `false` for right side)."""
function split_loop(loop)
	loop = copy(loop) # we will modify this array
	j, q = length(loop), last(loop)
	newloops = (typeof(loop)[], typeof(loop)[])
	loopstart = 0
	z = Int[] # indices of zero-crossings
	# we cannot iterate on pairs() since we grow the array
	i = 1; while i ≤ length(loop); p = loop[i]
		if (p[1] < 0 < q[1]) || (q[1] < 0 < p[1])
			p = SA[zero(p[1]), (p[1]*q[2]-p[2]*q[1])/(p[1]-q[1])]
			insert!(loop, i, p)
		end
		iszero(p[1]) && push!(z, i)
		j, q = i, p
		i+=1
	end
	sort!(z; by=i->loop[i][2])
	n = length(loop)
	b0, b1 = false, false # is left/right side of line in polygon?
	next0, next1 = Dict{Int,Int}(), Dict{Int,Int}()
	last0, last1 = 0, 0
	nextindex = @closure i->(i==n) ? 1 : i+1
	previndex = @closure i->(i==1) ? n : i-1
	# build correspondence tables on both sides
	for i2 in z#««
		s1, s3 = (Int(sign(loop[i][1])) for i in (previndex(i2), nextindex(i2)))
		if (s1 < 0) ⊻ (s3 < 0)
			if b0
				next0[last0] = i2; b0 = false
			else
				last0 = i2; b0 = true
			end
		end
		if (s1 > 0) ⊻ (s3 > 0)
			if b1
				next1[i2] = last1; b1 = false
			else
				last1 = i2; b1 = true
			end
		end
# 		(s1 < 0) ⊻ (s3 < 0) && (b0 = !b0; push!(c0, i2))
# 		(s1 > 0) ⊻ (s3 > 0) && (b1 = !b1; push!(c1, i2))
	end#»»
	done = falses(length(loop))
	newloops = Vector{Int}[]; sides = falses(0)
# 	newloops = (Vector{Int}[], Vector{Int}[])
	# we could also do something intelligent with a stack, but whatever:
	for startindex in eachindex(loop)#««
		i = startindex
		x0 = loop[i][1]
		(iszero(x0) || done[i]) && continue
		l = Int[]; push!(newloops, l); push!(sides, x0 < 0)
		while !done[i]
			done[i] = true
			push!(l, i)
			x = loop[i][1]
			if iszero(x)
				i = x0 > 0 ? next1[i] : next0[i]
				push!(l, i)
			end
			i = nextindex(i)
		end
	end#»»
	return ([loop[l] for l in newloops], sides)
end
# Separators (parametrized bisectors) ««1
# Segment positions and branches««2
struct Branch; sign::Int8; end
@inline CornerTables.int(b::Branch) = b.sign
const _BAD_BRANCH = Branch(Int8(-128)) # stable by unary minus
@inline isbad(b::Branch) = b == _BAD_BRANCH
@inline Base.convert(::Type{Branch}, x::Real) = Branch(x)
@inline branch(x::Real) = Branch(iszero(x) ? 0 : (x > 0) ? 1 : -1)

@inline Base.show(io::IO, b::Branch) = print(io, int(b))
@inline Base.:-(b::Branch) = Branch(-int(b))
@inline Base.:*(a::Integer, b::Branch) = Branch(a*int(b))
@inline Base.:<((b1,r1)::Tuple{Branch,Real}, (b2,r2)::Tuple{Branch,Real}) =
	int(b1)*r1 < int(b2)*r2
@inline sqrt(b::Branch, x::Real) = int(b)*sqrtp(x)

"""    segments_position(seg1, seg2)
Given two segments `seg1`, `seg2` with positive determinant,
returns the relative position of both segments, as indices in this matrix:
(drawings show seg1 as horizontal, seg2 vertical)

        │    │    │
     ──     ───     ───

     ── │   (*)   │ ──

     ──     ───     ───
        │    │    │

(*) the middle entry corresponds to crossing segments and thus throws
the corresponding `CrossingSegments` exception.
"""
function segments_position((p1,q1)::Segment, (p2,q2)::Segment)#««
	u1, u2 = q1-p1, q2-p2
	@assert det2(u1, u2) > 0
	Dp2, Dq2 = det2(u1, p2-p1), det2(u1, q2-p1)
	Dp1, Dq1 = det2(u2, p1-p2), det2(u2, q1-p2)
	pos2 = 2 - (Dp2 ≥ 0) + (Dq2 ≤ 0)
	pos1 = 2 - (Dq1 ≥ 0) + (Dp1 ≤ 0)
	pos1 == pos2 == 2 && throw(CrossingSegments())
	return pos2, pos1
end#»»
"""    segments_quadrants(seg1, seg2)
Returns a quadruple of signs identifying the quadrants for the - and + branches
of the separator with left=seg1 and right=seg2.

        │   ↘│↗   │
     ── ↘   ───   ↗ ───
        
       ↘           ↗
     ── │   (*)   │ ──
       ↙           ↖

     ── ↙   ───   ↖ ───
        │   ↙│↖   │
"""
@inline segments_quadrants(seg1::Segment, seg2::Segment) =#««
	position_quadrants[segments_position(seg1, seg2)...]

const position_quadrants = SMatrix{3,3,NTuple{2,NTuple{2,Int8}}}([
	((-1,+1),(+1,-1)) ((-1,+1),(+1,+1)) ((-1,-1),(+1,+1));
	((-1,+1),(-1,-1)) (( 0, 0),( 0, 0)) ((+1,-1),(+1,+1));
	((+1,+1),(-1,-1)) ((+1,-1),(-1,-1)) ((+1,-1),(-1,+1))])
#»»
"identifies which branch passes through this quadrant; 0 if no branch"
@inline segments_whichbranch(seg1::Segment, seg2::Segment, quadrant) =#««
	position_whichbranch(segments_position(seg1,seg2), quadrant)
@inline function position_whichbranch(pos, quadrant)
	xym, xyp = position_quadrants[pos...]
	quadrant == xym && return Branch(-1)
	quadrant == xyp && return Branch(+1)
	return _BAD_BRANCH
end
#»»
# Data structure ««2
"""
    Separator

A structure holding the parametrization for the bisector between two sites
(both of them either a point or a vector),
The separator between `a` and `b` is parametrized by the distance to the sites,
and represented as two branches: the `+` branch sees the site `a` on its right.

This is either:
 - the bisector of two points: the line parametrized as
   origin ± √(r²-rmin²)\\*tangent  (with normal == 0);
 - the bisector of a line and a point outside the line:
   the parabola parametrized as
   origin ± √(r-rmin)\\*tangent + (r-rmin)\\*normal;
 - the bisector of a segment and a point on the line supporting the segment:
   this is a line, and an error if the point is in the interior of the segment;
 - the bisector of two non-crossing segments on secant lines:
   the union of two half-lines parametrized as
   origin + r\\*tangent, origin + r\\*normal (with rmin == 0);
 - the bisector of two touching, non-parallel segments is a straight line;
 - the bisector of two parallel segments: the central line,
   described as origin + r\\*tangent, with normal = [NaN, NaN].

The separator of two sites A,B is in all cases the union of two
infinite branches: the + branch sees A on its right (along increasing r),
while the - branch sees B on its right.
"""
struct Separator{T}
	origin::SVector{2,T}
	tangent::SVector{2,T}
	normal::SVector{2,T}
	rmin::T
end

# This avoids problems when comparing -0. to 0. (not bitwise identical):
@inline Base.:(==)(s1::Separator, s2::Separator) =
	(s1.origin == s2.origin) && (s1.tangent == s2.tangent) &&
		(s1.normal == s2.normal) && (s1.rmin == s2.rmin)

# predicates
@inline isparallel(sep::Separator) = any(isnan, sep.normal)
@inline isstraight(sep::Separator) = iszero(sep.normal)
@inline ishalflines(sep::Separator)= iszero(sep.rmin) && !iszero(sep.normal)

@inline Base.show(io::IO, sep::Separator) =
	@printf(io, "sep %s(o=[%.3g,%.3g], r₀=%.3g, t=[%.3g,%.3g], n=[%.3g,%.3g])",
		isparallel(sep) ? "═" :
		isstraight(sep) ? "─" :
		ishalflines(sep) ? "⋁" : "◡",
		sep.origin..., sep.rmin, sep.tangent..., sep.normal...)

"""    reverse(separator)

Given `separator(a,b)`, returns `separator(b,a)`, i.e. such that
`evaluate(sep′,b,s) = evaluate(sep,-b,s)."""
@inline Base.reverse(s::Separator)=
	Separator(s.origin,-s.tangent, s.normal, s.rmin)

# Constructors ««2
function Separator(a::AbstractVector, b::AbstractVector)# two points««
	c = SVector{2}(a+b)/2
	d = √(distance²(a,b))
	u = quarterturn(a-b)/(d)
	# guarantee: ‖tangent‖ = 1
	return Separator(c, u, zero(u), d/2)
end#»»
function Separator((p1,q1)::Segment, p2::AbstractVector; k=1)#««
	p1q1, p1p2 = q1-p1, p2-p1
	x1, x2, y2 = norm²(p1q1), p1q1⋅p1p2, det2(p1q1, p1p2)
	f = √(x1) # scale factor
	v = quarterturn(p1q1)
	if iszero(y2) # point is on the line supporting the segment
		# By convention, the separator in this case passes through p2.
		# In the only practical case (i.e. p2 is one of the segment ends),
		# this gives the correct answer, and this is also the answer making
		# tripoint computations consistant.
		x2 ≤ 0  && return Separator(p2,  k/f*v, zero(p1q1), zero(f))
		x2 ≥ x1 && return Separator(p2, -k/f*v, zero(p1q1), zero(f))
		throw(PointInSegment())
	end
	rmin = y2/(2*f)
	return Separator(p2 - y2*v/x1,
		k*sign(y2)*√(2*abs(y2)/f)*p1q1/f,
		sign(y2)*v/f, abs(rmin))
end
@inline Separator(a::AbstractVector, b::Segment) = Separator(b, a; k=-1) #»»
function Separator((p1,q1)::Segment, (p2,q2)::Segment)#««
	p1q1, p2q2 = q1-p1, q2-p2
	d = det2(p1q1, p2q2)
	if iszero(d) # special case: parallel separator
		l = det2(p1q1, p2-p1)
		(l < 0) && ((p1, q1, p1q1, l) = (q1, p1, -p1q1, -l))
		u = unit(p1q1); l = det2(u, p2-p1)
		return Separator((p1+p2)/2, u, SA[oftype(l, NaN), oftype(l, NaN)], l/2)
	end
	# both segments are un-ordered, so we swap p2, q2 if needed:
	d < 0 && ((p2q2, p2, q2, d) = (-p2q2, q2, p2, -d))
	c = lineinter(p1, q1, p2, q2)
	u1, u2 = √(norm²(p2q2))*p1q1/d, √(norm²(p1q1))*p2q2/d
	((xm, ym), (xp, yp)) = segments_quadrants((p1,q1), (p2,q2))
	return Separator(lineinter(p1,q1,p2,q2),
		(xp-xm)/2 * u1 + (yp-ym)/2 * u2,
		(xp+xm)/2 * u1 + (yp+ym)/2 * u2,
		zero(d))
end#»»
# Evaluation, interpolation««2

"""    evaluate(separator, branch, r)

Returns the point on the separator situated at distance `r` from both
sites and on the branch given by sign `s` (either + or -).
The `+` branch sees `a` on its right and `b` on its left.
"""
@inline function evaluate(sep::Separator, b::Branch, r) # b is a sign
	ishalflines(sep) && return sep.origin + r*(sep.normal + int(b)*sep.tangent)
	isstraight(sep) &&
		return sep.origin + sqrt(b, r^2-sep.rmin^2)*sep.tangent
	# parabola arc # WARNING: sep.origin is *NOT* the parabola apex
	return sep.origin + r*sep.normal + sqrt(b, r-sep.rmin)*sep.tangent
end

"""    approximate(separator, r1, r2, atol)

Approximates the + branch of parabolic separator with precision `atol`
by a polygonal path. Returns vector of distance parameters."""
function approximate(sep::Separator, r1, r2, atol)
	iszero(sep.tangent) && return [r1, r2] # degenerate case
	isstraight(sep) && return [r1, r2]
	ishalflines(sep) && return [r1, r2]
	nt = √(norm²(sep.tangent))
	x1, x2 = nt*√(r1-sep.rmin), nt*√(r2-sep.rmin)
	x = approxparabola(sep.rmin, x1, x2, atol)
	y = sep.rmin .+ (x ./ nt) .^ 2
	y[begin] = r1; y[end] = r2
	return y
end

"""    vectors_angles(u, v)

Returns (angle(u), angle(v)-angle(u)); second angle always in [0,2π[."""
@inline vectors_angles(u::AbstractVector, v::AbstractVector) =
	(atan(u[2], u[1]), atan(det2(u,v), u⋅v))
# Tripoints ««1
@inline _BAD_TRIPOINT(x) = (oftype(x,NaN), _BAD_BRANCH,_BAD_BRANCH,_BAD_BRANCH)
# docstring ««2
"""    tripoint(c1,c2,c3)

This computes the tripoint (equidistant point) of a triple of cells,
each of which is either a point or an (unioriented) segment.
The cells are cyclically ordered:
tripoint(a,b,c) == tripoint(b,c,a) ≠ tripoint(c,b,a).
The data is returned as `(radius, branch1, branch2, branch3)`.
If no such tripoint exists, `nan, 2,2,2` is returned.

The branch positions are returned as `Int8`, encoded as:
branch1 = +1 iff the tripoint lies on branch seeing c1 on its left
and c2 on its right, i.e. **c1↑c2**
(the arrow marks the direction of increasing r on this branch).
Likewise, branch2 is +1 for c2↑c3, and branch3 is +1 for c3↑c1.

For example, the center of an equilateral triangle has branches +1,+1,+1.
"""
function tripoint end #»»
function tripoint(a::AbstractVector, b::AbstractVector, c::AbstractVector)#««
	ab, bc, ca = b-a, c-b, a-c
	det2(ca, ab) > 0 || return _BAD_TRIPOINT(ab[1])
	r = √(norm²(ab)*norm²(ca)*norm²(bc))/(2*abs(det2(ab,ca)))
	return r, -branch(bc ⋅ ca), -branch(ca ⋅ ab), -branch(ab ⋅ bc)
end#»»
function tripoint((p1,q1)::Segment, p2::AbstractVector, p3::AbstractVector)#««
# WLOG assume L1 has equation (y=0), P2=(0,a) and P3=(d,b).
# * If b < 0 then S12 and S23 do not intersect.
# * If b = 0 then S13 is a half-line XXX.
# * If b = a then the three separators meet at a single point,
# which is either H if d<0 or H' if d>0.
# * If b > 0, b≠a then the three separators meet in two points.
# H is the left-most point iff b > a; the right-most point iff a < 1.
#
#	H=(x,y) satisfies (b-a)x²+2adx+b+a(ab-b²-d²) = 0 and y=(x²+a²)/(2a).
# Let Δ=ab(d²+(b-a)²); then x=(ad+√Δ)/(a-b) whenever b≠1.
# (When b=1 and a<0, x=a/2; when b=1 and a>0, H does not exist).
#
# summary:
# branches are ++- for (b>1) and (a>-√(b-1))
#              -++ for (b>1) and (a>√(b(1-b)))
#              +-+ for (a<0) and (b< a^2/4)
#              +++ otherwise
	u1, u2, u3 = q1-p1, p2-p1, p3-p1
	x1, x2, y2, x3, y3 = norm²(u1), u1⋅u2, det2(u1,u2), u1⋅u3, det2(u1,u3)
	dx = x3-x2
	f = √(x1) # scale factor
	if iszero(y2)#««
		iszero(y3) && return _BAD_TRIPOINT(x1)
		if iszero(x2) # p2 is start of segment 1
			# possible positions for p3; only the lower ones make (123) tripoints:
			#  X╱│\ X  /
			#  ╱ │ `--'
			#   2┝━━━━1━━
			#  ╲ │ ,--.
			#  3↘↓↙  3 \
			#    ↓
			y3 > 0 && return _BAD_TRIPOINT(x1)
			s = branch(x3)
			return (x3^2+y3^2)/(-2*f*y3), Branch(1), -s, s
		elseif x2 == x1
			# possible positions for p3: only the top ones are (123) tripoints
			#        ↑ 
			#  \ 3  ↗↑↖ 3
			#   `--' │ ╲
			# ━━━1━━━┥2
			#   ,--. │ ╱
			#  /  X \│╱ X
			y3 < 0 && return _BAD_TRIPOINT(x1)
			s = branch(x3-x1)
			return ((x3-x1)^2+y3^2)/(2*f*y3), Branch(1), s, -s
		else
			throw(PointInSegment())
		end
	end#»»
	if iszero(y3)#««
		if iszero(x3) # p3 is start of segment 1
			y2 < 0 && return _BAD_TRIPOINT(x1)
			s = branch(x2)
			return (x2^2+y2^2)/(2*f*y2), s, -s, Branch(1)
		elseif x3 == x1
			y2 > 0 && return _BAD_TRIPOINT(x1)
			s = branch(x2-x1)
			return ((x2-x1)^2+y2^2)/(-2*f*y2), -s, s, Branch(1)
		else
			throw(PointInSegment())
		end
	end#»»
	# both points lie outside the segment; ensure they are on the same side,
	# and identify this side as the positive one:
	(y2 < 0) && ((y2,y3,dx) = (-y2,-y3,-dx))
	(y3 < 0) && return _BAD_TRIPOINT(x1)
	r = (y3 == y2) ?
		(4*y2^2+dx^2)/(8*f*y2) : # this limit formula is only valid when dx≤0
	let t = dx^2+(y2-y3)^2
		(2*dx*√(y2*y3*t) + (y2+y3)*t)/(2*f*(y2-y3)^2)
	end
	if dx ≥ 0
		(y3 == y2) && return _BAD_TRIPOINT(x1)
		return r, branch(y3-y2), branch(1), branch(y2-y3)
	end
	
	s0 = y3^2+dx^2-y2*y3
	s1 = 4*y2*y3 - dx^2
	s2 = y2^2+dx^2-y2*y3
	return r, branch(s0), branch(s1), branch(s2)
end#»»
function tripoint((p1,q1)::Segment, (p2,q2)::Segment, p3::AbstractVector)#««
	v1, v2 = q1-p1, q2-p2
	a1 = det2(v1, p3-p1)
	(a1 < 0) && ((a1, v1, p1, q1) = (-a1, -v1, q1, p1))
	c, s, l1, l2 = v1⋅v2, det2(v1, v2), √norm²(v1), √norm²(v2)
	# Reorient so that s ≥ 0 and a1 ≥ 0:
	(s < 0) && ((c, s, v2, p2, q2) = (-c, -s, -v2, q2, p2))
	a2 = det2(v2, p3-p2)

	if iszero(s) # parallel segments case««
		h1 = det2(v1, p2-p1)
		(h1 < 0) && ((h1, a1, p1, v1) = (-h1, -a1, q1, -v1)) # swap p1, q1
		(a1 < 0) || (a1 > h1) && return _BAD_TRIPOINT(c)
		if iszero(a1) # p3 ∈ line1; must be left of segment1
			x3 = (p3-p1)⋅v1 # x3/l1
			x3 > 0 && return _BAD_TRIPOINT(c)
			return h1/(2l1), Branch(0), Branch(-1), Branch(1)
		elseif a1 == h1
			xp3, xq3 = (p3-p2)⋅v1, (p3-q2)⋅v1
			(xp3 > 0 || xq3 > 0) && return _BAD_TRIPOINT(c)
			return h1/(2l1), Branch(0), Branch(1), Branch(-1)
		end
		return h1/(2l1), Branch(0), Branch(1), Branch(1)
	end#»»
	# Taking coordinates: u₁=(1,0), u₂=(c,s), p3=((c a₁-a₂)/s, a₁);
	# consider the branch of B₁₂ directed by (ηu₂-εu₁)/s = [(ηc-ε)/s, η].
	# (in general, η=sign(a1)=1 and ε=sign(a2)).
	#
	# This is equidistant from L1 and P3 whenever
	# ((1-εηc)r)² - 2(ηa₁+εa₂)(1-εηc)r + (a₁²+a₂²-2a₁ a₂ c) = 0;
	# Δ'= 2a₁a₂(εη+c) whence (1-εηc)r = ηa₁+εa₂±√Δ'; geometry forces sign to -εη
	#
	# ε,η represent the quadrant in which we look for the tripoint:
	e, f = Int(sign(a2)), Int(sign(a1))
	if iszero(e)
		iszero(f) && return a1, Branch(1), Branch(1), Branch(1) # trivial case
		e = (a1 ≥ det2(v1,q2-p1)) ? 1 : (a1 ≤ det2(v1,p2-p1)) ? -1 :
			throw(PointInSegment())
	elseif iszero(f)
		f = (a2 ≥ det2(v2,p1-p2)) ? 1 : (a2 ≤ det2(v2,q1-p2)) ? -1 :
			throw(PointInSegment())
	end
	b1 = -segments_whichbranch((p1,q1),(p2,q2), (-e, f)) # minus bc L=seg1, R=seg2
	isbad(b1) && return _BAD_TRIPOINT(c)
	d = √(2*a1*a2*(e*l1*l2+c)) # this has normalization factor l1*l2
	r = (f*a1*l2+e*(a2*l1-d))/(l1*l2-e*f*c)
	# compute position relative to both parabolas
	g1, g2 = 2*l2^2*a1-(l1*l2+c)*a2, 2*l1^2*a2-(l1*l2+c)*a1
	b2, b3 = (e != f) ? (Branch(1), Branch(1)) : (f*branch(g1), f*branch(g2))
	return r,b1,b2,b3
end#»»
function tripoint((p1,q1)::Segment, (p2,q2)::Segment, (p3,q3)::Segment)#««
	v1, v2 = q1-p1, q2-p2
	s12 = det2(v1, v2)
	if iszero(s12)
		error("parallel case: not implemented")
	end
	# we reorient so that (seg1,seg2) and (seg1,seg3) are positive angles
	(s12 < 0) && ((s12, v2, p2, q2) = (-s12, -v2, q2, p2))
	v3 = q3-p3
	s31 = det2(v3, v1)
	iszero(s31) && return tripoint((p3,q3),(p1,q1),(p2,q2)) # parallel case
	(s31 > 0) && ((s31, v3, p3, q3) = (-s31, -v3, q3, p3))

	s23 = det2(v2, v3)
	iszero(s23) && return tripoint((p2,q2),(p3,q3),(p1,q1)) # parallel case
	c12, c23, c31 = v1⋅v2, v2⋅v3, v3⋅v1
	l1, l2, l3 = √norm²(v1), √norm²(v2), √norm²(v3) # normalization factors
	# zij is the coordinate of line(i)∩line(j) on line(i):
	# let vᵢ=qᵢ-pᵢ, cᵢⱼ=vᵢ⋅vⱼ, sᵢⱼ=<vᵢ,vⱼ>,
	# x(A)=u₁.(A-p₁), y(A)=<u1,A-p₁>; xqᵢ = xpᵢ+c₁ᵢ, yqᵢ = ypᵢ+s₁ᵢ
	# then det([xpᵢ ypᵢ 1; xpᵢ+c₁ᵢ ypᵢ+s₁ᵢ 1;x(I₁ᵢ) 0 1]) = 0
	# so that z₁ᵢ = x(I₁ᵢ) = xpᵢ - c₁ᵢ/s₁ᵢ y₁ᵢ.
	x12, y12 = v1⋅(p2-p1), det2(v1, p2-p1)
	x13, y13 = v1⋅(p3-p1), det2(v1, p3-p1)
	z12 = x12 - c12/s12*y12
	z13 = x13 + c31/s31*y13
	a1 = abs(z12-z13)/l1

	x23, y23 = v2⋅(p3-p2), det2(v2, p3-p2)
	x21, y21 = v2⋅(p1-p2), det2(v2, p1-p2)
	z23 = x23 - c23/s23*y23
	z21 = x21 + c12/s12*y21
	a2 = abs(z23-z21)/l2

	x31, y31 = v3⋅(p1-p3), det2(v3, p1-p3)
	x32, y32 = v3⋅(p2-p3), det2(v3, p2-p3)
	z31 = x31 - c31/s31*y31
	z32 = x32 + c23/s23*y32
	a3 = abs(z31-z32)/l3

	e = Int(sign(z13-z12))
	iszero(e) && throw(ConcurrentLines())
	if s23 < 0 # the incenter
		r = sqrtp((a1+a2-a3)*(a2+a3-a1)*(a3+a1-a2)/(a1+a2+a3))/2
		# the situation must be this one (e=1) or its converse:
		# line2 ↖  ↗ line3
		#        ╲╱
		#        ╱╲
		#     q3╱  ╲q2
		#      ╱    ╲
		#   p3╱      ╲p2
		#    ╱        ╲
		# ───────────────→ line1
		#      p1   q1
		b1 = segments_whichbranch((p1,q1),(p2,q2), (-e, e))
		b3 = -segments_whichbranch((p1,q1),(p3,q3), (e, e))
		b2 = -segments_whichbranch((p3,q3),(p2,q2), (-e, -e))
		any(iszero, (b1,b2,b3)) && return _BAD_TRIPOINT(c12)
		return r, b1, b2, b3
	else # an excenter (depending on the segment positions)
		pos12 = segments_position((p1,q1), (p2,q2))
		pos13 = segments_position((p1,q1), (p3,q3))
		pos23 = segments_position((p2,q2), (p3,q3))
		# try excenters 1,2,3 in turn
		for (q1,q2,q3,s) in (
			((+e,-e),(-e,-e),(-e,-e),(-1,1,1)),
			((-e,+e),(-e,+e),(-e,+e),(1,-1,1)),
			((+e,+e),(+e,+e),(+e,-e),(1,1,-1)))
			b1 =-position_whichbranch(pos12, q1) # minus signs because L=1,R=2
			b2 =-position_whichbranch(pos23, q2)
			b3 = position_whichbranch(pos13, q3)
			any(iszero, (b1,b2,b3)) && continue
			r = sqrtp((s⋅(a2,a3,a1))*(s⋅(a3,a1,a2))*(a1+a2+a3)/(s⋅(a1,a2,a3)))/2
			return r, b1, b2, b3
		end
		return _BAD_TRIPOINT(c12)
	end
end#»»
# Triangulation««1
# Cell location««2
@inline geometrictriangle(t::AbstractTriangulation, points, q::Node) =
	(points[int(cell(t,q,1))], points[int(cell(t,q,2))], points[int(cell(t,q,3))])

"""    findnode(triangulation, points, point)

Returns the index of the node closest to this point.
(In a triangulation of points, this is the triangle containing the point).
"""
function findnode(v::AbstractTriangulation{J}, points, point) where{J}#««
	q = rand(eachnode(v))
	c = 0
	while true
		c+= 1; @assert c ≤ 1e3
		# this guarantees that v[i] = tail(side(q, i)) for i in 1:3:
		p1, p2, p3 = geometrictriangle(v, points, q)
		isleft(p1, point, p2) && (q = adjnode(v, q, 1); continue)
		isleft(p2, point, p3) && (q = adjnode(v, q, 2); continue)
		isleft(p3, point, p1) && (q = adjnode(v, q, 3); continue)
		return q
	end
end#»»

# Triangulation constructor ««2
"""    addpoint!(v, c): creates a cell for a point"""
function addpoint!(v::AbstractTriangulation, points, c, point)#««
	q0 = findnode(v, points, point)
	stack = [insert!(v, q0, c)...]
	while !isempty(stack)
		e = pop!(stack)
		@assert left(v, e) == c
		o = opposite(v, e)
		int(o) ≤ 3 && continue # this is the phony outer node
		q = node(o)
		isincircle(geometrictriangle(v, points, q)..., point) || continue
		ono, opo = opposite(v, next(o)), opposite(v, prev(o))
		if left(v, ono) == c
			error("closing cell to the right: ono=$ono")
		end
		if left(v, opo) == c
			error("closing cell to the left: opo=$opo")
		end
		# XXX check if we are closing a cell here!
		e1, e2 = flip!(v, o)
		push!(stack, e1, e2)
	end
end#»»

"""    triangulate(points)

Returns a triangulation of this set of points,
as a list of triples of integers."""
function triangulate(points; kw...)
	np = length(points)
	t = CornerTable{Int32}(points; kw...)
	# remove all superfluous nodes & cells ««
	# the nodes are sorted this way:
	# - inner nodes
	# - convex hull
	# - 1 backwards outer node
	k = nnodes(t)
	swapnodes!(t, Node(1), Node(k))
	k-= 1
	fakecells = np+1:np+3
	for i in nnodes(t)-1:-1:1; q = Node(i)
		w = int.(triangle(t,q))
		any(>(np), Int.(triangle(t, q))) || continue
		swapnodes!(t, q, Node(k))
		k-= 1
	end
	# »»
	resize!(points, np)
	nnodes!(t, k)
	return [(int(cell(t,q,1)),int(cell(t,q,2)),int(cell(t,q,3)))
		for q in eachnode(t)]
end#»»

function CornerTable{J}(points; extra = 0) where{J}
	# Builds a Delaunay triangulation of these points using Bowyer-Watson's
	# algorithm:
	T = float(eltype(eltype(points)))
	np = length(points)
	# build initial dihedron between points np+1, np+2, np+3 ««
	# node 1 is a fake node (“at infinity”);
	v = CornerTable{J}(J[4,6,5,1,3,2],J(np).+J[2,1,3,1,2,3], zeros(J, np+3))

	anyedge!(v, Cell(np+1), Edge(J(4)))
	anyedge!(v, Cell(np+2), Edge(J(5)))
	anyedge!(v, Cell(np+3), Edge(J(6)))
	m = maximum(abs(x) for p in points for x in p) + extra + 1
	append!(points, [SA[0,-3m], SA[3m,2m], SA[-3m,2m]])
  #»»
	# incrementally add all points ««
	Random.seed!(0)
	for c in Cell(J(1)):Cell(J(np)) # Random.randperm(np)
		addpoint!(v, points, c, points[int(c)])
	end #»»
	return v
end
# Voronoi diagram: topology ««1
# Data structure and accessor functions««2
abstract type AbstractVoronoi{J} <: AbstractTriangulation{J} end
"""    VoronoiDiagram{J,T}

Encodes the triangulation of a set of points and segments,
as well as the geometric structure of the Voronoi diagram.

Type parameters:
 - `J`: integer index type;
 - `T`: real distance type.
"""
struct VoronoiDiagram{J,T} <: AbstractVoronoi{J}
	# FIXME: make P = SVector{2,T}
	triangulation::CornerTable{J}
	points::Vector{SVector{2,T}}
	segments::Vector{NTuple{2,J}} # indices into points
	geomnode::Vector{SVector{2,T}} # indexed by nodes
	noderadius::Vector{T} # indexed by nodes
	separator::Vector{Separator{T}} # indexed by edges
	branch::Vector{Branch} # indexed by edges
	neighbours::Vector{J} # indexed by points

	@inline VoronoiDiagram{J,T}(triangulation::CornerTable, points, segments
		) where{J,T} =
		new{J,T}(triangulation, points, segments,
			Vector{SVector{2,T}}(undef, nnodes(triangulation)),
			Vector{T}(undef, nnodes(triangulation)),
			Vector{Separator{T}}(undef, nedges(triangulation)),
			Vector{Branch}(undef, nedges(triangulation)),
			zeros(J, length(points)), # neighbours
			)
end

@inline CornerTables.triangulation(v::VoronoiDiagram) = v.triangulation

@inline npoints(v::VoronoiDiagram) = length(v.points)
@inline nsegments(v::VoronoiDiagram) = length(v.segments)
@inline ispoint(v::AbstractVoronoi, c::Cell) = int(c) ≤ npoints(v)
@inline issegment(v::AbstractVoronoi, c::Cell) = int(c) > npoints(v)
@inline point(v::VoronoiDiagram, c::Cell) = v.points[int(c)]
@inline cellsegment(v::VoronoiDiagram, c::Cell) =
	Cell.(v.segments[int(c)-npoints(v)])
@inline segment(v::VoronoiDiagram, c::Cell) =
	tuple((point(v,i) for i in cellsegment(v,c))...)
@inline geometricnode(v::VoronoiDiagram, q::Node) = v.geomnode[int(q)]
@inline geometricnode!(v::VoronoiDiagram, l::Pair{<:Node,<:AbstractVector}...)=
	for (q, p) in l; v.geomnode[int(q)] = p; end
@inline noderadius(v::VoronoiDiagram, q::Node) = v.noderadius[int(q)]
@inline noderadius!(v::VoronoiDiagram, l::Pair{<:Node,<:Real}...) =
	for (q, r) in l; v.noderadius[int(q)] = r; end
@inline branch(v::VoronoiDiagram, e::Edge) = v.branch[int(e)]
@inline branch!(v::VoronoiDiagram, l::Pair{<:Edge,Branch}...) =
	for (e, b) in l; v.branch[int(e)] = b; end
@inline separator(v::VoronoiDiagram, e::Edge) = v.separator[int(e)]
@inline separator!(v::VoronoiDiagram, l::Pair{<:Edge,<:Separator}...) =
	for (e, s) in l; v.separator[int(e)] = s; end
@inline separators!(v::VoronoiDiagram, e::Edge, o::Edge, s::Separator) =
	separator!(v, e=>s, o=>reverse(s))


@inline edgedata(v::VoronoiDiagram, e::Edge) = separator(v,e), branch(v,e)
@inline edgedata!(v::VoronoiDiagram, e::Edge, (s, b)) =
	v.separator[int(e)], v.branch[int(e)] = s, b
# right-side segments if side=false, left-side if side=true:
@inline sidesegments(v::VoronoiDiagram{J}, side) where{J} =
	J(npoints(v)+1+side):J(2):ncells(v)

# Updates of `CornerTables` functions preserving geometry ««2
@inline geometrictriangle(v::VoronoiDiagram, q::Node) =
	geometrictriangle(v.triangulation, v.points, q)

function CornerTables.nnodes!(v::VoronoiDiagram, n)#««
	nnodes!(CornerTables.triangulation(v), n)
	ne = 3n
	resize!(v.geomnode, n)
	resize!(v.noderadius, n)
	resize!(v.separator, 3n)
	resize!(v.branch, 3n)
end#»»
function CornerTables.flip!(v::VoronoiDiagram, e::Edge)#««
	n, p, o = next(e), prev(e), opposite(v, e)
	no, po = next(o), prev(o)
	r = invoke(flip!, Tuple{AbstractTriangulation, Edge}, v, e)
	# recompute node radius for the two modified nodes
	d = edgedata(v, p)
	edgedata!(v, p, edgedata(v, n))
	edgedata!(v, n, edgedata(v, po))
	edgedata!(v, po, edgedata(v, no))
	edgedata!(v, no, d)
	edgedata!(v, e)
	edgedata!(v, o)
	nodedata!(v, node(e))
	nodedata!(v, node(opposite(v,e)))
	return r
end#»»
function CornerTables.insert!(v::VoronoiDiagram, q::Node, c::Cell)#««
	# call the “parent” method for basic edge manipulations
	e, n, p = sides(q)
	e0, e1, e2 = invoke(Base.insert!,
		Tuple{AbstractTriangulation, Node, Cell}, v, q, c)
	# then update all geometric information:
	edgedata!(v, e1, edgedata(v, n))
	edgedata!(v, e2, edgedata(v, p))
	edgedata!(v, n)
	edgedata!(v, p)
	edgedata!(v, next(e1)) # facing prev(e2)
	nodedata!(v, node(e))
	nodedata!(v, node(e1))
	nodedata!(v, node(e2))
end#»»

# Cell location functions ««2
"""    influences(v, a, b, point)
Returns true iff segment [a,b] sees the given point. """
function influences(v::AbstractVoronoi, i, j, q)
	a,b,g = point(v, i), point(v,j), geometricnode(v,q)
	ab, ag = b-a, g-a
	return (0 < dot(ab,ag) < dot(ab,ab))
end

"""    findrootnode(v, a, b)

Returns the node at which the segment [a,b] is inserted."""
function findrootnode(v::AbstractVoronoi, a,b)
	p, q = point(v,a), point(v,b)
	emin, dmin = nothing, nothing
# 	println("\e[35;7m finding root node for segment ($a,$b)\e[m")
# 	display((v,a))
	for e in star(v,a) # the ends of the segments are already cells...
# 		display((v,node(e)))
		influences(v,a,b,node(e)) || continue
# 		println("  $e is influenced by ($a,$b)")
		d = segdistance²(p,q,geometricnode(v, node(e)))
		(emin == nothing || d < dmin) && ((emin, dmin) = (e, d))
	end
	@assert emin ≠ nothing
# 	println("\e[35m return $(node(emin))\e[m:"); display((v,node(emin)))
	return node(emin)
end

# Segment insertion ««2
"""    edgecapture(v,e)

Returns a positive value if this edge must be flipped, i.e. if all its points
are closer to right(e) than to tail(e)/head(e)."""
function edgecapture(v::VoronoiDiagram, e::Edge)#««
	# edge is captured if for all r, d(p(r), L) < r
	# i.e. min(r² - d²(p(r), L)) > 0
# 	println("\e[32;7medgecapture($(right(v,e)) = $(cellsegment(v,right(v,e))) -> $e = $(tail(v,e))⋯$(head(v,e))\e[m")
	o = opposite(v,e)
	be, bo = branch(v, e), branch(v, o)
	bo == _BAD_BRANCH && return true
	@assert issegment(v, right(v, e))
	(a, b) = cellsegment(v, right(v,e))
	q = node(e)
	influences(v, a, b, q) || return false
	sep = separator(v, e)
	be, bo = branch(v,e), branch(v,o)
	re, ro = noderadius(v,node(e)), noderadius(v,node(o))
# 	z = capture(sep, be, re, -bo, ro, l)
# 	println("  \e[32mcapture($e) = $z\e[m")
	f = (be, re) < (-bo, ro)
# 	f = (-bo > be) || (-bo == be && ro > re)
# 	if f != (z > 0)
# 		println("\e[35;7m  (trying to connect $(right(v,e))→$(left(v,e)), disconnect $(tail(v,e))/$(head(v,e))) found f=$f, z=$z\e[m\e[35m")
# 	println("separator $e/$o = $(tail(v,e))/$(tail(v,o)) is $sep")
# 	println("  node($e)=$(node(e)) (in $(tail(v,e))) has parameter ($be, $re)")
# 	println("   geometricnode[$e] = $(geometricnode(v,node(e)))")
# 	println("  node($o)=$(node(o)) (in $(tail(v,o))) has parameter ($bo, $ro)")
# 	println("   geometricnode[$o] = $(geometricnode(v,node(o)))")
# 	println(evaluate(separator(v,e), be, re))
# 	println(evaluate(separator(v,o), bo, ro))
# 	println(evaluate(separator(v,e), -be, re))
# 	display((v,node(e)))
# 	display((v,node(o)))
# 
# 	println("\e[35;1m view from $(tail(v,e)): $(-bo),$ro  increases to $be,$re\e[m")
# 		error("stop")
# 	end
	return f
# 	return z > 0
end#»»
function addsegment!(v::VoronoiDiagram, c::Cell)#««
	a,b = cellsegment(v, c)
# 	println("\e[31;1;7minserting segment $c = ($a,$b)\e[m")
# 	display(v)
	q0 = findrootnode(v, a, b)
# 	print("\e[31mroot node is $q0:\e[m "); display((v,q0))
	stack = [opposite(v, e) for e in sides(q0)]
	insert!(v, q0, c)
# 	println("\e[31;7m after insert!($q0, $c) =$stack:\e[m") # ; display(v)
	# now grow the cell by repeatedly flipping edges as needed
	while !isempty(stack)
# 		println("\e[36;1m stack: $([left(v,e) for e in stack])\e[m")
		e = pop!(stack)
# 		println("\e[36;7m current status of graph (e=$e, stack=$stack)\e[m\n")
# 		display(v)
		tail(v,e) == c && continue # closing loop around $(head(v,e))
# 		tail(v,e) == c && error("error: closing loop around $(head(v,e))")
		@assert right(v, e) == c
# 		println("examining outgoing edge $e: l=$(left(v,e)) h=$(head(v,e)) t=$(tail(v,e))")
		o = opposite(v, e)
# 		println("  branches:  at e: $(branch(v,e)) at o=$o: $(branch(v,o))")
# 		int(o) ≤ 3 && continue # this is the phony outer node
		q = node(e)
# 		influences(v, a, b, q) || println("   segment($a,$b) does not see node $q =$(geometricnode(v,q))")
# 		influences(v, a, b, q) || continue
		edgecapture(v, e) || continue
# 		println("   \e[7m flipping edge $e: connect $(right(v,e))->$(left(v,e)), disconnect $(head(v,e)) - $(tail(v,e))\e[m")
		if left(v,e) == c
			error("  closing loop around cell $(tail(v,e))")
		end
		@assert int(e) ≠ 7
		e1, e2 = opposite(v, next(e)), opposite(v, prev(e))
		flip!(v, e)
# 		println("  flip done")
# 		println("  now e has h=$(head(v,e)) t=$(tail(v,e)) l=$(left(v,e)) r=$(right(v,e))")
# 		println("  now right($e1) = $(right(v,e1)); right($e2) = $(right(v,e2))")
		# XXX fixme: move separators for all concerned edges
		# & compute new separators
		push!(stack, e1, e2)
	end
end#»»
# Constructor ««2
@inline VoronoiDiagram(points::AbstractVector{P}, segments=[];kw...) where{P} =
	VoronoiDiagram{Int32,float(eltype(P))}(points, segments; kw...)

function VoronoiDiagram{J,T}(points, segments; extra=0) where{J,T}#««
	np, ns = length(points), length(segments)
	v = VoronoiDiagram{J,T}(CornerTable{J}(points), points, segments)

# 	println("\e[1;7m after triangulating all points:\e[m")
# 	global V=v
# 	display(v)
# 	gnuplot(v)
	# update geometric information ««
	for e in eachedge(v)
		e < opposite(v, e) && edgedata!(v, e)
	end
	for q in eachnode(v)
		nodedata!(v, q)
	end # »»

	ncells!(v, ncells(v) + ns)
	triangulation(v).anyedge[np+4:end] .= 0
	# incrementally add all segments ««
	for c in Cell(J(np+4)):Cell(J(np+ns+3)) # Random.randperm(ns)
		addsegment!(v, c)
	end
	#»»
# 	# remove all superfluous nodes & cells ««
# 	# the nodes are sorted this way:
# 	# - inner nodes
# 	# - convex hull
# 	# - 1 backwards outer node
# 	k = nnodes(v)
# 	swapnodes!(v, Node(1), Node(k))
# 	k-= 1
# 	fakecells = np+1:np+3
# 	for i in nnodes(t)-1:-1:1; q = Node{J}(i)
# 		w = Int.(triangle(t,q))
# 		w[1] ∈ fakecells || w[2] ∈ fakecells || w[3] ∈ fakecells || continue
# # 		any(>(ntotal), Int.(triangle(t, q))) || continue
# 		swapnodes!(v, q, Node(k))
# 		k-= 1
# 	end
# # 	nnodes!(t, k)
# 	resize!(points, np)
# 	# »»

# 	println("\e[1;7m before splitting segments:\e[m"); display(v)
#
	# split segments in two
	for s in segments, a in s; v.neighbours[a]+=one(J); end
	splitsegments!(v)

# 	println("segment = $(v.segments)")
# 	# connect each point to its right-side segment, if any:
# 	nc = zeros(J, np)
# 	for i in 1:2:length(v.segments); (a,b) = v.segments[i]
# 		nc[a] = iszero(nc[a]) ? np+3+i : -1
# 	end
# 	nc .= max.(0, nc)
# 	# and use this to build a doubly-linked list of cells
# 	for i in 1:2:ns; (a1,a2) = v.segments[i]
# 		j = nc[a2]; iszero(j) && continue
# 		c12, c21, c23, c32 = Cell.(J.((np+3+i, np+4+i, j, j+1)))
# 		# two configurations are now possible:
# 		# turn-left:              turn-right:
# 		#            |             _c21_|   c2    |_c32_
# 		#          ,' `.            c12  `.     .'  c23
# 		#   _c21_,'     `._c32_            `. ,'
# 		#    c12 |  c2   | c23               |
# 		# (the anyedge() data for segment cells is always the segment-split edge)
# 		b12, b32 = (before(v, anyedge(v, c)) for c in (c12, c32))
# 		a21, a23 = (after(v, anyedge(v, c)) for c in (c21, c23))
# 		bb12, bb32 = before(v, b12), before(v, b32)
# 		if head(v, bb32) == c21 # turn-left configuration
# 			nextedge!(v, b12, bb32)
# 			prevedge!(v, after(v, anyedge(v, c23)))
# 		else # turn-right
# 			@assert head(v, bb12) == c23
# 			nextedge!(v, bb12, b32)
# 			prevedge!(v, after(v, anyedge(v, c21)))
# 		end
# 	end
	return v
end#»»
# Split segments ««2
"""
    splitsegments!(voronoidiagram)

Splits segment cells in two parts depending on the orientation of the segment:
(right of segment, left of segment).
"""
function splitsegments!(v::VoronoiDiagram{J}) where{J}#««
	np = length(v.points)
	ns = length(v.segments)
	# sanity check:
	for (a,b) in v.segments; @assert (b,a) ∉ v.segments; end
# 	display(v)
	ncells!(v, np+2ns)
	# rename cells to alternate segments and their opposites;
	# since we are swapping with undefined cells,
	# there is no need to worry about side-effects.
	for i in ns:-1:1
		movecell!(v, Cell(np+i), Cell(np+2i-1))
	end
	origsegments = copy(v.segments)
	sizehint!(empty!(v.segments), 2ns)
	for (i, (a,b)) in pairs(origsegments)
		push!(v.segments, (a,b), (b,a))
	end
	# now split each cell in two
	for i in 1:ns
		c12, c21 = Cell(np+2i-1), Cell(np+2i)
		(c1,c2) = cellsegment(v, c12)
		(p1,p2) = point(v, c1), point(v, c2)
# 		println("\e[7m splitting cell $c12 = ($c1,$c2)\e[m")
# 		showcell(stdout, v, c12)
		e2 = anyedge(v, c12)
		while head(v, e2) ≠ c2; e2 = after(v, e2); end
		e1 = e2
		while head(v, e1) ≠ c1; tail!(v, e1, c21); e1 = after(v, e1); end
		# split the cell by inserting two new nodes
		# (= split vertex c12 of the dual triangulation)
		#
		#   c1  o1|e1  c12 e2|o2 c2 becomes:
		#
		#   o1│e12  c21 e2│
		#     │    e21    │e23
		#  c1 e1─────────e2 c2
		#  e13│    e11    │
		#     │e1      e22│o2
		#
		o1, o2 = opposite(v, e1), opposite(v, e2)
		q1, q2 = newnodes!(v, 2)
		e11, e12, e13 = side(q1,1), side(q1,2), side(q1,3)
		e21, e22, e23 = side(q2,1), side(q2,2), side(q2,3)
		tail!(v, e11=>c12, e12=>c21, e13=>c1, e21=>c21, e22=>c12, e23=>c2)
		opposites!(v, e11=>e21, e12=>o1, e13=>e1, e22=>o2, e23=>e2)
		anyedge!(v, c12, e11); anyedge!(v, c21, e21)
# 		display((v, q1)); display((v,q2)); display((v, c12)); display((v, c21))
		# fix geometric information:
		seg12, seg21 = segment(v, c12), segment(v, c21)
		separators!(v, e11, e21, Separator(seg12, seg21))
		separator!(v, e12 => separator(v, e1), e22 => separator(v, e2),
			e13 => separator(v, o1))
		edgedata!(v, o1) # this also fixes e12
		edgedata!(v, e2) # this also fixes e23
		for e in star(v, e21); edgedata!(v, e); end
		for e in star(v, e21); nodedata!(v, node(e)); end
		# two branches now have a bad orientation, fix this
		branch!(v, e23=>Branch(-1), e13=>Branch(-1))
		@assert iszero(noderadius(v,q1))
		@assert iszero(noderadius(v,q2))
	end
	return v
end#»»

# Voronoi diagram: geometry««1
# Geometric branch computation ««2
@inline function tripoint_ppp(v::VoronoiDiagram, c1,c2,c3)
	@assert ispoint(v,c1)
	@assert ispoint(v,c2)
	@assert ispoint(v,c3)
	return tripoint(point(v,c1), point(v,c2), point(v,c3))
end
@inline function tripoint_lpp(v::VoronoiDiagram, c1,c2,c3)
	@assert issegment(v,c1)
	@assert ispoint(v,c2)
	@assert ispoint(v,c3)
	return tripoint(segment(v, c1), point(v,c2), point(v,c3))
end
@inline function tripoint_llp(v::VoronoiDiagram, c1,c2,c3, s1)
	@assert issegment(v,c1)
	@assert issegment(v,c2)
	@assert ispoint(v,c3)
	return tripoint(segment(v,c1), segment(v,c2), point(v,c3))
end
@inline function tripoint_lll(v::VoronoiDiagram, c1,c2,c3)
	return tripoint(segment(v,c1), segment(v,c2), segment(v,c3))
end

@inline rot3l((r, a,b,c)) = (r, b,c,a) # rotate left
@inline rot3r((r, a,b,c)) = (r, c,a,b) # rotate right
function tripoint(v::VoronoiDiagram, q::Node)#««
	c1, c2, c3 = triangle(v, q)
	s1, s2, s3 = (separator(v, e) for e in sides(q))
	if issegment(v, c1)
		if issegment(v, c2)
			issegment(v, c3) && return       tripoint_lll(v, c1, c2, c3)
			                    return       tripoint_llp(v, c1, c2, c3, s1)
		else
			issegment(v, c3) && return rot3l(tripoint_llp(v, c3, c1, c2, s3))
			                    return       tripoint_lpp(v, c1, c2, c3)
		end
	else
		if issegment(v, c2)
			issegment(v, c3) && return rot3r(tripoint_llp(v, c2, c3, c1, s2))
			                    return rot3r(tripoint_lpp(v, c2, c3, c1))
		else
			issegment(v, c3) && return rot3l(tripoint_lpp(v, c3, c1, c2))
			                    return       tripoint_ppp(v, c1, c2, c3)
		end
	end
end#»»
# Edge updating ««2
@inline function edgedata!(v::VoronoiDiagram, e::Edge)
	o = opposite(v, e)
# 	println("  \e[34m edgedata!($e, $o)\e[m")
	# the separator is oriented with tail(o) = head(e) on its right,
	# i.e. pointing to the *left* of a:
	#       ╲n  head    ╱ 
	#  left  ╲____o____╱  right
	# +<⋯⋯   ╱    e    ╲  ⋯⋯>-
	#      p╱   tail    ╲
	#
	# branch[e] = +1 iff node(e) lies on the + branch of the separator
	# branch[e] = 0 iff this separator is a parallel bisector
	t, h, l, r = tail(v,e), tail(v,o), left(v,e), left(v,o)
	separators!(v, e, o, separator(v, h, t))
	# compute branches for e and o
end

# Node updating ««2
function nodedata!(v::VoronoiDiagram, q::Node)#««
	e1, e2, e3 = sides(q)
	c1, c2, c3 = tail(v,e1), tail(v,e2), tail(v,e3)
# 	println("\e[34;7m nodedata!($q = $c1,$c2,$c3):\e[m")
	r, b1, b2, b3 = tripoint(v, q)
	@assert (r ≥ 0) || (isnan(r))
	branch!(v, e1=>b1, e2=>b2, e3=>b3)
# 	# special case: two consecutive segments ««
# 	if issegment(v, c1) && issegment(v, c2) &&
# 		c3 ∈ cellsegment(v, c1) && c3 ∈ cellsegment(v, c2)
# 		v.noderadius[int(q)], v.geomnode[int(q)] = 0, point(v,c3)
# 	end
# 	if issegment(v, c2) && issegment(v, c3) &&
# 		c1 ∈ cellsegment(v, c2) && c1 ∈ cellsegment(v, c3)
# 		v.noderadius[int(q)], v.geomnode[int(q)] = 0, point(v,c1)
# 	end
# 	if issegment(v, c3) && issegment(v, c1) &&
# 		c2 ∈ cellsegment(v, c3) && c2 ∈ cellsegment(v, c1)
# 		v.noderadius[int(q)], v.geomnode[int(q)] = 0, point(v,c2)
# 	end
# 	#»»

	s1, s2, s3 = separator(v, e1), separator(v, e2), separator(v, e3)
	noderadius!(v, q=>r)
	p = evaluate(s1,b1,r)
# 	println("sep $c1/$c2: ", s1, b1, evaluate(s1,b1,r))
# 	println("sep $c2/$c3: ", s2, b2, evaluate(s2,b2,r))
# 	println("sep $c3/$c1: ", s3, (b3,r), evaluate(s3,b3,r))
	any(isnan, p) && return
	@assert evaluate(s1,b1,r) ≈ evaluate(s2,b2,r)
	@assert evaluate(s2,b2,r) ≈ evaluate(s3,b3,r)
	geometricnode!(v, q=>
		isstraight(s1) ? evaluate(s1, b1, r) :
		isstraight(s2) ? evaluate(s2, b2, r) :
		evaluate(s3, b3, r))
	@assert !isinf(geometricnode(v, q)[1])
# 	if int(q) == 11
# 	println("  tripoint($q = $c1,$c2,$c3) = $r, $b1, $b2, $b3")
# 	println(isstraight.((s1,s2,s3)))
# 	println("  e1 = $c1/$(head(v,e1))")
# 	println(s1)
# 	display((v,q))
# 	end
end#»»
# function CornerTables.swapnodes!(v::VoronoiDiagram, q1::Node, q2::Node)
# 	swapnodes!(CornerTables.triangulation(v), q1, q2)
# 	v.geomnode[SA[int(q1),int(q2)]] = v.geomnode[SA[int(q2),int(q1)]]
# 	v.noderadius[SA[int(q1),int(q2)]] = v.noderadius[SA[int(q2),int(q1)]]
# end

function separator(v::AbstractVoronoi, c1, c2)#««
	if issegment(v, c1)
		i1, j1 = cellsegment(v, c1)
		a1, b1 = point(v, i1), point(v, j1)
		if issegment(v, c2)
			return Separator(segment(v, c1), segment(v, c2))
		else # generic parabola separator
			return Separator(segment(v, c1), point(v,c2))
		end
	elseif issegment(v, c2)
		i2, j2 = cellsegment(v, c2)
			return Separator(point(v, c1), segment(v, c2))
	else # point-point separator
		return Separator(point(v, c1), point(v, c2))
	end
end#»»

# Offset ««1
# Single edge offset ««2
"""
    edgecross(v, e::Edge, r)

Given an edge bounding a cell and the r-offset region R,
returns the following booleans:
 - does the + branch of the separator cross the boundary ∂R?
 - does the - branch of the separator cross ∂R?
 - does the separator intersect R at all?
(The third boolean is *not* the disjunction of the two first ones:
namely, the separator may be fully included in the interior of R).
"""
@inline function edgecross(v::VoronoiDiagram, e::Edge, radius)#««
	o = opposite(v, e)
	q0, q1 = node(e), node(o)
	r0, r1 = noderadius(v, q0), noderadius(v, q1)
	# fi = does qi belong to offset region R?
	f0, f1 = (radius ≥ r0), (radius ≥ r1)
	f1 ≠ f0 && return (f1, f0, true)
	# we now know that q0, q1 are on the same side of ∂R
	sep = separator(v, e)
	b0, b1 = branch(v, e), branch(v, o)
# 	println("\e[32medgecross($e/$o, $radius):\e[m")
# 	println("    left node $q0 $b0 $(geometricnode(v, q0)) $r0 $b0 $f0")
# 	println("   right node $q1 $b1 $(geometricnode(v, q1)) $r1 $b1 $f1")

	if b0 == Branch(0) # this is a parallel bisector
		@assert b1 == Branch(0) "iszero(b0) -> iszero(b1)"
		@assert r0 == r1
		# TODO: find what to do with parallel bisectors
		return (false, false, f0 && (r1 > 0) && (r2 > 0))
	end
	@assert b1 ≠ Branch(0) "!iszero(b0) -> !iszero(b1)"
	# depending on (b0,b1), the edge is oriented this way:
	# ++ +- -+ --
	# <> << >> ><
	if b0 == Branch(+1) && b1 == Branch(+1) # non-monotonic edge:
		@assert r0 ≥ sep.rmin
		@assert r1 ≥ sep.rmin
		# case 1: the perigee of the separator lies outside R
		radius < sep.rmin && return (false, false, false)
		# case 2: the perigee lies inside, with zero or two nodes
		return (!f0, !f0, f0)
	end
	# this edge is monotonic with two nodes on the same side of ∂R:
	return (false, false, f0 && (r0 > 0) && (r1 > 0))
end#»»
"""    prevedge(v, e, r)
Given an edge bounding a cell c, return the previous edge where
the offset segment at distance r enters the cell."""
function prevedge(v::VoronoiDiagram, e0::Edge, r)#««
	@assert r ≥ 0
	for e in star(v, e0)
		(bl, _) = edgecross(v, e, r)
		bl && return e
	end
	return zero(e0)
end#»»
"""    nextedge(v, e, r)
Given an edge bounding a cell c, return the next edge where
the offset segment at distance r exits the cell."""
function nextedge(v::VoronoiDiagram, e0::Edge, r)#««
	@assert r ≥ 0
	for e in reverse(star(v, e0))
		(_, br) = edgecross(v, e, r)
		br && return e
	end
	return zero(e0)
end#»»
"""    firstedge(v, c, r)
Returns the first edge for an offset segment in cell `c`."""
@inline firstedge(v::VoronoiDiagram, c::Cell, r) = prevedge(v, anyedge(v, c), r)
"""    edgeinter(v, e, r)
Returns the status of 
"""

# Offset chain ««2
"""    offsetchains(v::VoronoiDiagram, radius, side)

Returns the set of chains encoding the offset curves for this radius.
Each chain is represented as a list of edges. Each edge correspond
to one cell traversed by the offset curve; it is the edge where
the curve enters the cell. The last edge in the chain represents either
the closure of the chain (if identical to the first) or the opposite edge to
the endpoint of the curve (if the chain is open).
"""
function offsetchains(v::VoronoiDiagram{J}, radius, side) where{J}#««
	# The last segment encodes the endpoint of the chain; it is either
	# identical to the first segment (if closed loop) or the opposite edge
	# of the true last point of the chain.
	#
	# At any point during this algorithm, the active chain is the last one.
	chains = Vector{Edge{J}}[]
	done = falses(nedges(v))
	@assert radius ≥ 0
	for startcell in sidesegments(v, side)
		c = Cell(startcell); e = firstedge(v, c, radius)
		iszero(e) && continue # this cell is not traversed by the offset curve
		done[int(e)] && continue # we already visited this curve segment

		# if this edge is not already done, then it lies on a new chain:
		l = [e]; push!(chains, l)
		while true
			e = opposite(v, nextedge(v, last(l), radius)); c = tail(v, e)
			push!(l, e); done[int(e)] = true
			!issegment(v, c) && v.neighbours[c] ≠ 2 && break
			e == first(l) && break
		end
		# if the chain is an open loop, we need to extend it to the left:
		first(l) == last(l) && continue
		while true
			e = prevedge(v, opposite(v, first(l)), radius); c = tail(v, e)
			!issegment(v, c) && v.neighbours[c] ≠ 2 && break
			pushfirst!(l, e); done[int(e)] = true
		end
	end
	return chains
end#»»

# Offset ««2
"""    interpolate(v, chain, radius, atol, start=1)

Interpolates an arc of ∂R as a polygonal path with absolute precision `atol`.
Returns (P = list of points, L = list of indices),
so that chain[i] corresponds to the points P[L[i]:L[i+1]].
"""
function interpolate(v::VoronoiDiagram{J,T}, chain, radius, atol) where{J,T}
	r = abs(radius)
	δ = √(r/(8*atol)) # used for discretizing circle arcs
	e0 = first(chain); sep0 = separator(v, e0)
	p0 = evaluate(sep0, Branch(+1), r)
	plist, llist = [p0], [1]
	for e1 in chain[2:end]
		sep1 = separator(v, e1)
		p1 = evaluate(sep1, Branch(+1), r)
		c = tail(v, e0)
		if !issegment(v, c) # circular arc
			p = point(v, c)
			a0, a01 = vectors_angles(p0-p, p1-p)
			n = ceil(Int, a01*δ); θ = a01/n
			for i in 1:n-1
				a = a0+i*θ
				push!(plist, SA[p[1]+cos(a)*r, p[2]+sin(a)*r])
			end
		end
		# this is either the single point (for a straight segment), or the
		# closure of a circular arc
		push!(plist, p1); push!(llist, length(plist))
		e0, sep0, p0 = e1, sep1, p1
	end
	(radius < 0) && reverse!(plist)
	return (plist, llist)
end
"""    offset(points, segments, radius; atol)

Returns the offset polygon(s) at distance `radius` from the polygons
defined by `points` and `segments`. Positive distance is right side.

Optional parameter `atol` is the maximal distance of consecutive
interpolated points on a circular arc.
"""
function offset(points, segments, radius::Real; atol=DEFAULT_ATOL)
	v = VoronoiDiagram(points, segments)
	r = abs(radius)
	chains = offsetchains(v, r, radius < 0)
	return [ interpolate(v, l, r, atol)[1] for l in chains ]
end
"""    offset(points, segments, radii; atol)

Equivalent to `[offset(points, segments, r) for r in radii]`,
but all offset paths are computed simultaneously.
"""
function offset(points, segments, radii::AbstractVector{<:Real};
		atol=DEFAULT_ATOL)
	v = VoronoiDiagram(points, segments)
	chains = [ offsetchains(v, abs(r), r < 0) for r in radii ]
	[[ interpolate(v, l, abs(r), atol)[1] for l in chains ] for r in radii ]
end
# Extrusion««1
# These functions compute a triangulation of the difference of two offset
# regions R(r2)∖R(r1), where r2 > r1.
# Mesh structure ««2
"""    Mesh

Points are doubly-indexed (by integers and points) to prevent repeats.
Faces are triangles of point indices."""
struct Mesh{J,P} <: AbstractVector{P}
	points::Vector{P}
	index::Dict{P,J}
	triangles::Vector{NTuple{3,J}}
end
@inline Mesh{J,P}() where{J,P} = Mesh{J,P}(P[], Dict{P,J}(), [])
@inline Base.size(mesh::Mesh) = (length(mesh.points),)
@inline Base.getindex(mesh::Mesh, i::Integer) = mesh.points[i]

function Base.push!(mesh::Mesh{J,P}, p) where{J,P} # returns index
	k = get(mesh.index, p, zero(J))
	!iszero(k) && return k
	push!(mesh.points, p)
	k = J(length(mesh))
	mesh.index[p] = k
	return k
end
@inline Base.append!(mesh::Mesh{J}, p) where{J} =
	J[ push!(mesh, x) for x in p ]

function Base.show(io::IO, mesh::Mesh)#««
	for (i, p) in pairs(mesh.points)
		println(io, " ",i, ": ", p)
		if mesh.index[p] ≠ i
			println(io, "\e[31;7m bad index[$p] = $(mesh.index[p]), should be $i")
			error("bad index")
		end
	end
end#»»

# Affine map ««2
"""    Affine3

An affine map mapping r1↦z1 and r2↦z2."""
struct Affine3{T}
	a::T
	b::T
	r1::T # we keep ri, zi to make this exact at these two points.
	r2::T
	z1::T
	z2::T
end
function Affine3((r1,z1)::Pair{T}, (r2,z2)::Pair{T}) where{T}
	a = (z2-z1)/(r2-r1)
	b = (z1*r2-r1*z2)/(r2-r1)
	return Affine3{T}(a, b, r1, r2, z1, z2)
end
# @inline Base.:-(aff::Affine3) =
# 	Affine3(-aff.a, aff.b, -aff.r1, -aff.r2, aff.z1, aff.z2)
function evaluate(aff::Affine3, r)
	# exact cases:
	r == aff.r1 && return aff.z1
	r == aff.r2 && return aff.z2
	return aff.a*r + aff.b
end

# Axial extrusion of a single point ««2
"""    AxialExtrude{J}

Represents the offset of a single point along the trajectory.

 - `radius`: the abscissa of this point
 - `chains`: the sequences of edge-crossings
 - `indices[i]`: the indices of points generated by each edge crossing
as points indexed by edge crossings."""
struct AxialExtrude{J,T}
	r::T
	z::T
	chains::Vector{Vector{Edge{J}}}
	indices::Dict{J,Vector{J}} # indexed by edges
# 	indices::Dict{J,NTuple{2,J}} # if the indices are always intervals
end
		
@inline npoints(a::AxialExtrude) = sum(length.(values(a.indices)))
@inline indices(a::AxialExtrude, e::Edge) = a.indices[int(e)]

function Base.show(io::IO, a::AxialExtrude)#««
	println(io, "axial offset of [", a.r, ",", a.z, "] has ",
		join(length.(a.chains),"+"), " crossings, $(npoints(a)) points:")
	for c in a.chains
		println(io, "  chain of $(length(c)) points: ")
		for e in c[1:end-1]
			println(io, "   edge $e -> points $(indices(a, e)) ->")
		end
		println(io, "   last edge $(last(c)) ",
			last(c) == first(c) ? "(closed)" : "(open)")
	end
end#»»
function AxialExtrude(v::VoronoiDiagram{J}, mesh::Mesh,#««
		rp::T, zp, side, atol) where{J,T}
	# `mesh`: `Mesh` collecting all mesh for this offset
	# `p`: the single point we are extruding
	indices = Dict{J,Vector{J}}()
	chains = offsetchains(v, rp, side)
	for ch in chains
		np = J(length(mesh))
		(newpoints, idx) = interpolate(v, ch, rp, atol)
		for i in 1:length(idx)-1
			j = idx[i]:idx[i+1]
			ind = append!(mesh, [[q; zp] for q in newpoints[j]])
			indices[int(ch[i])] = ind
		end
	end
	return AxialExtrude{J,T}(rp, zp, chains, indices)
end#»»

# Region between chains««2
"""    cell_band(v, e1in, e1out, c2next)
Returns a symbolic description of the contour of c ∩ (R(r2)∖R(r1)), where:
 - R(r1) enters c at `e1in` and exits at `e1out`;
 - c2next is the next-edge map for ∂R(r2).
The description is encoded as (edge, type, r), where `type` is:
1 for ∂R1 (in reverse direction), 2 for ∂R2,
3 for + branch of edge, 4 for - branch of edge;
`r` is the distance to site at starting point for this edge.

For example, the following cell band (delimited by the dashed lines
showing ∂R1 and ∂R2) is encoded as (e2,1) (e2,3) (e2,2) (e'4,4)
(assuming that e2, e4 are both the + branches of the separators).
     ________
    ↑   e3   ↑
   ┄│┄┄┄┄┄┄┄┄│┄┄┄┄ (∂R2)
 e'4│e4    e2│e'2
   ┄│┄┄┄┄┄┄┄┄│┄┄┄┄ (∂R1)
    │___e1___│
"""
function cell_band(v::VoronoiDiagram, e1in, e1out, c2next)
	elist, etype = [e1in], [Int8(1)]
	# orbit around the cell (fragment) c until we find the start point
	# for this edge
	e = e1in
	while true
		o = opposite(v, e)
		# is a decreasing part of edge `e` involved?
		int(branch(v, o)) > 0 && (push!(elist, o); push!(etype, Int8(4)))
		# is an increasing part of edge `e` involved?
		int(branch(v, e)) ≥ 0 && (push!(elist, e); push!(etype, Int8(3)))
		e == e1out && break
		c2n = c2next[int(e)]
		if !iszero(c2n)
			push!(elist, e); push!(etype, Int8(2))
			e = opposite(v, c2n)
		else
			e = after(v, e)
		end
	end
	return (elist, etype)
end
"""    chain_contour(v, chains1, chains2)

Given chains1 and chains2 enclosing R(r1) ⊂ R(r2),
returns all edges delimiting cells in R(r2)∖R(r1), as a list of tuples
(cell, edges, edgetyp)es), with:
 - edges = a (cyclic) list of edges around this cell,
 - edgetypes = list of types matching each edge, encoded as:
 1=∂R1, 2=∂R2, 3=positive edge, 4=negative edge.
"""
function chain_contour(v::VoronoiDiagram{J}, chains1, chains2) where{J}#««
	# build next-edge map for chains2 ««
	println("chains1 = $chains1; chains2=$chains2")
	c2next = zeros(Edge{J}, nedges(v))
	for ch in chains2
		e1 = first(ch)
		for e2 in ch[2:end]
			c2next[int(e1)] = e2
			e1 = e2
		end
	end #»»
	r = Tuple{Cell{J}, Vector{Edge{J}}, Vector{Int8}}[]
	for ch in chains1
		e1in = first(ch)
		for e1new in ch[2:end]
			c = tail(v, e1in)
			elist, etype = cell_band(v, e1in, opposite(v, e1new), c2next)
			println("  cell_contour($c; $e1in→opp($e1new)) is ($elist, $etype)")
			# ∂R1 enters the cell c at e1 and exits at o = opposite(e2)
			push!(r, (c, elist, etype))
			e1in = e1new
		end
		# if this chain is closed then the last examined cell did the
		# closing;
		# otherwise we still need to produce info about the “outer” boundary
		# (TODO)
	end
	return r
end#»»
"""    approximate(v, e, r1, r2, aff, atol)
Produces a list of points approximating the open interval [r1,r2] on edge e
(in this order), on the positive branch.
Returns list of [point; aff(r)] where r is distance to trajectory,
as well as the r parameters (rstart, rstop) delimitating this interpolation.
"""
function approximate(v::VoronoiDiagram{J,T}, e, r1, r2, aff, atol#««
		) where{J,T}
	o = opposite(v, e)
	qe, qo = node(e), node(o)
	ge, go = geometricnode(v, qe), geometricnode(v, qo)
	re, ro = noderadius(v, qe), noderadius(v, qo)
	be, bo = branch(v, e), branch(v, o)
	sep = separator(v, e)
	pt3 = @closure (p, r) -> SA[p[1], p[2], evaluate(aff, r)]
	isparallel(sep) && return [pt3(g, r) for (g,r) in ((ge,re), (go,ro))]

	# eliminate degenerate & parallel separators:
# 	ge == go && return [SA[ge[1], ge[2], evaluate(aff, re)]]
	re == ro && return [pt3(ge, re), pt3(go, re)]

	if int(be) ≥ 0
		r2 = min(r2, re)
		r1 = max(r1, int(bo) ≥ 0 ? sep.rmin : ro)
	else
		@assert int(bo) > 0
		r2 = min(r2, ro)
		r1 = max(r1, re)
	end
	r1 > r2 && return [pt3(ge, re)]

	rlist = approximate(sep, r1, r2, atol)

	return [pt3(evaluate(sep, Branch(+1), r), r) for r in rlist]
end#»»
"""    edgepoints(v, points, edge, edgetype, tlist, ax1, ax2)

Interpolates points on this edge, according to edgetype,
counter-clockwise around the cell. To avoid repeated points around the
cell, the common point between two edges is counted as belonging to the
firt of those two edges.

The points are pushed on the point-set and the corresponding indices are
returned.
"""
function edgepoints(v::VoronoiDiagram, mesh::Mesh,
		edge, edgetype, ax1, ax2, aff, atol)#««
# 	println("\e[34mdraw ($edge, $edgetype); r1=$(ax1.r) r2=$(ax2.r)\e[m")
	if edgetype == 1 # use segment from ∂R1, backwards
		return reverse(indices(ax1, edge)[begin+1:end])
	elseif edgetype == 2 # segment from ∂R2, forwards
		return indices(ax2, edge)[begin:end-1]
	end # interpolate along an edge:
	# (increasing if edgetype == 4, decreasing if edgetype == 3)
	seg = approximate(v, edge, ax1.r, ax2.r, aff, atol)
	edgetype == 4 && reverse!(seg)
	pop!(seg)
	return append!(mesh, seg)
end#»»

# Extrusion of a polygonal loop««2
function extrude_vertical(v::VoronoiDiagram, mesh, p, q, atol)
	r, side, zp, zq = abs(p[1]), p[1] < 0, p[2], q[2]
	println("\e[48;5;88mextrude a vertical face at r=$r: z=($zp, $zq)\e[m")
	chains = offsetchains(v, r, side)
	for ch in chains
		l = interpolate(v, ch, r, atol)[1]
		a = first(l);
		pa, qa = append!(mesh, (SA[a...; zp], SA[a...; zq]))
		for b in l[2:end]
			pb, qb = append!(mesh, (SA[b...; zp], SA[b...; zq]))
			push!(mesh.triangles, (pa,pb,qa), (pb,qb,qa))
			a, pa, qa = b, pb, qb
		end
	end
end
function extrude_sloped(v::VoronoiDiagram{J}, mesh, p, q, atol) where{J}
	println("\e[7mextrude a sloped face: $p->$q\e[m")
	side = p[1] < 0 || q[1] < 0
	rp, rq = abs(p[1]), abs(q[1])
	(r1,r2,z1,z2,rev) = rp < rq ? (rp,rq,p[2],q[2],side) : (rq,rp,q[2],p[2],!side)
	aff = Affine3(r1=>z1, r2=>z2)
	ax1 = AxialExtrude(v, mesh, r1, z1, side, atol)
	ax2 = AxialExtrude(v, mesh, r2, z2, side, atol)
	for (c, elist, tlist) in chain_contour(v, ax1.chains, ax2.chains)
		cellpoints = Int64[] # this will be passed to LibTriangle
		println("\e[35;7min cell $c: $elist, $tlist\e[m")
		for (e, t) in zip(elist, tlist) #(edge, edgetype)
			ep = edgepoints(v, mesh, e, t, ax1, ax2, aff, atol)
			println("  \e[1m($e,$t) contributes $ep=$(mesh[ep])\e[m")
			append!(cellpoints, ep)
		end
		unique!(cellpoints)
		if length(cellpoints) ≥ 3
			tri = triangulate_loop(mesh, cellpoints)
			append!(mesh.triangles, (rev ? (a,b,c) : (a,c,b) for (a,b,c) in tri))
		end
	end
end
"""    extrude_loop(v, loop)

Extrudes a loop of points [xi, yi] along the polygonal path(s);
returns (points, triangulation).
"""
function extrude_loop(v::VoronoiDiagram{J,T}, loop, atol) where{J,T}
	# axial paths: extrusions of individual points of the loop««
	mesh = Mesh{J,SVector{3,T}}()
	p = last(loop)
	for q in loop
		if p[1]*q[1] < 0
			s = SA[zero(p[1]), (p[1]*q[2]-p[2]*q[1])/(p[1]-q[1])]
			extrude_sloped(v, mesh, p, s, atol)
			extrude_sloped(v, mesh, s, q, atol)
		elseif p[1] == q[1]
			extrude_vertical(v, mesh, p, q, atol)
		else
			extrude_sloped(v, mesh, p, q, atol)
		end
		p = q
	end
	return (mesh.points, mesh.triangles)
end
"""    extrude(trajectory, profile, atol)

 - `trajectory`: vector of paths, either open or closed
 - `profile`: vector of open loops

Returns a vector of extrusions of each profile loop along the trajectory.
Each extrusion is a (points, triangles) pair.
"""
function extrude(trajectory, profile, atol)
	# decompose trajectory to (points, segments):
	plist = empty(first(trajectory))#««
	slist = NTuple{2,Int}[]
	for path in trajectory
		closed = last(path) == first(path)
		n = length(plist)
		append!(plist, path[begin:end-closed])
		for i in 1:length(path)-1
			push!(slist, (n+i, n+i+1))
		end
		closed && push!(slist, (n+length(path), n+1))
	end#»»
	extra = maximum(maximum(p[1] for p in loop) for loop in profile)
	println("plist=$plist\nslist=$slist\nextra=$extra\n")
	v = VoronoiDiagram(plist, slist; extra)
	return [ extrude_loop(v, loop, atol) for loop in profile ]
end
#»»1
# Displaying and debugging ««1
function gnuplot(io::IO, v::AbstractVoronoi; scale=10.)
# 	f1 = identity
	f1 = @closure x -> scale*sign(x)*log(1+abs(x)/scale)
	pt = @closure c -> f1.(point(v,c))
	# index 0: points (x y label)
	# index 1: segments (x y dx dy label)
	println(io, "# index 0: points (x y label)")
	for c in 1:npoints(v); c = Cell(c)
		g = pt(c)
		println(io, g[1], "\t", g[2], "\t", c)
	end
	println(io, "\n\n# index 1: segments (x y dx dy lx ly label)")
	for c in 1:nsegments(v); c = Cell(c+npoints(v))
		(c1, c2) = cellsegment(v, c)
		(g1, g2) = pt(c1), pt(c2)
		dg = g2 - g1; u = g1 + .5*dg - .1*quarterturn(dg)
		println(io, g1[1], "\t", g1[2], "\t", dg[1], "\t", dg[2], "\t",
			u[1], "\t", u[2], "\t", c, "\t# = ", (c1, c2))
	end
	nodepos = Vector{SVector{2,Float64}}(undef, nnodes(v))
	cellpos = Vector{SVector{2,Float64}}(undef, ncells(v))
	for q in eachnode(v)
		isone(int(q)) && continue
		g = geometricnode(v, q)
		d = zero(g)
		for i in 1:3
			c = tail(v, side(q, i))
			if issegment(v, c)
				(c1, c2) = cellsegment(v, c)
				d += .5*(pt(c1)+pt(c2))
			else
				d += pt(c)
			end
		end
		nodepos[int(q)] = f1.(g + .3*(d - g))
	end
	for c in eachcell(v)
		iszero(anyedge(v, c)) && continue
		n = 0; g = zero(first(nodepos))
		for e in star(v, c)
			n+= 1; g += nodepos[int(node(e))]
		end
		cellpos[int(c)] = g/n
	end
	println(io, "\n\n# index 2: nodes (x y label)")
	for q in eachnode(v)
		g = nodepos[int(q)]
		println(io, g[1], "\t", g[2], "\t", q, "\t# ", triangle(v,q))
	end
	println(io, "\n\n# index 3: edges (x y dx dy x1 y1 label1 x2 y2 label2)")
	for e in eachedge(v)
		int(e) ≤ 3 && continue
		o = opposite(v, e); o < e && continue
		q1, q2 = node(e), node(o)
		g1, g2 = nodepos[int(q1)], nodepos[int(q2)]
		dg = g2 - g1
		l1 = g1 + .5*dg + .1*f1.(quarterturn(dg))
		l2 = g1 + .5*dg - .1*f1.(quarterturn(dg))
		println(io, g1[1], "\t", g1[2], "\t", dg[1], "\t", dg[2],
			"\t", l1[1], "\t", l1[2], "\t", e,
			"\t", l2[1], "\t", l2[2], "\t", o,
			"\t# ", e, o, q1, q2)
	end
	println(io, "\n\n# index 4: cells (x y label)")
	for c in eachcell(v)
		g = cellpos[int(c)]
		println(io, g[1], "\t", g[2], "\t", c)
	end
end
function gnuplot(v::AbstractVoronoi; scale=.8, f_png="/tmp/a.png")
	f_dat = "/tmp/a.dat"
	open(f_dat, "w") do io gnuplot(io, v; scale); end
	f_gpi = "/tmp/a.gpi"
	open(f_gpi, "w") do io
		println(io, """
set style textbox 1 opaque border lc "blue"
set style textbox 2 opaque border lc "red"
f='$f_dat'
set terminal png fontscale .5 size 1000,800
set output '$f_png'
plot \\
  f index 1 u 1:2:3:4 w vectors lc "blue" lw 2, \\
	f index 3 u 1:2:3:4 w vectors lc "red", \\
	f index 4 u 1:2:3 w labels center textcolor "forest-green", \\
	f index 3 u 5:6:7 w labels center textcolor "red", \\
	f index 3 u 8:9:10 w labels center textcolor "red", \\
	f index 2 u 1:2:3 w labels center boxed bs 2 textcolor "red", \\
  f index 0 u 1:2:3 w labels center boxed bs 1
		""")
	end
	run(`gnuplot $f_gpi`)
end
function Base.show(io::IO, ::MIME"text/plain",
	(v,q)::Tuple{AbstractVoronoi,Node})
	bc = @closure i->"X-0+!X"[clamp(3+int(branch(v,side(q,i))), 1:6)]
	print(io, "\e[33m", q, triangle(v,q), "\e[m: ",
		@sprintf("(%.3g,%.3g) r=%.3g", geometricnode(v, q)..., noderadius(v, q)),
		" ", bc(1), bc(2), bc(3),
	)
# 	for i in (1,2,3)
# 		e = side(q, i); o = opposite(v, e); oo = opposite(v,o)
# 		oo ≠ e && println(io, "  \e[31;7m opposite($o) = $oo, should be $e\e[m")
# 	end
end
# »»1
end

V=Voronoi
using FastClosures

using StaticArrays
l = [[-2,0.],[1,0],[-1,2],[0,2],[0,3],[-1,3],[0,4],[-1,4],[1,5],[-2,7]]
l1 = [[-1,0.],[1,0],[1,1],[-1,1]]
l2 = [[1,1.],[-1,1],[-1,0],[1,0]]
V.split_loop(l1)
V.split_loop(l2)
# # TODO: tests
# # V = Voronoi
# # t=V.triangulate([[-10,0],[10,0],[0,10.]])
# # v=V.voronoi([(0,0),(10,0),(11,3),(6,2),(5,9),(1,10)],[(6,1),(1,2),(2,3),(3,4),(4,5),(5,6)])
# # v=V.VoronoiDiagram([[0.,0],[10,0],[0,10],[10,10],[5,9],[5,1]],[(3,4),(5,6)])
# 
# # v=V.OffsetDiagram([[0.,0],[10,0]],[])
# # println(V.point(v,V.Cell(3)))
# # println(V.point(v,V.Cell(4)))
# # println(V.point(v,V.Cell(5)))

# HERE:
s1a = (SA[-5,0.],SA[-3,0.])
s1b = (SA[-2,0.],SA[2,0.])
s1c = (SA[3,0.],SA[5,0.])
us2 = SA[3,4.]
s2a, s2b, s2c = (-2us2, -us2), (-us2, us2), (us2, 2us2)
s3a = (SA[-5,2.],SA[-2,2.])
s3b = (SA[-2,2.],SA[2,2.])
s3c = (SA[2,2.],SA[5,2.])
s4a = (SA[6,-3.],SA[6,-1.])
s4b = (SA[6,1.],SA[6,2.])

c1 = SA[0,0.]
c2 = SA[10,0.]
c3 = SA[5,1.]
c4 = SA[5,9.]
c8 = (c1, c2)
c9 = (c2, c3)
c10 = (c3, c4)



v = V.VoronoiDiagram([[0.,0],[10,0]],[(1,2)])
l = [[0,0.],[3,0],[0,4]]
v=V.VoronoiDiagram([[0.,0],[10,0],[5,1],[5,9]],[(1,2),(2,3),(3,4)];extra=0)
# v=V.VoronoiDiagram([[0.,0],[10.,0],[10,10.]],[(1,2),(2,3)];extra=5)

el = V.extrude_loop(v, [[-.4,-.6],[.4,-.44],[.4,.6],[-.6,.4]], .1)
# el = V.extrude_loop(v, [[.5,0],[2,0],[2,1]], .1)
# el = V.extrude_loop(v, [[-.5,-1],[1,-.5],[.5,1],[-1,.5]], .1)

# v=V.OffsetDiagram([[0.,0],[10,0],[0,10],[10,10],[5,9],[5,1]],[(1,2),(2,6),(6,5),(5,4),(3,1)])
# l=V.offsetchains(v, 1., false)
# o=V.offset([[0.,0],[10,0],[0,10],[10,10],[5,9],[5,1]],[(1,2),(2,6),(6,5),(5,4),(3,1)], 1.)
# ci = V.chain_contour(v, V.offsetchains(v, 1., false), V.offsetchains(v, 10., false))
# el=V.extrude_loop(v, [SA[1.,1],SA[-1.,0],SA[1.,-1]], .1)
