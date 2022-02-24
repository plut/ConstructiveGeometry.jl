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
# - remove Line:
#   - tripoint(SSP)
#   - tripoint(SSS)
# - use separator orientation to determine edge capture (remove `capture`)
#
"""    Voronoi

Computation of Voronoi diagrams for planar polygons and polygonal paths,
and of offset paths using these diagrams."""
module Voronoi
using StaticArrays
using FastClosures
using LinearAlgebra
using LazyArrays
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

@inline det2(u,v) = u[1]*v[2]-u[2]*v[1]
@inline det2(u,v,w) = det2(v-u, w-u)
@inline norm²(v) = v[1]^2+v[2]^2
@inline unit(v) = v/√(norm²(v))
@inline distance²(a,b) = norm²(a-b)
@inline quarterturn(v) = SA[-v[2], v[1]]

const Segment{T} = NTuple{2,<:AbstractVector{T}}

function lineinter(a,b,c,d)
	D = det2(a-b, c-d)
	t = det2(a-c, c-d)
	z = a+(t/D)*(b-a)
	return z
end
function linedistance²(a,b,c)
	# distance of c from segment (a,b)
	ab,ac = b-a,c-a
# 	d = dot(b,c); n = norm²(b)
# 	d < 0 && return norm²(c)
# 	d > n && return norm²(c-b)
	return det2(ab,ac)^2/norm²(ab)
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

# Lines ««2
Point{T} = SVector{2,T}
"the equation for a (normalized, oriented) straight line in the plane."
struct Line{T}
	# orientation: the normal vector points to the *left* of the line
	normal::SVector{2,T} # normalized to ‖u‖=1
	offset::T # line equation is normal⋅z + offset == 0
end
@inline Base.:*(a::Real, l::Line) = Line(a*l.normal, a*l.offset)
@inline normalize(l::Line) = (1/√(norm²(l.normal)))*l
@inline Base.:-(l::Line) = (-1)*l

@inline direction(l::Line) = quarterturn(l.normal)
# one arbitrary point on the line (actually the projection of the origin)
# FIXME: since this line is normalized, we should not need to divide by norm
@inline point(l::Line) = -l.normal*l.offset/(l.normal[1]^2+l.normal[2]^2)

@inline Line(a::AbstractVector, b::AbstractVector) = # oriented from a to b
	normalize(Line(SA[a[2]-b[2], b[1]-a[1]], a[1]*b[2]-a[2]*b[1]))

# signed distance from l to a; positive sign corresponds to left side of l.
@inline evaluate(l::Line, a::AbstractVector) = dot(l.normal, a) + l.offset

"returns either the intersection point, or `nothing`."
function Base.intersect(l1::Line, l2::Line)
	(a1, b1), c1 = l1.normal, l1.offset
	(a2, b2), c2 = l2.normal, l2.offset
	d = a1*b2-a2*b1
	iszero(d) && return nothing
	return SA[b1*c2-c1*b2,c1*a2-a1*c2]/d
end

"    det2(line1, line2): determinant of direction of these lines."
@inline det2(l1::Line, l2::Line) = det2(l1.normal, l2.normal)
"    dot(line1, line2): cosine of the angle formed by these lines."
@inline LinearAlgebra.dot(l1::Line, l2::Line) = dot(l1.normal, l2.normal)
"""    det2(line1, line2, line3):
Returns (twice) the oriented area of the triangle bounded by these
three lines (i.e. the intersection of the three half-planes).
This independent of the order of the lines,
but depends on the orientation of each line."""
function det2(l1::Line, l2::Line, l3::Line)#««
	D = det(SA[
		l1.normal[1] l1.normal[2] l1.offset
		l2.normal[1] l2.normal[2] l2.offset
		l3.normal[1] l3.normal[2] l3.offset
		])
	d1, d2, d3 = det2(l2, l3), det2(l3,l1), det2(l1, l2)
	return (l1.offset+l2.offset+l3.offset)*D/(d1*d2*d3)
end#»»

# # Circumscribed circle ««2
# function circumcenter(a,b,c)
# 	ab,ac = b-a,c-a
# 	m = SA[norm²(ab) ab[1] ab[2];norm²(ac) ac[1] ac[2]]
# 	kn = det(m[:,SA[2,3]])
# 	kx = det(m[:,SA[1,3]])/(2kn)
# 	ky = det(m[:,SA[2,1]])/(2kn)
# 	return a + SA[kx, ky]
# end
# 
# # """    circumcenter_orientation(a,b,c)
# # 
# # Returns +1 iff a views the circumcenter to the left of b."""
# # function circumcenter_orientation(a,b,c)
# # 	ab,ac = b-a,c-a
# # 	d = det2(ab, ac)
# # 	@assert !iszero(d)
# # 	# the power of c relative to the circle with diameter (ab) is
# # 	# cx² + cy² -(ax+bx)cx -(ay+by)cy + (ax bx + ay by)
# # 	u = c[1]^2+c[2]^2 - (a[1]+b[1])*c[1] - (a[2]+b[2])*c[2] + (a[1]*b[1]+a[2]*b[2])
# # 	iszero(u) && return Int8(0)
# # 	u > 0 && return (d > 0 ? Int8(+1) : Int8(-1))
# # 	return (d > 0 ? Int8(-1) : Int8(+1))
# # end
# # 	
# # 
# # function circumradius(a,b,c)
# # 	ab,ac,bc = b-a, c-a, c-b
# # 	return sqrt(norm²(ab)*norm²(ac)*norm²(bc))/(2*abs(det2(ab,ac)))
# # end
# 
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
# 
# """
#     isincircle(a,b,c,p,q)
# 
# Returns true iff open segment ]p,q[ intersects circumcircle of triangle (a,b,c).
# """
# function isincircle(a,b,c,p,q)
# 	a,b,c,q = a-p,b-p,c-p,q-p
# 	na, nb, nc = norm²(a), norm²(b), norm²(c)
# 	# equation of circumcircle is kn(x²+y²) - kx x - y y + k0 = 0
# 	m = SA[na a[1] a[2] 1;nb b[1] b[2] 1;nc c[1] c[2] 1]
# 	kn = det(m[:,SA[2,3,4]])
# 	kx = det(m[:,SA[1,3,4]])/kn
# 	ky = det(m[:,SA[1,4,2]])/kn
# 	C = det(m[:,SA[3,2,1]])/kn # equation for (t*q) is At² - Bt + C = 0
# 	A = norm²(q)
# 	B = kx*q[1] + ky*q[2]
# 	return (B ≥ 0) && (B ≤ 2A) && (B*B ≥ 4*A*C)
# end
# 
# function incenter(a,b,c)
# 	la, lb, lc = distance²(b,c), distance²(c,a), distance²(a,b)
# 	p = la+lb+lc
# 	ra, rb, rc = la/p, lb/p, lc/p
# 	return ra*a + rb*b + rc*c
# end
# 
# # Equidistant points ««2
# # p = point, s = segment, x = segment starting on previous point
# 
# @inline function equidistant_pps(a, b, p, q)#««
# 	# returns the point equidistant from a, b, (pq)
# 	# chosen so that the cells are oriented (a, b, pq)
# 	# imagine our segment as an x-axis; both points are on the same side,
# 	# we reorient the segment so they both have y > 0
# 	pqa = det2(p,q,a)
# 	@assert !iszero(pqa)
# 	(pqa < 0) && ((p,q) = (q,p); pqa = -pqa)
# 	pqb = det2(p,q,b)
# 	@assert pqb > 0
# 	ab2 = distance²(a,b)
# 	pq2 = distance²(p,q)
# 	pqab = dot(q-p, b-a)
# 	if pqa == pqb
# 		# special case if (a,b) and (p,q) are collinear:
# 		@assert false
# 	end
# 	# let z = (a+b)/2 + t*I*(b-a); then
# 	# d²(z,a) = d²(z,b) = ab2(t^2 + 1/4)
# 	# d²(z,pq) = <pqz>²/pq2 = (pqa/2+pqb/2+t*pqab)^2/pq2, so the eq. is:
# 	# (using the identity (pq2*ab2 = pqab^2 + (pqa-pqb)^2):
# 	# (pqa-pqb)^2 * t^2 - (pqa+pqb) pqab t + (pqab^2/4 - pqa*pqb) = 0
# 	#
# 	# We find Δ = 4*pqa*pqb*ab2*pq2
# 	# and geometry implies that the +√ is the correct sign
# 	Δ = 4*pqa*pqb*ab2*pq2
# 	t = ((pqa+pqb)*pqab + √(Δ)) / (2*(pqa-pqb)^2)
# 	return (a+b)/2 + t*quarterturn(b-a)
# end#»»
# 
# @inline function equidistant_pxs(a,b,p,q, ε)
# 	# return point equidistant from: a, (ab), (pq)
# 	# (in this order if ε=+1, opposite order if ε=-1)
# 	# z = a + tI(b-a) satisfies
# 	# d²(z,a) = d²(z,ab) = ab² t²
# 	# d²(z,pq) = <pqz>²/pq² = (<pqa> + t pq.ab)^2 / pq2
# 	# or: (pqa-pqb)^2*t^2 - 2 pqa*pqab*x - pqa^2 = 0
# 	# Δ = 4 pqa^2*pq2*ab2
# 	pqa = det2(p,q,a)
# 	pqb = det2(p,q,b)
# 	ab2 = distance²(a,b)
# 	if pqa == pqb # special case: both lines are parallel
# 		abp = det2(a,b,p)
# 		return a + abp/(2*ab2)*quarterturn(b-a)
# 	end
# 	pq2 = distance²(p,q)
# 	pqab = dot(q-p, b-a)
# 	Δ = 4*pqa^2*pq2*ab2
# 	t = (2*pqa*pqab+ ε*√(Δ)) / (2*(pqa-pqb)^2)
# 	z = a + t*quarterturn(b-a)
# 	return z
# end
# 
# "returns the point equidistant from point a and segments (pq), (rs)"
# function equidistant_pss(a, p, q, r, s)
# 	pqrs = det2(q-p, s-r)
# 	if iszero(pqrs) # special case: parallel lines
# 		pqa = det2(p,q,a)
# 		pqa < 0 && ((p,q, pqa) = (q,p, -pqa)) # ensure orientation
# 		pq = q-p
# 		pq2 = norm²(pq)
# 	pqr, pqs = det2(p,q,r), det2(p,q,s)
# # 	# `a` is not allowed to be on any of the two lines pq, rs
# # 	if pqr == pqs # special case: parallel lines
# 		# let v = quarterturn(q-p)
# 		λ = (pqr - 2*pqa)/(2*pq2)
# 		# c = a + λv is the projection of a on the middle line
# 		# this line is parametrized as z = c + t pq, then
# 		# az² = ac²+t² pq² = (λ²+t²) pq² must be pqr^2/(4 pq^2), hence
# 		# t^2 = (pqr^2 - (pqr-2pqa)^2)/(4 pq^4)
# 		#     = pqa(pqr-pqa)/pq^4
# 		# Geometry imposes the positive square root
# 		t = √(pqa*(pqr-pqa)) / pq2
# 		z = a+λ*quarterturn(pq)+t*pq
# 		return a + λ*quarterturn(pq) + t*pq
# 	end
# 	c = lineinter(p,q,r,s)
# 	ca = a-c
# 	pq = q-p; pq2 = norm²(pq); upq = pq/√(pq2)
# 	rs = s-r; urs = rs/√(norm²(rs))
# 	dot(upq, ca) < 0 && (upq = -upq)
# 	dot(urs, ca) < 0 && (urs = -urs)
# 	# parametrization of the inner angle bisector: z = c + t u
# 	c = lineinter(p, q, r, s)
# 	u = urs + upq
# 	# for z = c+tu: d²(z,pq) = t² pqu²/pq², while
# 	# d²(z,a) = ‖a-c-tu‖² = t² u²- 2t ca⋅u + ca²
# 	# the equation is thus t²(u²-pqu²/pq²) -2t ca⋅u + ca² = 0, or
# 	# t²(pq⋅u)²/pq² -2t ca⋅u + ca² = 0
# 	A = dot(pq, u)^2 / pq2
# 	B = dot(ca, u)
# 	C = norm²(ca)
# 	Δ = B^2-A*C
# 	ε = sign(det2(upq, urs))
# 	t = (B + ε*√(Δ))/A
# 	z = c+t*u
# 	return c + t*u
# end
# 
# function equidistant_sss(a1,b1,a2,b2,a3,b3)
# 	# returns the point equidistant from the *oriented* lines (ai, bi)
# 	# (i.e. either the incenter, or an excenter, of the triangle, according
# 	# to the orientations).
# 	u1 = quarterturn(b1 - a1); u1 /= √(norm²(u1))
# 	u2 = quarterturn(b2 - a2); u2 /= √(norm²(u2))
# 	u3 = quarterturn(b3 - a3); u3 /= √(norm²(u3))
# 	p1 = lineinter(a2, b2, a3, b3)
# 	p2 = lineinter(a3, b3, a1, b1)
# 	return lineinter(p1, p1+u2+u3, p2, p2+u1+u3)
# end
# 
# function equidistant_sss_parallel(a, u, b, p, pq)
# 	# returns the point equidistant from lines (a, a+u), (b,b+u), (p, pq)
# 	# (in this order)
# 	# parametrized as z = c+tu (with c = (a+b)/2)
# 	c, ab = (a+b)/2, b-a
# 	# d²(z, (a,a+u)) = l²/4 where l is the distance between the parallel lines
# 	# so (<pqc> + t <pqu>)² = |ab|⋅|pq|/2, or t=-<pqc>±|ab|.|pq|/(2<pqu>)
# 	# if <u,ab> > 0 then a lies to the right of the line (c, c+u)
# 	# and we take the + sign:
# 	pqc, pqu = det2(pq, c-p), det2(pq, u)
# 	l2, pq2 = det2(u,ab)^2/norm²(u), norm²(pq)
# 	t = (-pqc + sign(det2(u,ab))*sqrt(l2*pq2)/2)/pqu
# 	z = c + t*u
# 	return z
# end
# 
# # Bisectors ««2
# "mediator line of segment (ab)"
# function mediator(a,b)
# 	return SA[2*(a[1]-b[1]), 2*(a[2]-b[2]), norm²(b)-norm²(a)]
# end
# 
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
	vmat = [ points[i][j] for i in idx, j in 1:2 ]
	elist = [ idx[mod1(i+j,n)] for j in 1:n, i in 0:1]
	println("vmat=",vmat)
	println("idx=",idx)
	println("elist=",elist)
	return LibTriangle.constrained_triangulation(vmat, idx, elist)
end

# Minimizing a quadratic function
"""    min_quadratic((a,b,c), (x1,x2))
Returns the minimal value of (a x²+2bx+c) on the interval [x1,x2]."""
function min_quadratic((a,b,c), (x1,x2))
	x1,x2 = minmax(x1,x2)
	iszero(a) && return min(b*x1,b*x2)+c
	@assert a > 0
	xc = -b/a
	x1 < xc < x2 && return c-b^2/a
	return min(x1*(b+a*x1), x2*(b+a*x2))+c
end
"""    min_quartic(f, (x1,x2))
Returns the minimal value of the quartic f on the interval [x1,x2].
f is given as a list of coefficients (a0,a1,a2,a3,a4)."""
function min_quartic(f, (x1,x2))
	x1,x2 = minmax(x1,x2)
	r0 = (x1+x2)/2
	for _ in 1:5
		r1 = ((8*f[1]*r0+3*f[2])*r0^2-f[4])/(2*f[3]+r0*(6*f[2]+12*f[1]*r0))
		abs(r1 - r0) ≤ 1e-6*abs(r0) && break
		r0 = r1
	end
	r0 = clamp(r0, x1, x2)
	return f[5]+r0*(f[4]+r0*(f[3]+r0*(f[2]+r0*f[1])))
end
# Separators (parametrized bisectors) ««1
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

# predicates
@inline isparallel(sep::Separator) = any(isnan, sep.normal)
@inline isstraight(sep::Separator) = iszero(sep.normal)
@inline ishalflines(sep::Separator)= iszero(sep.rmin) && !iszero(sep.normal)
# @inline isparabola(sep::Separator) = !iszero(sep.normal) # default case

@inline Base.show(io::IO, sep::Separator) =
	@printf(io, "sep %s(o=[%.3g,%.3g], r₀=%.3g, t=[%.3g,%.3g], n=[%.3g,%.3g])",
		isparallel(sep) ? "═" :
		isstraight(sep) ? "─" :
		ishalflines(sep) ? "⋁" : "◡",
		sep.origin..., sep.rmin, sep.tangent..., sep.normal...)

"""    reverse(separator)

Given `separator(a,b)`, returns `separator(b,a)`, i.e. such that
`evaluate(sep′,b,s) = evaluate(sep,-b,s)."""
@inline Base.reverse(s::Separator)= ishalflines(s) ?
	Separator(s.origin, s.normal, s.tangent, s.rmin) :
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
# 		x2 ≤ 0 && return Separator((p1+p2)/2, k/f*v, zero(p1q1), -x2/(2*f))
# 		x2 ≥ x1 && return Separator((q1+p2)/2, -k/f*v, zero(p1q1), (x2-x1)/(2*f))
		throw(PointInSegment())
	end
	rmin = y2/(2*f)
	return Separator(p2 - y2/2*v/x1,
		k*sign(y2)*√(2*abs(y2)/f)*p1q1/f,
		sign(y2)*v/f, abs(rmin))
end
@inline Separator(a::AbstractVector, b::Segment) = Separator(b, a; k=-1) #»»
function Separator((p1,q1)::Segment, (p2,q2)::Segment)#««
	p1q1, p2q2 = q1-p1, q2-p2
	d = det2(p1q1, p2q2)
	if iszero(d)
		l = det2(p1q1, p2-p1)
		(l < 0) && ((p1, q1, p1q1, l) = (q1, p1, -p1q1, -l))
		u = unit(p1q1); l = det2(u, p2-p1)
		return Separator((p1+p2)/2, u, SA[oftype(l, NaN), oftype(l, NaN)], l/2)
	end
	# both segments are un-ordered, so we swap p2, q2 if needed:
	d < 0 && ((p2q2, p2, q2, d) = (-p2q2, q2, p2, -d))
	p1p2, p1q2 = p2-p1, q2-p1
# 	println("p1p2, p1q2, p2q2, d = ", (p1p2,p1q2,p2q2,d))
	c = lineinter(p1, q1, p2, q2)
	Dp2, Dq2, Dp1, Dq1 =
		det2(p1q1, p1p2), det2(p1q1, p1q2), det2(p1p2, p2q2), det2(p2q2, q1-p2)
# 	println("seg2 is ",
# 		Dp2 ≥ 0 ? "left of" : Dq2 ≤ 0 ? "right of" : "across", " line1; ",
# 	" seg1 is ",
# 		Dq1 ≥ 0 ? "left of" : Dp1 ≤ 0 ? "right of" : "across", " line2")
	@assert Dq2 ≈ Dp2+d
	@assert Dq1 ≈ Dp1-d
	@assert Dp2 < Dq2
	@assert Dp1 > Dq1
	u1, u2, z = √(norm²(p2q2))*p1q1/d, √(norm²(p1q1))*p2q2/d, zero(d)

	if Dp2 ≥ 0 # seg2 is left of line1
		Dq1 ≥ 0 && return Separator(c,  u1-u2, zero(u1), z)
		Dp1 ≤ 0 && return Separator(c,  u1+u2, zero(u1), z)
		return Separator(c, u1+u2, -u1+u2, z)
	elseif Dq2 ≤ 0 # seg2 is right of line1
		Dq1 ≥ 0 && return Separator(c, -u1-u2, zero(u1), z)
		Dp1 ≤ 0 && return Separator(c, -u1+u2, zero(u1), z)
		return Separator(c, -u1-u2, +u1-u2, z)
	else # seg2 is across line1
		Dq1 ≥ 0 && return Separator(c, -u1-u2, -u1+u2, z)
		Dp1 ≤ 0 && return Separator(c,  u1+u2,  u1-u2, z)
		throw(CrossingSegments())
	end
end#»»
# Evaluation, interpolation««2

"""    evaluate(separator, r, sign)

Returns the point on the separator situated at distance `r` from both
sites and on the branch given by sign `s` (either + or -).
The `+` branch sees `a` on its right and `b` on its left.
"""
@inline function evaluate(sep::Separator, b, r) # b is a sign
	ishalflines(sep) &&
		return sep.origin + r * (b ≥ 0 ? sep.tangent : sep.normal)
	isstraight(sep) &&
		return sep.origin + b*sqrt(max(r^2-sep.rmin^2, 0))*sep.tangent
	u = max(r-sep.rmin, 0)
	# parabola arc
	return sep.origin + u*sep.normal + b*sqrt(u)*sep.tangent
end

"""    approximate(separator, r1, r2, atol)

Approximates the + branch of parabolic separator with precision `atol`
by a polygonal path. Returns vector of distance parameters."""
function approximate(sep::Separator, r1, r2, atol)
	iszero(sep.tangent) && return [r1, r2] # degenerate case
	isstraight(sep) && return [r1, r2]
	nt = √(norm²(sep.tangent))
	x1, x2 = nt*√(r1-sep.rmin), nt*√(r2-sep.rmin)
	x = approxparabola(sep.rmin, x1, x2, atol)
	y = sep.rmin .+ (x ./ nt) .^ 2
	y[begin] = r1; y[end] = r2
	return y
end

"""    atan(separator)
Returns the angle of the initial normal of this separator."""
@inline Base.atan(sep::Separator) = atan(sep.normal[2], sep.normal[1])
# Capture ««
	# edge is captured if for all r, d(p(r), L) < r
	# i.e. min(r² - d²(p(r), L)) > 0
"""    capture(separator, b1, r1, b2, r2, line)

Returns a positive value iff the separator arc [r1,r2] is entirely closer
to the line than to either of its sides."""
function capture(sep::Separator, b1, r1, b2, r2, l::Line)
	# this returns a positive value iff, for all r: d²(eval(sep, r), L) < r
	# i.e. min(r²-d²(eval(sep, r), L)) > 0
	# Let p(r) = eval(sep, r), N=l.normal, A=l.offset;
	# we return min(f(r)), f(r) = r²-(N.p(r)+A)².
	if b1 == 2 || b2 == 2
		error("not implemented: separator from incorrect node")
	end
	if isparallel(sep)
		error("not implemented")
	end
	if isstraight(sep)
		# in this case: p(r) = o+ht, where h=±√(r²-r₀²)
		# so that f(r) = r₀²+h²-(N.o+A+ h N.t)²
		#              = h²(1-(N.t)²) - 2(N.t)(N.o+A) h + r₀²-(N.o+A)²
		# we return the minimum of f on [h1, h2]:
		Nt = l.normal ⋅ sep.tangent
		Noa= l.normal ⋅ sep.origin + l.offset
		r0 = sep.rmin
		a = 1 - Nt^2
		b = -Nt*Noa
		c = r0^2 - Noa^2
		return min_quadratic((a, b, c), (b1*sqrt(r1^2-r0^2), b2*sqrt(r2^2-r0^2)))
	end
	# in this case, p(r) = o+ht + h²n, where h=±√(r-r₀)  (so that r=r₀+h²)
	# f(r) = r₀+h²-(N.o+A + hN.t + h²N.n)²
	#      = -(N.n)²h⁴ -2(N.t)(N.n) h³ + (1-(N.t)²-(N.o+A)²) h²
	#        - 2(N.o-A)(N.t) h + r₀ - (N.o+A)²
	Nn = l.normal ⋅ sep.normal
	Nt = l.normal ⋅ sep.tangent
	Noa= l.normal ⋅ sep.origin + l.offset
	r0 = sep.rmin
	a4 = -Nn^2
	a3 = -2*Nt*Nn
	a2 = 1-Nt^2-Noa^2
	a1 = -2*Noa*Nt
	a0 = r0 - Noa^2
	return min_quartic((a0,a1,a2,a3,a4), (b1*sqrt(r1-r0), b2*sqrt(r2-r0)))
end
# Tripoints ««1
# Slight extension of `Base.sign`:
@inline Base.sign(T::Type{<:Integer}, x) =
	iszero(x) ? zero(T) : (x > 0) ? one(T) : -one(T)

@inline _BAD_TRIPOINT(x) = (oftype(x, NaN), Int8(2), Int8(2), Int8(2))
# docstring ««
"""    tripoint(a,b,c)

This computes the tripoint (equidistant point) of a triple of cells,
each of which is either a point or an (oriented) line.
The cells are cyclically ordered:
tripoint(a,b,c) == tripoint(b,c,a) ≠ tripoint(c,b,a).
The data is returned as `(radius, branch1, branch2, branch3)`.
If no such tripoint exists, `nan, 2,2,2` is returned.

The branch positions are returned as `Int8`, encoded as:
branch1 = +1 iff the tripoint lies on branch seeing a1 on its **left**
and a2 on its **right**, i.e. a1↑a2
(the arrow marks the direction of increasing r on this branch).
Likewise, branch2 is +1 for a2↑a3, and branch3 is +1 for a3↑a1.

For example, the center of an equilateral triangle has branches +1,+1,+1.
"""
function tripoint end #»»
function tripoint(a::AbstractVector, b::AbstractVector, c::AbstractVector)#««
	ab, bc, ca = b-a, c-b, a-c
	det2(ca, ab) > 0 || return _BAD_TRIPOINT(ab[1])
	r = √(norm²(ab)*norm²(ca)*norm²(bc))/(2*abs(det2(ab,ca)))
	return r, -sign(Int8, bc ⋅ ca), -sign(Int8, ca ⋅ ab), -sign(Int8, ab ⋅ bc)
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
			s = sign(Int8, x3)
			return (x3^2+y3^2)/(-2*f*y3), Int8(1), -s, s
		elseif x2 == x1
			# possible positions for p3: only the top ones are (123) tripoints
			#        ↑ 
			#  \ 3  ↗↑↖ 3
			#   `--' │ ╲
			# ━━━1━━━┥2
			#   ,--. │ ╱
			#  /  X \│╱ X
			y3 < 0 && return _BAD_TRIPOINT(x1)
			s = sign(Int8, x3-x1)
			return ((x3-x1)^2+y3^2)/(2*f*y3), Int8(1), s, -s
		else
			throw(PointInSegment())
		end
	end#»»
	if iszero(y3)#««
		if iszero(x3) # p3 is start of segment 1
			y2 < 0 && return _BAD_TRIPOINT(x1)
			s = sign(Int8, x2)
			return (x2^2+y2^2)/(2*f*y2), s, -s, Int8(1)
		elseif x3 == x1
			y2 > 0 && return _BAD_TRIPOINT(x1)
			s = sign(Int8, x2-x1)
			return ((x2-x1)^2+y2^2)/(-2*f*y2), -s, s, Int8(1)
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
		return r, sign(Int8, y3-y2), Int8(1), sign(Int8, y2-y3)
	end
	
	s0 = y3^2+dx^2-y2*y3
	s1 = 4*y2*y3 - dx^2
	s2 = y2^2+dx^2-y2*y3
	return r, sign(Int8, s0), sign(Int8, s1), sign(Int8, s2)
end#»»
function tripoint((p1,q1)::Segment, (p2,q2)::Segment, p3::AbstractVector)#««
# 	println("\e[36;7mtripoint(($p1,$q1), ($p2,$q2), $p3)\e[m")
	u1, u2 = q1-p1, q2-p2
	l1, l2 = √norm²(u1), √norm²(u2)
	a1, a2 = det2(u1, p3-p1), det2(u2, p3-p2)
	# Since this does not depend on the orientation of both lines,
	# we reorient them so that p3 lies on the positive side of both lines:
	(a1 < 0) && ((a1, u1, p1, q1) = (-a1, -u1, q1, p1))
	(a2 < 0) && ((a2, u2, p2, q2) = (-a2, -u2, q2, p2))
	c, s = u1⋅u2, det2(u1, u2)
	if iszero(s) # parallel segments case
		h1 = det2(u1, p2-p1)
		(h1 < 0) && ((h1, a1) = (-h1, -a1)) # swap p1, q1
		(a1 < 0) || (a1 > h1) && return _BAD_TRIPOINT(c)
		if iszero(a1) # p3 ∈ line1; must be left of segment1
			x3 = (p3-p1)⋅u1 # x3/l1
			x3 > 0 && return _BAD_TRIPOINT(c)
			return h1/(2l1), Int8(0), Int8(-1), Int8(1)
		elseif a1 == h1
			xp3, xq3 = (p3-p2)⋅u1, (p3-q2)⋅u1
			(xp3 > 0) || (xq3 > 0) && return _BAD_TRIPOINT(c)
			return h1/(2l1), Int8(0), Int8(1), Int8(-1)
		end
		return h1/(2l1), Int8(0), Int8(1), Int8(-1)
	end
	# Taking coordinates: u=(1,0), v=(c,s), p3=((ac-b)/s, a);
	# bisector B₁₂ is parametrized as r ↦ r*((c-1)/s,1).
	# This is equidistant from L1 and P3 whenever
	# (1-c)²r² - 2(a+b)(1-c)r + (a²+b²-2abc) = 0;
	# Δ'=2ab(1+c) hence (1-c)r = a+b±√Δ'.
	# Geometry gives the smaller root (-√Δ') when s>0, +√Δ' when s<0.
	#
	# Orientation of the branches: the tripoint is the apex of the parabola P₁₃
	# iff (c-1)r/s == (ac-b)/s i.e. a(1+c) = 2b; apex of P₂₃ iff b(1+c)=2a.
	# The branch on the line/line bisector depends on the relative position
	# of both segments...
	p1p2, p1q2 = p2-p1, q2-p1
	Dp2, Dq2, Dp1, Dq1 =
		det2(u1, p1p2), det2(u1, p1q2), det2(p1p2, u2), det2(u2, q1-p2)
# 	@assert Dq2 ≈ Dp2+det2(p1q1, p2q2)
# 	@assert Dq1 ≈ Dp1-det2(p1q1, p2q2)
# 	@assert Dp2 < Dq2
# 	@assert Dp1 > Dq1
	if iszero(a1) # special case: point is on one of both segments««
		iszero(a2) && return a1, Int8(1), Int8(1), Int8(1)
		# exchange p1↔q1 (a1 ← -a1 is a no-op)
		(s < 0) && ((s, c, Dp1, Dq1, Dp2, Dq2) = (-s, -c, Dq1, Dp1, -Dp2, -Dq2))
# 		println("  point is on line 1; (Dp1,Dq1,a2) = ", (Dp1,Dq1,a2))
# 		println("  both radii are: ", (a2*l1/(l1*l2+c), a2*l1/(l1*l2-c)))
		if a2 ≥ Dp1
# 			println("point is *left of* segment 1")
			Dq2 ≤ 0 && return _BAD_TRIPOINT(c)
			Dp1 ≤ 0 && return a2*l1/(l1*l2+c), Int8(-1), Int8(1), Int8(1)
			return a2*l1/(l1*l2-c), Int8(+1), Int8(-1), Int8(1)
		elseif a2 ≤ Dq1
# 			println("point is *right of* segment 1")
			Dp2 ≥ 0 && return _BAD_TRIPOINT(c)
			Dq1 ≥ 0 && return a2*l1/(l1*l2+c), Int8(-1), Int8(1), Int8(1)
			return a2*l1/(l1*l2-c), Int8(1), Int8(-1), Int8(1)
		else
			throw(PointInSegment())
		end
	elseif iszero(a2)
		# exchange p2, q2 (b ← -b is a no-op)
		(s < 0) && ((s, c, Dp1, Dq1, Dp2, Dq2) = (-s, -c, -Dp1, -Dq1, Dq2, Dp2))
		# FIXME: block some cases where the (12) branch is invalid
# 		println("  point is on line 2; (Dp2,Dq2,a1) = ", (Dp2,Dq2,a1))
# 		println("two radii: ", (a1*l2/(l1*l2+c), a1*l2/(l1*l2-c)))
		if a1 ≥ Dq2
# 			println("point is *after* segment 2")
			Dp1 < 0 && return _BAD_TRIPOINT(c)
			Dq2 ≤ 0 && return a1*l2/(l1*l2+c), Int8(-1), Int8(1), Int8(1)
			return a1*l2/(l1*l2-c), Int8(1), Int8(1), Int8(-1)
		elseif a1 ≤ Dp2
# 			println("point is *before* segment 2")
			Dq1 > 0 && return _BAD_TRIPOINT(c)
			Dp2 ≥ 0 && return a1*l2/(l1*l2+c), Int8(-1), Int8(1), Int8(1)
			return a1*l2/(l1*l2-c), Int8(1), Int8(1), Int8(-1)
		else
			throw(PointInSegment())
		end
	end#»»
	d = √(2*a1*a2*(l1*l2+c)) # normalization factor l1*l2
	# the factors indicating the position of the apex are: 2a-(1+c)b, 2b-(1+c)a
	# after renormalization: 2l2^2*a1-(l1*l2+c)*a2, 2l1^2*a2-(l1*l2+c)*a1
	if s > 0
		b1 = (Dp2 ≥ 0) ? (Dp1 < 0 ? Int8(2) : Int8(1)) : # seg2 left of line1
		     (Dq2 ≤ 0) ? (Dp1 ≤ 0 ? Int8(-1) : Int8(2)) :
				 Dq1 ≥ 0 ? Int8(1) : Dp1 ≤ 0 ? Int8(2) : throw(CrossingSegments())
		b1 == Int8(2) && return _BAD_TRIPOINT(c)
		b2 = sign(Int8, 2*l2^2*a1-(l1*l2+c)*a2)
		b3 = sign(Int8, 2*l1^2*a2-(l1*l2+c)*a1)
		r = (a1*l2+a2*l1-d)/(l1*l2-c)
		return r, b1, b2, b3
	else
		b1 = (Dq2 ≥ 0) ? (Dq1 < 0 ? Int8(2) : Int8(-1)) :
		     (Dp2 ≤ 0) ? (Dq1 ≤ 0 ? Int8(1) : Int8(2)) :
				 Dp1 ≥ 0 ? Int8(-1) : Dq1≤0 ? Int8(2) : throw(CrossingSegments())
		b1 == Int8(2) && return _BAD_TRIPOINT(c)
		r = (a1*l2+a2*l1+d)/(l1*l2-c)
		return r, b1, Int8(1), Int8(1)
	end
end#»»
function tripoint(l1::Line, l2::Line, l3::Line)#««
	# let lᵢ: aᵢ x + bᵢ y + cᵢ = 0
	# also define cᵢⱼ = aᵢ aⱼ + bᵢ bⱼ, sᵢⱼ = aᵢ bⱼ - bᵢ aⱼ
	# coordinate change: [u v]=[b1 -a1;a1 b1][x y]+[0 c1], then
	# [x y]=[b1 a1;-a1 b1][u v-c1] and the equation of l2 is now
	# a2*(b1 u+a1 (v-c1)) + b2*(-a1 u+b1 (v-c1)) + c2 = 0
	# (a2 b1 - a1 b2) u + (a1 a2 + b1 b2) (v-c1) + c2 = 0
	# -s12 u + c12 (v-c1) + c2 = 0
	# This has v=0 for u12= (c2-c12*c1)/s12.
	# Likewise, one finds u13 = (c3-c13*c1)/s13 = (-c3+c13*c1)/s31
	# and u12-u13 = (-(c12 s31 + s12 c31) c1 +s31c2 + s12 c3)/(s12 s31)
	#  = (s31 c2 + s12 c3 + s23 c1)/(s12 s13).
	#
	# Law of sines: d1/s23 = (∑ c1 s23)/(s12 s13 s23) = 2R
	# so the numerator is 2R*product of all sines
	c1, c2, c3 = l1.offset, l2.offset, l3.offset
# 	c12, c23, c31 = l1⋅l2, l2⋅l3, l3⋅l1
	s12, s23, s31 = det2(l1,l2), det2(l2,l3), det2(l3,l1)
# 	h1 = c1 + (s12*c3+s31*c2)/s23
# 	h2 = c2 + (s23*c2+s12*c3)/s31
# 	h3 = c3 + (s31*c1+s23*c1)/s12
	q = (c1*s23 + c2*s31 + c3*s12) / (s12*s23*s31)
	d1, d2, d3 = s23*q, s31*q, s12*q
	r = √(abs((d2+d3-d1)*(d3+d1-d2)*(d1+d2-d3)/(d1+d2+d3)))/2
# 	println("  triangle has sides $((d1,d2,d3))")
	z1, z2, z3 = sign(Int8, d1), sign(Int8, d2), sign(Int8, d3)
	z = sign(Int8, z1+z2+z3)
# 	println("  triangle has altitudes $((h1,h2,h3))")
	# if any altitude is zero, the three lines intersect in a single point:
# 	any(iszero, (h1, h2, h3)) && return _BAD_TRIPOINT(l1.offset)
	return r, z1*z, z2*z, z3*z
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
# 	println("\e[31;7m after insert!($q0, $c) = $stack:\e[m")
	while !isempty(stack)
		e = pop!(stack)
		@assert left(v, e) == c
# 		println("examining outgoing edge $e in cell $(tail(v,e)); right=$(right(v,e))")
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
# 	println(points)
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
	branch::Vector{Int8} # indexed by edges
	neighbours::Vector{J} # indexed by points

	@inline VoronoiDiagram{J,T}(triangulation::CornerTable, points, segments
		) where{J,T} =
		new{J,T}(triangulation, points, segments,
			Vector{SVector{2,T}}(undef, nnodes(triangulation)),
			Vector{T}(undef, nnodes(triangulation)),
			Vector{Separator{T}}(undef, nedges(triangulation)),
			Vector{Int8}(undef, nedges(triangulation)),
			Vector{J}(undef, length(points)))
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
@inline line(v::VoronoiDiagram, c::Cell) =
	let (a,b) = cellsegment(v, c); Line(point(v,a), point(v,b)) end
@inline geometricnode(v::VoronoiDiagram, q::Node) = v.geomnode[int(q)]
@inline geometricnode!(v::VoronoiDiagram, l::Pair{<:Node,<:AbstractVector}...)=
	for (q, p) in l; v.geomnode[int(q)] = p; end
@inline noderadius(v::VoronoiDiagram, q::Node) = v.noderadius[int(q)]
@inline noderadius!(v::VoronoiDiagram, l::Pair{<:Node,<:Real}...) =
	for (q, r) in l; v.noderadius[int(q)] = r; end
# @inline noderadius!(v::VoronoiDiagram, q::Node, r) = v.noderadius[int(q)] = r
@inline branch(v::VoronoiDiagram, e::Edge) = v.branch[int(e)]
@inline branch!(v::VoronoiDiagram, l::Pair{<:Edge,<:Integer}...) =
	for (e, b) in l; v.branch[int(e)] = b; end
@inline separator(v::VoronoiDiagram, e::Edge) = v.separator[int(e)]
@inline separator!(v::VoronoiDiagram, l::Pair{<:Edge,<:Separator}...) =
	for (e, s) in l; v.separator[int(e)] = s; end
@inline separators!(v::VoronoiDiagram, e::Edge, o::Edge, s::Separator) =
	separator!(v, e=>s, o=>reverse(s))


@inline edgedata(v::VoronoiDiagram, e::Edge) = separator(v,e), branch(v,e)
@inline edgedata!(v::VoronoiDiagram, e::Edge, (s, b)) =
	v.separator[int(e)], v.branch[int(e)] = s, b

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
 # (det2(ab,ag) ≥ 0) &&
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

# function meetscircle(v::AbstractVoronoi, q::Node, i, j)
# 	g, r = geometricnode(v, q), noderadius(v, q)
# 	a, b = point(v, i), point(v, j)
# 	ab2 = norm²(b-a)
# 	d = dot(g-a,b-a)
# # 	println("   seg($i,$j)-node($q=$(triangle(v,q))) distance:\n   a=$a\n   b=$b\n   g=$g\n   $(segdistance²(a,b,g))<? r^2=$r^2")
# 	return segdistance²(a, b, g) < r^2
# end
# Segment insertion ««2
"""    edgecapture(v,e)

Returns a positive value if this edge must be flipped, i.e. if all its points
are closer to right(e) than to tail(e)/head(e)."""
function edgecapture(v::VoronoiDiagram, e::Edge)#««
	# edge is captured if for all r, d(p(r), L) < r
	# i.e. min(r² - d²(p(r), L)) > 0
	println("\e[32;7medgecapture($(right(v,e)) = $(cellsegment(v,right(v,e))) -> $e = $(tail(v,e))⋯$(head(v,e))\e[m")
	o = opposite(v,e)
	be, bo = branch(v, e), branch(v, o)
	if bo == Int8(2)
		println("\e[32;1m right node $(node(o))=$(triangle(v,node(o))) is bad, must flip edge\e[m")
		return true
	end
	@assert issegment(v, right(v, e))
	(a, b) = cellsegment(v, right(v,e))
	q = node(e)
# 	if !influences(v, a, b, node(o))
# 		println("\e[32;3m segment ($a, $b) does not see new node $(node(o)), \e[1mmust flip\e[m")
# 		display((v, node(o)))
# 		return true
# 	end
	if !influences(v, a, b, q)
		println("\e[32m segment ($a,$b) does not see node $q=$(geometricnode(v,q))\e[m")
		return false
	end
	l = line(v, right(v,e))
	sep = separator(v, e)
	be, bo = branch(v,e), branch(v,o)
	re, ro = noderadius(v,node(e)), noderadius(v,node(o))
	z = capture(sep, be, re, -bo, ro, l)
# 	println("  \e[32mcapture($e) = $z\e[m")
	f = (-bo > be) || (-bo == be && ro > re)
	if f != (z > 0)
		println("\e[35;7m  (trying to connect $(right(v,e))→$(left(v,e)), disconnect $(tail(v,e))/$(head(v,e))) found f=$f, z=$z\e[m\e[35m")
	println("separator $e/$o = $(tail(v,e))/$(tail(v,o)) is $sep")
	println("  node($e)=$(node(e)) (in $(tail(v,e))) has parameter ($be, $re)")
	println("  node($o)=$(node(o)) (in $(tail(v,o))) has parameter ($bo, $ro)")
	println(evaluate(separator(v,e), be, re))
	println(evaluate(separator(v,o), bo, ro))
	println(evaluate(separator(v,e), -be, re))
	display((v,node(e)))
	display((v,node(o)))

	println("\e[35;1m view from $(tail(v,e)): $(-bo),$ro  increases to $be,$re\e[m")
	end
	return f
# 	return z > 0
end#»»
function addsegment!(v::VoronoiDiagram, c::Cell)#««
	a,b = cellsegment(v, c)
	println("\e[31;1;7minserting segment $c = ($a,$b)\e[m")
	display(v)
	q0 = findrootnode(v, a, b)
	print("\e[31mroot node is $q0:\e[m "); display((v,q0))
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
		println("   \e[7m flipping edge $e: connect $(right(v,e))->$(left(v,e)), disconnect $(head(v,e)) - $(tail(v,e))\e[m")
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

	# split segments in two
# 	splitsegments!(v)

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
# 	println("tripoint_llp($c1,$c2,$c3)")
# 	if c3 ∈ cellsegment(v, c1)
# 		c3 ∈ cellsegment(v, c2) &&
# 			return (zero(first(v.noderadius)), Int8(0), Int8(0), Int8(0))
# 		b1 = sum(cellsegment(v,c1)) - c3
# 		return tripoint_seg((point(v,c3), point(v,b1)), segment(v, c2), Int8(1))
# 		error("p3 ∈ s1: not implemented; other end is $b1")
# 	elseif c3 ∈ cellsegment(v, c2)
# 		b2 = sum(cellsegment(v,c2)) - c3
# 		return tripoint_seg((point(v,c3), point(v,b2)), segment(v,c1), Int8(-1))
# 		error("p3=$c3 ∈ s2=$c2: not implemented; other end is $b2")
# 	end
# 	# TODO: check if point c3 lies on either line; in this case, call
# 	# specialized code
	return tripoint(segment(v,c1), segment(v,c2), point(v,c3))
end
@inline function tripoint_lll(v::VoronoiDiagram, c1,c2,c3)
	return tripoint(line(v,c1), line(v,c2), line(v,c3))
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
	println("\e[34;7m nodedata!($q = $c1,$c2,$c3):\e[m")
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
	println((r, b1,b2, b3))
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

# @inline function geometricnode!(v::VoronoiDiagram, q::Node, g=equidistant(v,q))
# # 	println("geometricnode!($q)/$(nnodes(v))")
# 	v.geomnode[int(q)] = g
# 	s = cell(v, q, 1)
# 	v.noderadius[int(q)] = if issegment(v, s)
# 		(c1, c2) = cellsegment(v, s)
# 		segdistance²(point(v, c1), point(v, c2), g)
# 	else
# 		distance²(point(v, s), g)
# 	end
# end

function separator(v::AbstractVoronoi, c1, c2)#««
	if issegment(v, c1)
		i1, j1 = cellsegment(v, c1)
		a1, b1 = point(v, i1), point(v, j1)
		if issegment(v, c2)
			return Separator(segment(v, c1), segment(v, c2))
# 		elseif (c2 == i1 || c2 == j1) # degenerate parabola (exact)
# 			return degenerate_separator(line(v, c1), point(v, c2))
		else # generic parabola separator
			return Separator(segment(v, c1), point(v,c2))
		end
	elseif issegment(v, c2)
		i2, j2 = cellsegment(v, c2)
# 		if c1 == i2 || c1 == j2
# 			return reverse(degenerate_separator(line(v, c2), point(v, c1)))
# 		else # generic parabola
			return Separator(point(v, c1), segment(v, c2))
# 		end
	else # point-point separator
		return Separator(point(v, c1), point(v, c2))
	end
end#»»

# Computing Voronoi diagram ««1
# # Finding double edges««2
# # Double edge for a site s' = an edge between two nodes that are captured
# # by s', but which contains at least one non-captured point
# 
# """
# edgeexit(v, e, i, j)
# given the edge e and a new segment (i,j),
# returns the minimum for z∈e of d²(z,c1)-d²(z,(ij))
# (c1 being one of the two cells defining e)
# """
# function edgeexit(v::AbstractVoronoi, e, i, j)
# 	o = opposite(v, e)
# 	println("edgeexit($e/$o, $i, $j)")
# 	q0, q1 = node(e), node(o)
# 	g0, g1 = geometricnode(v, q0), geometricnode(v, q1)
# 	c1, c2 = tail(v,e), tail(v, o)
# 	g0 == g1 && return 0 # catch non-ternary nodes
# 	p,q = point(v,i), point(v,j)
# 	if issegment(v, c1)
# 		if issegment(v, c2)
# 			(i1, j1) = cellsegment(v, c1)
# 			(i2, j2) = cellsegment(v, c2)
# 			a1, b1, a2, b2 = point(v,i1), point(v,j1), point(v,i2), point(v,j2)
# 			return edgeexit_ss(a1,b1,a2,b2,g0,g1, p,q)
# 		end
# 		c1,c2 = c2,c1 # this replaces _sp case by equivalent _ps case just below:
# 	end
# 	if issegment(v, c2) # _ps case
# 		(i2, j2) = cellsegment(v, c2)
# 		a, a2, b2 = point(v,c1), point(v,i2), point(v,j2)
# 		return edgeexit_ps(a, a2,b2,g0,g1,p,q)
# 	else
# 		a, b = point(v,c1), point(v,c2)
# 		return edgeexit_pp(a,b,g0,g1, p,q)
# 	end
# end
# 
# "returns min on [0,1] of a*x^2+b*x+c"
# @inline function quadratic_min01(v)
# 	a,b,c = v
# 	a ≤ 0 && return min(c, a+b+c)
# 	(b ≥ 0) && return c
# 	(b ≤ -2a) && return (a+b+c)
# 	return (4*a*c-b*b)/(4*a)
# end
# 	
# function edgeexit_pp(a,b,g1,g2,p,q)
# 	# on the segment (g1,g2) of the mediator of (a,b),
# 	# compute min(d²(z,a)-d²(z,pq))
# 	
# 	# z = g1+t(g2-g1) = g1 + tv
# 	# d²(z,a) = ‖g1-a+tv‖² = ‖v‖² t² + 2(g1-a)⋅v t + ‖g1-a‖²
# 	# d²(z,pq) = <pqz>²/pq² = (<pqg1> + t <pqv>)²/pq²
# 	#  = <pqv>²/pq² t² + 2 <pqg₁><pqv>/pq² t² + <pqg₁>²/pq²
# 	v = g2-g1
# 	ag1 = g1 - a
# 	pq = q-p
# 	pqv = det2(pq,v)
# 	pqg1 = det2(pq,g1-p)
# 	da_quadratic = SA[norm²(v), 2*dot(ag1,v), norm²(ag1)]
# 	dpq_quadratic = SA[pqv^2, 2*pqg1*pqv, pqg1^2]/norm²(pq)
# 	return quadratic_min01(da_quadratic - dpq_quadratic)
# end
# 
# function edgeexit_ss(a,b, a2,b2, g1,g2, p,q)
# 	# on the segment (g1,g2) of the mediator of ((ab), (a2b2)),
# 	# compute min(d²(z,(ab)) - d²(z,(pq)))
# 	
# 	# z = g1 + tv just as in the pp case
# 	# d²(z,(ab)) = <abz>²/ab² = (<abg>+t<abv>)²/ab²
# 	# d²(z,(pq)) = <pqz>²/pq² = (<pqg>+t<pqv>)²/pq²
# 	v = g2-g1
# 	ag, ab = g1-a, b-a
# 	pg, pq = g1-p, q-p
# 	abg, abv, ab2 = det2(ab, ag), det2(ab, v), norm²(ab)
# 	pqg, pqv, pq2 = det2(pq, pg), det2(pq, v), norm²(pq)
# 	dab_quadratic = SA[abv^2, 2*abv*abg, abg^2]/ab2
# 	dpq_quadratic = SA[pqv^2, 2*pqv*pqg, pqg^2]/pq2
# 	return quadratic_min01(dab_quadratic - dpq_quadratic)
# end
# 
# function edgeexit_ps(a, b,c, g1,g2, p,q)
# 	# the intrinsic formulas here are a bit large...
# 	# we isometrically transform everything to a standard position instead,
# 	# by the transformation z ↦ (u⋅az, <u,bz>)
# 	# This maps the (bc) line to (y=0) and a to (0,d(bc,a))
# 	bc = c-b; u = bc/√(norm²(bc))
# 	xg1 = dot(u, g1-a) # no need to compute y for g1 and g2,
# 	xg2 = dot(u, g2-a) # which lie on the parabola
# 	xg1, xg2 = minmax(xg1, xg2)
# 	ya = det2(u,a-b)
# 	# now the parabola equation is 2y = ya + x^2/ya; x ∈ [xg1, xg2]
# 	# so that d²(z,bc) = y² = (ya+x²/ya)²/4
# 	fpx, fpy = (dot(u,p-a), det2(u,p-b))
# 	fqx, fqy = (dot(u,q-a), det2(u,q-b))
# 	pq2 = distance²(p,q)
# 	# now return the minimum for x ∈ [xg1, xg2] of
# 	# (ya + x²/ya)/2 - <fp, fq, z>/pq², where
# 	dpq = fpx*fqy-fqx*fpy
# 	dx = fpx-fqx
# 	dy = fpy-fqy
# 	#
# 	# This is a quartic in x:
# 	#= \\ gp verification script
# F(x) = (ya+x^2/ya)^2/4;
# G(x) = matdet([fpx,fpy,1;fqx,fqy,1;x,ya+x^2/ya,1])^2/pq^2;
# H = (F(x)-G(x))*4*ya^2*pq^2;
# dpq = fpx*fqy-fpy*fqx;
# dx = fpx-fqx;
# dy = fpy-fqy;
# 
# polcoeff(H,4,x)==pq^2-4*dx^2
# polcoeff(H,3,x)==8*ya*dx*dy
# polcoeff(H,2,x)==8*ya*dx*dpq+2*ya^2*(pq^2-4*dx^2-2*dy^2)
# polcoeff(H,1,x)==8*ya^2*dy*(dx*ya-dpq)
# polcoeff(H,0,x)==ya^4*pq^2-4*ya^2*(ya*dx-dpq)^2
# 	=#
# 	f = SA[pq2-4*dx^2, 8ya*dx*dy,
# 		8ya*dx*dpq+2*ya^2*(pq2-4*dx^2-2*dy^2),
# 		8ya^2*dy*(dx*ya-dpq),
# 		ya^4*pq2-4ya^2*(ya*dx-dpq)^2]/(4*pq2*ya^2)
# 	# this is guaranteed by geometry to have a unique minimum on ℝ
# 	# compute the root of the derivative by Newton's method (bounded iterations):
# 	r0 = 0.
# 	for _ in 1:10
# 		r1 = ((8*f[1]*r0+3*f[2])*r0^2-f[4])/(2*f[3]+r0*(6*f[2]+12*f[1]*r0))
# 		abs(r1 - r0) ≤ 1e-6*abs(r0) && break
# 		r0 = r1
# 	end
# 	(r0 < xg1) && (r0 = xg1)
# 	(r0 > xg2) && (r0 = xg2)
# 	return f[5]+r0*(f[4]+r0*(f[3]+r0*(f[2]+r0*f[1])))
# end
# 
# 
# Offset ««1
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
		s12, s21 = Cell(np+2i-1), Cell(np+2i)
		(c1,c2) = cellsegment(v, s12)
		(p1,p2) = point(v, c1), point(v, c2)
		println("\e[7m splitting cell $s12 = ($c1,$c2)\e[m")
# 		showcell(stdout, v, s12)
		e2 = anyedge(v, s12)
		while head(v, e2) ≠ c2; e2 = after(v, e2); end
		e1 = e2
		while head(v, e1) ≠ c1; tail!(v, e1, s21); e1 = after(v, e1); end
		println("found e1=$e1, e2=$e2")
		# split the cell by inserting two new nodes
		# (= split vertex s12 of the dual triangulation)
		#
		#   c1  o1|e1  s12 e2|o2 c2 becomes:
		#
		#   o1│q12  s21 e2│
		#     │    q21    │q23
		#  c1 q1─────────q2 c2
		#  q13│    q11    │
		#     │e1      q22│o2
		#
		o1, o2 = opposite(v, e1), opposite(v, e2)
		q1, q2 = newnodes!(v, 2)
		q11, q12, q13 = side(q1,1), side(q1,2), side(q1,3)
		q21, q22, q23 = side(q2,1), side(q2,2), side(q2,3)
		tail!(v, q11=>s12, q12=>s21, q13=>c1, q21=>s21, q22=>s12, q23=>c2)
		opposites!(v, q11=>q21, q12=>o1, q13=>e1, q22=>o2, q23=>e2)
		anyedge!(v, s12, q11); anyedge!(v, s21, q21)
		display((v, q1)); display((v,q2)); display((v, s12)); display((v, s21))
		# fix geometric information:
		l12, l21 = line(v, s12), line(v, s21)
		separators!(v, q11, q21, Separator(l12, l21))
		separator!(v, q12 => separator(v, e1), q22 => separator(v, e2),
			q13 => separator(v, o1))
		edgedata!(v, o1) # this also fixes q12
		edgedata!(v, e2) # this also fixes q23
		for e in star(v, q21); edgedata!(v, e); end
		for e in star(v, q21); nodedata!(v, node(e)); end
		@assert iszero(noderadius(v,q1))
		@assert iszero(noderadius(v,q2))
# 		branch!(v, q11=>0, q21=>0,
# 		)
# 		geometricnode!(v, q1, p1)
# 		geometricnode!(v, q2, p2)
	end
	return v
end#»»

# OffsetDiagram type and constructor ««2
struct OffsetDiagram{J,P,T} <: AbstractVoronoi{J}
	triangulation::CornerTable{J}
	voronoi::VoronoiDiagram{J,T}
	# parametrizations of edges by distance to sites:
	separator::Vector{Separator{T}} # indexed by edges
	# this encodes the direction of increasing distance of the cells
	# at the tail node of this edge; +1 = to the node, -1 = away from node
	# and 0 for a parallel bisector branch
	branch::Vector{Int8} # indexed by edges
	# number of neighbours of each point:
	# (this is used for deciding where to stop offsetting)
	neighbours::Vector{J}
end

@inline voronoi(v::OffsetDiagram) = v.voronoi
for f in (:(CornerTables.triangulation), :npoints, :nsegments, :point,
	:cellsegment, :geometricnode, :noderadius)
	@eval @inline $f(t::OffsetDiagram, args...;kwargs...) =
		$f(voronoi(t), args...;kwargs...)
end


@inline OffsetDiagram(points, segments; kw...) =
	OffsetDiagram{Int32}(points, segments; kw...)
function OffsetDiagram{J}(points::AbstractVector{P}, segments;
		extra = 0) where{J,P}
	v = VoronoiDiagram{J}(points, segments; extra)
# 	println("\e[1;7m in OffsetDiagram, got the following Voronoi diagram:\e[m")
# 	showall(v)
	splitsegments!(v)
println("\e[1;7m after split segments!:\e[m")
print(stdout, v.points)
showall(v)
gnuplot(v)
	# replace all node radii (squared distance) by their square roots:
	map!(√, v.noderadius, v.noderadius)
	# compute number of neighbours of each point:
	neighbours = zeros(J, npoints(v))
	for s in segments, a in s
		neighbours[a]+= one(J)
	end

	seplist = Vector{Separator{eltype(P)}}(undef, nedges(v))
	branch = Vector{Int8}(undef, nedges(v))
	for e in eachedge(v)
		o = opposite(v, e)
		o < e && continue
		# the separator is oriented with tail(o) = head(e) on its right,
		# i.e. pointing to the *left* of a:
		#       ╲n  head    ╱ 
		#  left  ╲____o____╱  right
		# +<⋯⋯   ╱    e    ╲  ⋯⋯>-
		#      p╱   tail    ╲
		#
		# branch[e] = +1 iff node(e) lies on the + branch of the separator
		# branch[e] = 0 iff this separator is a parallel bisector
		sep = separator(v, tail(v,o), tail(v,e))
		seplist[int(e)] = sep
		seplist[int(o)] = reverse(sep)

		if any(isnan, sep.normal)
			branch[int(e)] = branch[int(o)] = Int8(0)
			continue
		end

		q0, q1 = node(e), node(o)
		g0, g1 = geometricnode(v, q0), geometricnode(v, q1)
		r0, r1 = noderadius(v, q0), noderadius(v, q1)

		g0p, g0m = evaluate(sep, +1, r0), evaluate(sep, -1, r0)
		g1p, g1m = evaluate(sep, +1, r1), evaluate(sep, -1, r1)
		@debug """
edge $e (opposite $o): identifying node branches
  $q0 = $g0 at $r0
     (branch + $(g0≈g0p) $g0p
             - $(g0≈g0m) $g0m)
  $q1 = $g1 at $r1
     (branch + $(g1≈g1p) $g1p
             - $(g1≈g1m) $g1m)
  sep = $sep
"""

# 		println(" g0=$g0 $r0\n g1=$g1 $r1\n")
# 		println("  g0+=$g0p\n  g0-=$g0m\n  g1+=$g1p\n  g1-=$g1m\n")
		if r0 == sep.rmin
			if r1 == sep.rmin
				branch[int(e)] = branch[int(o)] = Int8(0)
			else
				@assert g1 ≈ g1m
				branch[int(e)] = Int8(-1); branch[int(o)] = Int8(+1)
			end
		elseif r1 == sep.rmin
			@assert g0 ≈ g0p
			branch[int(e)] = Int8(+1); branch[int(o)] = Int8(-1)
		else
# 			branch[int(e)] = iscloser(g0, g0p, g0m) ? one(Int8) : -one(Int8)
# 			branch[int(o)] = iscloser(g1, g1p, g1m) ? one(Int8) : -one(Int8)
			    if g0 ≈ g0p; branch[int(e)] = Int8(+1)
			elseif g0 ≈ g0m; branch[int(e)] = Int8(-1)
			else; error("no branch found for g0"); end
			    if g1 ≈ g1p; branch[int(o)] = Int8(-1)
			elseif g1 ≈ g1m; branch[int(o)] = Int8(+1)
			else; error("no branch found for g1"); end
		end
		@debug " branches = $(branch[int(e)]), $(branch[int(o)]) "
# 		println("branch[$e, $o] = $(branch[int(e)]), $(branch[int(o)])")
	end

	return OffsetDiagram{J,P,eltype(P)}(v, seplist, branch, neighbours)
end

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
@inline function edgecross(v::OffsetDiagram, e::Edge, radius)#««
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

	if iszero(b0) # this is a parallel bisector
		@assert iszero(b1) "iszero(b0) -> iszero(b1)"
		@assert r0 == r1
		# TODO: find what to do with parallel bisectors
		return (false, false, f0 && (r1 > 0) && (r2 > 0))
	end
	@assert !iszero(b1) "!iszero(b0) -> !iszero(b1)"
	# depending on (b0,b1), the edge is oriented this way:
	# ++ +- -+ --
	# <> << >> ><
	if b0 > 0 && b1 > 0 # the edge uses the two separator branches:
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
function prevedge(v::OffsetDiagram, e0::Edge, r)#««
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
function nextedge(v::OffsetDiagram, e0::Edge, r)#««
	@assert r ≥ 0
	for e in reverse(star(v, e0))
		(_, br) = edgecross(v, e, r)
		br && return e
	end
	return zero(e0)
end#»»
"""    firstedge(v, c, r)
Returns the first edge for an offset segment in cell `c`."""
@inline firstedge(v::OffsetDiagram, c::Cell, r) = prevedge(v, anyedge(v, c), r)
"""    edgeinter(v, e, r)
Returns the status of 
"""

# Offset chain ««2
"    finds a segment starting at a and *not* going to b"
function nextsegment(v::OffsetDiagram, a::Cell, b::Cell)
	v.neighbours[int(a)] ≠ 2 && return zero(a)
# 	showcell(stdout, v, a)
	for e in star(v, a)
		c = head(v, e)
# 		println("  $e => $c: $(issegment(v,c))")
		issegment(v, c) || continue
# 		println("     $(cellsegment(v,c))")
		cellsegment(v, c)[2] ≠ b && return c
	end
	return zero(a)
end
function zerochains_plus(v::OffsetDiagram{J,P}) where{J,P}#««
	chains = Vector{Edge{J}}[]
	points = Vector{P}[]
	done = falses(ncells(v))
	for startcell in J(npoints(v)+1):J(2):ncells(v)
		done[startcell] && continue; c = Cell(startcell)
		l = Edge{J}[]; push!(chains, l)
		p = P[]; push!(points, p)
		while !iszero(c) && !done[int(c)]
			e0 = anyedge(v, c) # guaranteed to be the segment-split edge
			e1, e2 = after(v, e0), opposite(v, before(v, e0))
			a, b = cellsegment(v, c)
			done[int(c)] = done[int(b)] = true
			push!(l, e1, e2)
			push!(p, point(v, b))
			c = nextsegment(v, b, a)
		end
		# now complete the chain to the left (if open)
		(a, b) = cellsegment(v, Cell(startcell))
		while !done[int(a)]
			done[int(a)] = true
			pushfirst!(p, point(v, a))
			c = nextsegment(v, a, b)
			iszero(c) && break
			done[int(c)] && break; done[int(c)] = true
			e0 = anyedge(v, c)
			e1, e2 = after(v, e0), opposite(v, before(v, e0))
			a, b = cellsegment(v, c)
			pushfirst!(l, e1, e2)
		end
	end
	return (chains, points)
end#»»
function zerochains_reverse(v::OffsetDiagram{J}, zplus) where{J}#««
	chains = Vector{Edge{J}}[]
	(cplus, pplus) = zplus
	for c in cplus
		l = Edge{J}[]; push!(chains, l)
		for i in length(c):-2:1
			e1, e2 = c[i], c[i-1]
			push!(l, opposite(v, next(opposite(v, e1))),
				opposite(v, prev(opposite(v, e2))))
		end
	end
	return (chains, reverse.(pplus))
end#»»
"""    zerochains(v::OffsetDiagram, reversed)
Returns the canonical chain corresponding to zero offset
on either side of the trajectory.
"""
function zerochains(v::OffsetDiagram)
	zplus = zerochains_plus(v)
	zminus = zerochains_reverse(v, zplus)
	(zplus, zminus)
end


"""    offsetchains(v::OffsetDiagram, radius, reversed)

Returns the set of chains encoding the offset curves for this radius.
Each chain is represented as a list of edges. Each edge correspond
to one cell traversed by the offset curve; it is the edge where
the curve enters the cell. The last edge in the chain represents either
the closure of the chain (if identical to the first) or the opposite edge to
the endpoint of the curve (if the chain is open).
"""
function offsetchains(v::OffsetDiagram{J}, radius, reversed) where{J}#««
	# The last segment encodes the endpoint of the chain; it is either
	# identical to the first segment (if closed loop) or the opposite edge
	# of the true last point of the chain.
	#
	# At any point during this algorithm, the active chain is the last one.
	chains = Vector{Edge{J}}[]
	done = falses(nedges(v))
	@assert radius ≥ 0
	for startcell in J(npoints(v)+1+reversed):J(2):ncells(v)
		c = Cell(startcell); e = firstedge(v, c, radius)
		iszero(e) && continue # this cell is not traversed by the offset curve
		done[int(e)] && continue # we already visited this curve segment

		# if this edge is not already done, then it lies on a new chain:
		l = [e]; push!(chains, l)
		while true
			e = opposite(v, nextedge(v, last(l), radius)); c = tail(v, e)
			push!(l, e); done[int(e)] = true
			!issegment(v, c) && (v.neighbours[c] ≠ 2) && break
			e == first(l) && break
		end
		# if the chain is an open loop, we need to extend it to the left:
		first(l) == last(l) && continue
		while true
			e = prevedge(v, opposite(v, first(l)), radius); c = tail(v, e)
			!issegment(v, c) && (v.neighbours[c] ≠ 2) && break
			pushfirst!(l, e); done[int(e)] = true
		end
	end
# 	if radius < 0 # correct the orientation of all chains
# 		for l in chains
# 			reverse!(l)
# 			for (i, e) in pairs(l); l[i] = opposite(v, e); end
# 		end
# 	end
	return chains
end#»»

# Offset ««2
"""    interpolate(v, chain, radius, atol, start=1)

Interpolates an arc of ∂R as a polygonal pathwith absolute precision `atol`.
Returns (P = list of points, L = list of indices),
so that chain[i] corresponds to the points P[L[i]:L[i+1]].
"""
function interpolate(v::OffsetDiagram{J,P,T}, chain, radius, atol) where{J,P,T}
	e0 = first(chain); sep0 = separator(v, e0)
	r = abs(radius)
	δ = √(r/(8*atol)) # used for discretizing circle arcs
	plist = [evaluate(sep0, +1, r)]
	llist = [1]
	for e1 in chain[2:end]
		sep1 = separator(v, e1)
		c = tail(v, e0)
		if !issegment(v, c) # circular arc
			@assert issegment(v, tail(v, e1))
			a0, a1 = atan(sep0), atan(sep1)
			a1 < a0 && (a1+= 2π)
			n = ceil(Int, (a1-a0)*δ)
			p = point(v, c); θ = (a1-a0)/n
			for i in 1:n-1
				a = a0+i*θ
				push!(plist, SA[p[1] + cos(a)*r, p[2]+sin(a)*r])
			end
		end
		# this is either the single point (for a straight segment), or the
		# closure of a circular arc
		push!(plist, evaluate(sep1, +1, r))
		push!(llist, length(plist))
		e0, sep0 = e1, sep1
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
	v = OffsetDiagram(points, segments)
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
	v = OffsetDiagram(points, segments)
	chains = [ offsetchains(v, abs(r), r < 0) for r in radii ]
	[[ interpolate(v, l, abs(r), atol)[1] for l in chains ] for r in radii ]
end
# Extrusion««1
# These functions compute a triangulation of the difference of two offset
# regions R(r2)∖R(r1), where r2 > r1.
# Point collection ««2
struct PointList{J,P} <: AbstractVector{P}
	points::Vector{P}
	index::Dict{P,J}
end
@inline PointList{J,P}() where{J,P} = PointList{J,P}(P[], Dict{P,J}())
@inline Base.size(plist::PointList) = (length(plist.points),)
@inline Base.getindex(plist::PointList, i::Integer) = plist.points[i]

function Base.push!(plist::PointList{J,P}, p) where{J,P} # returns index
	k = get(plist.index, p, zero(J))
	!iszero(k) && return k
	push!(plist.points, p)
	k = J(length(plist))
	plist.index[p] = k
	return k
end
@inline Base.append!(plist::PointList{J}, p) where{J} =
	J[ push!(plist, x) for x in p ]

function Base.show(io::IO, plist::PointList)
	for (i, p) in pairs(plist.points)
		println(io, " ",i, ": ", p)
		if plist.index[p] ≠ i
			println(io, "\e[31;7m bad index[$p] = $(plist.index[p]), should be $i")
			error("bad index")
		end
	end
# 	for (p, i) in pairs(plist.index)
# 		@assert plist.points[i] == p
# 	end
end

struct Affine3{T}
	a::T
	b::T
	r1::T
	r2::T
	z1::T
	z2::T
end
function Affine3((r1,z1)::Pair{T}, (r2,z2)::Pair{T}) where{T}
	a = (z2-z1)/(r2-r1)
	b = (z1*r2-r1*z2)/(r2-r1)
	return Affine3{T}(a, b, r1, r2, z1, z2)
end
@inline Base.:-(aff::Affine3) =
	Affine3(-aff.a, aff.b, -aff.r1, -aff.r2, aff.z1, aff.z2)
function evaluate(aff::Affine3, r)
	# exact cases:
	r == aff.r1 && return aff.z1
	r == aff.r2 && return aff.z2
	return aff.a*r + aff.b
end


# Axial extrusion of a single point ««2
"""    AxialExtrude{J}

Represents the offset of a single point along the trajectory,
as points indexed by edge crossings."""
struct AxialExtrude{J}
	chains::Vector{Vector{Edge{J}}}
	indices::Dict{J,Vector{J}}
# 	indices::Dict{J,NTuple{2,J}}
end
		
@inline npoints(a::AxialExtrude) = sum(length.(values(a.indices)))
@inline indices(a::AxialExtrude, e::Edge) = a.indices[int(e)]
function Base.reverse(a::AxialExtrude{J}, chains = reverse.(a.chains)) where{J}
	# assumes that each chain is reversed in place
	# c'1 => reverse(indices[cn])
	# c'2 => reverse(indices[c(n-1)])
	indices = empty(a.indices)
	println("\e[35;1m", (reverse, a.chains, chains), "\e[m")
	for (c1, c2) in zip(a.chains, chains)
		println("reverse $c1 => $c2")
		n = length(c1)
		@assert length(c2) == n
		for i in 1:n-1
			indices[int(c2[i])] = reverse(a.indices[int(c1[n-i])])
		end
	end
	return AxialExtrude(chains, indices)
end

function Base.show(io::IO, a::AxialExtrude)
	println(io, "axial offset of ", join(length.(a.chains),"+"),
		" crossings, $(npoints(a)) points:")
	for c in a.chains
		println(io, "  chain of $(length(c)) points: ")
		for e in c[1:end-1]
			println(io, "   edge $e -> points $(indices(a, e)) ->")
		end
		println(io, "   last edge $(last(c)) ",
			last(c) == first(c) ? "(closed)" : "(open)")
	end
end
function AxialExtrude(v::OffsetDiagram{J}, points, p, atol, zchains) where{J}
	rp = abs(p[1])
	indices = Dict{J,Vector{J}}()
	if iszero(rp)
		chains = zchains[1][1]
		for (ch, pts) in zip(chains, zchains[1][2])
			idx = append!(points, [ [q;p[2]] for q in pts ])
			for i in 1:length(ch)
				indices[int(ch[i])] =
					isodd(i) ? [idx[(i+1)>>1], idx[(i+3)>>1]] : [idx[(i+2)>>1]]
			end
		end
		return AxialExtrude{J}(chains, indices)
	end
	chains = offsetchains(v, rp, p[1] < 0)
	for ch in chains
		np = J(length(points))
		(newpoints, idx) = interpolate(v, ch, rp, atol)
		for i in 1:length(idx)-1
			j = idx[i]:idx[i+1]
			ind = append!(points, [[q; p[2]] for q in newpoints[j]])
			indices[int(ch[i])] = ind
		end
	end
	return AxialExtrude{J}(chains, indices)
end

@inline reverse2(x) = [reverse.(y) for y in x]
"""    standard_orientation(p, axp, q, axq, zchains)

Sort points (p,q) into (p1=closest to trajectory, p2=farthest),
and return corresponding axial extrusions, together with
the affine map (z = a*r+b).
"""
@inline function standard_orientation(p, axp, q, axq, zchains)
	rp, rq = p[1], q[1]
	@assert rp ≠ rq # we should already know that the face is not vertical
	aff = Affine3(p[1] => p[2], q[1] => q[2])
	println("\e[35;3m standard orientation ← rp=$rp, rq=$rq, aff=$aff\e[m")
	println("-aff = $(-aff)")
	if rp < rq
		if 0 ≤ rp # 0 ≤ rp < rq
			return (rp, axp, rq, axq, aff, false)
		else      # rp < rq ≤ 0
			iszero(rq) && (axq = reverse(axq, zchains[2][1]))
			return (-rq, axq, -rp, axp, -aff, false)
		end
	else
		if 0 ≤ rq # 0 ≤ rq < rp
			return (rq, axq, rp, axp, aff, true)
		else      # rq < rp ≤ 0
			iszero(rp) && (axp = reverse(axp, zchains[2][1]))
			return (-rp, axp, -rq, axq, -aff, true)
		end
	end
end

# Region between chains««2
"""    cell_contour(v, e1, e2, c2next)
Returns a description of the contour of c ∩ (R(r2)∖R(r1)), where:
 - R(r1) enters c at e1 and exits at e2
 - c2next is the next-edge map for ∂R(r2)
The description is encoded as (edge, type), where type is:
1 for ∂R1 (in reverse direction), 2 for ∂R2,
3 for + branch of edge, 4 for - branch of edge.
"""
function cell_contour(v::OffsetDiagram, e1, e2, c2next)
	elist, etype = [e1], [Int8(1)]
	lastedge = opposite(v, e2)
	# orbit around the cell (fragment) c until we find the start point
	# for this edge
	e = e1
	while true
		o = opposite(v, e)
		branch(v, o) > 0 && (push!(elist, o); push!(etype, Int8(4)))
		branch(v, e) ≥ 0 && (push!(elist, e); push!(etype, Int8(3)))
		e == lastedge && break
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
returns all edges delimiting cells in R(r2)∖R(r1), as

OLD: a list of tuples
(cell, edges, edgetypes), with
 - edges = a (cyclic) list of edges around this cell,
 - edgetypes = list of types matching each edge, encoded as:
 1=∂R1, 2=∂R2, 3=positive edge, 4=negative edge.
"""
function chain_contour(v::OffsetDiagram{J}, chains1, chains2) where{J}#««
	# build next-edge map for chains2 ««
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
		e1 = first(ch)
		for e2 in ch[2:end]
			c = tail(v, e1)
			elist, etype = cell_contour(v, e1, e2, c2next)
			# ∂R1 enters the cell c at e1 and exits at o = opposite(e2)
			push!(r, (c, elist, etype))
			e1 = e2
		end
		# if this chain is closed then the last examined cell did the
		# closing;
		# otherwise we still need to produce info about the “outer” boundary
		# (TODO)
	end
	return r
end#»»
"""    edge_transverse(v, e, r1, r2, atol)
Produces a list of points approximating the open interval [r1,r2] on edge e
(in this order), on the positive branch.
Returns list of [point; r], where r is distance to trajectory.
"""
function edge_transverse(v::OffsetDiagram{J,P,T}, e, r1, r2, aff, atol#««
		) where{J,P,T}
	o = opposite(v, e)
	q0, q1 = node(e), node(o)
	g0, g1 = geometricnode(v, q0), geometricnode(v, q1)
	d0, d1 = noderadius(v, q0), noderadius(v, q1)
	sep = separator(v, e)
# 	println("edge_transverse($e, $r1:$r2)")

	# eliminate degenerate & parallel separators:
	g0 == g1 && return [SA[g0[1], g0[2], evaluate(aff, d0)]]
	d0 == d1 && return [SA[g0[1], g0[2], evaluate(aff, d0)],
		SA[g1[1], g1[2], evaluate(aff, d0)]]

	t1, t2 = max(sep.rmin, r1), min(d0, r2)
	tlist = approximate(sep, t1, t2, atol)

	# remove points with exact distance r1 or r2, since they are already
	# taken account of as points of ∂R1 or ∂R2:
	first(tlist) == r1 && popfirst!(tlist)
	last(tlist) == r2 && pop!(tlist)
# 	println("  tlist=$tlist")
	return [[evaluate(sep, t, +1); evaluate(aff, t)] for t in tlist]
end#»»
"""    edgepoints(v, points, edge, edgetype, tlist, ax1, ax2)

Returns an interpolated list of points along this edge.
"""
function edgepoints(v::OffsetDiagram, points, edge, edgetype,
		r1,ax1,r2,ax2, aff, atol)#««
	if edgetype == 1 # use segment from ∂R1, backwards
# 		println("use segment $edge from ∂R1 (backwards)")
		return reverse(indices(ax1, edge))
	elseif edgetype == 2 # segment from ∂R2, forwards
# 		println("use segment $edge from ∂R2 (forwards)")
		return indices(ax2, edge)
	elseif edgetype == 3 # edge e transversally, forwards
# 		println("use edge $edge transversally (forwards)")
		seg = edge_transverse(v, edge, r1, r2, aff, atol)
		x = append!(points, seg)
# 		println("  added $(length(x)) points: $x => $seg")
		return x
# 		return append!(points, seg)
	elseif edgetype == 4
# 		println("use edge $edge transversally (backwards)")
		seg = edge_transverse(v, edge, r1, r2, aff, atol)
		x = append!(points, seg)
# 		println("  added $(length(x)) points: $x => $seg")
		return reverse(x)
# 		return reverse(append!(points, seg))
	end
end#»»

# Extrusion of a polygonal loop««2
"""    extrude_loop(v, loop)

Extrudes a loop of points [xi, yi] along the polygonal path(s);
returns (points, triangulation).
"""
function extrude_loop(v::OffsetDiagram{J,P,T}, loop, atol) where{J,P,T}
	# insert new points in the loop when it crosses [x=0] ««
	p = last(loop)
	newloop = []
	for q in loop
		(p[1]*q[1] < 0) && push!(newloop, SA[0, (p[1]*q[2]-p[2]*q[1])/(p[1]-q[1])])
		push!(newloop, q)
		p = q
	end
	loop = newloop
	#»»
	# axial paths: extrusions of individual points of the loop««
	points = PointList{J,SVector{3,T}}()
	zchains = zerochains(v)
	axial = [ AxialExtrude(v, points, p, atol, zchains) for p in loop ]
	println("\e[34m$points\e[m")
	for (p, ax) in zip(loop, axial)
		println("\e[36;1m extrusion of $p is:\e[m\n$ax")
	end
	# zero offset on the negative side:
	triangles = NTuple{3,Int}[]
# 	for (p, axp) in zip(loop, axial); println(p => axp); end
	# triangulate between consecutive axial paths««
	p, axp = last(loop), last(axial)
	for (q, axq) in zip(loop, axial)
		println("\e[1;7mtriangulate face: $p -> $q\e[m\n  axp=$axp\n  axq=$axq")
		if p[1] == q[1] # vertical face: easy case««
# 			println("\e[1mface is vertical\e[m")
			@assert axp.chains == axq.chains
			# axp and axq are composed of matched chains
			for c in axp.chains
				if first(c) == last(c)
# 					println("  build torus: $c")
				else
# 					println("  build tube: $c")
				end
# 					ip, iq = first(ap), first(aq)
# 					for i in 2:length(ap)
# 						jp, jq = ap[i], aq[i]
# # 						push!(triangles, (ip, jp, iq), (jp, jq, q))
# 						ip, iq = jp, jq
# 					end
# 					# TODO: end caps!
# 				end
			end
			#»»
		else # oblique face ««
			println("\e[1mface is oblique\e[m")
			r1, ax1, r2, ax2, aff, reversed =
				standard_orientation(p, axp, q, axq, zchains)
			println("  $r1:$r2 (reversed $reversed); z=$(aff.a)*r + $(aff.b)\n$r1: $ax1\n$r2: $ax2")
			println(" with aff = $aff, aff.b=$(aff.b)")
			# the surface between the axial paths for p1 and p2 is
			# split along the cells traversed by axial(p1) (= the closest one)
			# and each fragment is triangulated separately
			ct = chain_contour(v, ax1.chains, ax2.chains)
			println("contour:");for (c, el, et) in ct; println("   $c: $el, $et"); end
			for (c, elist, tlist) in ct
				println("\e[1m in cell $c: $elist, $tlist\e[m")
				cellpoints = Int[]
				for (edge, edgetype) in zip(elist, tlist)
					epoints = edgepoints(v, points, edge, edgetype,
						r1,ax1, r2,ax2, aff, atol)
					println("  ($edge, $edgetype) contributes $epoints = $(points[epoints])")
					append!(cellpoints, epoints)
				end
				println("  before unique!: $cellpoints")
				unique!(cellpoints)
				# build a loop for this cell fragment
				println("\e[36mcellpoints for $c = $cellpoints\e[m")
				for c in cellpoints
					println("  \e[36mpoint[$c] = $(points[c])\e[m")
				end
				if length(cellpoints) ≥ 3
					tri = triangulate_loop(points, cellpoints)
					for (a,b,c) in tri
						reversed && ((b,c) = (c,b))
						push!(triangles, (a,b,c))
					end
				end
			end
		end#»»
		p, axp = q, axq
	end
	(points, triangles)
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
	v = OffsetDiagram(plist, slist; extra)
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
	bc = @closure i->"X-0+!X"[clamp(3+branch(v,side(q,i)), 1:6)]
	print(io, "\e[33m", q, triangle(v,q), "\e[m: ",
		@sprintf("(%.3g,%.3g) r=%.3g", geometricnode(v, q)..., noderadius(v, q)),
		" ", bc(1), bc(2), bc(3),
	)
# 	for i in (1,2,3)
# 		e = side(q, i); o = opposite(v, e); oo = opposite(v,o)
# 		oo ≠ e && println(io, "  \e[31;7m opposite($o) = $oo, should be $e\e[m")
# 	end
end
# function showedge(io::IO, v::OffsetDiagram, a::Edge)
# 	if !iszero(branch(v, a))
# 		q = node(a)
# 		sep = separator(v, edge(v, a))
# 		g = geometricnode(v, q)
# 		r = √(noderadius(v, q))
# 		b = branch(v, a)
# 		h = evaluate(sep, r, b)
# 		println("evaluate(sep, $r, $b) = $h ≈ $g?")
# 		g ≈ h ||
# 			println(io, "  \e[31;7mgeometricnode($q) = $g; evaluate($r, $b) = $h\e[m")
# 	end
# end
# function Base.show(io::IO, ::MIME"text/plain", v::AbstractVoronoi)
# 	println(io, "\e[1m begin Voronoi diagram with $(nnodes(v)) nodes and $(ncells(v)) cells:\e[m")
# 	for q in eachnode(v); shownode(io, v, q); end
# 	for c in eachcell(v); showcell(io, v, c); end
# 	println(io, "\e[1m end Voronoi diagram\e[m")
# end
# function Base.show(io::IO, ::MIME"text/plain", v::OffsetDiagram)
# 	println(io, "\e[1m begin offset diagram with $(nnodes(v)) nodes and $(ncells(v)) cells:\e[m")
# 	for q in eachnode(v); shownode(io, v, q); end
# 	for c in eachcell(v); showcell(io, v, c); end
# 	println(io, "\e[1m end offset diagram\e[m")
# end
# end »»1
end

V=Voronoi
using StaticArrays
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
s1a = (SA[-5,0.],SA[-2,0.])
s1b = (SA[-2,0.],SA[2,0.])
s1c = (SA[2,0.],SA[5,0.])
s2a = (SA[-3,-6.],SA[-1,-2.])
s2b = (SA[-1,-2.],SA[1,2.])
s2c = (SA[1,2.],SA[3,6.])
s3a = (SA[-5,2.],SA[-2,2.])
s3b = (SA[-2,2.],SA[2,2.])
s3c = (SA[2,2.],SA[5,2.])

c1 = SA[0,0.]
c2 = SA[10,0.]
c3 = SA[5,1.]
c4 = SA[5,9.]
c8 = (c1, c2)
c9 = (c2, c3)
c10 = (c3, c4)



# v=V.VoronoiDiagram([[0.,0],[10,0],[5,1],[5,9]],[(1,2),(2,3),(3,4)];extra=0)

#
# v=V.OffsetDiagram([[0.,0],[10.,0],[10,10.]],[(1,2),(2,3)];extra=5)
# z = V.zerochains(v)
# el = V.extrude_loop(v, [[-.5,-1],[1,-.5],[.5,1],[-1,.5]], .1)

# v=V.OffsetDiagram([[0.,0],[10,0],[0,10],[10,10],[5,9],[5,1]],[(1,2),(2,6),(6,5),(5,4),(3,1)])
# l=V.offsetchains(v, 1., false)
# o=V.offset([[0.,0],[10,0],[0,10],[10,10],[5,9],[5,1]],[(1,2),(2,6),(6,5),(5,4),(3,1)], 1.)
# ci = V.chain_contour(v, V.offsetchains(v, 1., false), V.offsetchains(v, 10., false))
# el=V.extrude_loop(v, [SA[1.,1],SA[-1.,0],SA[1.,-1]], .1)
