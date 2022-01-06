# With some ideas taken from VRONI as described by [Huber 2008]:
# https://www.sciencedirect.com/science/article/pii/S0925772101000037
# https://www.sthu.org/research/publications/files/mscthesis.pdf
#
# Offset algorithm: [Kim 1998]
# https://www.sciencedirect.com/science/article/abs/pii/S0010448598000633
module Voronoi
using StaticArrays
using FastClosures
using LinearAlgebra
using LazyArrays
using Random

include("CornerTables.jl")
using .CornerTables
import .CornerTables: showall, showcell, shownode

# Geometry ««1
# Elementary geometry ««2

@inline det2(u,v) = u[1]*v[2]-u[2]*v[1]
@inline det2(u,v,w) = det2(v-u, w-u)
@inline norm²(v) = v[1]^2+v[2]^2
@inline distance²(a,b) = norm²(a-b)
@inline quarterturn(v) = SA[-v[2], v[1]]

@inline linethrough(a,b) =
	SA[a[2]-b[2], b[1]-a[1], det2(a,b)]

function lineinter(a,b,c,d)
	t = det2(a-c, c-d)
	D = det2(a-b, c-d)
	z = a+(t/D)*(b-a)
	return a+(t/D)*(b-a)
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

@inline iscloser(a,b,c) = dot(2a-b-c, b-c) ≥ 0
@inline isleft(a,b,c) = det2(a,b,c) > 0
@inline isleft(u,v) = det2(u,v) > 0

# Circumscribed circle ««2
function circumcenter(a,b,c)
	ab,ac = b-a,c-a
	m = SA[norm²(ab) ab[1] ab[2];norm²(ac) ac[1] ac[2]]
	kn = det(m[:,SA[2,3]])
	kx = det(m[:,SA[1,3]])/(2kn)
	ky = det(m[:,SA[2,1]])/(2kn)
	return a + SA[kx, ky]
end

function circumradius(a,b,c)
	ab,ac,bc = b-a, c-a, c-b
	return sqrt(norm²(ab)*norm²(ac)*norm²(bc))/(2*abs(det2(ab,ac)))
end

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

"""
    isincircle(a,b,c,p,q)

Returns true iff open segment ]p,q[ intersects circumcircle of triangle (a,b,c).
"""
function isincircle(a,b,c,p,q)
	a,b,c,q = a-p,b-p,c-p,q-p
	na, nb, nc = norm²(a), norm²(b), norm²(c)
	# equation of circumcircle is kn(x²+y²) - kx x - y y + k0 = 0
	m = SA[na a[1] a[2] 1;nb b[1] b[2] 1;nc c[1] c[2] 1]
	kn = det(m[:,SA[2,3,4]])
	kx = det(m[:,SA[1,3,4]])/kn
	ky = det(m[:,SA[1,4,2]])/kn
	C = det(m[:,SA[3,2,1]])/kn # equation for (t*q) is At² - Bt + C = 0
	A = norm²(q)
	B = kx*q[1] + ky*q[2]
	return (B ≥ 0) && (B ≤ 2A) && (B*B ≥ 4*A*C)
end

function incenter(a,b,c)
	la, lb, lc = distance²(b,c), distance²(c,a), distance²(a,b)
	p = la+lb+lc
	ra, rb, rc = la/p, lb/p, lc/p
	return ra*a + rb*b + rc*c
end

# Equidistant points ««2
# p = point, s = segment, x = segment starting on previous point

@inline function equidistant_pps(a, b, p, q)#««
	# returns the point equidistant from a, b, (pq)
	# chosen so that the cells are oriented (a, b, pq)
	# imagine our segment as an x-axis; both points are on the same side,
	# we reorient the segment so they both have y > 0
	pqa = det2(p,q,a)
	@assert !iszero(pqa)
	(pqa < 0) && ((p,q) = (q,p); pqa = -pqa)
	pqb = det2(p,q,b)
	@assert pqb > 0
	ab2 = distance²(a,b)
	pq2 = distance²(p,q)
	pqab = dot(q-p, b-a)
	if pqa == pqb
		# special case if (a,b) and (p,q) are collinear:
		@assert false
	end
	# let z = (a+b)/2 + t*I*(b-a); then
	# d²(z,a) = d²(z,b) = ab2(t^2 + 1/4)
	# d²(z,pq) = <pqz>²/pq2 = (pqa/2+pqb/2+t*pqab)^2/pq2, so the eq. is:
	# (using the identity (pq2*ab2 = pqab^2 + (pqa-pqb)^2):
	# (pqa-pqb)^2 * t^2 - (pqa+pqb) pqab t + (pqab^2/4 - pqa*pqb) = 0
	#
	# We find Δ = 4*pqa*pqb*ab2*pq2
	# and geometry implies that the +√ is the correct sign
	Δ = 4*pqa*pqb*ab2*pq2
	t = ((pqa+pqb)*pqab + √(Δ)) / (2*(pqa-pqb)^2)
	return (a+b)/2 + t*quarterturn(b-a)
end#»»

@inline function equidistant_pxs(a,b,p,q, ε)
	# return point equidistant from: a, (ab), (pq)
	# (in this order if ε=+1, opposite order if ε=-1)
	# z = a + tI(b-a) satisfies
	# d²(z,a) = d²(z,ab) = ab² t²
	# d²(z,pq) = <pqz>²/pq² = (<pqa> + t pq.ab)^2 / pq2
	# or: (pqa-pqb)^2*t^2 - 2 pqa*pqab*x - pqa^2 = 0
	# Δ = 4 pqa^2*pq2*ab2
	pqa = det2(p,q,a)
	pqb = det2(p,q,b)
	ab2 = distance²(a,b)
	if pqa == pqb # special case: both lines are parallel
		abp = det2(a,b,p)
		return a + abp/(2*ab2)*quarterturn(b-a)
	end
	pq2 = distance²(p,q)
	pqab = dot(q-p, b-a)
	Δ = 4*pqa^2*pq2*ab2
	t = (2*pqa*pqab+ ε*√(Δ)) / (2*(pqa-pqb)^2)
	z = a + t*quarterturn(b-a)
	return z
end

"returns the point equidistant from point a and segments (pq), (rs)"
function equidistant_pss(a, p, q, r, s)
	pqrs = det2(q-p, s-r)
	if iszero(pqrs) # special case: parallel lines
		pqa = det2(p,q,a)
		pqa < 0 && ((p,q, pqa) = (q,p, -pqa)) # ensure orientation
		pq = q-p
		pq2 = norm²(pq)
	pqr, pqs = det2(p,q,r), det2(p,q,s)
# 	# `a` is not allowed to be on any of the two lines pq, rs
# 	if pqr == pqs # special case: parallel lines
		# let v = quarterturn(q-p)
		λ = (pqr - 2*pqa)/(2*pq2)
		# c = a + λv is the projection of a on the middle line
		# this line is parametrized as z = c + t pq, then
		# az² = ac²+t² pq² = (λ²+t²) pq² must be pqr^2/(4 pq^2), hence
		# t^2 = (pqr^2 - (pqr-2pqa)^2)/(4 pq^4)
		#     = pqa(pqr-pqa)/pq^4
		# Geometry imposes the positive square root
		t = √(pqa*(pqr-pqa)) / pq2
		z = a+λ*quarterturn(pq)+t*pq
		return a + λ*quarterturn(pq) + t*pq
	end
	c = lineinter(p,q,r,s)
	ca = a-c
	pq = q-p; pq2 = norm²(pq); upq = pq/√(pq2)
	rs = s-r; urs = rs/√(norm²(rs))
	dot(upq, ca) < 0 && (upq = -upq)
	dot(urs, ca) < 0 && (urs = -urs)
	# parametrization of the inner angle bisector: z = c + t u
	c = lineinter(p, q, r, s)
	u = urs + upq
	# for z = c+tu: d²(z,pq) = t² pqu²/pq², while
	# d²(z,a) = ‖a-c-tu‖² = t² u²- 2t ca⋅u + ca²
	# the equation is thus t²(u²-pqu²/pq²) -2t ca⋅u + ca² = 0, or
	# t²(pq⋅u)²/pq² -2t ca⋅u + ca² = 0
	A = dot(pq, u)^2 / pq2
	B = dot(ca, u)
	C = norm²(ca)
	Δ = B^2-A*C
	ε = sign(det2(upq, urs))
	t = (B + ε*√(Δ))/A
	z = c+t*u
	return c + t*u
end

function equidistant_sss(a1,b1,a2,b2,a3,b3)
	# returns the point equidistant from the *oriented* lines (ai, bi)
	# (i.e. either the incenter, or an excenter, of the triangle, according
	# to the orientations).
	u1 = quarterturn(b1 - a1); u1 /= √(norm²(u1))
	u2 = quarterturn(b2 - a2); u2 /= √(norm²(u2))
	u3 = quarterturn(b3 - a3); u3 /= √(norm²(u3))
	p1 = lineinter(a2, b2, a3, b3)
	p2 = lineinter(a3, b3, a1, b1)
	return lineinter(p1, p1+u2+u3, p2, p2+u1+u3)
end

function equidistant_sss_parallel(a, u, b, p, pq)
	# returns the point equidistant from lines (a, a+u), (b,b+u), (p, pq)
	# (in this order)
	# parametrized as z = c+tu (with c = (a+b)/2)
	c, ab = (a+b)/2, b-a
	# d²(z, (a,a+u)) = l²/4 where l is the distance between the parallel lines
	# so (<pqc> + t <pqu>)² = |ab|⋅|pq|/2, or t=-<pqc>±|ab|.|pq|/(2<pqu>)
	# if <u,ab> > 0 then a lies to the right of the line (c, c+u)
	# and we take the + sign:
	pqc, pqu = det2(pq, c-p), det2(pq, u)
	l2, pq2 = det2(u,ab)^2/norm²(u), norm²(pq)
	t = (-pqc + sign(det2(u,ab))*sqrt(l2*pq2)/2)/pqu
	z = c + t*u
	return z
end

# Bisectors ««2
"mediator line of segment (ab)"
function mediator(a,b)
	return SA[2*(a[1]-b[1]), 2*(a[2]-b[2]), norm²(b)-norm²(a)]
end

# Lines ««2
#
Point{T} = SVector{2,T}
"the equation for a (normalized, oriented) straight line in the plane."
struct Line{T}
	# orientation: the normal vector points to the *right* of the line
	normal::SVector{2,T} # normalized to ‖u‖=1
	offset::T # line equation is normal⋅z + offset == 0
end
@inline Base.:*(a::Real, l::Line) = Line(a*l.normal, a*l.offset)
@inline normalize(l::Line) = (1/√(norm²(l.normal)))*l
@inline Base.:-(l::Line) = (-1)*l

@inline direction(l::Line) = quarterturn(l.normal)
# one arbitrary point on the line (actually the projection of the origin)
@inline point(l::Line) = -l.normal*l.offset/(l.normal[1]^2+l.normal[2]^2)

@inline Line(a::AbstractVector, b::AbstractVector) = # oriented from a to b
	normalize(Line(SA[b[2]-a[2], a[1]-b[1]], a[2]*b[1]-a[1]*b[2]))

# signed distance from l to a; positive sign corresponds to right side of l.
@inline evaluate(l::Line, a::AbstractVector) = dot(l.normal, a) + l.offset

"returns either the intersection point, or `nothing`."
function intersection(l1::Line, l2::Line)
	(a1, b1), c1 = l1.normal, l1.offset
	(a2, b2), c2 = l2.normal, l2.offset
	d = a1*b2-a2*b1
	iszero(d) && return nothing
	return SA[b1*c2-c1*b2,c1*a2-a1*c2]/d
end

# Parametrized bisectors ««2
# This is either:
#  - a parabola: origin + tangent*(t-t0) ± normal*sqrt(t-t0);
#  - a line: origin ± tangent*|t-t0|  (this is indicated by normal == 0)
# This could be rewritten as a union, but we can use one of the fields to
# store the type.
"""
    Separator

A structure holding the parametrization for the bisector between two sites
(both of them either a point or a vector),
The separator between `a` and `b` is parametrized by the distance to the sites,
and represented as two branches: the `+` branch sees the site `a` on its right.
"""
struct Separator{T}
	origin::SVector{2,T}
	tangent::SVector{2,T}
	normal::SVector{2,T}
	t0::T
end

"""
    evaluate(separator, t, sign)

Returns the point on the separator situated at distance `t` from both
sites and on the branch given by sign `s` (either + or -).
The `+` branch sees `a` on its right and `b` on its left.
"""
@inline function evaluate(b::Separator, t, s) # s is a sign
	if iszero(b.normal) # straight line
		return b.origin + s*sqrt(max(t^2-b.t0^2, 0))*b.tangent
	else # parabola arc
		u = max(t-b.t0, 0)
		return b.origin + u*b.normal + s*sqrt(u)*b.tangent
	end
end

function Separator(a::AbstractVector, b::AbstractVector)# bisector of two points
	c = SVector{2}(a+b)/2
	d = √(distance²(a,b))
	u = quarterturn(a-b)/(d)
	return Separator(c, u, zero(u), d/2)
end

function Separator(a::AbstractVector, l::Line; k=1)
	t0 = evaluate(l, a)/2
	if iszero(t0)
	# the degenerate case (t0==0) also makes sense (segment/endpoint boundary)
		return Separator(SVector{2}(a), zero(direction(l)), l.normal, zero(a[1]))
	end
	c = a-t0*l.normal
	# if t0 > 0: a is on the right side of l; parabola is oriented as l
	# if t0 < 0: a is on left side, parabola is oriented opposite to l
	return Separator(c, 2k*sign(t0)*√(abs(t0))*direction(l), sign(t0)*l.normal,
		abs(t0))
end
@inline Separator(l::Line, a::AbstractVector) = Separator(a, l; k=-1)

function Separator(l1::Line, l2::Line)
	d = det2(l1.normal, l2.normal)
	if iszero(d) # special case: two parallel lines
		c = (point(l1) + point(l2))/2
		return Separator(c, direction(l1), NaN16*l1.normal, abs(evaluate(l1, c)))
	end
	c = intersection(l1, l2)
	w = (l1.normal + l2.normal); w/= dot(w, l1.normal)
	return Separator(c, w, zero(w), zero(w[1]))
end

# Delaunay triangulation ««1
# Cell location functions ««2
"""
    findcell(cornertable, points, point)

In a Voronoi diagram,
finds the index of the cell closest to the given point.
"""
function findcell(t::AbstractTriangulation{J}, points, point) where{J}
	# Returns the cell index closest to this point
	c0 = tail(t, Arrow(1)); p0 = points[int(c0)]
	while true
		c1 = c0
		for s in neighbours(t, c0)
			p = points[int(s)]
			iscloser(point, p0, p) && continue
			c0, p0 = s, p
			break
		end
		c0 == c1 && return c0
	end
end

"""
    findnode(cornertable, points, point)

In a triangulation, returns the index of the triangle containing the point.
"""
function findnode(t::AbstractTriangulation{J}, points, i) where{J}
	q = rand(eachnode(t))
	point = points[i]
	c = 0
	while true
		c+= 1; @assert c ≤ 1e3
		# this guarantees that v[i] = tail(side(q, i)) for i in 1:3:
		v = [points[int(tail(t, e))] for e in arrows(q)]
		isleft(v[1], point, v[2]) && (q = adjnode(t, q, 1); continue)
		isleft(v[2], point, v[3]) && (q = adjnode(t, q, 2); continue)
		isleft(v[3], point, v[1]) && (q = adjnode(t, q, 3); continue)
		return q
	end
end
# Tree marking ««2
"""
    meetscircle(t, q, points, idx...)

Returns `true` iff site `points[idx]` meets circumcircle of node `q`.
"""
function meetscircle(t::AbstractTriangulation, q::Node, points::AbstractVector, i)#««
	c1, c2, c3 = Int.(triangle(t, q))
	return isincircle(points[SA[c1,c2,c3, i]]...)
end#»»
function badnodes(t::AbstractTriangulation, points, i)#««
	stack = [findnode(t, points, i)]
	tree = empty(stack)
	while !isempty(stack)
		q = pop!(stack)
		isone(int(q)) && continue # this is the reversed triangle
		q ∈ tree && continue
		meetscircle(t, q, points, i) || continue
		push!(tree, q)
		push!(stack, adjnodes(t, q)...)
	end
	return tree # the tree is usually small, no need to sort it
end#»»
function tree_boundary(t::AbstractTriangulation{J}, tree, doubleedges=()) where{J}#««
	# returns the list of half-edges running *outside* the tree,
	# cyclically ordered around the tree (random start).
	boundary = sizehint!(Arrow{J}[], length(tree)+2)
	@assert !isempty(tree)
	c = 0
	for e in arrows(first(tree))
		stack = [e]
		while !isempty(stack)
			c+=1; @assert c ≤ 1e2
			e = pop!(stack)
			o = opposite(t, e)
			if e ∈ doubleedges
				push!(boundary, o)
			elseif node(o) ∈ tree
				push!(stack, prev(o), next(o))
			else
				push!(boundary, o)
			end
		end
	end
	return boundary
	@assert false
end#»»
function star!(t::AbstractTriangulation{J}, s, tree, doubleedges=()) where{J}#««
	# replaces all nodes in `tree` by a star shape around `cell`,
	# adding two new nodes in the process.
# 	println("\e[1mstar!(", s, ", ", tree, ", double=", doubleedges, ")\e[m")
	boundary = tree_boundary(t, tree, doubleedges)
	push!(tree, newnodes!(t, 2)...)
	# rewrite boundary edges for double-edges, by substituting the
	# appropriate edge of the renamed triangle
	c2list = [ tail(t, e) for e in boundary ]
	j = 0
	for (i, o) in pairs(boundary)
		o ∈ doubleedges || continue
		if iszero(j)
			j = i
		else
			qi, qj = tree[i], tree[j]
			ei, ej = side(qi,1), side(qj, 1)
			boundary[i] = ej; boundary[j] = ei
			j = 0
		end
	end
	c1 = last(c2list)
	n = side(last(tree), 2)
	for (q, o, c2) in zip(tree, boundary, c2list)
		# set node `q` to triangle
		#      s
		#    ╱p n╲
		#   ╱  e  ╲
		#  c₁——————c₂
		#      o
# 		println("($q, $o, edge is $c1--$c2):")
		e = side(q, 1)
		p = side(q, 3)
		opposites!(t, e=>o, p=>n)
		n = side(q, 2)
		arrowsfrom!(t, s=>p, c1=>e, c2=>n)
		c1 = c2
	end
end#»»
# Triangulation ««2
"""
    triangulate(points, segments)

Returns a `CornerTable` whose cells are the Voronoi diagram of the input sites
and whose nodes are its Delaunay triangulation.

The sites are encoded as follows:
 - single points (given by their coordinates as a vector);
 - open segments (given by indices of their two end points).
"""
@inline triangulate(points; trim=true) = CornerTable{Int32}(points; trim)
function CornerTable{J}(points; trim=true) where{J}#««
	np = J(length(points))
# 	np < 3 && return CornerTable{J}(undef, 0)
	t = CornerTable{J}(J[4,6,5,1,3,2],J(np).+J[2,1,3,1,2,3],
		Vector{J}(undef, np+3))
	anyarrow!(t, Cell(np+1), Arrow(J(1)))
	anyarrow!(t, Cell(np+2), Arrow(J(2)))
	anyarrow!(t, Cell(np+3), Arrow(J(3)))

	m = maximum(abs(x) for p in points for x in p)
	push!(points, [0,-3m], [3m,2m], [-3m,2m])
  #»»
	# incrementally add all points ««
	for i in 1:np # Random.randperm(np)
		tree = badnodes(t, points, i)
		star!(t, Cell(i), tree)
	end
	#»»
	if trim
	# remove all superfluous nodes & cells ««
	# the nodes are sorted this way:
	# - inner nodes
	# - convex hull
	# - 1 backwards outer node
	k = nnodes(t)
	swapnodes!(t, Node(1), Node(k))
	k-= 1
	fakecells = np+1:np+3
	for i in nnodes(t)-1:-1:1; q = Node{J}(i)
		w = Int.(triangle(t,q))
		any(>(np), Int.(triangle(t, q))) || continue
		swapnodes!(t, q, Node(k))
		k-= 1
	end
	# »»
		resize!(points, np)
	end
	return t
end#»»
# Voronoi diagram ««1
# Data structure ««2
abstract type AbstractVoronoi{J} <: AbstractTriangulation{J} end
struct VoronoiDiagram{J,P,T} <: AbstractVoronoi{J}
	# J: integer index type
	# P: geometric point type
	# T: real distance type
	triangulation::CornerTable{J}
	points::Vector{P}
	segments::Vector{NTuple{2,J}} # indices into points
	geomnode::Vector{P}
	noderadius::Vector{T}
end

# @inline function VoronoiDiagram(t::CornerTable{J},
# 		points::AbstractVector{P}, segments::AbstractVector) where{J,P}
# 	v = VoronoiDiagram{J,P,eltype(P)}(t, points, NTuple{2,J}.(segments),
# 		Vector{P}(undef, nnodes(t)), vector{eltype(P)}(undef, nnodes(t))
# 	for q in allnodes(t); geometricnode!(v, q); end
# 	return v
# end
# 	edge = Vector{J}(undef, narrows(t))
# 	ne = narrows(t) ÷ J(2)
# 	firstarrow = sizehint!(J[], ne)
# 	for a in eacharrow(t)
# 		b = opposite(t, a)
# 		(b < a) && continue
# 		push!(firstarrow, a)
# 		edge[a] = edge[b] = length(firstarrow)
# 	end
# 	separator = [ Separator(points[int(right(t,Arrow(a)))],
# 		points[int(left(t,Arrow(a)))]) for a in firstarrow ]
# # 	noderadius = [ circumradius(points[SA[int.(triangle(t,q))...]]...)
# # 		for q in eachnode(t)]
# 	noderadius = sizehint!(eltype(P)[], nnodes(t))
# 	arrowtype = [ UInt8(2) for _ in 1:narrows(t) ]
# 	for q in eachnode(t)
# 		tri = triangle(t,q)
# 		a,b,c = points[tri[1]], points[tri[2]], points[tri[3]]
# 		ab, bc, ca = b-a, c-b, a-c
# 		push!(noderadius, circumradius(a,b,c))
# 		# if there is an obtuse angle, the corresponding edge
# 		# has the shape •→→
# 		if     dot(ab, ca) > 0 # a is obtuse
# 			arrowtype[3*int(q)-2] = UInt8(1)
# 		elseif dot(bc, ab) > 0 # b is obtuse
# 			arrowtype[3*int(q)-1] = UInt8(1)
# 		elseif dot(ca, bc) > 0 # c is obtuse
# 			arrowtype[3*int(q)-0] = UInt8(1)
# 		else  # all angles are acute; 1 is already the correct value here
# 		end
# 		
# 	end
# 	return VoronoiDiagram{J,P,eltype(P)}(t, firstarrow, edge,
# 		points, segments, noderadius, separator, arrowtype)
# end
@inline CornerTables.triangulation(v::VoronoiDiagram) = v.triangulation

@inline npoints(v::VoronoiDiagram) = length(v.points)
@inline nsegments(v::VoronoiDiagram) = length(v.segments)
@inline point(v::VoronoiDiagram, c::Cell) = v.points[int(c)]
@inline issegment(v::AbstractVoronoi, c::Cell) = int(c) > npoints(v)
@inline cellsegment(v::VoronoiDiagram, c::Cell) =
	Cell.(v.segments[int(c)-npoints(v)])
@inline geometricnode(v::VoronoiDiagram, q::Node) = v.geomnode[int(q)]
@inline noderadius(v::VoronoiDiagram, q::Node) = v.noderadius[int(q)]

function CornerTables.nnodes!(v::VoronoiDiagram, n)
	nnodes!(CornerTables.triangulation(v), n)
	resize!(v.geomnode, n)
	resize!(v.noderadius, n)
end
function CornerTables.swapnodes!(v::VoronoiDiagram, q1::Node, q2::Node)
	swapnodes!(CornerTables.triangulation(v), q1, q2)
	v.geomnode[SA[int(q1),int(q2)]] = v.geomnode[SA[int(q2),int(q1)]]
	v.noderadius[SA[int(q1),int(q2)]] = v.noderadius[SA[int(q2),int(q1)]]
end

# Geometric node computation ««2
@inline function geometricnode!(v::VoronoiDiagram, q::Node, g=equidistant(v,q))
# 	println("geometricnode!($q)/$(nnodes(v))")
	v.geomnode[int(q)] = g
	s = cell(v, q, 1)
	v.noderadius[int(q)] = if issegment(v, s)
		(c1, c2) = cellsegment(v, s)
		segdistance²(point(v, c1), point(v, c2), g)
	else
		distance²(point(v, s), g)
	end
end

function equidistant(v::AbstractVoronoi, q::Node)#««
	s = triangle(v, q)
	# orientation is important! some curves have two intersection points...
	if issegment(v, s[1])
		if issegment(v, s[2])
			issegment(v, s[3]) && return equidistant_sss(v, s[1], s[2], s[3])
			return equidistant_pss(v, s[3], s[1], s[2]) # points before segments
		else
			issegment(v, s[3]) && return equidistant_pss(v, s[2], s[3], s[1])
			return equidistant_pps(v, s[2], s[3], s[1])
		end
	else
		if issegment(v, s[2])
			issegment(v, s[3]) && return equidistant_pss(v, s[1], s[2], s[3])
			return equidistant_pps(v, s[3], s[1], s[2])
		else
			issegment(v, s[3]) && return equidistant_pps(v, s[1], s[2], s[3])
			return equidistant_ppp(v, s[1], s[2], s[3])
		end
	end
end#»»
function equidistant_sss(v::AbstractVoronoi, c1, c2, c3)#««
	(i1,j1),(i2,j2),(i3,j3)=cellsegment(v,c1),cellsegment(v,c2),cellsegment(v,c3)
	a1,b1,a2,b2,a3,b3 =
		point(v,i1),point(v,j1),point(v,i2),point(v,j2),point(v,i3),point(v,j3)
	u1,u2,u3 = b1-a1, b2-a2, b3-a3
	d1,d2,d3 = det2(u2,u3), det2(u3,u1), det2(u1,u2)

	iszero(d1) && return equidistant_sss_parallel(a2, u2, a3, a1, u1)
	iszero(d2) && return equidistant_sss_parallel(a3, u3, a1, a2, u2)
	iszero(d3) && return equidistant_sss_parallel(a1, u1, a2, a3, u3)

	return equidistant_sss(a1,b1,a2,b2,a3,b3)
# 	p1 = lineinter(a2,b2, a3,b3)
# 	p2 = lineinter(a3,b3, a1,b1)
# 	p3 = lineinter(a1,b1, a2,b2)
# # 	a = lineinter(v,i2,j2,i3,j3)
# # 	b = lineinter(v,i3,j3,i1,j1)
# # 	c = lineinter(v,i1,j1,i2,j2)
# 	@assert p1 ≠ nothing
# 	@assert p2 ≠ nothing
# 	@assert p3 ≠ nothing
# # 	a == nothing && return incenter_parallel(v, b, c, i2, j2, i3, j3)
# # 	b == nothing && return incenter_parallel(v, c, a, i3, j3, i1, j1)
# # 	c == nothing && return incenter_parallel(v, a, b, i1, j1, i2, j2)
# 	return incenter(p1,p2,p3)
end#»»

function lineinter(v,i1,j1,i2,j2)#««
	i1 ∈ (i2,j2) && return point(v,i1)
	j1 ∈ (i2,j2) && return point(v,j1)
	return lineinter(point(v,i1), point(v,j1), point(v,i2), point(v,j2))
end#»»
@inline equidistant_ppp(v::AbstractVoronoi, c1, c2, c3) =
	# all three cells are single points
	circumcenter(point(v,c1), point(v,c2), point(v,c3))

function equidistant_pps(v::AbstractVoronoi, c1, c2, c3)#««
	(i,j) = cellsegment(v, c3)
	# treat cases where one of the points is one end of the segment
	c1 == i && return equidistant_pxs(v, i, j, c2)
	c2 == i && return equidistant_pxs(v, i, j, c1)
	c1 == j && return equidistant_pxs(v, j, i, c2)
	c2 == j && return equidistant_pxs(v, j, i, c1)
	return equidistant_pps(point(v, c1), point(v, c2), point(v, i), point(v, j))
end#»»
function equidistant_pxs(v::AbstractVoronoi, c1, c2, c3)#««
	# node equidistant from: segment (c1,c2), points c1 & c3
	a, b, c = point(v,c1), point(v,c2), point(v,c3)
	ab, ac = b-a, c-a
	t = norm²(ac)/(2det2(ab,ac))
	return a+t*quarterturn(ab)
end#»»
function equidistant_pss(v::AbstractVoronoi, s, c2,c3)#««
	(i1,j1),(i2,j2) = cellsegment(v,c2), cellsegment(v, c3)
	# equidistant from s, (i1,j1), (i2,j2)
	s == i1 && return equidistant_pxs(v, i1, j1, i2, j2, 1)
	s == j1 && return equidistant_pxs(v, j1, i1, i2, j2, 1)
	s == i2 && return equidistant_pxs(v, i2, j2, i1, j1, -1)
	s == j2 && return equidistant_pxs(v, j2, i2, i1, j1, -1)
	return equidistant_pss(point(v, s), point(v,i1), point(v,j1),
		point(v, i2), point(v, j2))
end#»»
function equidistant_pxs(v::AbstractVoronoi, i1, j1, i2, j2, ε)#««
	i1 ∈ (i2,j2) && return point(v, i1)
	return equidistant_pxs(point(v,i1), point(v,j1), point(v,i2), point(v,j2), ε)
end#»»

function separator(v::AbstractVoronoi, c1, c2)
	a1 = issegment(v, c1) ? let (i1, j1) = cellsegment(v, c1)
		Line(point(v, i1), point(v, j1)) end : point(v, c1)
	a2 = issegment(v, c2) ? let (i2, j2) = cellsegment(v, c2)
		Line(point(v, i2), point(v, j2)) end : point(v, c2)
	return Separator(a1, a2)
end

# Cell location functions ««2
function findrootnode(v::AbstractVoronoi, i,j)
	# here the table is assumed built (for points)
	# so we can search in the nodes around the cell for point a
	# which one is closest to segment ab
# 	t = triangulation(v)
	a, b = point(v,i), point(v,j)
	emin, dmin = nothing, nothing
	global V=v
	for e in star(v,i)
		influences(v,i,j,node(e)) || continue
		d = segdistance²(a,b,geometricnode(v, node(e)))
		(emin == nothing || d < dmin) && ((emin, dmin) = (e, d))
	end
	@assert emin ≠ nothing
	return node(emin)
end

function meetscircle(v::AbstractVoronoi, q::Node, i, j)
	g, r = geometricnode(v, q), noderadius(v, q)
	a, b = point(v, i), point(v, j)
	ab2 = norm²(b-a)
	d = dot(g-a,b-a)
	return segdistance²(a, b, g) < r
end
function influences(v::AbstractVoronoi, i, j, q)
	a,b,g = point(v, i), point(v,j), geometricnode(v,q)
	ab, ag = b-a, g-a
	return (0 < dot(ab,ag) < dot(ab,ab))
 # (det2(ab,ag) ≥ 0) &&
end

# Computing Voronoi diagram ««1
# Finding double edges««2
# Double edge for a site s' = an edge between two nodes that are captured
# by s', but which contains at least one non-captured point

"""
edgeexit(v, e, i, j)
given the edge e and a new segment (i,j),
returns the minimum for z∈e of d²(z,s1)-d²(z,(ij))
(s1 being one of the two sites defining e)
"""
function edgeexit(v::AbstractVoronoi, e, i, j)
	c1, c2 = left(v,e), right(v,e)
	c2 = right(v, e)
	g1 = geometricnode(v, node(e))
	g2 = geometricnode(v, node(opposite(v,e)))
	g1 == g2 && return 0 # catch non-ternary nodes
	p,q = point(v,i), point(v,j)
	if issegment(v, c1)
		if issegment(v, c2)
			(i1, j1) = cellsegment(v, c1)
			(i2, j2) = cellsegment(v, c2)
			a1, b1, a2, b2 = point(v,i1), point(v,j1), point(v,i2), point(v,j2)
			return edgeexit_ss(a1,b1,a2,b2,g1,g2, p,q)
		end
		c1,c2 = c2,c1 # this replaces _sp case by equivalent _ps case just below:
	end
	if issegment(v, c2) # _ps case
		(i2, j2) = cellsegment(v, c2)
		a, a2, b2 = point(v,c1), point(v,i2), point(v,j2)
		return edgeexit_ps(a, a2,b2,g1,g2,p,q)
	else
		a, b = point(v,c1), point(v,c2)
		return edgeexit_pp(a,b,g1,g2, p,q)
	end
end

"returns min on [0,1] of a*x^2+b*x+c"
@inline function quadratic_min01(v)
	a,b,c = v
	a ≤ 0 && return min(c, a+b+c)
	(b ≥ 0) && return c
	(b ≤ -2a) && return (a+b+c)
	return (4*a*c-b*b)/(4*a)
end
	
function edgeexit_pp(a,b,g1,g2,p,q)
	# on the segment (g1,g2) of the mediator of (a,b),
	# compute min(d²(z,a)-d²(z,pq))
	
	# z = g1+t(g2-g1) = g1 + tv
	# d²(z,a) = ‖g1-a+tv‖² = ‖v‖² t² + 2(g1-a)⋅v t + ‖g1-a‖²
	# d²(z,pq) = <pqz>²/pq² = (<pqg1> + t <pqv>)²/pq²
	#  = <pqv>²/pq² t² + 2 <pqg₁><pqv>/pq² t² + <pqg₁>²/pq²
	v = g2-g1
	ag1 = g1 - a
	pq = q-p
	pqv = det2(pq,v)
	pqg1 = det2(pq,g1-p)
	da_quadratic = SA[norm²(v), 2*dot(ag1,v), norm²(ag1)]
	dpq_quadratic = SA[pqv^2, 2*pqg1*pqv, pqg1^2]/norm²(pq)
	return quadratic_min01(da_quadratic - dpq_quadratic)
end

function edgeexit_ss(a,b, a2,b2, g1,g2, p,q)
	# on the segment (g1,g2) of the mediator of ((ab), (a2b2)),
	# compute min(d²(z,(ab)) - d²(z,(pq)))
	
	# z = g1 + tv just as in the pp case
	# d²(z,(ab)) = <abz>²/ab² = (<abg>+t<abv>)²/ab²
	# d²(z,(pq)) = <pqz>²/pq² = (<pqg>+t<pqv>)²/pq²
	v = g2-g1
	ag, ab = g1-a, b-a
	pg, pq = g1-p, q-p
	abg, abv, ab2 = det2(ab, ag), det2(ab, v), norm²(ab)
	pqg, pqv, pq2 = det2(pq, pg), det2(pq, v), norm²(pq)
	dab_quadratic = SA[abv^2, 2*abv*abg, abg^2]/ab2
	dpq_quadratic = SA[pqv^2, 2*pqv*pqg, pqg^2]/pq2
	return quadratic_min01(dab_quadratic - dpq_quadratic)
end

function edgeexit_ps(a, b,c, g1,g2, p,q)
	# the intrinsic formulas here are a bit large...
	# we isometrically transform everything to a standard position instead,
	# by the transformation z ↦ (u⋅az, <u,bz>)
	# This maps the (bc) line to (y=0) and a to (0,d(bc,a))
	bc = c-b; u = bc/√(norm²(bc))
	xg1 = dot(u, g1-a) # no need to compute y for g1 and g2,
	xg2 = dot(u, g2-a) # which lie on the parabola
	xg1, xg2 = minmax(xg1, xg2)
	ya = det2(u,a-b)
	# now the parabola equation is 2y = ya + x^2/ya; x ∈ [xg1, xg2]
	# so that d²(z,bc) = y² = (ya+x²/ya)²/4
	fpx, fpy = (dot(u,p-a), det2(u,p-b))
	fqx, fqy = (dot(u,q-a), det2(u,q-b))
	pq2 = distance²(p,q)
	# now return the minimum for x ∈ [xg1, xg2] of
	# (ya + x²/ya)/2 - <fp, fq, z>/pq², where
	dpq = fpx*fqy-fqx*fpy
	dx = fpx-fqx
	dy = fpy-fqy
	#
	# This is a quartic in x:
	#= \\ gp verification script
F(x) = (ya+x^2/ya)^2/4;
G(x) = matdet([fpx,fpy,1;fqx,fqy,1;x,ya+x^2/ya,1])^2/pq^2;
H = (F(x)-G(x))*4*ya^2*pq^2;
dpq = fpx*fqy-fpy*fqx;
dx = fpx-fqx;
dy = fpy-fqy;

polcoeff(H,4,x)==pq^2-4*dx^2
polcoeff(H,3,x)==8*ya*dx*dy
polcoeff(H,2,x)==8*ya*dx*dpq+2*ya^2*(pq^2-4*dx^2-2*dy^2)
polcoeff(H,1,x)==8*ya^2*dy*(dx*ya-dpq)
polcoeff(H,0,x)==ya^4*pq^2-4*ya^2*(ya*dx-dpq)^2
	=#
	f = SA[pq2-4*dx^2, 8ya*dx*dy,
		8ya*dx*dpq+2*ya^2*(pq2-4*dx^2-2*dy^2),
		8ya^2*dy*(dx*ya-dpq),
		ya^4*pq2-4ya^2*(ya*dx-dpq)^2]/(4*pq2*ya^2)
	# this is guaranteed by geometry to have a unique minimum on ℝ
	# compute the root of the derivative by Newton's method (bounded iterations):
	r0 = 0.
	for _ in 1:10
		r1 = ((8*f[1]*r0+3*f[2])*r0^2-f[4])/(2*f[3]+r0*(6*f[2]+12*f[1]*r0))
		abs(r1 - r0) ≤ 1e-6*abs(r0) && break
		r0 = r1
	end
	(r0 < xg1) && (r0 = xg1)
	(r0 > xg2) && (r0 = xg2)
	return f[5]+r0*(f[4]+r0*(f[3]+r0*(f[2]+r0*f[1])))
end


# Tree traversal ««2
"""
    badnodes(v,i,j)

Returns a list of those nodes which are closer to segment (i,j)
than to their defining sites.

TODO: return a closed loop of edges which are fully closer to (i,j)
than to their defining sites, i.e. the boundary of the new cell for (i,j).
"""
function badnodes(v::AbstractVoronoi{J}, i, j) where{J} # segment case
	rootnode = findrootnode(v,i,j)
	# each entry in the stack is a pair (edge opposite which we are
	# exploring, cell around which we are turning):
	stack = [(e, Cell{J}(0)) for e in arrows(rootnode)]
	tree = [rootnode]
	loops = Cell{J}[]
	n = 0
	while !isempty(stack)
		n+=1
		@assert n ≤ 50
		(e, s) = pop!(stack)
		o = opposite(v, e)
		q = node(o)
		isone(int(q)) && continue # this is the reversed triangle
		if q ∈ tree # one loop
			!iszero(int(s)) && push!(loops, s)
			continue
		end
		(influences(v,i,j,q) && meetscircle(v, q, i, j)) || continue
		push!(tree, q)
		push!(stack, (prev(o), head(v,e)), (next(o), tail(v,e)))
	end
# 	@assert int(i) ≠ 5
	# break loops in tree by identifying double edges
# 	@assert isempty(loops)
	doubleedges = Arrow{J}[]
	for s in loops
		println("  breaking loop at $s")
		if s ∈ (i,j) # is this one of the ends of the segment we are inserting?
			u = point(v,j) - point(v,i); (s == j) && (u = -u)
			p = point(v,s)
			elist = collect(star(v,s))
			qlist = [ node(e) for e in elist ]
			k = det2(u, geometricnode(v, last(qlist))-p)
			for (e, q) in zip(elist, qlist)
				# q == right(e)
				k1 = k
				k = det2(u, geometricnode(v, q)-p)
				if (k1 > 0) && (k < 0)
					push!(doubleedges, prev(e), opposite(v, prev(e)))
					break
				end
			end
		else
			# s is a cell
			# look for an edge from s which (despite having its two end nodes
			# covered by the new segments) exits the cell of the new segment
			elist = collect(star(v, s))
			println("\e[32;1m elist=$elist\e[m")
			# groumpf
# 			e = argmin(edgeexit(v,e,i,j) for e in elist)
			e = first(elist); eemin = edgeexit(v,e,i,j)
			for e1 in elist[2:end]
				eemin1 = edgeexit(v, e1, i, j)
				(eemin1 < eemin) && ((e, eemin) = (e1, eemin1))
			end
			push!(doubleedges, e, opposite(v,e))
		end
	end
	@assert length(doubleedges) == 2*length(loops)
	return (tree, doubleedges) # the tree is usually small, no need to sort it
end
# end#»»
function star!(v::AbstractVoronoi, s, tree, doubleedges)#««
	star!(CornerTables.triangulation(v), s, tree, doubleedges)
	global V = v
	resize!(v.geomnode, length(v.geomnode)+2)
	resize!(v.noderadius, length(v.geomnode))
	for q in tree
		geometricnode!(v, q)
	end
end#»»
# Compute Voronoi diagram ««2
@inline VoronoiDiagram(points, segments=[]) =
	VoronoiDiagram{Int32}(points, segments)
# @inline function VoronoiDiagram(t::CornerTable{J},
# 		points::AbstractVector{P}, segments::AbstractVector) where{J,P}
# 	v = VoronoiDiagram{J,P,eltype(P)}(t, points, NTuple{2,J}.(segments),
# 		Vector{P}(undef, nnodes(t)), vector{eltype(P)}(undef, nnodes(t))
# 	for q in allnodes(t); geometricnode!(v, q); end
# 	return v
# end
function VoronoiDiagram{J}(points::AbstractVector{P}, segments) where{J,P}#««
	Random.seed!(0)
	np = length(points)
	ns= length(segments)
	ntotal = np+ ns
	points = SVector{2,Float64}.(points)
	println("generating corner table...")
	t = CornerTable{J}(points; trim=false)
	showall(stdout, t)
	println("expanding corner table for segments...")
	ncells!(t, ncells(t) + ns)
	v = VoronoiDiagram{J,P,eltype(P)}(t, points, NTuple{2,J}.(segments),
		similar(points, nnodes(t)), Vector{eltype(P)}(undef, nnodes(t)))
	for q in eachnode(t); geometricnode!(v, q); end

	# incrementally add all segments ««
	ns = length(segments)
	for i in 1:ns # Random.randperm(ns)
		s = Cell(np+3+i); (a,b) = cellsegment(v, s)
		println("\e[31;1;7minserting segment $s = ($a,$b)\e[m")
		tree, doubleedges = badnodes(v, a, b)
		star!(v, s, tree, doubleedges)
		showcell(stdout, v, s)
	end
	#»»
# 	# for each segment, identify the delimitation between both sides of the
# 	# cell (this will be useful for offsetting) #««
# 	# i.e. the first node 
# 	for i in 1:ns
# 		s = Cell(np+3+i); (a,b) = cellsegment(v,s)
# 		# 
# 		println("looking for first node for $s = ($a,$b)...")
# 		for e in star(v, s)
# 			s1 = left(v, e)
# 			println("  e = $e ($(tail(v,e)), $(right(v,e)), $(left(v,e)))")
# 			if (s1 == a) || (issegment(v, s1) && a ∈ cellsegment(v, s1))
# 				println("  cell $s1 is attached to $a")
# 			end
# 		end
# 	end
# 	#»»
	# remove all superfluous nodes & cells ««
	# the nodes are sorted this way:
	# - inner nodes
	# - convex hull
	# - 1 backwards outer node
	k = nnodes(t)
	swapnodes!(v, Node(1), Node(k))
	k-= 1
	fakecells = np+1:np+3
	for i in nnodes(t)-1:-1:1; q = Node{J}(i)
		w = Int.(triangle(t,q))
		w[1] ∈ fakecells || w[2] ∈ fakecells || w[3] ∈ fakecells || continue
# 		any(>(ntotal), Int.(triangle(t, q))) || continue
		swapnodes!(v, q, Node(k))
		k-= 1
	end
# 	nnodes!(t, k)
	# »»
# 	resize!(points, np)
	return v
end#»»
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
		showcell(stdout, v, s12)
		e2 = anyarrow(v, s12)
		while head(v, e2) ≠ c2; e2 = after(v, e2); end
		e1 = e2
		while head(v, e1) ≠ c1; tail!(v, e1, s21); e1 = after(v, e1); end
		# split the cell by inserting two new nodes
		# (= split vertex s12 of the dual triangulation)
		o1, o2 = opposite(v, e1), opposite(v, e2)
		q1, q2 = newnodes!(v, 2)
		q11, q12, q13 = side(q1,1), side(q1,2), side(q1,3)
		q21, q22, q23 = side(q2,1), side(q2,2), side(q2,3)
		tail!(v, q11=>s12, q12=>s21, q13=>c1, q21=>s21, q22=>s12, q23=>c2)
		opposites!(v, q11=>q21, q12=>o1, q13=>e1, q22=>o2, q23=>e2)
# 		pe1, pe2,  = prev(e1), prev(e2)
# 		ope1, ope2 = opposite(v, prev(e1)), opposite(v, prev(e2))
# 		tail!(v,
# 			side(q1,1)=>c1, side(q1,2)=>s12, side(q1,3)=>s21,
# 			side(q2,1)=>c2, side(q2,2)=>s21, side(q2,3)=>s12)
# 		opposites!(v, side(q1,1) => side(q2,1),
# 			pe1=>side(q1,3), ope1=>side(q1,2),
# 			pe2=>side(q2,3), ope2=>side(q2,2))
		geometricnode!(v, q1, p1)
		geometricnode!(v, q2, p2)
# 		shownode(stdout, v, q1)
# 		shownode(stdout, v, q2)
		anyarrow!(v, s12, q11); anyarrow!(v, s21, q21)
# 		showcell(stdout, v, s12)
# 		showcell(stdout, v, s21)
	end
	return v
end#»»

# Data type ««2
struct OffsetDiagram{J,P,T} <: AbstractVoronoi{J}
	voronoi::VoronoiDiagram{J,P,T}
	# arrows ↔ edges:
	edge::Vector{J}       # arrow → edge
	firstarrow::Vector{J} # edge → arrow
	# parametrizations of edges by distance to sites:
	separator::Vector{Separator{T}} # indexed by edges
	branch::Vector{Int8} # indexed by arrows;
	# this encodes the branch of the separator on which the node lies,
	# as either +1, -1, or 0 (for a parallel bisector)
end

@inline voronoi(v::OffsetDiagram) = v.voronoi
for f in (:(CornerTables.triangulation), :npoints, :nsegments, :point,
	:cellsegment, :geometricnode, :noderadius)
	@eval @inline $f(t::OffsetDiagram, args...;kwargs...) =
		$f(voronoi(t), args...;kwargs...)
end
@inline branch(v::OffsetDiagram, e::Arrow) = v.branch[int(e)]
@inline edge(v::OffsetDiagram, e::Arrow) = Edge(v.edge[int(e)])
@inline firstarrow(v::OffsetDiagram, e::Edge) = Arrow(v.firstarrow[int(e)])
@inline separator(v::OffsetDiagram, e::Edge) = v.separator[int(e)]


@inline OffsetDiagram(points, segments) = OffsetDiagram{Int32}(points, segments)
function OffsetDiagram{J}(points::AbstractVector{P}, segments) where{J,P}
	v = VoronoiDiagram{J}(points, segments)
	println("\e[1;7m in OffsetDiagram, got the following Voronoi diagram:\e[m")
	showall(stdout, v)
	splitsegments!(v)
	println("\e[1;7m after split segments!:\e[m")
	showall(stdout, v)
	gnuplot(v)
	ne = narrows(v) ÷ J(2)
	edge = Vector{J}(undef, narrows(v))
	sep = sizehint!(Separator{eltype(P)}[], ne)
	firstarrow = sizehint!(J[], ne)
	for a in eacharrow(v)
		b = opposite(v, a)
		b < a && continue
		push!(firstarrow, a)
		edge[a] = edge[b] = length(firstarrow)
		println("separator($(length(firstarrow))) = $(right(v,a))/ $(left(v,a))")
		push!(sep, separator(v, right(v, a), left(v, a)))
	end
	branch = sizehint!(Int8[], narrows(v))
	for q in eachnode(v), i in 1:3
		a = side(q, i)
		ea = edge[int(a)]
		s = sep[ea]
		println("computing geometric branch for arrow $a:$(tail(v,a))->$(head(v,a)), edge $ea, separator $sep")
		g, t = geometricnode(v, q), √(noderadius(v, q))
		orient = (firstarrow[ea] == int(a))
		# if orient == true then the + branch lies at the head of a, - at its tail
		sgplus, sgminus = evaluate(s, t, +1), evaluate(s, t, -1)
		# fun fact: there are exactly three “inverted” arrows (i.e. signs (-1,1));
		# these are the three arrows of the outer triangle.
		if t == s.t0; push!(branch, Int8(0))
		elseif any(isnan, s.normal); push!(branch, Int8(0))
		elseif sgplus  ≈ g; push!(branch, Int8(1))
		elseif sgminus ≈ g; push!(branch, Int8(-1))
		else
			println("""
\e[7m arrow $a (edge=$(edge[int(a)]), $(left(v,a))/$(right(v,a))):\e[m
  node $q = $g, r=$t
	separator $s
	candidates $sgplus, $sgminus
\e[31;1m no good candidate!\e[m""")
		end
	end
	for a in eacharrow(v)
		# fix branches of type (1,0) and (-1,0); replace by (1,1) and (-1,-1)
		iszero(branch[int(a)]) && (branch[int(a)] = branch[int(opposite(v, a))])
	end

	return OffsetDiagram{J,P,eltype(P)}(v, edge, firstarrow, sep, branch)
end
# Single cell offset ««2

function cut_parabola(a, b, c, g1, g2, r)
	# given a parabola (focus a, directrix (bc)),
	# returns the two points (if any) at distance²=r from the focus
	# (in order along the direction bc), or `nothing`.
	# let a' = proj(a), z = a'+(x+iy) * bc
	# the parabola equation is y = h/2 + x²/(2h), where h=d(a,bc)/bc
	# smallest y² value is y=h/2 for x=0
	# if r² > h²/4 then there are two solutions, both of them with y=r
	# or x = ±√(h(2r-h))
	bc, ba = c-b, a-b; ibc = quarterturn(bc)
	bc2 = norm²(bc); bc1 = √(bc2)
	h = det2(bc, ba)/bc2 # (signed) distance from bc to a
	a1 = a - h*ibc # orthogonal projection of a on bc
	ρ = r/bc1
	D = h*(2ρ-h)
	D ≤ 0 && return nothing
	s = sqrt(D)
	dot(g1-a, bc) < -s*bc2 || return nothing
	dot(g2-a, bc) >  s*bc2 || return nothing
	z1, z2 = (a1-s*bc+ρ*ibc, a1+s*bc+ρ*ibc)
	return (z1,z2)
end
@inline cut_parabola(v::AbstractVoronoi, a::Cell, b::Cell, c::Cell,
		q1::Node, q2::Node, r) =
	cut_parabola(point(v,a), point(v,b), point(v,c),
		geometricnode(v, q1), geometricnode(v, q2), r)

function edgeoffset(v::OffsetDiagram, e::Arrow, radius)
	# given an edge bounding a cell,
	# prints all the points on this edge lying at distance `r` from the cell
	# (edges are open at their far end and closed at their near end)
	r = radius^2
	o = opposite(v, e)
	ee = edge(v, e)
	orient = (e == firstarrow(v, ee))
	b0, b1 = branch(v, e), branch(v, o)
	b = b0 + 2b1

	sep = separator(v, edge(v, e))
	q0, q1 = node(e), node(o)
	r0, r1 = noderadius(v, q0), noderadius(v, q1)
	println("  edge $e (opp. $o; orient=$orient, branches=$b0, $b1): $q0($r0) -> $q1($r1)")
# 	if r0 ≤ r < r1
# 		println("intersection (IN) at $(evaluate(sep, radius, 1))")
# 	elseif r1 ≤ r < r0
# 		println("intersection (OUT) at $(evaluate(sep, radius, -1))")
# 	elseif (t0, t1) == (1, 1)
# 		println("monotonic edge, no intersection")
# 	elseif (t0, t1) == (1,2)
# 		println("  <-  -> ")
# 	elseif (t0, t1) == (2,1)
# 		println("  ->  <- ")
# 	elseif (t0, t1) == (0,0)
# 		println("parallel separator at distance $r0 == $r1")
# 	else
# 		error("bad edge types: ($t0, $t1)")
# 	end
end
function segmentoffset(v::OffsetDiagram, c::Cell, ::typeof(+), radius)
	# compute the boundary of the offset for a segment-cell (right side)
	# i.e. all intersection points of the edges of the cell and the offset line
	println(" walking the boundary of cell $c:")
	# is the previous node inside the offset line?
	is_in = true
	(a,b) = cellsegment(v, c)
	for e in ring(v, c)
		edgeoffset(v, e, radius)
	end
	for e in ring(v,c)
		c2 = right(v,e)
		o = opposite(v, e)
		q0, q1 = node(e), node(o)
		is_in1 = is_in
		is_in = noderadius(v, q1) ≤ radius^2
		println("edge $e ($(branch(v,e)), $(branch(v,o))):")
		println("  from $q0 to $q1; on our right is $c2; ($is_in1, $is_in)")
		println("  arrival $q1: $(geometricnode(v,q1)), d²=$(noderadius(v,q1))")
		if is_in1 ≠ is_in
			println("  \e[32;1m offset line crosses this edge!\e[m")
		elseif c2 ∈ (a,b) || issegment(v, c2)
			println("  this is a straight line, no crossing")
		else
			println("  this is a parabola arc! does it cross twice?")
			# compute points of intersection
			pts = cut_parabola(v, c2, a, b, q0, q1, -radius)
			if pts == nothing
				println("   no!")
			else
				println("   \e[32myes! at $pts\e[m")
			end
		end
	end
end
function offset(points, segments, radius)
	v = OffsetDiagram(points, segments)
	return offset(v, radius)
end
function offset(v::OffsetDiagram, radius)
	# do a walk on the set of segments by connectivity
	#   (we only start on segments but we can walk to points if an edge
	#   points us to them)
	# for each segment, determine the offset of the segment by the radius
	# and in particular all intersection points with the edges
	# then glue along the edges
	# as long as there are segments remaining...
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
	for e in eacharrow(v)
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
function gnuplot(v::AbstractVoronoi; scale=10., f_png="/tmp/a.png")
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
function shownode(io::IO, v::AbstractVoronoi, q::Node)
	println(io, "\e[33m", q, triangle(v,q), "\e[m: ",
		geometricnode(v, q), " r²=", noderadius(v,q))
# 	for i in (1,2,3)
# 		e = side(q, i); o = opposite(v, e); oo = opposite(v,o)
# 		oo ≠ e && println(io, "  \e[31;7m opposite($o) = $oo, should be $e\e[m")
# 	end
end
function showarrow(io::IO, v::OffsetDiagram, a::Arrow)
	if !iszero(branch(v, a))
		q = node(a)
		sep = separator(v, edge(v, a))
		g = geometricnode(v, q)
		r = √(noderadius(v, q))
		b = branch(v, a)
		h = evaluate(sep, r, b)
		g ≈ h ||
			println(io, "  \e[31;7mgeometricnode($q) = $g; evaluate($r, $b) = $h\e[m")
	end
end
function showedge(io::IO, v::OffsetDiagram, e::Edge)
	a = Arrow(v.firstarrow[int(e)])
	o = opposite(v, a)
	sep = separator(v, e)
	print(io, "\e[32m", e, " = (", a, ",", o, ")\e[m  ")
	print(io, left(v,a), "|", right(v,a), " ")
	print(io, node(a), "  (", branch(v, a), ") -> (",
	  branch(v, o), ") ", node(o), "  ",
		branch(v, a)+2*branch(v, o))
	println(io)
	ea, eo = Edge(v.edge[int(a)]), Edge(v.edge[int(o)])
	ea ≠ e && println(io, "  \e[31;7m edge($a) = $ea, should be $e\e[m")
	eo ≠ e && println(io, "  \e[31;7m edge($o) = $eo, should be $e\e[m")
	showarrow(io, v, a)
	showarrow(io, v, o)
end
function Base.display(io::IO, v::AbstractVoronoi)
	println(io, "\e[1m begin Voronoi diagram with $(nnodes(v)) nodes and $(ncells(v)) cells:\e[m")
	for q in eachnode(v); shownode(io, v, q); end
	for c in eachcell(v); showcell(io, v, c); end
# 	for e in 1:length(v.firstarrow); showedge(io, v, Edge(e)); end
	println(io, "\e[1m end Voronoi diagram\e[m")
end
function Base.display(io::IO, v::OffsetDiagram)
	println(io, "\e[1m begin offset diagram with $(nnodes(v)) nodes and $(ncells(v)) cells:\e[m")
	for q in eachnode(v); shownode(io, v, q); end
	for c in eachcell(v); showcell(io, v, c); end
	for e in 1:length(v.firstarrow); showedge(io, v, Edge(e)); end
	println(io, "\e[1m end offset diagram\e[m")
end
	
@inline Base.display(v::AbstractVoronoi) = display(stdout, v)
# end »»1
end

V=Voronoi
# using StaticArrays
# TODO: tests
# V = Voronoi
# t=V.triangulate([[-10,0],[10,0],[0,10.]])
# v=V.voronoi([(0,0),(10,0),(11,3),(6,2),(5,9),(1,10)],[(6,1),(1,2),(2,3),(3,4),(4,5),(5,6)])
# v=V.VoronoiDiagram([[0.,0],[10,0],[0,10],[10,10],[5,9],[5,1]],[(3,4),(5,6)])

v=V.OffsetDiagram([[0.,0],[10,0],[0,10],[10,10],[5,9],[5,1]],[(1,2),(2,6),(6,5),(5,4),(4,3),(3,1)])
# o = V.segmentoffset(v, V.Cell(10), +, .5)

# m = V.edgeexit_pp([5,9],[5,1],[1.6,5],[8.4,5],[0,10],[10,10])
# V.equidistant_pss([0,0],[5,9],[5,1],[0,10],[10,10]) # [5,5]

# p=NTuple{2,Float64}[(0,0),(10,0),(11,3),(6,2),(5,9),(1,10)]
# V.offset(p,[(6,1),(1,2),(2,3),(3,4),(4,5),(5,6)],0.5)

# collect(V.eachnode(v))
