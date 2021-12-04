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

import Base: Int

# Tools ««1
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
	b,c = b-a,c-a
	m = SA[norm²(b) b[1] b[2];norm²(c) c[1] c[2]]
	kn = det(m[:,SA[2,3]])
	kx = det(m[:,SA[1,3]])/(2kn)
	ky = det(m[:,SA[2,1]])/(2kn)
	return a + SA[kx, ky]
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
	@assert pqa ≠ 0
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
	if pqrs == 0 # special case: parallel lines
		pqa = det2(p,q,a)
		pqa < 0 && ((p,q, pqa) = (q,p, -pqa)) # ensure orientation
		pq = q-p
		pq2 = norm²(pq)
# 	pqr, pqs = det2(p,q,r), det2(p,q,s)
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

# Bisectors ««2
"mediator line of segment (ab)"
function mediator(a,b)
	return SA[2*(a[1]-b[1]), 2*(a[2]-b[2]), norm²(b)-norm²(a)]
end

"intersection point of segments (ab) and (cd)"
function intersection(a,b,c,d)
end

"angle bisector of the non-crossing segments (ab) and (cd)"
function bisector(a,b,c,d)
	# since these segments are non-crossing, at least one of them lies
	# fully on one side of their intersection point.
end
	

# Named indices ««2
struct NamedIndex{S,T<:Signed} i::T; end
@inline (T::Type{<:Integer})(i::NamedIndex) = T(i.i)
# @inline Base.convert(T::Type{<:NamedIndex},i::Integer) = T(i)
# @inline Base.zero(T::Type{<:NamedIndex}) = T(0)
# @inline Base.zero(x::T) where{T<:NamedIndex} = zero(T)
@inline Base.show(io::IO, i::NamedIndex{S}) where{S} = print(io, S, Int(i))
@inline Base.:+(a::T,b::T) where{T<:NamedIndex} = T(Int(a)+Int(b))
@inline Base.:-(a::T,b::T) where{T<:NamedIndex} = T(Int(a)-Int(b))
@inline NamedIndex{S}(x) where{S} = NamedIndex{S,typeof(x)}(x)
@inline Base.isless(x::NamedIndex{S}, y::NamedIndex{S}) where{S} = Int(x)<Int(y)
@inline Base.:(==)(x::NamedIndex{S},y::NamedIndex{S}) where{S} =(x.i==y.i)
macro NamedIndex(name,S)
	quote
	const $(esc(name)) = NamedIndex{$(QuoteNode(S))}
	Base.show(io::IO, ::Type{<:$(esc(name))}) = print(io, $(string(name)))
end end
@NamedIndex Corner e
@NamedIndex Cell s
@NamedIndex Node q

# CornerTable ««1
# Nodes and corners ««2
@inline next3(i::J) where{J<:Integer} = iszero(i%3) ? i-2 : i+1
@inline prev3(i::J) where{J<:Integer} = isone(i%3) ? i+2 : i-1
@inline next(e::Corner{J}) where{J} = Corner{J}(next3(Int(e)))
@inline prev(e::Corner{J}) where{J} = Corner{J}(prev3(Int(e)))

@inline node(e::Corner{J}) where{J} = Node{J}(fld1(Int(e), 3))
@inline side(q::Node{J}, i) where{J} = Corner{J}(3*Int(q)-3+i)
@inline corners(q::Node{J}) where{J} = Corner{J}.(3*Int(q) .- (2,1,0))
# AbstractTriangulation««2
abstract type AbstractTriangulation{J} end
# Objects of this type must define the following methods:
#  - triangulation() - returning a CornerTable, or something with the
# following methods:
#  - _apex, _apex!
for f in (:apex, :apex!, :opposite, :opposite!, :nnodes, :nnodes!,
	:anycorner, :anycorner!, :ncells, :ncells!)
	@eval @inline $f(t::AbstractTriangulation, args...;kwargs...) =
		$f(triangulation(t), args...;kwargs...)
end

# Data type ««2

"""
    CornerTable - generic table for a triangulation (or its dual).

The triangulation is represented as corners.
Each corner uniquely represents a corner of one of the triangles
around a given site.

The table contains the following:
 - nodes (triangles of three consecutive corners in the table),
 - cells (duals of nodes, made of adjacent corners)

The triangulation has nodes as faces and cells as vertices;
the opposite of a corner is pictured as <corner | opposite >.
The dual Voronoi diagram has nodes as vertices and cells as faces;
here a corner represents a half-edge from one node to another.

This structure contains only combinatorial information.
Geometric information, when needed, is passed as a separate array of points.
"""
struct CornerTable{J <: Integer} <: AbstractTriangulation{J}
	opposite::Vector{J} # corner->corner
	apex::Vector{J}     # corner->cell
	anycorner::Vector{J}# cell -> corner
end

@inline CornerTable{J}(::UndefInitializer, ncells::Integer) where{J} =
	CornerTable{J}(J[], J[], zeros(J, ncells))
@inline CornerTable{J}() where{J} = CornerTable{J}(undef, 0)
@inline triangulation(t::CornerTable) = t

# Accessors ««2
# Corner functions ««3
@inline apex(t::CornerTable, e::Corner) = Cell(t.apex[Int(e)])
@inline apex!(t::CornerTable, e::Corner, c::Cell) = t.apex[Int(e)] = Int(c)

@inline opposite(t::CornerTable, e::Corner) = Corner(t.opposite[Int(e)])
@inline opposite!(t::CornerTable, e::Corner, x::Corner) =
	t.opposite[Int(e)] = Int(x)

# Node functions ««3
@inline nnodes(t::CornerTable{J}) where{J} = J(length(t.opposite)) ÷ J(3)
@inline nnodes!(t::CornerTable, n) = 
	(resize!(t.opposite, 3n); resize!(t.apex, 3n); t)

# Cell functions ««3
@inline anycorner(t::CornerTable, c::Cell) = Corner(t.anycorner[Int(c)])
@inline anycorner!(t::CornerTable, c::Cell, e::Corner) =
	t.anycorner[Int(c)] = Int(e)
@inline ncells(t::CornerTable{J}) where{J} = J(length(t.anycorner))
@inline ncells!(t::CornerTable, n) = resize!(t.anycorner, n)

# AbstractTriangulation methods ««2
@inline apex!(t::AbstractTriangulation, l::Pair{<:Corner, <:Cell}...) =
	for (e,s) in l; apex!(t, e, s); end
@inline opposites!(t::AbstractTriangulation, l::Pair{<:Corner,<:Corner}...) =
	for(e1, e2) in l; opposite!(t, e1, e2); opposite!(t, e2, e1); end
@inline anycorner!(t::AbstractTriangulation, l::Pair{<:Cell, <:Corner}...) =
	for (s,e) in l; anycorner!(t, s, e); end
@inline cellcorner!(t::AbstractTriangulation, l::Pair{<:Cell, <:Corner}...) =
	for(s,e) in l; anycorner!(t, s, e); apex!(t, e, s); end

@inline right(t::AbstractTriangulation, e::Corner) = apex(t, next(e))
@inline left(t::AbstractTriangulation, e::Corner) = apex(t, prev(e))
@inline base(t::AbstractTriangulation, e::Corner) = (left(t,e),right(t,e))
@inline after(t::AbstractTriangulation, e::Corner) = next(opposite(t, next(e)))
@inline before(t::AbstractTriangulation, e::Corner) = prev(opposite(t, prev(e)))

@inline cell(t::AbstractTriangulation, q::Node, i) = apex(t, side(q,i))
@inline triangle(t::AbstractTriangulation, q::Node) =
	(cell(t,q,1), cell(t,q,2), cell(t,q,3))
@inline allcells(t::AbstractTriangulation) =
	(Cell(i) for i in Base.OneTo(ncells(t)))
@inline allnodes(t::AbstractTriangulation) =
	(Node(i) for i in Base.OneTo(nnodes(t)))
@inline alltriangles(t::AbstractTriangulation) =
	(triangle(t,q) for q in allnodes(t))
@inline adjnode(t::AbstractTriangulation, q::Node, i) =
	node(opposite(t, side(q, i)))
@inline adjnodes(t::AbstractTriangulation, q::Node) =
	(adjnode(t,q,1), adjnode(t,q,2), adjnode(t,q,3))

function newnodes!(t::AbstractTriangulation{J}, k) where{J}
	n = nnodes(t)
	nnodes!(t, n+k)
	return Node{J}.(n+1:n+k)
end

# Iterators ««2
struct CTIterator{S,X,J,T,A<:AbstractTriangulation}
	# S is an identifier
	# X is the return type
	# J is the index type
	# A, T are field types
	table::A
	start::T
	@inline CTIterator{S,X}(t::A, s::T) where{S,X,J,T,A<:AbstractTriangulation{J}}=
		new{S,X,J,T,A}(t,s)
end
@inline Base.eltype(::CTIterator{S,X}) where{S,X} = X
@inline Base.keys(g::Base.Generator{<:CTIterator}) = g.iter
@inline Base.keys(g::Base.Generator{<:Base.Generator{<:CTIterator}}) =
	g.iter.iter
@inline Base.IteratorSize(::CTIterator) = Base.SizeUnknown()
@inline Base.iterate(it::CTIterator) = (it.start, it.start)
@inline function Base.iterate(it::CTIterator, s)
	s = next(it, it.table, s)
	stop(it, it.table, s) && return nothing
	return (s, s)
end
@inline stop(it::CTIterator,_,e) = (e==it.start)
"""
    star(cornertable, cell)

Iterates all corners from the same cell.
"""
@inline star(t::AbstractTriangulation{J}, e::Corner) where{J} =
	CTIterator{star,Corner{J}}(t, e)
@inline star(t::AbstractTriangulation, c::Cell) = star(t, anycorner(t,c))
@inline next(::CTIterator{star}, t, e) = after(t, e)

# same as `star`, but in reverse order:
@inline revstar(t::AbstractTriangulation{J}, e::Corner) where{J} =
	CTIterator{revstar,Corner{J}}(t, e)
@inline revstar(t::AbstractTriangulation, c::Cell) = revstar(t, anycorner(t,c))
@inline next(::CTIterator{revstar},t,e) = before(t,e)

@inline Base.reverse(it::CTIterator{star,T}) where{T} =
	CTIterator{revstar,T}(it.table, it.start)
@inline Base.reverse(it::CTIterator{revstar,T}) where{T} =
	CTIterator{star,T}(it.table, it.start)

"""
    ring(cornertable, corner)

Iterates all corners left-adjacent to the indicated cell.
"""
@inline ring(t::AbstractTriangulation{J}, e::Corner) where{J} =
	CTIterator{ring,Corner{J}}(t, next(e))
@inline ring(t::AbstractTriangulation, c::Cell) = ring(t, anycorner(t, c))
@inline next(::CTIterator{ring}, t, e) = prev(opposite(t, e))
"""
   neighbours(cornertable, cell)

Iterates all cells adjacent to the indicated cell.
"""
@inline neighbours(t::AbstractTriangulation{J}, c::Cell) where{J} =
	(apex(t, x) for x in ring(t, anycorner(t, c)))

# Constructor from triangle list ««2
function CornerTable{J}(triangles) where{J}
	nt = length(triangles); nc = 3nt
	ns = 0
	slist = sizehint!(J[], nc)
	olist = Vector{J}(undef, nc)
	for q in triangles, s in q[1:3]
		(s > ns) && (ns = s)
		push!(slist, s)
	end
	clist = Vector{J}(undef, ns)
	# list of all corners seeing `s` on their right
	edge = [ sizehint!(J[], 6) for _ in 1:ns ]
# 	clist[slist] .= 1:nc
	for (e, s) in pairs(slist)
		clist[s] = e
		n = next3(e); sn = slist[n]
		push!(edge[sn], e)
	end
	for (e, s) in pairs(slist)
		sn, sp = slist[next3(e)], slist[prev3(e)]
		for e1 in edge[sp]
			slist[prev3(e1)] == sn || continue
			olist[e] = e1; olist[e1] = e
		end
	end
	return AbstractTriangulation{J}(olist, slist, clist)
end

# Elementary modifications: flip!, insert! ««2
"""
    flip!(cornertable, corner)

Flips the quadrilateral delimited by `corner` and its opposite.
Returns the two corners (in order) replacing `corner`.
"""
function flip!(t::AbstractTriangulation, e::Corner)
	n, p, o = next(e), prev(e), opposite(t, e)
	s, sn, so, sp = apex(t, e), apex(t, n), apex(t, o), apex(t, p)
	no, po = next(o), prev(o)
	# we leave (n ↔ on) and (no ↔ ono) edges unchanged:
	op, opo = opposite(t, p), opposite(t, po)
	opposites!(t, p=>po, e=>opo, o=>op)
	apex!(t, n=>so, no=>s, po=>sn)
	anycorner!(t, so=>n, s=>no, sn=>po)
	return (no, e)
end

"""
    insert!(cornertable, node, cell)

Inserts a new cell by splitting the given node.
Returns the three corners formed at the new cell.
"""
function insert!(t::AbstractTriangulation, q::Node, c::Cell)
	e0, n0, p0 = corners(q)
	c0, c1, c2 = apex(t, e0), apex(t, n0), apex(t, p0)
	# e0 and its opposite corner are unchanged
	on, op = opposite(t, n0), opposite(t, p0)
	q1, q2 = newnodes!(t, 2)
	e1, n1, p1 = corners(q1)
	e2, n2, p2 = corners(q2)
	opposites!(t, e1=>on, e2=>op, n0=>p1, n1=>p2, n2=>p0)
	apex!(t, e0=>c, e1=>c, e2=>c, n1=>c2, p1=>c0, n2=>c0, p2=>c1)
	anycorner!(t, c=>e0, c1=>n0, c2=>p0, c0=>n1)
	return (e0, e1, e2)
end

"""
    swapcells!(cornertable, cell, cell)

Swap cell indices for those two cells.
"""
function swapcells!(t::AbstractTriangulation, c1::Cell, c2::Cell)#««
	println("swapping cells $c1 and $c2")
	c1 == c2 && return
	for e in star(t, c1); apex!(t, e, c2); end
	for e in star(t, c2); apex!(t, e, c1); end
	e1 = anycorner(t, c1)
	anycorner!(t, c1, anycorner(t, c2))
	anycorner!(t, c2, e1)
end#»»
"moves cell c1 to (undefined) position c2, overwriting c2 in the process"
function movecell!(t::AbstractTriangulation, c1::Cell, c2::Cell)#««
	c1 == c2 && return
	for e in star(t, c1); apex!(t, e, c2); end
	anycorner!(t, c2, anycorner(t, c1))
end
function swapnodes!(t::AbstractTriangulation, q1::Node, q2::Node)#««
	q1 == q2 && return
	e1,n1,p1 = side(q1,1), side(q1,2), side(q1,3)
	e2,n2,p2 = side(q2,1), side(q2,2), side(q2,3)
	oc1,on1,op1 = opposite(t,e1), opposite(t,n1), opposite(t,p1)
	oc2,on2,op2 = opposite(t,e2), opposite(t,n2), opposite(t,p2)
	sc1,sn1,sp1 = apex(t,e1), apex(t,n1), apex(t,p1)
	sc2,sn2,sp2 = apex(t,e2), apex(t,n2), apex(t,p2)
	dc = typeof(e1)(3*(Int(q2)-Int(q1)))
	for (x,y) in (e2=>oc1, n2=>on1, p2=>op1, e1=>oc2, n1=>on2, p1=>op2)
		q = node(y)
		(q == q1) && (y+= dc); (q == q2) && (y-= dc)
		opposites!(t, x=>y)
# 	opposites!(t, e2=>oc1, n2=>on1, p2=>op1, e1=>oc2, n1=>on2, p1=>op2)
	end
	cellcorner!(t, sc1=>e2, sn1=>n2, sp1=>p2, sc2=>e1, sn2=>n1, sp2=>p1)
end#»»

# Delaunay triangulation ««1
# Cell location functions ««2
"""
    findcell(cornertable, points, point)

In a Voronoi diagram,
finds the index of the cell closest to the given point.
"""
function findcell(t::AbstractTriangulation{J}, points, point) where{J}
	# Returns the cell index closest to this point
	c0 = apex(t, Corner(1)); p0 = points[Int(c0)]
	while true
		c1 = c0
		for s in neighbours(t, c0)
			p = points[Int(s)]
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
	q = Node(J(rand(1:nnodes(t))))
	point = points[i]
	c = 0
	while true
		c+= 1; @assert c ≤ 1e3
		v = [points[Int(e)] for e in triangle(t,q)]
		isleft(v[1], point, v[2]) && (q = adjnode(t, q, 3); continue)
		isleft(v[2], point, v[3]) && (q = adjnode(t, q, 1); continue)
		isleft(v[3], point, v[1]) && (q = adjnode(t, q, 2); continue)
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
		isone(Int(q)) && continue # this is the reversed triangle
		q ∈ tree && continue
		meetscircle(t, q, points, i) || continue
		push!(tree, q)
		push!(stack, adjnodes(t, q)...)
	end
	return tree # the tree is usually small, no need to sort it
end#»»
function tree_boundary(t::AbstractTriangulation{J}, tree, doubleedges=()) where{J}#««
	# returns the list of half-edges pointing *into* the tree,
	# cyclically ordered around the tree (random start).
	boundary = sizehint!(Corner{J}[], length(tree)+2)
	for e in corners(first(tree))
		stack = [e]
		while !isempty(stack)
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
	boundary = tree_boundary(t, tree, doubleedges)
	push!(tree, newnodes!(t, 2)...)
	n = side(last(tree), 2)
	# rewrite boundary edges for double-edges, by substituting the
	# appropriate edge of the renamed triangle
	c2list = [ right(t, e) for e in boundary ]
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
	for (q, o, c2) in zip(tree, boundary, c2list)
		e = side(q, 1)
		p = side(q, 3)
		opposites!(t, e=>o, p=>n)
		n = side(q, 2)
		cellcorner!(t, s=>e, c1=>n, c2=>p)
		c1 = c2
	end
	check(t)
end#»»
# Triangulation
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
	npoints = length(points)
# 	npoints < 3 && return CornerTable{J}(undef, 0)
	t = CornerTable{J}(J[4,6,5,1,3,2],J(npoints).+J[1,3,2,1,2,3],
		Vector{J}(undef, npoints+3))
	anycorner!(t, Cell(npoints+1), Corner(J(1)))
	anycorner!(t, Cell(npoints+2), Corner(J(2)))
	anycorner!(t, Cell(npoints+3), Corner(J(3)))

	m = maximum(abs(x) for p in points for x in p)
	push!(points, [0,-3m], [3m,2m], [-3m,2m])
  #»»
	# incrementally add all points ««
	for i in 1:npoints # Random.randperm(npoints)
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
	fakecells = npoints+1:npoints+3
	for i in nnodes(t)-1:-1:1; q = Node{J}(i)
		w = Int.(triangle(t,q))
		any(>(npoints), Int.(triangle(t, q))) || continue
		swapnodes!(t, q, Node(k))
		k-= 1
	end
	# »»
		resize!(points, npoints)
	end
	return t
end#»»
# Voronoi diagram ««1
# Data structure ««2
struct VoronoiDiagram{J,P,T} <: AbstractTriangulation{J}
	# J: integer index type
	# P: geometric point type
	# T: real distance type
	triangulation::CornerTable{J}
	points::Vector{P}
	segments::Vector{NTuple{2,J}}
	geomnode::Vector{P}
	noderadius::Vector{T}
end

@inline VoronoiDiagram(t::CornerTable{J}, p::AbstractVector{P},
	s::AbstractVector, pos::AbstractVector, r::AbstractVector{T}) where{J,P,T} =
	VoronoiDiagram{J,P,T}(t, p, NTuple{2,J}.(s), P.(pos), r)
@inline VoronoiDiagram(t::CornerTable{J}, p::AbstractVector{P},
	s::AbstractVector) where{J,P} =
	VoronoiDiagram{J,P,eltype(P)}(t, p, NTuple{2,J}.(s),
	Vector{P}(undef, nnodes(t)), Vector{eltype(P)}(undef, nnodes(t)))
@inline triangulation(v::VoronoiDiagram) = v.triangulation

@inline npoints(v::VoronoiDiagram) = length(v.points)
@inline point(v::VoronoiDiagram, c::Cell) = v.points[Int(c)]
@inline issegment(v::VoronoiDiagram, c::Cell) = Int(c) > npoints(v)
@inline cellsegment(v::VoronoiDiagram, c::Cell) =
	Cell.(v.segments[Int(c)-npoints(v)])

function nnodes!(v::VoronoiDiagram, n)
	nnodes!(triangulation(v), n)
	resize!(v.geomnode, n)
	resize!(v.noderadius, n)
end
function swapnodes!(v::VoronoiDiagram, q1::Node, q2::Node)
	swapnodes!(triangulation(v), q1, q2)
	v.geomnode[SA[Int(q1),Int(q2)]] = v.geomnode[SA[Int(q2),Int(q1)]]
	v.noderadius[SA[Int(q1),Int(q2)]] = v.noderadius[SA[Int(q2),Int(q1)]]
end

# Geometric node computation ««2
@inline function geometricnode!(v::VoronoiDiagram, q::Node, g=equidistant(v,q))
	v.geomnode[Int(q)] = g
	s = cell(triangulation(v), q, 1)
	v.noderadius[Int(q)] = if issegment(v, s)
		(c1, c2) = cellsegment(v, s)
		segdistance²(point(v, c1), point(v, c2), g)
	else
		distance²(point(v, s), g)
	end
end
@inline geometricnode(v::VoronoiDiagram, q::Node) = v.geomnode[Int(q)]
@inline noderadius(v::VoronoiDiagram, q::Node) = v.noderadius[Int(q)]


function equidistant(v::VoronoiDiagram, q::Node)
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
end

function equidistant_sss(v::VoronoiDiagram, c1, c2, c3)
	(i1,j1),(i2,j2),(i3,j3)=cellsegment(v,c1),cellsegment(v,c2),cellsegment(v,c3)
	a = lineinter(v,i2,j2,i3,j3)
	b = lineinter(v,i3,j3,i1,j1)
	c = lineinter(v,i1,j1,i2,j2)
	# for parallel lines, there is no real distinction incenter/excenter...
	@assert a ≠ nothing
	@assert b ≠ nothing
	@assert c ≠ nothing
# 	a == nothing && return incenter_parallel(v, b, c, i2, j2, i3, j3)
# 	b == nothing && return incenter_parallel(v, c, a, i3, j3, i1, j1)
# 	c == nothing && return incenter_parallel(v, a, b, i1, j1, i2, j2)
	return incenter(a,b,c)
end

function lineinter(v,i1,j1,i2,j2)
	i1 ∈ (i2,j2) && return point(v,i1)
	j1 ∈ (i2,j2) && return point(v,j1)
	return lineinter(point(v,i1), point(v,j1), point(v,i2), point(v,j2))
end
@inline equidistant_ppp(v::VoronoiDiagram, c1, c2, c3) =
	# all three cells are single points
	circumcenter(point(v,c1), point(v,c2), point(v,c3))

function equidistant_pps(v::VoronoiDiagram, c1, c2, c3)
	(i,j) = cellsegment(v, c3)
	# treat cases where one of the points is one end of the segment
	c1 == i && return equidistant_pxs(v, i, j, c2)
	c2 == i && return equidistant_pxs(v, i, j, c1)
	c1 == j && return equidistant_pxs(v, j, i, c2)
	c2 == j && return equidistant_pxs(v, j, i, c1)
	return equidistant_pps(point(v, c1), point(v, c2), point(v, i), point(v, j))
end

function equidistant_pxs(v::VoronoiDiagram, c1, c2, c3)
	# node equidistant from: segment (c1,c2), points c1 & c3
	a, b, c = point(v,c1), point(v,c2), point(v,c3)
	ab, ac = b-a, c-a
	t = norm²(ac)/(2det2(ab,ac))
	return a+t*quarterturn(ab)
end

function equidistant_pss(v::VoronoiDiagram, s, c2,c3)
	(i1,j1),(i2,j2) = cellsegment(v,c2), cellsegment(v, c3)
	# equidistant from s, (i1,j1), (i2,j2)
	s == i1 && return equidistant_pxs(v, i1, j1, i2, j2, 1)
	s == j1 && return equidistant_pxs(v, j1, i1, i2, j2, 1)
	s == i2 && return equidistant_pxs(v, i2, j2, i1, j1, -1)
	s == j2 && return equidistant_pxs(v, j2, i2, i1, j1, -1)
	return equidistant_pss(point(v, s), point(v,i1), point(v,j1),
		point(v, i2), point(v, j2))
end

function equidistant_pxs(v::VoronoiDiagram, i1, j1, i2, j2, ε)
	i1 ∈ (i2,j2) && return point(v, i1)
	return equidistant_pxs(point(v,i1), point(v,j1), point(v,i2), point(v,j2), ε)
end

# Cell location functions ««2
function findrootnode(v::VoronoiDiagram, i,j)
	# here the table is assumed built (for points)
	# so we can search in the nodes around the cell for point a
	# which one is closest to segment ab
# 	t = triangulation(v)
	a, b = point(v,i), point(v,j)
	emin, dmin = nothing, nothing
	global V=v
	println("in star($i): $(collect(star(v,i)))")
	for e in star(v,i)
		influences(v,i,j,node(e)) || continue
		println("   ($i,$j) influences node $(node(e))")
		d = segdistance²(a,b,geometricnode(v, node(e)))
		(emin == nothing || d < dmin) && ((emin, dmin) = (e, d))
	end
	@assert emin ≠ nothing
	return node(emin)
end

function meetscircle(v::VoronoiDiagram, q::Node, i, j)
	g, r = geometricnode(v, q), noderadius(v, q)
	a, b = point(v, i), point(v, j)
	ab2 = norm²(b-a)
	d = dot(g-a,b-a)
	return segdistance²(a, b, g) < r
end
function influences(v::VoronoiDiagram, i, j, q)
	a,b,g = point(v, i), point(v,j), geometricnode(v,q)
	println("influences: ($i,$j,$q) points ($a,$b,$g)")
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
function edgeexit(v::VoronoiDiagram, e, i, j)
	# running between the geometric nodes g1 and g2,
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
	elseif issegment(v, c2) # _ps case
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
	@assert a >0
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
	# z = g1 + tv just as in the pp case
	# d²(z,(ab)) = <abz>²/ab² = (<abg>+t<abv>)²/ab²
	# d²(z,(pq)) = <pqz>²/pq² = (<pqg>+t<pqv>)²/pq²
	v = g2-g1
	ag, ab = g1-a, b-a
	pg, pq = g1-p, q-p
	abg, abv, ab2 = det2(ab, ag), det2(ab, v-a), norm²(ab)
	pqg, pqv, pq2 = det2(pq, pg), det2(pq, v-p), norm²(pq)
	dab_quadratic = SA[abv^2, 2*abv*abg, abg^2]/ab^2
	dpq_quadratic = SA[pqv^2, 2*pqv*pqg, pqg^2]/pq^2
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
function badnodes(v::VoronoiDiagram{J}, i, j) where{J} # segment case
	println("\e[34;1;7m badnodes($i, $j)\e[m")
	rootnode = findrootnode(v,i,j)
	println("\e[34;1m root node = $rootnode\e[m")
	stack = [(e, Cell{J}(0)) for e in corners(rootnode)]
	tree = [rootnode]
	loops = Cell{J}[]
	n = 0
	while !isempty(stack)
		n+=1
		@assert n ≤ 50
		(e, s) = pop!(stack)
		o = opposite(v, e)
		q = node(o)
		isone(Int(q)) && continue # this is the reversed triangle
		if q ∈ tree # one loop
			Int(s) ≠ 0 && push!(loops, s)
			continue
		end
		(influences(v,i,j,q) && meetscircle(v, q, i, j)) || continue
		push!(tree, q)
		push!(stack, (prev(o), left(v,e)))
		push!(stack, (next(o), right(v, e)))
	end
# 	@assert Int(i) ≠ 5
	# break loops in tree by identifying double edges
# 	@assert isempty(loops)
	doubleedges = Corner{J}[]
	for s in loops
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
			elist = collect(ring(v, s))
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
function star!(v::VoronoiDiagram, s, tree, doubleedges)#««
	star!(triangulation(v), s, tree, doubleedges)
	global V = v
	resize!(v.geomnode, length(v.geomnode)+2)
	resize!(v.noderadius, length(v.geomnode))
	for q in tree
		geometricnode!(v, q)
	end
end#»»
# Triangulation ««2
@inline voronoi(points, segments=[]) = VoronoiDiagram{Int32}(points, segments)
function VoronoiDiagram{J}(points, segments) where{J}#««
	Random.seed!(0)
	np = length(points)
	ns= length(segments)
	ntotal = np+ ns
	points = SVector{2,Float64}.(points)
	t = CornerTable{J}(points; trim=false)
	ncells!(t, ncells(t) + ns)
	v = VoronoiDiagram(t, points, segments)
	for q in allnodes(t); geometricnode!(v, q); end

	# incrementally add all segments ««
	ns = length(segments)
	for i in 1:ns # Random.randperm(ns)
		s = Cell(np+3+i); (a,b) = cellsegment(v, s)
		println("\e[31;1;7minserting segment $s = ($a,$b)\e[m")
		tree, doubleedges = badnodes(v, a, b)
		star!(v, s, tree, doubleedges)
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
# 			println("  e = $e ($(apex(v,e)), $(right(v,e)), $(left(v,e)))")
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
"""
    splitsegments!(voronoidiagram)

Splits segment cells in two parts depending on the orientation of the segment:
(right of segment, left of segment).
"""
function splitsegments!(v::VoronoiDiagram{J}) where{J}#««
	np = length(v.points)
	ns = length(v.segments)
# 	display(v)
	ncells!(v, np+2ns)
	# rename cells to alternate segments and their opposites
	# since we are swapping with undefined cells,
	# there is no need to worry about side-effects
	for i in ns:-1:1
		movecell!(v, Cell(np+i), Cell(np+2i-1))
	end
	origsegments = copy(v.segments)
	sizehint!(empty!(v.segments), 2ns)
	for (i, (a,b)) in pairs(origsegments)
		push!(v.segments, (a,b), (b,a))
	end
	println("segments(v) = $(v.segments)")
	# now split each cell in two
	for i in 1:ns
		c1, c2 = Cell(np+2i-1), Cell(np+2i)
		(a,b) = cellsegment(v, c1)
		(pa, pb) = point(v, a), point(v, b); ab = pb - pa
		println("\e[7m splitting cell $c1 = ($a,$b)\e[m")
# 		println("two points $pa->$pb; vector $ab")
		for e in star(v,c1)
			g = geometricnode(v, node(e))
# 			println("  edge $e ↔ node $(node(e)) ↔ cells $(triangle(v,node(e))):")
# 			println("  geom. $g: ", sign(det2(g-pa, ab)))
		end
		elist = star(v,c1)
		(e, status) = iterate(elist)
		sidechange=MVector(e,e)
		g = geometricnode(v, node(e))
		sd = (det2(g-pa, ab) ≥ 0)
		while true
			next = iterate(elist, status)
			next == nothing && break
			sd1 = sd
			(e, status) = next
			g = geometricnode(v, node(e))
			sd = (det2(g-pa, ab) ≥ 0)
# 			sd || apex!(v, e, c2)
			sd1 ≠ sd && (sidechange[1+sd] = e)
		end
		println(" \e[34msidechange $sidechange\e[m = nodes $(node.(sidechange))")
		# split the cell by inserting two new nodes
		e1, e2 = sidechange
		pe1, pe2 = prev(e1), prev(e2)
		ope1, ope2 = opposite(v, prev(e1)), opposite(v, prev(e2))
		q1, q2 = newnodes!(v, 2)
		apex!(v,
			side(q1,1)=>right(v,e1), side(q1,2)=>c2, side(q1,3)=>c1,
			side(q2,1)=>right(v,e2), side(q2,2)=>c1, side(q2,3)=>c2)
		opposites!(v, side(q1,1) => side(q2,1),
			pe1=>side(q1,3), ope1=>side(q1,2),
			pe2=>side(q2,3), ope2=>side(q2,2))
		geometricnode!(v, q1, pb)
		geometricnode!(v, q2, pa)
# 		shownode(stdout, v, q1)
# 		shownode(stdout, v, q2)
		for e in star(v, e1) # iterate on left side to change apex
			apex!(v, e, c2)
		end
		anycorner!(v, c1, side(q2,2)); anycorner!(v, c2, side(q1,2))
# 		showcell(stdout, v, c1)
# 		showcell(stdout, v, c2)
	end
	return v
end#»»
# function flip_rec(t::CornerTable, points, c::Corner)#««
# 	# c is the corner (p,new cell, a, b)
# 	# o is the opposite corner
# 	# if o is in circle of (a,b,p) then flip the corner c
# 	sc, sa, sb, so = apex(t, c), apex(t, next(c)), apex(t, prev(c)),
# 		apex(t, opposite(t, c))
# 	isincircle(points[Int.(SA[sc,sa,sb,so])]...) || return
# 	(c1, c2) = flip!(t, c)
# 	flip_rec(t, points, c1)
# 	flip_rec(t, points, c2)
# end#»»

function check(t::AbstractTriangulation)
	for i in 1:3*nnodes(t)
		e = Corner(i)
		@assert Int(opposite(t,opposite(t,e))) == Int(e) "opposite($e)=$(opposite(t,e)); opposite($(opposite(t,e)))=$(opposite(t,opposite(t,e)))"
	end
# 	for i in 1:ncells(t)
# 		s = Cell(i)
# 		@assert Int(apex(t, anycorner(t,s))) == Int(s) "apex(anycorner($s))"
# 	end
end
function shownode(io::IO, v::VoronoiDiagram, q::Node)
	println(io, "\e[33m", q, triangle(v,q), "\e[m: ", geometricnode(v,q), " d=",
		noderadius(v,q))
	for i in (1,2,3)
		e = side(q, i); o = opposite(v, e); oo = opposite(v,o)
		oo ≠ e && println(io, "  \e[31;7m opposite($o) = $oo, should be $e\e[m")
	end
end
function showcell(io::IO, v::AbstractTriangulation, c::Cell)
	println(io, "\e[34m", c, "\e[m: ", [node(e) for e in star(v,c)])
	for e in star(v,c)
		c1 = apex(v, e)
		c1 ≠ c && println(io, "  \e[31;7m apex($e) = $c1, should be $c\e[m")
	end
end
function Base.display(io::IO, v::VoronoiDiagram)
	println(io, "\e[1m begin Voronoi diagram with $(nnodes(v)) nodes and $(ncells(v)) cells:\e[m")
	for q in allnodes(v); shownode(io, v, q); end
	for c in allcells(v); showcell(io, v, c); end
	println(io, "\e[1m end Voronoi diagram\e[m")
end
@inline Base.display(v::VoronoiDiagram) = display(stdout, v)
# end »»1
# Offset ««1
function offset(points, segments, radius)
	v = voronoi(points, segments)
	return offset(v, radius)
end
function offset(v::VoronoiDiagram, radius)
	# do a walk on the set of segments by connectivity
	#   (we only start on segments but we can walk to points if an edge
	#   points us to them)
	# for each segment, determine the offset of the segment by the radius
	# and in particular all intersection points with the edges
	# then glue along the edges
	# as long as there are segments remaining...
end



#»»1
end

V=Voronoi
# using StaticArrays
# TODO: tests
# V = Voronoi
# t=V.triangulate([[-10,0],[10,0],[0,10.]])
# v=V.voronoi([(0,0),(10,0),(11,3),(6,2),(5,9),(1,10)],[(6,1),(1,2),(2,3),(3,4),(4,5),(5,6)])
v=V.voronoi([[0.,0],[10,0],[0,10],[10,10],[5,9],[5,1]],[(3,4),(5,6)])
# V.splitsegments!(v)
# m = V.edgeexit_pp([5,9],[5,1],[1.6,5],[8.4,5],[0,10],[10,10])
# V.equidistant_pss([0,0],[5,9],[5,1],[0,10],[10,10]) # [5,5]

# p=NTuple{2,Float64}[(0,0),(10,0),(11,3),(6,2),(5,9),(1,10)]
# V.offset(p,[(6,1),(1,2),(2,3),(3,4),(4,5),(5,6)],0.5)

# collect(V.allnodes(v))
