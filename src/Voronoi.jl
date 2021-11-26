# With some ideas taken from VRONI as described by [Huber 2008]:
# https://www.sthu.org/research/publications/files/mscthesis.pdf
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
	return a + t*quarterturn(b-a)
end

"returns the point equidistant from point a and segments (pq), (rs)"
function equidistant_pss(a, p, q, r, s)
	pqa = det2(p,q,a)
	pqa < 0 && ((p,q) = (q,p); pqa = -pqa) # ensure orientation
	pq = q-p
	pq2 = norm²(pq)
	pqr, pqs = det2(p,q,r), det2(p,q,s)
	# `a` is not allowed to be on any of the two lines pq, rs
	if pqr == pqs # special case: parallel lines
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
	rsa = det2(r,s,a)
	rsa > 0 && ((r,s) = (s,r); rsa = -rsa)
	# parametrization of the inner angle bisector: z = c + t u
	c = lineinter(p, q, r, s)
	u = (q-p)/√(pq2) + (s-r)/√(distance²(r,s))
	# for z = c+tu: d²(z,pq) = t^2 pqu²/pq², while
	# d²(z,a) = ‖a-c-tu‖² = t² u²- 2t ca⋅u + ca²
	# the equation is thus t²(u²-pqu²/pq²) -2t ca⋅u + ca² = 0, or
	# t²(pq⋅u)²/pq² -2t ca⋅u + ca² = 0
	A = dot(q-p, u)^2 / pq2
	B = dot(a-c, u)
	C = distance²(c, a)
	Δ = B^2-A*C
	ε = sign(pqs - pqr)
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
# Data type ««2
abstract type AbstractTriangulation{J} end
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

# Accessors ««2
# Corner functions ««3

# LazyRow accepts only `Int` indexes:
@inline apex(t::CornerTable, e::Corner) = Cell(t.apex[Int(e)])
@inline apex!(t::CornerTable, e::Corner, s::Cell) = t.apex[Int(e)] = Int(s)
@inline apex!(t::CornerTable, l::Pair{<:Corner, <:Cell}...) =
	for (e,s) in l; apex!(t, e, s); end

@inline opposite(t::CornerTable, e::Corner) = Corner(t.opposite[Int(e)])
@inline opposite!(t::CornerTable, e::Corner, x::Corner) =
	t.opposite[Int(e)] = Int(x)
@inline opposites!(t::AbstractTriangulation, l::Pair{<:Corner,<:Corner}...) =
	for(e1, e2) in l; opposite!(t, e1, e2); opposite!(t, e2, e1); end

# Node functions ««3
@inline nnodes(t::CornerTable) = length(t.opposite) ÷ 3
@inline nnodes!(t::CornerTable, n) = 
	(resize!(t.opposite, 3n); resize!(t.apex, 3n))

# Cell functions ««3
@inline anycorner(t::CornerTable, s::Cell) = Corner(t.anycorner[Int(s)])
@inline anycorner!(t::CornerTable, s::Cell, e::Corner) =
	t.anycorner[Int(s)] = Int(e)
@inline ncells(t::CornerTable) = length(t.anycorner)
@inline ncells!(t::CornerTable, n) = resize!(t.anycorner, n)

# AbstractTriangulation methods ««2
@inline right(t::AbstractTriangulation, e::Corner) = apex(t, next(e))
@inline left(t::AbstractTriangulation, e::Corner) = apex(t, prev(e))
@inline after(t::AbstractTriangulation, e::Corner) = next(opposite(t, next(e)))
@inline before(t::AbstractTriangulation, e::Corner) = prev(opposite(t, prev(e)))

@inline cell(t::AbstractTriangulation, q::Node, i) = apex(t, side(q,i))
@inline cells(t::AbstractTriangulation, q::Node) =
	(cell(t,q,1), cell(t,q,2), cell(t,q,3))
@inline allnodes(t::AbstractTriangulation{J}) where{J} =
	(Node(J(i)) for i in 1:nnodes(t))
@inline alltriangles(t::AbstractTriangulation) =
	(cells(t,q) for q in allnodes(t))
@inline adjnode(t::AbstractTriangulation, q::Node, i) =
	node(opposite(t, side(q, i)))
@inline adjnodes(t::AbstractTriangulation, q::Node) =
	(adjnode(t,q,1), adjnode(t,q,2), adjnode(t,q,3))

function newnodes!(t::AbstractTriangulation{J}, k) where{J}
	n = nnodes(t)
	nnodes!(t, n+k)
	return Node{J}.(n+1:n+k)
end

@inline anycorner!(t::AbstractTriangulation, l::Pair{<:Cell, <:Corner}...) =
	for (s,e) in l; anycorner!(t, s, e); end
@inline cellcorner!(t::AbstractTriangulation, l::Pair{<:Cell, <:Corner}...) =
	for(s,e) in l; anycorner!(t, s, e); apex!(t, e, s); end
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
"""
    star(cornertable, cell)

Iterates all corners from the same cell.
"""
@inline star(t::AbstractTriangulation{J}, s::Cell) where{J} =
	CTIterator{star,Corner{J}}(t, anycorner(t, s))
@inline next(::CTIterator{star}, t, e) = after(t, e)
@inline stop(it::CTIterator{star}, t, e) = (e == it.start)

"""
    ring(cornertable, corner)

Iterates all corners left-adjacent to the indicated cell.
"""
@inline ring(t::AbstractTriangulation{J}, e::Corner) where{J} =
	CTIterator{ring,Corner{J}}(t, next(e))
@inline ring(t::AbstractTriangulation, s::Cell) = ring(t, anycorner(t, s))
@inline next(::CTIterator{ring}, t, e) = prev(opposite(t, e))
@inline stop(it::CTIterator{ring}, t, e) = (e == it.start)
"""
   neighbours(cornertable, cell)

Iterates all cells adjacent to the indicated cell.
"""
@inline neighbours(t::AbstractTriangulation{J}, s::Cell) where{J} =
	(apex(t, x) for x in ring(t, anycorner(t, s)))

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
function insert!(t::AbstractTriangulation, q::Node, s::Cell)
	e0, n0, p0 = corners(q)
	s0, s1, s2 = apex(t, e0), apex(t, n0), apex(t, p0)
	# e0 and its opposite corner are unchanged
	on, op = opposite(t, n0), opposite(t, p0)
	q1, q2 = newnodes!(t, 2)
	e1, n1, p1 = corners(q1)
	e2, n2, p2 = corners(q2)
	opposites!(t, e1=>on, e2=>op, n0=>p1, n1=>p2, n2=>p0)
	apex!(t, e0=>s, e1=>s, e2=>s, n1=>s2, p1=>s0, n2=>s0, p2=>s1)
	anycorner!(t, s=>e0, s1=>n0, s2=>p0, s0=>n1)
	return (e0, e1, e2)
end

"""
    swapcells!(cornertable, cell, cell)

Swap cell indices for those two cells.
"""
function swapcells!(t::AbstractTriangulation, s1::Cell, s2::Cell)#««
	s1 == s2 && return
	for e in star(t, s1); apex!(t, e, s2); end
	for e in star(t, s2); apex!(t, e, s1); end
	e1 = anycorner(t, s1)
	anycorner!(t, s1, anycorner(t, s2))
	anycorner!(t, s2, e1)
end#»»
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
function findcell(t::CornerTable{J}, points, point) where{J}
	# Returns the cell index closest to this point
	s0 = apex(t, Corner(1)); p0 = points[Int(s0)]
	while true
		s1 = s0
		for s in neighbours(t, s0)
			p = points[Int(s)]
			iscloser(point, p0, p) && continue
			s0, p0 = s, p
			break
		end
		s0 == s1 && return s0
	end
end

"""
    findnode(cornertable, points, point)

In a triangulation, returns the index of the triangle containing the point.
"""
function findnode(t::CornerTable{J}, points, i) where{J}
	q = Node(J(rand(1:nnodes(t))))
	point = points[i]
	c = 0
	while true
		c+= 1; @assert c ≤ 1e3
		v = [points[Int(e)] for e in cells(t,q)]
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
function meetscircle(t::CornerTable, q::Node, points, i)#««
	s1, s2, s3 = Int.(cells(t, q))
	return isincircle(points[SA[s1,s2,s3, i]]...)
end#»»
function badnodes(t::CornerTable, points, i)#««
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
function tree_boundary(t::CornerTable{J}, tree) where{J}#««
	# returns the list of half-edges pointing *into* the tree,
	# cyclically ordered around the tree (random start).
	boundary = sizehint!(Corner{J}[], length(tree)+2)
	for e in corners(first(tree))
		stack = [e]
		while !isempty(stack)
			e = pop!(stack)
			o = opposite(t, e)
			if node(o) ∈ tree
				push!(stack, prev(o), next(o))
			else
				push!(boundary, o)
			end
		end
	end
	return boundary
end#»»
function star!(t::CornerTable, s, tree)#««
	# replaces all nodes in `tree` by a star shape around `cell`,
	# adding two new nodes in the process.
	boundary = tree_boundary(t, tree)
	push!(tree, newnodes!(t, 2)...)
	n = side(last(tree), 2)
	for (q, o) in zip(tree, boundary)
		e = side(q, 1)
		p = side(q, 3)
		opposites!(t, e=>o, p=>n)
		n = side(q, 2)
		cellcorner!(t, s=>e, right(t,o)=>p, left(t,o)=>n)
	end
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
function triangulate(points; trim=true)#««
	J = Int32
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
	for i in Random.randperm(npoints)
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
		w = Int.(cells(t,q))
		any(>(npoints), Int.(cells(t, q))) || continue
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
@inline triangulation(v::VoronoiDiagram) = v.triangulation
@inline apex(v::VoronoiDiagram, e::Corner) = apex(triangulation(v), e)
@inline opposite(v::VoronoiDiagram, e::Corner) = opposite(triangulation(v), e)
@inline anycorner(v::VoronoiDiagram, q::Cell) = anycorner(triangulation(v), q)
@inline nnodes(v::VoronoiDiagram) = nnodes(triangulation(v))
@inline npoints(v::VoronoiDiagram) = length(v.points)
@inline point(v::VoronoiDiagram, s::Cell) = v.points[Int(s)]
@inline issegment(v::VoronoiDiagram, s::Cell) = Int(s) > npoints(v)
@inline cellsegment(v::VoronoiDiagram, s::Cell) =
	Cell.(v.segments[Int(s)-npoints(v)])
function swapnodes!(v::VoronoiDiagram, q1::Node, q2::Node)
	swapnodes!(triangulation(v), q1, q2)
	v.geomnode[SA[Int(q1),Int(q2)]] = v.geomnode[SA[Int(q2),Int(q1)]]
	v.noderadius[SA[Int(q1),Int(q2)]] = v.noderadius[SA[Int(q2),Int(q1)]]
end

# Geometric node computation ««2
@inline function geometricnode!(v::VoronoiDiagram, q::Node)
	g = equidistant(v, q)
	println("  geometricnode($(cells(v,q))): found $g")
	v.geomnode[Int(q)] = g
	s = cell(triangulation(v), q, 1)
	v.noderadius[Int(q)] = if issegment(v, s)
		(s1, s2) = cellsegment(v, s)
		segdistance²(point(v, s1), point(v, s2), g)
	else
		distance²(point(v, s), g)
	end
end
@inline geometricnode(v::VoronoiDiagram, q::Node) = v.geomnode[Int(q)]
@inline noderadius(v::VoronoiDiagram, q::Node) = v.noderadius[Int(q)]

@inline allnodes(v::VoronoiDiagram) =
	((cells(triangulation(v), Node(q)), v.geomnode[q], v.noderadius[q])
	for q in 1:nnodes(v))

function equidistant(v::VoronoiDiagram, q::Node)
	s = cells(v, q)
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

function equidistant_sss(v::VoronoiDiagram, s1, s2, s3)
	(i1,j1),(i2,j2),(i3,j3)=cellsegment(v,s1),cellsegment(v,s2),cellsegment(v,s3)
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
@inline equidistant_ppp(v::VoronoiDiagram, s1, s2, s3) =
	# all three cells are single points
	circumcenter(point(v,s1), point(v,s2), point(v,s3))

function equidistant_pps(v::VoronoiDiagram, s1, s2, s3)
	(i,j) = cellsegment(v, s3)
	# treat cases where one of the points is one end of the segment
	s1 == i && return equidistant_pxs(v, i, j, s2)
	s2 == i && return equidistant_pxs(v, i, j, s1)
	s1 == j && return equidistant_pxs(v, j, i, s2)
	s2 == j && return equidistant_pxs(v, j, i, s1)
	return equidistant_pps(point(v, s1), point(v, s2), point(v, i), point(v, q))
end

function equidistant_pxs(v::VoronoiDiagram, s1, s2, s3)
	# node equidistant from: segment (s1,s2), points s1 & s3
	a, b, c = point(v,s1), point(v,s2), point(v,s3)
	ab, ac = b-a, c-a
	t = norm²(ac)/(2det2(ab,ac))
	return a+t*quarterturn(ab)
end

function equidistant_pss(v::VoronoiDiagram, s, s2,s3)
	(i1,j1),(i2,j2) = cellsegment(v,s2), cellsegment(v, s3)
	# equidistant from s, (i1,j1), (i2,j2)
# 	println("  node2: equidistant from $s, ($i1,$j1), ($i2,$j2)")
	s == i1 && return equidistant_pxs(v, i1, j1, i2, j2, 1)
	s == j1 && return equidistant_pxs(v, j1, i1, i2, j2, 1)
	s == i2 && return equidistant_pxs(v, i2, j2, i1, j1, -1)
	s == j2 && return equidistant_pxs(v, j2, i2, i1, j1, -1)
	return equidistant_pss(point(v, s), point(v,i1), point(v,j1),
		point(v, i2), point(v, j2))
end

function equidistant_pxs(v::VoronoiDiagram, i1, j1, i2, j2, ε)
# 	println("  node2a: equidistant from $i1, ($i1,$j1), ($i2, $j2)")
	i1 ∈ (i2,j2) && return point(v, i1)
	return equidistant_pxs(point(v,i1), point(v,j1), point(v,i2), point(v,j2), ε)
end

# Cell location functions ««2
function findnode(v::VoronoiDiagram, a, b)
	# here the table is assumed built (for points)
	# so we can search in the nodes around the cell for point a
	# which one is closest to segment ab
	t = triangulation(v)
	pa, pb = point(v,a), point(v,b)
	e = argmin(segdistance²(pa, pb, geometricnode(v, node(e)))
		for e in star(t, a))
	return node(e)
end

function meetscircle(v::VoronoiDiagram, q, i, j)
	g, r = geometricnode(v, q), noderadius(v, q)
	a, b = point(v, Cell(i)), point(v, Cell(j))
	ab2 = norm²(b-a)
	d = dot(g-a,b-a)
	return segdistance²(a, b, g) < r
end

# Computing Voronoi diagram ««1
# Tree traversal ««2
function badnodes(v::VoronoiDiagram{J}, i, j) where{J} # segment case
	t = triangulation(v)
	stack = [findnode(v, Cell(i),Cell(j))]
	tree = empty(stack)
	loops = Cell{J}[]
	while !isempty(stack)
		q = pop!(stack)
		isone(Int(q)) && continue # this is the reversed triangle
		@assert q ∉ tree # loops should have already been prevented
		meetscircle(v, q, i, j) || continue
		push!(tree, q)
		# continue walking the tree, avoiding loops««
		a = adjnodes(t,q)
# 		println("    this segments crosses the node circle, continuing... (adj=$a)")
		if a[1] in tree
			if a[2] in tree
				@assert a[3] ∉ tree
				push!(stack, a[1]); push!(loops, cell(t, q, 1))
			else
				if a[3] in tree push!(stack, a[2]); push!(loops, cell(t, q, 2))
				else push!(stack, a[2], a[3])
				end
			end
		else
			if a[2] in tree
				if a[3] in tree push!(stack, a[1]); push!(loops, cell(t, q, 1))
				else push!(stack, a[3], a[1])
				end
			else
				if a[3] in tree push!(stack, a[1], a[2])
				else push!(stack, a[1], a[2], a[3])
				end
			end
		end#»»
	end
	# break loops in tree by identifying double edges
	@assert isempty(loops)
	for q in loops
		# find an edge in which we are sure that *some point* is closer to
		# its defining sites (i.e. left(e), right(e)) than to the site we are
		# inserting
		for e in star(t, q)
		end
	end
	return tree # the tree is usually small, no need to sort it
end
# end#»»
function star!(v::VoronoiDiagram, s, tree)#««
	star!(triangulation(v), s, tree)
	resize!(v.geomnode, length(v.geomnode)+2)
	resize!(v.noderadius, length(v.geomnode))
	for q in tree
		geometricnode!(v, q)
	end
end#»»
# Triangulation ««2
function voronoi(points, segments = [])#««
	J = Int32
	Random.seed!(0)
	npoints = length(points)
	nsegments = length(segments)
	ntotal = npoints + nsegments
	t = triangulate(points; trim=false)
	ncells!(t, ncells(t) + nsegments)
	println(collect(alltriangles(t)))
	# compute all circumcenters ««
	# (not doing this during pointwise triangulation saves a bit of time)
	geomnodes = [ circumcenter(points[Int(cell(t,q,1))],
		points[Int(cell(t,q,2))], points[Int(cell(t,q,3))]) for q in allnodes(t) ]
	noderadius = [ distance²(g, points[Int(apex(t, Corner(3i)))])
		for (i,g) in pairs(geomnodes) ]
	v = VoronoiDiagram(t, points, segments, geomnodes, noderadius)
	#»»
	# incrementally add all segments ««
	nsegments = length(segments)
	for i in Random.randperm(nsegments)
		(a,b) = segments[i]
		tree = badnodes(v, a, b)
		println("  tree for cell $i = ($a,$b) is $tree")
		star!(v, Cell(npoints+3+i), tree)
	end
	#»»
	# remove all superfluous nodes & cells ««
	# the nodes are sorted this way:
	# - inner nodes
	# - convex hull
	# - 1 backwards outer node
	k = nnodes(t)
	swapnodes!(v, Node(1), Node(k))
	k-= 1
	fakecells = npoints+1:npoints+3
	for i in nnodes(t)-1:-1:1; q = Node{J}(i)
		w = Int.(cells(t,q))
		w[1] ∈ fakecells || w[2] ∈ fakecells || w[3] ∈ fakecells || continue
# 		any(>(ntotal), Int.(cells(t, q))) || continue
		swapnodes!(v, q, Node(k))
		k-= 1
	end
# 	nnodes!(t, k)
	# »»
	resize!(points, npoints)
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

# end »»1

vert = [SA[0,0], SA[10,0], SA[0,10], SA[11,11]]

end

using StaticArrays
V = Voronoi
# t=V.triangulate([[-10,0],[10,0],[0,10.]])
# v=V.voronoi([[-10,0],[10,0],[0,10.]],[(1,3),(3,2),(2,1)])
# collect(V.allnodes(v))
