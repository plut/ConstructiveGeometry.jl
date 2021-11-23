# [guibas stolfi]
module Voronoi
using StaticArrays
using StructArrays
using LazyArrays
using FastClosures
using LinearAlgebra

import Base: Int

# Tools ««1
# Geometry ««2

"""
    iscloser(a,b,c)

Returns true iff d(a,b) ≤ d(a,c).
"""
@inline iscloser(a,b,c) = dot(2a-b-c, b-c) ≥ 0

@inline det2(u,v) = u[1]*v[2]-u[2]*v[1]
@inline det2(u,v,w) = det2(v-u, w-u)
"""
    isleft(a,b,c)

Returns true iff c is to the left of ab.
"""
@inline isleft(a,b,c) = det2(a,b,c) > 0
@inline isleft(u,v) = det2(u,v) > 0

@inline isincircle(a,b,c,x) = isincircle(a-x, b-x, c-x)
function isincircle(a,b,c) # is point zero in this circle (assumed oriented)
	@assert !isleft(a,c,b) "incircle: triangle ($a,$b,$c) has wrong orientation"
	m = SA[a[1] a[2] a[1]^2+a[2]^2; b[1] b[2] b[1]^2+b[2]^2; c[1] c[2] c[1]^2+c[2]^2]
	return det(m) > 0
end

# Named indices ««2
struct NamedIndex{S,T<:Signed} i::T; end
@inline Int(i::NamedIndex) = i.i
@inline Base.convert(T::Type{<:Integer}, i::NamedIndex) = T(Int(i))
@inline Base.show(io::IO, i::NamedIndex{S}) where{S} = print(io, S, Int(i))
@inline NamedIndex{S}(x) where{S} = NamedIndex{S,typeof(x)}(x)
@inline Base.isless(x::NamedIndex{S}, y::NamedIndex{S}) where{S} = Int(x)<Int(y)
macro NamedIndex(name,S)
	quote
	const $(esc(name)) = NamedIndex{$(QuoteNode(S))}
	Base.show(io::IO, ::Type{<:$(esc(name))}) = print(io, $(string(name)))
end end
@NamedIndex Corner c
@NamedIndex Site s
@NamedIndex Triangle t

# CornerTable ««1
# Data type ««2
"""
    CornerTable - generic table for a triangulation (or its dual)
    of a list of sites.

This table encodes two types of objects:
 - triangles (made of three consecutive corners in the table),
 - sites (duals of triangles, made of adjacent corners)

For triangulations, triangles=faces and sites=vertices;
the opposite of a corner is pictured as <corner | opposite >.
For Voronoi cells, triangles=vertices and sites=faces;
the opposite of a corner is pictured as corner>———<opposite.
"""
struct CornerTable{J,CS<:StructArray}
	# the exact CS type is computed at constructor time
	cornerdata::CS
	anycorner::Vector{J} # site -> corner map
	@inline (T::Type{<:CornerTable{J}})(cdata::CS, a::AbstractVector
		) where{J,CS<:StructArray} = new{J,CS}(cdata, a)
	@inline (T::Type{<:CornerTable{J}})(undef, nsites) where{J} =
		T(StructArray((opposite=J[], site=J[])), zeros(J, nsites))
	@inline (T::Type{<:CornerTable{J}})() where{J} = T(undef, 0)
end

# for broadcast:
@inline Base.Broadcast.broadcastable(x::CornerTable) = (x,)


# Accessors ««2
@inline next3(i::J) where{J<:Integer} = iszero(i%3) ? i-2 : i+1
@inline prev3(i::J) where{J<:Integer} = isone(i%3) ? i+2 : i-1
@inline next(c::Corner{J}) where{J} = Corner{J}(next3(Int(c)))
@inline prev(c::Corner{J}) where{J} = Corner{J}(prev3(Int(c)))
@inline triangle(c::Corner{J}) where{J} = Triangle{J}(fld1(Int(c), 3))
@inline side(q::Triangle{J}, i) where{J} = Corner{J}(3*Int(q)-3+i)
@inline corners(q::Triangle{J}) where{J} = Corner{J}.(3*Int(q) .- (2,1,0))
# LazyRow accepts only `Int` indexes:
@inline getc(t::CornerTable, c::Corner) = LazyRow(t.cornerdata, Int(Int(c)))
@inline opposite(t::CornerTable, c::Corner) = Corner(getc(t, c).opposite)
@inline opposite!(t::CornerTable, c::Corner, x::Corner) =
	getc(t, c).opposite = Int(x)
@inline opposites!(t::CornerTable, l::Pair{<:Corner,<:Corner}...) =
	for(c1, c2) in l; opposite!(t, c1, c2); opposite!(t, c2, c1); end
@inline site(t::CornerTable, c::Corner) = Site(getc(t, c).site)
@inline site!(t::CornerTable, c::Corner, s::Site) =
	getc(t, c).site = Int(s)
@inline site!(t::CornerTable, l::Pair{<:Corner, <:Site}...) =
	for (c,s) in l; site!(t, c, s); end
@inline sites(t::CornerTable, q::Triangle) =
	(site(t, side(q,1)), site(t, side(q,2)), site(t, side(q,3)))
# 	(site(t, side(q, i)) for i in (1,2,3))
@inline after(t::CornerTable, c::Corner) = next(opposite(t, next(c)))
@inline before(t::CornerTable, c::Corner) = prev(opposite(t, prev(c)))

@inline anycorner(t::CornerTable, s::Site) = Corner(t.anycorner[Int(s)])
@inline anycorner!(t::CornerTable, s::Site, c::Corner) =
	t.anycorner[Int(s)] = Int(c)
@inline anycorner!(t::CornerTable, l::Pair{<:Site, <:Corner}...) =
	for (s,c) in l; anycorner!(t, s, c); end
@inline sitecorner!(t::CornerTable, l::Pair{<:Site, <:Corner}...) =
	for(s,c) in l; anycorner!(t, s, c); site!(t, c, s); end
@inline ntriangles(t::CornerTable) = length(t.cornerdata) ÷ 3
@inline ntriangles!(t::CornerTable, n) = resize!(t.cornerdata, 3n)
@inline triangles(t::CornerTable) =
	(((Int(site(t, c)) for c in corners(q))...,)
		for q in Triangle.(1:ntriangles(t)))
@inline nsites(t::CornerTable) = length(t.anycorner)
@inline nsites!(t::CornerTable, n) = resize!(t.anycorner, n)
function new_triangles!(t::CornerTable{J}, k) where{J}
	n = ntriangles(t)
	ntriangles!(t, n+k)
	return Triangle{J}.(n+1:n+k)
end

# Iterators ««2
struct CTIterator{S,X,J,T}
	table::CornerTable{J}
	start::T
	@inline CTIterator{S,X}(t::CornerTable{J}, s::T) where{S,X,J,T}=
		new{S,X,J,T}(t,s)
end
@inline Base.eltype(::CTIterator{S,X}) where{S,X} = X
@inline Base.IteratorSize(::CTIterator) = Base.SizeUnknown()
@inline Base.iterate(it::CTIterator) = (it.start, it.start)
@inline function Base.iterate(it::CTIterator, s)
	s = next(it, it.table, s)
	stop(it, it.table, s) && return nothing
	return (s, s)
end
"""
    star(cornertable, site)

Iterates all corners from the same site.
"""
@inline star(t::CornerTable{J}, s::Site) where{J} =
	CTIterator{star,Corner{J}}(t, anycorner(t, s))
@inline next(::CTIterator{star}, t, c) = after(t, c)
@inline stop(it::CTIterator{star}, t, c) = (c == it.start)

"""
    ring(cornertable, corner)

Iterates all corners left-adjacent to the indicated site.
"""
@inline ring(t::CornerTable{J}, c::Corner) where{J} =
	CTIterator{ring,Corner{J}}(t, next(c))
@inline next(::CTIterator{ring}, t, c) = prev(opposite(t, c))
@inline stop(it::CTIterator{ring}, t, c) = (c == it.start)
"""
   neighbours(cornertable, site)

Iterates all sites adjacent to the indicated site.
"""
@inline neighbours(t::CornerTable{J}, s::Site) where{J} =
	(site(t, x) for x in ring(t, anycorner(t, s)))

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
	for (c, s) in pairs(slist)
		clist[s] = c
		n = next3(c); sn = slist[n]
		push!(edge[sn], c)
	end
	for (c, s) in pairs(slist)
		sn, sp = slist[next3(c)], slist[prev3(c)]
		for c1 in edge[sp]
			slist[prev3(c1)] == sn || continue
			olist[c] = c1; olist[c1] = c
		end
	end
	return CornerTable{J}(StructArray((opposite=olist, site=slist)), clist)
end

@inline CornerTable{J}(::Val{:tetrahedron}) where{J} =
	# initial Voronoi diagram for 4 sites (3 points + 1 at infinity):
	CornerTable{J}(
		StructArray((opposite=J[4,7,10,1,12,8,2,6,11,3,9,5],
			site=J[1,2,3, 4,3,2, 4,1,3, 4,2,1])),
			J[1,2,3,4])

# Elementary modifications: flip!, insert! ««2
"""
    flip!(cornertable, corner)

Flips the quadrilateral delimited by `corner` and its opposite.
Returns the two corners (in order) replacing `corner`.
"""
function flip!(t::CornerTable, c::Corner)
	n, p, o = next(c), prev(c), opposite(t, c)
	s, sn, so, sp = site(t, c), site(t, n), site(t, o), site(t, p)
	no, po = next(o), prev(o)
	# we leave (n ↔ on) and (no ↔ ono) edges unchanged:
	op, opo = opposite(t, p), opposite(t, po)
	opposites!(t, p=>po, c=>opo, o=>op)
	site!(t, n=>so, no=>s, po=>sn)
	anycorner!(t, so=>n, s=>no, sn=>po)
	return (no, c)
end

"""
    insert!(cornertable, triangle, site)

Inserts a new site by splitting the given triangle.
Returns the three corners formed at the new site.
"""
function insert!(t::CornerTable, q::Triangle, s::Site)
	c0, n0, p0 = corners(q)
	s0, s1, s2 = site(t, c0), site(t, n0), site(t, p0)
	# c0 and its opposite corner are unchanged
	on, op = opposite(t, n0), opposite(t, p0)
	q1, q2 = new_triangles!(t, 2)
	c1, n1, p1 = corners(q1)
	c2, n2, p2 = corners(q2)
	opposites!(t, c1=>on, c2=>op, n0=>p1, n1=>p2, n2=>p0)
	site!(t, c0=>s, c1=>s, c2=>s, n1=>s2, p1=>s0, n2=>s0, p2=>s1)
	anycorner!(t, s=>c0, s1=>n0, s2=>p0, s0=>n1)
	return (c0, c1, c2)
end

"""
    swapsites!(cornertable, site, site)

Swap site indices for those two sites.
"""
function swapsites!(t::CornerTable, s1::Site, s2::Site)
	s1 == s2 && return
	for c in star(t, s1); site!(t, c, s2); end
	for c in star(t, s2); site!(t, c, s1); end
	c1 = anycorner(t, s1)
	anycorner!(t, s1, anycorner(t, s2))
	anycorner!(t, s2, c1)
end

function swaptriangles!(t::CornerTable, q1::Triangle, q2::Triangle)
	q1 == q2 && return
	c1,n1,p1 = side(q1,1), side(q1,2), side(q1,3)
	c2,n2,p2 = side(q2,1), side(q2,2), side(q2,3)
	oc1,on1,op1 = opposite(t,c1), opposite(t,n1), opposite(t,p1)
	oc2,on2,op2 = opposite(t,c2), opposite(t,n2), opposite(t,p2)
	sc1,sn1,sp1 = site(t,c1), site(t,n1), site(t,p1)
	sc2,sn2,sp2 = site(t,c2), site(t,n2), site(t,p2)
	opposites!(t, c2=>oc1, n2=>on1, p2=>op1, c1=>oc2, n1=>on2, p1=>op2)
	sitecorner!(t, sc1=>c2, sn1=>n2, sp1=>p2, sc2=>c1, sn2=>n1, sp2=>p1)
end


# Site location functions ««2

"""
    locate(cornertable, points, point)

In a Voronoi diagram,
finds the index of the site closest to the given point.
"""
function locate_site(t::CornerTable{J}, points, point) where{J}
	# Returns the site index closest to this point
	s0 = Site(J(1)); p0 = points[Int(s0)]
	while true
		s1, p1 = s0, p0
		# find a closer neighbour if it exists
		for s in neighbours(t, s0)
			p = points[Int(s)]
			iscloser(point, p0, p) && continue
			s0, p0 = s, p
			break
		end
		s0 == s1 && return s0
	end
end

function locate_triangle(t::CornerTable{J}, points, point) where{J}
	q0 = Triangle(J(1));
	c = 0
	while true
		c+= 1
		@assert c ≤ 5
		v0 = ((points[Int(site(t, s))] for s in corners(q0))...,)
		if isleft(v0[1], point, v0[2])
			q0 = triangle(opposite(t, side(q0, 3)))
			continue
		elseif isleft(v0[2], point, v0[3])
			q0 = triangle(opposite(t, side(q0, 1)))
			continue
		elseif isleft(v0[3], point, v0[1])
			q0 = triangle(opposite(t, side(q0, 2)))
			continue
		else
			return q0
		end
	end
end


# Delaunay triangulation / Voronoi tessellation ««1
"""
    triangulate(points)

Returns a `CornerTable` whose sites are the Voronoi diagram of the input points
and whose triangles are its Delaunay triangulation.
"""
function triangulate(points)
	J = Int32
	N = length(points)
	N < 3 && return CornerTable{J}(undef, 0)
	# add  three extra points around the whole picture««
	# and build an initial dihedron
	m = maximum(abs(x) for p in points for x in p)
	push!(points, [0,-3m], [3m,2m], [-3m,2m])
	t = CornerTable{J}(
		StructArray((opposite=J[4,6,5,1,3,2], site=J(N).+J[1,2,3,1,3,2])),
		Vector{J}(undef, N+3))
	anycorner!(t, Site(N+1), Corner(J(1)))
	anycorner!(t, Site(N+2), Corner(J(2)))
	anycorner!(t, Site(N+3), Corner(J(3)))
  #»»
	# incrementally add all points ««
	for i in N:-1:1
		point = points[i]
		q0 = locate_triangle(t, points, point)
		(c0, c1, c2) = insert!(t, q0, Site(i))
		flip_rec(t, points, c0)
		flip_rec(t, points, c1)
		flip_rec(t, points, c2)
	end#»»
	# remove all superfluous triangles & sites ««
	k = ntriangles(t)
	for i in ntriangles(t):-1:1; q = Triangle{J}(i)
		any(>(N), Int.(sites(t, q))) || continue
		swaptriangles!(t, q, Triangle(k))
		k-= 1
	end
	ntriangles!(t, k)
	resize!(points, N)
	# »»
	return t
end

function flip_rec(t::CornerTable, points, c::Corner)
	# c is the corner (p,new site, a, b)
	# o is the opposite corner
	# if o is in circle of (a,b,p) then flip the corner c
	sc, sa, sb, so = site(t, c), site(t, next(c)), site(t, prev(c)),
		site(t, opposite(t, c))
	isincircle(points[Int.(SA[sc,sa,sb,so])]...) || return
	(c1, c2) = flip!(t, c)
	flip_rec(t, points, c1)
	flip_rec(t, points, c2)
end

# end »»1

c = Corner.([1,2,3])
t = CornerTable{Int}(Val(:tetrahedron))
t2 = CornerTable{Int}(StructArray((opposite=[4,6,5,1,3,2],site=[1,2,3,1,3,2])),
	[1,2,3])
vert = [SA[0,-3000],SA[3000,2000],SA[-3000,2000],SA[0,0],SA[10,0],SA[0,10]]
vert = [SA[0,0], SA[10,0], SA[0,10], SA[11,11]]

end
