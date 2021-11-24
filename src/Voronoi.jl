# [guibas stolfi]
module Voronoi
using StaticArrays
using LazyArrays
using FastClosures
using LinearAlgebra

import Base: Int

# Tools ««1
# Geometry ««2

@inline norm²(v) = v[1]^2+v[2]^2
@inline distance²(a,b) = norm²(a-b)
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
	m = SA[a[1] a[2] norm²(a); b[1] b[2] norm²(b); c[1] c[2] norm²(c)]
	return det(m) > 0
end

# Named indices ««2
struct NamedIndex{S,T<:Signed} i::T; end
@inline Int(i::NamedIndex) = i.i
@inline Base.convert(T::Type{<:Integer}, i::NamedIndex) = T(Int(i))
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
@NamedIndex Corner c
@NamedIndex Cell s
@NamedIndex Node t

# CornerTable ««1
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
"""
struct CornerTable{J <: Integer}
	opposite::Vector{J}
	apex::Vector{J}
	anycorner::Vector{J} # cell -> corner map
end

@inline CornerTable{J}(::UndefInitializer, ncells::Integer) where{J} =
	CornerTable{J}(J[], J[], zeros(J, ncells))
@inline CornerTable{J}() where{J} = CornerTable{J}(undef, 0)

# for broadcast:
@inline Base.Broadcast.broadcastable(x::CornerTable) = (x,)


# Accessors ««2
# Corner functions ««3
@inline next3(i::J) where{J<:Integer} = iszero(i%3) ? i-2 : i+1
@inline prev3(i::J) where{J<:Integer} = isone(i%3) ? i+2 : i-1
@inline next(c::Corner{J}) where{J} = Corner{J}(next3(Int(c)))
@inline prev(c::Corner{J}) where{J} = Corner{J}(prev3(Int(c)))

# LazyRow accepts only `Int` indexes:
@inline apex(t::CornerTable, c::Corner) = Cell(t.apex[Int(c)])
@inline apex!(t::CornerTable, c::Corner, s::Cell) = t.apex[Int(c)] = Int(s)
@inline apex!(t::CornerTable, l::Pair{<:Corner, <:Cell}...) =
	for (c,s) in l; apex!(t, c, s); end

@inline opposite(t::CornerTable, c::Corner) = Corner(t.opposite[Int(c)])
@inline opposite!(t::CornerTable, c::Corner, x::Corner) =
	t.opposite[Int(c)] = Int(x)
@inline opposites!(t::CornerTable, l::Pair{<:Corner,<:Corner}...) =
	for(c1, c2) in l; opposite!(t, c1, c2); opposite!(t, c2, c1); end

@inline right(t::CornerTable, c::Corner) = apex(t, next(c))
@inline left(t::CornerTable, c::Corner) = apex(t, prev(c))
@inline after(t::CornerTable, c::Corner) = next(opposite(t, next(c)))
@inline before(t::CornerTable, c::Corner) = prev(opposite(t, prev(c)))

# Node functions ««3
@inline nnodes(t::CornerTable) = length(t.opposite) ÷ 3
@inline nnodes!(t::CornerTable, n) = 
	(resize!(t.opposite, 3n); resize!(t.apex, 3n))
@inline node(c::Corner{J}) where{J} = Node{J}(fld1(Int(c), 3))

@inline side(q::Node{J}, i) where{J} = Corner{J}(3*Int(q)-3+i)
@inline corners(q::Node{J}) where{J} = Corner{J}.(3*Int(q) .- (2,1,0))
@inline cells(t::CornerTable, q::Node) =
	(apex(t, side(q,1)), apex(t, side(q,2)), apex(t, side(q,3)))
@inline triangles(t::CornerTable) =
	(((Int(apex(t, c)) for c in corners(q))...,)
		for q in Node.(1:nnodes(t)))
@inline adjacent(t::CornerTable, q::Node, i) =
	node(opposite(t, side(q, i)))
@inline adjacents(t::CornerTable, q::Node) =
	(adjacent(t,q,1), adjacent(t,q,2), adjacent(t,q,3))

function newnodes!(t::CornerTable{J}, k) where{J}
	n = nnodes(t)
	nnodes!(t, n+k)
	return Node{J}.(n+1:n+k)
end

# Cell functions ««3
@inline anycorner(t::CornerTable, s::Cell) = Corner(t.anycorner[Int(s)])
@inline anycorner!(t::CornerTable, s::Cell, c::Corner) =
	t.anycorner[Int(s)] = Int(c)
@inline anycorner!(t::CornerTable, l::Pair{<:Cell, <:Corner}...) =
	for (s,c) in l; anycorner!(t, s, c); end
@inline cellcorner!(t::CornerTable, l::Pair{<:Cell, <:Corner}...) =
	for(s,c) in l; anycorner!(t, s, c); apex!(t, c, s); end
@inline ncells(t::CornerTable) = length(t.anycorner)
@inline ncells!(t::CornerTable, n) = resize!(t.anycorner, n)
# Iterators ««2
struct CTIterator{S,X,J,T}
	table::CornerTable{J}
	start::T
	@inline CTIterator{S,X}(t::CornerTable{J}, s::T) where{S,X,J,T}=
		new{S,X,J,T}(t,s)
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
@inline star(t::CornerTable{J}, s::Cell) where{J} =
	CTIterator{star,Corner{J}}(t, anycorner(t, s))
@inline next(::CTIterator{star}, t, c) = after(t, c)
@inline stop(it::CTIterator{star}, t, c) = (c == it.start)

"""
    ring(cornertable, corner)

Iterates all corners left-adjacent to the indicated cell.
"""
@inline ring(t::CornerTable{J}, c::Corner) where{J} =
	CTIterator{ring,Corner{J}}(t, next(c))
@inline ring(t::CornerTable, s::Cell) = ring(t, anycorner(t, s))
@inline next(::CTIterator{ring}, t, c) = prev(opposite(t, c))
@inline stop(it::CTIterator{ring}, t, c) = (c == it.start)
"""
   neighbours(cornertable, cell)

Iterates all cells adjacent to the indicated cell.
"""
@inline neighbours(t::CornerTable{J}, s::Cell) where{J} =
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
	return CornerTable{J}(olist, slist, clist)
end

# Elementary modifications: flip!, insert! ««2
"""
    flip!(cornertable, corner)

Flips the quadrilateral delimited by `corner` and its opposite.
Returns the two corners (in order) replacing `corner`.
"""
function flip!(t::CornerTable, c::Corner)
	n, p, o = next(c), prev(c), opposite(t, c)
	s, sn, so, sp = apex(t, c), apex(t, n), apex(t, o), apex(t, p)
	no, po = next(o), prev(o)
	# we leave (n ↔ on) and (no ↔ ono) edges unchanged:
	op, opo = opposite(t, p), opposite(t, po)
	opposites!(t, p=>po, c=>opo, o=>op)
	apex!(t, n=>so, no=>s, po=>sn)
	anycorner!(t, so=>n, s=>no, sn=>po)
	return (no, c)
end

"""
    insert!(cornertable, node, cell)

Inserts a new cell by splitting the given node.
Returns the three corners formed at the new cell.
"""
function insert!(t::CornerTable, q::Node, s::Cell)
	c0, n0, p0 = corners(q)
	s0, s1, s2 = apex(t, c0), apex(t, n0), apex(t, p0)
	# c0 and its opposite corner are unchanged
	on, op = opposite(t, n0), opposite(t, p0)
	q1, q2 = newnodes!(t, 2)
	c1, n1, p1 = corners(q1)
	c2, n2, p2 = corners(q2)
	opposites!(t, c1=>on, c2=>op, n0=>p1, n1=>p2, n2=>p0)
	apex!(t, c0=>s, c1=>s, c2=>s, n1=>s2, p1=>s0, n2=>s0, p2=>s1)
	anycorner!(t, s=>c0, s1=>n0, s2=>p0, s0=>n1)
	return (c0, c1, c2)
end

"""
    swapcells!(cornertable, cell, cell)

Swap cell indices for those two cells.
"""
function swapcells!(t::CornerTable, s1::Cell, s2::Cell)#««
	s1 == s2 && return
	for c in star(t, s1); apex!(t, c, s2); end
	for c in star(t, s2); apex!(t, c, s1); end
	c1 = anycorner(t, s1)
	anycorner!(t, s1, anycorner(t, s2))
	anycorner!(t, s2, c1)
end#»»
function swapnodes!(t::CornerTable, q1::Node, q2::Node)#««
	q1 == q2 && return
	c1,n1,p1 = side(q1,1), side(q1,2), side(q1,3)
	c2,n2,p2 = side(q2,1), side(q2,2), side(q2,3)
	oc1,on1,op1 = opposite(t,c1), opposite(t,n1), opposite(t,p1)
	oc2,on2,op2 = opposite(t,c2), opposite(t,n2), opposite(t,p2)
	sc1,sn1,sp1 = apex(t,c1), apex(t,n1), apex(t,p1)
	sc2,sn2,sp2 = apex(t,c2), apex(t,n2), apex(t,p2)
	dc = typeof(c1)(3*(Int(q2)-Int(q1)))
	for (x,y) in (c2=>oc1, n2=>on1, p2=>op1, c1=>oc2, n1=>on2, p1=>op2)
		q = node(y)
		(q == q1) && (y+= dc); (q == q2) && (y-= dc)
		opposites!(t, x=>y)
# 	opposites!(t, c2=>oc1, n2=>on1, p2=>op1, c1=>oc2, n1=>on2, p1=>op2)
	end
	cellcorner!(t, sc1=>c2, sn1=>n2, sp1=>p2, sc2=>c1, sn2=>n1, sp2=>p1)
end#»»

# Delaunay triangulation / Voronoi tessellation ««1
# Tree traversal ««2
function badnodes(t::CornerTable{J}, points, point) where{J}
	stack = [locate_node(t, points, point)]
	tree = empty(stack)
	while !isempty(stack)
		q = pop!(stack)
		isone(Int(q)) && continue # this is the reversed triangle
		q ∈ tree && continue
		isincircle(t, q, points, point) || continue
		push!(tree, q)
		push!(stack, adjacents(t, q)...)
	end
	return tree # the tree is usually small, no need to sort it
end
# end#»»
function tree_boundary(t::CornerTable{J}, tree) where{J}#««
	# returns the list of half-edges pointing *into* the tree,
	# cyclically ordered around the tree (random start).
	boundary = sizehint!(Corner{J}[], length(tree)+2)
	for c in corners(first(tree))
		stack = [c]
		while !isempty(stack)
			c = pop!(stack)
			o = opposite(t, c)
			if node(o) ∈ tree
				push!(stack, prev(o), next(o))
			else
				push!(boundary, o)
			end
		end
	end
	return boundary
end#»»
function star_transform!(t::CornerTable, tree, s)#««
	# replaces all nodes in `tree` by a star shape around `cell`,
	# adding two new nodes in the process.
	boundary = tree_boundary(t, tree)
	push!(tree, newnodes!(t, 2)...)
	n = side(last(tree), 2)
	for (q, o) in zip(tree, boundary)
		c = side(q, 1)
		p = side(q, 3)
		opposites!(t, c=>o, p=>n)
		n = side(q, 2)
		cellcorner!(t, s=>c, right(t,o)=>p, left(t,o)=>n)
	end
end#»»

# Cell location functions ««2
"""
    locate_cell(cornertable, points, point)

In a Voronoi diagram,
finds the index of the cell closest to the given point.
"""
function locate_cell(t::CornerTable{J}, points, point) where{J}
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

function locate_node(t::CornerTable{J}, points, point) where{J}
	q0 = Node(J(1));
	c = 0
	while true
		c+= 1
		@assert c ≤ 5
		v0 = ((points[Int(apex(t, s))] for s in corners(q0))...,)
		if isleft(v0[1], point, v0[2])
			q0 = node(opposite(t, side(q0, 3)))
			continue
		elseif isleft(v0[2], point, v0[3])
			q0 = node(opposite(t, side(q0, 1)))
			continue
		elseif isleft(v0[3], point, v0[1])
			q0 = node(opposite(t, side(q0, 2)))
			continue
		else
			return q0
		end
	end
end

# Triangulation ««2
"""
    triangulate(points)

Returns a `CornerTable` whose cells are the Voronoi diagram of the input points
and whose nodes are its Delaunay triangulation.
"""
function triangulate(points)
	J = Int32
	N = length(points)
	N < 3 && return CornerTable{J}(undef, 0)
	# add  three extra points around the whole picture««
	# and build an initial dihedron
	m = maximum(abs(x) for p in points for x in p)
	push!(points, [0,-3m], [3m,2m], [-3m,2m])
	t = CornerTable{J}(J[4,6,5,1,3,2],J(N).+J[1,3,2,1,2,3], Vector{J}(undef, N+3))
	anycorner!(t, Cell(N+1), Corner(J(1)))
	anycorner!(t, Cell(N+2), Corner(J(2)))
	anycorner!(t, Cell(N+3), Corner(J(3)))
  #»»
	# incrementally add all points ««
	for i in N:-1:1
		point = points[i]
		tree = badnodes(t, points, point)
		star_transform!(t, tree, Cell(i))
	end#»»
	# remove all superfluous nodes & cells ««
	# the nodes are sorted this way:
	# - inner nodes
	# - convex hull
	# - 1 backwards outer node
	k = nnodes(t)
	swapnodes!(t, Node(1), Node(k))
	k-= 1
	for i in nnodes(t)-1:-1:1; q = Node{J}(i)
		any(>(N), Int.(cells(t, q))) || continue
		swapnodes!(t, q, Node(k))
		k-= 1
	end
# 	nnodes!(t, k)
	resize!(points, N)
	# »»
	return t
end

function isincircle(t::CornerTable, q::Node, points, point)
	s1, s2, s3 = Int.(cells(t, q))
# 	println("$s1, $s2, $s3: $(points[[s1,s2,s3]])")
	return isincircle(points[s1], points[s2], points[s3], point)
end

function flip_rec(t::CornerTable, points, c::Corner)
	# c is the corner (p,new cell, a, b)
	# o is the opposite corner
	# if o is in circle of (a,b,p) then flip the corner c
	sc, sa, sb, so = apex(t, c), apex(t, next(c)), apex(t, prev(c)),
		apex(t, opposite(t, c))
	isincircle(points[Int.(SA[sc,sa,sb,so])]...) || return
	(c1, c2) = flip!(t, c)
	flip_rec(t, points, c1)
	flip_rec(t, points, c2)
end

# end »»1

vert = [SA[0,0], SA[10,0], SA[0,10], SA[11,11]]

end
