module CornerTables
# Tools««1
# Named indices ««2
# These are “strongly typed” integers, used to prevent confusion when
# writing accessor functions for meshes.
struct NamedIndex{S,T<:Signed}<:Integer i::T; end
@inline int(i::NamedIndex) = i.i
@inline Integer(i::NamedIndex) = int(i)
@inline (T::Type{<:Integer})(i::NamedIndex) = T(Integer(i))
@inline NamedIndex{S,T}(i::NamedIndex{S,T}) where{S,T<:Signed}=
	NamedIndex{S,T}(int(i))
@inline NamedIndex{S}(x) where{S} = NamedIndex{S,typeof(x)}(x)

# @inline Base.convert(T::Type{<:NamedIndex},i::Integer) = T(i)
@inline Base.zero(A::Type{NamedIndex{S,T}}) where{S,T} = A(zero(T))
@inline Base.one(A::Type{NamedIndex{S,T}}) where{S,T} = A(one(T))
@inline Base.zero(x::T) where{T<:NamedIndex} = zero(T)
@inline Base.one(x::T) where{T<:NamedIndex} = one(T)
for op in (:+, :-)
	@eval @inline Base.$op(a::T,b::T) where{T<:NamedIndex} = T($op(int(a),int(b)))
end
for op in (:<, :(<=), :(==))
	@eval @inline Base.$op(x::NamedIndex{S}, y::NamedIndex{S}) where{S} =
		$op(int(x), int(y))
end

macro NamedIndex(name,S) quote
	const $(esc(name)) = NamedIndex{$(QuoteNode(name))}
	Base.show(io::IO, ::Type{<:$(esc(name))}) = print(io, $(string(name)))
@inline Base.show(io::IO, i::NamedIndex{$(QuoteNode(name))}) =
	print(io, $(QuoteNode(S)), int(i))
end end

# AbstractTriangulation ««1
# Nodes and arrows ««2
@NamedIndex Arrow a
@NamedIndex Edge e
@NamedIndex Cell c
@NamedIndex Node q

@inline next3(i::J) where{J<:Integer} = iszero(i%3) ? i-J(2) : i+J(1)
@inline prev3(i::J) where{J<:Integer} = isone(i%3) ? i+J(2) : i-J(1)
@inline next(e::Arrow{J}) where{J} = Arrow{J}(next3(int(e)))
@inline prev(e::Arrow{J}) where{J} = Arrow{J}(prev3(int(e)))

@inline node(e::Arrow{J}) where{J} = Node{J}(fld1(int(e), 3))
@inline side(q::Node{J}, i) where{J} = Arrow{J}(3*int(q)-3+i)
@inline arrows(q::Node{J}) where{J} = let a=J(3)*int(q)
	Arrow(a-J(2)):Arrow(a); end
# Data type««2
abstract type AbstractTriangulation{J} end
# Objects of this type must define the following methods:
#  - triangulation() - returning a CornerTable, or something with the
for f in (:tail, :tail!, :opposite, :opposite!, :anyarrow, :anyarrow!,
	:nnodes, :nnodes!, :narrows, :ncells, :ncells!)
	@eval @inline $f(t::AbstractTriangulation, args...;kwargs...) =
		$f(triangulation(t), args...;kwargs...)
end

# Methods ««2
@inline head(t::AbstractTriangulation, e::Arrow) = tail(t, opposite(t, e))
@inline tail!(t::AbstractTriangulation, l::Pair{<:Arrow, <:Cell}...) =
	for (e,s) in l; tail!(t, e, s); end
@inline opposites!(t::AbstractTriangulation, l::Pair{<:Arrow,<:Arrow}...) =
	for(e1, e2) in l; opposite!(t, e1, e2); opposite!(t, e2, e1); end
@inline anyarrow!(t::AbstractTriangulation, l::Pair{<:Cell, <:Arrow}...) =
	for (s,e) in l; anyarrow!(t, s, e); end
@inline arrowsfrom!(t::AbstractTriangulation, l::Pair{<:Cell, <:Arrow}...) =
	for(s,e) in l; anyarrow!(t, s, e); tail!(t, e, s); end

@inline left(t::AbstractTriangulation, e::Arrow) = tail(t, prev(e))
@inline right(t::AbstractTriangulation, e::Arrow) = left(t, opposite(t, e))
# @inline right(t::AbstractTriangulation, e::Arrow) = tail(t, next(e))
# @inline left(t::AbstractTriangulation, e::Arrow) = tail(t, prev(e))
# @inline base(t::AbstractTriangulation, e::Arrow) = (left(t,e),right(t,e))
# `after`, `before`: next/previous arrows with same tail
# @inline after(t::AbstractTriangulation, e::Arrow) = next(opposite(t, next(e)))
# @inline before(t::AbstractTriangulation, e::Arrow) = prev(opposite(t, prev(e)))
@inline after(t::AbstractTriangulation, e::Arrow) = opposite(t, prev(e))
@inline before(t::AbstractTriangulation, e::Arrow) = next(opposite(t, e))

@inline cell(t::AbstractTriangulation, q::Node, i) = tail(t, side(q,i))
@inline triangle(t::AbstractTriangulation, q::Node) =
	(cell(t,q,1), cell(t,q,2), cell(t,q,3))
@inline lastcell(t::AbstractTriangulation) = Cell(ncells(t))
@inline eachcell(t::AbstractTriangulation) = Base.OneTo(lastcell(t))
@inline lastnode(t::AbstractTriangulation) = Node(nnodes(t))
@inline eachnode(t::AbstractTriangulation) = Base.OneTo(lastnode(t))
@inline lastarrow(t::AbstractTriangulation) = Arrow(narrows(t))
@inline eacharrow(t::AbstractTriangulation) = Base.OneTo(lastarrow(t))
@inline alltriangles(t::AbstractTriangulation) =
	(triangle(t,q) for q in eachnode(t))
@inline adjnode(t::AbstractTriangulation, q::Node, i) =
	node(opposite(t, side(q, i)))
@inline adjnodes(t::AbstractTriangulation, q::Node) =
	(adjnode(t,q,1), adjnode(t,q,2), adjnode(t,q,3))

@inline function newnodes!(t::AbstractTriangulation{J}, k) where{J}
	n = nnodes(t)
	nnodes!(t, n+k)
	return Node(J(n+1)):Node(J(n+k))
end
@inline function newcells!(t::AbstractTriangulation{J}, k) where{J}
	n = ncells(t)
	ncells!(t, n+k)
	return Cell(J(n+1)):Cell(J(n+k))
end

# Iterators ««2
struct ATIterator{S,X,J,T,A<:AbstractTriangulation}
	# S is an identifier
	# X is the return type
	# J is the index type
	# A, T are field types
	table::A
	start::T
	@inline ATIterator{S,X}(t::A, s::T) where{S,X,J,T,A<:AbstractTriangulation{J}}=
		new{S,X,J,T,A}(t,s)
end
@inline Base.eltype(::ATIterator{S,X}) where{S,X} = X
@inline Base.keys(g::Base.Generator{<:ATIterator}) = g.iter
@inline Base.keys(g::Base.Generator{<:Base.Generator{<:ATIterator}}) =
	g.iter.iter
@inline Base.IteratorSize(::ATIterator) = Base.SizeUnknown()
@inline Base.iterate(it::ATIterator) = (it.start, it.start)
@inline function Base.iterate(it::ATIterator, s)
	s = next(it, it.table, s)
	stop(it, it.table, s) && return nothing
	return (s, s)
end
@inline stop(it::ATIterator,_,e) = (e==it.start)
"""
    star(cornertable, cell)

Iterates all arrows sharing the same tail cell.
"""
@inline star(t::AbstractTriangulation{J}, e::Arrow) where{J} =
	ATIterator{star,Arrow{J}}(t, e)
@inline star(t::AbstractTriangulation, c::Cell) = star(t, anyarrow(t,c))
@inline next(::ATIterator{star}, t, e) = after(t, e)

# same as `star`, but in reverse order:
@inline revstar(t::AbstractTriangulation{J}, e::Arrow) where{J} =
	ATIterator{revstar,Arrow{J}}(t, e)
@inline revstar(t::AbstractTriangulation, c::Cell) = revstar(t, anyarrow(t,c))
@inline next(::ATIterator{revstar},t,e) = before(t,e)

@inline Base.reverse(it::ATIterator{star,T}) where{T} =
	ATIterator{revstar,T}(it.table, it.start)
@inline Base.reverse(it::ATIterator{revstar,T}) where{T} =
	ATIterator{star,T}(it.table, it.start)
function arrow(t::AbstractTriangulation, c1::Cell, c2::Cell)
	for e in star(t, c1)
		head(t, e) == c2 && return e
	end
	return nothing
end

"""
    ring(cornertable, cell)

Iterates all arrows left-adjacent to the indicated cell.
Equivalent to next.(star(...)).
"""
@inline ring(t::AbstractTriangulation{J}, e::Arrow) where{J} =
	ATIterator{ring,Arrow{J}}(t, next(e))
@inline ring(t::AbstractTriangulation, c::Cell) = ring(t, anyarrow(t, c))
@inline next(::ATIterator{ring}, t, e) = prev(opposite(t, e))
"""
   neighbours(cornertable, cell)

Iterates all cells adjacent to the indicated cell.
Equivalent to right.(ring(...)).
"""
@inline neighbours(t::AbstractTriangulation{J}, c::Cell) where{J} =
	(tail(t, x) for x in ring(t, anyarrow(t, c)))

# Elementary modifications ««2
"""
    split!(triangulation, edge)

Splits an edge in two, creating two new triangles.
Returns (newtriangle1, newtriangle2, newcell).
"""
function split!(t::AbstractTriangulation, e::Arrow)
	o = opposite(t, e)
	c1, c2, cl, cr = tail(t, e), tail(t, o),
		tail(t, prev(e)), tail(t, prev(o))
	(c3,) = newcells!(t, 1)
# 	c3 = pushpoint!(t, (point(t, c1)+point(t, c2))/2)
	x, y = newnodes!(t, 2)
	opposites!(t, side(x,1) => side(y,1),
		side(x,2) => opposite(t, next(e)), side(x,3) => next(e),
		side(y,2) => prev(o), side(y,3) => opposite(t, prev(o)))
	tail!(t, o=>c3, next(e)=>c3,
		side(x,1) => c3, side(x,2) => c2, side(x,3) => cl,
		side(y,1) => c2, side(y,2) => c3, side(y,3) => cr)
	anyarrow!(t, c2=>side(y,1), c3=>o)
	return (x, y, c3)
end
# FIXME: check that flip! is still up-to-date
"""
    flip!(triangulation, corner)

Flips the quadrilateral delimited by `corner` and its opposite.
Returns the two arrows (in order) replacing `corner`.
"""
function flip!(t::AbstractTriangulation, e::Arrow)#««
	n, p, o = next(e), prev(e), opposite(t, e)
	s, sn, so, sp = tail(t, e), tail(t, n), tail(t, o), tail(t, p)
	no, po = next(o), prev(o)
	# we leave (n ↔ on) and (no ↔ ono) edges unchanged:
	op, opo = opposite(t, p), opposite(t, po)
	opposites!(t, p=>po, e=>opo, o=>op)
	tail!(t, n=>so, no=>s, po=>sn)
	anyarrow!(t, so=>n, s=>no, sn=>po)
	return (no, e)
end#»»
"""
    insert!(triangulation, node, cell)

Inserts a new cell by splitting the given node.
Returns the three corners formed at the new cell.
"""
function insert!(t::AbstractTriangulation, q::Node, c::Cell)#««
	e0, n0, p0 = arrows(q)
	c0, c1, c2 = tail(t, e0), tail(t, n0), tail(t, p0)
	# e0 and its opposite corner are unchanged
	on, op = opposite(t, n0), opposite(t, p0)
	q1, q2 = newnodes!(t, 2)
	e1, n1, p1 = arrows(q1)
	e2, n2, p2 = arrows(q2)
	opposites!(t, e1=>on, e2=>op, n0=>p1, n1=>p2, n2=>p0)
	tail!(t, e0=>c, e1=>c, e2=>c, n1=>c2, p1=>c0, n2=>c0, p2=>c1)
	anyarrow!(t, c=>e0, c1=>n0, c2=>p0, c0=>n1)
	return (e0, e1, e2)
end#»»
"swap cell indices for those two cells."
function swapcells!(t::AbstractTriangulation, c1::Cell, c2::Cell)#««
	c1 == c2 && return
	for e in star(t, c1); tail!(t, e, c2); end
	for e in star(t, c2); tail!(t, e, c1); end
	e1 = anyarrow(t, c1)
	anyarrow!(t, c1, anyarrow(t, c2))
	anyarrow!(t, c2, e1)
end#»»
"moves cell c1 to (undefined) position c2, overwriting c2 in the process"
function movecell!(t::AbstractTriangulation, c1::Cell, c2::Cell)#««
	c1 == c2 && return
	for e in star(t, c1); tail!(t, e, c2); end
	anyarrow!(t, c2, anyarrow(t, c1))
end#»»
"swaps two nodes in the triangulation"
function swapnodes!(t::AbstractTriangulation, q1::Node, q2::Node)#««
	q1 == q2 && return
	e1,n1,p1 = side(q1,1), side(q1,2), side(q1,3)
	e2,n2,p2 = side(q2,1), side(q2,2), side(q2,3)
	oc1,on1,op1 = opposite(t,e1), opposite(t,n1), opposite(t,p1)
	oc2,on2,op2 = opposite(t,e2), opposite(t,n2), opposite(t,p2)
	sc1,sn1,sp1 = tail(t,e1), tail(t,n1), tail(t,p1)
	sc2,sn2,sp2 = tail(t,e2), tail(t,n2), tail(t,p2)
	dc = typeof(e1)(3*(int(q2)-int(q1)))
	for (x,y) in (e2=>oc1, n2=>on1, p2=>op1, e1=>oc2, n1=>on2, p1=>op2)
		q = node(y)
		(q == q1) && (y+= dc); (q == q2) && (y-= dc)
		opposites!(t, x=>y)
# 	opposites!(t, e2=>oc1, n2=>on1, p2=>op1, e1=>oc2, n1=>on2, p1=>op2)
	end
	arrowsfrom!(t, sc1=>e2, sn1=>n2, sp1=>p2, sc2=>e1, sn2=>n1, sp2=>p1)
end#»»

# CornerTable««1
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
	# corner ≈ 6V, cell ≈ 2V
	# total ≈ 14V
	opposite::Vector{J} # corner->corner
	tail::Vector{J}     # corner->cell
	anyarrow::Vector{J}# cell -> corner
end

@inline CornerTable{J}(::UndefInitializer, ncells::Integer) where{J} =
	CornerTable{J}(J[], J[], zeros(J, ncells))
@inline CornerTable{J}() where{J} = CornerTable{J}(undef, 0)
@inline triangulation(t::CornerTable) = t

# Accessors ««2
# Arrow functions ««3
@inline tail(t::CornerTable, e::Arrow) = Cell(t.tail[int(e)])
@inline tail!(t::CornerTable, e::Arrow, c::Cell) = t.tail[int(e)] = int(c)

@inline opposite(t::CornerTable, e::Arrow) = Arrow(t.opposite[int(e)])
@inline opposite!(t::CornerTable, e::Arrow, x::Arrow) =
	t.opposite[int(e)] = int(x)

# Node functions ««3
@inline narrows(t::CornerTable{J}) where{J} = J(length(t.opposite))
@inline nnodes(t::CornerTable{J}) where{J} = narrows(t) ÷ J(3)
@inline nnodes!(t::CornerTable, n) = 
	(resize!(t.opposite, 3n); resize!(t.tail, 3n); t)

# Cell functions ««3
@inline anyarrow(t::CornerTable, c::Cell) = Arrow(t.anyarrow[int(c)])
@inline anyarrow!(t::CornerTable, c::Cell, e::Arrow) =
	t.anyarrow[int(c)] = int(e)
@inline ncells(t::CornerTable{J}) where{J} = J(length(t.anyarrow))
@inline ncells!(t::CornerTable, n) = resize!(t.anyarrow, n)

# Constructor from triangle list ««2
function CornerTable{J}(triangles) where{J}
	nf = length(triangles); nc = 3nf
	nv = 0
	vlist = sizehint!(J[], nc)
	olist = Vector{J}(undef, nc)
	for q in triangles, v in q[1:3]
		(v > nv) && (nv = v)
		push!(vlist, v)
	end
	clist = Vector{J}(undef, nv)
	# edgeto[v] = list of all edges pointing to `v`
	edgeto = [ sizehint!(J[], 6) for _ in 1:nv ]
	for (e, v) in pairs(vlist)
		clist[v] = e
		push!(edgeto[v], prev3(e))
	end
	for (e, v1) in pairs(vlist)
		v2 = vlist[next3(e)]
		for o in edgeto[v1]
			vlist[o] == v2 || continue
			olist[e], olist[o] = o, e
		end
	end
	return CornerTable{J}(olist, vlist, clist)
end

# Displaying & debugging ««1
function showall(io::IO, t::AbstractTriangulation)
	println("\e[33;1m", nnodes(t), " nodes:\e[m")
	for q in eachnode(t); shownode(io, t, q, "\n"); end
	println("\e[31;1m", ncells(t), " cells:\e[m")
	for c in eachcell(t); showcell(io, t, c, "\n"); end
end

function shownode(io::IO, t::AbstractTriangulation, q::Node, s = "\n")
	print(io, "\e[33m", q, triangle(t,q), "\e[m ", s)
	for i in (1,2,3)
		e = side(q, i); o = opposite(t, e); oo = opposite(t,o)
		oo ≠ e && println(io, "  \e[31;7m opposite($o) = $oo, should be $e\e[m")
		e1 = next(e); t1 = tail(t, e1); to = tail(t, o)
		t1 ≠ to && println(io, "  \e[31;7m tail($o) = $to, should be $t1\e[m")
	end
end
function showcell(io::IO, v::AbstractTriangulation, c::Cell, s = "\n")
	print(io, "\e[31m star(", c, ")\e[m:");
	for e in star(v, c)
		print(io, " ", e, "(", node(e), ")→", head(v,e))
	end
	print(io, s)
	for e in star(v,c)
		c1 = tail(v, e)
		c1 ≠ c && println(io, "  \e[31;7m tail($e) = $c1, should be $c\e[m")
	end
end

#»»1
export AbstractTriangulation, CornerTable, Arrow, Cell, Edge, Node, int
export tail, tail!, head, opposite, opposite!, anyarrow, anyarrow!
export narrows, nnodes, nnodes!, ncells, ncells!
export next, prev, node, side, arrows
export opposites!, arrowsfrom!
export right, left, after, before, cell, triangle
export lastcell, eachcell, lastnode, eachnode, lastarrow, eacharrow
export alltriangles, adjnode, adjnodes, newnodes!
export star, ring
export swapcells!, movecell!, swapnodes!
end
