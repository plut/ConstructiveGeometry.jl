# FIXME: upon deleting a segment, newinter between new neighbours
# FIXME: due to float inexactitude issues, the ordering of segments in the
# active line is not completely consistent; it may happen that the
# intersection point `p` of two consecutive segments `s1, s2` is
# approximated, so that the ordering of `s1, s2` at `p` is wrong.
#
# For now we use a lazy solution, which works only for small graphs: we
# remove segments using a O(n) loop iterating on all segments.
# Two possible 'clean' solutions:
# (1) use double precision (e.g. `DoubleFloats.jl`) for all points and
# hope that the ordering becomes consistent;
# (2) use an indexed structure for storing active segments, making
# segment removal O(log n) given its tree pointer
#
# (1) was tried but, at least “out-of-the-box”, does not work.
# (2) is much more complicated and probably more expensive for the short
# structures we handle anyway.
#
# [Melhorn, Näher 1994] _Implementation of a sweep line algorithm
# for the straight line segment intersection problem_
# http://people.scs.carleton.ca/~michiel/lecturenotes/ALGGEOM/bentley-ottmann.pdf
# J. Hobby, Practical segment intersection with finite precision output
# https://www.sciencedirect.com/science/article/pii/S0925772199000218/pdf?isDTMRedir=true&download=true
# J. Hersberger, Stable snap rounding
# https://www.sciencedirect.com/science/article/pii/S0925772112001551
# J. Hershberger, Improved output-sensitive snap rounding
# https://link.springer.com/content/pdf/10.1007/s00454-007-9015-0.pdf
# P. Guigue, thèse
# https://tel.archives-ouvertes.fr/tel-00471447/document
#
# Must return enough data to
#  - tweak existing points (changing coordinates + lift is enough)
#  - merge corresponding points (?)
#  - link all of these
module SegmentIntersections
using StructArrays
using DataStructures
using StaticArrays
using FastClosures
# Geometry tools««1
@inline det2(a,b)=a[1]*b[2]-a[2]*b[1]
@inline det2(a,b,c)=det2(b-a,c-a)
"returns the point (p1p2) ∩ (p3p4), or `nothing` if no intersection"
function inter(p1, p2, p3, p4)#««
	@assert p1 < p2
	println("inter: $p1 $p2 $p3 $p4")
	# special case: (p1p2) is vertical, we *need* the intersection point
	# to be on that same vertical line.
	(v12, v13, v14) = (p2-p1, p3-p1, p4-p1)
	d2 = det2(v13, v14)
	d3 = det2(v12, v14)
	d4 = det2(v12, v13)
	println("v2,v3,v4=$v12 $v13 $v14")
	println("d2=$d2, d3=$d3, d4=$d4")
	iszero(d3) && iszero(d4) && return min(p2,p4)
	sign(d3) == sign(d4) && return nothing
	d = d3 - d4
# 	if d == 0 # two segments are parallel or confounded
# 		d2 == 0 && println("d2 is zero")
# 		d2 == 0 && return min(p2, p4)
# 		return nothing
# 	end
# 	d == 0 && return nothing
	d1 = d2 - d
	println("d2=$d2, d1=$d1")
	sign(d1) == sign(d2) && return nothing
# 	abs(d) ≤ ε && return nothing
	z= p1+d2/d*v12
		return z
	return nothing
end#»»
"""
compares two segments at their intersection with the sweepline either
before of after n-the point.
"""
function cmpseg(p, (a1, b1), (a2, b2), before=false)#««
	(u1, v1, u2, v2) = (x - p for x in (a1, b1, a2, b2))
	d1 = det2(v1, u1) # d1 > 0 iff (a1,b1) crosses *above* p
	d2 = det2(v2, u2) # ditto
	dx1 = v1[1]-u1[1]
	dx2 = v2[1]-u2[1]
# 	println("\e[31;1mp=$p, a1=$a1\e[m, u1=$u1")
# 	print(""" cmpseg at $p: ($a1,$b1) and ($a2,$b2)
#   seg1 crosses at $d1/$dx1, seg2 at $d2/$dx2
# """)
# 	println(u1, v1, u2, v2)
	if dx1 == 0 # segment 1 is vertical through
		dx2 == 0 && return cmp(b1[2], b2[2]) # some arbitrary, consistent test...
		d2 == 0 && return before ? -1 : 1
		return cmp(0, d2)
	elseif dx2 == 0
		d1 == 0 && return before ? 1 : -1
		return cmp(d1, 0)
	end
	# segment (a1, b1) crosses at d1/dx1, etc.
	c = cmp(d1*dx2, d2*dx1)
	c ≠ 0 && return c
	# both segments cross at the same point: check directions
	d = det2(v1-u1, v2-u2)
	s = before ? (d1 < 0) : (d1 ≤ 0)
	return (-1)^s*Int(sign(d))
end#»»

# L∞ pixels ««1
@inline roundreal(x, ε) = iszero(ε) ? x : ε*floor(Int,x/ε+1/2)
@inline roundLinf(point, ε) = roundreal.(point, ε)
	
"returns the point where segment (p,q) exits the pixel, or `nothing`."
function pixel_exit(pixel, ε, p, q)#««
	c = pixel .+ ε/2 # upper-right corner
	(cp, cq) = (p-c, q-c)
	@assert any(cp .< 0)
	all(cq .< 0) && return nothing
	d2 = det2(cp, cq)
	if d2 ≤ 0 # exits through the top part
		@assert cp[2] < 0
		@assert cq[2] ≥ 0
		# d = d3-d4 = ε*(cp[2]-cq[2])
		# z = c + d2/d*(-ε,0) = c+d2/ε*(cp[2]-cq[2])*(-ε,0)
		x = d2/(cp[2]-cq[2])
		@assert x ≤ ε
		return c .- (x, 0)
	end
	# test for exit through the right part:
	r = ε*(cq[1]-cp[1])
	d2 ≤ r && return c .- (0, d2/(cq[1]-cp[1]))
	# exit through the bottom part:
	# new origin c' = c-(0,ε)
	# d'2 = det(c',p,q) = d2-r
	# v12 = (-ε, 0)
	# d3 = det(v12,v14+(0,ε)) = |-ε 0;cq[1] cq[2]+ε| = -ε*cq[2]-ε²
	# likewise d4=-ε*cp[2]-ε², so that d=d3-d4=ε*(cq[2]-cp[2])
# 	d = ε*(cq[2]-cp[2])
	x = (d2-r)/(cp[2]-cq[2])
	@assert x ≤ ε
	return c .- (x,ε)
end#»»
# SnapRounding structure««1
"""
 - `active`: all segments crossing sweepline, in order
 - `segmentsfrom`: all segments starting from this point
   (≈ events associated with this point in the queue)
 - `queue`: heap of future points to examine (sorted lex.)
"""
struct SnapRounding{P,I<:Integer,E,EV}
	# P: point type
	# E: ε type
	# I: index type
	# EV: StructVector (constructor computed) for events
	# active: (start, end, last visited point)
	events::EV
	segments::Vector{MVector{3,I}} # (a,b,last-active-point)
	ε::E
	queue::Vector{I}
	active::Vector{I}
	pointid::Dict{P,I}
	pixel::Dict{P,I}
	links::Set{NTuple{2,I}}
	@inline SnapRounding{P,I}(events::EV, segments, ε::E, queue
			) where{P,I,EV,E} =
		new{P,I,E,EV}(events, segments, ε, queue, [], Dict(), Dict(), Set())
	@inline SnapRounding{P,I}(points, segments, ε) where{P,I} =
		SnapRounding{P,I}(StructVector((
			# start: segments starting at this point (original segment start)
			# pass: segments passing through this intersection point
			# reinsert: segments reinserted at this point
				point=points,
				start = [collect(searchsorted(segments, i;by=first))
					for i in eachindex(points)],
				pass = [I[] for _ in points ],
				reinsert=[ I[] for _ in points ],
			)),
			[(s[1], s[2], s[1]) for s in segments], # segments
			ε,  # ε
			heapify!([1:length(points)...], Base.Order.By(@closure i->points[i])),
			)
end

@inline event(snap::SnapRounding, i) = LazyRow(snap.events, i)
@inline point(snap::SnapRounding, i) = event(snap, i).point
@inline point(snap::SnapRounding, i, j) = point(snap, snap.segments[i][j])
@inline points(snap::SnapRounding, l...)= (x.point for x in snap.events[[l...]])
@inline ord(snap::SnapRounding) = Base.Order.By(@closure i->point(snap, i))
@inline Base.round(snap::SnapRounding, point) = roundLinf(point, snap.ε)

function SnapRounding(points, segments, ε)#««
# 	println("\e[48;5;57m****************************\e[m"*"\n"^30)
	# make all segments left-to-right
	points .= roundLinf.(points, ε)
	for (i, s) in pairs(segments)
		(points[s[1]] >  points[s[2]]) && (segments[i] = (s[2], s[1]))
	end
	unique!(sort!(segments))
	snap = SnapRounding{eltype(points),Int}(points, segments, ε)
	for i in reverse(eachindex(points))
		p = point(snap, i)
		snap.pointid[p] = i
		snap.pixel[round(snap, p)] = i
	end
	for (s, (a,b)) in pairs(snap.segments)
		push!(event(snap, b).pass, s)
	end
	return snap
end#»»

function point!(snap::SnapRounding, point)
	@assert !isnan(point[1])
	n = get!(snap.pointid, point, length(snap.events)+1)
	if n == length(snap.events)+1
		println("\e[32;7;1m create event $n=$point\e[m")
		push!(snap.events, (point=point, start=[], reinsert=[], pass=[]))
		heappush!(snap.queue, n, ord(snap))
	end
	return n
end
@inline points!(snap::SnapRounding, plist...)= [point!(snap, p) for p in plist]

function pixel!(snap::SnapRounding, point)
	center = round(snap, point)
	n = get!(snap.pixel, center, length(snap.pixel)+1)
	return n
end


# Adding an intersection point««1
"""
intersects active[j] and active[j+1] if needed.
pushes the intersection point in the queue (with appropriate segmentsfrom[]).
"""
function newinter(snap::SnapRounding, n, j)#««
	s1, s2 = snap.active[j:j+1]
	(a1, b1) = snap.segments[s1]
	(a2, b2) = snap.segments[s2]
	println("\e[33;1m intersect active[$j]=$s1=($a1,$b1) and [$(j+1)]=$s2=($a2,$b2)\e[m")
	println(point(snap,a1))
	println(point(snap,b1))
	println(point(snap,a2))
	println(point(snap,b2))
	@assert point(snap, a1) < point(snap, b1)
	@assert point(snap, a2) < point(snap, b2)
	a1 == a2 && return false
	b1 == b2 && return false
	println("$(point(snap, a1)) $(point(snap, b1))")
	p = inter(points(snap, a1, b1, a2, b2)...)
	println(collect(points(snap, a1,b1,a2,b2)))
	println("  found point p=$p")
# 	@assert (a1,b1,a2,b2) ≠ (5,7,17,16)
	p == nothing && return false
	p > point(snap, n) || return false
	c = point!(snap, p)
	println("   \e[33;7;1m ($a1,$b1) ∩ ($a2,$b2) found point $c=$p\e[m")
	println(" pushing $s1 and/or $s2 on point $c")
	ev = event(snap, c)
	println(" before: pass($c)=$(Tuple.(snap.segments[ev.pass]))")
	c ≠ a1 && c ≠ b1 && push!(ev.pass, s1)
	c ≠ a2 && c ≠ b2 && push!(ev.pass, s2)
	println(" after: pass($c)=$(Tuple.(snap.segments[ev.pass]))")
end#»»

# Inserting new segments ««1
@inline cmpf(snap::SnapRounding, n, before) =
	@closure (s1, s2) -> cmpseg(point(snap, n),
		(point(snap, s1, 1), point(snap, s1, 2)),
		(point(snap, s2, 1), point(snap, s2, 2)), before)

function cmp1(snap::SnapRounding, n, before, u, v)
	c = cmpf(snap, n, before)(u, v)
# 	println("\e[36mchecking ($n, $before)($u, $v) = $c\e[m")
	return c
end
function check(snap::SnapRounding, n, before)
	for i in 1:length(snap.active), j in 1:length(snap.active)
# 		println("\e[36m($n,$before) checking ($i, $j) = $(Tuple.(snap.segments[snap.active[[i,j]]]))\e[m")
		c = cmp1(snap, n, before, snap.active[[i,j]]...)
# 		c = cmpf(snap, n, before)(snap.active[i], snap.active[j])
# 		println("    \e[1mreturned $c\e[m")
		@assert c == cmp(i,j)
	end
end
function locate(snap::SnapRounding, n, s, before)
	(lo, hi) = (0, 1+length(snap.active)); while lo+1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = cmpf(snap, n, before)(s, snap.active[m])
# 		println("comparing @$n segments $(snap.segments[s]) and $m=$(snap.segments[snap.active[m]]) returned $c")
		c == 0 && return m
		c ≥ 0 ? (lo = m) : (hi = m)
	end
	return lo+1
end

function Base.deleteat!(snap::SnapRounding, n, i)
	deleteat!(snap.active, i)
	# intersect i-1 and i if needed
	2 ≤ i ≤ length(snap.active) && newinter(snap, n, i-1)
end
"""
deletes segment `s` *before* point `n`
"""
function Base.delete!(snap::SnapRounding, n, s)
	println("\e[31mdelete!\e[m($n, $s=$(snap.segments[s])) $(Tuple.(snap.segments[snap.active]))")
	i = locate(snap, n, s, true)
	println("  found i=$i / $(Tuple.(snap.segments[snap.active]))")
# 	@assert snap.active[i] == s
# 	deleteat!(snap.active, s)
	for (i, t) in pairs(snap.active)
		t ≠ s && continue
		deleteat!(snap, n, i)
		return true
	end
	return false
end
"""
inserts segment `s` in the correct place *after* point `n`
"""
function Base.insert!(snap::SnapRounding, n, s)
	println("\e[32minsert!\e[m($n, $s=$(snap.segments[s])) // $(snap.segments[snap.active])")
	@assert s ∉ snap.active
	i = locate(snap, n, s, false)
	insert!(snap.active, i, s)
	i > 1 && newinter(snap, n, i-1)
	i < length(snap.active) && newinter(snap, n, i)
# 	println("  after insertion of $(snap.segments[s]): $(snap.segments[snap.active])")
end

# Treating a pixel««1
"""
treats the n-th event in the structure
heats the pixel containing `point`, if needed,
clips outgoing segments
"""
function process(snap::SnapRounding, n)
# 	check(snap, n, true)
	p = point(snap, n)
	e = event(snap, n)
	ε = snap.ε
	unique!(sort!(e.pass))
	center = round(snap, p)
	if !isempty(e.pass)
		pix2 = pixel!(snap, p)
		for s in e.pass
			c = point(snap, s, 3)
			pix1 = pixel!(snap, c)
			pix1 ≠ pix2 && push!(snap.links, (pix1, pix2))
			println("\e[34;7;1m link pixels $pix1 => $pix2\e[m on segment $(snap.segments[s])")
# 			@assert pix2 < 18
# 			@assert (pix1, pix2) ≠ (17,9)
			found = delete!(snap, n, s)
			# if no segment was found
			found || continue
			snap.segments[s][2] == n && continue # n is the end point
			q = pixel_exit(center, ε, point(snap, s, 1), point(snap, s, 2))
			q == nothing && continue
			# mark for reinsertion at q:
			println("\e[35m mark $(snap.segments[s]) for reinsertion at point $q\e[m")
			c = point!(snap, q)
			snap.segments[s][3] = n
			push!(event(snap, c).reinsert, s)
# 		println("reinserting $(snap.segments[s]) at point $c=$q")
	# 		push!(snap.links, (snap.segments[s][3], n))
# 			println("\e[34;1m $(snap.segments[s]) =>,$n\e[m")
		end
	end
		

	for s in e.start
# 		println("at $n: segments $s=$(snap.segments[s]) starts")
		insert!(snap, n, s)
	end
	for s in e.reinsert
		insert!(snap, n, s)
	end
# 	check(snap, n, false)
end
# Main function««1
"""
Snap-rounding after the [Hershberger 1996] algorithm.
Input: geometric points, symbolic segments, precision ε.
Output:
 * new geometric points
 * list of point merges (if any)
 * all segments between the new points
"""
function snapround(points, segments, ε)#««
	println("\e[48;5;87m****************************\e[m"*"\n"^30)
	println(points)
	println(segments)
	io = open("/tmp/points", "w")#««
	println(io, "# plist=$points\n# elist=$segments")
	for (i,p) in pairs(points)
		println(io,"$(p[1])\t$(p[2]) # $i")
	end
	println(io,"\n\n")
	for (i1,i2) in segments
		(x1,y1) = points[i1]; (x2,y2) = points[i2]
		println(io,"$x1\t$y1\t$(x2-x1)\t$(y2-y1)\t# $i1--$i2")
	end
	close(io)#»»
	snap = SnapRounding(points, segments, ε)

# 	count = 0
	while !isempty(snap.queue)
# 		count += 1
# 		@assert count < 10
		n = heappop!(snap.queue, ord(snap)); p = point(snap, n)
		print("""
\e[1;7m event $n = $p \e[m => $(Tuple.(snap.segments[event(snap,n).start]))
  \e[1mactive=$(snap.active)=$(Tuple.(snap.segments[snap.active]))\e[m
  pass=$(event(snap,n).pass)=$(Tuple.(snap.segments[event(snap,n).pass]))
  reinsert=$(event(snap,n).reinsert) $(Tuple.(snap.segments[event(snap,n).reinsert]))
  queue=$(sort(snap.queue; order=ord(snap)))
  segment[9]=$(length(snap.segments) ≥ 9 ? snap.segments[9] : "")
""")
	process(snap, n)
# 	@assert n ≠ 15
	end
	pix = similar(points, length(snap.pixel))
	for (k, v) in snap.pixel; pix[v] = k; end
	return (pix, sort(collect(snap.links)))
end#»»

function snapround!(plist, elist, ε)
	(pixels, links) = snapround(plist, elist, ε)
	println(links)
	@assert length(links)≠25
	println("returned $(length(plist)) => $(length(pixels)) points
	and $(length(elist))=>$(length(links)) segments")
	for i in length(plist)+1:length(pixels)
		push!(plist, pixels[i])
	end
	empty!(elist); push!(elist, links...)
	return nothing
end

# ««1
end
