# After [Melhorn, Näher 1994] _Implementation of a sweep line algorithm
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

@inline det2(a,b)=a[1]*b[2]-a[2]*b[1]
@inline det2(a,b,c)=det2(b-a,c-a)
"(p1p2) ∩ (p3p4), or `nothing` if no intersection"
function inter(p1, p2, p3, p4, ε)
	@assert p1 < p2
	# special case: (p1p2) is vertical, we *need* the intersection point
	# to be on that same vertical line.
	(v12, v13, v14) = (p2-p1, p3-p1, p4-p1)
	d2 = det2(v13, v14)
	d3 = det2(v12, v14)
	d4 = det2(v12, v13)
	sign(d3) == sign(d4) && return nothing
	d = d3 - d4
	d1 = d2 - d
	sign(d1) == sign(d2) && return nothing
	abs(d) ≤ ε && return nothing
	z= p1+d2/d*v12
# 	println("p1=$p1\np2=$p2\np3=$p3\np4=$p4\n \e[31mz=$z\e[m\n")
# 	(p1[1] ≤ z[1] ≤ p2[1]) && (p3[1] ≤ z[1] ≤ p4[1]) &&
# 	(p1[2] ≤ z[2] ≤ p2[2]) && (p3[2] ≤ z[2] ≤ p4[2]) &&
		return z
	return nothing
end
"""
Returns > 0 if point p is above the line (p1p2),
0 if it on the line (up to ε),
and < 0 if it is below this line.
"""
function point_cmp(p1,p2,p)#««
	return sign(det2(p1,p2,p))
end#»»
"""
Returns the insertion point for `p` in the list of active segments `active`.
"""
function findpoint(points, active, p)#««
	# compare with one segment:
	# points[SA[segments[front[i]]...]] = pair of geometric points
	
	# we go for the simpler version returning a range, i.e. doing this in
	# two passes [http://www.tbray.org/ongoing/When/200x/2003/03/22/Binary]
# 	println("  *** findpoint($p)")
	(lo, hi) = (0, 1+length(active))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		(a,b) = active[m]
# 		println("looking in $lo:$hi; m=$m, (a,b)=($a,$b)")
		c = point_cmp(points[a], points[b], p)
# 		println("point_cmp(points[$a],points[$b],p) = $c")
		if c > 0
			lo = m
		else
			hi = m
		end
	end
	a1 = lo
	(lo, hi) = (0, 1+length(active))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		(a,b) = active[m]
		c = point_cmp(points[a], points[b], p)
# 		p[1]==8 && println("point_cmp(points[$a],points[$b],p) = $c")
		if c ≥ 0
			lo = m
		else
			hi = m
		end
	end
	a2 = hi
	return (a1+1:a2-1)
end#»»

# SNAP-ROUNDING, after [Hershberger 2006] ««1
# L¹ pixels ««2
"rounds in L¹ norm to an *even* multiple of ε, e.g. (0,0), (0,2), (1,1)…"
function roundL1(point, ε)
	iszero(ε) && return point
	u = round(Int,(point[1]+point[2])/(2ε))
	v = round(Int,(point[1]-point[2])/(2ε))
	return SA[(u+v)*ε,(u-v)*ε]
end

"returns the point where (p,q) exits the diagonal pixel with width ε"
function diagonal_exit(pixel, ε, p, q)#««
	c = pixel .+ (ε, 0) # right corner
	(cp, cq) = (p-c, q-c)
	@assert (cp[1]+cp[2] < 0 || cp[1]-cp[2] < 0)
	# point `q` is inside the pixel:
	(cq[1]+cq[2] < 0) && (cq[1]-cq[2] < 0) && return nothing
	d2 = det2(cp, cq)
	if d2 ≥ 0 # exits through the bottom
		x = d2/(cq[1]-cp[1]+cp[2]-cq[2])
		return c .- (x, x)
	else
		x = d2/(cp[1]+cp[2]-cq[1]-cq[2])
		return c .+ (-x,x)
	end
end#»»

# L∞ pixels««2

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

# Events (x-structure) ««2

@enum EventType begin
	SegmentStart
	SegmentStop
	SegmentReinsert
	SegmentIntersect
end

"""
Encodes either a segment intersection, or another event,
according to the constants below:
"""
struct Event{P}
	type::EventType
	point::P
	dest1::P
	dest2::P
end

# function trim(pixel, segment, ε)#««
# 	println("\e[1mtrim($sweepline, $segment)\e[m")
# 	i = findfirst(==(segment), sweepline) # FIXME sorted search!
# 	println("   found at $i")
# 	@assert i ≠ nothing "segment $segment found in sweepline"
# 	deleteat!(sweepline, i)
# 	(a,b) = segment
# 	q = pixel .+ ε/2
# 	all(points[b] .< q) && return
# 	# compute the re-entering point `r`
# 	r = pixel_exit(pixel, ε, points[a], points[b])
# 	println("  exit point is $r")
# 	r == nothing && return nothing
# 	return Event(r, segment, SEGMENT_REINSERT)
# end#»»


# Active segments (y-structure) ««2

struct SnapRounding{P,I<:Integer,E} # P is point type
	points::Vector{P}
	segmentsfrom::Vector{Vector{I}}
	ε::E
	queue::Vector{I} # should be a tree instead
	active::Vector{NTuple{2,I}} # should be a tree instead
	pixelid::Dict{P,Int}
# 	nextpoint::Vector{I}
	@inline SnapRounding{P,I}(points, segmentsfrom, ε::E) where{P,I,E} =
		new{P,I,E}(points, segmentsfrom, ε, eachindex(points), [], Dict{P,Int}())
end

@inline ord(snap::SnapRounding) = Base.Order.By(@closure i->snap.points[i])
@inline roundL1(snap::SnapRounding, point) = roundL1(point, snap.ε)
byslope(points,i) = (j,k)->det2(points[[i,j,k]]...) > 0

function SnapRounding(points, segments, ε)#««
	println("\e[48;5;57m****************************\e[m"*"\n"^30)
# 	z = inter(points[10],points[8],points[2],points[3],ε)
# 	println(sort([10,8,2,3],by=i->points[i]))
# 	println(z)
# 	return z
	# make all segments left-to-right
	for (i, s) in pairs(segments)
		(points[s[1]] >  points[s[2]]) && (segments[i] = (s[2], s[1]))
	end
	unique!(sort!(segments))

	segmentsfrom = [ sort(last.(segments[searchsorted(segments,i;by=first)]);
		lt=byslope(points, i)) for i in eachindex(points)]
	println("\e[36m",collect(pairs(segmentsfrom)),"\e[m")
	snap = SnapRounding{eltype(points),Int}(points, segmentsfrom, ε)
	for i in reverse(eachindex(points))
		snap.pixelid[roundL1(points[i], ε)] = i
	end
	heapify!(snap.queue, ord(snap))
	return snap
end#»»

function identify(snap::SnapRounding, point)
	if point[1] == -8.0
		println(collect(pairs(snap.pixelid)))
		println(point)
		println(point==snap.points[15])
	end
	n = get!(snap.pixelid, roundL1(snap, point), length(snap.points)+1)
	println("n = $n")
	if n == length(snap.points)+1
	println("\e[36mcreate point $n=$point\e[m")
		push!(snap.points, point)
		push!(snap.segmentsfrom, [])
		heappush!(snap.queue, n, ord(snap))
	end
	return n
end
function newsegmentfrom(snap::SnapRounding, p, q)#««
	println("newsegmentfrom($p, $q)")
	@assert snap.points[p] < snap.points[q] "points[$p] < points[$q]"
	slist = snap.segmentsfrom[p]
	pos = searchsorted(slist, q; lt=byslope(snap.points, p))
	if !isempty(pos) # keep only shorter vector
		println("segment ($p, $q) is collinear to segments[$pos] = ($p, $(slist[pos]))")
		i = first(pos)
		r = slist[i]
		@assert r ≠ p
		q == r && return # do nothing if q already present
		snap.points[r] < snap.points[q] && ((q,r) = (r,q))
		println("  points aligned in order $p--$q--$r")
		println("  replacing $(slist[i]) by $q")
		slist[i] = q
		println("  and calling newsegmentfrom($q, $r)")
		q ≠ r && newsegmentfrom(snap, q, r)
	else
		insert!(slist, first(pos), q)
	end
end#»»
"returns the insertion point for new segment (i,j) in the list of active segments"
function locate(snap::SnapRounding, cmpf)#««
# 	p = snap.points[i]
	(lo, hi) = (0, 1+length(snap.active)); while lo+1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = cmpf(snap.points[[snap.active[m]...]])
		c > 0 ? (lo = m) : (hi = m)
	end
	a1 = lo
	(lo, hi) = (0, 1+length(snap.active)); while lo+1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = cmpf(snap.points[[snap.active[m]...]])
		c ≥ 0 ? (lo = m) : (hi = m)
	end
	a2 = hi
	return (a1+1:a2-1)
end#»»
@inline cmppoint(p) = s-> det2(s[1]-p, s[2]-p)
@inline cmpleft(p, q) = s -> begin
	(a,b) = s
	d = det2(a-p, b-p)
	!iszero(d) && return d
	# FIXME
end
"returns the insertion point for this event in the list of active segments"
function locate(list, p, q)#««
	cmpf = @closure s -> begin # returns 0 if point p is *above* s
		(a,b) = s
		d = det2(a-p, b-p)
		!iszero(d) && return d
		# FIXME
	end
	(lo, hi) = (0, 1+length(list)); while lo+1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = cmpf(list[m].segment)
		c > 0 ? (lo = m) : (hi = m)
	end
	a1 = lo
	(lo, hi) = (0, 1+length(list)); while lo+1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = cmpf(list[m].segment)
		c ≥ 0 ? (lo = m) : (hi = m)
	end
	a2 = hi
	return (a1+1:a2-1)
end#»»

"""
intersects active[j] and active[j+1] if needed.
"""
function newinter(snap::SnapRounding, j)#««
	println("newinter(snap, $j)\n$(snap.active)")
	(a1, b1) = snap.active[j]
	(a2, b2) = snap.active[j+1]
	println("\e[34;1m intersect active[$j]=($a1,$b1) and [$(j+1)]=($a2,$b2)\e[m")
	b1 == b2 && return false
	point = inter(snap.points[[a1,b1,a2,b2]]..., 0)
	point == nothing && return false
	@assert point ≥ snap.points[a1]
	@assert point ≥ snap.points[a2]
	n = identify(snap, point)
	n ∈ (b1,b2) && return true
	newsegmentfrom(snap, n, b1)
	newsegmentfrom(snap, n, b2)
	replace!(snap.segmentsfrom[a1], b1 => n)
	replace!(snap.segmentsfrom[a2], b2 => n)
	snap.active[j]   = (a1, n)
	snap.active[j+1] = (a2, n)
	println("  \e[34m point $n=$point ($(snap.segmentsfrom[n]))\e[m")
	println("  segmentsfrom[$a1]=$(snap.segmentsfrom[a1])")
	println("  segmentsfrom[$a2]=$(snap.segmentsfrom[a2])")
	return true
end#»»»»

"""
inserts new segments, computes intersection for top of range
and bottom of range.
no intersection is computed between the inserted segments themselves
(they are assumed divergent).
"""
function Base.splice!(snap::SnapRounding, range, segments)#««
	(bot, top) = (first(range), last(range))
	# show changes in active list««
	print("\e[1m〈\e[m\e[38;5;8m")
	join(stdout, snap.active[begin:bot-1], ",")
	print("\e[m\e[1m|\e[m\e[31m")
	join(stdout, snap.active[range],",")
	print("\e[m\e[1m→\e[m\e[32m")
	join(stdout,segments, ",")
	print("\e[m\e[1m|\e[m\e[38;5;8m")
	join(stdout, snap.active[top+1:end], ",")
	println("\e[m\e[1m〉\e[m")#»»
	# insert new snap.active segments
	splice!(snap.active, range, segments)
	println(" now active=$(snap.active) bot=$bot")
	if isempty(segments)
		# TODO: intersect bottom and top of range
		bot < length(segments) && newinter(snap, bot-1)
	else
		bot > 1 && newinter(snap, bot-1)
		top = bot + length(segments)-1
		top < length(snap.active) && newinter(snap, top)
	end
	println("  after inter: active=$(snap.active)")
end#»»

"""
Snap-rounding after the [Hershberger 1996] algorithm.
Input: geometric points, symbolic segments, precision ε.
Output:
 * new geometric points
 * list of point merges (if any)
 * all segments between the new points
"""
function snapround(points, segments, ε)#««
# 	println("\e[48;5;87m****************************\e[m"*"\n"^30)
	snap = SnapRounding(points, segments, ε)

	# pixel rounding:
	hotgraph = empty(segments)
	hotpixel = empty(points)
# 	roundcoord = iszero(ε) ? identity : @closure x->round(Int,x/ε)*ε
# 	pixelid = Dict(roundcoord.(points[i])=>i for i in reverse(eachindex(points)))
# 	hotpixel = @closure p -> get!(pixelid, roundcoord.(p), length(points)+1)

	count = 0
	while !isempty(snap.queue)
		count += 1
# 		@assert count < 6
		p = heappop!(snap.queue, ord(snap)); point = snap.points[p]
		print("""
\e[1;7m point $p = $point\e[m  $(sort(snap.queue; order=ord(snap)))
  \e[3m$(snap.active)\e[m
  sfrom=\e[32m$(snap.segmentsfrom[p])\e[m
""")
		# first remove all segments ending at this point
		range = locate(snap, cmppoint(point))
		(bot, top) = (first(range), last(range))
		println("\e[31m found range = $range\e[m")
		slist = snap.segmentsfrom[p]
		splice!(snap, range, [(p,i) for i in slist])
	end

# 	unique!(sort!(hotpixel))
	return hotpixel
end#»»

# ««1
end
