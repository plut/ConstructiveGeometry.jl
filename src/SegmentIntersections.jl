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
# Geometry tools««1
@inline det2(a,b)=a[1]*b[2]-a[2]*b[1]
@inline det2(a,b,c)=det2(b-a,c-a)
"(p1p2) ∩ (p3p4), or `nothing` if no intersection"
function inter(p1, p2, p3, p4, ε)#««
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
end#»»
"""
Returns > 0 if point p is above the line (p1p2),
0 if it on the line (up to ε),
and < 0 if it is below this line.
"""
@inline function point_cmp(p1,p2,p)#««
	return sign(det2(p1,p2,p))
end#»»

# L¹ pixels ««1
"rounds in L¹ norm to an *even* multiple of ε, e.g. (0,0), (0,2), (1,1)…"
function roundL1(point, ε)#««
	# we need diamond-pixels closed on the right side:
	# we use ceil() instead of round(), which has strange behaviour
	# for half-integers (rounding to closest *even* integer).
	iszero(ε) && return point
# 	println(((point[1]+point[2])/(2ε),(point[1]-point[2])/(2ε)))
	u = ceil(Int,(point[1]+point[2]-ε)/(2ε))
	v = ceil(Int,(point[1]-point[2]-ε)/(2ε))
	return SA[(u+v)*ε,(u-v)*ε]
end#»»
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
# L∞ pixels ««1
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
# SnapRounding structure««1
"""
 - `active`: all segments crossing sweepline, in order
 - `segmentsfrom`: all segments starting from this point
   (≈ events associated with this point in the queue)
 - `queue`: heap of future points to examine (sorted lex.)
"""
struct SnapRounding{P,I<:Integer,E} # P is point type
	points::Vector{P}
	segmentsfrom::Vector{Vector{I}}
	ε::E
	queue::Vector{I} # should be a tree instead
	active::Vector{NTuple{2,I}} # should be a tree instead
	isreinsert::Vector{Bool} # should be collected with other point properties as a StructVector
	pointid::Dict{P,I}
	@inline SnapRounding{P,I}(points, segmentsfrom, ε::E) where{P,I,E} =
		new{P,I,E}(points, segmentsfrom, ε, eachindex(points), [],
			[false for _ in points],
			Dict{P,I}(points[i]=>i for i in reverse(eachindex(points))),
			)
end

@inline ord(snap::SnapRounding) = Base.Order.By(@closure i->snap.points[i])
@inline Base.round(snap::SnapRounding, point) = roundL1(point, snap.ε)
@inline byslope(points,i) = (j,k)->det2(points[[i,j,k]]...) > 0
@inline DataStructures.heappush!(snap::SnapRounding, plist...) =
	for p in plist; p ∉ snap.queue && heappush!(snap.queue, p, ord(snap)); end

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
	heapify!(snap.queue, ord(snap))
	return snap
end#»»

function point!(snap::SnapRounding, point)
	@assert !isnan(point[1])
	n = get!(snap.pointid, point, length(snap.points)+1)
	if n == length(snap.points)+1
		println("\e[32;7;1m create point $n=$point\e[m")
		push!(snap.points, point)
		push!(snap.segmentsfrom, [])
		push!(snap.isreinsert, false)
	end
	return n
end
@inline points!(snap::SnapRounding, plist...)= [point!(snap, p) for p in plist]
# Heating a pixel««1
"""
heats the pixel containing `point`, if needed,
clips outgoing segments
"""
function pixel!(snap::SnapRounding, n)
	point = snap.points[n]
	ε = snap.ε
	if point[1] == -8.0
		println(point)
		println(point==snap.points[15])
	end
	center = round(snap, point)
	cright = center .+ (ε, 0)
	cbelow = center .- (0, ε)
	cabove = center .+ (0, ε)
	(nr, nb, na) = points!(snap, cright, cbelow, cabove)
# 	newsegmentfrom(snap, nb, nr)
# 	newsegmentfrom(snap, na, nr)
	if cright == point # end of hot pixel:
		println("\e[31;1m end of hot pixel at $cright\e[m")
		range = intercept(snap, snap.points[nb], snap.points[na])
		splice!(snap, range, [])
		return
	end
	heappush!(snap, nr)
	println("\e[36mpixel for $center: $nr=$cright $nb=$cbelow $na=$cabove\e[m")
	# the guards are added
	range = intercept(snap, snap.points[nb], snap.points[na])
	println("  \e[43msegments intercepted:\e[m $(snap.active[range])")
# 	range = min(j, first(range)):max(j+1,last(range))
# 	println("    after correction: $(snap.active[range])")
	newseg = empty(snap.active)
	# trim incoming segments:
	for i in range
		(a, b) = snap.active[i]
		b == nr && continue
		p = diagonal_exit(center, snap.ε, snap.points[[a,b]]...)
		p == nothing && continue
		n = point!(snap, p)
		snap.isreinsert[n] = true
		println("segment ($a,$b) exits at $n=$p")
		newsegmentfrom(snap, n, b)
		heappush!(snap, n)
		push!(newseg, (n, b))
		println("   segment $(snap.active[i]) exits at $n=$p")
	end
	# insert (truncated) segments from point in list
	# segments going down are inserted between nb and nr
	# segments going up are inserted between na and nr (reverse order)
	segdown = [nb]
	segup = [na]
	println("\e[35m insert truncated segments: $n -- $(snap.segmentsfrom[n])\e[m")
	for b in snap.segmentsfrom[n]
		c = point!(snap, diagonal_exit(center, snap.ε, snap.points[[n, b]]...))
		println("  segment ($n, $b) exits pixel at $c")
		y = snap.points[c][2] - center[2]
		snap.isreinsert[c] = true
		newsegmentfrom(snap, c, b)
		println("  segmentsfrom[$c] = $(snap.segmentsfrom[c])")
		(x, y) = snap.points[b] .- center
		y == 0 && continue
		heappush!(snap, c)
		x ≠ 0 && push!(y < 0 ? segdown : segup, c)
	end
	newsegs = sizehint!(NTuple{2,Int}[], length(segup)+length(segdown))
	push!(segdown, nr)
	reverse!(segup); push!(segup, nr)
	for s in (segdown, segup), i in 1:length(s)-1
		push!(newsegs, (s[i], s[i+1]))
	end
# 	newsegs = [
# 		(segdown[i], segdown[i+1]) for i in 1:length(segdown)-1
# 		(segup[i+1], segup[i]) for i in length(segup)-1:-1:1
# 	]
	println("\e[1mreplace segments by $newsegs\e[m ($segdown, $segup)")
	splice!(snap, range, newsegs)
end
function newsegmentfrom(snap::SnapRounding, p, q)#««
	println("newsegmentfrom($p, $q)")
	@assert snap.points[p] < snap.points[q] "points[$p] < points[$q]"
	plist = snap.segmentsfrom[p]
	pos = searchsorted(plist, q; lt=byslope(snap.points, p))
	if !isempty(pos) # keep only shorter vector
		println("segment ($p, $q) is collinear to segments[$pos] = ($p, $(plist[pos]))")
		i = first(pos)
		r = plist[i]
		@assert r ≠ p
		q == r && return # do nothing if q already present
		snap.points[r] < snap.points[q] && ((q,r) = (r,q))
		println("  points aligned in order $p--$q--$r")
		println("  replacing $(plist[i]) by $q")
		plist[i] = q
		println("  and calling newsegmentfrom($q, $r)")
		q ≠ r && newsegmentfrom(snap, q, r)
	else
		insert!(plist, first(pos), q)
	end
end#»»
"returns the range of active segments intercepted by points p and q"
function intercept(snap::SnapRounding, p, q)#««
	# two passes [http://www.tbray.org/ongoing/When/200x/2003/03/22/Binary]
# 	p = snap.points[i]
	(lo, hi) = (0, 1+length(snap.active)); while lo+1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = det2(p, snap.points[[snap.active[m]...]]...)
# 		c = cmpf(snap.points[[snap.active[m]...]])
		c > 0 ? (lo = m) : (hi = m)
	end
	a1 = lo
	(lo, hi) = (0, 1+length(snap.active)); while lo+1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = det2(q, snap.points[[snap.active[m]...]]...)
# 		c = cmpf(snap.points[[snap.active[m]...]])
		c ≥ 0 ? (lo = m) : (hi = m)
	end
	a2 = hi
	return (a1+1:a2-1)
end#»»
# @inline cmppoint(p) = s-> det2(s[1]-p, s[2]-p)
# @inline cmpleft(p, q) = s -> begin
# 	(a,b) = s
# 	d = det2(a-p, b-p)
# 	!iszero(d) && return d
# 	# FIXME
# end

# Adding an intersection point««1
# "returns the insertion point for this event in the list of active segments"
# function locate(list, p, q)#««
# 	cmpf = @closure s -> begin # returns >0 if point p is *above* s
# 		(a,b) = s
# 		d = det2(a-p, b-p)
# 		!iszero(d) && return d
# 		# FIXME
# 	end
# 	(lo, hi) = (0, 1+length(list)); while lo+1 < hi
# 		m = Base.Sort.midpoint(lo, hi)
# 		c = cmpf(list[m].segment)
# 		c > 0 ? (lo = m) : (hi = m)
# 	end
# 	a1 = lo
# 	(lo, hi) = (0, 1+length(list)); while lo+1 < hi
# 		m = Base.Sort.midpoint(lo, hi)
# 		c = cmpf(list[m].segment)
# 		c ≥ 0 ? (lo = m) : (hi = m)
# 	end
# 	a2 = hi
# 	return (a1+1:a2-1)
# end#»»
"""
intersects active[j] and active[j+1] if needed.
pushes the intersection point in the queue (with appropriate segmentsfrom[]).
"""
function newinter(snap::SnapRounding, j)#««
	(a1, b1) = snap.active[j]
	(a2, b2) = snap.active[j+1]
	a1 == a2 && return false
	b1 == b2 && return false
	println("newinter(snap, $j)\n$(snap.active)")
	println("\e[34;1m intersect active[$j]=($a1,$b1) and [$(j+1)]=($a2,$b2)\e[m")
	p = inter(snap.points[[a1,b1,a2,b2]]..., 0)
	p == nothing && return false
	n = point!(snap, p)
	println("   \e[34m ($a1,$b1) ∩ ($a2,$b2) found point $n\e[m")
	if n == a1
		snap.active[j+1] = (a2, a1)
		newsegmentfrom(snap, a1, b2)
	elseif n == a2
		println("\e[35;1;7m case a2=$a2:\e[m")
		snap.active[j] = (a1, a2)
		println("\e[35mchange active[$j] to ($a1,$a2)\e[m")
		newsegmentfrom(snap, a2, b1)
	elseif n == b1
		snap.active[j+1] = (a2, b1)
		newsegmentfrom(snap, b1, b2)
	elseif n == b2
		snap.active[j] = (a1, b2)
		newsegmentfrom(snap, b2, b1)
	else
		snap.active[j] = (a1, n)
		snap.active[j+1] = (a2, n)
		newsegmentfrom(snap, n, b1)
		newsegmentfrom(snap, n, b2)
		heappush!(snap, n)
	end
end#»»

# Inserting new segments ««1
"""
inserts new segments, computes intersection for top of range
and bottom of range.
no intersection is computed between the inserted segments themselves
(they are assumed divergent).
"""
function Base.splice!(snap::SnapRounding, range, segments)#««
	println("splice: $range => $segments")
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
	println(" now active=$(snap.active)")
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
# 	println("\e[48;5;87m****************************\e[m"*"\n"^30)
	io = open("/tmp/points", "w")#««
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

	# pixel rounding:
	hotgraph = empty(segments)
	hotpixel = empty(points)
# 	roundcoord = iszero(ε) ? identity : @closure x->round(Int,x/ε)*ε

	count = 0
	while !isempty(snap.queue)
		count += 1
# 		@assert count < 6
		n = heappop!(snap.queue, ord(snap)); p = snap.points[n]
		print("""
\e[1;7m p $n = $p $(snap.isreinsert[n])\e[m  $(sort(snap.queue; order=ord(snap)))
  active=\e[3m$(snap.active)\e[m
  segfrom=\e[32m$(snap.segmentsfrom[n])\e[m
""")
		if snap.isreinsert[n]
			println("\e[36;7;1mpoint $n=$p is reinsert: $(snap.segmentsfrom[n])\e[m")
			range = intercept(snap, p, p)
			println("\e[36m range = $range\e[m")
			splice!(snap, range, [(n, i) for i in snap.segmentsfrom[n]])
		else
			pixel!(snap, n)
		end
		# first remove all segments ending at this p
# 		range = intercept(snap, p, p)
# 		(bot, top) = (first(range), last(range))
# 		println("\e[31m found range = $range\e[m")
# 		slist = snap.segmentsfrom[n]
# 		splice!(snap, range, [(n,i) for i in slist])
	end

# 	unique!(sort!(hotpixel))
	return hotpixel
end#»»

# ««1
end
