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
	z= p1+det2(v13,v14)/d*v12
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
"compares segments `s1` and `s2` right of starting point of `s1`"
function seg_cmp_a(points, (a1,b1), (a2,b2))
	c = sign(det2(points[[a1,a2,b2]]...))
	c ≠ 0 && return c
	return sign(det2(points[[a1,b2,b1]]...))
end
"compares segments `s1` and `s2` left of ending point of `s1`"
function seg_cmp_b(points, (a1,b1), (a2,b2))
	c = sign(det2(points[[b1,a2,b2]]...))
	c ≠ 0 && return c
	return sign(det2(points[[b1,a1,a2]]...))
end
"""
Returns the intersection point for `seg` in the `list`
given the comparison function `cmpf`.
"""
function findsegment(points, list, seg, cmpf)
	(lo, hi) = (0, 1+length(list))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		if cmpf(points, seg, list[m]) > 0; lo = m; else; hi = m; end
	end
	a1 = lo
	(lo, hi) = (0, 1+length(list))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		if cmpf(points, seg, list[m]) ≥ 0; lo = m; else; hi = m; end
	end
	a2 = hi
	return (a1+1:a2-1)
end#»»
"""
Returns the intersection point for `seg` in the `list`
at left of ending point of `seg`.
"""
function findsegment_left(points, list, seg)
	(u,v) = seg
	(lo, hi) = (0, 1+length(list))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		(a,b) = list[m]
		c = point_cmp(points[a], points[b], points[v])
		(c == 0) && (c = point_cmp(points[u], points[b], points[v]))
		if c > 0; lo = m; else; hi = m; end
	end
	a1 = lo
	(lo, hi) = (0, 1+length(list))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		(a,b) = list[m]
		c = point_cmp(points[a], points[b], points[u])
		(c == 0) && (c = point_cmp(points[u], points[b], points[v]))
		if c ≥ 0; lo = m; else; hi = m; end
	end
	a2 = hi
	return (a1+1:a2-1)
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


"""
Encodes either s1 ∩ s2, start of s1 (as `(s1,0)`),
or end of s1 (as `(s1,-1)`).
"""
struct Event{P,I<:Integer}
	point::P
	seg1::NTuple{2,I}
	seg2::NTuple{2,I} # seg2 == 0 for start point, -1 for end point
end
@inline function Base.isless(e1::Event, e2::Event)
	e1.point < e2.point && return true
	e1.point > e2.point && return false
	# (re-insertions) < (end points) < (start points):
	return e1.seg2[1] > e2.seg2[1]
end

# function trim(points, sweepline, pixel, segment, ε)#««
# 	println("\e[33m trim($sweepline, $segment)\e[m")
# 	i = findfirst(==(segment), sweepline) # FIXME sorted search!
# 	println("   found at $i")
# 	deleteat!(sweepline, i)
# 	(a,b) = segment
# 	q = pixel .+ ε/2
# 	all(points[b] .< q) && return
# 	# compute the re-entering point `r`
# 	r = inter(points[a], points[b], pixel .+ (0,ε/2), q, 0)
# 	if r == nothing
# 		r = inter(points[a], points[b], pixel .+ (ε/2,0), q, 0)
# 		r == nothing && return
# 	end
# 	return Event(r, segment, (-2,0))
# end#»»
# function simplify_collinear!(slist, p, points)#««
# 	isempty(slist) && return
# 	j = 1
# 	for i in 2:length(slist)
# 		d = det2(points[p], points[slist[j]], points[slist[i]])
# 		if d > 0 # these are different points
# 			j+= 1
# 			slist[j] = slist[i]
# 		elseif points[j] < points[i] # collinear points, shorter vector
# 			slist[j] = slist[i]
# 		else # collinear points, longer vector: do nothing
# 		end
# 	end
# 	resize!(slist, j)
# 	slist
# end#»»

"""
creates (if it exists) a new intersection between (p,i) and active[j].
returns `true` iff an intersection was found.
"""
function newinter(points, active, segmentsfrom, p, i, j; above=true)#««
	(a,b) = active[j]
	println("\e[34;1m intersect ($p, $i) and active[$j]=($a,$b)\e[m")
	b == i && return false
	q = inter(points[[p,i,a,b]]..., ε)
	(q == nothing || q < points[p] || q < points[a]) && return false
	n = get!(pixel_id, q, length(points)+1)
	if n == length(points)+1
		push!(points, q)
		push!(segmentsfrom, above ? [b2,b1] : [b1,b2])
		heappush!(queue, n, point_ord)
	end
	replace!(segmentsfrom[a], b => n)
	replace!(segmentsfrom[p], i => n)
	active[j] = (a,n)
	println("  \e[34m point $n=$q ($(segmentsfrom[n]))\e[m")
	println("  segmentsfrom[$a1]=$(segmentsfrom[a1])")
	println("  segmentsfrom[$a2]=$(segmentsfrom[a2])")
	return true
end#»»
"""
Inputs
 - `points`: list of points (as vectors)
 - `segments`: list of pairs `(a, b)` of indices of linked points
 - `ε`: diameter of pixels for snap rounding.

Returns
 - the list of points (old points + new points (appended at the end));
 - the list of segments between all those points.
"""
function intersections!(points, segments, ε = 0)
	roundcoord = iszero(ε) ? identity : @closure x->round(Int,x/ε)*ε
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

	# snap-rounding [Hershberger 2006]:
# 	queue = sizehint!(Event{eltype(points),Int}[],2*length(segments))#««
# 	for (a,b) in (segments)
# 		push!(queue, Event(points[a], (a,b), (0,0)))
# 		push!(queue, Event(points[b], (a,b), (-1,0)))
# 	end
# 	heapify!(queue)
# 	hotpixel = similar(points, 0)
# 	sweepline = similar(segments, 0)
# 
# 
# 	while !isempty(queue)
# 		event = heappop!(queue)
# 		println("\e[1;7m event at $(event.point): $(event.seg1) $(event.seg2)\e[m")
# 		println("  sweepline = $sweepline")
# 		push!(hotpixel, roundcoord.(event.point))
# 		if event.seg2[1] ≠ -2
# 			println("NOT A RE-INSERTION, hot pixel")
# 		end
# 		if event.seg2[1] == -1 # end point
# 			seg = event.seg1
# 			println("  \e[31;1m at point $(seg[2])=$(event.point)\e[m")
# 			println("  remove $seg from sweepline")
# 			range = findsegment(points, sweepline, seg, seg_cmp_b)
# 			deleteat!(sweepline, [ i for i in range if sweepline[i] == seg ])
# 		elseif event.seg2[1] == 0 # start of segment
# 			seg = event.seg1
# 			println("  insert $seg in sweepline")
# 			range = findsegment(points, sweepline, seg, seg_cmp_a)
# 			println("    range=$range")
# 			if first(range) < length(sweepline)
# 				seg2 = sweepline[first(range)]
# 				if seg2[1] ≠ seg[1]
# 					z = inter(points[[seg...]]..., points[[seg2...]]..., ε)
# 					z ≠ nothing && heappush!(queue, Event(z, seg, seg2))
# 				end
# 			end
# 			if last(range) ≥ 1
# 				seg1 = sweepline[last(range)]
# 				if seg1[1] ≠ seg[1]
# 					z = inter(points[[seg...]]..., points[[seg1...]]..., ε)
# 					z ≠ nothing && heappush!(queue, Event(z, seg1, seg))
# 				end
# 			end
# 			insert!(sweepline, first(range), seg)
# 			println("  now sweepline = $sweepline")
# 		else # s1 ∩ s2
# 			pixel = roundcoord.(event.point)
# 			println("\e[35;1m at intersection: $(event.seg1) ∩ $(event.seg2) = $(event.point)\e[m")
# 			println("  point=$(event.point)")
# 			println("  seg1: $(points[[event.seg1...]])")
# 			println("  seg2: $(points[[event.seg2...]])")
# 			ev1 = trim(points, sweepline, pixel, event.seg1, ε)
# 			ev1 ≠ nothing && heappush!(queue, ev1)
# 			ev2 = trim(points, sweepline, pixel, event.seg2, ε)
# 			ev2 ≠ nothing && heappush!(queue, ev2)
# 		end
# 	end
# 
# 
# 	unique!(sort!(hotpixel))
# 	return hotpixel#»»
	# “exact” Bentley-Ottmann
	# After [Melhorn, Näher 1994]

	byslope = let points=points;
		i ->(j,k)->det2(points[[i,j,k]]...) > 0
	end
	segmentsfrom = [ sort(last.(segments[searchsorted(segments,i;by=first)]);
		lt=byslope(i)) for i in eachindex(points)]
	println("\e[36m",collect(pairs(segmentsfrom)),"\e[m")
	newsegmentfrom = @closure (p,q) -> begin
		println("  calling newsegmentfrom($p, $q); list=$(segmentsfrom[p])")
		@assert points[p] < points[q] "points[$p] < points[$q]"
		slist = segmentsfrom[p]
		pos = searchsorted(slist, q; lt=byslope(p))
		if !isempty(pos) # keep only shorter vector
			println("segment ($p, $q) is collinear to segments[$pos] = ($p, $(slist[pos]))")
			i = first(pos)
			r = slist[i]
			@assert r ≠ p
			q == r && return # do nothing if q already present
			points[r] < points[q] && ((q,r) = (r,q))
			println("  points aligned in order $p--$q--$r")
			println("  replacing $(slist[i]) by $q")
			slist[i] = q
			println("  and calling newsegmentfrom($q, $r)")
			q ≠ r && newsegmentfrom(q, r)
		else
			insert!(slist, first(pos), q)
		end
	end
	# FIXME: this should insert the segment at correct position *after*
	# point i, i.e. (11,14) is *after* (12,3) in the cross picture.
# 	addsegment = @closure (i,q)-> begin
# 		pos = searchsortedfirst(segmentsfrom[i], q;lt=byslope(i))
# 		insert!(segmentsfrom[i], pos, q)
# 	end
# 	end
# 	rsegments = reverse.(segments); sort!(rsegments)
# 	segmentsto = [ sort(last.(rsegments[searchsorted(rsegments,i;by=first)]);
# 		lt=(j,k)->det2(points[[i,j,k]]...) < 0)
# 		for i in eachindex(points)]

	# build a queue of all future events:
	point_ord = Base.Order.By(@closure i -> points[i])
	queue = heapify(eachindex(points), point_ord)
	println("queue = $queue")

	# list of segments crossing sweepline,
	# represented as (i,j) where
	# `i` is the last point before the sweepline
	# `j` is the first point after the sweepline
	# FIXME: replace vector by a searchable list (skip list/balanced tree)
	active = NTuple{2,Int}[]
	# for all active segments, store either the intersection with the
	# previous active segment, or 0 if no such intersection exists:
	nextinter = Int[]

	# debug points to /tmp/points««
	io = open("/tmp/points", "w")
	for (i,p) in pairs(points)
		println(io,"$(p[1])\t$(p[2]) # $i")
	end
	println(io,"\n\n")
	for (i1,i2) in segments
		(x1,y1) = points[i1]; (x2,y2) = points[i2]
		println(io,"$x1\t$y1\t$(x2-x1)\t$(y2-y1)\t# $i1--$i2")
	end
	close(io)
#»»

	# we ue the incoming segments list to store the result:
	empty!(segments)

# 	# hot pixels are used to merge segments as needed:
# 	# the `reverse` operation ensures that only the smallest index is stored
# 	# rounding to `Int` prevents the dreaded `-0.0`
	pixel_id = Dict(roundcoord.(points[i])=>i for i in reverse(eachindex(points)))
# 	hotpixel = @closure p -> get!(pixel_id, roundcoord.(p), length(points)+1)
# 	# renaming of points (if needed):
# 	newname = Pair{Int,Int}[]
# 	for (i,p) in pairs(points)
# 		j = hotpixel(p)
# 		j ≠ i && push!(newname, i=>j)
# 	end


	while !isempty(queue)
		p = heappop!(queue, point_ord)
		# STEP 1. output vertex p (nothing to do)
		println("\e[1;7m p = $p $(points[p]) \e[m\n  queue: $(sort(queue;order=point_ord))")
		println("  segments from $p = $(segmentsfrom[p])")
		println("  current sweepline list=$active")
# 		new_inter = @closure n -> begin#««
# 			(n < 1 || n > length(active)) && return
# 			(a,b) = active[n]
# 			for q in segmentsfrom[p]
# 				println("    \e[3mtrying intersection ($p,$q) ∩ sweepline ($a,$b)…\e[m")
# 				(a==p||a==q||b==p||b==q) && continue
# 				z = inter(points[p],points[q],points[a],points[b],ε)
# 				z == nothing && continue
# 				# h is the point index for the new point; either an existing one,
# 				# or a new one, created for the occasion:
# 				h = get!(pixel_id, roundcoord.(z), length(points)+1)
# 				# two cases here:
# 				# 1. this is a new hot pixel: create a point
# 				# (and push to segmentsfrom[] as appropriate)
# 				println("\e[32;7m * \e[mcreate point $h = $z for ($a,$b)∩($p,$q)")
# 				println("   round = $(roundcoord.(z))")
# 				if h == length(points)+1
# 					z = roundcoord.(z)
# 					push!(points, z)
# 					if det2(points[[h,b,q]]...) > 0
# 						push!(segmentsfrom, [b,q])
# 					else
# 						push!(segmentsfrom, [q,b])
# 					end
# 					heappush!(queue, h, point_ord)
# 					if h == 19
# 						println("\e[35m")
# 						println(sort([10,19,15,16]; order=point_ord))
# 						println(points[[10,19,15,16]])
# 						println("\e[m")
# 					end
# 					println("  pushing $h on queue: $(sort(queue;order=point_ord))")
# 				# 2. this is a pre-existing pixel:
# 				# break segments (a,b) and (p,q) as (a,h)+(h,b) and (p,h)+(h,q)
# 				# i.e. insert h in segmentsfrom[a] and segmentsfrom[b]
# 				else
# 					h ∈ (a,b,p,q) && continue
# 					println("  inserting segments $h => $b and $q")
# 					println("     before: $(segmentsfrom[h])")
# 					addsegment(h, b)
# 					addsegment(h, q)
# 					println("     after: $(segmentsfrom[h])")
# 				end
# 				replace!(segmentsfrom[a], b=>h)
# 				replace!(segmentsfrom[p], q=>h)
# 				# replace n-th sweepline segment (a,b) by (a,h)
# 				active[n] = (a,h)
# 				println("  segmentfrom[$a] = $(segmentsfrom[a])")
# 				println("  segmentfrom[$p] = $(segmentsfrom[p])")
# 				println("  segmentfrom[$h] = $(segmentsfrom[h])")
# 				println("  now queue = $(sort(queue;order=point_ord))")
# 			end
# 		end#»»
		# STEP 2. determine range of segments meeting p
		range = findpoint(points, active, points[p])
		println("  segments meeting $p: $range")
		# STEP 3. output segments ending at p, and
		# STEP 4. mark segments not ending at p for re-insertion.
		slist = segmentsfrom[p]
		println("    at start: segmentsfrom[$p]=$slist")
		for (a,b) in active[range]
			println("\e[32;1;7m($a, $p)\e[m (-- $b)")
			println(sort([a,b,p]; order=point_ord))
			push!(segments, (a,p))
			# STEP 5. reinsert segments starting from p, and
			# STEP 6. sort in order, removing colinear segments.
			b ≠ p && newsegmentfrom(p, b)
		end
		println("    after insertions: segmentsfrom[$p]=$slist")

		# STEP 7. update Y-structure and X-structure (compute new intersections)
		(bot, top) = (first(range), last(range))

		# computes the intersection point for two (newly) consecutive segments
		newinter = @closure (s, j, above) -> begin #««
			(a,b) = active[j]
			println("\e[34;1m intersect ($p, $s) and active[$j]=($a,$b) ($above)\e[m")
			b == s && return true # nothing to do
			q = inter(points[[p,s,a,b]]..., ε)
			(q == nothing || q < points[p] || q < points[a]) && return false
			n = get!(pixel_id, q, length(points)+1)
			if n == length(points)+1
				push!(points, q)
				push!(segmentsfrom, Int[])
				heappush!(queue, n, point_ord)
			end
			if n == p
				# line (a,b) was supposed to not pass through p,
				# but we found (due to inexactitude) (a,b) ∩ (p,s) = p.
				# * replace (a,b) by (a,p) + (p,b)
				# * leave (p,s) alone
				println("\e[35m inexact case ($a,$b) ∩ ($p, $s)\e[m")
				println("should replace active segment ($a, $b) by ($a, $p)")
				newsegmentfrom(p, b)
				replace!(segmentsfrom[a], b => p)
				return true
			end
			n ≠ b && newsegmentfrom(n, b)
			n ≠ s && newsegmentfrom(n, s)
			n ≠ a && replace!(segmentsfrom[a], b => n)
			n ≠ p && replace!(segmentsfrom[p], s => n)
			active[j] = (a,n)
			println("  \e[34m point $n=$q ($(segmentsfrom[n]))\e[m")
			println("  segmentsfrom[$a]=$(segmentsfrom[a])")
			println("  segmentsfrom[$p]=$(segmentsfrom[p])")
			return true
		end#»»
		# (p, i) and active[j]
		if isempty(slist)
		else
			# find intersection(s) above:
			if top < length(active)
				for s in reverse(slist)
					newinter(s, top+1, true) || break
				end
			end
			if bot > 1
				for s in slist
					newinter(s, bot-1, false) || break
				end
			end
		end
		print("  \e[1m〈\e[m\e[38;5;8m")
		join(stdout, active[begin:first(range)-1], ",")
		print("\e[m\e[1m|\e[m\e[31m")
		join(stdout, active[range],",")
		print("\e[m\e[1m→\e[m\e[32m")
		join(stdout,[(p,i) for i in slist], ",")
		print("\e[m\e[1m|\e[m\e[38;5;8m")
		join(stdout, active[last(range)+1:end], ",")
		println("\e[m\e[1m〉\e[m")
		# insert new active segments
		active = [active[begin:bot-1]; [(p,i) for i in slist]; active[top+1:end]]
		p == 11 && return
	end
end



end
