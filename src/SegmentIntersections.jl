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
	d = det2(p1-p2, p3-p4)
	abs(d) ≤ ε && return nothing
	c = 1/d
	d12 = c*det2(p1, p2)
	d34 = c*det2(p3, p4)
	z = d12*(p3-p4) - d34*(p1-p2)
	(p1 ≤ z ≤ p2) && (p3 ≤ z ≤ p4) && return z
	return nothing
end
"""
Returns > 0 if point p is above the line (p1p2),
0 if it on the line (up to ε),
and < 0 if it is below this line.
"""
function point_cmp(p1,p2,p)#««
	println(det2(p1,p2,p))
	return sign(det2(p1,p2,p))
end#»»
"""
Returns the insertion point for `p` in the list of active segments `active`.
"""
function segfind(points, active, p)#««
	# compare with one segment:
	# points[SA[segments[front[i]]...]] = pair of geometric points
	
	# we go for the simpler version returning a range, i.e. doing this in
	# two passes [http://www.tbray.org/ongoing/When/200x/2003/03/22/Binary]
	p[1] == 8.0 && println("  *** segfind($p)")
	(lo, hi) = (0, 1+length(active))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		(a,b) = active[m]
		c = point_cmp(points[a], points[b], p)
		p[1] == 8 && println("point_cmp(points[$a],points[$b],p) = $c")
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
		p[1]==8 && println("point_cmp(points[$a],points[$b],p) = $c")
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
Inputs
 - `points`: list of points (as vectors)
 - `segments`: list of pairs `(a, b)` of indices of linked points
 - `ε`: diameter of pixels for snap rounding.

Returns
 - the list of points (old points + new points (appended at the end));
 - the list of segments between all those points.
"""
function intersections!(points, segments, ε = 0)
	# make all segments left-to-right
	for (i, s) in pairs(segments)
		points[s[1]] ≤ points[s[2]] && continue
		segments[i] = (s[2], s[1])
	end
	# sorted by increasing slope:
	unique!(sort!(segments))
	# segmentsto[i] = origins for segments to i, by increasing slope:

	byslope = let points=points;
		i ->(j,k)->det2(points[[i,j,k]]...) > 0
	end
	segmentsfrom = [ sort(last.(segments[searchsorted(segments,i;by=first)]);
		lt=byslope(i))
		for i in eachindex(points)]
	# FIXME: this should insert the segment at correct position *after*
	# point i, i.e. (11,14) is *after* (12,3) in the cross picture.
	addsegment = @closure (i,q)-> begin
		pos = searchsortedfirst(segmentsfrom[i], q;lt=byslope(i))
		insert!(segmentsfrom[i], pos, q)
	end
# 	end
# 	rsegments = reverse.(segments); sort!(rsegments)
# 	segmentsto = [ sort(last.(rsegments[searchsorted(rsegments,i;by=first)]);
# 		lt=(j,k)->det2(points[[i,j,k]]...) < 0)
# 		for i in eachindex(points)]

	# build a queue of all future events:
	point_lt = @closure (i,j) -> points[i] < points[j]
	point_ord = Base.Order.Lt(point_lt)
	queue = heapify(eachindex(points), point_ord)
	println("queue = $queue")

	# list of active segments, sorted by y-intercept
	# FIXME: use an appropriate data structure instead of a vector
	active = NTuple{2,Int}[]
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

	# hot pixels are used to merge segments as needed:
	# the `reverse` operation ensures that only the smallest index is stored
	# rounding to `Int` prevents the dreaded `-0.0`
	roundcoord = iszero(ε) ? identity : @closure x->round(Int,x/ε)*ε
	pixel_id = Dict(roundcoord.(points[i])=>i
		for i in reverse(eachindex(points)))
	hotpixel = @closure p -> get!(pixel_id, roundcoord.(p), length(points)+1)
	# renaming of points (if needed):
	newname = Pair{Int,Int}[]
	for (i,p) in pairs(points)
		j = hotpixel(p)
		j ≠ i && push!(newname, i=>j)
	end


	while !isempty(queue)
		p = heappop!(queue, point_ord)
		println("\e[1;7m p = $p \e[m queue: $queue")
		new_inter = @closure n -> begin#««
			(n < 1 || n > length(active)) && return
			(a,b) = active[n]
			for q in segmentsfrom[p]
				(a==p||a==q||b==p||b==q) && continue
				z = inter(points[a],points[b],points[p],points[q],ε)
				z == nothing && continue
				# h is the point index for the new point; either an existing one,
				# or a new one, created for the occasion:
				h = get!(pixel_id, roundcoord.(z), length(points)+1)
				# two cases here:
				# 1. this is a new hot pixel: create a point
				# (and push to segmentsfrom[] as appropriate)
				println("\e[32;7m * \e[mcreate point $h = $z for ($a,$b)∩($p,$q)")
				println("   round = $(roundcoord.(z))")
				h == 18 && println("   round(12) = $(roundcoord.(points[12]))")
				if h == length(points)+1
					z = roundcoord.(z)
					push!(points, z)
					if det2(points[[h,b,q]]...) > 0
						push!(segmentsfrom, [b,q])
					else
						push!(segmentsfrom, [q,b])
					end
					heappush!(queue, h, point_ord)
				# 2. this is a pre-existing pixel:
				# break segments (a,b) and (p,q) as (a,h)+(h,b) and (p,h)+(h,q)
				# i.e. insert h in segmentsfrom[a] and segmentsfrom[b]
				else
					h ∈ (a,b,p,q) && continue
					println("  inserting segments $h => $b and $q")
					println("     before: $(segmentsfrom[h])")
					addsegment(h, b)
					addsegment(h, q)
					println("     after: $(segmentsfrom[h])")
				end
				replace!(segmentsfrom[a], b=>h)
				replace!(segmentsfrom[p], q=>h)
				# replace n-th active segment (a,b) by (a,h)
				active[n] = (a,h)
				println("  segmentfrom[$a] = $(segmentsfrom[a])")
				println("  segmentfrom[$p] = $(segmentsfrom[p])")
				println("  segmentfrom[$h] = $(segmentsfrom[h])")
				println("  now queue = $queue")
			end
		end#»»
		println("  segments from $p = $(segmentsfrom[p])")
		println("  current active list=$active")
# 		println("  segments to $p = $(segmentsto[p])")
		# find insert position for this point in `active`
		println("  ### segfind($(points[p]))")
		range = segfind(points, active, points[p])
		println("  insert at position $range")
		# if some segment(s) match the given point position:
		# replace (a,b) by (p,b) and output segment (a,b)
			# compute intersection of segment first(range)-1 and
			# all segments starting at p
		first(range) > 1 && new_inter(first(range)-1)
		last(range) < length(points) && new_inter(last(range)+1)
		for i in range
			(a,b) = active[i]
			push!(segments, (a,p))
# 			b == p && continue
# 			active[i] = (p,b)
			println("  \e[34;7moutput: segment ($a,$p)\e[m")
		end
		print("  \e[1m〈\e[m\e[38;5;8m")
		join(stdout, active[begin:first(range)-1], ",")
		print("\e[m\e[1m|\e[m\e[31m")
		join(stdout, active[range],",")
		print("\e[m\e[1m→\e[m\e[32m")
		join(stdout,[(p,i) for i in segmentsfrom[p]], ",")
		print("\e[m\e[1m|\e[m\e[38;5;8m")
		join(stdout, active[last(range)+1:end], ",")
		println("\e[m\e[1m〉\e[m")
# 		@assert length(points) ≤ 7
# 		println("  \e[38;5;8m$(active[begin:first(range)-1])\e[m\e[1m|\e[m\e[32m$([(p,i) for i in segmentsfrom[p]])\e[m\e[1m|\e[m\e[38;5;8m$(active[last(range)+1:end])\e[m")
		# insert new active segments
		active = [active[begin:first(range)-1];
			[(p,i) for i in segmentsfrom[p]];
			active[last(range)+1:end]]
		for j in segmentsfrom[p]
		end
	end
end



end
