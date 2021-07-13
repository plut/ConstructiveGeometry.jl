# J. Hobby, Practical segment intersection with finite precision output
# https://www.sciencedirect.com/science/article/pii/S0925772199000218/pdf?isDTMRedir=true&download=true
# J. Hersberger, Stable snap rounding
# https://www.sciencedirect.com/science/article/pii/S0925772112001551
# J. Hershberger, Improved output-sensitive snap rounding
# https://link.springer.com/content/pdf/10.1007/s00454-007-9015-0.pdf
# P. Guigue, thèse
# https://tel.archives-ouvertes.fr/tel-00471447/document
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
	(lo, hi) = (0, 1+length(active))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		(a,b) = active[m]
		c = point_cmp(points[a], points[b], p)
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
 - `points`: list of points (as vectors)
 - `segments`: list of pairs `(a, b)` of indices of linked points
 - `names`: list of name of each point (indices used by default)

Returns
 - the list of points (old points + new points (appended at the end))
 - the list of segments between all those points.
"""
function intersections!(points, segments, ε = 0)
	# http://people.scs.carleton.ca/~michiel/lecturenotes/ALGGEOM/bentley-ottmann.pdf
	# make all segments left-to-right
	for (i, s) in pairs(segments)
		points[s[1]] ≤ points[s[2]] && continue
		segments[i] = (s[2], s[1])
	end
	# segmentsfrom[i] = destination points for segments from i,
	# sorted by increasing slope:
	unique!(sort!(segments))
	segmentsfrom = [ sort(last.(segments[searchsorted(segments,i;by=first)]);
		lt=(j,k)->det2(points[[i,j,k]]...) > 0)
		for i in eachindex(points)]
	# segmentsto[i] = origins for segments to i, by increasing slope:
	rsegments = reverse.(segments); sort!(rsegments)
	segmentsto = [ sort(last.(rsegments[searchsorted(rsegments,i;by=first)]);
		lt=(j,k)->det2(points[[i,j,k]]...) < 0)
		for i in eachindex(points)]

	# build a queue of all future events:
	point_lt = @closure (i,j) -> points[i] < points[j]
	point_ord = Base.Order.Lt(point_lt)
	queue = heapify(eachindex(points), point_ord)
	println("queue = $queue")

	# list of active segments, sorted by y-intercept
	# FIXME: use an appropriate data structure instead of a vector
	active = NTuple{2,Int}[]

	# we ue the incoming segments list to store the result:
	empty!(segments)

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
				# FIXME: search if z is an existing point
				# replace segments (a,b) by (a,z), (p,q) by (p,z)
				push!(points, z); n = length(points)
				replace!(segmentsfrom[a], b=>n)
				replace!(segmentsfrom[p], q=>n)
				if det2(points[[n,b,q]]...) > 0
					push!(segmentsfrom, [b,q])
				else
					push!(segmentsfrom, [q,b])
				end
				println("\e[32;7m * \e[ecreate point $n = $z for ($a,$b)∩($p,$q)")
				println("  segmentfrom[$a] = $(segmentsfrom[a])")
				println("  segmentfrom[$p] = $(segmentsfrom[p])")
				println("  segmentfrom[$n] = $(segmentsfrom[n])")
				heappush!(queue, n, point_ord)
				println("  now queue = $queue")
			end
		end#»»
		println("  segments from $p = $(segmentsfrom[p])")
		println("  current active list=$active")
# 		println("  segments to $p = $(segmentsto[p])")
		# find insert position for this point in `active`
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
