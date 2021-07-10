module SegmentIntersections
using DataStructures
using StaticArrays
using FastClosures

"""
Sorts points (in lex (x,y) order);
returns the pair (sort permutation, dict from point name to new index).
"""
function rename_points(points, names)
	# sort (lex) points
	σ = sortperm(points)
	σinv = inverse_permutation(σ)
	dict = Dict(n => σinv[i] for (i, n) in pairs(names))
	return (σ, dict)
end

function inverse_permutation(σ)#««
	σinv = similar(σ)
	for (i,j) in pairs(σ)
		σinv[j] = i
	end
	return σinv
end#»»
function sortperm_inv(vect)#««
	σ = sortperm(vect)
	return (vect[σ], σ, inverse_permutation(σ))
end#»»
@inline det2(a,b)=a[1]*b[2]-a[2]*b[1]
@inline det3(a,b,c)=det2(b-a,c-a)

"""
Returns > 0 if point p is above the line (p1p2),
0 if it on the line (up to ε),
and < 0 if it is below this line.
"""
function point_cmp(p1,p2,p)#««
	return sign(det3(p1,p2,p))
end#»»
"""
Returns the insertion point for `p` in the list of segments `front`.
"""
function segfind(points, segments, front, p)#««
	# compare with one segment:
	# front[i] = index of segment
	# segments[front[i]] = pair of point indices
	# points[SA[segments[front[i]]...]] = pair of geometric points
	
	# we go for the simpler version returning a range, i.e. doing this in
	# two passes [http://www.tbray.org/ongoing/When/200x/2003/03/22/Binary]
	(lo, hi) = (0, 1+length(front))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = point_cmp(points[SA[segments[front[m]]...]]..., p)
		if c > 0
			lo = m
		else
			hi = m
		end
	end
	a1 = lo
	(lo, hi) = (0, 1+length(front))
	while lo + 1 < hi
		m = Base.Sort.midpoint(lo, hi)
		c = point_cmp(points[SA[segments[front[m]]...]]..., p)
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
Inserts a new segment in its correct position in the `front` list.
Deletes intersection points of former neighbors
"""
function insert_segment!(front, queue, segment)
end


struct Event{P}
	point::P
	segments::Vector{Int}
end
"""
 - `points`: list of points (as vectors)
 - `segments`: list of pairs `(a, b)` of indices of linked points
 - `names`: list of name of each point (indices used by default)

Returns
 - the list of points (old points + new points (appended at the end))
 - the list of segments between all those points.
"""
function intersections(points, segments, ε = 0)
	# first thing is to sort points (lex)
	(points, ϖ, ϖinv) = sortperm_inv(points)
	segments = [minmax(ϖinv[i], ϖinv[j]) for (i,j) in segments]
	# remove duplicate segments:
	unique!(sort!(segments))
	(segments_rev, ρ, ρinv) = sortperm_inv(reverse.(segments))
	# segments *starting* at j are searchsorted(segments,j;by=first)
	# segments *ending* at j are ρinv[searchsorted(segments_rev,j;by=first)]

	println(points)
	println(segments)
	println(segments_rev)
# 	println("searching for all segments arriving at 3:")
# 	println(ρinv[searchsorted(segments_rev, 3; by=first)])

	# we need to compute
	#  - all points which belong (up to ε) to another segment,
	#  - all (strict) intersection points of two segments
	# An event is either:
	#  - (segment, segment) for a strict intersection,
	#  - (-point, segment) for an endpoint on a segment.
	Event = NTuple{2,Int}
	front = SortedSet{Int}() # segment indices
	queue = SortedDict{eltype(eltype(points)),Event}()
	newsegments = NTuple{2,Int}()

	sweep = beforestartsemitoken(queue)
	for (i, p) in pairs(points)
		# 1. process any events happening before `x`
		while status((queue, sweep)) ==1 && deref_key((queue, sweep)) < p[1]
			# this event must be an intersection point at x=t:
			(t, event) = deref((queue, sweep))
			sweep = advance((queue, sweep))
		end
		# 2. process point i
		position = segfind(points, segments, front, p)
		# we look first for all segments *starting* at i:
		# insert them in the front structure, sorted by ascending slope
		println("\e[1mProcess point $i=$p\e[m")
		sseg = collect(searchsorted(segments, i; by=first))
		sort!(sseg; lt=(j,k)->
			det3(p,points[segments[j][2]],points[segments[k][2]])>0)
		println("   new segments: $(segments[sseg])")


		eseg = ρinv[searchsorted(segments_rev, i; by=first)]
		sort!(eseg; lt=(j,k)->
			det3(p,points[segments[j][1]],points[segments[k][1]])<0)
		println("   old segments: $(segments[eseg])")
		println("insertion of point $p in $front:")
		println(segfind(points, segments, front, p))
	end
	# 3. process any events remaining in the queue
end



end
