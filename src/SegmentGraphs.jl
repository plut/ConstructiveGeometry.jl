"""
    SegmentsGraphs

Computes self-intersection and dissection for a set of segments.

Input data type: (points, segments), where points is a 2d vector type
and segments is a list of pairs of point indices.

Interface: simplify!(points, segments, ε).
"""
module SegmentGraphs
using LinearAlgebra
using StaticArrays
include("SpatialSorting.jl")

@inline boundingbox(v::V...) where{V<:AbstractVector} =
	SpatialSorting.Box{V}(min.(v...), max.(v...))
@inline det2(u,v) = u[1]*v[2]-u[2]*v[1]

@inline segboxes(p,l) = [ boundingbox(p[s[1]], p[s[2]]) for s in l ]
@inline ptboxes(l, ε) =  [ SpatialSorting.Box(p, p .+ ε) for p in l ]

"computes all intersection points in the set of segments"
function intersections(points, segments, ε = 0;
		sboxes = segboxes(points, segments), pboxes = ptboxes(points, ε))
	r = similar(points, 0)
	for (s1, s2) in SpatialSorting.intersections(sboxes)
		p1 = points[segments[s1][1]]
		p2 = points[segments[s1][2]]
		p3 = points[segments[s2][1]]
		p4 = points[segments[s2][2]]
		d = det2(p1-p2, p3-p4)
		abs(d) ≤ ε && continue # collinear segments: no new point!
		c = 1/d
		d12 = c*det2(p1, p2)
		d34 = c*det2(p3, p4)
		z = d12*(p3-p4) - d34*(p1-p2)
		z ∈ sboxes[s1] || continue
		z ∈ sboxes[s2] || continue
		push!(r, z)
	end
	return r
end

"returns adjacency matrix, as a list of pairs of indices (segment, point)."
function seg_point_adjacency(points, segments, ε = 0;
		pboxes = ptboxes(points, ε), sboxes = segboxes(points, segments))
	boxes = [ sboxes; pboxes]
	ns = length(segments)
	r = NTuple{2,Int}[]

	for (i1, i2) in SpatialSorting.intersections(boxes)
		(i1, i2) = minmax(i1, i2)
		# only segment-point intersections:
		(i2 ≤ ns || i1 > ns) && continue
		is = i1; ip = i2-ns
		p1 = points[segments[is][1]]
		p2 = points[segments[is][2]]
		p  = points[ip]
		d = det2(p2-p1, p-p1)
		norm(p-p1, Inf) ≤ ε && continue
		norm(p-p2, Inf) ≤ ε && continue
		abs(d) > ε && continue
		push!(r, (is, ip))
	end
	return sort!(r)
end

"""
    simplify!(points, segments, ε)

Splits all segments at intersection points and points in segments.
Data is returned by modifying the input datasets
(`points` is modified only by appending new points).
"""
function simplify!(points, segments, ε = 0)
	sboxes = segboxes(points, segments)
	pboxes = ptboxes(points, ε)
	newpoints = intersections(points, segments, ε; pboxes, sboxes)
	newboxes = [SpatialSorting.Box(p, p .+ ε) for p in newpoints]
	allboxes = [pboxes; newboxes]

	redundant_points =
		filter(>(length(points)), maximum.(SpatialSorting.intersections(allboxes)))
	if !isempty(redundant_points)
		unique!(sort!(redundant_points))
		# merge non-redundant points:
		n=length(points); j = 1; rp = redundant_points[j]-n
		for i in eachindex(newpoints)
			while rp < i && j < length(redundant_points)
				j+= 1
				rp = redundant_points[j]-n
			end
			rp == i && continue
			push!(points, newpoints[i])
			push!(pboxes, newboxes[i])
		end
	end
	adjacency = seg_point_adjacency(points, segments, ε; pboxes, sboxes)
	# group by segment:
	start = 1
	while start ≤ length(adjacency)
		stop = start
		s = adjacency[start][1]
		while stop≤ length(adjacency) && adjacency[stop][1] == s; stop+=1; end
		# find best coordinate and sort
		direction = points[segments[s][2]]-points[segments[s][1]]
		a1 = abs(direction[1]); a2 = abs(direction[2])
		plist = last.(adjacency[start:stop-1])
		if a1 > a2
			ins = (direction[1] > 0) ?
				[ points[i][1] for i in plist] : [-points[i][1] for i in plist]
		else
			ins = (direction[2] > 0) ?
			  [ points[i][2] for i in plist] : [-points[i][2] for i in plist]
		end
		perm = sortperm(ins)
		seg = segments[s]
		u = seg[1]
		v = adjacency[start-1+perm[1]][2]
		segments[s] = minmax(u,v)
		for i in perm[2:end]
			u = adjacency[start-1+i][2]
			push!(segments, minmax(v, u))
			v = u
		end
		push!(segments, minmax(v, seg[2]))
		start = stop
	end
	unique!(sort!(segments))
end

end # module
