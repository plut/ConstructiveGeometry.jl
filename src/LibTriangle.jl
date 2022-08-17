module LibTriangle
using Triangulate: Triangulate

# constrained_triangulation
# basic_triangulation
function triangulation(vertices::AbstractVector, names = [],
		edges::AbstractVector = [], holes = []; reverse = false)
	pointlist = Float64[v[i] for i in 1:2, v in vertices]
	segmentlist = Int32[e[i] for i in 1:2, e in edges]
	pointmarkerlist = Vector{Int32}(names)
	holelist = Float64[v[i] for i in 1:2, v in holes]
	graph = Triangulate.TriangulateIO(;
		pointlist, segmentlist, pointmarkerlist, holelist)
	(tri, vor) = Triangulate.triangulate("-Q", graph)
	reverse && for t in eachcol(tri.trianglelist)
		t[1], t[2] = t[2], t[1]
	end
	return ((t[1], t[2], t[3]) for t in eachcol(tri.trianglelist))
end
end

using .LibTriangle
t = LibTriangle.triangulation([[0,0],[1,0],[0,1],[2,2]])
