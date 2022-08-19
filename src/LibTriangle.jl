module LibTriangle
using Triangulate: Triangulate

# constrained_triangulation
# basic_triangulation
function triangulation(vertices::AbstractVector, names = [],
		edges::AbstractVector = [], boundary = [], holes = []; reverse = false)
	pointlist = Float64[v[i] for i in 1:2, v in vertices]
	segmentlist = Int32[e[i] for i in 1:2, e in edges]
	segmentmarkerlist = isempty(boundary) ? Int32[0 for e in edges] :
		Int32.(boundary)
	pointmarkerlist = Vector{Int32}(names)
	holelist = Float64[v[i] for i in 1:2, v in holes]
	graph = Triangulate.TriangulateIO(;
		pointlist, segmentlist, pointmarkerlist, segmentmarkerlist, holelist)
	(tri, vor) = Triangulate.triangulate("pQ", graph)
	reverse && for t in eachcol(tri.trianglelist)
		t[1], t[2] = t[2], t[1]
	end
	return ((t[1], t[2], t[3]) for t in eachcol(tri.trianglelist))
end
end

# # isdefined(Main,:Triangle) ||
# (include("../tri/src/Triangle.jl"); using .Triangle)
# L = LibTriangle
# T = Triangle
# N=T.Triangulate.NativeInterface
# p=[[10.0, 1.0], [1.0, 1.0], [1.0, 10.0], [0.0, 10.0], [0.0, 0.0], [10.0, 0.0]]
# s=NTuple{2,Int}[(1,2),(2,3),(3,4),(4,5),(5,6),(6,1)]
# 
# v=Float64[transpose.(p)...;]
# vm=Int[1,2,3,4,5,6]
# e=Int[1 2; 2 3; 3 4; 4 5; 5 6; 6 1]
# em=Bool[true for _ in 1:6]
# 
# t2=T.constrained_triangulation(v, vm, e, em)
# t1=L.triangulation(p,[], s, em).trianglelist
