using ConstructiveGeometry
using Test
using StaticArrays
using LinearAlgebra

G = ConstructiveGeometry

@testset "2d" begin#««
nva(s::G.Shapes.PolygonXor) = (length.(s.paths), G.Shapes.area(s))
r1 = G.Square(2.,1)
r2 = G.Square(1.,2)
r3 = G.Square(3.,3)
r4 = G.AffineTransform(G.TranslationMap([4,4]), G.Square(1,1))
@test nva(mesh(union(r1,r2))) == ([6], 3)
@test nva(mesh(union(r1,r3))) == ([4], 9)
@test nva(mesh(union(r1,r4))) == ([4,4], 3)
@test nva(mesh(intersect(r1,r2))) == ([4], 1)
@test nva(mesh(intersect(r2,r1))) == ([4], 1)
@test nva(mesh(intersect(r1,r3))) == ([4], 2)
@test nva(mesh(intersect(r3,r1))) == ([4], 2)
@test nva(mesh(intersect(r1,r4))) == ([], 0)
@test nva(mesh(setdiff(r1,r2))) == ([4], 1)
@test nva(mesh(setdiff(r3,r1))) == ([6], 7)
end#»»
@testset "Spatial sorting" begin#««
include("../src/SpatialSorting.jl")
# Bounding box««
struct BoundingBox{N,T}
	min::SVector{N,T}
	max::SVector{N,T}
end

BoundingBox{N}(min::AbstractVector, max::AbstractVector) where{N} =
	BoundingBox{N,promote_type(eltype(min), eltype(max))}(min, max)
BoundingBox(min::StaticVector{N}, max::StaticVector{N}) where{N} =
	BoundingBox{N}(min, max)
function BoundingBox(min::AbstractVector, max::AbstractVector)
	length(min) == length(max) || throw(DimensionMismatch(
	"min and max must have same length ($(length(min)),$(length(max)))"))
	return BoundingBox{length(min)}(min, max)
end

@inline Base.merge(box1::BoundingBox{N}, box2::BoundingBox{N}) where{N} =
	BoundingBox{N}(min.(box1.min, box2.min), max.(box1.max, box2.max))
@inline Base.:∩(box1::BoundingBox{N}, box2::BoundingBox{N}) where{N} =
	BoundingBox{N}(max.(box1.min, box2.min), min.(box1.max, box2.max))
@inline Base.isempty(box::BoundingBox) = any(box.min .> box.max)
@inline SpatialSorting.position(box::BoundingBox) = box.min + box.max

function rbb(int1, int2, n)
	v = rand(int1, n)
	w = rand(int2, n)
	return BoundingBox(v, v+w)
end
boxes=[ rbb(1:500, 1:10, 3) for _ in 1:3000 ]

function validate(boxes)
	int = extrema.(SpatialSorting.intersections(boxes))
	for i = 1:length(boxes), j = 1:i-1
		a = SpatialSorting.intersects(boxes[i], boxes[j])
		b = (j,i) ∈ int
		(a ≠ b) && return false
	end
	return true
end

@test validate(boxes)

#»»
end#»»
@testset "Cubes" begin#««
volume(x) = ConstructiveGeometry.CornerTables.volume(x)
@test volume(mesh(cube(1))) ≈ 1
@test volume(mesh(union(cube([5,1,1]),cube([1,5,1])))) ≈ 9
@test volume(mesh(setdiff(cube(2),cube(1)))) ≈ 7
end#»»
@testset "Extrusion" begin #««1
# using ConstructiveGeometry: nvertices, nfaces
# C = vertices(circle(3.),(precision=.01,accuracy=1,symmetry=1))
# c = [Point(20*cos(i),20*sin(i)) for i in 0:.1:π]; c=[c;[Point(0.,-1.)]]
# @test (path_extrude(c, C)) != 0
# 
d=setdiff(square(15.), translate([1,1])*square(8.))
# e=linear_extrude(10.)*d
# m=mesh(e)
# @test nvertices(m) == 16
# @test nfaces(m) == 32
end
# @testset "Types" begin #««1
# @testset "Basic types" begin #««2
# V = Point(1,2)
# @test V.coords === SA[1,2]
# end
# @testset "Conversion to/from Clipper.jl" begin #««2
# V = Point(1,2)
# test_from_to(Int, 3)
# test_from_to(Int, V)
# test_from_to(Int, [V, V])
# test_from_to(_FIXED, 3)
# test_from_to(_FIXED, Point{2,_FIXED}(1,2))
# test_from_to(Float64, 3)
# end
# end
# 
# @testset "Handling of objects" begin #««1
# s = square(1)
# @testset "Primitives" begin #««2
# @test s == square([1,1])
# end
# @testset "Operations" begin #««2
# @test union(s, union()) === s
# @test length(children(union(s, union(s)))) == 2
# end
# @testset "Transforms" begin #««2
# @test 2s == scale(2, s)
# @test scale(2)*s == scale(2, s)
# @test scale(2)*[s] == scale(2, s)
# @test color("red")*s == color("red", s)
# end
# end
# @testset "Clipper" begin #««1
# s = square(1)
# # FIXME
# end
# @testset "Convex hull" begin #««1
# using ConstructiveGeometry: convex_hull, convex_hull_list
# P(x...) = Point(Float64.(x)...)
# CH = convex_hull([P(0,0,0),P(0,0,10),P(10,0,0),P(0,10,0),P(1,1,1),P(1,0,0),])
# @test Set(CH[1]) == Set([P(0,0,0),P(0,0,10),P(10,0,0),P(0,10,0),])
# @test length(CH[2]) ==4
# @test convex_hull_list([ P(0., 0), P(1., 1), P(-1., 1), P(.2, 2.), P(0, .8), P(-.2, .8), ]) == [1,2,4,3]
# # @test convex_hull_list(ConstructiveGeometry.rows(SA[0. 0;1. 1;-1. 1;.2 2.;0 .8;-.2 .8])) == [1,2,4,3]
# @test convex_hull_list([
# 	P(-2.627798062316817, 1.075268817204301),
# 	P(-0.5030257403564974, -1.720430107526882),
# 	P(0.7927283156659947, 2.7956989247311825),
# 	P(0.0, 2.396978520135108),
# 	P(0.0, 0.03278003249806397),
# 	]) == [2,5,3,1]
# @test convex_hull_list([
# 	P(-2.150537634408602, 1.3494700327417308),
# 	P(-0.4301075268817205, -2.097580910437773),
# 	P(2.3655913978494625, 0.04739817471346019),
# 	P(0.0, 0.0),
# 	P(2.3038140933536018, 0.0),
# 	P(0.0, 0.7294358146330306),
# 	])== [2,5,3,6,1]
# end
# @testset "Triangulation" begin #««1
# using ConstructiveGeometry: LibTriangle, triangulate, triangulate_loop, identify_polygons, PolygonXor
# @test length(triangulate_loop(Point{2,Float64}.([(0,0),(1,0),(1,1),(0,1)]))) == 2
# square_with_hole = PolygonXor([[5.0, 5.0], [0.0, 5.0], [0.0, 0.0], [5.0, 0.0]],[[1.0, 1.0], [1.0, 3.0], [3.0, 3.0], [3.0, 1.0]])
# @test identify_polygons(square_with_hole) == [1, -1]
# @test length(triangulate(square_with_hole)) == 8
# @test CG.ladder_triangles(5,4,0,10) ==
# 	[(0,10,1),(1,10,11),(1,11,2),(2,11,12),(2,12,3),(3,12,13),(3,13,4)]
# end
# @testset "Basic geometry" begin#««1
# v = [[-1,0],[0,-1],[1,0],[0,1]]
# m = [CG.circular_sign(i,j) for i in v, j in v]
# @test sign.(m) == [0 1 1 1; -1 0 1 1; -1 -1 0 1; -1 -1 -1 0]
# end
# @testset "Triangle intersection" begin#««1
# TI=CG.TriangleIntersections
# Pos=TI.Constants
# pts(m) = ([Float64.(m[i,:]) for i in 1:size(m,1)]...,)
# function check_inter1(s, t, n)#««
# 	it = inter(s,t)
# 	return length(it)==n &&
# 		all(u[1] ==  pttype(p, s) && u[2] == pttype(p,t) for (p,u) in it)
# end#»»
# function check_inter(s1,t1,n)#««
# 	perm = (([1],),
# 		([1,2],[2,1]),
# 		([1,2,3],[2,3,1],[3,1,2],[1,3,2],[2,1,3],[3,2,1]))
# 	s = (s1[p] for p in perm[length(s1)])
# 	if length(s1) == 2
# 		t = (t1[p] for p in perm[3][1:3])
# 	else
# 		t = (t1[p] for p in perm[length(t1)])
# 	end
# 	return all(check_inter1(s2,t2,n) for s2 in s, t2 in t)
# end#»»
# function pttype(a,(u,v)::NTuple{2},ε=1e-8)#««
# 	abs(TI.det2(u,v,a)) > ε && return -1
# 	t = dot(a-u,v-u)/dot(v-u,v-u)
# 	t < 0-ε && return Pos.invalid
# 	t > 1+ε && return Pos.invalid
# 	t ≤ 0+ε && return Pos.vertex1
# 	t ≥ 1-ε && return Pos.vertex2
# 	return Pos.edge12
# end#»»
# function pttype(a,(p,q,r)::NTuple{3,<:StaticVector{2}},ε=1e-8)#««
# 	(d1,d2,d3) = (TI.det2(q,r,a),TI.det2(r,p,a),TI.det2(p,q,a))
# 	s(x) = (x > ε) ? '+' : (x<-ε) ? '-' : '0'
# 	m = s(d1)s(d2)s(d3)
# 	@assert '+' ∈ m "bad point type ($m)for $a in ($p, $q, $r)"
# 	m == "+++" && return Pos.interior
# 	m == "0++" && return Pos.edge23
# 	m == "+0+" && return Pos.edge31
# 	m == "++0" && return Pos.edge12
# 	m == "+00" && return Pos.vertex1
# 	m == "0+0" && return Pos.vertex2
# 	m == "00+" && return Pos.vertex3
# 	return Pos.invalid
# end#»»
# @inline inter(s::NTuple{2}, t::NTuple{3}) = TI.inter_segment2_triangle2(s,t)
# @inline inter(s::NTuple{3,<:StaticVector{3}}, t::NTuple{3}) = TI.inter(s,t)
# t2=pts(SA[0 0;10 0;0 10])
# 
# # vertices:
# @test check_inter(pts(SA[-1 1;10 0]), t2, 2)
# @test check_inter(pts(SA[-1 1;0 10]), t2, 1)
# @test check_inter(pts(SA[-5 5;0 0]), t2, 1)
# # edges:
# @test check_inter(pts(SA[-1 1;2 0]), t2, 2)
# @test check_inter(pts(SA[-1 1;0 2]), t2, 1)
# @test check_inter(pts(SA[-1 1;7 3]), t2, 2)
# # contact with vertices:
# @test check_inter(pts(SA[-5 5;1 -1]), t2, 1)
# # inside:
# @test check_inter(pts(SA[-1 1;1 1]), t2, 2)
# # outside:
# @test check_inter(pts(SA[-1 1;2 10]), t2, 2)
# @test check_inter(pts(SA[-1 5;2 20]), t2, 1)
# @test check_inter(pts(SA[-1 5;2 30]), t2, 0)
# @test check_inter(pts(SA[-1 5;10 -3]), t2, 2)
# function pttype(a,(p,q,r)::NTuple{3,<:StaticVector{3}},ε=1e-8)#««
# 	normal2 = cross(q-p,r-p)
# 	if abs(dot(normal2, p-a)) > ε
# 		return -1
# 	end
# 	proj = TI.Projector(normal2)
# 	(a1, p1, q1, r1) = proj.((a,p,q,r))
# 	return pttype(a1, (p1,q1,r1), ε)
# end#»»
# t3 = pts(SA[0 0 0;5 0 0;0 5 0])
# function dotest(c, a, b, n)#««
# 	# triangle on plane [x+y=c],
# 	# crossing [z=0] at (3-a,a,0) and (3-b,b,0)
# 	# (edges 3 and 2)
# 	t1 = pts(SA[c-a a 1; c-a a -1; c-2b+a 2b-a  -1])
# 
# 	check_inter(t3,t1,n)
# 	check_inter(t1,t3,n)
# end#»»
# @test check_inter(pts(SA[2. 0 0;2 2 0;1 1 1]),
# 		pts(SA[3. 2 0;1 2 0;2 1 1]), 2)
# 	# touch configuration:
# @test dotest(0, -2, -1, 0)
# @test dotest(0, -2, 0,  1)
# @test dotest(0, -2, 1, 1)
# @test dotest(0, 0, 1, 1)
# @test dotest(0, 1, 2, 0)
# 
# 	# general position:
# @test dotest(3, -2,-1, 0)
# @test dotest(3, -1, 0, 1)
# @test dotest(3, -1, 2, 2)
# @test dotest(3, -1, 3, 2)
# @test dotest(3, -1, 4, 2)
# 
# @test dotest(3, 0, 2, 2)
# @test dotest(3, 0, 3, 2)
# @test dotest(3, 0, 4, 2)
# 
# @test dotest(3, 1, 2, 2)
# @test dotest(3, 1, 3, 2)
# @test dotest(3, 1, 4, 2)
# 
# @test dotest(3, 3, 4, 1)
# @test dotest(3, 5, 4, 0)
# 
# 	# arrow configuration:
# @test check_inter(pts(SA[-1 1 1;-1 1 -1;2 2 0]), t3, 2)
# 
# 	# border configuration:
# @test check_inter(pts(SA[0 0 0;5 0 0;0 0 1]), t3, 2)
# @test check_inter(pts(SA[-1 1 0; 2 2 0;0 0 1]), t3, 2)
# @test check_inter(pts(SA[1 0 0;3 0 0;0 0 1]), t3, 2)
# 
# 	# coplanar
# @test check_inter(pts(SA[1 1 0;2 1 0;1 2 0]), t3, 3)
# @test check_inter(pts(SA[-1 2 0;2 -1 0;6 6 0]), t3, 6)
# @test check_inter(pts(SA[28. 0 0;12 0 0;20 0 -8]),
# 		pts(SA[25. 0 0;20 0 -5;15 0 0]), 2)
# end
# @testset "Surfaces" begin #««1
# using ConstructiveGeometry: Surface, merge, mesh
# # using ConstructiveGeometry: nvertices, nfaces
# v=[[0,-1],[1,0],[0,1],[-1,0]]
# for j in eachindex(v), i in 1:j-1
# 	@test ConstructiveGeometry.circular_lt(v[i], v[j])
# end
# nvf(s) = (CG.nvertices(s), CG.nfaces(s))
# p=Surface([[-1.,0,0],[1,0,0],[0,1,0],[0,0,1]],
#          [[3,2,1],[4,1,2],[4,3,1],[4,2,3]])
# @test nvf(mesh(2p \ p)) == (8, 12)
# # @test connected_components([:a, :b, :c, :d, :e], [[1,2],[1,3],[4,5]]) ==
# # 	[([:a, :b, :c], [[1,2],[1,3]]),
# # 	 ([:d, :e], [[1,2]]) ]
# function pyramid(t=[0,0,0], n=0)
# 	points = ([ t, t+[2.,0,0], t+[2,2,0], t+[0,2,0], t+[1,1,1]])
#   faces = [[4,3,2],[2,1,4],[1,2,5],[2,3,5],[3,4,5],[4,1,5]]
# 	faces1 = [ f .+ n for f in faces ]
# 	return surface(points, faces1)
# end
# p1 = pyramid()
# p2 = pyramid([1,0,0])
# u12 = mesh(p1 ∪ p2)
# i12 = mesh(p1 ∩ p2)
# d12 = mesh(p1 \ p2)
# @test nvf(u12) == (14, 24)
# @test nvf(i12) == (8, 12)
# @test nvf(d12) == (8, 12)
# 
# function tetrahedron(;center=[0,0,0], radius=1)
# 	return surface([center+radius*[1,0,-.5], center+radius*[-1,0,-.5],
# 		center+radius*[0,1,.5],center+radius*[0,-1,.5]],
# 		[[1,3,4],[1,2,3],[1,4,2],[2,4,3]])
# end
# 
# t1 = tetrahedron()
# t2 = 2t1
# t3 = [3,0,0]+t1
# @test nvf(mesh(t1∪t2)) == (4,4)
# @test nvf(mesh(t1∪t3)) == (8,8)
# @test nvf(mesh(t1∩t2)) == (4,4)
# @test nvf(mesh(t1∩t3)) == (0,0)
# @test nvf(mesh(t1\t1)) == (0,0)
# @test nvf(mesh(t1\t3)) == (4,4)
# @test nvf(mesh(t1\t2)) == (0,0) # t1 ⊂ t2
# @test nvf(mesh(t2\t1)) == (8,8)
# 
# end
# @testset "Difference of rotate_extrude()" begin# ««1
# m1=mesh(rotate_extrude(10,[2,0]+circle(1)))
# m2=mesh(rotate_extrude(11,[2,0]+circle(.5)))
# @debug "###### difference ######"
# m3=mesh(m1\m2)
# @test CG.nfaces(m3) == 2*CG.nvertices(m3)
# end
# #»»1
# vim: noet ts=2 fmr=««,»»
