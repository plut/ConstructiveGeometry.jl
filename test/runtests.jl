using ConstructiveGeometry
using Test
using StaticArrays
using ConstructiveGeometry: _FIXED, Vec, Point, Path
using ConstructiveGeometry: from_clipper, to_clipper
using ConstructiveGeometry: children, vertices
CG = ConstructiveGeometry

function test_from_to(T, x)
	return @test from_clipper(T, to_clipper(T, x)) == x
end

@testset "Types" begin #««1
@testset "Basic types" begin #««2
V = Point(1,2)
@test V.coords === SA[1,2]
end
@testset "Conversion to/from Clipper.jl" begin #««2
V = Point(1,2)
test_from_to(Int, 3)
test_from_to(Int, V)
test_from_to(Int, [V, V])
test_from_to(_FIXED, 3)
test_from_to(_FIXED, Point{2,_FIXED}(1,2))
test_from_to(Float64, 3)
end
end

@testset "Handling of objects" begin #««1
s = square(1)
@testset "Primitives" begin #««2
@test s == square([1,1])
end
@testset "Operations" begin #««2
@test union(s, union()) === s
@test length(children(union(s, union(s)))) == 2
end
@testset "Transforms" begin #««2
@test 2s == scale(2, s)
@test scale(2)*s == scale(2, s)
@test scale(2)*[s] == scale(2, s)
@test color("red")*s == color("red", s)
end
end
@testset "Clipper" begin #««1
s = square(1)
# FIXME
end
@testset "Extrusion" begin #««1
using ConstructiveGeometry: nvertices, nfaces
C = vertices(circle(3.),(precision=.01,accuracy=1,symmetry=1))
c = [Point(20*cos(i),20*sin(i)) for i in 0:.1:π]; c=[c;[Point(0.,-1.)]]
@test (path_extrude(c, C)) != 0

d=difference(square(15), translate([1,1])*square(8))
e=linear_extrude(10)*d
m=mesh(e)
@test nvertices(m) == 16
@test nfaces(m) == 32
end
@testset "Convex hull" begin #««1
using ConstructiveGeometry: convex_hull, convex_hull_list
P(x...) = Point(Float64.(x)...)
CH = convex_hull([P(0,0,0),P(0,0,10),P(10,0,0),P(0,10,0),P(1,1,1),P(1,0,0),])
@test Set(CH[1]) == Set([P(0,0,0),P(0,0,10),P(10,0,0),P(0,10,0),])
@test length(CH[2]) ==4
@test convex_hull_list([ P(0., 0), P(1., 1), P(-1., 1), P(.2, 2.), P(0, .8), P(-.2, .8), ]) == [1,2,4,3]
# @test convex_hull_list(ConstructiveGeometry.rows(SA[0. 0;1. 1;-1. 1;.2 2.;0 .8;-.2 .8])) == [1,2,4,3]
@test convex_hull_list([
	P(-2.627798062316817, 1.075268817204301),
	P(-0.5030257403564974, -1.720430107526882),
	P(0.7927283156659947, 2.7956989247311825),
	P(0.0, 2.396978520135108),
	P(0.0, 0.03278003249806397),
	]) == [2,5,3,1]
@test convex_hull_list([
	P(-2.150537634408602, 1.3494700327417308),
	P(-0.4301075268817205, -2.097580910437773),
	P(2.3655913978494625, 0.04739817471346019),
	P(0.0, 0.0),
	P(2.3038140933536018, 0.0),
	P(0.0, 0.7294358146330306),
	])== [2,5,3,6,1]
end
@testset "Triangulation" begin #««1
using ConstructiveGeometry: LibTriangle, triangulate, triangulate_loop, identify_polygons, PolygonXor
@test length(triangulate_loop(Point{2,Float64}.([(0,0),(1,0),(1,1),(0,1)]))) == 2
square_with_hole = PolygonXor([[5.0, 5.0], [0.0, 5.0], [0.0, 0.0], [5.0, 0.0]],[[1.0, 1.0], [1.0, 3.0], [3.0, 3.0], [3.0, 1.0]])
@test identify_polygons(square_with_hole) == [1, -1]
@test length(triangulate(square_with_hole)) == 8
@test CG.ladder_triangles(5,4,0,10) ==
	[(0,10,1),(1,10,11),(1,11,2),(2,11,12),(2,12,3),(3,12,13),(3,13,4)]
end
@testset "Basic geometry" begin#««1
v = [[-1,0],[0,-1],[1,0],[0,1]]
m = [CG.circular_sign(i,j) for i in v, j in v]
@test sign.(m) == [0 1 1 1; -1 0 1 1; -1 -1 0 1; -1 -1 -1 0]
end
@testset "Intersection" begin#««1
seg1 = CG.Segment(Point(0,0,2), Point(0,2,0))
seg2 = CG.Segment(Point(0,1,1), Point(1,0,1))
@test CG.inter(seg1, seg2) == Point(0,1,1)
end
@testset "Surfaces" begin #««1
using ConstructiveGeometry: Surface, merge, select_faces, mesh
using ConstructiveGeometry: nvertices, nfaces
v=[[0,-1],[1,0],[0,1],[-1,0]]
for j in eachindex(v), i in 1:j-1
	@test ConstructiveGeometry.circular_lt(v[i], v[j])
end
# @test connected_components([:a, :b, :c, :d, :e], [[1,2],[1,3],[4,5]]) ==
# 	[([:a, :b, :c], [[1,2],[1,3]]),
# 	 ([:d, :e], [[1,2]]) ]
function pyramid(t=[0,0,0], n=0)
	points = ([ t, t+[2.,0,0], t+[2,2,0], t+[0,2,0], t+[1,1,1]])
  faces = [[4,3,2],[2,1,4],[1,2,5],[2,3,5],[3,4,5],[4,1,5]]
	faces1 = [ f .+ n for f in faces ]
	return surface(points, faces1)
end
nvf(s) = (CG.nvertices(s), CG.nfaces(s))
p1 = pyramid()
p2 = pyramid([1,0,0])
u12 = mesh(p1 ∪ p2)
i12 = mesh(p1 ∩ p2)
d12 = mesh(p1 \ p2)
@test nvf(u12) == (14, 24)
@test nvf(i12) == (8, 12)
@test nvf(d12) == (8, 12)

function tetrahedron(;center=[0,0,0], radius=1)
	return surface([center+radius*[1,0,-.5], center+radius*[-1,0,-.5],
		center+radius*[0,1,.5],center+radius*[0,-1,.5]],
		[[1,3,4],[1,2,3],[1,4,2],[2,4,3]])
end

t1 = tetrahedron()
t2 = 2t1
t3 = [3,0,0]+t1
@test nvf(mesh(t1∪t2)) == (4,4)
@test nvf(mesh(t1∪t3)) == (8,8)
@test nvf(mesh(t1∩t2)) == (4,4)
@test nvf(mesh(t1∩t3)) == (0,0)
@test nvf(mesh(t1\t1)) == (0,0)
@test nvf(mesh(t1\t3)) == (4,4)
@test nvf(mesh(t1\t2)) == (0,0) # t1 ⊂ t2
@test nvf(mesh(t2\t1)) == (8,8)

end
#»»1

# vim: noet ts=2 fmr=««,»»
