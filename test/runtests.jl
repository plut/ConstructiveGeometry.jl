using Test, StaticArrays, Colors

# push!(LOAD_PATH, "./src", "src")
using Solids: _FIXED, Vec, Path
using Solids: from_clipper, to_clipper

function test_from_to(T, x)
	return @test from_clipper(T, to_clipper(T, x)) == x
end

@testset "Types" begin #<<<1
@testset "Basic types" begin #<<<2
V = Vec{2}(1,2)
@test V === SA[1,2]
end
@testset "Conversion to/from Clipper.jl" begin #<<<2
V = Vec{2}(1,2)
test_from_to(Int, 3)
test_from_to(Int, V)
test_from_to(Int, [V, V])
test_from_to(_FIXED, 3)
test_from_to(_FIXED, V)
test_from_to(Float64, 3)
end
end

@testset "Handling of objects" begin #<<<1
using Solids: Square, Circle
using Solids: children
using Solids: mult_matrix, translate, scale
using Solids: color
s = Square(1)
@testset "Primitives" begin #<<<2
@test s == Square([1,1])
end
@testset "Operations" begin #<<<2
@test union(s, union()) === s
@test length(children(union(s, union(s)))) == 2
end
@testset "Transforms" begin #<<<2
@test 2s == scale(2, s)
@test scale(2)*s == scale(2, s)
@test scale(2)*[s] == scale(2, s)
@test color("red", s) == color(parse(Colorant, "red"), s)
end
end
@testset "Clipper" begin #<<<1
s = Square(1)
# FIXME
end
@testset "Extrusion" begin #<<<1
using Solids: path_extrude, points
C = points(Circle(3.),(precision=.01,accuracy=1))
c = [SA[cos(i),sin(i)] for i in 0:.1:π]; c=[c;[SA[0,-1]]]
@test (path_extrude(20c, C)) != 0
end
@testset "Convex hull" begin #<<<1
using Solids: convex_hull
CH = convex_hull([SA[0,0,0],SA[0,0,10],SA[10,0,0],SA[0,10,0],SA[1,1,1],SA[1,0,0]])
@test Set(CH[1]) == Set([SA[0,0,0],SA[0,0,10],SA[10,0,0],SA[0,10,0]])
@test length(CH[2]) ==4
end
@testset "Surfaces" begin #<<<1
using Solids: connected_components
@test connected_components([:a, :b, :c, :d, :e], [[1,2],[1,3],[4,5]]) ==
	[([:a, :b, :c], [[1,2],[1,3]]),
	 ([:d, :e], [[1,2]]) ]
end
#>>>1

# vim: noet ts=2 fmr=<<<,>>>

