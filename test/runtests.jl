using Test
using StaticArrays
using Colors

push!(LOAD_PATH, "./src", "src")
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
@testset "Clipper" begin
s = Square(1)
end

#>>>1
# vim: noet ts=2 fmr=<<<,>>>

