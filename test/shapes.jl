include("../src/Shapes.jl")
using Test
using StaticArrays

PX = Shapes.PolygonXor{Float64}
rect(a,b,c=0,d=0) = PX([SA[c,d],SA[c+a,d],SA[c+a,d+b],SA[c,d+b]])
r1 = rect(2,1)
r2 = rect(1,2)
r3 = rect(3,3)
r4 = rect(1,1,4,4)
nva(s::Shapes.PolygonXor) = (length.(s.paths), Shapes.area(s))

@testset "Union" begin#««
@test nva(Shapes.clip(:union,r1,r2)) == ([6], 3)
@test nva(Shapes.clip(:union,r1,r3)) == ([4], 9)
@test nva(Shapes.clip(:union,r1,r4)) == ([4,4], 3)
end#»»
@testset "Intersection" begin#««
@test nva(Shapes.clip(:intersection,r1,r2)) == ([4], 1)
@test nva(Shapes.clip(:intersection,r2,r1)) == ([4], 1)
@test nva(Shapes.clip(:intersection,r1,r3)) == ([4], 2)
@test nva(Shapes.clip(:intersection,r3,r1)) == ([4], 2)
@test nva(Shapes.clip(:intersection,r1,r4)) == ([], 0)
end#»»
@testset "Difference" begin#««, 0)
@test nva(Shapes.clip(:difference,r1,r2)) == ([4], 1)
@test nva(Shapes.clip(:difference,r3,r1)) == ([6], 7)
end#»»
@testset "Convex hull" begin#««
@test nva(Shapes.convex_hull(r1,r2)) == ([5], 3.5)
@test nva(Shapes.convex_hull(r1,r4)) == ([6], 11)
end#»»

nothing
