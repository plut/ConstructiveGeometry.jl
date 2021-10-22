# using ConstructiveGeometry
include("../src/ConstructiveGeometry.jl")
C = ConstructiveGeometry

r=20.
s1 = C.set_parameters(precision=.0001)*C.sphere(20)
s2 = C.cube(20)
k=50

m=C.mesh(union(
  [ k, 0, 0]+C.color("cyan")*(s1 ∪ s2),
	[ k, k, 0]+C.color("pink")*(s1 ∩ s2),
	[ 0, 0, 0]+C.color("lime")*(s1 \ s2),
	[ 0, k, 0]+C.color("gold")*C.hull(s1, s2)
))
