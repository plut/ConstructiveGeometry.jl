# using ConstructiveGeometry
include("../src/ConstructiveGeometry.jl")
C = ConstructiveGeometry

r=20.
s1 = C.set_parameters(precision=.0001)*C.sphere(r)
s2 = C.cube(r)
k=50

m=C.mesh(union(
  [ k, 0, 0]+C.color("cyan")*(s1 ∪ s2),
	[ k, k, 0]+C.color("pink")*(s1 ∩ s2),
	[ 0, 0, 0]+C.color("lime")*(s1 \ s2),
	[ 0, k, 0]+C.color("gold")*C.hull(s1, s2)
))

# T.explain(mesh(s1∪s2), "/tmp/u.scad", scale=sc, name=:mu)
# T.explain(mesh(s1∩s2), "/tmp/i.scad", scale=sc, name=:mi)
# T.explain(mesh(s1\s2), "/tmp/d.scad", scale=sc, name=:md)
# T.explain(mesh(hull(s1,s2)), "/tmp/c.scad", scale=sc, name=:mc)

#=
use <u.scad>
use <i.scad>
use <d.scad>
use <c.scad>

k=50;

mu(pos=[k,0,0], c="cyan");
mi(pos=[k,k,0], c="pink");
md(pos=[0,0,0], c="lime");
mc(pos=[0,k,0], c="gold");
=#
