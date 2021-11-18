using ConstructiveGeometry
using GLMakie
# using CairoMakie

hexagon = polygon([[cospi(t/3),sinpi(t/3)] for t in 0:5])
c1, c2, c3 = parse.(ConstructiveGeometry.Colorant, ( "#cb3c33", "#9558b2", "#389826"))

bolt = union(linear_extrude(5)*(8*hexagon),
	cylinder(15,4), rotate_extrude(7*360, slide=14)*translate([1,0])*square(4,1))

m = union(c1*bolt, [20,0,0]+c2*bolt, [10,17,0]+c3*bolt)

s = plot(m)
zoom!(s, Vec3f0(0,0,0), 5., false)
s.center = false
save("logo.png", s)
