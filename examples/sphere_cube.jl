using ConstructiveGeometry

r=20.
s1 = sphere(20)
s2 = cube(20)
k=50

u =union(
  [ k, 0, 0]+color("cyan")*(s1 ∪ s2),
	[ k, k, 0]+color("pink")*(s1 ∩ s2),
	[ 0, 0, 0]+color("lime")*(s1 \ s2),
	[ 0, k, 0]+color("gold")*hull(s1, s2)
)
save("sphere_cube.png", u)
