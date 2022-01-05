using ConstructiveGeometry
using StaticArrays
using FastClosures
using GLMakie
using Colors

b1(p0,u0,t0,p1,u1,t1,n) = ConstructiveGeometry.bezier(p0,p0+t0*u0,p1+t1*u1,p1,n)

# first define a half-fleurdelis shape using Bézier curves
demi_lis(n=15) = [0, 17.4] + .3*[1 0;0 -1]*([-271,-277]+polygon([
	b1([271,165],[20,18],0.5,[288,234],[10,-14],1.5,n)...,
	b1([288,234],[10,-9],1.,[342,249],[0,-10],3.,n)...,
	b1([342,249],[0,10],1.,[320,274],[15,0],1.,n)...,
	b1([320,274],[-15,-15],2.,[287,277],[0,-10],1.,n)...,
	[287,290],
	[360,290],
	[360,340],
	[271,340],
	]))
lis_w = 22.6
lis_t = 1.6
lis(n=15) = (demi_lis(n) ∪ [-1 0;0 1]*demi_lis(n))
# 15 points seems enough, compile the polygon for the fleurdelis:
lis_poly = polygon(lis(15))

tol = 0.2;
lozenge = polygon([[0,-4],[6,0],[0,4],[-6,0]]);
lozenge1 = offset(tol)*lozenge
# define the flat 3d shape of a section of the crown:
lis3d = surface([0 0 1;1 0 0;0 1 0]*(set_parameters(atol=1e-2,rtol=1e-3)*
  intersect(union(
		linear_extrude(lis_t)*lis_poly,
		sweep(lis_poly; miter_limit=10.)*half(:top)*circle(2),
# 	linear_extrude(1.6)*offset(-5)*lis
		),
	# reserve space for inset stones:
	~([0,7.5,.8]+cylinder(10,4+tol)),
	~([lis_w,7.5,.8]+linear_extrude(10)*lozenge1),
	~([-lis_w,7.5,.8]+linear_extrude(10)*lozenge1),
	translate([-lis_w, 0, 0])*cube([2*lis_w, 1000, 100]))))

# and deform it by wrapping; ConstructiveGeometry wraps on a cylinder,
# here we use a one-sheet revolution hyperboloid:
function hyperboloid(r, (x,y,z))
	θ = y/r
	k = 1+(z+17.4)^2/2e4
	return SA[(x+r)*cos(θ)*k, (x+r)*sin(θ)*k, z]
end

N = 8
R = N*lis_w/π-.01
def = deform(@closure(p->hyperboloid(R, p));
	distance2=@closure((p,q)->ConstructiveGeometry.wrapsagitta(R, p,q)))
lis3d_def = def*lis3d
lis3d_def1 = surface(lis3d_def)
crown=union((rotate(45*i)*lis3d_def1 for i in 1:8)...)
crown1 = surface(crown)


pyramid = [0 0 1;1 0 0;0 1 0]*([lis_w,7.5,.8]+cone(3)*lozenge)
pyramid_stone = def*pyramid

oval = [0 0 1;1 0 0;0 1 0]*([0,7.5,.8]+half(:top)*sphere(4))
oval_stone = def*oval


# display a nice color image:
gold = colorant"gold"
ruby = coloralpha(colorant"crimson", .5)
emerald = coloralpha(colorant"chartreuse2", .5)

fullcrown = union(gold*crown1,
	ruby*union((rotate(45*i)*(def*oval) for i in 1:8)...),
	emerald*union((rotate(45*i)*(def*pyramid) for i in 1:8)...))

save("crown.png", plot(fullcrown))

# save("crown-main.stl", crown1)
# save("crown-pyramid.stl", rotate(90,axis=[0,-1,0])*pyramid_stone)
# save("crown-oval.stl", rotate(90,axis=[0,-1,0])*oval_stone)

# save("crown-main.ply", crown1)
# save("crown-pyramid.ply", rotate(90,axis=[0,-1,0])*pyramid_stone)
# save("crown-oval.ply", rotate(90,axis=[0,-1,0])*oval_stone)

nothing
