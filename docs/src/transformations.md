# [Transformations](@id transformations)
```@meta
CurrentModule = ConstructiveGeometry
```
```@setup 0
using ConstructiveGeometry
using Makie
using Colors
using CairoMakie
png(name, s) = save(name*".png", Makie.plot(s));
```

All single-object transformations accept two possible syntaxes:
```julia
    transform(parameters, solid1, solid2, ...)
    transform(parameters) * solid1
```
The second, multiplicative form allows easy chaining of transformations:
```julia
    transform1(param1) * transform2(param2) * solid
```
This form may also be applied to several solids by either wrapping them in a
`union`, or equivalently, by applying it to a `Vector` of such objects:
```julia
    transform(parameters) * [ solid1, solid2, ... ]
```

## Affine transformations
```@docs
mult_matrix
```
```@repl 0
s = mult_matrix([1 0 0;0 1 0;0 .5 1])*cube(10);
png("mult_matrix1", s); # hide
```
![a skewed cube](mult_matrix1.png)

Only invertible affine transformations are supported.
Transformations with a negative determinant reverse the object
(either reverse the polygon loops, or reverse the triangular faces of meshes)
to preserve orientation.

For non-invertible transformations, see [`project`](@ref).

### [Three-dimensional embeddings of two-dimensional objects](@id embed)
As an exception, it is allowed to apply a (2d -> 3d) transformation
to any three-dimensional object.
The result of such a transformation is still two-dimensional
(and will accordingly be rendered as a polygon),
but the information about the embedding will be used when computing
convex hull or Minkowski sum with a three-dimensional object.

```@repl 0
s = hull([30,0,0]+[1 0 0;0 1 0;.5 0 0]*circle(20), [0,0,30]);
png("embed_2d_3d", s); # hide
```
![convex hull of a non-canonically embedded circle and a point](embed_2d_3d.png)



```@docs
translate
```
```@docs
scale
```
```@repl 0
s = [1,1.5,2]*sphere(50);
png("scaled_sphere", s); # hide
```
![a diagonally scaled sphere produces an ellipsoid](scaled_sphere.png)
```@docs
rotate
```

Angles are in degrees by default.
Angles in radians are supported through the use of `Unitful.rad`.

```@repl 0
s = rotate(30)*square(20);
png("rotate", s); # hide
```
![a rotated square](rotate.png)
```@docs
reflect
```
```@docs
raise
lower
```
## Overloaded operators

The following operators are overloaded.

 - `matrix * solid` is a linear transformation.
 - `vector * solid` is a multiplication by a diagonal matrix.
 - `vector + solid` is a translation.
 - `real * solid` is a scaling.
 - `complex * 2dshape` is a similitude.
 - `color * solid` is a `color` operation.
 - `color % solid` is a `highlight` operation.

## Two-dimensional drawing
```@docs
offset
```
The grid size used for offsetting
is derived from the `atol` and `rtol` parameters,
and upper bounded by the optional `maxgrid` parameter
(if this is different from zero).

```@repl 0
s1 = offset(10)*[square(100,50), square(50,100)];
s2 = offset(3)*cube(30);
png("offset_L", s1); # hide
png("offset_cube", s2); # hide
```
![example: an offset L-shape](offset_L.png)
![example: an offset cube](offset_cube.png)
```@docs
opening
```
```@repl 0
s = opening(10)*[square(100,50), square(50,100)];
png("opening", s); # hide
```
![example: the opening of the L-shape](opening.png)
```@docs
closing
```
```@repl 0
s = closing(10)*[square(100,50), square(50,100)];
png("closing", s); # hide
```
![example: the closing of the L-shape](closing.png)

## Extrusion
### Linear extrusion
```@docs
linear_extrude
```
```@repl 0
s1 = linear_extrude(10)*[square(10,5), square(5,15)];
png("linear_extrude", s1); # hide
s2 = linear_extrude(20, twist=45, scale=.8)*[square(10,5), square(5,15)];
png("linear_extrude_twist", s2); # hide
```
![example: linear extrusion of a L-shape](linear_extrude.png)
![example: linear extrusion with twist](linear_extrude_twist.png)
```@docs
rotate_extrude
```
```@repl 0
s1 = rotate_extrude(245)*[square(10,5), square(5,15)];
png("rotate_extrude", s1); # hide
s2 = rotate_extrude(720, slide=30)*translate([10,0])*square(5);
png("slide", s2); # hide
```
![example: rotation extrusion of a L-shape](rotate_extrude.png)
![example: rotational extrusion with slide](slide.png)

The `cone` function may also be used as an operator
to build a cone out of an arbitrary shape:
```@repl 0
s = cone([1,2,3])*square(5);
png("cone_pyramid", s); # hide
```
![example: using cone to build a pyramid](cone_pyramid.png)


### Surface sweep
```@docs
sweep
```

A swept surface is similar to a (closed) path extrusion:

```@repl 0
s = sweep(square(50))*circle(5);
png("swept_circle", s); # hide
```
![example: a circle swept along a square](swept_circle.png)

```@repl 0
f(t) =([ cospi(t) -sinpi(t) 0;sinpi(t) cospi(t) 0;0 0 1],[0 0 10*t]);
s = sweep(f; nsteps=100,maxgrid=100)*cube(20);
png("swept_cube", s); # hide

```
![example: a cube swept along a helix](swept_cube.png)
### Path extrusion
```@docs
path_extrude
```

## User-defined volume deformations
```@docs
deform
```
```@repl 0
s1 = deform(p->p/(3+p[1]))*cube(5);
s2 = deform(p->p/sqrt(1+p[1]))*cube(5);
png("homographic_cube", s1); # hide
png("nonlinear_cube", s2); # hide
```
![example: homographic image of a cube](homographic_cube.png)
![example: non-linear image of a cube](nonlinear_cube.png)

### Cylindrical wrapping (experimental)

```@docs
wrap
```
```@repl 0
s = wrap(3)*cube(5);
png("wrapped_cube", s); # hide
```
![example: a cube wrapped around a cylinder](wrapped_cube.png)

## Decimation

These operations either reduce or increase the number of faces in
a three-dimensional object.

```@docs
decimate
```
```@docs
loop_subdivide
```
```@repl 0
s = loop_subdivide(4)*cube(20);
png("loop_subdivide", s); # hide
```
![example: loop subdivision of a cube](loop_subdivide.png)

```@docs
refine
```

## Coloring objects

```@docs
color
```
```@repl 0
green, red = colorant"green", colorant"red"
s = union(green * cube(10), [10,0,0]+red*sphere(10));
png("color", s); # hide
```
![example: union of a sphere and a cube](color.png)

```@docs
highlight
```
```@repl 0
s = intersect(green % cube(10), red % ([10,0,0]+sphere(10)));
png("highlight", s); # hide
```
![example: intersection of a highlighted sphere and a highlighted cube](highlight.png)

Highlighted parts of objects are shown only
when the object is represented as an image via the `plot` method.
For SVG and STL output, all highlighted parts are ignored.

Highlighted objects are preserved only by CSG operations
and (invertible) affine transformations. For other transformations:
 - convex hull and Minkowski sum are generally increasing
   transformations, and would cover highlighted parts anyway;
 - projection, slicing and extrusion modify the dimension of object,
   making it impossible to preserve highlighted parts.

```@docs
randomcolor
```
```@repl 0
s = randomcolor()*sphere(5);
png("randomcolor", s); # hide
```
![example: sphere with random-colored triangular faces](randomcolor.png)

## Modifying meshing parameters
```@docs
set_parameters
```

The `set_parameters` transformation allows attaching arbitrary metadata.
This is on purpose (although there currently exists no easy way
for an user to recover these metadata while meshing an object).

The values for these parameters are explained in [`atol` and `rtol`](@ref atol_rtol).

