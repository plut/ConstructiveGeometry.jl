# [Transformations](@id transformations)
```@meta
CurrentModule = ConstructiveGeometry
```
```@setup 0
using ConstructiveGeometry
using Makie
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
Transformations with a negative determinant will reverse the object
(reverse the polygons, or reverse the faces of meshes)
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
```@docs
rotate
```
```@repl 0
s = rotate(30)*square(20);
png("rotate", s); # hide
```
![a rotated square](rotate.png)
```@docs
mirror
```
```@docs
raise
lower
```
## Overloaded operators

The following operators are overloaded.

 - `vector + solid` is a translation.
 - `real * solid` is a scaling.
 - `complex * 2dshape` is a similitude.
 - `matrix * solid` is a linear transformation.
 - `vector * solid` is a multiplication by a diagonal matrix.
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
s = linear_extrude(10)*[square(10,5), square(5,15)];
png("linear_extrude", s); # hide
```
![example: linear extrusion of a L-shape](linear_extrude.png)
```@docs
rotate_extrude
```
```@repl 0
s = rotate_extrude(245)*[square(10,5), square(5,15)];
png("rotate_extrude", s); # hide
```
![example: rotation extrusion of a L-shape](rotate_extrude.png)
```@docs
sweep
```

A swept surface is similar to a (closed) path extrusion:

```@repl 0
s = sweep(square(50))*circle(5);
png("swept_circle", s); # hide
```
![example: a circle swept along a square](swept_circle.png)
!!! warning "Swept surfaces"

    A surface may only be swept along a closed loop
    (or the union of several closed loops) for now;
    this is a limitation of the `clipper` library,
    which [does not support single-path extrusion](http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Types/EndType.htm)
    for now (and this is unlikely to change in the near future).

```@repl 0
f(t) =([ cospi(t) -sinpi(t) 0;sinpi(t) cospi(t) 0;0 0 1],[0 0 10*t]);
s = sweep(f; nsteps=100,maxgrid=100)*cube(20);
png("swept_cube", s); # hide

```
![example: a cube swept along a helix](swept_cube.png)

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

## Coloring objects

```@docs
color
```
```@repl 0
green, red = parse.(ConstructiveGeometry.Colorant, ("green", "red"))
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

## Modifying meshing parameters
```@docs
set_parameters
```

The `set_parameters` transformation allows attaching arbitrary metadata.
This is on purpose (although there currently exists no easy way
for an user to recover these metadata while meshing an object).

The values for these parameters are explained in [Accuracy and
precision](@ref atol_rtol).

