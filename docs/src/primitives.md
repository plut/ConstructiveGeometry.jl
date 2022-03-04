# [Primitive solids](@id primitives)
```@meta
CurrentModule = ConstructiveGeometry
```
```@setup 0
using ConstructiveGeometry
using Makie
using CairoMakie
png(name, s) = save(name*".png", Makie.plot(s));
```

`ConstructiveGeometry.jl` supports two basic families of objects:
two-dimensional shapes and three-dimensional volumes.

## Two-dimensional shapes
### Square
```@docs
square
```
```@repl 0
s = square(20,15);
png("square", s); # hide
```

![example: a square](square.png)
`square([20,15])` also works; `square(20)` produces a real square.

### Circle
```@docs
circle
```
```@repl 0
s = circle(20);
png("circle", s); # hide
```
![example: a circle](circle.png)

When a circle (and circle-like objects such as cylinders, cones,
or rotational extrusions) is converted to an explicit mesh,
it is instantiated as a regular polygon,
which number of sides depends on the `atol` and `rtol` parameters;
see [Number of vertices of circles](@ref circle_vertices).

### Stroke path
```@docs
stroke
```
```@repl 0
s = stroke([[0,0], [100,0],[100,100],[50,150],[0,100]],10);
s1 = [120,0]+ stroke([[0,0], [100,0],[100,100],[50,150],[0,100]],10;ends=:loop,join=:square);
png("stroke", s ∪ s1); # hide
```
![example: a stroked path](stroke.png)
### Polygon
```@docs
polygon
```
```@repl 0
s = polygon([[0,0], [100,0],[100,100],[50,150],[0,100]]);
png("polygon", s); # hide
```
![example: a polygon](polygon.png)

## Three-dimensional volumes

### Cube
```@docs
cube
```
```@repl 0
s = cube(10,20,30);
png("cube", s); # hide
```
![example: a cube](cube.png)

### Cone
```@docs
cone
```
```@repl 0
s = cone(50,10);
png("cone", s); # hide
```
![example: a cone](cone.png)

### Cylinder
```@docs
cylinder
```
```@repl 0
s = cylinder(50,10);
png("cylinder", s); # hide
```
![example: a cylinder](cylinder.png)

### Sphere
```@docs
sphere
```
```@repl 0
s = sphere(50);
png("sphere", s); # hide
```
![example: a sphere](sphere.png)

When a sphere is converted to an explicit mesh,
a Fibonacci sphere is used, with a number of vertices
depending on the `atol` and `rtol` parameters,
as explained in [Number of vertices of spheres](@ref sphere_vertices).

TODO: allow various other models of spheres (icosphere, etc.).

### Surface
```@docs
surface(::Any,::Any)
```
```@repl 0
s = surface([[0,0,0],[10,0,0],[10,10,0],[0,10,0],[5,5,2]],
  [(1,2,5),(2,3,5),(3,4,5),(4,1,5),(4,3,2,1)]);
png("surface",s); # hide
```
![example: a pyramidal surface](surface.png)

