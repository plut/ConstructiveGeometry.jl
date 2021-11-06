# [CSG operations](@id operations)
```@meta
CurrentModule = ConstructiveGeometry
```
```@setup 0
using ConstructiveGeometry
using Makie
using CairoMakie
png(name, s) = save(name*".png", Makie.plot(s));
```

## Boolean operations
```@docs
union(::AbstractGeometry,::AbstractGeometry)
```
```@repl 0
s = union(cube(50), sphere(50));
png("union", s); # hide
```
![example: union of a sphere and a cube](union.png)

```@docs
intersect(::AbstractGeometry{D},::AbstractGeometry{D}) where{D}
```
```@repl 0
s = intersect(cube(50), sphere(50));
png("intersection", s); # hide
```
![example: intersection of a sphere and a cube](intersection.png)
```@docs
setdiff(::AbstractGeometry{D},::AbstractGeometry{D}) where{D}
```
```@repl 0
s = setdiff(cube(50), [50,0,0]+sphere(50));
png("setdiff", s); # hide
```
![example: difference of a sphere and a cube](setdiff.png)
```@docs
complement
```
Complements are provided as a shortcut to simplify “subtractive”
operations, *i.e.* chains of intersections and differences.
See [Three-dimensional embeddings of two-dimensional objects](@ref embed).

## Convex hull
```@docs
hull
```
```@repl 0
s = hull(cube(50), [50,0,0]+sphere(50));
png("hull", s); # hide
```
![example: convex hull of a sphere and a cube](hull.png)

In the case of mixed dimensions, two-dimensional objects are understood
as included in the horizontal plane, unless they have been subjected
to a three-dimensional affine transformation; in that case,
this transformation is applied to their vertices.
See [Three-dimensional embeddings of two-dimensional objects](@ref embed).

## Minkowski sum
```@docs
minkowski
```
```@repl 0
c = cube(10);
s1 = minkowski(square(50), circle(20));
s2 = minkowski(c, cube(20,1,1));
s3 = minkowski(c, polygon([[0,0],[0,30],[30,0]]));
png("minkowski_square_circle",s1); # hide
png("minkowski_cube_cube",s2); # hide
png("minkowski_cube_polygon", s3); # hide
```

The Minkowski sum between a polygon and a circle of radius `r`
is the same as the offset of the polygon by this radius:
![example: Minkowski sum of a square and a circle](minkowski_square_circle.png)

Minkowski sum between volumes is allowed; e.g. the Minkowski sum of two
axis-aligned parallelepipeds is again a parallelepiped:
![example: Minkowski sum of two axis-aligned cubes](minkowski_cube_cube.png)

Minkowski sum between a volume and a polygon is also allowed;
here the polygon is a triangle in the horizontal plane:
![example: Minkowski sum of a cube and a polygon](minkowski_cube_polygon.png)

## Slicing and projection

Slicing and projection convert a volume to a shape.
These transformations are only defined with respect to horizontal planes,
since these are the only planes in which canonical `(x,y)` coordinates
are defined.

To use another plane, say the image of the horizontal plane by a rotation
`R`, apply the inverse rotation of `R` to the object to bring the
situation back to the horizontal plane.

### `slice`
```@docs
slice
```
```@repl 0
s = slice()*setdiff(sphere(20),sphere(18));
png("slice", s); # hide
```
![example: slicing a hollow sphere](slice.png)

### `project`
```@docs
project
```
```@repl 0
s = project()*setdiff(sphere(20),sphere(18));
png("project", s); # hide
```
![example: projecting a hollow sphere](project.png)

## Intersection with half-space
```@docs
halfspace
```
```@repl 0
s = halfspace([0,0,-1],[0,0,0])*setdiff(sphere(20),sphere(18));
png("halfspace", s); # hide
```
![example: one half of a hollow sphere](halfspace.png)

