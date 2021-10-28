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
```@docs
intersect(::AbstractGeometry{D},::AbstractGeometry{D}) where{D}
```
```@docs
setdiff(::AbstractGeometry{D},::AbstractGeometry{D}) where{D}
```
## Convex hull
```@docs
hull
```

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

