```@meta
CurrentModule = ConstructiveGeometry
```
```@setup 0
using ConstructiveGeometry
using Makie
using CairoMakie
png(name, s) = save(name*".png", Makie.plot(s));
```
# [`ConstructiveGeometry.jl`](https://github.com/plut/ConstructiveGeometry.jl) Documentation

This package provides tools for describing 3d objects in Julia.
For example, this is the Julia code used to draw the logo of this
page:
```julia
using ConstructiveGeometry
using CairoMakie
using Colors

hexagon = polygon([[cospi(t/3),sinpi(t/3)] for t in 0:5])
c1, c2, c3 = colorant"#cb3c33", colorant"#9558b2", colorant"#389826"

bolt = linear_extrude(5)*(8*hexagon) ∪ cylinder(15,4) ∪
	rotate_extrude(7*360, slide=14)*translate([1,0])*square(4,1)

m = union(c1*bolt, [20,0,0]+c2*bolt, [10,17,0]+c3*bolt)

save("logo.png", Makie.plot(m))
```

# Overview

Drawing an object is generally done in two steps.
First, an abstract, “ideal” geometric object
(in either two or three dimensions)
is defined using the following constructions:

 - [primitive geometric objects](@ref primitives), such as cubes,
   spheres, explicit polygons, etc.;
 - [geometric transformations](@ref transformations) acting on one
   object, such as (invertible) affine transformations, extrusions,
   projections, color-change, etc.;
 - [CSG operations](@ref operations) combining several objects,
   such as boolean operations or Minkowski sum.

Such an abstract object is stored as a tree with primitive objects as
leaves:
```@repl 0
hexagon = polygon([[cospi(t/3),sinpi(t/3)] for t in 0:5]); # hide
bolt = linear_extrude(5)*(8*hexagon) ∪ cylinder(15,4) ∪ rotate_extrude(7*360, slide=14)*translate([1,0])*square(4,1); # hide
display(bolt)
```


Any geometric object defined in this way can then be
[instantiated as an explicit mesh](@ref meshing).
The mesh may be visualized directly within Julia (using `Makie`)
or [exported](@ref io) in several formats,
including as an STL (for 3d objects) or SVG (for 2d objects) file.
