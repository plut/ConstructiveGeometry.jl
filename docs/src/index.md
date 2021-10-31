# `ConstructiveGeometry.jl` Documentation

This package provides tools for describing 3d objects in Julia.
For example, this is the Julia code used to draw the logo of this
page:
```julia
using ConstructiveGeometry
using CairoMakie

hexagon = polygon([[cospi(t/3),sinpi(t/3)] for t in 0:5])
c1, c2, c3 = parse.(ConstructiveGeometry.Colorant, ( "#cb3c33", "#9558b2", "#389826"))

bolt = linear_extrude(5)*(8*hexagon) ∪ cylinder(15,4)

m = union(c1*bolt, [20,0,0]+c2*bolt, [10,17,0]+c3*bolt)

save("logo.png", Makie.plot(m))
```

# Overview

This package defines three kinds of abstract geometric objects
(in either two or three dimensions):

 - [primitive geometric objects](@ref primitives), such as cubes,
   spheres, etc.;
 - [geometric transformations](@ref transformations) acting on one
   object, such as (invertible) affine transformations, extrusions,
   projections, color-change, etc.;
 - [CSG operations](@ref operations) combining several objects,
   such as boolean operations or Minkowski sum.

Any geometric object defined in this way can then be
[instantiated as an explicit mesh](@ref meshing).
The mesh may be visualized directly within Julia (using `Makie`)
or [exported](@ref io) as an STL (for 3d objects) or SVG (for 2d objects) file.
