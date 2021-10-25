# `ConstructiveGeometry.jl` Documentation

This package provides tools for describing 3d objects in Julia.
For example, this is the Julia code used to draw the logo of this
package:
````@eval
Markdown.parse("""
```julia
$(readstring("../logo.jl"))
```
""")
````

# Quick-start

This package defines three kinds of abstract geometric objects
(in either two or three dimensions):

 - [primitive geometric objects](@ref primitives), such as cubes,
   spheres, etc.;
 - [geometric transformations](@ref transformations) acting on one
   object, such as (invertible) affine transformations, extrusions,
   projections, color-change, etc.;
 - [CSG operations](@ref operations) combining several objects,
   such as boolean operations or Minkowski sum.

Any geometric object defined in this way may then be instantiated as an
explicit mesh. The mesh can be visualized directly within Julia (using
`Makie`) or exported as an STL (for 3d objects) or SVG (for 2d objects)
file.
