# `ConstructiveGeometry.jl` Documentation

This package provides tools for describing 3d objects in Julia.
This includes mainly a syntax for building a CSG tree
and functions for representing the objects.
This syntax is inspired by OpenSCAD, but is actual Julia code:
```julia
using ConstructiveGeometry

square(20)

linear_extrude(30) * [
  intersection(
    translate([10,0]) * circle(3),
    translate([13,0]) * circle(3),
  ),
  color("pink") * scale(2) * square(1),
]

```

# Quick-start

## Basic example
```julia
using ConstructiveGeometry

s1 = union(
  color("pink")*
  translate([3,0])*
  scale([2,1])*
  circle(3),

  color("cyan")*
  translate([0,5])*
  square([2,3])
)

mesh(s1)
```
