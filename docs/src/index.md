# Solids.jl Documentation

This package provides tools for a simple syntax
for constructive solid geometry objects in Julia.
This syntax is inspired by OpenSCAD, but is actual Julia code:
```julia
using Solids

Square(20)

linear_extrude(30) * [
  intersection(
    translate([10,0]) * Circle(3),
    translate([13,0]) * Circle(3),
  ),
  color("pink") * scale(2) * Square(1),
]

```

This package is made of two parts:

 - a system for describing geometric objects;
 - back-ends for converting these objects to useful formats.

As of 2021-01, the only useable output format is an OpenSCAD file. (Since
`Solids.jl` has some constructions which do not exist in OpenSCAD, this
is available only for those objects which do not use these
constructions).

Other planned output formats include:
 - (**TODO**) represented graphically using one of the Julia plotting packages;
 - (**TODO**) converted to either a mesh or a signed distance field;
 - (**TODO**) directly exported as a 2d (`.svg`) or 3d file (`.stl` or `.ply`).


# Quick-start

## Basic example
```julia
using Solids
import Solids: Square, Circle, mult_matrix, translate, scale, color

union(
  color("pink")*
  translate([3,0])*
  scale([2,1])*
  Circle(3),

  color("cyan")*
  translate([0,5])*
  Square([2,3])
)
```

## I/O

```@docs
Solids.include
Solids.scad
```
