# Solids.jl Documentation

!!! warning

This package is very much a work-in-progress. Right now only very basic
functionality is available (describing geometries, and some cases for
meshing). The code is absolutely not stable, and any contributions are welcome.


This package provides tools for describing 3d objects
in Julia.
This includes both geometry functions
and a syntax describing constructive geometry.
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

The package is made of two parts:

 - a system for describing geometric objects;
 - back-ends for converting these objects to useful formats.

As of 2021-01, the only useable output format is an OpenSCAD file.

Other planned output formats include:
 - (**TODO**) represented graphically using one of the Julia plotting packages;
 - (**TODO**) converted to either a mesh or a signed distance field;
 - (**IN PROGRESS**) directly exported as a 2d (`.svg`)
 or 3d file (`.stl` or `.ply`).


# Quick-start

## Basic example
```julia
using Solids
import Solids: Square, Circle, Surface, mult_matrix, translate, scale, color

s1 = union(
  color("pink")*
  translate([3,0])*
  scale([2,1])*
  Circle(3),

  color("cyan")*
  translate([0,5])*
  Square([2,3])
)

Surface(s1)
```

## I/O

```@docs
Solids.include
Solids.scad
```
