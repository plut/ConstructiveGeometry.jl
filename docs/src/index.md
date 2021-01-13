# Solids.jl Documentation

This package presents a set of types and syntax providing tools for a
simple description of constructive solid geometry objects in Julia.
Syntax is inspired by OpenSCAD, but is actual Julia code:
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

Structures defined by this package can be used in one of the following
way:

 - converted to an OpenSCAD file;
 - (**TODO**) represented graphically using one of the Julia plotting packages;
 - (**TODO**) directly exported as a 2d (`.svg`) or 3d file (`.stl` or `.ply`).

```@docs
Solids.AbstractSolid
```
