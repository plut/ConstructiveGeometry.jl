# Quick-start

## Basic example
```julia
using Solid
import Solid: Square, Circle, mult_matrix, translate, scale, color

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
Solid.include
Solid.scad
```

