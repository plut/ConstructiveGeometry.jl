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

