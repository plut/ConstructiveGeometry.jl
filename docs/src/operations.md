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

