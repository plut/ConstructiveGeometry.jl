# CSG operations
```@meta
CurrentModule = ConstructiveGeometry
```
```@setup 0
using ConstructiveGeometry
using Makie
using CairoMakie
png(name, s) = save(name*".png", Makie.plot(s));
```

## CSG operations
```@docs
union(::AbstractGeometry,::AbstractGeometry)
intersect(::AbstractGeometry{D},::AbstractGeometry{D}) where{D}
setdiff(::AbstractGeometry{D},::AbstractGeometry{D}) where{D}
hull
minkowski
```

