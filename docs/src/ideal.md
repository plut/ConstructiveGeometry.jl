# Ideal objects
```@meta
CurrentModule = ConstructiveGeometry
```

## Primitive solids

### 2d
```@docs
ConstructiveGeometry.square
ConstructiveGeometry.circle
ConstructiveGeometry.stroke
ConstructiveGeometry.polygon
```

### 3d
```@docs
ConstructiveGeometry.cube
ConstructiveGeometry.cone
ConstructiveGeometry.cylinder
ConstructiveGeometry.sphere
ConstructiveGeometry.surface
```

## Transformations

All transformations accept two possible syntaxes:
```julia
    transform(parameters, solid1, solid2, ...)
    transform(parameters) * solid1
```
The second, multiplicative form allows easy chaining of transformations:
```julia
    transform1(param1) * transform2(param2) * solid
```
This form may also be applied to several solids by either wrapping them in a
`union`, or equivalently, by applying it to a `Vector` of ConstructiveGeometry:
```julia
    transform(parameters) * [ solid1, solid2, ... ]
```

### Affine transformations
```@docs
ConstructiveGeometry.mult_matrix
ConstructiveGeometry.translate
ConstructiveGeometry.scale
ConstructiveGeometry.rotate
ConstructiveGeometry.mirror
ConstructiveGeometry.raise
ConstructiveGeometry.lower
```

### 2d drawing
```@docs
ConstructiveGeometry.offset
ConstructiveGeometry.opening
ConstructiveGeometry.closing
```
## Operations

### CSG operations
```@docs
union(::AbstractGeometry,::AbstractGeometry)
intersect(::AbstractGeometry{D},::AbstractGeometry{D}) where{D}
setdiff(::AbstractGeometry{D},::AbstractGeometry{D}) where{D}
ConstructiveGeometry.hull
ConstructiveGeometry.minkowski
```

### Extrusion
```@docs
ConstructiveGeometry.linear_extrude
ConstructiveGeometry.rotate_extrude
ConstructiveGeometry.path_extrude
```
### Slicing
```@docs
ConstructiveGeometry.slice
ConstructiveGeometry.project
ConstructiveGeometry.half_space
```

### Decimation
```@docs
decimate
```

### Inserting metadata

A couple of transformations attach metadata to objects.
These are defined using the same base types as affine transforms
and can therefore be applied using the same syntax,
i.e. either as `transform(parameters, s...)`
or as a product `transform(parameters) * s`.

```@docs
ConstructiveGeometry.color
ConstructiveGeometry.set_parameters
```

The `set_parameters` transformation allows attaching arbitrary metadata.
This is on purpose (although there currently exists no easy way
for an user to recover these metadata while meshing an object).

### Defining a custom transformation

