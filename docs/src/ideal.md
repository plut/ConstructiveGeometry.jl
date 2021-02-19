# Ideal objects
```@meta
CurrentModule = ConstructiveGeometry
```

## Primitive solids

```@docs
ConstructiveGeometry.Square
ConstructiveGeometry.Circle
ConstructiveGeometry.Cube
ConstructiveGeometry.Cylinder
ConstructiveGeometry.Sphere
ConstructiveGeometry.Polygon
ConstructiveGeometry.Surface
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
```

TODO: `ConstructiveGeometry.project`, `ConstructiveGeometry.cut`.

### Extrusion
```@docs
ConstructiveGeometry.linear_extrude
ConstructiveGeometry.rotate_extrude
ConstructiveGeometry.path_extrude
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

### Defining a custom transformation

## Operations
```@docs
union(::AbstractGeometry,::AbstractGeometry)
intersect(::AbstractGeometry,::AbstractGeometry)
ConstructiveGeometry.difference
ConstructiveGeometry.hull
ConstructiveGeometry.minkowski
```
