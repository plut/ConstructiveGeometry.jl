# Ideal `Solids` objects
```@meta
CurrentModule = Solids
```

```@docs
Solids.AbstractSolid
```
## Primitive solids

```@docs
Solids.Square
Solids.Circle
Solids.Cube
Solids.Cylinder
Solids.Sphere
Solids.Polygon
Solids.Surface
Solids.NeutralSolid
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
`union`, or equivalently, by applying it to a `Vector` of Solids:
```julia
    transform(parameters) * [ solid1, solid2, ... ]
```

### Affine transformations
```@docs
Solids.mult_matrix
Solids.translate
Solids.scale
Solids.rotate
Solids.mirror
```

TODO: `Solids.project`, `Solids.cut`.

### Extrusion
```@docs
Solids.linear_extrude
Solids.rotate_extrude
Solids.path_extrude
```

### Inserting metadata

A couple of transformations attach metadata to objects.
These are defined using the same base types as affine transforms
and can therefore be applied using the same syntax,
i.e. either as `transform(parameters, s...)`
or as a product `transform(parameters) * s`.

```@docs
Solids.color
Solids.set_parameters
```

### Defining a custom transformation

## Operations
```@docs
Solids.union
Solids.intersect
Solids.difference
Solids.hull
Solids.minkowski
```
