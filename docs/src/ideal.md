# Ideal `Solid` objects
```@meta
CurrentModule = Solid
```

```@docs
Solid.AbstractSolid
```
## Primitive solids

```@docs
Solid.Square
Solid.Circle
Solid.Cube
Solid.Cylinder
Solid.Sphere
Solid.Polygon
Solid.Surface
Solid.NeutralSolid
```

## Transformations

Transformations accept two possible syntaxes:
```julia
    transform(parameters, solid1, solid2, ...)
    transform(parameters) * solid1
```
The second, multiplicative form allows easy chaining of transformations:
```julia
    transform1(param1) * transform2(param2) * solid
```
It may also be applied to several solids by either wrapping them in a
`union`, or equivalently, by applying it to a `Vector` of Solids:
```julia
    transform(parameters) * [ solid1, solid2, ... ]
```

### Affine transformations
```@docs
Solid.mult_matrix
Solid.translate
Solid.scale
Solid.rotate
Solid.mirror
```
TODO: Solid.project

### Extrusion
```@docs
Solid.linear_extrude
Solid.rotate_extrude
Solid.path_extrude
```

### Inserting metadata

A couple of transformations attach metadata to objects.
These are defined using the same base types as affine transforms
and can therefore be applied using the same syntax,
i.e. either as `transform(parameters, s...)`
or as a product `transform(parameters) * s`.

```@docs
Solid.color
Solid.set_parameters
```

### Defining a custom transformation

## CSG operations
```@docs
Solid.union(::Solid.AbstractSolid...)
Solid.intersect
Solid.difference
Solid.hull
Solid.minkowski
```
