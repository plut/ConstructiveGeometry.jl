# Transformations
```@meta
CurrentModule = ConstructiveGeometry
```
```@setup 0
using ConstructiveGeometry
using Makie
using CairoMakie
png(name, s) = save(name*".png", Makie.plot(s));
```

All single-object transformations accept two possible syntaxes:
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

## Affine transformations
```@docs
mult_matrix
```
```@repl 0
s = mult_matrix([1 0 0;0 1 0;0 .5 1])*cube(10);
png("mult_matrix1", s); # hide
```
![a skewed cube](mult_matrix1.png)
```@docs
translate
```
```@docs
scale
```
```@docs
rotate
```
```@repl 0
s = rotate(30)*square(20)
png("rotate", s); # hide
```
![a rotated square](rotate.png)
```@docs
mirror
```
```@docs
raise
lower
```

## Two-dimensional drawing
```@docs
offset
```
```@repl 0
s = offset(10)*[square(100,50), square(50,100)];
png("offset", s); # hide
```
![example: an offset L-shape](offset.png)
```@docs
opening
```
```@repl 0
s = opening(10)*[square(100,50), square(50,100)];
png("opening", s); # hide
```
![example: the opening of the L-shape](opening.png)
```@docs
closing
```
```@repl 0
s = closing(10)*[square(100,50), square(50,100)];
png("closing", s); # hide
```
![example: the closing of the L-shape](closing.png)

## Extrusion
```@docs
linear_extrude
```
```@docs
rotate_extrude
```
```@docs
path_extrude
```
## Slicing
```@docs
slice
```
```@docs
project
```
```@docs
half_space
```

## Decimation
```@docs
decimate
```

## Inserting metadata

A couple of transformations attach metadata to objects.
These are defined using the same base types as affine transforms
and can therefore be applied using the same syntax,
i.e. either as `transform(parameters, s...)`
or as a product `transform(parameters) * s`.

```@docs
color
set_parameters
```

The `set_parameters` transformation allows attaching arbitrary metadata.
This is on purpose (although there currently exists no easy way
for an user to recover these metadata while meshing an object).
