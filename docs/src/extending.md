# [Extending `ConstructiveGeometry`](@id extending)
```@meta
CurrentModule = ConstructiveGeometry
```
```@setup 0
using ConstructiveGeometry
using Makie
using CairoMakie
png(name, s) = save(name*".png", Makie.plot(s));
```

## Defining a new primitive object

Defining a new primitive object has two main parts:

1. define a new concrete subtype of `AbstractGeometry{D}`
(where `D` = 2 or 3 is the dimension of the object);

2. implement a meshing method for this type, by extending
the `ConstructiveGeometry.mesh` method.

## Defining a new transformation

1. If this transformation applies to a single object,
then define a new concrete subtype of `AbstractTransform{D}`
with a single `child` field;

2. otherwise, define a new concrete subtype of `AbstractGeometry{D}`
and extend the `children` method to the new type,
returning a proper list (vector, tuple etc; it must only be iterable)
of the object's children.

3. extend the `mesh` method to this type.

4. (optional) add some syntactic sugar to make it easier for the user
to call the transformation, e.g. by using the `operator()` function.
