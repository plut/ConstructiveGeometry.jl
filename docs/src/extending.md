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

## The `mesh` method

The main work of converting ideal objects to concrete
meshes is done by the following `mesh` method:

    mesh(opt::MeshOptions, object, children_meshes)

This method should in principle not be directly called by the user.
It is invoked by the `ConstructiveGeometry` module with the following
parameters:

 - `opt` is a `MeshOptions{T}` structure holding
   the [parameters for computing the mesh](@ref Meshing-parameters)
   (the type `T` is the coordinate type used for mesh vertices);
 - `object` is the object being meshed itself;
 - `children_meshes` is a list (vector or tuple) of meshes computed
   for the contents of `children(object)`.

This method must return a `FullMesh` object, which in turns contains:
 - a main mesh: this is the mesh which will be used when outputting e.g. a STL file;
 - a vector of auxiliary meshes, representing any [highlighted](@ref Highlight) parts of the object.

The main mesh can be either a `ShapeMesh` (i.e. an exclusive union of
closed polygonal loops) or a `VolumeMesh` (i.e. the volume delimited by a
closed triangulated surface). These types are also the types
returned by the [`polygon()`](@ref Polygon) and [`surface()`](@ref Surface)
primitives.

Any two-dimensional object (viz. three-dimensional) object `x`
may be instantiated as a mesh by invoking `polygon(x)` (viz. `surface(x)`).

## Defining a new primitive object

Defining a new primitive object has two main parts:

1. define a new concrete subtype of `AbstractGeometryLeaf{D}`
   (where `D` = 2 or 3 is the dimension of the object)
   holding the geometric information for the “ideal” object;

2. implement a meshing method for this type, by extending
   the `ConstructiveGeometry.mesh` method; more specifically,
   `mesh(opt::MeshOptions{T}, s::NewObjectType, _)`.
   (the last parameter is the list of children meshes,
   which will always be empty for a leaf type).

On top of this, several methods may be added as syntactic sugar
to help the user define objects of this new type.

## Defining a new transformation

1. If this transformation applies to a single object,
   then define a new concrete subtype of `AbstractTransform{D}`
   with a single `child` field;

2. otherwise, define a new concrete subtype of `AbstractGeometry{D}`
   and extend the `children` method to the new type,
   returning a proper list (vector, tuple etc; it must only be iterable)
   of the object's children.

3. extend the `mesh` method to this type.

As in the case of leaf objects,
it is also possible to add some methods to make it easier for the user
to call the transformation, e.g. by using the `operator()` function.

### The `operator()` function

This function is a shortcut to define the multiplicative syntax
for an object transformation. It has two main methods:
 - `operator(Constructor, parameters, solid)` (where `parameters` is a tuple)
   is the same as `Constructor(parameters..., solid)`;
 - `operator(Constructor, parameters)` returns an encapsulation
   such that `operator(Constructor, parameters) * solid`
   is again interpreted as `Constructor(parameters..., solid)`.

For example, multiplicative syntax can be given to a transformation
`Frobnicate` by defining
`frobnicate(param1, param2, s...) = operator(Frobnicate, (param1, param2,), s...)`.
