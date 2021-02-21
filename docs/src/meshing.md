# Meshing

## Interface

`Surface(objects...)`


## Accuracy and precision

The meshing of objects is governed by two parameters:
`accuracy` and `precision`.

 - `accuracy` is the maximum absolute deviation allowed when meshing an object.
 This is the maximum distance between the mesh and the ideal shape.
 Its dimensionality is the same as basic length units for the object
 (*i.e.* it will be understood as millimeters by most 3d slicers).

 - `precision` is the maximum relative deviation allowed when meshing.
 This is a dimensionless number.

When meshing an object, the minimum value will be used
between those given by these two definitions.
This means that `precision` gives an absolute maximum
on the number of vertices for large objects,
while `accuracy` governs 

### Default values

The default values are
`accuracy = 0.2` and `precision = 1/200`.
The latter value corresponds to the fact
that large circles have 32 sides (see below).

### Modifying the values

To set values other than the defaults for an object,
apply the `set_parameters` transform to that object:

```julia
set_parameters(accuracy=0.2)*
Circle(2)
```

### Circles

A circle of radius ``r`` is replaced by an inscribed ``n``-gon.
The deviation between the ideal circle and the ``n``-gon
is the sagitta of the [circular
segment](https://en.wikipedia.org/wiki/Circular_segment)
with radius ``r`` and central angle ``2π/n``;
its value is hence ``s = r(1-\cos(π/n)) ≈ \frac{π^2 r}{2 n^2}``.

By definition, ``\texttt{accuracy} = s``
while ``\texttt{precision} = s/r \approx \frac{π^2}{2 n^2}``.
This gives

``n = \min(π √{r/(2\texttt{accuracy})}, π √{1/(2\texttt{precision})}).``

In addition, the number of sides is bounded below to always be at least 4.
The number of sides thus increases as the square root of the radius,
with an upper bound. With the default parameters, that upper bound is
``n=32``.
