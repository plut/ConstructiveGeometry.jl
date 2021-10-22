# Meshing

## Interface

`mesh(object, parameters...)`

The meshing of objects is governed by a few parameters:
 - `accuracy` and `precision` determine the number of faces inserted in the mesh;
 - `symmetry` allows to impose a given rotational symmetry to circles;
 - `type` dictates the coordinate type of the returned mesh (e.g.
	 `Float64` or `Rational{Int}`);
 - `ε` (experimental) is a value of thickness of planes, i.e. any point
	 closer than `ε` from a plane is considered as belonging to the plane.

To set values other than the defaults for an object,
apply the `set_parameters` transform to that object:

```julia
set_parameters(accuracy=1)*
circle(2)
```


## Accuracy and precision

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
while `accuracy` governs the number of vertices for small objects.

### Default values

The default values are
`accuracy = 0.1` and `precision = 1/200`.
The first value means that a circle will deviate by at most 0.1mm from
the ideal circle, and 
the latter value corresponds to the fact
that large circles have 32 sides (see below).

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

``n = \min(π √{r/(2\texttt{accuracy})}, π/ √{\texttt{precision})}).``

In addition, the number of sides is bounded below to always be at least 4.
The number of sides thus increases as the square root of the radius,
with an upper bound.
With the default parameters, one has
``n ≈ \min(7√r, 32)``.

The corresponding value for OpenSCAD is
``n = \min(2πr/\texttt{\textdollar fs},360/\texttt{\textdollar fa})``;
with the default values ``\texttt{\textdollar fa}=12``
and ``\texttt{\textdollar fs=2}``, this gives
``n ≈ \min(π r, 30)``.

### Spheres

Spheres are rendered as [Fibonacci
spheres](http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/).
This produces a more regular mesh than latitude-longitude grids.

```@docs
ConstructiveGeometry.sphere_nvertices
```



## Symmetry

In addition to `accuracy` and `precision`,
the `symmetry` parameter allows forcing the number of vertices
of a circle to be a multiple of a defined value
(by rounding up, if needed, to a multiple of `symmetry`).

# Mesh type for 2d objects

2d objects are represented as the exclusive union (XOR)
of simple-loop polygons.

# Mesh type for 3d objects

3d objects are represented as a triangle mesh,
in a way compatible with LibIGL's functions.
