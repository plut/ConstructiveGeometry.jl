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

A circle of radius ``r`` is approximated by an inscribed ``n``-gon.
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
This produces a more regular mesh than latitude-longitude grids
(in particular, the grid does not have singularities at the poles).


A sphere is approximated by an inscribed polyhedron with ``n`` vertices.
Such a polyhedron has ``2n-4`` triangular faces;
the average area of a face is ``\frac{4π r^2}{2n-4} = \frac{2π r^2}{n-2}``,
thus the average (squared) edge length is
``d² ≈ (8π/√3) \frac{r^2}{n-2}``
(according to the unit equilateral triangle area ``√3/4``).

The sagitta for a chord of length ``d`` is given by
``s/r = 1 - √{1-d^2/4r^2} ≈ (1-(1-d^2/8 r^2)) ≈ (π/√3)/(n-2)``.
Hence we find
``n ≈ 2 + (π/√3)/(\textt{max}(\texttt{precision},\textt{accuracy}/r))``.

With the default values for `accuracy` and `precision`:
 - small spheres have approximately ``2+18r`` vertices (and always at least 6 vertices);
 - large spheres have 365 vertices.


## Symmetry

In addition to `accuracy` and `precision`,
the `symmetry` parameter allows forcing the number of vertices
of a circle to be a multiple of a defined value
(by rounding up, if needed, to a multiple of `symmetry`).

# Mesh types and orientation

## Two-dimensional shapes

Two-dimensional objects are represented as the exclusive union (XOR)
of simple-loop polygons (using even-odd rule).
Internally, direct loops (counter-clockwise) represent polygons,
and retrograde loops (clockwise) represent holes.

## Three-dimensional volumes

Three-dimensional objects are represented as a triangle mesh,
in a way compatible with LibIGL's functions.
The triangles are oriented so that their normal vector points outside the
volume of the object.
