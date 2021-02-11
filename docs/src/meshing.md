# Meshing

## Interface

`Surface(objects...)`


## Accuracy and precision

 - `accuracy` is the maximum absolute deviation allowed when meshing an object.
The default value is `2.0`, corresponding to the default value in
OpenSCAD for `$fs`.

 - `precision` is the maximum relative deviation allowed when meshing.
The default value is `0.02`, corresponding to the default value of `$fa`
in OpenSCAD (roughly, `1-cos(180°/$fa)`).

To set values other than the defaults for an object,
apply the `set_parameters` transform to that object:

```julia
set_parameters(accuracy=0.2)*
Circle(2)
```


### Circles

A circle of radius ``r`` is replaced by an inscribed ``n``-gon,
where ``n`` is determined such that:

 - each side has a length ≈``2π r/n``, which must not be smaller than
	 `accuracy`
	 (hence ``n \leq 2πr/\texttt{accuracy}``);
 - deviation from the ideal circle (sagitta) is ``1/cos(π/n)``,
	which must not be smaller than `precision`
	(hence ``n \leq π/\sqrt{2\texttt{precision}}``).

In all cases, the number of sides is at least 4.
