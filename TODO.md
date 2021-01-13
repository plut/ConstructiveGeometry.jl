# Immediate work
 - distinguish ideal solid and elements
 + replace include() by a function that reads a `Solid` from a `JuliaFile`
 - Minkowski difference
 - what to do for polygons with holes?
 - replace minkowski with circle by an offset
 + implement polygon convolution
 - finish grouping all Clipper stuff in one section
 - fix `Offset` for `Region` values
  + add a NamedTuple to convert symbols to ClipperJoinType
 - choose a correct value for `Clipper` precision
 * check `convex_hull`
# Basic types
 - add something for fake-3d objects (embedded in a subspace):
   this would represent both `mult_matrix` with zero determinant,
   `mult_matrix` of 2d object with 3d matrix,
   and `project()` (or `cut()`).
   It would also allow, say, convex hull with a translated 2d object.
 - this needs a plane object type (which could be the image by a 2x3
   multmatrix of `:full`).
 - think of using `LabelledArrays.jl` (`SLArray`) as an alternative to
   `StaticVector` for `Vec` types
 - add a 1d type (points, segments; paths) for minkowski (/ extrusions)?
   - this makes sense; `Clipper.jl` seems happy to do Minkowski with a path
   *no*, feature request written
 - abstract directions (`up`, `left`) etc., interpreted differently
   depending on the dimension.
# Primitives
 + decide whether to use Square or square as a name
 suggestion: `Square` is the raw constructor;
 `square` is the convenience user function
  - (which might return e.g. a rounded square)
  - also do cylinder, sphere, cube
 - add convenience constructors for rounded square, cone, …
 - simple syntax for making conditionals (⇒ use those empty objects)
   - or also allow `Nothing` in vectors of objects
 - import `.stl` and `.ply`
# Transformations
 + make `hull`, `minkowski` formal and `Region`-able
 * a move system (= a representation of abstract affine rotations)
   - allow `NamedTuple` for this
 - possible via `move(origin, s...; direction, spin)`
 + anchor/attachment system
    anchor(square(…), [-1,0])
    anchor(square(…), :left)
    square(…, anchor=:left)
 - make difference() a strictly binary operation?
 + add a reduce() operator that multiplies all the matrices
 - check that it is easy for the user to define arbitrary `Transform`s.
 - rewrite `attach` using `Transform`
  - and allow:
    attach(X) * [
      :left => Y, :right => Z,
    ] # as array *or* tuple
# Issues in other packages
 + `HybridArrays.jl`: fix the type of static * hybrid matrix
   multiplication (currently fully dynamic...)
 - `StaticArrays.jl`: SDiagonal is currently *not* a static matrix?
    julia> SDiagonal(1,2,3) isa StaticMatrix
    false
 + `Clipper.jl`: IntPoint should really be a `StaticVector{2,Fixed}`
   (even better: an abstract type, or a labelled vector...).
 + `Clipper.jl`: Minkowski sum
 - `Rotations.jl`: using the same type for angles and coordinates is not
   terribly useful (in particular with angles in radians).
 - *Julia*: add `cossin` to `sincos` (helps with complex units).
# Syntax
 - think of replacing parameters by kwargs
 +  ∪, ∩, \
 - `+ ⊕` Minkowski sum; translation
 - `- ⊖` Minkowski difference
 - `:` hull ?
 - `×` linear_extrude, rotate_extrude
 + `*` multmatrix; scaling
 - think of overloading `{...}` or`[...]` (either `hcat` or `vcat`,
   and `braces`, `bracescat`).
    dump(:({a b})) => :bracescat
   **no**: will not work (but could in a macro...)
 - really really stupid idea: *n*-dimensional matrix actually arranges
   objects in a matrix...
# Transformations
 * call Clipper to provide offset algorithm
	- orientation, area, pointinpolygon
	- this provides polygon intersection, difference, …
	- also offset and `get_bounds`
 + draw(path, width) (using Clipper.offset)
 + convex hull (in dim 2)
 + convex hull (in dim 3)
 - convex hull in mixed dimensions makes sense
    (also: image of 2d object in another plane?).
 - : overload extrude() (for paths, angles, numbers)
 - minkowski has a convexity parameter
  - `convexity`'s place is in `SetParameters`
# Packaging
 - write a minimal regression test
 * make this a proper package
 - distinguish between core and sub-packages (implementing BOSL2 stuff)
# Future
[https://www.researchgate.net/publication/220184531_Efficient_Clipping_of_Arbitrary_Polygons/link/0912f510a5ac9191e9000000/download]()
 - Minkowski with Clipper.jl
 - interface with CGAL for computing CSG in 3d
 - add some visualization (`Makie`?)
 - export to SVG/STL/PLY
# Extras
 - improve `unit_n_gon` to take advantage of symmetries
 + Color
 - Annotations in 2d
 - Annotations in 3d (this might depend on the visualizer though)
 * rewrite Annotations in terms of `Transform`
 + (more generally, metadata)
 - add an Annotation type, which passes through all transformations
 - *(Obsolete)*: Offset using OpenSCAD `offset()`
 * things from BOSL2 to look at:
 - transforms, distributors, mutators,
 + attachments,
 - primitives, shapes, shapes2d, masks
 - math, vectors, arrays, quaternions, affine, coords
geometry, edges, vnf, paths, regions, debug
common, strings, constants, errors,
bezier, threading, rounding, partitions, knurling, skin, hull,
triangulation
polyhedra, screws, metric\_screws
