# Immediate work
 - a common name for meshing objects (either `Mesh` or `elements` ?)
   => **realize**
 - Boolean operations work even if non-connected triangulation
   -> remove conn.comp. split from `Triangulation()` converter
 - add complement (of triangulation, and symbolic)
 - add difference (of triangulations)
 * intersection and union: ray-shooting (it is possible even within the
	 merged structure)
   - force-remove opposite faces (neither ∪ nor ∩)
	 - check if it works for complement
 + use `import` for modules used only a few times once (`Color`, all geometry)
   to avoid polluting the namespace
 * test suite
 - what to do for polygons with holes?
   - look in `Makie`
	 - in `BasicGeometry`: a list of polygons + list of holes
	 - parity is probably simplest bet (easily allows both non-connected
		 polys and holes)
 - replace minkowski with circle by an offset
 - Minkowski difference
 + finish grouping all Clipper stuff in one section
 - fix `Offset` for polygon unions
 - choose a correct value for `Clipper` precision
 * check `convex_hull` works
# Basic types
 - clarify `Path`: maybe add `points` iterator and `matrix` abstract
	 conversion. Or make `Path` an actual struct type and add accessors.
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
 - abstract directions (`up`, `left`) etc., interpreted differently
   depending on the dimension.
 - add a LineNode reference to constructors
   (i.e. first thing in call stack outside module).
# 2d vs 3d
 - before deciding what *“should”* be done, write a set of examples of
	 **how it should work** (for various operations, e.g. linear maps,
	 Minkowski, hull) and then decide how to best implement this behavior
 - Objects should really have *two* dimensions: intrinsic and embedding.
 E.g. a `square(1)` has dimensions `(2,2)`, while its translation by
 `[0,0,1]` has dimensions `(2,3)` and an embedding given by the
 corresponding matrix.
 - CSG operations may be performed either
   * on objects of the same dimension, same embedding
   * `hull`: use embedding to push all objects to same space if possible
   * `minkowski`: ditto
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
 - STL export
# Syntax
 - using `:` for transforms is very tempting:
```
   color(red):
	 translate([-1,1,1-]):
	 square(3)
```
 - find something like OpenSCAD' # ! % operators.
   (a) prefix multiplication by integer constant
	 (b) unary operators: + - ! √ ~ ¬
	 (c) ad-hoc `Transform` with one-letter name, e.g. `H*square(1)`
 - think of replacing parameters by kwargs
 +  `∪`, `∩`, `\`
 - `+ ⊕` Minkowski sum; translation
 - `- ⊖` Minkowski difference
 - `:` hull ?
 - `¬` complement
 - `×` linear_extrude, rotate_extrude
 + `*` multmatrix; scaling
 - think of overloading `{...}` or`[...]` (either `hcat` or `vcat`,
   and `braces`, `bracescat`).
    dump(:({a b})) => :bracescat
   **no**: will not work (but could in a macro...)
 - really really stupid idea: *n*-dimensional matrix actually arranges
   objects in a matrix...
# Transformations
 - make transformations even lazier, so that they are evaluated only once
   their subjects (and more importantly, their dimension) are known
 + call Clipper to provide offset algorithm
	+ orientation, area, pointinpolygon
	- this provides polygon intersection, difference, …
	+ also offset and `get_bounds`
 + draw(path, width) (using Clipper.offset)
 + convex hull (in dim 2)
 + convex hull (in dim 3)
 - convex hull in mixed dimensions makes sense
    (also: image of 2d object in another plane?).
 - : overload extrude() (for paths, angles, numbers)
 - minkowski has a convexity parameter
  - `convexity`'s place is in `SetParameters`
 - Complex * 2d object
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
 - 3d Minkowski
# Issues in other packages
 - `StaticArrays.jl`: SDiagonal is currently *not* a static matrix?
    julia> SDiagonal(1,2,3) isa StaticMatrix
    false
 - `Rotations.jl`: using the same type for angles and coordinates is not
   terribly useful (in particular with angles in radians).
 - *Julia*: add `cossin` to `sincos` (helps with complex units).
# Packaging
 - write a full doc about how to define a new transform
 - complete the list of exports
 * write a minimal regression test
 * make this a proper package
 - distinguish between core and sub-packages (implementing BOSL2 stuff)?
# Future
 * [https://www.researchgate.net/publication/220184531_Efficient_Clipping_of_Arbitrary_Polygons/link/0912f510a5ac9191e9000000/download]()
 - add some visualization (`Makie`?)
 - export to SVG/STL/PLY
   - `MeshIO`
# Extras
 - improve `unit_n_gon` to take advantage of symmetries
 - investigate Fibonacci spheres
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
