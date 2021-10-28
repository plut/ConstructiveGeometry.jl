# By order of priority
 - [x] rewrite affine transforms
   either use 3x3+3 matrices internally in all places (ugly)
   or use any types of transform (Julia-esque) and instantiate on meshing
    * decide size for 2d meshes
 - [x] make a nice logo (threaded bolt? some variant of Julia logo?)
 - [x] allow 2d->3d transforms (and take account of these for Minkowski,
   convex hull)
 - [x] Minkowski sum in mixed dimensions
 - [x] allow complement (for intersection)
 - [ ] make meshing type-stable
 - [ ] allow self-union (for fixing meshes)
 - [x] triangulate faces of `Surface`
 - [x] update `iglwrap` and use latest version
   - [x] `loop`
   - [x] `minkowski_sum` (might need tetrahedralize)
   - [ ] `upsample` (only useful once we have deformations)
   - [ ] `swept_volume`
 - [ ] add a parameter to circles, spheres and cylinders to mesh them as
   circumscribed
 - [x] Minkowski sum for holed polygons: slice a connected, holed polygon
   as an almost-simple loop (by connecting outer + inner paths)
   and call binary Minkowski sum
 - [ ] `TriangleMeshes`: have a way to detect non-pwm meshes *and explain
   why*
   - [ ] fix meshes on stl import
 - [ ] also in 2d, regularize polygons (e.g. one backwards loop)
 - [x] libigl contains `offset_surface`
 - [x] `intersect_with_half_space`
 - [ ] add some way of marking individual objects
   - [ ] this requires expanding the mesh types to include marked objects
   - [ ] as well as new syntax, e.g. `!object` or `mark()*object`
 - [x] half-space intersection
  - [ ] clarify parameters for `halfspace`
  - [ ] overload `left_half` etc. for 2d children
 - [x] plane intersection: `slice`
 - [x] projection
 - [ ] define a path type (for Minkowski + stroke) ?
 - [ ] overload `color*object` for `color::Colorant`
 - [x] add examples (with images) in documentation
 - [x] clarify priority: `linear_extrude(8)*(5*object)`
 - [ ] `text`
   - [ ] use `Pango` for text and `FreeType` for fonts
 - [ ] move doc examples to `WGLMakie`
 - [ ] swept surfaces
 - [ ] swung surfaces (`path_extrude`)
 - [ ] wrapped volumes
 - [ ] projection of hollow sphere does not work: replace the temporary
   fix by something better (but this is likely Clipper's fault)
 - [ ] Bézier curves (used as path for `path_extrude`, `stroke`, `polygon`)
 - [ ] get a better syntax for transforms, e.g.
 `symmetrize = transform(axis,s->s ∪ mirror(axis,s))` ??
 - [x] libigl contains `minkowki_sum`
 - [ ] and even `convex_hull`
 - [ ] 2d Minkowski difference
 - [ ] `linear_extrude` with twist and scale
 - [ ] find a way to fix path extrusion? either
   - [ ] cut “by hand” the result of a “butt” extrusion;
   - [ ] intersect the result of a custom “fill” extrusion;
 - [ ] bring back some `Angle` types (with degrees) to allow overloading
   of `extrude`, e.g. `extrude(90°)`; likewise, use `cospi` and `sinpi`
   - [ ] `using Unitful.°` is probably (almost) enough
 - [ ] compute center of gravity (and use it for scaling etc.)
 - [ ] `color`: add more properties (e.g. shininess) to be able to show
   e.g. metal parts
 - [ ] find some way of referencing parts of objects,
 e.g. `cylinder().edge(:top)` references the top edge,
 and have it accessible through CSG hierarchy
 (thus, maybe later providing a way to e.g. fillet it?)
  - [ ] do something for keyword argument dispatch; e.g. the following should be equivalent:
    circle(r=1)
    circle(d=2)
    circle(radius=1)
    circle(diameter=2)
    circle(1)
  and it should be “open” for adding new keywords, e.g.
    circle(1, center=[0,0]) # => translate([0,0])* circle(1)
    square(3, round=1) # => calls rounded_square(3, 1)
  - [ ] path wrapping
 - **Definitions of geometric objects**
 - [ ] Annotations:
    [1,0,0] + annotate("blah", (.5,.5,5))* sphere(3);
    # when meshing, produces
    Annotation("blah",(1.5,.5,.5), mesh(sphere(3)))
 - [ ] import `.stl` and `.ply`
# Basic types
 - [ ] find a way to access `.x`, `.y` and `.z` for `Point` and `Vec`
   types
    [ ] probably requires making `Vec` a separate type from `SArray`
    (to avoid piracy)
 - [x] add a `symmetry` parameter for circles
   - [ ] (and spheres?)
 - [ ] this needs a plane object type (which could be the image by a 2x3
   multmatrix of `:full`), for intersections etc.
 - [ ] think of using `LabelledArrays.jl` (`SLArray`) as an alternative to
   `StaticVector` for `Vec` types
 - [ ] add a 1d type (points, segments; paths) for minkowski (/ extrusions)?
   - [ ] this makes sense; `Clipper.jl` seems happy to do Minkowski with a path
 - [ ] abstract directions (`up`, `left`) etc., interpreted differently
   depending on the dimension.
 - [ ] add a LineNode reference to constructors
   (i.e. first thing in call stack outside module).
# Primitives
 - [ ] add convenience constructors for rounded square, cone, …
 - [ ] simple syntax for making conditionals (⇒ use those empty objects)
   - [ ] or also allow `Nothing` in vectors of objects
   - [ ] even better, `EmptyUnion`
# Syntax
 - [ ] using `:` for transforms is very tempting:
```
   color(red):
   translate([-1,1,1-]):
   square(3)
```
 - [ ] `symmetrize(m, s...) = union(s..., mirror(m, s...))`
   find some easy syntax allowing also `symmetrize(m)*s`
 - [ ] find something like OpenSCAD' # ! % operators.
  - [ ] integer%solid (with various highlight colors)
  - [ ] and find a way to propagate this through the hierarchy
   (a) prefix multiplication by integer constant
   (b) unary operators: + - ! √ ~ ¬
   (c) ad-hoc `Transform` with one-letter name, e.g. `H*square(1)`
 - [ ] think of replacing parameters by kwargs
 - [x]  `∪`, `∩`, `\`: booleans
 - [x] `+ ⊕` Minkowski sum; translation
 - [ ] `object ± real` = offset
 - [ ] `- ⊖` Minkowski difference
 - [ ] `:` convex hull ?
 - [ ] `~` complement
 - [ ] `×` extrusion (scalar => linear_extrude; interval =>
   rotate_extrude; path => path_extrude)
 - [x] `*` multmatrix; scaling
 - [x] Complex * 2d object
 - [ ] make kwargs open for user extension (e.g. rounded squares)
# Transformations
 - [ ] make transformations even lazier, so that they are evaluated only once
   their subjects (and more importantly, their dimension) are known
 - [ ] : overload extrude() (for paths, angles, numbers)
 - [?] a move system (= a representation of abstract affine rotations)
   - [ ] allow `NamedTuple` for this
 - [ ] possible via `move(origin, s...; direction, spin)`
 - [ ] anchor/attachment system
    anchor(square(…), [-1,0])
    anchor(square(…), :left)
    square(…, anchor=:left)
 - [ ] check that it is easy for the user to define arbitrary `Transform`s.
 - [ ] rewrite `attach` using `Transform`
  - [ ] and allow:
    attach(X) * [
      :left => Y, :right => Z,
    ] # as array *or* tuple
 - [ ] 3d Minkowski
 - [ ] replace minkowski with circle by an offset
# Issues in other packages
 - [ ] `Rotations.jl`: using the same type for angles and coordinates is not
   terribly useful (in particular with angles in radians).
# Packaging
 - [?] write a full doc about how to define a new transform
 - [ ] complete the list of exports
 - [x] write a minimal regression test
   - [ ] the doc is the test
   - [ ] add some more complicated examples in the doc to expand the tests
 - [ ] distinguish between core and sub-packages (implementing BOSL2 stuff)?
# Future
 - [?] [https://www.researchgate.net/publication/220184531_Efficient_Clipping_of_Arbitrary_Polygons/link/0912f510a5ac9191e9000000/download]()
 - [x] export to SVG/STL/PLY
 - [ ] splines (and enclosed area)
 - [ ] NURBS
# Extras
 - [ ] deformation https://libigl.github.io/tutorial/#chapter-4-shape-deformation
 - [ ] https://student.cs.uwaterloo.ca/~cs779/Gallery/Winter2018/anietoro/doc/
 - [ ] add a display method that shows the tree
 - [ ] icosphere (from Blender)
 - [ ] sphere made from extruding a half-circle
 - [ ] improve `unit_n_gon` to take advantage of symmetries
 - [ ] Annotations in 2d
 - [ ] Annotations in 3d (this might depend on the visualizer though)
 - [?] rewrite Annotations in terms of `Transform`
 - [x] (more generally, metadata)
 - [ ] add an Annotation type, which passes through all transformations
 - [?] things from BOSL2 to look at:
 - [ ] transforms, distributors, mutators,
 - [x] attachments,
 - [ ] primitives, shapes, shapes2d, masks
 - [ ] math, vectors, arrays, quaternions, affine, coords
geometry, edges, vnf, paths, regions, debug
common, strings, constants, errors,
bezier, threading, rounding, partitions, knurling, skin, hull,
triangulation
polyhedra, screws, metric\_screws

# vim: et:
