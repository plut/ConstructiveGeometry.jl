`rotate(90)*cube`
# Misc.
  - [ ] find a way to fix path extrusion? either
    - [ ] cut “by hand” the result of a “butt” extrusion;
    - [ ] intersect the result of a custom “fill” extrusion;
  - [ ] do something for keyword argument dispatch; e.g. the following should be equivalent:
    circle(r=1)
    circle(d=2)
    circle(radius=1)
    circle(diameter=2)
    circle(1)
  and it should be “open” for adding new keywords, e.g.
    circle(1, center=[0,0]) # => translate([0,0])* circle(1)
    square(3, round=1) # => calls rounded_square(3, 1)
 - **2d subsystem**
  - [x] reconstruct polygonal shape by list of segments (for projections)
 - [ ] Minkowski sum in 2d
   - [ ] and Minkowski difference
 - **3d subsystem**
  - [ ] (re-)add half-space and
  - [ ] plane intersections
 - **Definitions of geometric objects**
  - [x] Some 3d primitives are accessible via extrusion:
   - cube, cylinder, frustum, cone
   - simplicity vs. efficiency?
  - [x]_Geometry transforms_:
   - [x] invertible affine transform
   - [x] projection,
   - [x] plane intersection,
   - [x] linear extrusion, rotational extrusion
   - [x] path extrusion
   - [ ] path wrapping
   - [ ] rename `Path` to `Stroke`
   todo: enlarge parameters for rotate_Extrude
  - [x] _Non-geometry transforms_: `set_parameters`, color
  - [ ] _CG operations_:
   - [x] union, inter, diff
   - [x] hull
   - [ ] 3d offset
   - [ ] minkowski+, minkowski-
 - **Syntactic sugar**
 - **Import/export**
  - [x] SCAD
  - [x] STL
 - **Visualization** (TODO)
 - [x] **Convex hull** (in its own file)
# Visualization
 - [x] Annotations:
    [1,0,0] + annotate("blah", (.5,.5,5))* sphere(3);
    # when meshing, produces
    Annotation("blah",(1.5,.5,.5), mesh(sphere(3)))

 - [x] 3d visualization with `Makie` `mesh`
   - [x] add some per-face visualization data (colors?).
 - [ ] 2d with `polygon`
# For version 0.2 (performance update)
 - [x] what to do for polygons with holes? find a representation that
   must be useable for extrusion + (makie) drawing + Clipper + openscad conversion
   possibilities include:

|           | U(poly with hole) | (polys) + (holes) | xor list         |
|-----------|-------------------|-------------------|------------------|
|Clipper in | easy (as xor)     | easy (as xor)     | easy             |
|Clipper out| **hard**          |                   | easy             |
|extrude    | easy              |                   |xor of extrusions |
|minkowski  | easy-ish          | ?                 |quite hard (split)|
|Makie      | easy              |                   |**hard**?         |
|svg        | easy              |                   |easy (`evenodd`)  |
|openscad   | easy              | just works        | just works       |

Convolution: parse polys + holes; add to polys, subtract from holes.
=> must implement Minkowski difference first!
http://acg.cs.tau.ac.il/copy_of_projects/in-house-projects/exact-and-efficient-construction-of-planar/mink-conv.pdf

Clipper out: must parse xor list as a list of polys and holes.
Same hardness for ⋃(poly+hole) and (polys)+(holes).

Extrude of ⋃(p+h): triangulate faces and build manually.

    1. union of polygon with holes: trivial to convert
    (2) list of polys + list of holes,
    (3) flat xor of polygons]
   - [ ] xor is doable for extrusions but might be harder to convert to
     openscad
   - in `BasicGeometry`: a list of polygons + list of holes
   - this is simplest (it works as a xor polygon, whereas converting any
     xor to this is harder)
# Performance
 - [ ] try to use SIMD (e.g. bbox computation)
 - [?] `@inbounds` wherever possible
# Dependencies
 - [ ] `CircularArrays.jl` ?
 - [ ] `Dictionaries.jl` ?
 - [ ] `KeywordDispatch.jl` ?
 - [ ] `Reexport.jl` ?
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
# 2d vs 3d
 - [ ] add something for fake-3d objects (embedded in a subspace):
   this would represent both `mult_matrix` with zero determinant,
   `mult_matrix` of 2d object with 3d matrix,
   and `project()` (or `cut()`).
   It would also allow, say, convex hull with a translated 2d object.
 - [ ] before deciding what *“should”* be done, write a set of examples of
   **what it should do** (for various operations, e.g. linear maps,
   Minkowski, hull) and then decide how to best implement this behavior
 - [ ] Objects should really have *two* dimensions: intrinsic and embedding.
 E.g. a `square(1)` has dimensions `(2,2)`, while its translation by
 `[0,0,1]` has dimensions `(2,3)` and an embedding given by the
 corresponding matrix.
 - [ ] CSG operations may be performed either
   - [?] on objects of the same dimension, same embedding
   - [?] `hull`: use embedding to push all objects to same space if possible
   - [?] `minkowski`: ditto
# Primitives
 - [x] decide whether primitives always have origin=0
   - [x] this is probably simpler for generating points (translating later)
   - [x] and transparent for the user: the constructor can do the translation
 - [x] add trivial types for EmptyUnion and EmptyIntersect
 - [x] decide whether to use Square or square as a name
 suggestion: `Square` is the raw constructor;
 `square` is the convenience user function
  - [x] (which might return e.g. a rounded square)
  - [x] also do cylinder, sphere, cube
 - [ ] add convenience constructors for rounded square, cone, …
 - [ ] simple syntax for making conditionals (⇒ use those empty objects)
   - [ ] or also allow `Nothing` in vectors of objects
 - [ ] import `.stl` and `.ply`
 - [ ] STL export
# Syntax
 - [ ] using `:` for transforms is very tempting:
```
   color(red):
   translate([-1,1,1-]):
   square(3)
```
 - [ ] find something like OpenSCAD' # ! % operators.
  - [ ] integer%solid (with various highlight colors)
  - [ ] and find a way to propagate this through the hierarchy
   (a) prefix multiplication by integer constant
   (b) unary operators: + - ! √ ~ ¬
   (c) ad-hoc `Transform` with one-letter name, e.g. `H*square(1)`
 - [ ] think of replacing parameters by kwargs
 - [x]  `∪`, `∩`, `\`: booleans
 - [ ] `+ ⊕` Minkowski sum; translation
 - [ ] `object ± real` = offset
 - [ ] `- ⊖` Minkowski difference
 - [ ] `:` convex hull ?
 - [ ] `~` complement
 - [ ] `×` extrusion (scalar => linear_extrude; interval =>
   rotate_extrude; path => path_extrude)
 - [x] `*` multmatrix; scaling
 - [ ] Complex * 2d object
 - [ ] make kwargs open for user extension (e.g. rounded squares)
# Transformations
 - [ ] make transformations even lazier, so that they are evaluated only once
   their subjects (and more importantly, their dimension) are known
 - [x] call Clipper to provide offset algorithm
  + orientation, area, pointinpolygon
  + this provides polygon intersection, difference, …
  + also offset and `get_bounds`
 - [x] draw(path, width) (using Clipper.offset)
 - [x] convex hull (in dim 2)
 - [x] convex hull (in dim 3)
 - [ ] convex hull in mixed dimensions makes sense
    (also: image of 2d object in another plane?).
 - [ ] : overload extrude() (for paths, angles, numbers)
 - [?] a move system (= a representation of abstract affine rotations)
   - [ ] allow `NamedTuple` for this
 - [ ] possible via `move(origin, s...; direction, spin)`
 - [x] anchor/attachment system
    anchor(square(…), [-1,0])
    anchor(square(…), :left)
    square(…, anchor=:left)
 - [x] make difference() a strictly binary operation?
 - [x] add a reduce() operator that multiplies all the matrices
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
 - [ ] write another package with `BinaryBuilder.jl` for IGL wrapper(s)
 - [x] complete the list of exports
 - [x] write a minimal regression test
 - [x] make this a proper package
 - [ ] distinguish between core and sub-packages (implementing BOSL2 stuff)?
# Future
 - [?] [https://www.researchgate.net/publication/220184531_Efficient_Clipping_of_Arbitrary_Polygons/link/0912f510a5ac9191e9000000/download]()
 - [ ] add some visualization (`Makie`?)
 - [ ] export to SVG/STL/PLY
   - [ ] `MeshIO`
# Extras
 - [?] [surface decimation](https://caffeineviking.net/papers/stima.pdf)
 - [ ] decimation (`igl::decimate`)
 - [ ] [Loop subdivision](https://github.com/cmu462/Scotty3D/wiki/Loop-Subdivision)
   - [ ] `igl::upsample`, `igl::loop`
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
