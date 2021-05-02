# MAKE mesh() TYPE-STABLE
  - by defining a Mesh{T} pseudo-type (only a constructor)
# Split/rewrite of module
 - **2d subsystem**
  - [x] interface via simple `Vector{SVector{2}}` objects and short list of
    functions
  - [ ] maybe even better, interface as parametric types
  - [ ] rename `PolygonXor`? `PolygonalShape`?
  - [x] add half-plane intersections in this module
  - [ ] reconstruct polygonal shape by list of segments (for projections)
 - **3d subsystem**
  - [x] interface via `CornerTables` and short list of functions
   - [x] document this list of functions
 - **Definitions of geometric objects**
  - [x] _2d primitives_: circle, square, polygon
  - [x] _3d primitives_: sphere, cube, cylinder, surface (?)
  - [x] Some 3d primitives are accessible via extrusion:
   - cube, cylinder, frustum, cone
   - simplicity vs. efficiency?
  - [ ]_Geometry transforms_:
   - [x] invertible affine transform
   - [ ] projection,
   - [ ] plane intersection,
   - [x] linear extrusion, rotational extrusion
  - [x] _Non-geometry transforms_: `set_parameters`, color
  - [ ] _CG operations_:
   - [x] union, inter, diff
   - [ ] hull
   - [ ] offset
   - [ ] minkowski+, minkowski-
 - **Syntactic sugar**
 - **Import/export**
  - SCAD
  - STL
 - **Visualization** (TODO)
 - [ ] **Convex hull** (in its own file)
# For version 0.2 (performance update)
 # corner table http://graphics.stanford.edu/courses/cs468-12-spring/LectureSlides/02_Mesh_Data_Structures.pdf
 - [Aleardi, Devillers] https://hal.inria.fr/inria-00623762/document
 operators: LeftBack, LeftFront, RightBack, RightFront, Source, Target,
 Left, Right; Point, Edge
 - [+] write generic functions for face deletion/insertion in
   HalfEdgeMesh, and use that for (a) mesh construction, (b)
   subtriangulation, (c) simplification (retriangulate faces, removing
   vertices).
 - [+] use a singular half-edge data structure such as:
 https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.11.4605&rep=rep1&type=pdf
 - [x] write set of tests for 2d subsystem: number of vertices and paths
   for known shapes, etc.
 - [x] group together what belongs to `Square`, etc. (in particular,
   include "scad.jl" at the top instead of bottom)
 - [x] use some form of explicit representation for edges
   - [x] may as well use half-edge structure while we're at it
 - [x] document the changes
 - [x] convert Polygon{Int} to Polygon{Float}
 - [?] clean the type system once and for all:
  - [x] suggestion: let `AbstractGeometry` objects have a mess of types
    if they so want,
  - [x] and only decide on a coord. type at meshing type (as a parameter).
  - [ ] take advantage to allow exact (rational) arithmetic
 - [?] split in several packages:
  - [ ] `AbstractGeometry`: remove stale `Meshes` dependency
  - [x] `StrongIndices` -> `StrongArrays`, not used (for now at least)
  - [x] `AABBTree` -> `SpatialSorting` + `TriangleIntersections`
  - [x] `Meshing` -> `CornerTables`
 - [x] implement ZGZJ's cluster triangulation to prevent crossing faces
  - [x] “thicken” faces (and edges) wrt intersection (i.e. add small,
    well-defined tolerance). As a parameter in `set_parameters`?
 - [x] fix difference of two touching solids (e.g. noodle extrusion):
 inclusion (as defined on connected components) does not work in this
 case
   - [x] solution could be in correctly counting opposite faces...
 - [x] linear extrusion of PolygonXor
 - [x] rotate extrusion
 - [ ] intersection with half-plane and half-space
 - [ ] Minkowski sum in 2d
   - [ ] and Minkowski difference
 - [x] clear meshing parameter propagation
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
   - [ ] look in `Makie`
   - in `BasicGeometry`: a list of polygons + list of holes
   - this is simplest (it works as a xor polygon, whereas converting any
     xor to this is harder)
 - [x] reorganize inside the file
 - [ ] 3d => 2d projection
 - [ ] 3d => 2d intersection with plane
 - [x] a common name for meshing objects => **mesh**
 - [x] check using vs. import
 - [x] write a few examples
 - [x] fix `Offset` for polygon xor
 - [x] choose a correct value for `Clipper` precision
 - [x] check that `convex_hull` works
# For version 0.3 (visualization update)
 - [ ] add some per-face visualization data (colors?).
# Performance
 - [ ] try to prevent `TriangleIntersections.inter` from typing
   (e.g. by returning separately the intersection + an array of points,
   or passing the return type of points as a parameter?).
 - [ ] try to use SIMD (e.g. bbox computation)
 - [?] `@inbounds` wherever possible
 - [x] use `sizehint!` for arrays
 - [x] use an AABB tree (bounding box tree) for intersection detection
# Dependencies
 - [ ] `CircularArrays.jl` ?
 - [ ] `Dictionaries.jl` ?
 - [ ] `KeywordDispatch.jl` ?
 - [ ] `Reexport.jl` ?
decide `Meshes.jl`, `GeometryBasics.jl`, or nothing:
 - [x] `Meshes.jl`:
   - [ ] which basic useful algorithms are in this package? not many...
   - [x] seems to work with `Makie`
   - [x] simple base types (`Point` is good; has basic operations + rand)
 - [ ] `GeometryBasics.jl`:
   - [ ] a bit more basic algorithms (meshing of spheres, cylinders?)
     - only very basic stuff
   - [ ] `MeshIO.jl` works with `GeometryBasics`
   - [ ] `Point` is bad
   - [ ] mesh types are awfully long and depend on bad `Point` type
# Basic types
 - [ ] use abstract types wherever possible
 - [ ] find a way to access `.x`, `.y` and `.z` for `Point` and `Vec`
   types
    [ ] probably requires making `Vec` a separate type from `SArray`
    (to avoid piracy)
 - [x] add a `symmetry` parameter for circles
   - [ ] (and spheres?)
 - [ ] clarify `Path`: maybe add `points` iterator and `matrix` abstract
   conversion. Or make `Path` an actual struct type and add accessors.
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
 - [ ] decide whether primitives always have origin=0
   - this is probably simpler for generating points (translating later)
   - and transparent for the user: the constructor can do the translation
 - [x] add trivial types for EmptyUnion and EmptyIntersect
 - [?] decide whether to use Square or square as a name
 suggestion: `Square` is the raw constructor;
 `square` is the convenience user function
  - [ ] (which might return e.g. a rounded square)
  - [ ] also do cylinder, sphere, cube
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
 - [ ] minkowski has a convexity parameter
  - [ ] `convexity`'s place is in `SetParameters`
   - we don't need convexity
 - [ ] Complex * 2d object
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
 - [x] complete the list of exports
 - [x] write a minimal regression test
 - [?] make this a proper package
 - [ ] distinguish between core and sub-packages (implementing BOSL2 stuff)?
# Future
 - [ ] create incidence structure on triangulated surface creation?
   use [directed edges structure](https://core.ac.uk/download/pdf/190807228.pdf)
 - [?] [https://www.researchgate.net/publication/220184531_Efficient_Clipping_of_Arbitrary_Polygons/link/0912f510a5ac9191e9000000/download]()
 - [ ] add some visualization (`Makie`?)
 - [ ] export to SVG/STL/PLY
   - [ ] `MeshIO`
# Extras
 - [ ] [Loop subdivision](https://github.com/cmu462/Scotty3D/wiki/Loop-Subdivision)
 - [ ] https://student.cs.uwaterloo.ca/~cs779/Gallery/Winter2018/anietoro/doc/
 - [ ] add a display method that shows the tree
 - [ ] try 2 strategies for n-ary intersection/union:
   - [ ] merge all structures and compute multiplicity with ray-shooting,
   - [x] or reduce with binary op
 - [ ] icosphere (from Blender)
 - [ ] sphere made from extruding a half-circle
 - [?] [surface decimation](https://caffeineviking.net/papers/stima.pdf)
 - [ ] improve `unit_n_gon` to take advantage of symmetries
 - [x] investigate Fibonacci spheres
 - [x] Color
 - [ ] Annotations in 2d
 - [ ] Annotations in 3d (this might depend on the visualizer though)
 - [?] rewrite Annotations in terms of `Transform`
 - [x] (more generally, metadata)
 - [ ] add an Annotation type, which passes through all transformations
 - [ ] *(Obsolete)*: Offset using OpenSCAD `offset()`
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
