# Immediate work
 - [x] clear meshing parameter propagation
 - [-] `Region`: finish 2d subsystem
 - [x] what to do for polygons with holes? find a representation that
	 must be useable for extrusion + (makie) drawing + Clipper + openscad conversion
	 possibilities include:

|           | U(poly with hole) | (polys) + (holes) | xor list         |
|-----------|-------------------|-------------------|------------------|
|Clipper in | easy (as xor)     | easy (as xor)     | easy             |
|Clipper out| **hard**          |                   | easy             |
|extrude    | easy              |                   |xor of extrusions |
|minkowski  | easy-ish          | ?                 |quite hard (split)|
|Makie      | easy              |                   |**hard**          |
|svg        | easy              |                   |easy (`evenodd`)  |
|openscad   | easy              | just works        | just works       |

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
 - [?] reorganize inside the file
 - [ ] split in several files after all (will help in `include()` time)?
 - [x] a common name for meshing objects
   => **mesh**
 - [x] check using vs. import
 - [ ] replace `triangulate_between` by an actual Triangulate call
 - [x] write a few examples
 - [x] test suite
 - [ ] Minkowski difference of polygons
 - [ ] fix `Offset` for polygon xor
 - [ ] choose a correct value for `Clipper` precision
   - [ ] using the `sides` function...
 - [?] check that `convex_hull` works
# Dependencies
 - [ ] use `Dictionaries.jl`
 - [ ] use `Reexport.jl`
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
 - [ ] use an AABB tree (bounding box tree) for intersection detection
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
   (a) prefix multiplication by integer constant
	 (b) unary operators: + - ! √ ~ ¬
	 (c) ad-hoc `Transform` with one-letter name, e.g. `H*square(1)`
 - [ ] think of replacing parameters by kwargs
 - [x]  `∪`, `∩`, `\`
 - [ ] `+ ⊕` Minkowski sum; translation
 - [ ] `- ⊖` Minkowski difference
 - [ ] `:` hull ?
 - [ ] `¬` complement
 - [ ] `×` linear_extrude, rotate_extrude
 - [x] `*` multmatrix; scaling
 - [ ] think of overloading `{...}` or`[...]` (either `hcat` or `vcat`,
   and `braces`, `bracescat`).
    dump(:({a b})) => :bracescat
   - [?]*no**: will not work (but could in a macro...)
 - [ ] really really stupid idea: *n*-dimensional matrix actually arranges
   objects in a matrix...
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
 - [ ] Complex * 2d object
 - [?] a move system (= a representation of abstract affine rotations)
   - [ ] allow `NamedTuple` for this
 - [ ] possible via `move(origin, s...; direction, spin)`
 - [x] anchor/attachment system
    anchor(square(…), [-1,0])
    anchor(square(…), :left)
    square(…, anchor=:left)
 - [ ] make difference() a strictly binary operation?
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
 - [ ] complete the list of exports
 - [?] write a minimal regression test
 - [?] make this a proper package
 - [ ] distinguish between core and sub-packages (implementing BOSL2 stuff)?
# Future
 - [ ] create incidence structure on triangulated surface creation?
   use [directed edges structure](https://core.ac.uk/download/pdf/190807228.pdf)
 - [ ] maybe using ~ instead of `-` for reversed faces
 - [ ] and use AABBtree for computing face intersections
 - [?] [https://www.researchgate.net/publication/220184531_Efficient_Clipping_of_Arbitrary_Polygons/link/0912f510a5ac9191e9000000/download]()
 - [ ] add some visualization (`Makie`?)
 - [ ] export to SVG/STL/PLY
   - [ ] `MeshIO`
# Extras
 - [ ] add a display method that shows the tree
 - [ ] try 2 strategies for n-ary intersection/union:
   - [ ] merge all structures and compute multiplicity with ray-shooting,
	 - [x] or reduce with binary op
 - [ ] icosphere (from Blender)
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
