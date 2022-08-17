# 1d paths
 * use `AbstractGeometry{1}` for this!
 * concrete subtypes include segments, circular arcs, splines, loop closing
 * allow moves to be relative as well as absolute
 * circular arc: defined either by angle+center+radius, sagitta+endpoints...
 * meshing methods produce polylines
 * concrete type for polyline
 * operators 1d ↔ 2d: border, interior
 * `path_extrude(surface, profile) => path_extrude(∂surface)`
 * `path(element1, element2, ...)`
# Bug fixes
 - [ ] `path_extrude` of shape with holes: hole extrusions are reversed
  - [ ] instead of trying funny stuff with symdiff, we could just concatenate everything and do a self-union to regularize
 - [x] `matrix*set_parameters` etc.; `raise(.8)*mat*cone(3)*lozenge`...
 - [x] `raise` in 2d
 - [x] volume - volume - volume
 - [x] implement *n*-ary `setdiff`
 - [x] hull of 2d points: hull([1,1], ...)
 - [x] linear extrude of repositioned 2d mesh: use normal (this would
   also fix Minkowski)
 - [x] `polygon([[0,0], bezier(...)...])`: element type
 - [x] simplify crossing polygon
 - [x] rotation with axis
 - [x] highlighed sphere difference
 - [x] `minkowski(volume, polygon)` seems broken
 - [x] `rotate_extrude()` with slide and polygon touching y-axis
 - [x] `linear_extrude(scale=0)` is a cone: merge points on top
 - [x] fix problem of triangulating tangent polygons: separately
   triangulate each connected component
 - [x] clarify priority: `linear_extrude(8)*(5*object)`
 - [x] projection of hollow sphere does not work: replace the temporary
   fix by something better (but this is likely Clipper's fault)
 - [x] simple syntax for making conditionals (⇒ use those empty objects)
# Simple fixes
 - [ ] make default `children` return an error and add
   `AbstractGeometryLeaf`
 - [x] `atol`/`rtol` ?
 - [x] auto-compute `offset` npoints from meshing options
 - [x] likewise for `sweep` (volume case)
 - [x] check whether `translate(highlight()*s)` works
   * it should be enough to check that all transformations apply to highlights
 - [x] make a nice logo (threaded bolt? some variant of Julia logo?)
 - [x] bring back some `Angle` types (with degrees) to allow overloading
   of `extrude`, e.g. `extrude(90°)`; likewise, use `cospi` and `sinpi`
   - [x] `using Unitful: °` is probably (almost) enough
 - [x] add a `symmetry` parameter for circles
   - [ ] (and spheres?)
 - [ ] a move system (= a representation of abstract affine rotations)
  - [ ] allow `NamedTuple` for this
  - [ ] possible via `move(origin, s...; direction, spin)`
# Code cleaning
 - [x] remove `DataStructures` dependency
 - [x] `import`: stl, amf, ply, dxf (-> **load**)
   - [ ] stl and ply would follow from integrating `Meshes.jl`
For MeshIO: stl needs
    decompose(Point3f0, mesh)
    decompose(GLTriangleFace, mesh)
    decompose_normals(mesh)
   - [ ] this still needs to somehow merge points in stl (e.g. self-union) ?
 - [x] replace `include` by something using `FileIO` e.g. `.cg.jl`?
 - [x] https://www.usenix.org/legacy/event/usenix05/tech/freenix/full_papers/kirsch/kirsch.pdf : CSG tree normalization (only useful when using `convexity` rendering...)
 - [x] split `Offset` in two structures `OffsetShape` and `OffsetVolume`
 - [x] make transformations even lazier, so that they are evaluated only once
   their subjects (and more importantly, their dimension) are known
 - [x] check if operator associativity is still needed
 - [ ] overload `extrude()` (for paths, angles, numbers)
 - [ ] move doc examples to `WGLMakie`
 - [ ] remove `ConvexHull.jl` dependency (use method from `igl` instead?)
 - [ ] `linear_extrude` / `prism` ?
 - [ ] `rotate_extrude` / `revolution` ?
 - [ ] replace ad-hoc `plot` methods by correct `Makie` interface
 - [x] make meshing type-stable
 - [x] update `iglwrap` and use latest version
   - [x] `loop`
   - [x] `minkowski_sum` (might need tetrahedralize)
   - [ ] `upsample` (only useful once we have deformations)
   - [x] `swept_volume`
   - [x] `centroid`
   - [ ] `convex_hull`
   - [x] `offset_surface`
   - [x] `intersect_with_half_space`
 - [x] `cylinder(..., center)`
 - [x] `cylinder(h, r1, r2)`
 - [x] add a parameter to circles, spheres and cylinders to mesh them as
   circumscribed
 - [ ] document how to extend (e.g. new object type)
# New features
## Geometry
 - [ ] `multiply`: dispose copies of a mesh (better than `union` because
   we compute the child mesh only once)
 - [ ] replace `attributes` by a pointer to the original colored object; this could allow detecting edges etc.
 - [x] `refine`: shorten all edges until no longer than given length
  - [ ] 2d (just divide edges)
  - [x] 3d (split triangles)
 - [ ] path type: e.g. take a rectangular path, round corners, then create the shape defined by 1-directional offset
 - [ ] ellipse
 - [ ] geometry operations:
  - [ ] point in shape,
  - [ ] random sample of shape,
  - [ ] area/volume of shape
 - [ ] add some relational definitions: mutual tangent, intersection, etc.
 - [ ] path sweep: use straight skeleton/medial axis for disappearing vertices
  - [x] compute medial axis
  - [x] compute straight skeleton
 - [ ] volume \ surface := volume \ extrude(surface)
 - [ ] propagate `atol` and `rtol` through transformation matrices (use largest eigenvalue)
 - [x] equivalent of OpenSCAD's for loop?!
 - [x] `linear_extrude` with twist and scale
 - [x] `rotate_extrude` with slide (per-turn) along the axis
 - [ ] compute center of gravity (and use it for scaling etc.)
 - [ ] add a 1d type (points, segments; paths) for minkowski (/ extrusions)?
   - [ ] this makes sense; `Clipper.jl` seems happy to do Minkowski with a path
 - [ ] define a path type (for Minkowski + stroke) ?
  - [ ] probably useful for 3d paths at least
 - [ ] `TriangleMeshes`: have a way to detect non-pwm meshes *and explain
   why*
   - [ ] fix meshes on stl import: *a bit harder, IGL does not have a function for this*
 - [x] also in 2d, regularize polygons (e.g. one backwards loop)
 - [x] triangulate faces of `Surface`
 - [x] allow self-union (for fixing meshes)
 - [x] Minkowski sum in mixed dimensions
 - [x] half-plane intersection
 - [x] swept surfaces (`path_extrude`)
  - [ ] allow planar sweep too?
  - [x] volume sweep (use `swept_volume`) ?
  - [ ] fix Clipper's missing sweep? (e.g. adding a few extra points far
    away (preserving tangents) and removing anything close to those points)
    [ ] or write a patch for the C++ library?
 - [x] find a way to fix path extrusion? either
   - [x] write something from scratch...
   - [ ] cut “by hand” the result of a “butt” extrusion;
   - [ ] intersect the result of a custom “fill” extrusion;
   - [ ] patch the `clipper` library...
 - [ ] lofting/skinning
 - [x] allow complement (for intersection)
 - [x] half-space intersection
 - [x] plane intersection: `slice`
 - [x] projection
 - [ ] path wrapping
 - [x] wrapped volumes
 - [ ] `text`
   - OpenSCAD: `src/FreetypeRenderer.cc`, `::render` function
   - [ ] use `Pango` for text and `FreeType` for fonts
 - [x] Bézier curves (used as path for `path_extrude`, `stroke`, `polygon`)
 - [ ] `symmetrize(m, s...) = union(s..., mirror(m, s...))`
   find some easy syntax allowing also `symmetrize(m)*s`
 - [ ] 2d Minkowski difference
## Misc.
 - [x] `surface(volume)` = instantiate as a mesh (avoids recomputation)
 - [ ] `plot(...; zoom=...)`
## Attachments
 - [ ] anchor/attachment system
    attach(square(1),
    :left => circle(3),
    )
    anchor(square(…), [-1,0])
    anchor(square(…), :left)
    square(…, anchor=:left)
 - [ ] this needs allowing “children” for the primitives, e.g.
    square(5)*[
      position(:top) circle(5),
    ]
 - [ ] check that it is easy for the user to define arbitrary `Transform`s.
 - [ ] rewrite `attach` using `Transform`
  - [ ] and allow:
    attach(X) * [
      :left => Y, :right => Z,
    ] # as array *or* tuple
## Speed
 - [x] check `@code_warntype` everywhere for a start
 - [ ] still slow; see if we can resurrect `CornerTables` (using `Clipper` for the hard 2d part).
 - [ ] make objects mutable to store computed mesh
## Syntax
```
    [1,0,0] + annotate("blah", (.5,.5,5))* sphere(3);
    Annotation("blah",(1.5,.5,.5), mesh(sphere(3)))
```
 - [ ] hook them in existing highlight procedure
 - [ ] add a LineNode reference to constructors
   (i.e. first thing in call stack outside module).
 - [ ] `color`: add more properties (e.g. shininess) to be able to show
   e.g. metal parts
 - [ ] find some way of referencing parts of objects,
 e.g. `cylinder().edge(:top)` references the top edge,
 and have it accessible through CSG hierarchy
 (thus, maybe later providing a way to e.g. fillet it?)
 - [ ] get a better syntax for transforms, e.g.
 `symmetrize = transform(axis,s->s ∪ mirror(axis,s))` ??
 - [x] add examples (with images) in documentation
  - [ ] more sophisticated/real examples
 - [x] overload `color*object` for `color::Colorant`
 - [ ] find a way to access `.x`, `.y` and `.z` for `Point` and `Vec`
   types
    [ ] probably requires making `Vec` a separate type from `SArray`
    (to avoid piracy)
 - [ ] think of using `LabelledArrays.jl` (`SLArray`) as an alternative to
   `StaticVector` for `Vec` types
 - [ ] do something for keyword argument dispatch; e.g. the following should be equivalent:
    circle(r=1)
    circle(d=2)
    circle(radius=1)
    circle(diameter=2)
    circle(1)
  and it should be “open” for adding new keywords, e.g.
    circle(1, center=[0,0]) # => translate([0,0])* circle(1)
    square(3, closing=1) # => calls closing(1)*square(3)
# Syntax
 - [ ] using `:` for transforms is very tempting:
```
   color(red):
   translate([-1,1,1-]):
   square(3)
```
 - [x]  `∪`, `∩`, `\`: booleans
 - [x] `+` translation
 - [ ]  `⊕` Minkowski sum
 - [ ] `object ± real` = offset
 - [ ] `- ⊖` Minkowski difference
 - [ ] `:` convex hull ?
 - [x] `×` extrusion (scalar => linear_extrude; interval =>
   rotate_extrude; path => path_extrude)
 - [ ] make kwargs open for user extension (e.g. rounded squares)
 - [ ] replace minkowski with circle by an offset
# Issues in other packages
 - [x] `Rotations.jl`: using the same type for angles and coordinates is not
   terribly useful (in particular with angles in radians).
# Packaging
 - [ ] write a full doc about how to define a new transform
 - [ ] check the list of exports
 - [x] write a minimal regression test
   - [x] the doc is the test
   - [ ] add some more complicated examples in the doc to expand the tests
 - [ ] distinguish between core and sub-packages (implementing BOSL2 stuff)?
# Future
 - [?] [https://www.researchgate.net/publication/220184531_Efficient_Clipping_of_Arbitrary_Polygons/link/0912f510a5ac9191e9000000/download]()
 - [x] export to SVG/STL/PLY
 - [ ] splines (and enclosed area)
 - [ ] NURBS
 - [ ] deformation https://libigl.github.io/tutorial/#chapter-4-shape-deformation
 - [ ] https://student.cs.uwaterloo.ca/~cs779/Gallery/Winter2018/anietoro/doc/
 - [x] add a display method that shows the tree
 - [x] icosphere (from Blender)
 - [ ] improve `unit_n_gon` to take advantage of symmetries
 - [ ] Annotations in 2d
 - [ ] Annotations in 3d (this might depend on the visualizer though)
 - [ ] rewrite Annotations in terms of `Transform`
 - [ ] (more generally, metadata)
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
