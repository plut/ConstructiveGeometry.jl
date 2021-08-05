# ConstructiveGeometry


[![Documentation|Dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://plut.github.io/ConstructiveGeometry.jl/dev/)

Defining CSG objects from within Julia.

This module is currently a work in progress,
roughly at the “proof of concept” stage.
The following work at least in some cases:

 - an easy-to-use syntax for defining solids and CSG operations;
 - meshing of most 2d and 3d operations;
 - output as a STL or OpenSCAD file and viewing.

It is possible to use this module to define basic 3d models:
```julia
s1 = sphere(20)
s2 = cube(20)
mesh(s1 ∪ s2)
mesh(s1 ∩ s2)
mesh(s1 \ s2)
mesh(hull(s1, s2))
```
gives the following output:
![CSG operations on a sphere and a cube](examples/sphere_cube.png)

*(although the meshes were computed by the Julia module,
OpenSCAD was used for rendering the image.
This will change in a future release.)*

It is also possible to export a model as a `.svg` (for 2d shapes)
or `.stl` (for 3d volumes) file.

It would be dishonest to pretend that this project has nothing to do with
OpenSCAD. The goal is to succeed where OpenSCAD did (i.e. in having a
simple syntax for script-base CAD), while fixing its less-fun parts
(which include at least the esoteric programming language and the
“opaque” data types for geometric objects).
Julia is the right tool for this: it is a full-blown programming
language with all the capacity for defining a convenient syntax,
implementing efficient algorithms, and linking with pre-existing
libraries.

Currently (as of 2021-08), `ConstructiveGeometry.jl` uses the following
libraries:
 - `libigl` for boolean mesh operations;
 - `Clipper` for boolean polygon operations;
 - `Triangle` for triangulation;
 - `Polyhedra` and `GLPK` for convex hull;
 - `Makie` for visualization.

The roadmap now includes at least the following:
 - built-in visualization of models;
 - `.stl` file import, and I/O to more formats;
 - add more convenience constructors (e.g. for transforms);
 - implement missing operators: 3d Minkowski sum; 2d/3d Minkowski;
   path extrusion; `text()`;
 - add an annotation system to ease the design of complex models;
 - add an anchor system.

Any contributions are welcome!
