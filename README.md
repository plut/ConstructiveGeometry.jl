# ConstructiveGeometry


[![Documentation|Dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://plut.github.io/ConstructiveGeometry.jl/dev/)

A module providing a syntax and tool
to define CAD objects directly as Julia scripts.

Currently a work in progress. What should work, at least in part:

 - an easy-to-use syntax for defining solids and CSG operations;
 - meshing of most 2d and 3d operations;
 - output format: OpenSCAD file.

This means that it should already be possible to use this module to
define basic 3d models.


The roadmap now includes at least the following:
 - interfacing with `Makie.jl` for in-built visualization of models;
 - direct output to STL files (and input of STL files);
 - add more convenience constructors (e.g. for transforms);
 - see if 3d algorithms could be made faster with better data structures
	 (e.g. halfedge structures);
 - implement missing operators: 3d Minkowski sum; 2d/3d Minkowski
	 difference;
 - add an annotation system to ease the design of complex models.

Any contributions are welcome!
