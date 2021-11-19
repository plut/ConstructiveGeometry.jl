# [Mesh I/O](@id io)

Loading and saving of files can be done through the `FileIO` functions
`load` and `save`.

## File loading

Files can be loaded in this way:
```julia
using FileIO
bust = load("beethoven.ply")
```

The following formats should be supported: `.stl` (ASCII and binary)
`.ply` (ASCII and little-endian/big-endian binary).

TODO: `.dxf`

## File saving

Files can be saved in this way:
```julia
using FileIO
save("foo.stl", object)
```

The following formats are supported: `.stl` (ASCII)
and `.ply` (ASCII) for volumes;
`.svg` for shapes.

Image file formats (`.png`, `.pdf`) are also supported;
they are delegated to `Makie`'s `plot` function.

## OpenSCAD output

```@docs
ConstructiveGeometry.scad
```

!!! warning "Deprecated"

    Since the constructions from this package have started to diverge
    from those of OpenSCAD, (and also since the possibilities for export
    and visualization have improved), OpenSCAD output is being deprecated
    and may be removed in a future version.

## OpenSCAD to Julia conversion.

Might be possible for a limited subset of OpenSCAD language. TODO.
