# Corner tables

This describes the data type used for 3d meshing.

```@meta
CurrentModule = ConstructiveGeometry.CornerTables
```

```@docs
CornerTables
```

## The `CornerTable` data structure.
```@docs
CornerTable
```
## Iterators

 - `star`: iterates around a vertex.

## Elementary mesh modification functions

Vertices, fans and corners are deleted by moving the last entry of the
table in their place (and then shortening the table).
This means that, for example,
if a function manipulates two vertices ``v_1`` and ``v_2`` and deletes ``v_1``,
the value of ``v_2`` might need to be modified.
For this reason, all data-moving functions take, as an optional argument,
an array of values which might need to be modified.
For example, consider the following code:
```
    v = [some_vertex, some_other_vertex]
    delete_vertex!(table, v[1], v)
```
If `v[2]` was the last vertex before the `delete_vertex!` call,
then it is modified to point to the new location of this vertex
(i.e. the previous value of `v[1]`).

These functions update the array via the `Base.replace!` function;
thus, only small arrays (``O(1)`` length) should be passed.

This interface is used by the following functions:
`{move,delete}_{face,vertex,fan}!`.

`delete_vertices!` is passed a list of vertices and consistently deletes
these (renumbering them on-the-fly).
