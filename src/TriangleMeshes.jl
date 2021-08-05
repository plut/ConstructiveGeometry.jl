module TriangleMeshes
using StaticArrays

const libiglboolean_so = joinpath(@__DIR__, "libiglboolean.so")
# the bits types are hard-coded on the C side:
const Point=SVector{3,Cdouble}
const Face=NTuple{3,Cint}
@inline _face(a)=(Cint.(a[1:3])...,)

struct TriangleMesh{A}
	vertices::Vector{Point}
	faces::Vector{Face}
	attributes::Vector{A}
	@inline TriangleMesh{A}(v, f, a) where{A} = new{A}(v, _face.(f), a)
	@inline TriangleMesh(v, f, a::AbstractVector{A}) where{A} =
		TriangleMesh{A}(v, f, a)
end

@inline points(m::TriangleMesh) = m.vertices
@inline faces(m::TriangleMesh) = m.faces
@inline nvertices(m::TriangleMesh) = size(m.vertices, 1)
@inline nfaces(m::TriangleMesh) = size(m.faces, 1)
@inline shift(f, k) = f .+ Face((k,k,k))

function boolean(op, m1::TriangleMesh{A}, m2::TriangleMesh{A}) where{A}#««
	n = nfaces(m1)
	nv3 = Ref(Cint(0))
	nf3 = Ref(Cint(0))
	v3 = Ref(Ptr{Cdouble}(0))
	f3 = Ref(Ptr{Cint}(0))
	j = Ref(Ptr{Cint}(0));
	@inline vpointer(m) = convert(Ptr{Cdouble}, pointer(m.vertices))
	@inline fpointer(f) = convert(Ptr{Cint}, pointer(f))
	r = ccall((:igl_mesh_boolean, libiglboolean_so), Cint,
		(Cint,
		Cint, Cint, Ref{Cdouble}, Ref{Cint},
		Cint, Cint, Ref{Cdouble}, Ref{Cint},
		Ref{Cint}, Ref{Cint}, Ref{Ptr{Cdouble}}, Ref{Ptr{Cint}}, Ref{Ptr{Cint}}),
		op,
		nvertices(m1), nfaces(m1), vpointer(m1), fpointer(shift.(m1.faces, -1)),
		nvertices(m2), nfaces(m2), vpointer(m2), fpointer(shift.(m2.faces, -1)),
		nv3, nf3, v3, f3, j)
	rv3 = unsafe_wrap(Array, convert(Ptr{Point},v3[]), Int(nv3[]); own=true)
	rf3 = unsafe_wrap(Array, convert(Ptr{Face}, f3[]), Int(nf3[]); own=true)
	index = unsafe_wrap(Array, j[], (Int(nf3[]),); own=true) .+ 1
	a3 = [ i ≤ n ? m1.attributes[i] : m2.attributes[i-n] for i in index ]
	return TriangleMesh{A}(rv3, shift.(rf3, 1), a3)
end#»»

@inline Base.union(m1::TriangleMesh, m2::TriangleMesh) = boolean(0, m1, m2)
@inline Base.intersect(m1::TriangleMesh, m2::TriangleMesh) = boolean(1, m1, m2)
@inline Base.setdiff(m1::TriangleMesh, m2::TriangleMesh) = boolean(2, m1, m2)

@inline Base.union(m1::TriangleMesh, m2::TriangleMesh, m::TriangleMesh...) =
	union(union(m1, m2), m...)
@inline Base.intersect(m1::TriangleMesh, m2::TriangleMesh, m::TriangleMesh...)=
	intersect(intersect(m1, m2), m...)

export TriangleMesh
end
