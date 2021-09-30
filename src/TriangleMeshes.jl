module TriangleMeshes
using StaticArrays
using IGLWrap_jll

# the bits types are hard-coded on the C side:
const Point=SVector{3,Cdouble}
const Face=NTuple{3,Cint}
@inline _face(a)=(Cint.(a[1:3])...,)

struct TriangleMesh{T,A}
	vertices::Vector{SVector{3,T}}
	faces::Vector{Face}
	attributes::Vector{A}
	@inline TriangleMesh{T,A}(v, f, a) where{T,A} =
		new{T,A}(v, _face.(f), a)
	@inline TriangleMesh{T}(v, f, a::AbstractVector{A}) where{T,A} =
		TriangleMesh{T,A}(v, f, a)
end

const CTriangleMesh = TriangleMesh{Cdouble}

@inline vertices(m::TriangleMesh) = m.vertices
@inline faces(m::TriangleMesh) = m.faces
@inline attributes(m::TriangleMesh) = m.attributes
@inline nvertices(m::TriangleMesh) = size(m.vertices, 1)
@inline nfaces(m::TriangleMesh) = size(m.faces, 1)
@inline shift(f::Face, k) = f .+ Face((k,k,k))
@inline vpointer(m::TriangleMesh{Cdouble}) =
	convert(Ptr{Cdouble}, pointer(m.vertices))
@inline fpointer(m::TriangleMesh) = convert(Ptr{Cint}, pointer(m.faces))

function boolean(op, m1::CTriangleMesh{A}, m2::CTriangleMesh{A}) where{A}#««
	n = nfaces(m1)
	nv3 = Ref(Cint(0))
	nf3 = Ref(Cint(0))
	v3 = Ref(Ptr{Cdouble}(0))
	f3 = Ref(Ptr{Cint}(0))
	j = Ref(Ptr{Cint}(0));
	r = ccall((:igl_mesh_boolean, libiglwrap), Cint,
		(Cint,
		Cint, Cint, Ref{Cdouble}, Ref{Cint},
		Cint, Cint, Ref{Cdouble}, Ref{Cint},
		Ref{Cint}, Ref{Cint}, Ref{Ptr{Cdouble}}, Ref{Ptr{Cint}}, Ref{Ptr{Cint}}),
		op,
		nvertices(m1), nfaces(m1), vpointer(m1), fpointer(m1),
		nvertices(m2), nfaces(m2), vpointer(m2), fpointer(m2),
		nv3, nf3, v3, f3, j)
	rv3 = unsafe_wrap(Array, convert(Ptr{Point},v3[]), Int(nv3[]); own=true)
	rf3 = unsafe_wrap(Array, convert(Ptr{Face}, f3[]), Int(nf3[]); own=true)
	index = unsafe_wrap(Array, j[], (Int(nf3[]),); own=true);
	a3 = [ i ≤ n ? m1.attributes[i] : m2.attributes[i-n] for i in index ]
	return TriangleMesh{Cdouble,A}(rv3, rf3, a3)
end#»»

function ispwn(m::TriangleMesh)
	r = ccall((:igl_mesh_is_pwn, libiglwrap), Cint,
		(Cint, Cint, Ref{Cdouble}, Ref{Cint},),
		nvertices(m), nfaces(m), vpointer(m), fpointer(m),
		)
	return (r ≠ 0)
end

@inline Base.union(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(0, m1, m2)
@inline Base.intersect(m1::CTriangleMesh, m2::CTriangleMesh)= boolean(1, m1, m2)
@inline Base.setdiff(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(2, m1, m2)
@inline Base.xor(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(3, m1, m2)

# @inline Base.union(m1::TriangleMesh, m2::TriangleMesh, m::TriangleMesh...) =
# 	union(union(m1, m2), m...)
# @inline Base.intersect(m1::TriangleMesh, m2::TriangleMesh, m::TriangleMesh...)=
# 	intersect(intersect(m1, m2), m...)

export TriangleMesh
end
