module TriangleMeshes
using StaticArrays
using FastClosures
using DataStructures # heappush!, 
using IGLWrap_jll
# libiglwrap="/home/jerome/src/iglwrap/local/libiglwrap.so"
# the bits types are hard-coded on the C side:
const Point=SVector{3,Cdouble}
const Face=NTuple{3,Cint}
@inline _face(a)=(Cint.(a[1:3])...,)

const Vec3d = MVector{3,Cdouble}
const Mat3d = MMatrix{3,3,Cdouble,9}

include("CornerTables.jl")
using .CornerTables

# Geometry««1
@inline det(a::SVector{2}, b::SVector{2}) = a[1]*b[2]-a[2]*b[1]
@inline det(a::SVector{2}, b::SVector{2}, c::SVector{2}) = det(a-c, b-c)
@inline distance2(a::SVector{3}, b::SVector{3}) = norm2(a-b)
@inline norm2(a::SVector{3}) = a[1]^2+a[2]^2+a[3]^2
# Data type ««1

struct TriangleMesh{T,A}
	vertices::Vector{SVector{3,T}}
	faces::Vector{Face}
	attributes::Vector{A}
	@inline TriangleMesh{T,A}(v, f, a) where{T,A} =
		new{T,A}(v, _face.(f), a)
	@inline TriangleMesh(v::AbstractVector{<:AbstractVector{T}}, f,
		a::AbstractVector{A}) where{T,A} = TriangleMesh{T,A}(v, f, a)
# 	@inline TriangleMesh{T}(v, f, a::AbstractVector{A}) where{T,A} =
# 		TriangleMesh{T,A}(v, f, a)
end

const CTriangleMesh = TriangleMesh{Cdouble}

@inline vertices(m::TriangleMesh) = m.vertices
@inline faces(m::TriangleMesh) = m.faces
@inline attributes(m::TriangleMesh) = m.attributes
@inline nvertices(m::TriangleMesh) = size(m.vertices, 1)
@inline nfaces(m::TriangleMesh) = size(m.faces, 1)
@inline shift(f::Face, k) = f .+ Face((k,k,k))

# IGL interface ««1

# a type to capture what is returned by IGL functions ««
@inline unsafe_array(T, p, n) =
	unsafe_wrap(Array, convert(Ptr{T}, p), n; own=true)
mutable struct Cmesh
	nv::Cint
	nf::Cint
	vertices::Ptr{Cdouble}
	faces::Ptr{Cint}
	@inline Cmesh() = new()
end
@inline fieldpointer(m::T, name) where{T} =
	pointer_from_objref(m) + fieldoffset(T, findfirst(==(name), fieldnames(T)))
@inline nvertices(m::Cmesh) = fieldpointer(m, :nv)
@inline nfaces(m::Cmesh) = fieldpointer(m, :nf)
@inline vpointer(m::Cmesh) = fieldpointer(m, :vertices)
@inline fpointer(m::Cmesh) = fieldpointer(m, :faces)
CTriangleMesh(m::Cmesh, a::AbstractVector{A}) where{A} = CTriangleMesh{A}(
	unsafe_array(Point, m.vertices, m.nv),
	unsafe_array(Face , m.faces   , m.nf),
	a)

@inline vpointer(m::TriangleMesh{Cdouble}) =
	convert(Ptr{Cdouble}, pointer(m.vertices))
@inline fpointer(m::TriangleMesh) = convert(Ptr{Cint}, pointer(m.faces))

#»»

function boolean(op, m1::CTriangleMesh{A}, m2::CTriangleMesh{A}) where{A}#««
	n = nfaces(m1)
	cm = Cmesh()
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.mesh_boolean(op::Cint,
		nvertices(m1)::Cint, nfaces(m1)::Cint,
		vpointer(m1)::Ref{Cdouble}, fpointer(m1)::Ref{Cint},
		nvertices(m2)::Cint, nfaces(m2)::Cint,
		vpointer(m2)::Ref{Cdouble}, fpointer(m2)::Ref{Cint},
		nvertices(cm)::Ptr{Cint}, nfaces(cm)::Ptr{Cint},
		vpointer(cm)::Ptr{Ptr{Cdouble}}, fpointer(cm)::Ptr{Ptr{Cint}},
		j::Ref{Ptr{Cint}})::Cint
	index = unsafe_array(Cint, j[], cm.nf)
	ao = [ i ≤ n ? m1.attributes[i] : m2.attributes[i-n] for i in index ]
	return CTriangleMesh(cm, ao)
end#»»
function minkowski_sum(m1::CTriangleMesh{A}, m2::CTriangleMesh{A}) where{A}#««
	n = nfaces(m1)
	cm = Cmesh()
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.minkowski_sum(
		nvertices(m1)::Cint, nfaces(m1)::Cint,
		vpointer(m1)::Ref{Cdouble}, fpointer(m1)::Ref{Cint},
		nvertices(m2)::Cint, nfaces(m2)::Cint,
		vpointer(m2)::Ref{Cdouble}, fpointer(m2)::Ref{Cint},
		3::Cint,
		nvertices(cm)::Ptr{Cint}, nfaces(cm)::Ptr{Cint},
		vpointer(cm)::Ptr{Ptr{Cdouble}}, fpointer(cm)::Ptr{Ptr{Cint}},
		j::Ref{Ptr{Cint}})::Cint
# 	index = unsafe_array(Cint, j[], cm.nf)
# 	ao = fill(first(attributes(m1)), size(index))
	return CTriangleMesh(cm, fill(first(m1.attributes), cm.nf))
# 	ao = [ i ≤ n ? m1.attributes[i] : m2.attributes[i-n] for i in index ]
# 	return (CTriangleMesh{A}(rvo, rfo, ao), index)
end#»»
function minkowski_sum(m1::CTriangleMesh{A}, v2::Vector{Point},
		e2::Vector{NTuple{2,Cint}}) where{A}#««
	n = nfaces(m1)
	cm = Cmesh()
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.minkowski_sum(
		nvertices(m1)::Cint, nfaces(m1)::Cint,
		vpointer(m1)::Ref{Cdouble}, fpointer(m1)::Ref{Cint},
		length(v2)::Cint, length(e2)::Cint,
		pointer(v2)::Ptr{Cdouble}, pointer(e2)::Ptr{NTuple{2,Cint}},
		2::Cint,
		nvertices(cm)::Ptr{Cint}, nfaces(cm)::Ptr{Cint},
		vpointer(cm)::Ptr{Ptr{Cdouble}}, fpointer(cm)::Ptr{Ptr{Cint}},
		j::Ref{Ptr{Cint}})::Cint
	# index = unsafe_array(Cint, j[], cm.nf)
# 	ao = fill(first(attributes(m1)), size(index))
# 	return (rvo, rfo, index)
	return CTriangleMesh(cm, fill(first(attributes(m1)), cm.nf))
# 	ao = [ i ≤ n ? m1.attributes[i] : m2.attributes[i-n] for i in index ]
# 	return (CTriangleMesh{A}(rvo, rfo, ao), index)
end#»»
function ispwn(m::TriangleMesh)#««
	r = ccall((:mesh_is_pwn, libiglwrap), Cint,
		(Cint, Cint, Ref{Cdouble}, Ref{Cint},),
		nvertices(m), nfaces(m), vpointer(m), fpointer(m),
		)
	return (r ≠ 0)
end#»»
function offset(m::CTriangleMesh{A}, level::Real, grid::Integer) where{A}#««
	cm = Cmesh()
	r = @ccall libiglwrap.offset_surface(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		level::Cdouble, grid::Cint,
		nvertices(cm)::Ptr{Cint}, nfaces(cm)::Ptr{Cint},
		vpointer(cm)::Ptr{Ptr{Cdouble}}, fpointer(cm)::Ptr{Ptr{Cint}})::Cint
	@assert r == 0

	return CTriangleMesh(cm, fill(first(m.attributes), cm.nf))
end#»»
function decimate(m::CTriangleMesh{A}, max_faces::Integer) where{A}#««
	cm = Cmesh()
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.decimate(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		max_faces::Cint,
		nvertices(cm)::Ptr{Cint}, nfaces(cm)::Ptr{Cint},
		vpointer(cm)::Ptr{Ptr{Cdouble}}, fpointer(cm)::Ptr{Ptr{Cint}},
		j::Ref{Ptr{Cint}})::Cint
	@assert r == 0

	index = unsafe_array(Cint, j[], cm.nf)
	return CTriangleMesh(cm, m.attributes[index])
end#»»
function loop(m::CTriangleMesh{A}, count::Integer) where{A}#««
	cm = Cmesh()
	r = @ccall libiglwrap.loop(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		count::Cint,
		nvertices(cm)::Ptr{Cint}, nfaces(cm)::Ptr{Cint},
		vpointer(cm)::Ptr{Ptr{Cdouble}}, fpointer(cm)::Ptr{Ptr{Cint}})::Cint
	@assert r == 0

	return CTriangleMesh(cm, attributes(m)[fld1.(1:cm.nf, 4^count)])
end#»»
function centroid_and_volume(m::CTriangleMesh)#««
	c = Vec3d(undef)
	v = @ccall libiglwrap.centroid_and_volume(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		c::Ref{Vec3d})::Cdouble
	return (c, v)
end#»»
function halfspace(direction, origin, m::CTriangleMesh{A}, color::A) where{A}#««
	cm = Cmesh()
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.intersect_with_half_space(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		Vec3d(origin...)::Ref{Vec3d},
		Vec3d(direction...)::Ref{Vec3d},
		nvertices(cm)::Ptr{Cint}, nfaces(cm)::Ptr{Cint},
		vpointer(cm)::Ptr{Ptr{Cdouble}}, fpointer(cm)::Ptr{Ptr{Cint}},
		j::Ref{Ptr{Cint}})::Cint
	@assert r == 0
	index = unsafe_array(Cint, j[], cm.nf)

	return CTriangleMesh(cm,
		[(i <= nfaces(m) ? m.attributes[i] : color) for i in index])
end#»»
function swept_volume(m::CTriangleMesh, transform, #««
	nsteps, gridres, isolevel = 0)
	cm = Cmesh()
	f = (t,a,b) -> begin
		u, v = transform(t)
		unsafe_store!(a, Mat3d(u)) # column-major, same as Eigen
		unsafe_store!(b, Vec3d(v))
		return
	end
	ctransform = @cfunction($f, Cvoid, (Cdouble, Ptr{Mat3d}, Ptr{Vec3d}))
	r = @ccall libiglwrap.swept_volume(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		ctransform::Ptr{Cvoid}, nsteps::Csize_t, gridres::Csize_t,
		isolevel::Cdouble,
		nvertices(cm)::Ptr{Cint}, nfaces(cm)::Ptr{Cint},
		vpointer(cm)::Ptr{Ptr{Cdouble}}, fpointer(cm)::Ptr{Ptr{Cint}})::Cint
	@assert r == 0

	return CTriangleMesh(cm, fill(first(attributes(m)), cm.nf))
end#»»

@inline Base.union(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(0, m1, m2)
@inline Base.intersect(m1::CTriangleMesh, m2::CTriangleMesh)= boolean(1, m1, m2)
@inline Base.setdiff(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(2, m1, m2)
@inline Base.symdiff(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(3, m1, m2)

# Own functions ««1
# plane_splice ««2
"""
    plane_slice(m::TriangleMesh)

Returns the set of all edges formed by this mesh ∩ the horizontal plane,
as `(vertices, edges)`, where `vertices` are 2d points,
and `edges` are indices into `vertices`.
"""
function plane_slice(z::Real, m::TriangleMesh)
	# build a list of intersection points + connectivity
	# each intersection point is either:
	#  - a vertex v, represented as (v, 0)
	#  - in the edge (v1v2), represented as (v1,v2)
	points = Dict{NTuple{2,Int},Int}()
	elist = NTuple{2,Int}[]
	pindex = @closure v->get!(points, extrema(v), length(points)+1)
	edge! = @closure (v,w)-> push!(elist, minmax(pindex(v),pindex(w)))
	# build list of all edges:
	for (i1, i2, i3) in faces(m)
		(v1, v2, v3) = vertices(m)[[i1,i2,i3]]
		(z1, z2, z3) = (v1[3], v2[3], v3[3]) .- z
		if z1 == 0#««
			if z2 == 0
				z3 == 0 && continue # 000: ignore horizontal triangle
				edge!(i1=>0,i2=>0) # 00+, 00-
			elseif z2 > 0
				z3 == 0 && edge!(i1=>0,i3=>0)
				z3 < 0 && edge!(i1=>0,i2=>i3)
			else # z2 < 0
				z3 == 0 && edge!(i1=>0,i3=>0)
				z3 > 0 && edge!(i1=>0,i2=>i3)
			end
		elseif z1 > 0
			if z2 == 0
				z3 == 0 && edge!(i2=>0,i3=>0) # +00
				z3 < 0 && edge!(i2=>0,i1=>i3) # +0-
			elseif z2 > 0
				z3 < 0 && edge!(i1=>i3,i2=>i3) # ++-
			else # z2 < 0
				z3 == 0 && edge!(i3=>0,i1=>i2)
				z3 > 0 && edge!(i1=>i2,i2=>i3)#+-+
				z3 < 0 && edge!(i1=>i2,i1=>i3)#+--
			end
		else # z1 < 0
			if z2 == 0
				z3 == 0 && edge!(i2=>0,i3=>0) # -00
				z3 > 0 && edge!(i2=>0,i1=>i3) # -0+
			elseif z2 > 0
				z3 == 0 && edge!(i3=>0,i1=>i2)
				z3 > 0 && edge!(i1=>i2,i1=>i3)#-++
				z3 < 0 && edge!(i1=>i2,i2=>i3)#-+-
			else # z2 < 0
				z3 > 0 && edge!(i1=>i3,i2=>i3) # --+
			end
		end#»»
	end
	vlist = Vector{SVector{2,Float64}}(undef, length(points))
	for ((i1,i2),j) in pairs(points)
		if iszero(i1)
			v = vertices(m)[i2]
			vlist[j] = SA[v[1],v[2]]
		else
			(v1,v2) = vertices(m)[[i1,i2]]
			(z1,z2) = (v1[3]-z,v2[3]-z); f = 1/(z2-z1)
			# (z2 v1 - z1 v2)/(z2-z1)
			vlist[j] = SA[f*(z2*v1[1]-z1*v2[1]), f*(z2*v1[2]-z1*v2[2])]
		end
	end
	return (vlist, unique!(sort!(elist)))
end


# project««2
"""
    project(mesh)

Returns the list of triangles produced as the 2d projection of `mesh`.
"""
function project(m::TriangleMesh{T}) where{T}
	triangles = SVector{3,SVector{2,T}}[]
	for (i,j,k) in faces(m)
		a,b,c = vertices(m)[SA[i,j,k]]
		a1,b1,c1 = (SA[a[1],a[2]], SA[b[1],b[2]], SA[c[1],c[2]])
		d = det(a1,b1,c1)
		iszero(d) && continue
		push!(triangles, d > 0 ? SA[a1,b1,c1] : SA[a1,c1,b1])
	end
	return triangles
end

# Corner table structure ««2
struct CTMesh{J,T,A} <: AbstractTriangulation{J}
	triangulation::CornerTable{J}
	points::Vector{SVector{3,T}}
	attributes::Vector{A}
end

Base.show(io::IO, t::CTMesh) = CornerTables.showall(io, t)

@inline CornerTables.triangulation(t::CTMesh) = t.triangulation
@inline point(t::CTMesh, c::Cell) = t.points[int(c)]
@inline edgelength(t::CTMesh, a::Arrow) =
	distance2(point(t, CornerTables.head(t,a)), point(t, CornerTables.tail(t,a)))

# conversion to and from `TriangleMesh`

@inline CTMesh{J}(m::TriangleMesh{T,A}) where{J,T,A} =
	CTMesh{J,T,A}(CornerTable{J}(faces(m)), vertices(m), attributes(m))
@inline CTMesh(m::TriangleMesh) = CTMesh{Int32}(m)

@inline TriangleMesh(t::CTMesh{J,T,A}) where{J,T,A} =
	TriangleMesh(t.points, [ Int.(f) for f in alltriangles(t) ], t.attributes)

# splitedges ««2
# The following structure is used as a heap for the edges of the mesh,
# sorted by size; we will always remove (and split) the longest edge:
struct VecSortedSet{J,V}
	set::SortedSet{Tuple{V,J}}
	size::Vector{V}
end
@inline VecSortedSet{J}(s::AbstractVector{V}) where{J,V} =
	VecSortedSet{J,V}(SortedSet{Tuple{V,J}}((x,i) for (i,x) in pairs(s)), s)
@inline VecSortedSet(s::AbstractVector) = VecSortedSet{Int}(s)

@inline Base.last(s::VecSortedSet) = last(s.set)[2]
@inline function Base.setindex!(s::VecSortedSet, v, k)
	delete!(s.set, (s.size[k], k))
	s.size[k] = v
	push!(s.set, (s.size[k], k))
	return s
end
@inline function Base.push!(s::VecSortedSet, vlist...)
	for v in vlist
		push!(s.size, v)
		push!(s.set, (last(s.size), length(s.size)))
	end
	return s
end
"""
    splitedges!(m::TriangleMesh, maxlen)

Repeatedly split all edges of `m`, starting by the longest ones,
until no edge has length² > `maxlen`.
"""
@inline splitedges(m::TriangleMesh, maxlen) =
	TriangleMesh(splitedges!(CTMesh(m), maxlen))

function splitedges!(t::CTMesh{J,T,A}, maxlen) where{J,T,A}
	elist = VecSortedSet{J}(T[ a > opposite(t, a) ? zero(T) : edgelength(t, a)
		for a in eacharrow(t) ])
	while true
	e = Arrow(last(elist)); o = opposite(t, e)
	elist.size[e] ≤ maxlen && break
	ne, po = next(e), prev(o)
	o_ne, o_po = opposite(t, ne), opposite(t, po)
	
	ct, ch = CornerTables.tail(t, e), CornerTables.head(t, e)
	pm = (point(t,ct) + point(t,ch))/2
	push!(t.points, pm)
	push!(t.attributes, t.attributes[int(node(e))], t.attributes[int(node(o))])
	(q1, q2, cm) = CornerTables.split!(t, e)
	# update elist accordingly:
	# edge (ct, cm) is now (e, opp(e))  (with e < opp(e))
	# edge (cm, ch) is (side(q1,1) < side(q2,1))
	# edge (cm, cl) is (next(e) < side(q1,3))
	# edge (cm, cr) is (prev(o) < side(q2,2))
	# edge (ch, cl) is (opp(next(e)) < side(q1,2))
	# edge (ch, cr) is (opp(prev(o)) < side(q2,3))
	pt, ph, pl, pr = point(t, ct), point(t, ch),
		point(t, left(t,e)), point(t, right(t,e))
	elist[e] = distance2(ph, pm) # elist[o] is already 0
	push!(elist, distance2(pm, ph), 0, 0, 0, 0, 0)
	elist[ne] = distance2(pm, pl)
	elist[po] = distance2(pm, pr)
	elist[o_ne] = distance2(ph, pl)
	elist[o_po] = distance2(ph, pr)
	end
	return t
end

#  »»1
export TriangleMesh
end
