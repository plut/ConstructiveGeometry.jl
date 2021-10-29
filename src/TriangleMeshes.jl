module TriangleMeshes
using StaticArrays
using FastClosures
using IGLWrap_jll
libiglwrap="/home/jerome/src/iglwrap/local/libiglwrap.so"

# Data type ««1
@inline det(a::SVector{2}, b::SVector{2}) = a[1]*b[2]-a[2]*b[1]
@inline det(a::SVector{2}, b::SVector{2}, c::SVector{2}) = det(a-c, b-c)

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

# IGL interface ««1

@inline vpointer(m::TriangleMesh{Cdouble}) =
	convert(Ptr{Cdouble}, pointer(m.vertices))
@inline fpointer(m::TriangleMesh) = convert(Ptr{Cint}, pointer(m.faces))
function boolean(op, m1::CTriangleMesh{A}, m2::CTriangleMesh{A}) where{A}#««
	n = nfaces(m1)
	nvo = Ref(Cint(0))
	nfo = Ref(Cint(0))
	vo = Ref(Ptr{Cdouble}(0))
	fo = Ref(Ptr{Cint}(0))
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.mesh_boolean(op::Cint,
		nvertices(m1)::Cint, nfaces(m1)::Cint,
		vpointer(m1)::Ref{Cdouble}, fpointer(m1)::Ref{Cint},
		nvertices(m2)::Cint, nfaces(m2)::Cint,
		vpointer(m2)::Ref{Cdouble}, fpointer(m2)::Ref{Cint},
		nvo::Ref{Cint}, nfo::Ref{Cint},
		vo::Ref{Ptr{Cdouble}}, fo::Ref{Ptr{Cint}}, j::Ref{Ptr{Cint}})::Cint
	rvo = unsafe_wrap(Array, convert(Ptr{Point},vo[]), Int(nvo[]); own=true)
	rfo = unsafe_wrap(Array, convert(Ptr{Face}, fo[]), Int(nfo[]); own=true)
	index = unsafe_wrap(Array, j[], (Int(nfo[]),); own=true);
	ao = [ i ≤ n ? m1.attributes[i] : m2.attributes[i-n] for i in index ]
	return CTriangleMesh{A}(rvo, rfo, ao)
end#»»
function minkowski_sum(m1::CTriangleMesh{A}, m2::CTriangleMesh{A}) where{A}#««
	n = nfaces(m1)
	nvo = Ref(Cint(0))
	nfo = Ref(Cint(0))
	vo = Ref(Ptr{Cdouble}(0))
	fo = Ref(Ptr{Cint}(0))
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.minkowski_sum(
		nvertices(m1)::Cint, nfaces(m1)::Cint,
		vpointer(m1)::Ref{Cdouble}, fpointer(m1)::Ref{Cint},
		nvertices(m2)::Cint, nfaces(m2)::Cint,
		vpointer(m2)::Ref{Cdouble}, fpointer(m2)::Ref{Cint}, 3::Cint,
		nvo::Ref{Cint}, nfo::Ref{Cint},
		vo::Ref{Ptr{Cdouble}}, fo::Ref{Ptr{Cint}}, j::Ref{Ptr{Cint}})::Cint
	rvo = unsafe_wrap(Array, convert(Ptr{Point},vo[]), Int(nvo[]); own=true)
	rfo = unsafe_wrap(Array, convert(Ptr{Face}, fo[]), Int(nfo[]); own=true)
	index = unsafe_wrap(Array, j[], (Int(nfo[]),2); own=true);
	ao = fill(first(attributes(m1)), size(index))
# 	return (rvo, rfo, index)
	return CTriangleMesh{A}(rvo, rfo, fill(first(m1.attributes), length(rfo)))
# 	ao = [ i ≤ n ? m1.attributes[i] : m2.attributes[i-n] for i in index ]
# 	return (CTriangleMesh{A}(rvo, rfo, ao), index)
end#»»
function minkowski_sum(m1::CTriangleMesh{A}, v2::Vector{Point},
		e2::Vector{NTuple{2,Cint}}) where{A}#««
	n = nfaces(m1)
	nvo = Ref(Cint(0))
	nfo = Ref(Cint(0))
	vo = Ref(Ptr{Cdouble}(0))
	fo = Ref(Ptr{Cint}(0))
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.minkowski_sum(
		nvertices(m1)::Cint, nfaces(m1)::Cint,
		vpointer(m1)::Ref{Cdouble}, fpointer(m1)::Ref{Cint},
		length(v2)::Cint, length(e2)::Cint,
		pointer(v2)::Ptr{Cdouble}, pointer(e2)::Ptr{NTuple{2,Cint}},
		2::Cint,
		nvo::Ref{Cint}, nfo::Ref{Cint},
		vo::Ref{Ptr{Cdouble}}, fo::Ref{Ptr{Cint}}, j::Ref{Ptr{Cint}})::Cint
	rvo = unsafe_wrap(Array, convert(Ptr{Point},vo[]), Int(nvo[]); own=true)
	rfo = unsafe_wrap(Array, convert(Ptr{Face}, fo[]), Int(nfo[]); own=true)
	index = unsafe_wrap(Array, j[], (Int(nfo[]),2); own=true);
	ao = fill(first(attributes(m1)), size(index))
# 	return (rvo, rfo, index)
	return CTriangleMesh{A}(rvo, rfo, fill(first(m1.attributes), length(rfo)))
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
	nvo = Ref(Cint(0))
	nfo = Ref(Cint(0))
	vo = Ref(Ptr{Cdouble}(0))
	fo = Ref(Ptr{Cint}(0))
	r = @ccall libiglwrap.offset_surface(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		level::Cdouble, grid::Cint,
		nvo::Ref{Cint}, nfo::Ref{Cint}, vo::Ref{Ptr{Cdouble}}, fo::Ref{Ptr{Cint}}
		)::Cint
	@assert r == 0

	rvo = unsafe_wrap(Array, convert(Ptr{Point},vo[]), Int(nvo[]); own=true)
	rfo = unsafe_wrap(Array, convert(Ptr{Face}, fo[]), Int(nfo[]); own=true)
	return CTriangleMesh{A}(rvo, rfo,
		fill(first(m.attributes), nfo[]))
end#»»
function decimate(m::CTriangleMesh{A}, max_faces::Integer) where{A}#««
	nvo = Ref(Cint(0))
	nfo = Ref(Cint(0))
	vo = Ref(Ptr{Cdouble}(0))
	fo = Ref(Ptr{Cint}(0))
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.decimate(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		max_faces::Cint,
		nvo::Ref{Cint}, nfo::Ref{Cint}, vo::Ref{Ptr{Cdouble}}, fo::Ref{Ptr{Cint}},
		j::Ref{Ptr{Cint}})::Cint
	@assert r == 0

	rvo = unsafe_wrap(Array, convert(Ptr{Point},vo[]), Int(nvo[]); own=true)
	rfo = unsafe_wrap(Array, convert(Ptr{Face}, fo[]), Int(nfo[]); own=true)
	index = unsafe_wrap(Array, j[], (Int(nfo[]),); own=true);
	return CTriangleMesh{A}(rvo, rfo, [m.attributes[i] for i in index])
end#»»
function loop(m::CTriangleMesh{A}, count::Integer) where{A}#««
	nvo = Ref(Cint(0))
	nfo = Ref(Cint(0))
	vo = Ref(Ptr{Cdouble}(0))
	fo = Ref(Ptr{Cint}(0))
	r = @ccall libiglwrap.loop(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		count::Cint,
		nvo::Ref{Cint}, nfo::Ref{Cint}, vo::Ref{Ptr{Cdouble}}, fo::Ref{Ptr{Cint}}
		)::Cint
	@assert r == 0

	rvo = unsafe_wrap(Array, convert(Ptr{Point},vo[]), Int(nvo[]); own=true)
	rfo = unsafe_wrap(Array, convert(Ptr{Face}, fo[]), Int(nfo[]); own=true)
	return CTriangleMesh{A}(rvo, rfo,
		[ attributes(m)[fld1(i, 4^count)] for i in 1:length(rfo)])
end#»»
mutable struct Vec3d
	x::Float64
	y::Float64
	z::Float64
end
function centroid_and_volume(m::CTriangleMesh)
	nvo = Ref(Cint(0))
	nfo = Ref(Cint(0))
	vo = Ref(Ptr{Cdouble}(0))
	fo = Ref(Ptr{Cint}(0))
	c = Vec3d(0,0,0)
	v = @ccall libiglwrap.centroid_and_volume(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		c::Ref{Vec3d})::Cdouble
	return (SA[c.x,c.y,c.z], v)
end
function halfspace(direction, origin, m::CTriangleMesh{A}, color) where{A}
	nvo = Ref(Cint(0))
	nfo = Ref(Cint(0))
	vo = Ref(Ptr{Cdouble}(0))
	fo = Ref(Ptr{Cint}(0))
	j = Ref(Ptr{Cint}(0));
	r = @ccall libiglwrap.intersect_with_half_space(
		nvertices(m)::Cint, nfaces(m)::Cint,
		vpointer(m)::Ref{Cdouble}, fpointer(m)::Ref{Cint},
		Vec3d(origin...)::Ref{Vec3d},
		Vec3d(direction...)::Ref{Vec3d},
		nvo::Ref{Cint}, nfo::Ref{Cint}, vo::Ref{Ptr{Cdouble}}, fo::Ref{Ptr{Cint}},
		j::Ref{Ptr{Cint}})::Cint
	@assert r == 0

	rvo = unsafe_wrap(Array, convert(Ptr{Point},vo[]), Int(nvo[]); own=true)
	rfo = unsafe_wrap(Array, convert(Ptr{Face}, fo[]), Int(nfo[]); own=true)
	index = unsafe_wrap(Array, j[], (Int(nfo[]),); own=true);
	return CTriangleMesh{A}(rvo, rfo,
		[(i <= nfaces(m) ? m.attributes[i] : color) for i in index])
end

@inline Base.union(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(0, m1, m2)
@inline Base.intersect(m1::CTriangleMesh, m2::CTriangleMesh)= boolean(1, m1, m2)
@inline Base.setdiff(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(2, m1, m2)
@inline Base.symdiff(m1::CTriangleMesh, m2::CTriangleMesh) = boolean(3, m1, m2)

# Own functions ««1


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


# project««
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
end#»»

#  »»1

export TriangleMesh
end
