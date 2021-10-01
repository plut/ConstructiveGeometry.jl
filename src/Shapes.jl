"""
    Shapes

A module encapsulating types and methods for bidimensional shapes.

Exported types: `PolygonXor`

Useful functions:
"""
module Shapes

using LinearAlgebra
using StaticArrays
using FixedPointNumbers
using FastClosures
using Triangle
using DataStructures
import Clipper


norm²(x) = x[1]*x[1]+x[2]*x[2]
distance²(x,y) = norm²(x-y)
Path{D,T} = Vector{SVector{D,T}}
# @inline det(u::StaticVector{2}, v::StaticVector{2}) =
# 	@inbounds u[1]*v[2]-u[2]*v[1]
# @inline det(p::StaticVector{2}, q::StaticVector{2}, r::StaticVector{2}) =
# 	det(q-p, r-p)

# Clipper types and constants««1
# default number of bits for Clipper types
# this *must* be a 64-bit type:
const _CLIPPER_FIXED = Fixed{Int64,16}
const _CLIPPER_ONE = FixedPointNumbers.rawone(_CLIPPER_FIXED)
# constants
const _CLIPPER_ENUM = (#««
	clip=(
		union       =Clipper.ClipTypeUnion,
		intersection=Clipper.ClipTypeIntersection,
		difference  =Clipper.ClipTypeDifference,
		xor         =Clipper.ClipTypeXor,
	),
	# the enum order in `Clipper.jl` is wrong, we fix it here:
	ends=(
		fill  = Clipper.EndType(0), #Clipper.EndTypeClosedPolygon,
		loop  = Clipper.EndType(1), #Clipper.EndTypeClosedLine,
		butt  = Clipper.EndType(2), #Clipper.EndTypeOpenButt,
		square= Clipper.EndType(3), #Clipper.EndTypeOpenSquare,
		round = Clipper.EndType(4), #Clipper.EndTypeOpenRound,
	),
	join=(
		square=Clipper.JoinTypeSquare,
		round =Clipper.JoinTypeRound,
		miter =Clipper.JoinTypeMiter,
	),
	fill=(
		nonzero=Clipper.PolyFillTypeNonZero,
		evenodd=Clipper.PolyFillTypeEvenOdd,
		positive=Clipper.PolyFillTypePositive,
	),
)#»»
# conversions««2

@inline toc(x::Real) = reinterpret(convert(_CLIPPER_FIXED, x))
@inline toc(x::StaticVector{2,<:Real}) = Clipper.IntPoint(toc.(x)...)
@inline toc(x) = toc.(x) # recursively...

@inline fromc(x::Int64) = Float64(reinterpret(_CLIPPER_FIXED,x))
@inline fromc(x::Clipper.IntPoint) = SA[fromc(x.X), fromc(x.Y)]
@inline fromc(x::StaticVector) = SA[fromc.(x)...]
@inline fromc(x) = fromc.(x)

# this is a special case; the `delta` parameter wants a Float64,
# and we *might* use a different type for `delta` than for coordinates:
# @inline to_clipper_float(T::Type{<:Real}, x, n=1)::Float64 = x*_CLIPPER_ONE^n

# numbers...
# Wrappers for Clipper calls««1
# We wrap all Clipper objects in a NamedTuple with the original type
# struct Marked{T,X}
# 	data::X
# 	@inline Marked{T}(x::X) where{T,X} = new{T,X}(x)
# end

# @inline ClipperClip(T::Type) = (T, Clipper.Clip())
# @inline ClipperClip(T::Type) = Marked{T}(Clipper.Clip())
# @inline ClipperOffset(T::Type, miterLimit::Real, roundPrecision::Real) =
# 	Marked{T}(Clipper.ClipperOffset(Float64(miterLimit),
# 		to_clipper_float(T, roundPrecision)))

# @inline add_path!(Type{<:Real}, c, path, args...) =
# 	Clipper.add_path!(c, to_clipper(T, path), args...)
# @inline add_paths!(T::Type{<:Real}, c, paths, args...) =
# 	Clipper.add_paths!(c, to_clipper(T, paths), args...)

# @inline execute(T::Type{<:Real}, c::Clipper.Clip, args...) =
# 	from_clipper(T, Clipper.execute(c, args...)[2])
# @inline execute(T::Type{<:Real}, c::Clipper.ClipperOffset, delta::Real) =
# 	from_clipper(T, Clipper.execute(c, to_clipper_float(T, delta)))

# Clipper calls on Path values««1
@inline function clip(op::Symbol,
		v1::AbstractVector{Path{2,T}},
		v2::AbstractVector{Path{2,T}};
		fill = :evenodd)::Vector{Path{2,T}} where {T}
	c = Clipper.Clip()
	Clipper.add_paths!(c, toc(v1), Clipper.PolyTypeSubject, true) # closed=true
	Clipper.add_paths!(c, toc(v2), Clipper.PolyTypeClip, true)

	f = _CLIPPER_ENUM.fill[fill]
	return fromc(Clipper.execute(c, _CLIPPER_ENUM.clip[op], f, f)[2])
end
@inline function offset(v::AbstractVector{Path{2,T}}, r::Real;
		join = :round,
		ends = :fill,
		miter_limit = 2.,
		precision = 0.2
		)::Vector{Path{2,T}} where{T}
	c = Clipper.ClipperOffset(miter_limit, precision*_CLIPPER_ONE)
	add_paths!(T, c, v, _CLIPPER_ENUM.join[join], _CLIPPER_ENUM.ends[ends])
	return execute(c, r*_CLIPPER_ONE)
end
@inline function offset(v::AbstractVector{Path{2,T}}, r::AbstractVector{<:Real};
		join = :round,
		ends = :fill,
		miter_limit = 2.,
		precision = 0.2
		)::Vector{Vector{Path{2,T}}} where{T}
	# “Simultaneously” computes offset for several offset values.
	# Used by path_extrude().
	c = Clipper.ClipperOffset(miter_limit, precision*_CLIPPER_ONE)
	add_paths!(T, c, v, _CLIPPER_ENUM.join[join], _CLIPPER_ENUM.ends[ends])
	return [ execute(c, ρ*_CLIPPER_ONE) for ρ in r]
end
@inline simplify_paths(p::AbstractVector{Path{2,T}}; fill=:nonzero) where{T} =
	from_clipper(T,
		Clipper.simplify_polygons(to_clipper(T, p), _CLIPPER_ENUM.fill[fill]))
"""
    orientation(p::Path{2})

Returns `true` iff p is a direct loop (i.e. if area >= 0).
"""
@inline orientation(p::Path{2}) = Clipper.orientation(toc(p))
@inline area(p::Path{2}) = Clipper.area(toc(p)) / _CLIPPER_ONE^2

"""
    point_in_polygon(pt, poly::Path{2})

Returns 1 if point is in the interior, -1 on boundary, and 0 outside the
given polygon.

Polygon is assumed not self-intersecting.
"""
@inline point_in_polygon(point::StaticVector{2}, path) =
	Clipper.pointinpolygon(toc(point), toc(path))

# # convex hull of 2d points (monotone chain)««1
# # this is the version using MiniQhull: #««
# # convex_hull(points::AbstractVector{<:AnyVec(2)}) =
# #		# Delaunay triangulation:
# #		let T = MiniQhull.delaunay([p[i] for i in 1:2, p in points]),
# #				N = length(points),
# #				M = zeros(Bool, N, N) # ⚠ memory O(N²)
# #		# mark all edges:
# #		for (a,b,c) in eachcol(T)
# #			b1 = points[b] - points[a]
# #			c1 = points[c] - points[a]
# #			d = b1[1]*c1[2] - b1[2]*c1[1] # determinant === orientation of triangle
# #			if d < 0
# #				M[a,b] = M[b,c] = M[c,a] = true
# #			else
# #				M[b,a] = M[c,b] = M[a,c] = true
# #			end
# #		end
# #		# list of remaining edges (retrograde oriented)
# #		L= sort([(i,j) for i in 1:N, j in 1:N if M[i,j] && ! M[j,i]], by=v->v[1])
# #		next(i) = L[searchsorted(L, i, by=y->y[1])][1][2]
# #		R = zeros(Int, length(L))
# #		R[1:2] .= L[1] # initialize with first edge
# #		for i in 3:length(R)
# #			R[i] = next(R[i-1])
# #		end
# #		# returns in retrograde ordering (OpenSCAD convention):
# #		points[R]
# # end#»»
# """
#     convex_hull_list(points)
# 
# Returns the convex hull of the points, as a list of indexes (in direct
# order, starting at a reproducible index in the list of points).
# """
# function convex_hull_list(points::AbstractVector{<:StaticVector{2}})
# 	list = Int[]
# 	length(points) ≤ 2 && return points
# 	p = sortperm(points) # lex by x then y
# 	k = 1
# 	pdet = @closure (i,j,k) -> begin
# 		qi = points[i]; qj = points[j]; qk = points[k]
# 		vij = qj-qi; vik = qk-qi
# 		return vij[2]*vik[1]-vij[1]*vik[2]
# 	end
# 	for i in p # lower hull
# 		while length(list) ≥ 2 && pdet(list[end-1], list[end], i) ≥ 0
# 			pop!(list)
# 		end
# 		push!(list, i)
# 	end
# 	t = length(list)
# 	for i in reverse(p) # upper hull
# 		while length(list) ≥ t+1 && pdet(list[end-1], list[end], i) ≥ 0
# 			pop!(list)
# 		end
# 		push!(list, i)
# 	end
# 	pop!(list) # remove loop
# 	return list
# end
# """
#     convex_hull([vector of 2d points])
# 
# Returns the convex hull (as a vector of 2d points, ordered in direct
# order).
# """
# @inline convex_hull(points::AbstractVector{<:StaticVector{2}}) =
# 	points[convex_hull_list(points)]
# Convolution««1
function convolution(p::Path, q::Path)
	(np, nq) = (length(p), length(q))
	ep = [p[cyclindex(i, p)]-p[i] for i in eachindex(p)] # edges of p
	eq = [q[cyclindex(i, q)]-q[i] for i in eachindex(q)]
	j0 = 0
	newpoly = similar(p, 0)
	for ip in eachindex(p)
		for iq in eachindex(q)
			iq0 = cyclindex(iq, q, -1)
			if circularcmp(eq[iq0], ep[ip], eq[iq], Val(:offset))
				push!(newpoly, p[ip]+q[iq])
				push!(newpoly, p[cyclindex(ip, p)]+q[iq])
			end
		end
	end
	newpoly
end
p⋆q = convolution(p, q)
# Minkowski sum#««1
# function minkowski(p, q; fill=:nonzero)
# 	r = convolution(p, q)
# 	return simplify_paths([r]; fill)
# end
# # Convolution of polygons««2
# # http://acg.cs.tau.ac.il/tau-members-area/general%20publications/m.sc.-theses/thesis-lienchapter.pdf
# """
#     circularcmp(v1, v2, v3, [Val(:offset)])
# 
# Circular comparison predicate; returns true iff directions of vectors
# `v1`, `v2`, `v3` are arranged in a trigonometric ordering along the unit
# circle.
# 
# If `Val(:offset)` is passed then `v1`, `v3` are infinitesimally rotated
# in the positive direction compared to `v2`.
# """
# function circularcmp(v1, v2, v3)
# 	d1 = v2[1]*v3[2] ≥ v2[2]*v3[1]
# 	d2 = v3[1]*v1[2] ≥ v3[2]*v1[1]
# 	d3 = v1[1]*v2[2] ≥ v1[2]*v2[1]
# 	return (d1+d2+d3) ≥ 2
# end
# function circularcmp(v1, v2, v3, ::Val{:offset})
# 	d1 = v2[1]*v3[2] > v2[2]*v3[1]
# 	d2 = v3[1]*v1[2] ≥ v3[2]*v1[1]
# 	d3 = v1[1]*v2[2] ≥ v1[2]*v2[1]
# 	return (d1+d2+d3) ≥ 2
# end
# 
# Draw path««1
# # """
# #     draw(path, width; kwargs...)
# # 
# #     ends=:round|:square|:butt|:closed
# #     join=:round|:miter|:square
# # """
# # function draw(path::Path{2,T}, width::Real;
# # 		ends::Symbol = :round, join::Symbol = :round,
# # 		miter_limit::Float64 = 2.0, precision::Real = 0.2) where{T}
# # 	CT = clipper_type(T)
# # 	RT = clipper_rettype(T)
# # 	c = ClipperOffset(miter_limit, clipper_float(CT, precision))
# # 	println("join=$join, round=$round")
# # 	Clipper.add_path!(c, clipper_path(path),
# # 		JoinType(Val(join)), EndType(Val(ends)))
# # 	println("$(clipper_type(T)) $(CT(1.)); prec=$(Float64(CT(precision)))")
# # 	ret = clipper_unpath.(RT, Clipper.execute(c, clipper_float(CT, width)/2))
# # 	return PolyUnion(ret)
# # end
# # 
# PolygonXor««1
"""
    PolygonXor{T}(polygon...)
    PolygonXor{T}([points], [points]...)

Exclusive union of several polygons.
(This may be used to represent either the disjoint union of several polygons,
or a polygon with holes, etc.)
"""
struct PolygonXor{T}
	paths::Vector{Vector{SVector{2,T}}}
# 	@inline PolygonXor{T}(paths::AbstractVector{<:AbstractVector{<:StaticVector{2,<:Real}}}) where{T<:Real} =
# 		new{T}(Vector{SVector{2,T}}.(paths))
	@inline PolygonXor{T}(paths::AbstractVector) where{T} = new{T}(paths)
	@inline PolygonXor(paths::AbstractVector) =
		PolygonXor{eltype(eltype(paths))}(paths)
end
@inline PolygonXor{T}(paths::AbstractVector{<:StaticVector{2,<:Real}}...
	) where{T} = PolygonXor{T}([paths...])

@inline coordtype(::Type{PolygonXor{T}}) where{T} = T
@inline coordtype(s::PolygonXor) = coordtype(typeof(s))
@inline paths(p::PolygonXor) = p.paths
@inline vertices(p::PolygonXor) = reduce(vcat, paths(p))

"""
    perimeters(::PolygonXor)

Returns a list of perimeters and holes for this region.
Perimeters are oriented ↺ and holes ↻.
"""
function perimeters(p::PolygonXor)
	firstindex = zeros(Int, length(p.paths))
	firstindex[1] = 0
	for i in 1:length(p.paths)-1
		firstindex[i+1] = firstindex[i] + length(p.paths[i])
	end
	return [[firstindex[i]+1 : firstindex[i] + length(p.paths[i]);]
		for i in eachindex(p.paths)]
end

# Identify holes of PolygonXor
"""
    identify_polygons(s::PolygonXor)

Given a `PolygonXor` defined as the exclusive union of polygons and holes,
returns a vector mapping individual polygons of `s` to:
 - positive indices 1, 2, 3... for the exterior polygons of `s`;
 - negative indices for the holes, indicating in which exterior polygon is the corresponding hole.
"""
function identify_polygons(s::PolygonXor)
	n = length(paths(s))
	m = [ point_in_polygon(first(p), q) for p in paths(s), q in paths(s) ]
	index = zeros(Int, n)
	hole = zeros(Int, n)
	c = 0 # number of connected components
	for i in 1:n
		# if this is a hole, then it will have m[i,j] > 0 for some i
		k = findfirst(m[i,:] .> 0)
		if k == nothing
			index[i] = (c+= 1)
		else
			hole[i] = k
		end
	end
	for i in 1:n
		if hole[i] > 0
			index[i] = -index[hole[i]]
		end
	end
	return index
end

@inline clip(op, s::PolygonXor{T}...) where{T} =
	reduce((p,q)->PolygonXor{T}(clip(op, paths(p), paths(q), fill=:evenodd)), s)

# @inline Base.union(s::PolygonXor...) = clip(:union, s...)
# @inline Base.intersect(s::PolygonXor...) = clip(:intersection, s...)
# @inline Base.setdiff(s1::PolygonXor, s2::PolygonXor) = clip(:difference, s1, s2)

@inline area(shape::PolygonXor) = sum(area.(paths(shape)))
@inline convex_hull(shape::PolygonXor{T}...) where{T} =
	PolygonXor{T}(convex_hull([vertices.(shape)...;]))

# Minkowski sum of polygons and their unions ««2
# function minkowski(vp::AbstractVector{<:AbstractVector{<:StaticVector{2}}},
# 		vq::AbstractVector{<:AbstractVector{<:StaticVector{2}}}; fill=:nonzero)
# 	vr = vec([(convolution(p, q)) for p in vp, q in vq])
# 	return simplify_paths(vr; fill)
# end
# function minkowski(p, q)
# 	return PolygonXor(minkowski(p, q; fill=:evenodd)...)
# end
# function minkowski(p::PolygonXor, q::PolygonXor)
# 	cp = [ coordinates.(vertices(x)) for x in paths(p) ]
# 	cq = [ coordinates.(vertices(y)) for y in paths(q) ]
# 	return PolygonXor(minkowski(cp, cq, fill=:evenodd)...)
# end

# Triangulation and reconstruction from triangle««1
function triangulate(m::PolygonXor)#««
	v = vertices(m)
	id = identify_polygons(m)
	peri = perimeters(m)
	is_hole = falses(length(v))
	edges = Matrix{Int}(undef, sum(length.(peri)), 2)
	c = 0
	for (k, p) in pairs(peri)
		id[k] < 0 && (is_hole[p] .= true)
		n = length(p)
		for i in 1:n, j in 1:2
			edges[c+i,j] = p[mod1(i+j-1,n)]
		end
		c+= n
	end
	tri = constrained_triangulation(
		Matrix{Float64}([transpose.(v)...;]),
		collect(1:length(v)), edges)
	# remove triangle made entirely of hole vertices
	return [ (t[1], t[2], t[3]) for t in tri if !all(is_hole[t]) ]
# 	return tri[[!all(is_hole[t]) for t in tri]]
end#»»
# reconstruction from triangles
# this is used by 3d->2d projection:
function PolygonXor(points::AbstractVector{SVector{2,T}},
	faces::AbstractVector{<:Tuple{<:Integer,<:Integer,<:Integer}}) where{T}
	v=[PolygonXor{T}([points[f[1]], points[f[2]], points[f[3]]]) for f in faces ]
	return clip(:union, v...)
end

# this is used by 3d->2d slicing:
function glue_segments(points::AbstractVector, segments)
	n = maximum(x[2] for x in segments)
	# compute the set of all edges
	edges = [ SortedSet{eltype(eltype(segments))}() for _ in 1:n ]
	for (a, b) in segments
		push!(edges[a], b)
		push!(edges[b], a)
	end
	# decomposes the graph as a set of loops
	loops = Vector{eltype(points)}[]
	for (a0, e) in pairs(edges)
		isempty(e) && continue
		a = a0; l = [a]
		while !isempty(e)
			b = first(e)
			delete!(edges[a], b)
			delete!(edges[b], a)
			(b == a0) && (length(l) > 1) && break
			push!(l, b)
			e = edges[b]
			a = b
		end
		push!(loops, points[l])
	end
	return PolygonXor{eltype(eltype(points))}(loops)
end


# Intersections (2d)««1

struct HalfPlane{T<:Real}
	# equation a*x+b = 0
	a::SVector{2,T}
	b::T
end
@inline (h::HalfPlane)(v::AbstractVector) = dot(h.a, v) + h.b

function Base.intersect(halfplane::HalfPlane,
		path::AbstractVector{<:StaticVector{2}})
	newpath = similar(path, 0)
	n = length(path)
	s = halfplane.(path) # vector of signs
	for i in 1:n
		j = i+1; (j > n) && (j = 1)
		(si, sj) = (s[i], s[j])
		si >= 0 && push!(newpath, path[i])
		if si*sj < 0
			newpoint = (s[j]*path[i]-s[i]*path[j])/(s[j]-s[i])
			(isempty(newpath) || last(newpath) ≠ newpoint) && push!(newpath, newpoint)
		end
	end
	return newpath
end

@inline function bound(s::PolygonXor)
	r = maximum(norm(v, Inf) for v in vertices(s))
	return r + one(r)
end

@inline function bounding_square(s::PolygonXor)
	r = bound(s)
	return [SA[r,r], SA[-r,r], SA[-r,-r], SA[r,-r]]
end

@inline function Base.intersect(halfplane::HalfPlane, s::PolygonXor)
	c = intersect(halfplane, bounding_square(s))
	return clip(:intersection, PolygonXor{coordtype(s)}(c), s)
end

# Path extrusion««1
# triangulate_between: triangulate between two parallel paths««2
"""
		triangulate_between(poly1, poly2, start1, start2)

Given two polygons `poly1` and `poly2`, both of them represented as a
vector of paths, and produced as offsets from a common path,
find a triangulation for the region between the two polygons.

This functions returns a pair `(triangulation, edge)`, where:

 - the triangulation is a vector of triples of integers,
where each point is represented by its index. Indices in `poly1` start at
value `start1`, and in `poly2` at `start2`.

 - the edge is a pair `(lastidx1, lastidx2)` corresponding to the last
	 points visited on each polygon. (this will be useful for closing the
	 extrusion).

"""
function triangulate_between(poly1, poly2, start1=1, start2=1)
	Big = typemax(eltype(eltype(eltype(poly1))))
	triangles = NTuple{3,Int}[]
	# head is the marker of current leading edge
	# headpoint[i] is the point marked to by head[i]
	# headidx is the new index for this marked point
	# status[i][j] is the number of last used point in j-th path of i-th poly
	head = [(1,1), (1,1)]
	headpoint = [poly1[1][1], poly2[1][1]]
	headidx = [start1, start2]
	status = zeros.(Int,length.((poly1, poly2)))
	# so far we used exactly one point on each side:
	status[1][1] = status[2][1] = 1

	# we need a way to convert (poly, path, index) to integer index««
	function first_indices(start::Int, l::Vector{Int})::Vector{Int}
		f = zeros.(Int, length(l))
		f[1] = start
		for i in 1:length(l)-1
			@inbounds f[i+1] = f[i] + l[i]
		end
		f
	end
	# firstindex[poly][path] is the first index for this path
	# firstindex[1][1] = start1
	# firstindex[1][2] = start1 + len(poly1[1]) etc.
	firstindex = (first_indices(start1, length.(poly1)),
								first_indices(start2, length.(poly2)))
	newindex(poly::Int, path::Int, index::Int)::Int =
		firstindex[poly][path] + index - 1
#»»
	# computing diagonal distances to find the smallest one:««
	distance(pt, path, i) =
		i > length(path) ? Big : distance²(pt, path[i])

	closest(pt, poly, status) =
		findmin([distance(pt, poly[i], status[i]+1) for i in eachindex(poly)])
#»»

	while true
		d1, i1 = closest(headpoint[2], poly1, status[1])
		d2, i2 = closest(headpoint[1], poly2, status[2])
		# if no more points are left, we return:
		(d1 == d2 == Big) && break

		if d1 < d2 # we append a point from poly1
			# add the triangle: head1, head2, newpoint
			s = status[1][i1] += 1
			newidx = newindex(1, i1, s)
			push!(triangles, (headidx[1], headidx[2], newidx))
			# update head1 to point to new point
			headidx[1] = newidx
			head[1] = (i1, s)
			headpoint[1] = poly1[i1][s]
		else
			# add the triangle: head1, head2, newpoint
			s = status[2][i2] += 1
			newidx = newindex(2, i2, s)
			push!(triangles, (headidx[1], headidx[2], newidx))
			# update head1 to point to new point
			headidx[2] = newidx
			head[2] = (i2, s)
			headpoint[2] = poly2[i2][s]
		end
	end
	(triangles, (headidx[1], headidx[2]))
end#»»
# path_extrude««
"""
		path_extrude(path, poly, options...)

Extrudes the given polygon (a path of points forming a simple loop)
along the given path. Both arguments are provided as a
`Vector{SVector{2}}`.

Returns a `Surface` (defined by points and a triangulation).
"""
function path_extrude(path, poly;
		join=:round, closed=false, miter_limit=2.0, precision=0.2)
	# in kwargs: ends = closed ? :fill : :butt
	N = length(poly)
	# offset_path is a vector of vector of paths
	@assert closed == true "Open-path extrusion is currently not implemented by ClipperLib"
	offset_path = Shapes.offset([path], [pt[1] for pt in poly ];
		join, miter_limit, ends=closed ? :fill : :single)
	new_points = SVector{3,eltype(eltype(eltype(offset_path)))}[]
	first_face = Int[]; last_ff = 1
	for (a, paths) in zip(poly, offset_path)
		push!(first_face, last_ff)
		for p in paths
			pxyz = [SA[b[1], b[2], a[2]] for b in p]
			new_points = [new_points; pxyz]
			last_ff+= length(p)
		end
	end
	# new_points is a flat list of all 3d points produced
# 	new_points = [[
# 		[ SA[[pt[1], pt[2], poly[i][2]]] for pt in [p...;] ]
# 		for (i, p) in pairs(offset_path)
# 	]...;]

	# first index for each path
# 	first_face = cumsum([1; # initial
# 		map(p->sum(length.(p)), offset_path)])
# # 	println("first_face=$first_face")

	triangles = map(1:N) do i
		i1 = mod1(i+1, N)
		triangulate_between(offset_path[i], offset_path[i1],
			first_face[i], first_face[i1])
		# XXX keep the last edge for closing the poly
	end
	# this completes the set of triangles for the tube:
	tube_triangles = vcat([ t[1] for t in triangles ]...)
	last_face = [ t[2][1] for t in triangles ]
# 	println("last_face=$last_face")
	# here we decide if it is closed or open
	# if open, triangulate the two facets
	# if closed, join them together
	if closed
		more_triangles = vcat(map(1:N) do i
			j = (i%N)+1
			[ SA[first_face[i], last_face[i], first_face[j]],
				SA[first_face[j], last_face[i], last_face[j]] ]
		end...)
# 		println("more_triangles=$more_triangles")
		tube_triangles = [ tube_triangles; more_triangles ]
	else
	# TODO: triangulate the surface
	# or, for now, close with two non-triangular facets...
		more_triangles = [ reverse(first_face), last_face ]
# 		println("more_triangles=$more_triangles")
	end
	return (new_points, reverse.(tube_triangles))
end#»»
# Exports««1

# export Path
# export PolygonXor
end

