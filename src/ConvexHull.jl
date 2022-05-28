module ConvexHull
using LinearAlgebra
using StaticArrays

import Polyhedra # for convex hull
import GLPK

module LibTriangle
	using Triangle
end
include("Projectors.jl")
using .Projectors
# using MiniQhull

# tools««1
# findextrema««2
"""
    findextrema(itr; lt=isless)

Like `findmin`, except that
 - it returns both extrema, as a `NamedTuple`;
 - it accepts an `lt` parameter, like `sort`.
"""
function findextrema(itr; lt=isless)
  p = pairs(itr); y = iterate(p)
  if y == nothing
    throw(ArgumentError("collection must be non-empty"))
  end
  (mi, m), s = y; (Mi, M) = (mi, m); i = mi
  while true
    y = iterate(p, s)
    y == nothing && break
    (ai, a), s = y
    if lt(a, m) m = a; mi = ai; end
    if lt(M, a) M = a; Mi = ai; end
  end
  return (min=(m, mi), max=(M, Mi))
end
# determinants««2
@inline det2(x,y) = x[1]*y[2]-x[2]*y[1]
@inline det2(x,y,z) = det2(y-x, z-x)
# triangulations««2
"""
    triangulate_face(points; direction, map)

Returns a triangulation of the face (assumed convex; points in any order)
Optional keyword arguments:
 - `direction` is a normal vector (used for projecting to 2d).
 - `map` is a labeling of points (default is identity map).
 - `convex` is a Val(true) or Val(false).

The triangulation is returned as a vector of StaticVector{3,Int},
containing the labels of three points of each triangle.
"""
function triangulate_face(points::AbstractVector{<:StaticVector{3}};
		direction::AbstractVector = face_normal(points),
		map::AbstractVector{<:Integer} = [1:length(points)...],
		)
	axis = main_axis(direction)
	N = length(points)
	# this common case deserves a shortcut:
	N == 3 && return [(map[1], map[2], map[3])]

	points2d = Matrix{Float64}(undef, N, 2)
	for (i, p) in pairs(points)
		points2d[i,:] .= project2d(axis, p)
	end
	tri = LibTriangle.basic_triangulation(points2d, map)
	return [(t[1], t[2], t[3]) for t in tri]
end
# Polyhedra interface««1
@inline polyhedra_lib(T::Type{<:Real}) =
	Polyhedra.DefaultLibrary{T}(GLPK.Optimizer)

# fixing an oversight in Polyhedra.jl: it has multiplications but no
# divisions
@inline Base.:(/)(h::Polyhedra.HyperPlane, α::Real) =
	Polyhedra.HyperPlane(h.a / α, h.β / α)

# HRepElement is the supertype of HalfSpace and HyperPlane
@inline direction(h::Polyhedra.HRepElement) = h.a
@inline function normalize(h::Polyhedra.HRepElement)
	n = norm(direction(h))
	(n ≠ 0) ? (h / n) : h
end
@inline (h::Polyhedra.HRepElement)(p) = dot(h.a, p) - h.β
@inline ∈(p, h::Polyhedra.HyperPlane) = iszero(h(p))
@inline Base.convert(T::Type{<:Polyhedra.HRepElement},
		h::Polyhedra.HRepElement) = T(h.a, h.β)

# 2d convex hull ««1


"""
    convex_hull_list(points)

Returns the convex hull of the points, as a list of indexes (in direct
order, starting at a reproducible index in the list of points).
"""
function convex_hull_list(points)
  # Uses Graham scan
  # In practice this is faster than using `Polyhedra.jl`.
#   println("points=$points, length=$(length(points))")
  i0 = findextrema(points;
    lt=(p,q)->(p[2]<q[2])|| (p[2]==q[2] && p[1]>q[1])).min[2]
  @inline detp2(i,j,k) = det2(points[[i,j,k]]...)
	# 1024 is an arbitrary value for detecting “aligned” points (i.e. up to
	# representation errors), which should be fast for both Float and Fixed
	# division
  @inline function are_aligned(i,j,k)
    v1 = points[j]-points[i]
    v2 = points[k]-points[j]
    d = det2(v1, v2)
    c = dot(v1, v2)
    return abs(d) < abs(c)/1024
  end
  scan = sort(filter(!isequal(i0), eachindex(points)),
    lt=(p,q)->detp2(i0,p,q) > 0)
#   println("i0=$i0, scan=$scan")
  stack = [i0, scan[1]]
  for h in scan[2:end]
#     println("scanning: $stack + $h")
    v1 = points[stack[end]] - points[stack[end-1]]
    v2 = points[h] - points[stack[end]]
    s = det2(v1, v2)
    c = dot(v1, v2)
    if abs(s) < abs(c)/1024 && c < 0 # points are aligned and backwards
			# here we know that we can insert at (end)
			# look for an insertion point i:
			i = length(stack)
			while i > 2
# 				println(" try to insert at $i")
				v1 = points[stack[i]] - points[stack[i-1]]
				v2 = points[h] - points[stack[i]]
				s = det2(v1, v2)
				c = dot(v1, v2)
				if s < -1e-3*abs(c)
# 					println(" break at $i")
					break
				end
				i -= 1
			end
# 			println("  inserting at $i")
			insert!(stack, i, h)
			continue
# 			println("  now stack=$stack")
    end
    while detp2(last(stack), h, stack[end-1]) < 0
      pop!(stack)
    end
    push!(stack, h)
  end
  return stack
end

"""
    convex_hull([vector of 2d points])

Returns the convex hull (as a vector of 2d points, ordered in direct
order).
"""
function convex_hull(points::AbstractVector{<:StaticVector{2}})
	upoints = unique(sort(points))
	return upoints[convex_hull_list(upoints)]
end

# """
#     convex_hull(x::Geometry{3}...)
# 
# Returns the convex hull of the union of all the given solids, as a
# pair `(points, faces)`. `faces` is a list of triangles.
# """
# @inline convex_hull(x::AbstractGeometry{3}) =
# 	convex_hull(vcat([vertices(y) for y in x]...))

# this is the version using MiniQhull: #««
# convex_hull(points::AbstractVector{<:AnyVec(2)}) =
#		# Delaunay triangulation:
#		let T = MiniQhull.delaunay([p[i] for i in 1:2, p in points]),
#				N = length(points),
#				M = zeros(Bool, N, N) # ⚠ memory O(N²)
#		# mark all edges:
#		for (a,b,c) in eachcol(T)
#			b1 = points[b] - points[a]
#			c1 = points[c] - points[a]
#			d = b1[1]*c1[2] - b1[2]*c1[1] # determinant === orientation of triangle
#			if d < 0
#				M[a,b] = M[b,c] = M[c,a] = true
#			else
#				M[b,a] = M[c,b] = M[a,c] = true
#			end
#		end
#		# list of remaining edges (retrograde oriented)
#		L= sort([(i,j) for i in 1:N, j in 1:N if M[i,j] && ! M[j,i]], by=v->v[1])
#		next(i) = L[searchsorted(L, i, by=y->y[1])][1][2]
#		R = zeros(Int, length(L))
#		R[1:2] .= L[1] # initialize with first edge
#		for i in 3:length(R)
#			R[i] = next(R[i-1])
#		end
#		# returns in retrograde ordering (OpenSCAD convention):
#		points[R]
# end#»»
# this is the version using Polyhedra:««
# """
#     convex_hull([vector of 2d points])
# 
# Returns the convex hull (as a vector of 2d points).
# """
# function convex_hull(points::Union{AbstractVector{<:Vec{2}},
# 	AbstractMatrix{<:Number}})
# # M is a matrix with the points as *columns*, hence the transpose
# # below:
# 	PH = Polyhedra
# 	poly = vpoly(points)
# 	PH.removevredundancy!(poly)
# 	return Vec{2}.(collect(PH.points(poly)))
# end
# »»

# 3d convex hull««1
"""
    convex_hull(vector of 3d points)

Returns the convex hull of these points, as a pair `(points, faces)`.
All the faces are triangles.
"""
function convex_hull(p::AbstractVector{<:StaticVector{3,T}}) where{T}
	M = reduce(hcat, p)
	PH = Polyhedra
	poly = PH.polyhedron(PH.vrep(transpose(M)), polyhedra_lib(T))
	R = PH.removevredundancy!(poly)
	V = SVector{3,T}.(collect(PH.points(poly)))

	triangles = SVector{3,Int}[]
	for i in PH.eachindex(PH.halfspaces(poly)) # index of halfspace
		h = PH.get(poly, i)
		pts = PH.incidentpointindices(poly, i) # vector of indices of points
		vlist = [SVector{3,T}(PH.get(poly, j)) for j in pts]
		length(vlist) ≥ 3 || continue
		for t in triangulate_face( vlist;
				direction = h.a,
				map = [j.value for j in pts])
			(a,b,c) = (V[j] for j in t)
			k = det([b-a c-a h.a])
			push!(triangles, (k > 0) ? t : (t[1], t[3], t[2]))
		end
	end
	return (points=V, faces=triangles)
end

#»»1
export convex_hull
end # module
