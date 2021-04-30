# OpenSCAD output ««1
"""
    scad(io::IO, name::Symbol, parameters, children)
    scad(filename::AbstractString, s::AbstractGeometry...)
    scad(io::IO, s::AbstractGeometry)

Returns, in printable form (e.g. `Symbol` or `String`), the OpenSCAD name
of this object.

"""
function scad(io::IO, s)
	(name, parameters) = scad_info(s)
	indent(io)
	print(io, name, "(")
	f = true; for (k, v) in pairs(parameters)
		f || print(io, ", "); f = false
		print(io, k, "=", v)
	end
	print(io, ")")
	ch = children(s)
	isempty(ch) && (println(io, ";"); return)

	print(io, " {\n")
	io2 = IOContext(io, :indent => get(io, :indent, 0)+1)
	for c in ch
		scad(io2, c)
	end
	indent(io); print(io, "}")
end
@inline scad(filename::AbstractString, s...) =
	open(filename, "w") do f scad(f, s...) end

# scad_parameters(p::PolygonXor) =
# 	(points = vertices(p),
# 	paths = [ f .- 1 for f in perimeters(p) ])

# special case: Surface, with annotations for points
# function scad(io::IO, s::AbstractSurface)
# 	println(io, "polyhedron(points=[ // ", nvertices(s), " points:")
# 	for (i,p) in pairs(vertices(s))
# 		indent(io)
# 		print(io, " ", Vector{Float64}(coordinates(p)))
# 		if i < nvertices(s) print(io, ","); end
# 		println(io, " // ", i)
# 	end
# 	println(io, "], faces=[ // ", nfaces(s), " faces:")
# 	for (i,f) in pairs(faces(s))
# 		indent(io)
# 		print(io, " ", Vector(f .- 1))
# 		if i < nfaces(s) print(io, ","); end
# 		println(io, " // ", i, "=", Vector(f))
# 	end
# 	indent(io); println(io, "] );")
# end

# function scad(io::IO, s::SetParameters)
# 	indent(io); println(io, "{ // SetParameters")
# 	for (i,j) in pairs(s.data)
# 		indent(io); println(io, "// ", i, "=", j, ";")
# 	end
# 	scad(io, s.child)
# 	indent(io); println(io, "} // SetParameters")
# end

# @inline scad_parameters(io::IO, s::Offset) =
# 	scad_parameters(io, s, Val(parameters(s).join), parameters(s))
# @inline scad_parameters(io::IO, ::Offset, ::Val{:round}, param) =
# 	scad_parameters(io, (r=param.r,))
# @inline scad_parameters(io::IO, ::Offset, ::Val{:miter}, param) =
# 	scad_parameters(io, (delta=param.r, chamfer=false,))
# @inline scad_parameters(io::IO, ::Offset, ::Val{:square}, param) =
# 	scad_parameters(io, (delta=param.r, chamfer=true,))

function strscad(args...)
	b = IOBuffer()
	scad(b, args...)
	return String(take!(b))
end

@inline indent(io::IO) = print(io, " "^get(io, :indent, 0))
# @inline function Base.show(io::IO, l::AbstractGeometry...)
# 	for s in l
# 		indent(io); scad(io, s)
# 	end
# end

# SCAD (for mesh types)««1
function scad(io::IO, m::PolygonXor)
	t = Shapes.identify_polygons(m)
	println(io, "difference() { union() { ")
	for (s, p) in zip(t, Shapes.paths(m)); s > 0 || continue
		println(io, "polygon(", Vector{Float64}.(p), ");")
	end
	println(io, "} union() { ")
	for (s, p) in zip(t, Shapes.paths(m)); s < 0 || continue
		println(io, "polygon(", Vector{Float64}.(p), ");")
	end
	println(io, "} } ")
end
function scad(io::IO, m::CornerTable)
	println(io, "polyhedron(points=",
		Vector{Float64}.(CornerTables.points(m)), ",")
	println(io, "  faces=", [collect(f.-1) for f in CornerTables.faces(m)], ");")
end

# FIXME: mat44 for 4x4 matrices

export scad
