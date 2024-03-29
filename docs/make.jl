push!(LOAD_PATH, "src/", "../src/")
using Documenter
include("../src/ConstructiveGeometry.jl"); using .ConstructiveGeometry
makedocs(
	pages = [
		"index.md",
		"primitives.md",
		"transformations.md",
		"operations.md",
		"meshing.md",
		"io.md",
		"extending.md",
	],
# 	modules = [ConstructiveGeometry],
	sitename="ConstructiveGeometry.jl",
	format=Documenter.HTML(prettyurls=get(ENV,"CI",nothing) !=nothing),
)

deploydocs(
	repo = "github.com/plut/ConstructiveGeometry.jl.git",
)
