push!(LOAD_PATH, "src/", "../src/")
using Documenter
using ConstructiveGeometry
makedocs(
	pages = [
		"index.md",
		"primitives.md",
		"transformations.md",
		"operations.md",
		"meshing.md",
		"io.md",
	],
# 	modules = [ConstructiveGeometry],
	sitename="ConstructiveGeometry.jl",
	format=Documenter.HTML(prettyurls=get(ENV,"CI",false)),
)

# deploydocs(
# 	repo = "github.com/plut/ConstructiveGeometry.jl.git",
# )
