push!(LOAD_PATH, "src/")
push!(LOAD_PATH, "../src/")
using Documenter, ConstructiveGeometry
makedocs(
# 	modules = [ConstructiveGeometry],
	sitename="ConstructiveGeometry.jl",
)

deploydocs(
	repo = "github.com/plut/ConstructiveGeometry.jl.git",
)
