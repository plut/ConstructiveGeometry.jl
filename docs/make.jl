push!(LOAD_PATH, "src/")
push!(LOAD_PATH, "../src/")
using Documenter, ConstructiveGeometry
makedocs(sitename="ConstructiveGeometry.jl")

deploydocs(
    repo = "github.com/plut/ConstructiveGeometry.jl.git",
)
