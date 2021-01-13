push!(LOAD_PATH, "src/")
push!(LOAD_PATH, "../src/")
using Documenter, Solids
makedocs(sitename="Solids.jl")

deploydocs(
    repo = "github.com/plut/Solids.jl.git",
)
