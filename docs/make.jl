pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter
using KerrGeodesics

makedocs(
    sitename = "KerrGeodesics.jl",
    modules = [KerrGeodesics],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API Reference" => "APIs.md",
    ],
    checkdocs = :none,
)

deploydocs(
    repo = "github.com/CuberYyc808/KerrGeodesics.jl.git",
    devbranch = "main",
)
