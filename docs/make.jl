using Documenter
using KerrGeodesics

makedocs(
    sitename = "KerrGeodesics.jl",
    modules = [KerrGeodesics],
    remotes = nothing,
    format = Documenter.HTML(
        repolink = "https://github.com/CuberYyc808/KerrGeodesics.jl",
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
    versions = ["v0.4.0" => "v0.4.0"],
)
