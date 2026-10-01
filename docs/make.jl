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
        assets = ["assets/custom.css"],
        size_threshold_warn = 128 * 2^10,
        size_threshold = 400 * 2^10,
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "Orbit classes" => "classes.md",
            "Working with an orbit" => "orbits.md",
            "Conventions" => "conventions.md",
            "Numerics and accuracy" => "accuracy.md",
        ],
        "Examples" => "examples.md",
        "APEX and finite-window interfaces" => "interfaces.md",
        "API reference" => "api.md",
    ],
    checkdocs = :exports,
)

# The home page shows the animation of the 56 catalogue orbits from example/.
cp(joinpath(@__DIR__, "..", "example", "animations", "showcase_all.gif"),
    joinpath(@__DIR__, "build", "assets", "showcase_all.gif"); force = true)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/CuberYyc808/KerrGeodesics.jl.git",
        devbranch = "main",
        versions = ["stable" => "v^", "v#.#"],
    )
end
