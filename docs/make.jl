using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using Documenter
using LowLevelFEM

const DOC_PAGES = [
    "Introduction" => "index.md",
    "API Reference" => [
        "General" => "General.md",
        "Linear" => "Linear.md",
        "Heat" => "Heat.md",
        "Nonlinear" => "Nonlinear.md",
        "Operators" => "Operators.md",
    ],
    "Examples" => "Examples.md",
]

# HTML dokumentáció
makedocs(
    sitename = "LowLevelFEM",
    authors = "Balázs Pere",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1,
        sidebar_sitename = true,
        size_threshold_warn = 250.0,  # KiB
        size_threshold = 400.0,       # KiB
    ),
    pages = DOC_PAGES,
    doctest = false,
)

# LaTeX (PDF) dokumentáció
#@info "STARTING LATEX BUILD"
#makedocs(
#    sitename = "LowLevelFEM",
#    authors = "Balázs Pere",
#    format = Documenter.LaTeX(),
#    pages = DOC_PAGES,
#    doctest = false,
#)

deploydocs(
    repo = "github.com/perebalazs/LowLevelFEM.jl.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    push_preview = true
)
