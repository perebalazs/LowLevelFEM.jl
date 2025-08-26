using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using Documenter
using LowLevelFEM

const DOC_PAGES = [
    "Introduction" => "index.md",
    "Functions" => "Functions.md",
    "Examples" => "Examples.md"
]

# HTML dokumentáció
makedocs(
    sitename = "LowLevelFEM",
    authors = "Balázs Pere",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1
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