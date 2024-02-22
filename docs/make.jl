push!(LOAD_PATH,"../src/")
include("../src/LowLevelFEM.jl")
using Documenter, .LowLevelFEM

DocMeta.setdocmeta!(LowLevelFEM, :DocTestSetup, :(using LowLevelFEM); recursive=true)

#makedocs(sitename="LowLevelFEM", format=Documenter.LaTeX())
makedocs(;sitename="LowLevelFEM",
    authors="Balázs Pere",
    repo="https://github.com/perebalazs/LowLevelFEM.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://perebalazs.github.io/LowLevelFEM.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(; repo="github.com/perebalazs/LowLevelFEM.jl",
    devbranch = "main",
    devurl="dev",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ]
)
