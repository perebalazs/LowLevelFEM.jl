push!(LOAD_PATH,"../src/")
include("../src/LowLevelFEM.jl")
using Documenter, .LowLevelFEM

#DocMeta.setdocmeta!(LowLevelFEM, :DocTestSetup, :(using LowLevelFEM); recursive=true)

#makedocs(sitename="LowLevelFEM", format=Documenter.LaTeX())
makedocs(;sitename="LowLevelFEM",
    authors="Balázs Pere",
    doctest = false,
    remotes = nothing,
    #repo="https://github.com/perebalazs/LowLevelFEM.jl.git",
    pages=[
        "Home" => "index.md",
    ],
    format=Documenter.HTML(;
        repolink = "https://localhost:8000",
        canonical="https://perebalazs.juliahub.io/LowLevelFEM/stable/",
        #assets=String[],
    ),
)

deploydocs(; 
    repo="github.com/perebalazs/LowLevelFEM.jl.git",
    devbranch = "main",
    devurl="dev",
    target = "build",
    branch = "main",
    versions = ["stable" => "v^", "v#.#" ],
    push_preview = true,
)
