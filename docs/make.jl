push!(LOAD_PATH,"../src/")
include("../src/LowLevelFEM.jl")
using Documenter, .LowLevelFEM

DocMeta.setdocmeta!(LowLevelFEM, :DocTestSetup, :(using LowLevelFEM); recursive=true)

#makedocs(sitename="LowLevelFEM", format=Documenter.LaTeX())
makedocs(;sitename="LowLevelFEM",
    authors="Balázs Pere",
#    repo="https://github.com/perebalazs/LowLevelFEM.jl/stable/{commit}{path}#{line}",
    format=Documenter.HTML(;
        canonical="https://docs.juliahub.com/General/LowLevelFEM/stable/",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(; repo="github.com/perebalazs/LowLevelFEM.jl.git"#,
    #devbranch = "main",
    #devurl="dev",
    #target = "build",
    #branch = "gh-pages",
    #versions = ["stable" => "v^", "v#.#" ]
)
