push!(LOAD_PATH,"../src/")
using Documenter
using LowLevelFEM

makedocs(
    sitename = "LowLevelFEM",
    authors="Balázs Pere",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        #assets=String[],
        collapselevel=1
    ),
    pages = [
         "Introduction" => "Introduction.md"],
         "Functions" => "index.md"],
         "Examples" => "Examples.md"],
    doctest=false,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/perebalazs/LowLevelFEM.jl.git",
    branch = "gh-pages",
    target="build",
    devbranch="main",
    push_preview=true
)
