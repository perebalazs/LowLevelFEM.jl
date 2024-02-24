push!(LOAD_PATH,"../src/")
using Documenter
using LowLevelFEM

makedocs(
    sitename = "LowLevelFEM",
    authors="Balázs Pere",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets=String[],
    ),
    pages = [
         "Introduction" => "index.md"],
    repo="https://github.com/perebalazs/LowLevelFEM.jl/blob/{commit}{path}#L{line}",
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
