using Documenter
using LowLevelFEM

makedocs(
    sitename = "LowLevelFEM",
    format = Documenter.HTML(),
#    modules = [LowLevelFEM],
    pages = [
         "Introduction" => "index.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/perebalazs/LowLevelFEM.jl.git",
    target="build",
    push_preview=false
)
