push!(LOAD_PATH,"../src/")
include("../src/LowLevelFEM.jl")
using Documenter, .LowLevelFEM

#makedocs(sitename="LowLevelFEM", format=Documenter.LaTeX())
makedocs(sitename="LowLevelFEM")
