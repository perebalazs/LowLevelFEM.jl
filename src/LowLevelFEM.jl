module LowLevelFEM

using LinearAlgebra, SparseArrays
using StaticArrays
using Arpack
using JLD2
using Polyester
#using Base.Threads
import gmsh_jll
include(gmsh_jll.gmsh_api)
import .gmsh
export gmsh

include("general.jl")
include("operators.jl")
include("linear.jl")
include("heat.jl")
include("nonlinear.jl")

if Sys.CPU_THREADS != Threads.nthreads()
    @warn "Number of threads ≠ logical threads in CPU"
end

end #module
