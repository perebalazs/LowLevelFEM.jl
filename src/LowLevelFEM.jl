module LowLevelFEM

using LinearAlgebra, SparseArrays
using StaticArrays
using Arpack
using JLD2
using Polyester
using PrecompileTools
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

macro disp(expr)
    return :(display("$(string(expr)) = $($expr)"))
end

#=
# --- PRECOMPILE BLOKK ---
@compile_workload begin
    gmsh.initialize()

    # nagyon kicsi 3D próba-probléma (pl. egy 1×1×1 kocka)
    gmsh.open("cube.geo")
    #gmsh.model.add("cube")
    #b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
    #gmsh.model.occ.synchronize()
    #gmsh.model.addPhysicalGroup(3, [1], 10)
    #gmsh.model.setPhysicalName(3, 10, "body")
    ##gmsh.model.mesh.setSize([[0, -1]], 1.0)
    #gmsh.model.mesh.generate(3)

    mat = material("body")
    prob = Problem([mat], type = :Solid)
    openPostProcessor()

    stiffnessMatrix(prob)

    gmsh.finalize()
end
=#

end #module
