module LowLevelFEM

using LinearAlgebra, SparseArrays
using StaticArrays
using Arpack
using JLD2
using Polyester
using PrecompileTools: @setup_workload, @compile_workload
#using Base.Threads
import gmsh_jll
include(gmsh_jll.gmsh_api)
import .gmsh
export gmsh

macro disp(expr)
    s = string(expr)
    return :(display($s * " = " * string($(esc(expr)))))
end

include("general.jl")
include("operators.jl")
include("linear.jl")
include("heat.jl")
include("nonlinear.jl")

@setup_workload begin
    @compile_workload begin
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        
        # nagyon kicsi 3D próba-probléma (pl. egy 1×1×1 kocka)
        gmsh.model.add("cube")
        b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
        pts = gmsh.model.getEntities(0)
        gmsh.model.occ.mesh.setSize(pts, 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.model.addPhysicalGroup(3, [1], 10)
        gmsh.model.setPhysicalName(3, 10, "body")
        gmsh.model.addPhysicalGroup(2, [1], 11)
        gmsh.model.setPhysicalName(2, 11, "left")
        gmsh.model.addPhysicalGroup(2, [2], 12)
        gmsh.model.setPhysicalName(2, 12, "right")
        
        mat = material("body")
        prob = Problem([mat])
        
        stiffnessMatrix(prob)
        
        gmsh.finalize()
    end
end

end #module
