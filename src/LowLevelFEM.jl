module LowLevelFEM

using LinearAlgebra, SparseArrays
using IterativeSolvers
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
        mat = material("dummy")
        prob = Problem([mat], type = :dummy)
        stiffnessMatrix(prob)
        solveDisplacement(prob, [], [])
        solveDisplacement(prob, [], [], condensed=true)
        solveDisplacement(prob, [], [], iterative=true)
        solveDisplacement(prob, [], [], condensed=true, iterative=true)
        solveDisplacement(prob, [], [], iterative=true, reltol=0)
        solveDisplacement(prob, [], [], iterative=true, maxiter=0)
        solveDisplacement(prob, [], [], iterative=true, maxiter=0, reltol=0)
        solveDisplacement(prob, [], [], condensed=true, iterative=true, maxiter=0, reltol=0, preconditioner=nothing, ordering=true)
        solveDisplacement(prob, [], [], [])
        solveDisplacement(prob, [], [], [], condensed=true)
        solveDisplacement(prob, [], [], [], iterative=true)
        solveDisplacement(prob, [], [], [], condensed=true, iterative=true)
        solveDisplacement(prob, [], [], [], iterative=true, reltol=0)
        solveDisplacement(prob, [], [], [], iterative=true, maxiter=0)
        solveDisplacement(prob, [], [], [], iterative=true, maxiter=0, reltol=0)
        solveDisplacement(prob, [], [], [], condensed=true, iterative=true, maxiter=0, reltol=0, preconditioner=nothing, ordering=true)
        q = VectorField([], [;;], [], [], 0, :dummy, prob)
        T = ScalarField([], [;;], [], [], 0, :dummy, prob)
        solveStress(q)
        solveStress(q, DoFResults=true)
        solveStress(q, T=T, DoFResults=true)
    end
end

end #module
