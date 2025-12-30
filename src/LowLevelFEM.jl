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

include("debugtools.jl")
using .DebugTools
include("general.jl")
include("operators.jl")
include("linear.jl")
include("heat.jl")
include("nonlinear.jl")
include("poisson.jl")
#include("fieldtools.jl")
#using .FieldTools
include("extra.jl")

export @showfields, @showstruct, @showdef, @showtype, @showmem, @showmethods, @disp, @showsize
export probe_field

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
        solvePressure(prob, [], [], 0.0; cav=false, periodicSlave="", periodicMaster="")

        tsteps = [1.0]
        sfA = ScalarField([reshape(collect(1.0:3.0), 3, 1)], [;;], tsteps, [1], 1, :scalar, prob)
        sfB = ScalarField([reshape(fill(2.0, 3), 3, 1)], [;;], tsteps, [1], 1, :scalar, prob)
        vfA = VectorField([reshape(collect(1.0:3.0), 3, 1)], [;;], tsteps, [1], 1, :v3D, prob)
        vfB = VectorField([reshape(collect(2.0:4.0), 3, 1)], [;;], tsteps, [1], 1, :v3D, prob)
        tf_block_I = reshape([1.0, 0.0, 0.0,
                              0.0, 1.0, 0.0,
                              0.0, 0.0, 1.0], 9, 1)
        tf_block_diag = reshape([2.0, 0.0, 0.0,
                                 0.0, 3.0, 0.0,
                                 0.0, 0.0, 4.0], 9, 1)
        tfA = TensorField([tf_block_I], [;;], tsteps, [1], 1, :e, prob)
        tfB = TensorField([tf_block_diag], [;;], tsteps, [1], 1, :e, prob)

        _ = sfA + sfB
        _ = sfA - sfB
        _ = sfA * sfB
        _ = sfA / sfB
        _ = sfA * 3.0
        _ = 3.0 * sfA
        _ = sfA / 2.0
        _ = 2.0 / sfA
        _ = log(sfA)
        _ = sqrt(sfA)
        _ = cbrt(sfA)

        _ = vfA + vfB
        _ = vfA - vfB
        _ = vfA * sfA
        _ = vfA / sfA
        _ = dot(vfA, vfB)
        _ = ×(vfA, vfB)
        _ = norm(vfA)
        _ = diagm(vfA)
        _ = ∘(vfA, vfB)

        _ = tfA + tfB
        _ = tfA - tfB
        _ = tfA * tfB
        _ = ⋅(tfA, tfB)
        _ = tfA * vfA
        _ = vfA * tfA
        _ = tfA * sfA
        _ = tfA / sfA
        _ = transpose(tfA)
        _ = adjoint(tfA)
        _ = unitTensor(tfA)
        _ = trace(tfA)
        _ = det(tfA)
        _ = inv(tfA)
        (_, _) = eigen(tfA)
        _ = sqrt(tfA)
        _ = cbrt(tfA)
        _ = log(tfA)
    end
end

end #module
