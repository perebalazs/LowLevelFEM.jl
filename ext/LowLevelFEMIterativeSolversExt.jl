module LowLevelFEMIterativeSolversExt

using LowLevelFEM
using IterativeSolvers
using LinearAlgebra: norm

LowLevelFEM._require_iterative_solver(::Val{:cg}) = nothing
LowLevelFEM._require_iterative_solver(::Val{:gmres}) = nothing

function solve_with_history(algorithm, solver, A, b;
                            preconditioner=nothing, solveroptions=(;))
    policy = get(solveroptions, :on_failure, :error)
    options = Base.structdiff(solveroptions, NamedTuple{(:on_failure,)})
    pc = preconditioner === nothing ? (;) : (; Pl=preconditioner)
    x, history = algorithm(A, b; log=true, pc..., options...)

    if !history.isconverged
        failure = LowLevelFEM.LinearSolveConvergenceError(
            solver, history.iters, norm(A * x - b)
        )
        if policy === :warn
            @warn sprint(showerror, failure)
        else
            throw(failure)
        end
    end

    return x
end

function LowLevelFEM._solve_iterative(::Val{:cg}, A, b;
                                     preconditioner=nothing, solveroptions=(;))
    return solve_with_history(cg, :cg, A, b; preconditioner, solveroptions)
end

function LowLevelFEM._solve_iterative(::Val{:gmres}, A, b;
                                     preconditioner=nothing, solveroptions=(;))
    return solve_with_history(gmres, :gmres, A, b; preconditioner, solveroptions)
end

end # module
