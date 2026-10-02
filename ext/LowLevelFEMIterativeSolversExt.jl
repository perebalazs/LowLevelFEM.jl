module LowLevelFEMIterativeSolversExt

using LowLevelFEM
using IterativeSolvers


LowLevelFEM._require_iterative_solver(::Val{:cg}) = nothing
LowLevelFEM._require_iterative_solver(::Val{:gmres}) = nothing


function LowLevelFEM._solve_iterative(
    ::Val{:cg},
    A,
    b;
    preconditioner=nothing,
    solveroptions=(;)
)
    if preconditioner === nothing
        x, history = cg(
            A,
            b;
            log=true,
            solveroptions...
        )
    else
        x, history = cg(
            A,
            b;
            Pl=preconditioner,
            log=true,
            solveroptions...
        )
    end

    if !history.isconverged
        @warn "CG did not converge after $(history.iters) iterations."
    end

    return x
end


function LowLevelFEM._solve_iterative(
    ::Val{:gmres},
    A,
    b;
    preconditioner=nothing,
    solveroptions=(;)
)
    if preconditioner === nothing
        x, history = gmres(
            A,
            b;
            log=true,
            solveroptions...
        )
    else
        x, history = gmres(
            A,
            b;
            Pl=preconditioner,
            log=true,
            solveroptions...
        )
    end

    if !history.isconverged
        @warn "GMRES did not converge after $(history.iters) iterations."
    end

    return x
end


end # module