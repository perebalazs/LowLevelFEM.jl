module LowLevelFEMPreconditionersExt

using LowLevelFEM
using Preconditioners
using LinearAlgebra
using SparseArrays


function LowLevelFEM._build_preconditioner(
    ::Val{:ich},
    A;
    preconditioneroptions=(;)
)
    # Expand the authoritative triangle rather than exposing ignored storage.
    A0 = A isa Union{Symmetric,Hermitian} ? sparse(A) : A

    memory = get(preconditioneroptions, :memory, 2)

    options = Base.structdiff(
        preconditioneroptions,
        NamedTuple{(:memory,)}
    )

    return CholeskyPreconditioner(
        A0,
        memory;
        options...
    )
end


end # module
