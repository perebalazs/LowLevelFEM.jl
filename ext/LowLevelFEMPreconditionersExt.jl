module LowLevelFEMPreconditionersExt

using LowLevelFEM
using Preconditioners
using LinearAlgebra


function LowLevelFEM._build_preconditioner(
    ::Val{:ich},
    A;
    preconditioneroptions=(;)
)
    A0 = A isa Union{Symmetric,Hermitian} ? parent(A) : A

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