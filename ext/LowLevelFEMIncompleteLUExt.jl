module LowLevelFEMIncompleteLUExt

using LowLevelFEM
using IncompleteLU
using LinearAlgebra
using SparseArrays


function LowLevelFEM._build_preconditioner(
    ::Val{:ilu},
    A;
    preconditioneroptions=(;)
)
    A0 =
        A isa Union{Symmetric,Hermitian} ?
        parent(A) :
        A

    if !(A0 isa SparseMatrixCSC)
        A0 = sparse(A0)
    end

    return ilu(
        A0;
        preconditioneroptions...
    )
end


end # module