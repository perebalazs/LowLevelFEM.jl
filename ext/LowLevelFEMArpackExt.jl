module LowLevelFEMArpackExt

using LowLevelFEM
import Arpack

function LowLevelFEM._arpack_eigs(
    A::AbstractMatrix,
    B::AbstractMatrix;
    kwargs...
)
    return Arpack.eigs(A, B; kwargs...)
end

end # module