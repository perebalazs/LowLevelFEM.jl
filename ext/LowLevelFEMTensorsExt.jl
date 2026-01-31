module LowLevelFEMTensorsExt

using LowLevelFEM
using Tensors
using StaticArrays

# =========================
# Stress (PK2)
# =========================
function LowLevelFEM._stress_from_energy(
    ψ,
    Cmat::SMatrix{3,3,Float64},
    p
)
    # SMatrix -> Tensor
    C = symmetric(Tensor{2,3}(Cmat))

    # AD: first derivative only
    S = 2 * gradient(C -> ψ(C, p), C)

    # Tensor -> SMatrix
    return @SMatrix [
        S[1,1] S[1,2] S[1,3];
        S[2,1] S[2,2] S[2,3];
        S[3,1] S[3,2] S[3,3]
    ]
end


# =========================
# Tangent (material)
# =========================
function LowLevelFEM._tangent_from_energy(
    ψ,
    Cmat::SMatrix{3,3,Float64},
    p
)
    C = symmetric(Tensor{2,3}(Cmat))

    # AD: second derivative only
    C4 = 4 * hessian(C -> ψ(C, p), C)

    # ---- convert 4th order tensor to Voigt 6x6 SMatrix ----
    return tensor4_to_voigt(C4)
end

const voigt = (
    (1,1),
    (2,2),
    (3,3),
    (2,3),
    (1,3),
    (1,2)
)

@inline function tensor4_to_voigt(C4)
    @SMatrix [
        C4[i,j,k,l] for (i,j) in voigt, (k,l) in voigt
    ]
end


end # module