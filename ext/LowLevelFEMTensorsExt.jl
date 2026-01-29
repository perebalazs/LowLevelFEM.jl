module LowLevelFEMTensorsExt

using LowLevelFEM
using Tensors

# Ide később jönnek majd:
# - Tensor/Matrix konverziók
# - Gauss-ponti stress/tangent kernelfüggvények
# - hyperelastic_from_energy implementáció (ha akarod)

function neo_hooke_energy(C::SymmetricTensor{2,3}, μ, K)
    I1 = tr(C)
    J = sqrt(det(C))
    return μ / 2 * (I1 / J^(2 / 3) - 3) + K / 2 * (log(J))^2
end

function PK2_stress(C::SymmetricTensor{2,3}, μ, K)
    return 2 * gradient(C -> neo_hooke_energy(C, μ, K), C)
end

function material_tangent(C::SymmetricTensor{2,3}, μ, K)
    return 4 * hessian(C -> neo_hooke_energy(C, μ, K), C)
end

function LowLevelFEM._tensor_extension_test()
    F = Tensor{2,3}((1.2, 0.3, 0.0,
        0.0, 1.1, 0.0,
        0.0, 0.0, 1.0))

    C = symmetric(F' ⋅ F)
    return tr(C)
end

end # module
