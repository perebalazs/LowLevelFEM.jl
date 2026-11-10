@testset "Energy material: shared engineering Voigt convention" begin
    # Normal/shear coupling and unequal shear moduli expose both permutations
    # and factors of two that isotropic examples cannot distinguish.
    D0 = Matrix(Diagonal([11.0, 13.0, 17.0, 2.0, 3.0, 5.0]))
    D0[1, 4] = D0[4, 1] = 0.4
    D0[5, 6] = D0[6, 5] = -0.3
    function energy(C, p)
        e = ((C[1,1]-1)/2, (C[2,2]-1)/2, (C[3,3]-1)/2,
             (C[1,2]+C[2,1])/2, (C[2,3]+C[3,2])/2,
             (C[3,1]+C[1,3])/2)
        return p.scale * sum(e[i]*D0[i,j]*e[j] for i in 1:6, j in 1:6)/2
    end
    C0 = LowLevelFEM.SMatrix{3,3,Float64}(
        [1.2 0.04 0.06; 0.04 1.4 0.05; 0.06 0.05 1.6]
    )
    params = (; scale=1.0)
    stress(C) = LowLevelFEM._stress_from_energy(energy, C, params)
    voigt(S) = [S[1,1], S[2,2], S[3,3], S[1,2], S[2,3], S[3,1]]
    expected = D0 * [0.1, 0.2, 0.3, 0.04, 0.05, 0.06]
    @test voigt(stress(C0)) ≈ expected atol=1e-12
    tangent = LowLevelFEM._tangent_from_energy(energy, C0, params)
    @test tangent ≈ D0 atol=1e-12

    # Independent finite differences with respect to engineering Green strain.
    for (j, (a, b)) in enumerate(((1,1),(2,2),(3,3),(1,2),(2,3),(3,1)))
        dC = zeros(3, 3)
        dC[a,b] += 1
        dC[b,a] += 1
        step = 1e-6 * LowLevelFEM.SMatrix{3,3,Float64}(dC)
        numerical = (voigt(stress(C0 + step)) - voigt(stress(C0 - step))) / 2e-6
        @test numerical ≈ tangent[:, j] rtol=1e-8 atol=1e-9
    end

    with_test_gmsh() do
        structured_box_mesh(n=1, order=1)
        U = Field([Material("body")], type=:VectorField, dim=3, fieldName=:u)
        # The legacy total-Lagrangian builder shares the same public convention.
        # At F=I its material tangent must reduce to the modern linear DSL form.
        F = tensorField(U, "body", Matrix{Float64}(I,3,3))
        reference = ∫(SymGrad(U) ⋅ D0 ⋅ SymGrad(U))
        explicit = materialTangentMatrix(U; F, C=D0)
        automatic = materialTangentMatrix(U; F, energy, params)
        @test explicit.A ≈ reference.A rtol=1e-11 atol=1e-11
        @test automatic.A ≈ reference.A rtol=1e-11 atol=1e-11
        for C in (TensorField(U, "body", Matrix(C0)),
                  tensorField(U, "body", Matrix(C0)))
            response = energyMaterial(energy, C; params)
            for (i, (a, b)) in enumerate(((1,1),(2,2),(3,3),(1,2),(2,3),(3,1)))
                test_constant_field(response.S[a,b], expected[i])
                for j in 1:6
                    test_constant_field(response.D[i,j], D0[i,j])
                end
            end
            test_constant_field(IIPiolaKirchhoff(energy, C; params),
                                LowLevelFEM._stress_from_energy(energy, C0, params))
            tangent_only = tangentMatrix(energy, C; params)
            @test all(tangent_only[i,j].A ≈ response.D[i,j].A for i in 1:6, j in 1:6)

            # Use the returned tangent in the DSL without a permutation.
            u = vectorField(U, "body", [(x,y,z)->x+0.2y+0.3z,
                                       (x,y,z)->2y+0.4z, (x,y,z)->3z])
            K = ∫(SymGrad(U) ⋅ response.D ⋅ SymGrad(U))
            strain = [1.0, 2.0, 3.0, 0.2, 0.4, 0.3]
            @test dot(u.a[:,1], K.A*u.a[:,1]) ≈ dot(strain, D0*strain) atol=1e-9
        end

        C = TensorField(U, "body", Matrix(C0))
        scaled = energyMaterial(energy, C; params=(scale=ScalarField(U, "body", 2.0),))
        for i in 1:6, j in 1:6
            test_constant_field(scaled.D[i,j], 2D0[i,j])
        end
    end
end
