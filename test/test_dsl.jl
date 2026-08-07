@testset "Weak-form DSL" begin
    gmsh.initialize()
    try
        # ------------------------------------------------------
        # 3D linear elasticity:
        # DSL assembly must reproduce the legacy stiffness matrix
        # ------------------------------------------------------
        structured_box_mesh(n=2, order=1)

        mat = Material("body")
        Plegacy = Problem([mat])
        Pu = Problem([mat], type=:VectorField, dim=3, field=:u)

        Klegacy = stiffnessMatrix(Plegacy)

        μ = mat.μ
        λ = mat.λ

        D = [
            λ+2μ  λ     λ     0  0  0
            λ     λ+2μ  λ     0  0  0
            λ     λ     λ+2μ  0  0  0
            0     0     0     μ  0  0
            0     0     0     0  μ  0
            0     0     0     0  0  μ
        ]

        Kdsl = ∫(SymGrad(Pu) ⋅ D ⋅ SymGrad(Pu))

        @test size(Kdsl.A) == size(Klegacy.A)
        @test norm(Kdsl.A - Klegacy.A) / norm(Klegacy.A) < 1e-12

        # ------------------------------------------------------
        # Threaded and serial DSL assembly must agree
        # ------------------------------------------------------
        Kserial = ∫(
            SymGrad(Pu) ⋅ D ⋅ SymGrad(Pu);
            multithread=false
        )

        Kthreaded = ∫(
            SymGrad(Pu) ⋅ D ⋅ SymGrad(Pu);
            multithread=true
        )

        @test norm(Kthreaded.A - Kserial.A) / norm(Kserial.A) < 1e-12

    finally
        gmsh.finalize()
    end
end