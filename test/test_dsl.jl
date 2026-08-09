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

        # ------------------------------------------------------
        # Threaded and serial linear-form assembly must agree
        # ------------------------------------------------------
        fvolume_serial = ∫(
            Pu ⋅ [1.0, 1.0, 1.0];
            multithread=false
        )

        fvolume_threaded = ∫(
            Pu ⋅ [1.0, 1.0, 1.0];
            multithread=true
        )

        @test fvolume_threaded.a ≈ fvolume_serial.a rtol=1e-12 atol=1e-12

        volume_resultant = vec(sum(
            reshape(fvolume_serial.a[:, 1], Pu.pdim, :);
            dims=2
        ))
        @test volume_resultant ≈ [1.0, 1.0, 1.0] rtol=1e-12 atol=1e-12

        ftraction_serial = ∫(
            Pu ⋅ [1.0, 1.0, 1.0];
            Γ="right",
            multithread=false
        )

        ftraction_threaded = ∫(
            Pu ⋅ [1.0, 1.0, 1.0];
            Γ="right",
            multithread=true
        )

        @test ftraction_threaded.a ≈ ftraction_serial.a rtol=1e-12 atol=1e-12

        traction_resultant = vec(sum(
            reshape(ftraction_serial.a[:, 1], Pu.pdim, :);
            dims=2
        ))
        @test traction_resultant ≈ [1.0, 1.0, 1.0] rtol=1e-12 atol=1e-12

        # ------------------------------------------------------
        # Single-field and block solveField paths must agree
        # ------------------------------------------------------
        support = [
            BoundaryCondition(
                "left";
                problem=Pu,
                ux=0.0,
                uy=0.0,
                uz=0.0
            )
        ]

        u = solveField(Kserial, ftraction_serial; support=support)

        fixed = constrainedDoFs(Pu, support)
        free = freeDoFs(Pu, support)
        residual = Kserial.A * u.a[:, 1] - ftraction_serial.a[:, 1]

        @test norm(residual[free]) / norm(ftraction_serial.a[free, 1]) < 1e-10
        @test all(iszero, u.a[fixed, 1])

        Kblock = SystemMatrix(reshape([Kserial], 1, 1))
        Fblock = SystemVector([ftraction_serial])
        ublock = solveField(Kblock, Fblock; support=support)

        @test ublock.a ≈ u.a rtol=1e-12 atol=1e-12

    finally
        gmsh.finalize()
    end
end
