@testset "Integration and quadrature" begin
    with_test_gmsh() do
        structured_box_mesh(n=2, order=1)

        mat = Material("body")
        P = Field([mat], type=:ScalarField, dim=3, fieldName=:p)

        f(x, y, z) = 1.0 + x + 2.0y + 3.0z
        exact_volume_integral = 4.0

        nodal_field = scalarField(P, "body", f)
        elementwise_field = ScalarField(P, "body", f)

        @test isNodal(nodal_field)
        @test isElementwise(elementwise_field)

        I_function = ∫(P, "body", f; order=2)
        I_nodal = ∫(P, "body", nodal_field; order=2)
        I_elementwise = ∫(P, "body", elementwise_field; order=2)

        @test I_function ≈ exact_volume_integral rtol=1e-12 atol=1e-12
        @test I_nodal ≈ exact_volume_integral rtol=1e-12 atol=1e-12
        @test I_elementwise ≈ exact_volume_integral rtol=1e-12 atol=1e-12
        @test I_nodal ≈ I_function rtol=1e-12 atol=1e-12
        @test I_elementwise ≈ I_function rtol=1e-12 atol=1e-12

        # ------------------------------------------------------
        # The elementwise fallback path must support reordered
        # element tags and values
        # ------------------------------------------------------
        permutation = collect(reverse(eachindex(elementwise_field.numElem)))
        reordered_field = ScalarField(
            elementwise_field.A[permutation],
            [;;],
            elementwise_field.t,
            elementwise_field.numElem[permutation],
            elementwise_field.nsteps,
            elementwise_field.type,
            elementwise_field.model
        )

        I_reordered = ∫(P, "body", reordered_field; order=2)
        @test I_reordered ≈ I_elementwise rtol=1e-12 atol=1e-12

        # ------------------------------------------------------
        # Surface integration uses the same function and nodal
        # field paths on the x = 1 face
        # ------------------------------------------------------
        exact_surface_integral = 4.5

        I_surface_function = integrate(P, "right", f; order=2)
        I_surface_nodal = ∫(P, "right", nodal_field; order=2)

        @test I_surface_function ≈ exact_surface_integral rtol=1e-12 atol=1e-12
        @test I_surface_nodal ≈ exact_surface_integral rtol=1e-12 atol=1e-12

        # ------------------------------------------------------
        # A higher-order rule must be forwarded by the ∫ wrapper
        # ------------------------------------------------------
        quadratic(x, y, z) = x^2 + y^2 + z^2
        @test ∫(P, "body", quadratic; order=3) ≈ 1.0 rtol=1e-12 atol=1e-12

        # ------------------------------------------------------
        # Stored-column selection is a field-data operation; no transient
        # equation or time integrator is exercised here.
        # ------------------------------------------------------
        time_field = ScalarField(
            P,
            "body",
            (x, y, z, t) -> (1.0 + t) * f(x, y, z);
            steps=2,
            tmin=0.0,
            tmax=1.0
        )

        @test ∫(P, "body", time_field; step=1, order=2) ≈ 4.0 rtol=1e-12 atol=1e-12
        @test ∫(P, "body", time_field; step=2, order=2) ≈ 8.0 rtol=1e-12 atol=1e-12
    end
end

@testset "DSL element matrices and assembly storage" begin
    with_test_gmsh() do
        line_mesh(n=1, lx=2.0, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        K = ∫(Grad(P) ⋅ Grad(P); threads=1)
        M = ∫(P ⋅ P)
        f = ∫(P ⋅ 1.0)
        @test Matrix(K.A) ≈ [0.5 -0.5; -0.5 0.5] atol=1e-12
        @test Matrix(M.A) ≈ [2.0 1.0; 1.0 2.0] / 3 atol=1e-12
        @test vec(f.a) ≈ [1.0, 1.0] atol=1e-12
        @test ∫(AxialGrad(P) ⋅ AxialGrad(P)).A ≈ K.A atol=1e-12

        # Nodal coefficients are interpolated at quadrature points: mean(x)=1.
        x = scalarField(P, "body", (x,y,z)->x)
        @test ∫(Grad(P) ⋅ x ⋅ Grad(P)).A ≈ K.A atol=1e-12
        @test ∫(Grad(P) ⋅ 2.0 ⋅ Grad(P)).A ≈ 2K.A atol=1e-12
        @test ∫((Grad(P) ⋅ Grad(P)) * 3.0).A ≈ 3K.A atol=1e-12
        @test ∫(Grad(P) ⋅ Grad(P); assembly=:ijv).A ≈ K.A atol=1e-12

        pattern = build_csc_pattern(P, P; Ω="body")
        first = ∫(Grad(P) ⋅ Grad(P); csc_matrix=pattern)
        reference = copy(first.A)
        fill!(pattern.nzval, 0.0)
        second = ∫(Grad(P) ⋅ 2.0 ⋅ Grad(P); csc_matrix=pattern)
        @test second.A === pattern
        @test second.A ≈ 2reference atol=1e-12
        ∫(Grad(P) ⋅ Grad(P); csc_matrix=pattern)
        @test pattern ≈ 3reference atol=1e-12
    end
end

@testset "DSL elasticity, loads and parallel assembly" begin
    with_test_gmsh() do
        structured_box_mesh(n=2, order=1)
        mat = Material("body")
        U = Field([mat], type=:VectorField, dim=3, fieldName=:u)
        C = D(:Solid, mat)
        expr = SymGrad(U) ⋅ C ⋅ SymGrad(U)
        K = ∫(expr; threads=1)
        Kauto = ∫(expr; threads=:auto)
        @test relerr(Kauto.A, K.A) < 1e-12
        @test relerr(∫(expr; assembly=:ijv).A, K.A) < 1e-12
        @test relerr(K.A', K.A) < 1e-12

        u = vectorField(U, "body", [(x,y,z)->x+2y, (x,y,z)->3y, 0.0])
        strain = [1.0, 3.0, 0.0, 2.0, 0.0, 0.0]
        @test dot(u.a[:, 1], K.A * u.a[:, 1]) ≈ dot(strain, C * strain) rtol=1e-12
        translation = vectorField(U, "body", [1.0, -2.0, 3.0])
        @test norm(K.A * translation.a) / norm(K.A) < 1e-12

        for domain in (nothing, "right")
            fserial = domain === nothing ? ∫(U ⋅ [1.0, 1.0, 1.0]; threads=1) :
                                          ∫(U ⋅ [1.0, 1.0, 1.0]; Γ=domain, threads=1)
            fauto = domain === nothing ? ∫(U ⋅ [1.0, 1.0, 1.0]; threads=:auto) :
                                        ∫(U ⋅ [1.0, 1.0, 1.0]; Γ=domain, threads=:auto)
            @test fauto.a ≈ fserial.a rtol=1e-12 atol=1e-12
            @test vec(sum(reshape(fserial.a[:, 1], 3, :); dims=2)) ≈ ones(3) rtol=1e-12 atol=1e-12
        end

        S = Field([mat], type=:ScalarField, dim=3, fieldName=:s)
        sf = scalarField(S, "body", (x,y,z)->2y+3z)
        Ksurface = ∫(SurfaceGrad(S) ⋅ SurfaceGrad(S); Γ="right")
        @test dot(sf.a[:, 1], Ksurface.A * sf.a[:, 1]) ≈ 13.0 atol=1e-10
    end
end

@testset "CSC reuse rejects incomplete vector component patterns" begin
    with_test_gmsh() do
        structured_rect_mesh(nx=2, ny=1, order=1)
        U = Field([Material("body")]; type=:VectorField, dim=2, fieldName=:u)
        term = Grad(U)⋅Grad(U)
        B = [1.0 0 0 0; 0 1.0 0 0; 0 0 1.0 0; 0 0 0 1.0]
        compound = Grad(U)⋅B⋅Grad(U)
        reference = ∫(term)
        pattern = build_csc_pattern(U, U)
        actual = ∫(term; csc_matrix=copy(pattern))
        @test actual.A ≈ reference.A atol=1e-12
        actual = ∫(compound; csc_matrix=copy(pattern))
        @test actual.A ≈ reference.A atol=1e-12

        # Unequal column lengths inside one trial node's component block.
        unequal = copy(pattern)
        unequal.nzval .= 1.0
        unequal[2, 1] = 0.0
        dropzeros!(unequal)
        # Equal lengths but different row patterns in the component columns.
        different = sparse(ones(ndofs(U), ndofs(U)))
        different[3, 1] = different[4, 1] = 0.0
        different[5, 2] = different[6, 2] = 0.0
        dropzeros!(different)
        # Identical column patterns, but an incomplete test-component block.
        partial = sparse(ones(ndofs(U), ndofs(U)))
        partial[2, 1] = partial[2, 2] = 0.0
        dropzeros!(partial)
        for bad in (unequal, different, partial), expr in (term, compound)
            saved = copy(bad.nzval)
            @test_throws ArgumentError ∫(expr; csc_matrix=bad, threads=1)
            @test bad.nzval == saved
        end
    end
end
