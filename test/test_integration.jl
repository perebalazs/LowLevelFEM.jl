@testset "Integration" begin
    gmsh.initialize()
    try
        # ------------------------------------------------------
        # Integration over a unit cube
        # ------------------------------------------------------
        structured_box_mesh(n=2, order=1)

        mat = Material("body")
        P = Problem([mat])

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
        # The requested time step of an elementwise field must be
        # selected during integration
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

    finally
        gmsh.finalize()
    end
end
