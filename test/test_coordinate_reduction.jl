@testset "Reduced-order vector fields in varying local coordinates" begin
    with_test_gmsh() do
        structured_rect_mesh(n=2, order=2)
        mat = Material("body")
        U = Field([mat]; type=:VectorField, dim=2, field=:u, reducedOrder=true)
        angle(x, y, z) = 0.2 + 0.7x + 0.4y
        e1 = VectorField(U, "body", [(x,y,z) -> cos(angle(x,y,z)),
                                     (x,y,z) -> sin(angle(x,y,z)), 0.0])
        coordSys = [CoordinateSystem(e1)]
        Q = rotateNodes(U, only(coordSys)).T
        T, R = reductionMatrices(U)
        L = Q * T
        local_u = vectorField(U, "body", [(x,y,z) -> 1+x+2y,
                                          (x,y,z) -> 0.5-y, 0.0])
        exact = Q * local_u.a
        K = ∫(Grad(U) ⋅ Grad(U)) + 2.0 * ∫(U ⋅ U)
        M = ∫(U ⋅ U)
        support = [BoundaryCondition("left"; field=U, ux=(x,y,z) -> 1+2y)]
        f = K * LowLevelFEM._transient_field(U, exact, [0.0])

        @testset "Static single and coupled blocks" begin
            u = solveField(K, Global(f); support, coordSys)
            @test u.a ≈ exact rtol=1e-11
            @test rotateResult(u; toCoordSys=coordSys).a ≈ local_u.a rtol=1e-11
            @test norm(exact - T * (R * exact)) > 1e-3

            P = Field([mat]; type=:ScalarField, dim=2, field=:p)
            C = ∫(P ⋅ [0.1 -0.2] ⋅ U)
            A = SystemMatrix([K C'; C ∫(P ⋅ P)])
            p = scalarField(P, "body", (x,y,z) -> 2+x-y)
            F = SystemVector([f + C' * p, C * u + ∫(P ⋅ P) * p])
            mpc = [MPC(master="left", slave="right", field=U, ux=false, uy=true)]
            uc, pc = solveField(A, Global(F); support, coordSys, mpc)
            @test uc.a ≈ exact rtol=1e-11
            @test pc.a ≈ p.a rtol=1e-11
        end

        @testset "Eigenmodes" begin
            fixed = constrainedDoFs(U, support)
            fixed_r = findall(!iszero, vec(sum(abs.(R[:, fixed]); dims=2)))
            Lfree = L[:, setdiff(axes(L, 2), fixed_r)]
            Kr, Mr = Lfree' * K.A * Lfree, Lfree' * M.A * Lfree
            expected = eigvals(Symmetric(Matrix(Kr)), Symmetric(Matrix(Mr)))[1:3]
            modes = solveEigen(K, M; n=3, solver=:direct, support, coordSys)
            @test modes.t ≈ expected rtol=1e-10
            @test norm(Lfree' * (K.A * modes.a - M.A * modes.a * Diagonal(modes.t))) < 1e-9
        end

        @testset "Transient Cartesian input and local supports" begin
            n, dt = 9, 0.01
            t = collect(0:n-1) .* dt
            Uhist = vec(exact) * (1 .+ t)'
            Vhist = repeat(exact, 1, n)
            C = 0.3 * M
            U0 = LowLevelFEM._transient_field(U, exact, [0.0])
            load = LowLevelFEM._transient_field(U, K.A * Uhist + C.A * Vhist, t)
            bc = [BoundaryCondition("left"; field=U, ux=(x,y,z,t) -> (1+2y)*(1+t))]
            for method in (HHT, CDM)
                u, v = method(K, M, C, Global(load), Global(U0), Global(U0), n, dt;
                              support=bc, coordSys)
                @test u.a ≈ Uhist rtol=1e-10
                @test v.a ≈ Vhist rtol=1e-9
            end
            u = FDM(K, C, Global(load), Global(U0), n, dt; support=bc, coordSys)
            @test u.a ≈ Uhist rtol=1e-10
        end
    end
end
