@testset "Contact geometry and DSL operators" begin
    for dim in (2, 3)
        @testset "Flat interface in $(dim)D" begin
            with_test_gmsh() do
                (; U, C, zero_u, gap) = flat_contact_fixture(; dim=dim)
                @test !isempty(C.slave_nodes)
                @test C.gap_values ≈ fill(gap, length(C.slave_nodes)) atol=1e-10
                @test !any(C.active)
                @test (C.t2 === nothing) == (dim == 2)
                test_constant_field(C.n, [0.0, 1.0, 0.0])
                test_constant_field(norm(C.t1), 1.0)
                test_constant_field(dot(C.n, C.t1), 0.0)
                if dim == 3
                    test_constant_field(norm(C.t2), 1.0)
                    test_constant_field(dot(C.n, C.t2), 0.0)
                    test_constant_field(dot(C.t1, C.t2), 0.0)
                end

                # Open contact contributes no active normal stiffness.
                cn = 2.5
                Gn = ContactGap(C)
                Kopen = ∫(Gn ⋅ cn ⋅ Gn)
                @test norm(Kopen.A) < 1e-12

                Gall = ContactGap(C; active=:all)
                Kall = ∫(Gall ⋅ cn ⋅ Gall; threads=1)
                @test norm(Kall.A) > 0
                @test Kall.A ≈ Kall.A' atol=1e-12
                @test eigmin(Symmetric(Matrix(Kall.A))) > -1e-10
                @test ∫(Gall ⋅ cn ⋅ Gall; assembly=:ijv).A ≈ Kall.A atol=1e-12
                @test ∫(Gall ⋅ cn ⋅ Gall; threads=:auto).A ≈ Kall.A atol=1e-12

                translation = vectorField(U, "body", [0.2, -0.3, 0.1])
                @test norm(Kall.A * translation.a) < 1e-10
                updateContact!(C, translation)
                @test C.gap_values ≈ fill(gap, length(C.slave_nodes)) atol=1e-10
                @test !any(C.active)

                # Relative normal displacement of one on a unit interface has
                # penalty energy cn * area = cn; tangential motion has zero energy.
                normal_jump = vectorField(U, "body", [0.0, (x,y,z)->y>1.0 ? 1.0 : 0.0, 0.0])
                tangent_jump = vectorField(U, "body", [(x,y,z)->y>1.0 ? 1.0 : 0.0, 0.0, 0.0])
                @test dot(normal_jump.a[:, 1], Kall.A * normal_jump.a[:, 1]) ≈ cn atol=1e-10
                @test norm(Kall.A * tangent_jump.a) < 1e-10

                updateContact!(C, zero_u)
                Gfull = ContactGap(C; components=:all, active=:all)
                Kfull = ∫(Gfull ⋅ ContactStiffness(C; cn=cn, ct=0.0) ⋅ Gfull)
                @test Kfull.A ≈ Kall.A atol=1e-10

                if dim == 2
                    # Independent hand-assembled two-node interface matrix.
                    master = LowLevelFEM._nodesOnPhysicalGroup(U, "master")
                    slave = LowLevelFEM._nodesOnPhysicalGroup(U, "slave")
                    by_x(node) = gmsh.model.mesh.getNode(node)[1][1]
                    sort!(master; by=by_x)
                    sort!(slave; by=by_x)
                    dofs = 2 .* vcat(slave, master)
                    Me = [2.0 1.0; 1.0 2.0] / 6
                    expected = cn * [Me -Me; -Me Me]
                    @test Matrix(Kall.A[dofs, dofs]) ≈ expected atol=1e-10
                end

                # A fixed prescribed penetration is a static geometry update,
                # with an exact gap, penalty energy and normal force resultant.
                penetration = 0.02
                u = -(gap + penetration) * normal_jump
                @test updateContact!(C, u) === C
                @test all(C.active)
                @test C.gap_values ≈ fill(-penetration, length(C.slave_nodes)) atol=1e-10
                Kactive = ∫(Gn ⋅ cn ⋅ Gn)
                @test Kactive.A ≈ Kall.A atol=1e-10
                gap_field = ContactGap(C, u)
                gap_nodes = elementsToNodes(gap_field)
                @test gap_nodes.a[C.slave_nodes, 1] ≈ fill(-penetration, length(C.slave_nodes)) atol=1e-10
                r = vectorField(U, "body", [(x,y,z)->x, (x,y,z)->y, (x,y,z)->z])
                current = r.a[:, 1] + u.a[:, 1]
                @test dot(current, Kactive.A * current) ≈ cn * penetration^2 atol=1e-10
                force = Kactive.A * current
                slave_dofs = dim .* (C.slave_nodes .- 1) .+ 2
                @test sum(force[slave_dofs]) ≈ -cn * penetration atol=1e-10
                @test norm(vec(sum(reshape(force, dim, :); dims=2))) < 1e-10

                # Multiplier coupling enters the DSL directly, as a rectangular
                # test/trial matrix. No legacy contact matrices are accessed.
                Λ = Field([Material("body")], type=:ScalarField, dim=dim, fieldName=:λ)
                B = ∫(Λ ⋅ Gn)
                @test size(B.A) == (ndofs(Λ), ndofs(U))
                @test B.model === U
                @test B.test_model === Λ
                @test sum(B.A * current) ≈ -penetration atol=1e-10

                # updateFrom refreshes integration geometry without replacing
                # the Contact object captured by the already-built operator.
                reopened = ∫(Gn ⋅ cn ⋅ Gn; updateFrom=zero_u)
                @test norm(reopened.A) < 1e-10
                updateContact!(ContactSet(C), zero_u)
                @test !any(C.active)
            end
        end
    end
end

@testset "Static contact: analytical compression with finite penalty" begin
    with_test_gmsh() do
        (; U, C, gap) = flat_contact_fixture()
        cn, penetration = 4.0, 0.02
        # E=1, ν=0 gives a unit normal stiffness for the lower unit square.
        K = ∫(SymGrad(U) ⋅ Diagonal([1.0, 1.0, 0.5]) ⋅ SymGrad(U))
        G = ContactGap(C; active=:all)
        Kc = ∫(G ⋅ cn ⋅ G)
        r = vectorField(U, "body", [(x,y,z)->x, (x,y,z)->y, (x,y,z)->z])
        support = [BoundaryCondition("body"; field=U, ux=0.0),
                   BoundaryCondition("bottom"; field=U, uy=0.0),
                   BoundaryCondition("upper"; field=U, uy=-(gap+penetration))]
        u = solveField(K + Kc, -Kc * r; support=support)
        updateContact!(C, u)
        master = LowLevelFEM._nodesOnPhysicalGroup(U, "master")
        # Equilibrium: w = cn*(-penetration-w), so w=-cn*penetration/(1+cn).
        expected_w = -cn * penetration / (1 + cn)
        expected_gap = -penetration / (1 + cn)
        @test u.a[2 .* master, 1] ≈ fill(expected_w, length(master)) atol=1e-10
        @test C.gap_values ≈ fill(expected_gap, length(C.slave_nodes)) atol=1e-10
        @test all(C.active)
        free = freeDoFs(U, support)
        residual = (K.A + Kc.A) * u.a[:, 1] + Kc.A * r.a[:, 1]
        @test norm(residual[free]) < 1e-10
        gap_gp = ContactGap(C, u; gauss=:full)
        @test gap_gp.a[C.slave_nodes, 1] ≈ fill(expected_gap, length(C.slave_nodes)) atol=1e-10
    end
end
