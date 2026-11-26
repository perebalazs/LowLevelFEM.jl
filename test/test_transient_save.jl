# Compare all returned columns and metadata, including nested multifield outputs.
function check_transient_saved(actual, full, indices)
    if full isa Tuple
        @test length(actual) == length(full)
        for (a, f) in zip(actual, full)
            check_transient_saved(a, f, indices)
        end
    else
        @test actual.a == full.a[:, indices]
        @test actual.t == full.t[indices]
        @test actual.nsteps == length(indices) == size(actual.a, 2)
        @test actual.model === full.model
        @test typeof(actual) === typeof(full) && actual.kind === full.kind
    end
end

@testset "Transient: save_every, schemes, backends and callback state" begin
    with_test_gmsh() do
        line_mesh(n=2, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        N = ndofs(P)
        # Consistent matrices exercise factorized mass and boundary cross terms.
        M = ∫(P⋅P)
        K = ∫(Grad(P)⋅Grad(P)) + 2M
        C = 0.3M
        initial = transient_field(P, fill(0.2, N, 1), [0.0])
        velocity = transient_field(P, fill(0.1, N, 1), [0.0])
        for n in (2, 3, 12)
            dt, t0 = 0.005, 0.7
            t = [t0 + (j-1)*dt for j in 1:n]
            load = transient_field(P, repeat(reshape(sin.(3t), 1, :), N, 1), t)
            prescribed = transient_field(P, repeat(reshape(t.^3, 1, :), N, 1), t)
            support = [BoundaryCondition("left"; field=P, p=prescribed)]
            solvers = (:auto, :backslash, :lu, :lu_no_ordering, :cholesky, :ldlt, :qr, :cg, :gmres)
            runs = Any[]
            for solver in solvers
                options = solver in (:cg, :gmres) ?
                    (; solver, preconditioner=:none, solveroptions=(reltol=1e-12, maxiter=100)) : (; solver)
                push!(runs, (; run=(;kw...)->HHT(Symmetric(K), Symmetric(M), Symmetric(C),
                    load, initial, velocity, n, dt; support, t0, α=-0.1, options..., kw...)))
                for ϑ in (0.0, 0.5, 1.0)
                    push!(runs, (; run=(;kw...)->FDM(Symmetric(K), Symmetric(M), load,
                        initial, n, dt; support, t0, ϑ, options..., kw...)))
                end
            end
            for dampingScheme in (:lagged, :central), velocityScheme in (:backward, :central)
                push!(runs, (; run=(;kw...)->CDM(K, M, C, load, initial, velocity, n, dt;
                    support, t0, dampingScheme, velocityScheme, kw...)))
            end
            for entry in runs, modified in (false, true)
                seen = Any[]
                function update_rhs(b, state)
                    # A previous accepted state must be the previous integration
                    # step, even when it is absent from the returned history.
                    push!(seen, (state.step, state.iteration, state.phase, state.t,
                        state.fields[1].a, state.velocity === nothing ? nothing : state.velocity[1].a,
                        state.acceleration === nothing ? nothing : state.acceleration[1].a,
                        state.previous === nothing ? nothing :
                            (state.previous.step, state.previous.fields[1].a)))
                    out = copy(b)
                    out.a .+= 0.01 .* state.fields[1].a.^2
                    if state.previous !== nothing
                        out.a .+= 0.005 .* state.previous.fields[1].a
                    end
                    return out
                end
                hooks = modified ? (; modifyRHS=update_rhs, copyOperators=false) : (;)
                full = entry.run(; hooks...)
                full_seen = copy(seen)
                for (stride, indices) in ((1, collect(1:n)),
                        (4, n == 12 ? [1, 5, 9, 12] : [1, n]),
                        (n+5, [1, n]),
                        (big(typemax(Int))+1, [1, n]))
                    empty!(seen)
                    actual = entry.run(; save_every=stride, hooks...)
                    check_transient_saved(actual, full, indices)
                    @test seen == full_seen
                end
                if n == 12
                    # Divisible interval count: the endpoint must not be duplicated.
                    actual = entry.run(; save_every=11, hooks...)
                    check_transient_saved(actual, full, [1, 12])
                end
                for invalid in (0, -2, 1.5, 2.0, true, false, nothing, "2")
                    @test_throws ArgumentError entry.run(; save_every=invalid, hooks...)
                end
            end
        end
        # Compatibility calling forms and all-constrained systems.
        n, dt = 12, 0.005
        t = collect(0:n-1)*dt
        allbc = [BoundaryCondition("body"; field=P, p=(x,y,z,t)->t^3)]
        indices = [1, 5, 9, 12]
        for method in (HHT, CDM), hooks in ((;), (; modifyRHS=(b,s)->nothing))
            full = method(K, M, initial, allbc, initial, velocity, n, dt;
                          damping=C, hooks...)
            actual = method(K, M, initial, allbc, initial, velocity, n, dt;
                            damping=C, save_every=4, hooks...)
            check_transient_saved(actual, full, indices)
            full = method(K, M, C, initial, allbc, initial, velocity, n, dt; hooks...)
            actual = method(K, M, C, initial, allbc, initial, velocity, n, dt;
                            save_every=4, hooks...)
            check_transient_saved(actual, full, indices)
        end
        for hooks in ((;), (; modifyRHS=(b,s)->nothing))
            full = FDM(K, M, initial, allbc, initial, n, dt; hooks...)
            actual = FDM(K, M, initial, allbc, initial, n, dt; save_every=4, hooks...)
            check_transient_saved(actual, full, indices)
        end
    end
end

@testset "Transient: thinned multifield reconstruction with MPC" begin
    with_test_gmsh() do
        line_mesh(n=2, order=1)
        mat = Material("body")
        P = Field([mat], type=:ScalarField, dim=1, fieldName=:p)
        S = Field([mat], type=:ScalarField, dim=1, fieldName=:s)
        N = ndofs(P)
        Ip = SystemMatrix(spdiagm(0=>ones(N)), P)
        Is = SystemMatrix(spdiagm(0=>ones(N)), S)
        Z = SystemMatrix(spzeros(N, N), S, P)
        K = SystemMatrix([2Ip Z; Z' 3Is])
        M = SystemMatrix([Ip Z; Z' Is])
        C = 0.2M
        zeroP = transient_field(P, zeros(N, 1), [0.0])
        zeroS = transient_field(S, zeros(N, 1), [0.0])
        oneP = transient_field(P, ones(N, 1), [0.0])
        oneS = transient_field(S, ones(N, 1), [0.0])
        initial, velocity = (zeroP, zeroS), (oneP, oneS)
        mpc = [MPC(master="left", slave="right", field=P, p=true)]
        support = [BoundaryCondition("right"; field=P, p=(x,y,z,t)->sin(t))]
        for hooks in ((;), (; modifyRHS=(b,s)->nothing))
            for method in (HHT, CDM)
                full = method(K, M, C, initial, initial, velocity, 12, 0.01; mpc, support, hooks...)
                actual = method(K, M, C, initial, initial, velocity, 12, 0.01;
                                mpc, support, save_every=4, hooks...)
                check_transient_saved(actual, full, [1, 5, 9, 12])
            end
            full = FDM(K, M, initial, initial, 12, 0.01; mpc, support, hooks...)
            actual = FDM(K, M, initial, initial, 12, 0.01; mpc, support, save_every=4, hooks...)
            check_transient_saved(actual, full, [1, 5, 9, 12])
        end
    end
end

@testset "Transient: thinned local and reduced-space output" begin
    for basis in (:coordinates, :reduced, :spectral)
        with_test_gmsh() do
            if basis === :coordinates
                structured_rect_mesh(n=1, order=1)
                P = Field([Material("body")], type=:VectorField, dim=2, fieldName=:u)
                cs = [CoordinateSystem(VectorField(P, "body", [0.0, 1.0, 0.0]))]
                mpc = MPC[]
            else
                line_mesh(n=2, order=3)
                options = basis === :reduced ? (; reducedOrder=true) : (; basis=:spectral)
                P = Field([Material("body")]; type=:ScalarField, dim=1, fieldName=:p, options...)
                cs = CoordinateSystem[]
                mpc = basis === :reduced ? [MPC(master="left", slave="right", field=P, p=true)] : MPC[]
            end
            M = ∫(P⋅P)
            K = ∫(Grad(P)⋅Grad(P)) + 2M
            C = 0.2M
            initial = transient_field(P, zeros(ndofs(P), 1), [0.0])
            velocity = transient_field(P, reshape(collect(1.0:ndofs(P)), :, 1), [0.0])
            for hooks in ((;), (; modifyRHS=(b,s)->nothing))
                for method in (HHT, CDM)
                    full = method(K, M, C, initial, initial, velocity, 12, 0.001; coordSys=cs, mpc, hooks...)
                    actual = method(K, M, C, initial, initial, velocity, 12, 0.001;
                                    coordSys=cs, mpc, save_every=4, hooks...)
                    check_transient_saved(actual, full, [1, 5, 9, 12])
                end
                full = FDM(K, M, velocity, initial, 12, 0.001; coordSys=cs, mpc, hooks...)
                actual = FDM(K, M, velocity, initial, 12, 0.001;
                             coordSys=cs, mpc, save_every=4, hooks...)
                check_transient_saved(actual, full, [1, 5, 9, 12])
            end
        end
    end
end
