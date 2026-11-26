@testset "Parameter integration: numbers and arrays" begin
    # No mesh or Gmsh session is needed for numeric integration.
    @test ∫(z -> 3.0, -2, 4) ≈ 18.0
    @test integrate(z -> z^5 + 2z^4 - z^2, -1, 2) ≈
          (2^6-1)/6 + 2*(2^5+1)/5 - (2^3+1)/3
    @test ∫(z -> z, 2, -1; order=1) ≈ -1.5
    @test ∫(z -> z^2, -1, 1; order=2) ≈ 2/3
    @test ∫(z -> [1, z, z^2], -1, 1) ≈ [2, 0, 2/3]
    @test ∫(z -> [1 z; z^2 2], 0, 2) ≈ [2 2; 8/3 4]
    @test ∫(z -> [z 1; 2 z], 1, 1) == zeros(2, 2)
    for bad in (0, -1, true, 2.0, nothing)
        @test_throws ArgumentError ∫(identity, 0, 1; order=bad)
    end
    @test_throws ArgumentError ∫(identity, NaN, 1)
    @test_throws ArgumentError ∫(identity, 0, Inf)
    @test_throws ArgumentError ∫(z -> z < 0 ? [1, 2] : [1], -1, 1; order=3)
    @test_throws ArgumentError ∫(z -> "unsupported", 0, 1)
end

# Compare element tags/values independently of the source ordering.
function check_parameter_field(actual, expected, factor=1.0; offset=0.0)
    @test typeof(actual) === typeof(expected)
    @test isElementwise(actual)
    @test actual.t == expected.t && actual.nsteps == expected.nsteps
    @test actual.model === expected.model && actual.type === expected.type && actual.kind === expected.kind
    @test Set(actual.numElem) == Set(expected.numElem)
    index = Dict(e => i for (i, e) in enumerate(expected.numElem))
    for (i,e) in enumerate(actual.numElem)
        @test actual.A[i] ≈ factor .* expected.A[index[e]] .+ offset atol=1e-12
    end
end

@testset "Parameter integration: fields, metadata and mixed coefficients" begin
    with_test_gmsh() do
        structured_rect_mesh(n=2, order=1)
        mat = Material("body")
        P = Field([mat], type=:ScalarField, dim=2)
        U = Field([mat], type=:VectorField, dim=2)
        s = scalarField(P, "body", (x,y,z)->1+x+2y)
        se = nodesToElements(s)
        # Give the scalar field multiple time columns, preserving metadata.
        sn = ScalarField([], hcat(s.a, 2s.a), [0.4, 0.8], [], 2, s.type, s.model)
        en = nodesToElements(sn)
        v = vectorField(U, "body", [1.0, 2.0, 0.0])
        T = TensorField(U, "body", [1.0 0.2 0.0; 0.3 2.0 0.0; 0.0 0.0 3.0])
        for f in (s, se, sn, en, v, nodesToElements(v), T)
            expected = nodesToElements(f)
            saved = (deepcopy(f.A), copy(f.a), copy(f.t), copy(f.numElem))
            actual = ∫(z -> (1+z^2)*f, -1, 2)
            check_parameter_field(actual, expected, 6.0)
            check_parameter_field(∫(z -> f, 2, -1), expected, -3.0)
            check_parameter_field(∫(z -> f, 1, 1), expected, 0.0)
            @test (f.A, f.a, f.t, f.numElem) == saved
        end
        reversed = ScalarField(reverse(en.A), [;;], en.t, reverse(en.numElem), en.nsteps, en.type, en.model)
        # Same field with different storage/order at different quadrature points.
        check_parameter_field(∫(z -> z < 0 ? sn : reversed, -1, 1; order=3), en, 2.0)
        check_parameter_field(∫(z -> z < 0 ? 3.0 : s, -1, 1; order=3), se, 1.0; offset=3.0)
        check_parameter_field(∫(z -> z < 0 ? s : 3.0, -1, 1; order=3), se, 1.0; offset=3.0)
        mixed = ∫(z -> [z^2*s 0; 2 z*se], 0, 2)
        @test size(mixed) == (2,2)
        @test mixed[1,2] isa Number && mixed[1,2] == 0
        @test mixed[2,1] isa Number && mixed[2,1] ≈ 4
        check_parameter_field(mixed[1,1], se, 8/3)
        check_parameter_field(mixed[2,2], se, 2.0)
        vecresult = ∫(z -> [s, 2z], 0, 2)
        check_parameter_field(vecresult[1], se, 2.0)
        @test vecresult[2] ≈ 4

        # Different time grids, element supports and Problems must not be merged.
        shifted = ScalarField(en.A, [;;], en.t .+ 1, en.numElem, en.nsteps, en.type, en.model)
        partial = ScalarField(se.A[1:end-1], [;;], se.t, se.numElem[1:end-1], se.nsteps, se.type, se.model)
        other = ScalarField(U, "body", 1.0)
        @test_throws ArgumentError ∫(z -> z < 0 ? en : shifted, -1, 1; order=3)
        @test_throws ArgumentError ∫(z -> z < 0 ? se : partial, -1, 1; order=3)
        @test_throws ArgumentError ∫(z -> z < 0 ? se : other, -1, 1; order=3)
        @test_throws ArgumentError ∫(z -> partial, -1, 1; onPhysicalGroup="body")
        # Explicit surface restriction of nodal fields uses the target group.
        edge = nodesToElements(s; onPhysicalGroup="left")
        check_parameter_field(∫(z -> s, 0, 2; onPhysicalGroup="left"), edge, 2.0)
        @test_throws ArgumentError ∫(z -> se, 0, 2; onPhysicalGroup="left")
    end
end

@testset "Parameter integration: layered constitutive matrix and weak forms" begin
    with_test_gmsh() do
        structured_rect_mesh(n=2, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=2)
        c = ScalarField(P, "body", (x,y,z)->2+x+y)
        # Deliberately discontinuous values at shared nodes.
        for (i, A) in enumerate(c.A)
            A .+= i
        end
        Q = ([c 0; 0 4.0], [2c 0; 0 6.0])
        z = [-1.0, 0.0, 2.0]
        A = sum(∫(ζ -> Q[k], z[k], z[k+1]) for k in 1:2)
        B = sum(∫(ζ -> ζ*Q[k], z[k], z[k+1]) for k in 1:2)
        D = sum(∫(ζ -> ζ^2*Q[k], z[k], z[k+1]) for k in 1:2)
        # Analytic zeroth, first and second moment for the two layers.
        for (M, cscale, numeric) in ((A, 5, 16), (B, 3.5, 10), (D, 17/3, 52/3))
            check_parameter_field(M[1,1], c, cscale)
            @test M[1,2] == M[2,1] == 0
            @test M[2,2] ≈ numeric
        end
        section = sum(∫(ζ -> [Q[k] ζ*Q[k]; ζ*Q[k] ζ^2*Q[k]], z[k], z[k+1]) for k in 1:2)
        for i in eachindex(section)
            ref = [A B; B D][i]
            if ref isa ScalarField
                check_parameter_field(section[i], ref)
            else
                @test section[i] ≈ ref
            end
        end
        # Actual compound assembly uses elementwise coefficients directly.
        G = [1.0 0.0; 0.0 1.0; 0.1 0.0; 0.0 0.2] ⋅ Grad(P)
        exact = [5c 0 3.5c 0; 0 16 0 10; 3.5c 0 (17/3)*c 0; 0 10 0 52/3]
        for assembly in (:csc, :ijv)
            @test ∫(G' ⋅ section ⋅ G; Ω="body", assembly).A ≈
                  ∫(G' ⋅ exact ⋅ G; Ω="body", assembly).A atol=1e-11
        end
        term = LowLevelFEM.BilinearTerm(Grad(P), A, Grad(P))
        @test ∫(term + term; Ω="body").A ≈ 2*∫(term; Ω="body").A
        @test ∫(P ⋅ ∫(ζ -> c, 0, 2); Ω="body").a ≈ 2*∫(P ⋅ c; Ω="body").a
        # Nodal coefficients remain usable, and missing elementwise data must
        # never be silently averaged or zero-filled, on either assembly path.
        partial = ScalarField(c.A[1:end-1], [;;], c.t, c.numElem[1:end-1], c.nsteps, c.type, c.model)
        Cp = [partial 0; 0 1.0]
        for assembly in (:csc, :ijv)
            @test_throws ArgumentError ∫(Grad(P) ⋅ Cp ⋅ Grad(P); Ω="body", assembly)
            @test_throws ArgumentError ∫(Grad(P) ⋅ Grad(P); Ω="body", weight=partial, assembly)
            @test_throws ArgumentError ∫(LowLevelFEM.BilinearTerm(Grad(P), Cp, Grad(P)) + term; Ω="body", assembly)
        end
        @test_throws ArgumentError ∫(P ⋅ partial; Ω="body")
        empty = ScalarField([], [;;], c.t, [], c.nsteps, c.type, c.model)
        @test_throws ArgumentError ∫(P ⋅ empty; Ω="body")
        @test_throws ArgumentError ∫(P ⋅ 1; Ω="body", weight=partial)
        @test_throws ArgumentError ∫(P ⋅ c; Γ="left")
        @test ∫(P ⋅ elementsToNodes(c); Γ="left").a ≈
              ∫(P ⋅ nodesToElements(elementsToNodes(c); onPhysicalGroup="left"); Γ="left").a
    end
end

@testset "Parameter integration: embedded shell membrane and bending" begin
    with_test_gmsh() do
        structured_box_mesh(n=1, order=1)
        mat = Material("body")
        U = Field([mat], type=:VectorField, dim=3, fieldName=:u)
        R = Field([mat], type=:VectorField, dim=3, fieldName=:r)
        # A spatially varying constituent on an embedded surface, and a
        # homogeneous layer. The explicit group selects only surface elements.
        q = scalarField(U, "top", (x,y,z)->2+x+y)
        Q = [q 0 0; 0 2.0 0; 0 0 1.5]
        h = 0.2
        section = ∫(ζ -> [Q ζ*Q; ζ*Q ζ^2*Q], -h/2, h/2; onPhysicalGroup="top")
        qe = nodesToElements(q; onPhysicalGroup="top")
        Qe = [qe 0 0; 0 2.0 0; 0 0 1.5]
        reference = [h*Qe zeros(3,3); zeros(3,3) (h^3/12)*Qe]
        I3 = Matrix{Float64}(I, 3, 3)
        Bmb = [I3; zeros(3,3)] ⋅ SurfaceSymGrad(U) +
              [zeros(3,3); I3] ⋅ SurfaceSymGrad(R)
        K = ∫(Bmb' ⋅ section ⋅ Bmb; Γ="top")
        Kref = ∫(Bmb' ⋅ reference ⋅ Bmb; Γ="top")
        @test K.problems == Kref.problems
        @test K.A ≈ Kref.A atol=1e-12
        @test norm(Kref.A) > 0
        @test issymmetric(Kref.A) || Kref.A ≈ Kref.A'
    end
end
