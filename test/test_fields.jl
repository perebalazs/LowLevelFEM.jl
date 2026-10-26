@testset "Field algebra" begin
    P = Problem([material("dummy")], type=:dummy)

    t = [0.0]

    A = ScalarField(
        [reshape([1.0, 2.0, 3.0], 3, 1)],
        [;;],
        t,
        [1],
        1,
        :scalar,
        P
    )

    B = ScalarField(
        [reshape([2.0, 4.0, 6.0], 3, 1)],
        [;;],
        t,
        [1],
        1,
        :scalar,
        P
    )

    @test isElementwise(A)
    @test isElementwise(B)

    @test (A + B).A[1][:,1] ≈ [3.0, 6.0, 9.0]
    @test (B - A).A[1][:,1] ≈ [1.0, 2.0, 3.0]
    @test (A * B).A[1][:,1] ≈ [2.0, 8.0, 18.0]
    @test (B / A).A[1][:,1] ≈ [2.0, 2.0, 2.0]
    @test (2.0 * A).A[1][:,1] ≈ [2.0, 4.0, 6.0]
    @test (A + 1.0).A[1][:,1] ≈ [2.0, 3.0, 4.0]
end

@testset "Field constructors, conversion and restricted support" begin
    with_test_gmsh() do
        structured_box_mesh(n=1, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=3, fieldName=:p)
        f(x, y, z) = 1.0 + x + 2y + 3z
        nodal = scalarField(P, "body", f)
        elementwise = ScalarField(P, "body", f)
        @test isNodal(nodal)
        @test isElementwise(elementwise)
        @test elementsToNodes(elementwise).a ≈ nodal.a atol=1e-12
        converted = nodesToElements(nodal)
        @test converted.numElem == elementwise.numElem
        @test converted.A ≈ elementwise.A
        @test elementsToNodes(converted).a ≈ nodal.a atol=1e-12

        for (field, components) in ((ScalarField(P, "body", 2.0), 2.0),
                                    (VectorField(P, "body", [1.0, 2.0, 3.0]), [1.0, 2.0, 3.0]),
                                    (TensorField(P, "body", Matrix{Float64}(I, 3, 3)), Matrix{Float64}(I, 3, 3)))
            restricted = elementsToElements(field; onPhysicalGroup="right", fromPhysicalGroup="body")
            @test !isempty(restricted.numElem)
            @test restricted.numElem != field.numElem
            test_constant_field(restricted, components)
            test_constant_field(nodesToElements(elementsToNodes(field)), components)
        end
    end
end

@testset "Field representation and eigen metadata lifecycle" begin
    with_test_gmsh() do
        structured_rect_mesh(n=1)
        mat = Material("body")
        P = Problem([mat], type=:VectorField, dim=2, fieldName=:u)
        Ps = Problem([mat], type=:ScalarField, dim=2, fieldName=:s)
        Pt = Problem([mat], type=:TensorField, dim=2, fieldName=:T)

        t = [2.0, 5.0]
        nv, ns, nt = ndofs(P), ndofs(Ps), ndofs(Pt)

        u = VectorField(Matrix{Float64}[], ones(nv, 2), copy(t), Int[],
                        2, :v2D, P, :eigen_real)
        s = ScalarField(Matrix{Float64}[], ones(ns, 2), copy(t), Int[],
                        2, :scalar, Ps, :eigen_imag)
        T = TensorField(Matrix{Float64}[], ones(nt, 2), copy(t), Int[],
                        2, :tensor, Pt, :eigen_real)

        # A representation-only conversion retains eigenvalue metadata.
        ue = nodesToElements(u)
        @test ue.kind === :eigen_real
        @test elementsToNodes(ue).kind === :eigen_real
        @test ue.t == t

        # copy must be independent, including element data matrices.
        uc = copy(u)
        @test uc.kind === :eigen_real
        @test uc.t == u.t
        uc.a[1, 1] = 99.0
        @test u.a[1, 1] == 1.0
        uec = copy(ue)
        uec.A[1][1, 1] = 99.0
        @test ue.A[1][1, 1] == 1.0

        # Uniform real scalar scaling preserves the eigenvector identity.
        for f in (u, s, T)
            @test (2.0 * f).kind === f.kind
            @test (f * 2.0).kind === f.kind
            @test (f / 2.0).kind === f.kind
            @test (f * 2.0).t == t
        end
        @test (-u).kind === :eigen_real

        # The result of a general field operation is an ordinary field.
        @test (u + u).kind === :field
        @test (u - u).kind === :field
        @test (s + s).kind === :field
        @test (T + T).kind === :field
        @test u[1].kind === :field
        @test (3.0 / s).kind === :field

        # Concatenation as a time history must not shift eigenvalues.

        # JLD2 persistence preserves kind and eigenvalue labels.
        mktempdir() do tmp
            for (i, f) in enumerate((u, s, T))
                path = joinpath(tmp, "kind_$(i)")
                saveField(path, f)
                loaded = loadField(path)
                @test typeof(loaded) === typeof(f)
                @test loaded.kind === f.kind
                @test loaded.t == f.t
                @test loaded.a == f.a
            end
        end
    end
end
