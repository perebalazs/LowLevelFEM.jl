# Run in the LowLevelFEM development environment with the patched sources.
using Test, LowLevelFEM

@testset "Field kind lifecycle" begin
    structured_rect_mesh(n=2)
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
    @test_throws ArgumentError LowLevelFEM.eigenValue(u[1])

    # Concatenation as a time history must not shift eigenvalues.
    @test_throws ArgumentError mergeFields([u, u])
    @test mergeFields([u[1], u[1]]).kind === :field

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
