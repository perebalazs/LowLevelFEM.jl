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