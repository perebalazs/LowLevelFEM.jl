@testset "Operators" begin
    gmsh.initialize()
    try
        if get(ENV, "CI", "false") == "true"
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.Verbosity", 0)
        end

        structured_rect_mesh()

        mat = material("body")
        P = Problem([mat], type=:ScalarField, dim=2)

        c = ScalarField(P, "body", 3.0)
        g = grad(c)

        @test isElementwise(g)
        @test all(all(abs.(A) .< 1e-10) for A in g.A)

    finally
        gmsh.finalize()
    end
end