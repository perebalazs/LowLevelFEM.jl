@testset "Minimal cube problem" begin
    # Write your tests here.
    #@test LowLevelFEM.greet_your_package_name() == "Hello LowLevelFEM!"
    #@test LowLevelFEM.greet_your_package_name() != "Hello world!"
    #returntest = "test", 0, 0, 0
    #@test FEM.displacementConstraint("test", ux=0, uy=0, uz=0) == returntest
    #display(FEM.displacementConstraint("test", ux=0, uy=0, uz=0))
    try
        gmsh.initialize()
        if get(ENV, "CI", "false") == "true"
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.Verbosity", 0)
        end
        structured_box_mesh()

        mat = material("body")
        prob = Problem([mat])

        supp = displacementConstraint("left", ux=0, uy=0, uz=0)
        load1 = load("right", fx=1)
        u = solveDisplacement(prob, load=[load1], support=[supp])
        S = solveStress(u)
        σ = probe(S, 0.5, 0.5, 0.5)
        @test all(!isnan, σ)
    finally
        gmsh.finalize()
    end
end
