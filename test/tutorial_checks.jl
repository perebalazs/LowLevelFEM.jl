# Sidecar checks inspect the notebook's actual results after all its code cells.
# These complement execution checks; they do not replace model validation.
using LinearAlgebra

if _tutorial_name == "examples/lin_elast_2D_static.ipynb"
    @testset "Tutorial: analytical plane-stress solution" begin
        @test all(isfinite, u.a)
        @test all(isfinite, u2.a)
        @test probe(u[1], 1.0, 1.0, 0.0) ≈ 1 / mat.E rtol=1e-8
        @test probe(u2[1], 1.0, 1.0, 0.0) ≈ 1 / mat.E rtol=1e-8
        @test probe(u2[2], 1.0, 1.0, 0.0) ≈ -mat.ν / mat.E rtol=1e-8
        @test probe(σx, 0.5, 0.5, 0.0) ≈ 1.0 rtol=1e-8
        @test abs(probe(σy, 0.5, 0.5, 0.0)) < 1e-8
        @test abs(probe(τxy, 0.5, 0.5, 0.0)) < 1e-8
    end
elseif _tutorial_name == "examples/plane-beam-structure.ipynb"
    @testset "Tutorial: coupled beam equilibrium and constraints" begin
        # Postprocessing expands displacement to 3D and converts rotation to
        # element storage; restore the original solver representation here.
        un = projectTo2D(u)
        φn = elementsToNodes(φ)
        values = vcat(un.a[:, 1], φn.a[:, 1])
        @test all(isfinite, values)
        @test length(values) == size(K.A, 1)
        free, fixed, _ = LowLevelFEM.multifield_bc_data(K, [supp_u, supp_φ]; nsteps=1)
        residual = K.A * values - F.a[:, 1]
        scale = max(LinearAlgebra.opnorm(K.A, Inf) * maximum(abs, values) + maximum(abs, F.a), eps())
        @test maximum(abs, residual[free]) / scale < 1e-10
        @test maximum(abs, values[fixed]) < 1e-10
    end
else
    println("NUMERICAL_CHECKS none registered; cell execution was checked only")
end
