# Linear 3D elasticity benchmark for Gridap.jl.
#
# Run with, for example:
#   julia --threads=4 gridap_elasticity_benchmark.jl
#
# Requires Gridap and GridapGmsh:
#   julia -e 'using Pkg; Pkg.add(["Gridap", "GridapGmsh"])'
#
# The script solves a unit cube fixed at x = 0 and loaded by a unit traction
# in the x direction at x = 1. Timings exclude Julia compilation.

using Gridap
using GridapGmsh
using LinearAlgebra
using Statistics

const repeats = 5
const mesh_file = joinpath(@__DIR__, "elasticity_benchmark.msh")
const E = 2.0e5
const nu = 0.3
const lambda = E * nu / ((1 + nu) * (1 - 2 * nu))
const mu = E / (2 * (1 + nu))

epsilon(u) = 0.5 * (∇(u) + ∇(u)')
const identity_tensor = TensorValue(1.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0)
sigma(eps) = lambda * tr(eps) * identity_tensor + 2 * mu * eps

function median_timed(f; repeats = repeats)
    f() # Warm-up: compilation is deliberately excluded.
    samples = map(1:repeats) do _
        GC.gc()
        @timed f()
    end
    return sort(samples; by = s -> s.time)[cld(repeats, 2)]
end

function setup_problem()
    isfile(mesh_file) || error("Mesh file not found: $mesh_file")
    # The Gmsh physical groups "body", "left", and "right" are read directly.
    return GmshDiscreteModel(mesh_file)
end

function assemble_and_solve(model)
    reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, 1)
    V0 = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "left")
    U = TrialFESpace(V0, VectorValue(0.0, 0.0, 0.0))
    Omega = Triangulation(model)
    dOmega = Measure(Omega, 2)
    Gamma_right = BoundaryTriangulation(model; tags = "right")
    dGamma = Measure(Gamma_right, 2)
    traction = VectorValue(1.0, 0.0, 0.0)

    a(u, v) = ∫(epsilon(v) ⊙ sigma(epsilon(u)))dOmega
    l(v) = ∫(v ⋅ traction)dGamma
    op = AffineFEOperator(a, l, U, V0)
    return solve(op), dOmega
end

function stress_checksum(uh, dOmega)
    stress = sigma(epsilon(uh))
    return sum(∫(stress ⊙ stress)dOmega)
end

function complete_analysis(model)
    # AffineFEOperator constructs the algebraic operator and its sparse pattern.
    uh, dOmega = assemble_and_solve(model)
    return (; uh, stress = stress_checksum(uh, dOmega))
end

function main()
    BLAS.set_num_threads(1)
    mesh = @timed setup_problem()
    analysis_stats = median_timed(() -> complete_analysis(mesh.value))
    result = analysis_stats.value

    println("Gridap.jl linear elasticity benchmark")
    println("Threads: ", Threads.nthreads(), ", cells: ", num_cells(Triangulation(mesh.value)))
    println("Mesh import (not timed):        ", round(mesh.time; digits = 4), " s, ", round(mesh.bytes / 2^20; digits = 1), " MiB")
    println("Complete analysis (median):     ", round(analysis_stats.time; digits = 4), " s, ", round(analysis_stats.bytes / 2^20; digits = 1), " MiB")
    println("Includes: FE spaces, algebraic operator/pattern, assembly, solve, and stress postprocessing")
    println("Checksums: ", sum(get_free_values(result.uh)), ", ", result.stress)
end

main()
