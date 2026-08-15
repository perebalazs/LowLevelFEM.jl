# Linear 3D elasticity benchmark for LowLevelFEM.jl.
#
# Run with, for example:
#   julia --threads=4 lowlevelfem_elasticity_benchmark.jl
#
# The script solves a unit cube fixed at x = 0 and loaded by a unit traction
# in the x direction at x = 1. Timings exclude Julia compilation.

using LowLevelFEM
using LinearAlgebra
using Statistics

const n = 30
const repeats = 5

function median_timed(f; repeats = repeats)
    f() # Warm-up: compilation is deliberately excluded.
    samples = map(1:repeats) do _
        GC.gc()
        @timed f()
    end
    return sort(samples; by = s -> s.time)[cld(repeats, 2)]
end

function setup_problem(n)
    structured_box_mesh(n = n)
    gmsh.write("elasticity_benchmark.msh")
    mat = Material("body")
    P = Problem([mat], type = :VectorField, dim = 3, field = :u, bandwidth=:RCMK)
    mu = mat.μ
    lambda = mat.λ
    D = [
        lambda + 2mu  lambda         lambda         0.0  0.0  0.0
        lambda         lambda + 2mu  lambda         0.0  0.0  0.0
        lambda         lambda         lambda + 2mu  0.0  0.0  0.0
        0.0            0.0            0.0            mu   0.0  0.0
        0.0            0.0            0.0            0.0  mu   0.0
        0.0            0.0            0.0            0.0  0.0  mu
    ]
    support = [BoundaryCondition("left", ux = 0.0, uy = 0.0, uz = 0.0)]
    return (; mat, P, D, support)
end

function assemble_and_solve(p)
    K = ∫(SymGrad(p.P) ⋅ p.D ⋅ SymGrad(p.P); threads = :auto)
    f = ∫(p.P ⋅ [1.0, 0.0, 0.0], Γ = "right")
    u = solveField(Symmetric(K), f; support = p.support)
    return (; K, f, u)
end

function stress_checksum(p, u)
    strain = (u ∘ ∇ + ∇ ∘ u) / 2.0
    identity = TensorField(p.P, "body", [1 0 0; 0 1 0; 0 0 1])
    E = p.mat.E
    nu = p.mat.ν
    stress = E / (1 + nu) * (strain + nu / (1 - 2nu) * trace(strain) * identity)
    stress = elementsToNodes(stress)
    return sum(abs, stress.a)
end

function complete_analysis(p)
    # The weak-form integral builds the CSC pattern and assembles its numerical values.
    result = assemble_and_solve(p)
    return (; result..., stress = stress_checksum(p, result.u))
end

function main()
    BLAS.set_num_threads(1)
    mesh = @timed setup_problem(n)
    analysis_stats = median_timed(() -> complete_analysis(mesh.value))
    result = analysis_stats.value

    println("LowLevelFEM.jl linear elasticity benchmark")
    println("Threads: ", Threads.nthreads(), ", mesh parameter n: ", n)
    println("Mesh generation (not timed):    ", round(mesh.time; digits = 4), " s, ", round(mesh.bytes / 2^20; digits = 1), " MiB")
    println("Complete analysis (median):     ", round(analysis_stats.time; digits = 4), " s, ", round(analysis_stats.bytes / 2^20; digits = 1), " MiB")
    println("Includes: CSC pattern, assembly, solve, and stress postprocessing")
    println("Checksums: ", sum(abs, result.K.A), ", ", sum(result.u.a), ", ", result.stress)
end

main()
