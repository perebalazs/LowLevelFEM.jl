# Minimal standalone LowLevelFEM assembly benchmark, essentially the worker part of
# the upstream benchmark_multithread.jl (see ../README.md for the full setup).
# Usage: julia --threads=N llfem_benchmark.jl
# Requires: LowLevelFEM

using LowLevelFEM, SparseArrays, LinearAlgebra

const h = 0.025
const repeats = 5

function make_mesh(h)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("benchmark")
    volume = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()
    group = gmsh.model.addPhysicalGroup(3, [volume])
    gmsh.model.setPhysicalName(3, group, "body")
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)
    gmsh.model.mesh.generate(3)
    return length(gmsh.model.mesh.getNodes()[1]),
        sum(length, gmsh.model.mesh.getElements(3)[2])
end

function median_timed(fn, repeats)
    stats = map(1:repeats) do _
        GC.gc()
        @timed fn()
    end
    return sort(stats; by = s -> s.time)[cld(repeats, 2)]
end

function main()
    BLAS.set_num_threads(1)
    try
        nnodes, nelements = make_mesh(h)
        mat = Material("body") # E = 2.0e5, ν = 0.3 defaults
        problem = Problem([mat]; type = :VectorField, dim = 3, field = :u)
        μ, λ = mat.μ, mat.λ
        D = [
            λ + 2μ  λ       λ       0.0  0.0  0.0
            λ       λ + 2μ  λ       0.0  0.0  0.0
            λ       λ       λ + 2μ  0.0  0.0  0.0
            0.0     0.0     0.0     μ    0.0  0.0
            0.0     0.0     0.0     0.0  μ    0.0
            0.0     0.0     0.0     0.0  0.0  μ
        ]
        println("Threads: ", Threads.nthreads(), ", elements: ", nelements,
            ", nodes: ", nnodes)
        assemble_K = () -> ∫(SymGrad(problem) ⋅ D ⋅ SymGrad(problem))
        assemble_f = () -> ∫(problem ⋅ [1.0, 2.0, 3.0])
        K = assemble_K(); f = assemble_f() # warm-up (exclude compilation)
        for (label, fn) in ["Matrix assembly" => assemble_K, "Vector assembly" => assemble_f]
            s = median_timed(fn, repeats)
            println(label, ": ", round(s.time; digits = 4), " s, ",
                round(s.bytes / 2^20; digits = 1), " MiB")
        end
        println("Checksums: ", sum(abs, K.A), ", ", sum(f.a))
    finally
        gmsh.isInitialized() != 0 && gmsh.finalize()
    end
end

main()
