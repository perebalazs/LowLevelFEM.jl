using Pkg
using LowLevelFEM
using LinearAlgebra
using SparseArrays

const WORKER_FLAG = "--worker"
const REPEATS = parse(Int, get(ENV, "LLFEM_BENCH_REPEATS", "5"))
const MESH_SIZE = parse(Float64, get(ENV, "LLFEM_BENCH_H", "0.025"))

function thread_counts()
    nmax = Sys.CPU_THREADS
    counts = Int[1]

    n = 2
    while n < nmax
        push!(counts, n)
        n *= 2
    end

    nmax ∉ counts && push!(counts, nmax)
    return counts
end

function run_workers()
    script = abspath(@__FILE__)
    project = dirname(Base.active_project())
    julia = Base.julia_cmd()

    println("LowLevelFEM multithread benchmark")
    println("Launching independent Julia processes...")
    println()

    for nt in thread_counts()
        cmd = `$julia --project=$project --threads=$nt $script $WORKER_FLAG`
        run(cmd)
    end
end

function package_information()
    version = try
        string(Base.pkgversion(LowLevelFEM))
    catch
        "unknown"
    end

    source = try
        pathof(LowLevelFEM)
    catch
        "unknown"
    end

    return version, source
end

function median_timed(f, repeats)
    times = Float64[]
    bytes = Int[]
    result = nothing

    for _ in 1:repeats
        GC.gc()
        stats = @timed f()
        result = stats.value
        push!(times, stats.time)
        push!(bytes, stats.bytes)
    end

    order = sortperm(times)
    middle = order[cld(length(order), 2)]

    return result, times[middle], bytes[middle]
end

function make_mesh(h)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("multithread_benchmark")

    volume = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()

    group = gmsh.model.addPhysicalGroup(3, [volume])
    gmsh.model.setPhysicalName(3, group, "body")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)
    gmsh.model.mesh.generate(3)

    node_tags = gmsh.model.mesh.getNodes()[1]
    element_tags = gmsh.model.mesh.getElements(3)[2]

    return length(node_tags), sum(length, element_tags)
end

function worker()

    LinearAlgebra.BLAS.set_num_threads(1)

    version, tree = package_information()

    println("="^72)
    println("Threads:              ", Threads.nthreads())
    println("Julia:                ", VERSION)
    println("OS / architecture:    ", Sys.KERNEL, " / ", Sys.ARCH)
    println("CPU:                  ", Sys.cpu_info()[1].model)
    println("Logical CPU threads:  ", Sys.CPU_THREADS)
    println("Total memory:         ", round(Sys.total_memory() / 2.0^30; digits=1), " GiB")
    println("BLAS threads:         ", LinearAlgebra.BLAS.get_num_threads())
    println("LowLevelFEM version:  ", version)
    println("Package tree hash:    ", tree)
    println("Mesh size parameter:  ", MESH_SIZE)
    println("Repetitions:          ", REPEATS)

    try
        nnodes, nelements = make_mesh(MESH_SIZE)

        mat = Material("body")
        problem = Problem(
            [mat];
            type=:VectorField,
            dim=3,
            field=:u,
        )

        μ = mat.μ
        λ = mat.λ

        D = [
            λ + 2μ  λ       λ       0.0  0.0  0.0
            λ       λ + 2μ  λ       0.0  0.0  0.0
            λ       λ       λ + 2μ  0.0  0.0  0.0
            0.0     0.0     0.0     μ    0.0  0.0
            0.0     0.0     0.0     0.0  μ    0.0
            0.0     0.0     0.0     0.0  0.0  μ
        ]

        println("Nodes:                ", nnodes)
        println("Elements:             ", nelements)

        # Warm-up: compilation is intentionally excluded from the measurements.
        K = ∫(
            SymGrad(problem) ⋅ D ⋅ SymGrad(problem), assembly=:csc
        )

        f = ∫(
            problem ⋅ [1.0, 2.0, 3.0]
        )

        K, matrix_time, matrix_bytes = median_timed(REPEATS) do
            ∫(
                SymGrad(problem) ⋅ D ⋅ SymGrad(problem), assembly=:csc
            )
        end

        f, vector_time, vector_bytes = median_timed(REPEATS) do
            ∫(
                problem ⋅ [1.0, 2.0, 3.0]
            )
        end

        println("Matrix assembly:      ",
            round(matrix_time; digits=6), " s, ",
            round(matrix_bytes / 2.0^20; digits=1), " MiB")

        println("Vector assembly:      ",
            round(vector_time; digits=6), " s, ",
            round(vector_bytes / 2.0^20; digits=1), " MiB")

        println("Matrix size:          ", size(K.A))
        println("Matrix nonzeros:      ", nnz(K.A))
        println("Matrix checksum:      ", sum(K.A))
        println("Vector checksum:      ", sum(f.a))

    finally
        if LowLevelFEM.gmsh.isInitialized() != 0
            LowLevelFEM.gmsh.finalize()
        end
    end

    println()
end

if WORKER_FLAG in ARGS
    worker()
else
    run_workers()
end
