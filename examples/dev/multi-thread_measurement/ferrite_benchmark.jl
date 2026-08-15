# Minimal standalone Ferrite assembly benchmark (see ../README.md for the full setup).
# Usage: julia --threads=N ferrite_benchmark.jl
# Requires: Ferrite, FerriteGmsh, OhMyThreads

using Ferrite, FerriteGmsh, OhMyThreads, SparseArrays, LinearAlgebra

const h = 0.025
const repeats = 5

function make_grid(h)
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
    grid = FerriteGmsh.togrid()
    gmsh.finalize()
    return grid
end

function material_stiffness()
    E = 2.0e5; ν = 0.3
    λ = E * ν / ((1 + ν) * (1 - 2ν)); μ = E / (2(1 + ν))
    δ(i, j) = i == j ? 1.0 : 0.0
    return SymmetricTensor{4, 3}() do i, j, k, l
        λ * δ(i, j) * δ(k, l) + μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
    end
end

function assemble_cell!(Ke, fe, cv, C, b)
    for q in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q)
        for i in 1:getnbasefunctions(cv)
            fe === nothing || (fe[i] += (shape_value(cv, q, i) ⋅ b) * dΩ)
            if Ke !== nothing
                ∇δui = shape_symmetric_gradient(cv, q, i)
                for j in 1:getnbasefunctions(cv)
                    Ke[i, j] += (∇δui ⊡ C ⊡ shape_symmetric_gradient(cv, q, j)) * dΩ
                end
            end
        end
    end
end

# Grid coloring approach from the Ferrite howto on multithreaded assembly:
# parallelize the cell loop within each color; task local scratch via @local.
function assemble_global!(K, f, dh, colors, cvt, C, b; ntasks = Threads.nthreads())
    K === nothing || (nonzeros(K) .= 0)
    f === nothing || fill!(f, 0)
    n = ndofs_per_cell(dh)
    for color in colors
        scheduler = OhMyThreads.DynamicScheduler(; ntasks)
        OhMyThreads.@tasks for cellidx in color
            @set scheduler = scheduler
            @local scratch = (;
                cc = CellCache(dh), cv = copy(cvt),
                Ke = K === nothing ? nothing : zeros(n, n),
                fe = f === nothing ? nothing : zeros(n),
                asm = K === nothing ? nothing :
                    start_assemble(K, f === nothing ? Float64[] : f; fillzero = false),
            )
            (; cc, cv, Ke, fe, asm) = scratch
            reinit!(cc, cellidx); reinit!(cv, cc)
            Ke === nothing || fill!(Ke, 0)
            fe === nothing || fill!(fe, 0)
            assemble_cell!(Ke, fe, cv, C, b)
            if asm !== nothing
                fe === nothing ? assemble!(asm, celldofs(cc), Ke) :
                    assemble!(asm, celldofs(cc), Ke, fe)
            else
                assemble!(f, celldofs(cc), fe)  # safe: colors don't share dofs
            end
        end
    end
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
    grid = make_grid(h)
    ip = Lagrange{RefTetrahedron, 1}()^3
    cv = CellValues(QuadratureRule{RefTetrahedron}(3), ip)
    dh = close!(add!(DofHandler(grid), :u, ip))
    coloring = @timed create_coloring(grid)
    pattern = @timed allocate_matrix(dh)
    colors, K = coloring.value, pattern.value
    f = zeros(ndofs(dh))
    C, b = material_stiffness(), Vec{3}((1.0, 2.0, 3.0))

    println("Threads: ", Threads.nthreads(), ", elements: ", getncells(grid),
        ", dofs: ", ndofs(dh))
    println("Coloring:           ", round(coloring.time; digits = 3), " s")
    println("Pattern allocation: ", round(pattern.time; digits = 3), " s")
    for (label, fn) in [
            "Matrix assembly" => () -> assemble_global!(K, nothing, dh, colors, cv, C, b),
            "Vector assembly" => () -> assemble_global!(nothing, f, dh, colors, cv, C, b),
            "Matrix + vector" => () -> assemble_global!(K, f, dh, colors, cv, C, b),
        ]
        fn() # warm-up (exclude compilation)
        s = median_timed(fn, repeats)
        println(label, ":    ", round(s.time; digits = 4), " s, ",
            round(s.bytes / 2^20; digits = 1), " MiB")
    end
    println("Checksums: ", sum(abs, K), ", ", sum(f))
end

main()
