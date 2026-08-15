# Linear 3D elasticity benchmark for Ferrite.jl.
#
# Run with, for example:
#   julia --threads=4 ferrite_elasticity_benchmark.jl
#
# The script solves a unit cube fixed at x = 0 and loaded by a unit traction
# in the x direction at x = 1. Timings exclude Julia compilation.

using Ferrite
using LinearAlgebra
using OhMyThreads
using SparseArrays
using Statistics
using FerriteGmsh

const n = 30
const repeats = 5
const E = 2.0e5
const nu = 0.3
const lambda = E * nu / ((1 + nu) * (1 - 2 * nu))
const mu = E / (2 * (1 + nu))

function median_timed(f; repeats = repeats)
    f() # Warm-up: compilation is deliberately excluded.
    samples = map(1:repeats) do _
        GC.gc()
        @timed f()
    end
    return sort(samples; by = s -> s.time)[cld(repeats, 2)]
end

function elasticity_tensor()
    delta(i, j) = i == j ? 1.0 : 0.0
    return SymmetricTensor{4,3}() do i, j, k, l
        lambda * delta(i, j) * delta(k, l) +
        mu * (delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k))
    end
end

function assemble_volume!(K, f, dh, colors, cv_template, C)
    nonzeros(K) .= 0.0
    fill!(f, 0.0)
    ndofs_cell = ndofs_per_cell(dh)

    # Cells of one color do not share degrees of freedom.
    for color in colors
        scheduler = DynamicScheduler(; ntasks = Threads.nthreads())
        @tasks for cellidx in color
            @set scheduler = scheduler
            @local scratch = (
                cc = CellCache(dh),
                cv = copy(cv_template),
                Ke = zeros(ndofs_cell, ndofs_cell),
                fe = zeros(ndofs_cell),
                assembler = start_assemble(K, f; fillzero = false),
            )
            (; cc, cv, Ke, fe, assembler) = scratch
            reinit!(cc, cellidx)
            reinit!(cv, cc)
            fill!(Ke, 0.0)
            fill!(fe, 0.0)

            for q in 1:getnquadpoints(cv)
                dOmega = getdetJdV(cv, q)
                for i in 1:getnbasefunctions(cv)
                    eps_i = shape_symmetric_gradient(cv, q, i)
                    for j in 1:getnbasefunctions(cv)
                        eps_j = shape_symmetric_gradient(cv, q, j)
                        Ke[i, j] += (eps_i ⊡ C ⊡ eps_j) * dOmega
                    end
                end
            end
            assemble!(assembler, celldofs(cc), Ke, fe)
        end
    end
end

function assemble_traction!(f, dh, right_facets, fv_template)
    traction = Vec{3}((1.0, 0.0, 0.0))
    for facet in FacetIterator(dh, right_facets)
        fv = copy(fv_template)
        reinit!(fv, facet)
        fe = zeros(ndofs_per_cell(dh))
        for q in 1:getnquadpoints(fv), i in 1:getnbasefunctions(fv)
            fe[i] += (shape_value(fv, q, i) ⋅ traction) * getdetJdV(fv, q)
        end
        assemble!(f, celldofs(facet), fe)
    end
end

function make_grid(n)
    gmsh.initialize()
    gmsh.open("elasticity_benchmark.msh")
    grid = FerriteGmsh.togrid()
    gmsh.finalize()
    #grid = generate_grid(Hexahedron, (n, n, n), Vec((0.0, 0.0, 0.0)), Vec((1.0, 1.0, 1.0)))
    ## generate_grid provides the standard "left" and "right" facetsets.
    return grid
end

function setup_problem(grid)

    ip = Lagrange{RefHexahedron,1}()^3
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    ch = ConstraintHandler(dh)
    zero_displacement(x, t) = Vec{3}((0.0, 0.0, 0.0))
    add!(ch, Dirichlet(:u, getfacetset(grid, "left"), zero_displacement, [1, 2, 3]))
    close!(ch)
    update!(ch, 0.0)

    cv = CellValues(QuadratureRule{RefHexahedron}(2), ip)
    fv = FacetValues(FacetQuadratureRule{RefHexahedron}(2), ip)
    colors = create_coloring(grid)
    K = allocate_matrix(dh)
    f = zeros(ndofs(dh))
    return (; grid, dh, ch, cv, fv, colors, K, f, C = elasticity_tensor(),
        right_facets = getfacetset(grid, "right"))
end

function assemble_and_solve!(p)
    assemble_volume!(p.K, p.f, p.dh, p.colors, p.cv, p.C)
    assemble_traction!(p.f, p.dh, p.right_facets, p.fv)
    apply!(p.K, p.f, p.ch)
    return cholesky(Symmetric(p.K)) \ p.f
end

function stress_checksum(p, u)
    total = 0.0
    for cell in CellIterator(p.dh)
        cv = copy(p.cv)
        reinit!(cv, cell)
        ue = u[celldofs(cell)]
        for q in 1:getnquadpoints(cv)
            eps = function_symmetric_gradient(cv, q, ue)
            sigma = lambda * tr(eps) * one(eps) + 2 * mu * eps
            total += (sigma ⊡ sigma) * getdetJdV(cv, q)
        end
    end
    return total
end

function complete_analysis(grid)
    # This deliberately rebuilds the dof handler, coloring, and sparse pattern.
    p = setup_problem(grid)
    u = assemble_and_solve!(p)
    return (; p, u, stress = stress_checksum(p, u))
end

function main()
    BLAS.set_num_threads(1)
    mesh = @timed make_grid(n)
    analysis_stats = median_timed(() -> complete_analysis(mesh.value))
    result = analysis_stats.value

    println("Ferrite.jl linear elasticity benchmark")
    println("Threads: ", Threads.nthreads(), ", cells: ", getncells(mesh.value), ", dofs: ", ndofs(result.p.dh))
    println("Mesh generation (not timed):    ", round(mesh.time; digits = 4), " s, ", round(mesh.bytes / 2^20; digits = 1), " MiB")
    println("Complete analysis (median):     ", round(analysis_stats.time; digits = 4), " s, ", round(analysis_stats.bytes / 2^20; digits = 1), " MiB")
    println("Includes: dofs, coloring, pattern, assembly, solve, and stress postprocessing")
    println("Checksums: ", sum(abs, result.p.K), ", ", sum(result.u), ", ", result.stress)
end

main()
