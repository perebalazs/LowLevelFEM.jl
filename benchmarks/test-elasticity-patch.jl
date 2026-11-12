include("common.jl")

function elasticity_patch()
    gmsh.model.add("elasticity-patch")
    _, edges = unit_square()
    for (edge, name) in zip(edges, ("bottom", "right", "top", "left"))
        gmsh.model.addPhysicalGroup(1, [edge], -1, name)
    end
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.2)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)
    gmsh.model.mesh.generate(2)

    # Affine displacement with distinct normal strains and nonzero shear.
    ux(x, y, z) = 1e-3 * (1 + 2x + 3y)
    uy(x, y, z) = 1e-3 * (-2 + x - y)
    E, ν = 210e9, 0.3
    εxx, εyy, γxy = 2e-3, -1e-3, 4e-3
    σxx = E / (1 - ν^2) * (εxx + ν * εyy)
    σyy = E / (1 - ν^2) * (εyy + ν * εxx)
    τxy = E / (2 * (1 + ν)) * γxy

    problem = Problem([Material("body"; E, ν)]; type=:PlaneStress, thickness=0.2)
    support = [displacementConstraint("left"; ux, uy)]
    # Exact tractions σn on the other three sides.
    loads = [load("bottom"; fx=-τxy, fy=-σyy),
             load("right"; fx=σxx, fy=τxy),
             load("top"; fx=τxy, fy=σyy)]
    u = solveDisplacement(problem; support, load=loads)
    σ = solveStress(u; DoFResults=true)

    tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = eachcol(reshape(coordinates, 3, :))
    u_exact = [f(point...) for f in (ux, uy), point in points]
    u_fem = reshape(DoFs(u), 2, :)[:, tags]
    σ_exact = vec([σxx τxy 0.0; τxy σyy 0.0; 0.0 0.0 0.0])
    σ_fem = reshape(DoFs(σ), 9, :)[:, tags]
    println("  Linear triangles, $(length(tags)) nodes")
    check_limit("normalized max displacement error",
                maximum(abs, u_fem .- u_exact) / maximum(abs, u_exact), 1e-10)
    check_limit("normalized max nodal stress error",
                maximum(abs, σ_fem .- σ_exact) / maximum(abs, σ_exact), 1e-9)
end

run_benchmark(elasticity_patch, "Analytic plane-stress patch test")
