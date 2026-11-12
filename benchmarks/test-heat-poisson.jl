include("common.jl")

function heat_error(n)
    gmsh.model.add("heat-poisson-$n")
    surface, edges = unit_square()
    gmsh.model.addPhysicalGroup(1, edges, -1, "boundary")
    for edge in edges
        gmsh.model.mesh.setTransfiniteCurve(edge, n + 1)
    end
    gmsh.model.mesh.setTransfiniteSurface(surface)
    gmsh.model.mesh.setRecombine(2, surface)
    gmsh.model.mesh.generate(2)

    # -k ΔT = q on the unit square, with zero temperature on the boundary.
    k = 3.0
    exact(x, y, z) = sinpi(x) * sinpi(y)
    source(x, y, z) = 2π^2 * k * exact(x, y, z)
    problem = Problem([Material("body"; k)]; type=:PlaneHeatConduction)
    T = solveTemperature(problem;
                         heatFlux=[heatFlux("body"; qn=source)],
                         temperatureConstraint=[temperatureConstraint("boundary"; T=0.0)])

    tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = eachcol(reshape(coordinates, 3, :))
    err = maximum(abs, DoFs(T)[tags] .- [exact(point...) for point in points])
    @printf("  %d × %d Q1 elements: max nodal error %.8g K\n", n, n, err)
    return err
end

function heat_poisson()
    coarse = heat_error(8)
    fine = heat_error(16)
    isfinite(coarse) || error("FAIL: coarse-grid nodal error is non-finite")
    check_limit("fine-grid max nodal error [K]", fine, 0.004)
    check_limit("fine/coarse nodal error ratio", fine / coarse, 1 / 3)
end

run_benchmark(heat_poisson, "Analytic steady heat conduction")
