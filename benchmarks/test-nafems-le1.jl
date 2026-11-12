include("common.jl")

function le1_stress(n_arc, n_radial)
    gmsh.model.add("nafems-le1-$n_arc")
    geo = gmsh.model.geo
    outer_a, outer_b = 3.25, 2.75
    center = geo.addPoint(0.0, 0.0, 0.0)
    a = geo.addPoint(0.0, 1.0, 0.0)
    b = geo.addPoint(0.0, outer_b, 0.0)
    c = geo.addPoint(outer_a, 0.0, 0.0)
    d = geo.addPoint(2.0, 0.0, 0.0)
    outer = geo.addEllipseArc(c, center, c, b)
    inner = geo.addEllipseArc(d, center, d, a)
    x_axis = geo.addLine(d, c)
    y_axis = geo.addLine(b, a)
    surface = geo.addPlaneSurface([geo.addCurveLoop([x_axis, outer, y_axis, -inner])])
    geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [surface], -1, "body")
    gmsh.model.addPhysicalGroup(1, [x_axis], -1, "x-axis")
    gmsh.model.addPhysicalGroup(1, [y_axis], -1, "y-axis")
    gmsh.model.addPhysicalGroup(1, [outer], -1, "outer-ellipse")
    for edge in (inner, outer)
        gmsh.model.mesh.setTransfiniteCurve(edge, n_arc + 1)
    end
    for edge in (x_axis, y_axis)
        gmsh.model.mesh.setTransfiniteCurve(edge, n_radial + 1)
    end
    gmsh.model.mesh.setTransfiniteSurface(surface, "Left", [d, c, b, a])
    gmsh.model.mesh.setRecombine(2, surface)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)

    problem = Problem([Material("body"; E=210e9, ν=0.3)]; type=:PlaneStress, thickness=0.1)
    support = [displacementConstraint("x-axis"; uy=0.0),
               displacementConstraint("y-axis"; ux=0.0)]
    # Apply 10 MPa outward using the normalized gradient of the outer ellipse.
    pressure = 10e6
    normal_length(x, y) = hypot(x / outer_a^2, y / outer_b^2)
    fx(x, y, z) = pressure * (x / outer_a^2) / normal_length(x, y)
    fy(x, y, z) = pressure * (y / outer_b^2) / normal_length(x, y)
    u = solveDisplacement(problem; support, load=[load("outer-ellipse"; fx, fy)])
    σ = solveStress(u; DoFResults=true)
    σyy = probe(σ, 2.0, 0.0, 0.0)[2, 2] / 1e6
    @printf("  %d × %d Q2 elements: σyy(D) = %.8f MPa\n", n_arc, n_radial, σyy)
    return σyy
end

function nafems_le1()
    reference = 92.7 # MPa, NAFEMS R0006, LE1 (elliptic membrane).
    coarse = le1_stress(12, 8)
    fine = le1_stress(24, 16)
    println("  NAFEMS reference σyy(D) = $reference MPa")
    check_limit("coarse-grid relative stress error [%]", 100 * abs(coarse / reference - 1), 5.0)
    check_limit("fine-grid relative stress error [%]", 100 * abs(fine / reference - 1), 2.0)
end

run_benchmark(nafems_le1, "NAFEMS LE1 elliptic membrane")
