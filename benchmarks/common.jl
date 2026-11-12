using LowLevelFEM
using Printf

function run_benchmark(f, name)
    gmsh.initialize(String[], false) # Ignore local Gmsh configuration files.
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        println(name)
        println("Julia $VERSION | LowLevelFEM $(pkgversion(LowLevelFEM)) | Gmsh ",
                gmsh.option.getString("General.Version"))
        f()
    finally
        gmsh.finalize()
    end
    println("PASS ", name)
end

function check_limit(label, value, limit)
    @printf("  %s = %.8g (limit %.8g)\n", label, value, limit)
    isfinite(value) && value <= limit ||
        error("FAIL: $label = $value exceeds limit $limit or is non-finite")
end

function unit_square()
    geo = gmsh.model.geo
    corners = [geo.addPoint(x, y, 0.0) for (x, y) in ((0, 0), (1, 0), (1, 1), (0, 1))]
    edges = [geo.addLine(corners[i], corners[mod1(i + 1, 4)]) for i in 1:4]
    surface = geo.addPlaneSurface([geo.addCurveLoop(edges)])
    geo.synchronize()
    gmsh.model.addPhysicalGroup(2, [surface], -1, "body")
    return surface, edges
end
