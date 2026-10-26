# Run each scenario in a fresh process from the package root, for example:
# julia --startup-file=no --project=. --threads=2 test/precompile_latency.jl dsl
# Scenarios: dsl, mixed, eigen, contact, fields, legacy.
# MEASURE columns: scenario, operation, first wall time, first compile time,
# first recompile time, second wall time, second compile time (seconds).
# Per-stage timings are ordered: later stages share earlier stages' compiled code.
# Compilation counters may include parallel worker activity and exceed wall time.
# This developer diagnostic is deliberately separate from runtests.jl.
isempty(ARGS) && push!(ARGS, "dsl")
load_measurement = @timed @eval using LowLevelFEM, LinearAlgebra, SparseArrays
println("LOAD\t", ARGS[1], "\t", load_measurement.time, "\t", load_measurement.compile_time, "\t", load_measurement.recompile_time)

function sample(f, name)
    first = @timed Base.invokelatest(f)
    second = @timed Base.invokelatest(f)
    println("MEASURE\t", ARGS[1], "\t", name, "\t", first.time, "\t",
            first.compile_time, "\t", first.recompile_time, "\t", second.time,
            "\t", second.compile_time)
    flush(stdout)
    return first.value
end

function contact_fixture()
    gmsh.clear()
    gmsh.model.add("latency_contact")
    lower = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0)
    upper = gmsh.model.occ.addRectangle(0.0, 1.1, 0.0, 1.0, 1.0)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [lower, upper], -1, "body")
    for (name, entity, y) in (("master", lower, 1.0), ("slave", upper, 1.1))
        sides = gmsh.model.getBoundary([(2, entity)], false, false)
        tags = [tag for (d, tag) in sides if abs(gmsh.model.occ.getCenterOfMass(d, tag)[2] - y) < 1e-10]
        gmsh.model.addPhysicalGroup(1, tags, -1, name)
    end
    for (_, tag) in gmsh.model.getEntities(1)
        gmsh.model.mesh.setTransfiniteCurve(tag, 2)
    end
    for (_, tag) in gmsh.model.getEntities(2)
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)
    end
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.model.mesh.generate(2)
    return nothing
end

@assert gmsh.isInitialized() == 0 "Package loading must not initialize Gmsh."
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
try
    if ARGS[1] == "dsl"
        sample("mesh") do; structured_rect_mesh(n=1, order=1); end
        mat = Material("body")
        P = Field([mat], type=:ScalarField, dim=2, fieldName=:p)
        U = Field([mat], type=:VectorField, dim=2, fieldName=:u)
        C = D(:PlaneStress, mat)
        K = sample("elasticity") do; ∫(SymGrad(U) ⋅ C ⋅ SymGrad(U)); end
        f = sample("surface_load") do; ∫(U ⋅ [1.0, 0.0]; Γ="right"); end
        sample("scalar_diffusion") do; ∫(Grad(P) ⋅ Grad(P)); end
        sample("mass") do; ∫(P ⋅ P); end
        c = scalarField(P, "body", 2.0)
        sample("field_coefficient") do; ∫(Grad(P) ⋅ c ⋅ Grad(P)); end
        B = Grad(P) + [1.0; 0.0;;] ⋅ P
        sample("compound") do; ∫(B' ⋅ [2.0 0.0; 0.0 1.0] ⋅ B); end
        sample("ijv") do; ∫(SymGrad(U) ⋅ C ⋅ SymGrad(U); assembly=:ijv); end
        support = [BoundaryCondition("left"; field=U, ux=0.0, uy=0.0)]
        sample("solve") do; solveField(K, f; support=support); end
        sample("symmetric_solve") do; solveField(Symmetric(K), f; support=support); end
        mpc = MPC(master="rightbottom", slave="right", field=U, ux=true, uy=true)
        sample("mpc") do; solveField(K, f; support=support, mpc=[mpc]); end
        cs = CoordinateSystem(VectorField(U, "body", [0.0, 1.0, 0.0]))
        sample("coordinates") do; solveField(K, Global(f); support=support, coordSys=[cs]); end
    elseif ARGS[1] == "mixed"
        structured_rect_mesh(n=2, order=2)
        mat = Material("body")
        P = Field([mat], type=:ScalarField, dim=2, fieldName=:p, reducedOrder=true)
        V = Field([mat], type=:VectorField, dim=2, fieldName=:v)
        A = ∫((SymGrad(V) ⋅ SymGrad(V)) * 2.0)
        B = sample("coupling") do; ∫(Div(V) ⋅ P); end
        C = 0.0 * ∫(P ⋅ P)
        K = sample("block_matrix") do; SystemMatrix([A B; B' C]); end
        F = SystemVector([∫(V ⋅ [1.0, 0.0]), ∫(P ⋅ 0.0)])
        support = [BoundaryCondition("rightbottom"; field=P, p=0.0),
                   BoundaryCondition("top"; field=V, vx=0.0, vy=0.0),
                   BoundaryCondition("bottom"; field=V, vx=0.0, vy=0.0)]
        sample("block_solve_reduced") do; solveField(K, F; support=support); end
    elseif ARGS[1] == "eigen"
        line_mesh(n=2, order=1)
        mat = Material("body")
        P = Field([mat], type=:ScalarField, dim=1, fieldName=:p)
        K, M = ∫(Grad(P) ⋅ Grad(P)), ∫(P ⋅ P)
        bc = [BoundaryCondition("left"; field=P, p=0.0), BoundaryCondition("right"; field=P, p=0.0)]
        sample("direct_modes") do; solveEigen(K, M; n=1, solver=:direct, support=bc); end
        sample("mass_normalized") do; solveEigen(Symmetric(K), Symmetric(M); n=1, solver=:direct, support=bc, normalize=:mass); end
        sample("raw_modes") do; solveEigen(K, M; n=1, solver=:direct, support=bc, output=:raw); end
    elseif ARGS[1] == "contact"
        contact_fixture()
        U = Field([Material("body")], type=:VectorField, dim=2, fieldName=:u)
        u = vectorField(U, "body", [0.0, 0.0, 0.0])
        C = sample("geometry") do; contact(U; master="master", slave="slave", displacement=u); end
        G = ContactGap(C; active=:all)
        sample("penalty_assembly") do; ∫(G ⋅ 2.0 ⋅ G); end
        sample("gap_postprocess") do; ContactGap(C, u; gauss=:full); end
        sample("update") do; updateContact!(C, u); end
        P = Field([Material("body")], type=:ScalarField, dim=2, fieldName=:λ)
        sample("multiplier_coupling") do; ∫(P ⋅ G); end
    elseif ARGS[1] == "fields"
        structured_box_mesh(n=1, order=1)
        P = Field([Material("body")], type=:VectorField, dim=3, fieldName=:u)
        s = scalarField(P, "body", 2.0)
        S = ScalarField(P, "body", 2.0)
        V = VectorField(P, "body", [1.0, 2.0, 3.0])
        T = TensorField(P, "body", [2.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 4.0])
        sample("nodal_gradient") do; ∇(s); end
        sample("element_gradient") do; ∇(S); end
        sample("vector_gradient") do; ∇(V); end
        sample("curl") do; ∇(V; nabla=:curl); end
        sample("tensor_divergence") do; ∇(T; nabla=:div); end
        sample("restriction") do; elementsToElements(T; onPhysicalGroup="left"); end
        sample("scalar_algebra") do; sqrt(S * S + 1.0); end
        sample("tensor_algebra") do; (inv(T), eigen(T), sqrt(T), log(T)); end
    elseif ARGS[1] == "legacy"
        structured_rect_mesh(n=1, order=1)
        P = Problem([Material("body")], type=:PlaneStress)
        sample("stiffness") do; stiffnessMatrix(P); end
        support = [displacementConstraint("left"; ux=0.0), displacementConstraint("bottom"; uy=0.0)]
        u = sample("displacement") do; solveDisplacement(P; load=[load("right"; fx=1.0)], support=support); end
        sample("stress") do; solveStress(u); end
    else
        error("Unknown scenario")
    end
finally
    gmsh.clear()
    gmsh.finalize()
end
