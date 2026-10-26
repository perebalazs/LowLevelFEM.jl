# Helpers own Gmsh cleanup, even when an assertion or assembly throws.
# Mesh sizes are explicit: the public mesh helpers default to much larger models.
function with_test_gmsh(f)
    initialized_here = gmsh.isInitialized() == 0
    initialized_here && gmsh.initialize()
    terminal = gmsh.option.getNumber("General.Terminal")
    verbosity = gmsh.option.getNumber("General.Verbosity")
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.clear()
        return f()
    finally
        if gmsh.isInitialized() != 0
            gmsh.clear()
            gmsh.option.setNumber("General.Terminal", terminal)
            gmsh.option.setNumber("General.Verbosity", verbosity)
            initialized_here && gmsh.finalize()
        end
    end
end

relerr(a, b) = norm(a - b) / max(norm(b), eps(Float64))

function test_constant_field(f, expected; atol=1e-10)
    values = expected isa Number ? [expected] : vec(expected)
    blocks = isNodal(f) ? (f.a,) : f.A
    @test !isempty(blocks)
    for block in blocks
        samples = reshape(block, length(values), :)
        @test samples ≈ repeat(values, 1, size(samples, 2)) atol=atol rtol=1e-10
    end
end

function test_eigen_residual(A, B, values, vectors; transform=nothing, atol=1e-8)
    @test length(values) == size(vectors, 2)
    @test issorted(real.(values))
    for j in eachindex(values)
        v = vectors[:, j]
        Av, Bv = A * v, B * v
        residual = Av - values[j] * Bv
        # Essential/MPC constraints admit reactions. Check equilibrium only
        # against admissible variations, rather than at constrained full DOFs.
        residual = transform === nothing ? residual : transform' * residual
        scale = max(norm(Av) + abs(values[j]) * norm(Bv), eps(Float64))
        @test norm(residual) / scale < atol
    end
end

# Taylor-Hood-style mixed regression fixture, shared by block and BC tests.
function stokes_fixture()
    structured_rect_mesh(
        lx=2.0,
        n=2,
        order=2
    )

    material = Material("body")

    P = Field(
        [material],
        type=:ScalarField,
        fieldName=:p,
        rhsName=:fp,
        dim=2,
        reducedOrder=true
    )

    V = Field(
        [material],
        type=:VectorField,
        fieldName=:v,
        rhsName=:fv,
        dim=2
    )

    pres = BoundaryCondition(
        "rightbottom",
        field=P,
        p=0.0
    )

    supp_top = BoundaryCondition(
        "top",
        field=V,
        vx=0.0,
        vy=0.0
    )

    supp_bottom = BoundaryCondition(
        "bottom",
        field=V,
        vx=0.0,
        vy=0.0
    )

    support = [
        pres,
        supp_top,
        supp_bottom
    ]

    fv_x = ∫(V ⋅ [1.0, 0.0])

    fp = ∫(P ⋅ 0.0)
    F_x = SystemVector([fv_x, fp])

    μ = 1.0
    γ = 1e-1

    A = ∫((SymGrad(V) ⋅ SymGrad(V)) * 2μ)
    B = ∫(Div(V) ⋅ P)
    C = 0 * ∫(P ⋅ P)
    D = ∫(Div(V) ⋅ Div(V) * γ)

    K = SystemMatrix([
        A + D  B
        B'     C
    ])


    return (; P, V, A, B, C, D, K, fv_x, fp, F_x,
            pres, supp_top, supp_bottom, support)
end

# Two disconnected unit blocks with flat matching interfaces. The lower
# block's top is the master and the upper block's bottom is the slave.
function flat_contact_fixture(; dim=2, gap=0.1)
    gmsh.clear()
    gmsh.model.add("flat_contact_$(dim)d")
    if dim == 2
        lower = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0)
        upper = gmsh.model.occ.addRectangle(0.0, 1.0 + gap, 0.0, 1.0, 1.0)
    else
        lower = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        upper = gmsh.model.occ.addBox(0.0, 1.0 + gap, 0.0, 1.0, 1.0, 1.0)
    end
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(dim, [lower, upper], -1, "body")
    gmsh.model.addPhysicalGroup(dim, [lower], -1, "lower")
    gmsh.model.addPhysicalGroup(dim, [upper], -1, "upper")
    for (name, entity, y) in (("master", lower, 1.0),
                              ("slave", upper, 1.0 + gap),
                              ("bottom", lower, 0.0))
        sides = gmsh.model.getBoundary([(dim, entity)], false, false)
        tags = [tag for (d, tag) in sides
                if abs(gmsh.model.occ.getCenterOfMass(d, tag)[2] - y) < 1e-10]
        @assert length(tags) == 1
        gmsh.model.addPhysicalGroup(dim - 1, tags, -1, name)
    end
    for (_, curve) in gmsh.model.getEntities(1)
        gmsh.model.mesh.setTransfiniteCurve(curve, 2)
    end
    for (_, surface) in gmsh.model.getEntities(2)
        gmsh.model.mesh.setTransfiniteSurface(surface)
        gmsh.model.mesh.setRecombine(2, surface)
    end
    if dim == 3
        for volume in (lower, upper)
            gmsh.model.mesh.setTransfiniteVolume(volume)
        end
    end
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.model.mesh.generate(dim)
    U = Field([Material("body")], type=:VectorField, dim=dim, fieldName=:u)
    zero_u = vectorField(U, "body", [0.0, 0.0, 0.0])
    # The master manifold's orientation is mesher-dependent; choose the
    # documented open-positive sign explicitly for this geometric fixture.
    C = contact(U; master="master", slave="slave", displacement=zero_u)
    if sum(C.gap_values) < 0
        C = contact(U; master="master", slave="slave", displacement=zero_u,
                    normal_sign=-1.0)
    end
    return (; U, C, zero_u, gap)
end
