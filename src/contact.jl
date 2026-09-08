###############################################################################
#                                                                             #
#                              Contact                                        #
#                                                                             #
###############################################################################

export Contact, contact, updateContact!, contactPenaltyMatrix
export CONTACT_OPEN, CONTACT_STICK, CONTACT_SLIP

const CONTACT_OPEN  = UInt8(0)
const CONTACT_STICK = UInt8(1)
const CONTACT_SLIP  = UInt8(2)

# -----------------------------------------------------------------------------
# Internal geometric data
# -----------------------------------------------------------------------------

struct _ContactElement
    tag::Int
    etype::Int
    dim::Int
    name::String
    order::Int
    node_tags::Vector{Int}
    coords::Matrix{Float64}          # 3 × number of element nodes, current configuration
    local_nodes::Matrix{Float64}     # dim × number of element nodes
    bmin::Vector{Float64}
    bmax::Vector{Float64}
end

mutable struct _ContactAABBNode
    bmin::Vector{Float64}
    bmax::Vector{Float64}
    left::Union{Nothing,_ContactAABBNode}
    right::Union{Nothing,_ContactAABBNode}
    elements::Vector{Int}
end

struct _ContactProjection
    element_index::Int
    element_tag::Int
    ξ::Vector{Float64}
    x::Vector{Float64}
    N::Vector{Float64}
    normal::Vector{Float64}
    tangent1::Vector{Float64}
    tangent2::Union{Nothing,Vector{Float64}}
    gap::Float64
    distance2::Float64
end

# -----------------------------------------------------------------------------
# Public contact object
# -----------------------------------------------------------------------------

"""
    Contact

Container for one slave-master contact pair.

The contact kinematics are stored in a reduced local contact space: every slave
contact node contributes `pdim` rows to `G`. In 2D the local ordering is
`(normal, tangent)`, while in 3D it is `(normal, tangent1, tangent2)`.

Therefore, for `nc` slave contact nodes,

    size(G) == (nc * displacement.pdim, ndofs(displacement))

and the penalty contribution is represented directly by

    Kc = G' * C * G

where `C` is block diagonal in the local contact basis. For isotropic friction
its local blocks are

    [cn  0 ]                 # 2D
    [0   ct]

and

    [cn  0   0 ]             # 3D
    [0   ct  0 ]
    [0   0   ct]

on active contact points. Inactive blocks are zero.

If a Lagrange multiplier field is supplied, `multiplier_dofs` maps the reduced
rows of `G` to the corresponding global multiplier DoFs on the slave nodes.
This mapping is intentionally stored in the `Contact` object; the solver can
later eliminate inactive multiplier DoFs without changing the contact
kinematics.
"""
mutable struct Contact
    master::String
    slave::String

    displacement::Problem
    multiplier::Union{Nothing,Problem}
    U::VectorField

    slave_nodes::Vector{Int}
    master_element_tags::Vector{Int}
    master_local_coordinates::Vector{Vector{Float64}}
    master_points::Matrix{Float64}

    gap::ScalarField
    G::SparseMatrixCSC{Float64,Int}
    C::SparseMatrixCSC{Float64,Int}

    n::VectorField
    t1::VectorField
    t2::Union{Nothing,VectorField}

    active::BitVector
    state::Vector{UInt8}

    # Local contact-space vectors. They are initialized here and are intended
    # to hold the converged traction/slip history required by Coulomb friction.
    traction::Vector{Float64}
    slip::Vector{Float64}

    cn::Any
    ct::Any
    μ::Any
    cn_values::Vector{Float64}
    ct_values::Vector{Float64}
    μ_values::Vector{Float64}

    # Mapping from reduced contact rows to global multiplier DoFs.
    multiplier_dofs::Vector{Int}

    # Geometry/search options reused by updateContact!.
    options::NamedTuple
end

function Base.show(io::IO, c::Contact)
    nc = length(c.slave_nodes)
    na = count(c.active)
    nstick = count(==(CONTACT_STICK), c.state)
    nslip = count(==(CONTACT_SLIP), c.state)

    print(
        io,
        "Contact(\"$(c.slave)\" -> \"$(c.master)\", " *
        "$nc candidate nodes, $na active, " *
        "stick=$nstick, slip=$nslip, " *
        "G=$(size(c.G)), C=$(size(c.C)))"
    )
end

# -----------------------------------------------------------------------------
# Public interface
# -----------------------------------------------------------------------------

"""
    contact(U::VectorField; master, slave, displacement=U.model,
            LagrangeMultiplierField=nothing, cn=1.0, ct=0.0, μ=0.0,
            activation_tol=0.0, step=U.nsteps, kwargs...) -> Contact

Construct one node-to-manifold contact pair in the current deformed
configuration.

The mesh coordinates used by the contact search are

    x = X + U

for both the slave and master sides. Thus the closest-point projection, normal,
tangential basis, signed gap and kinematic matrix are all evaluated on the
current geometry represented by `U`.

# Arguments

- `U::VectorField`: current displacement field.
- `master::String`: master physical group.
- `slave::String`: slave physical group. Contact quantities are discretized on
  the slave nodes.
- `displacement::Problem=U.model`: displacement trial field associated with the
  columns of `G`.
- `LagrangeMultiplierField::Union{Nothing,Problem}=nothing`: optional vector
  multiplier field. In 2D it must have two components and in 3D three
  components. Only the slave-node DoFs are mapped into the reduced contact
  space.
- `cn`: normal penalty stiffness. May be a number, `ScalarField`, or function
  `f(x,y,z)`.
- `ct`: isotropic tangential penalty stiffness. May be a number,
  `ScalarField`, or function `f(x,y,z)`. Set `ct=0` for frictionless contact.
- `μ`: isotropic Coulomb friction coefficient. It is stored in the contact
  state for the later stick-slip update. May be a number, `ScalarField`, or
  function `f(x,y,z)`.
- `activation_tol::Real=0.0`: a point is active when `gap <= activation_tol`.
- `step::Int=U.nsteps`: displacement time/load step used to construct the
  current geometry.

# Geometry options

- `normal_sign::Real=1.0`: multiply the master-side contact normal by `+1` or
  `-1`. The signed gap is `(x_slave - x_master) ⋅ n`.
- `self_contact::Bool=(slave == master)`: enable local-topology exclusion.
- `self_exclusion_layers::Int=1`: number of additional node-connected master
  element layers excluded in self-contact.
- `aabb_padding::Real=0.05`: relative AABB padding.
- `leaf_size::Int=8`: maximum number of elements in an AABB leaf.
- `projection_tol::Real=1e-10`: closest-point solver tolerance.
- `projection_maxiter::Int=40`: maximum projected Gauss-Newton iterations.

# Notes

`G` contains all slave candidate nodes, including currently open ones. Contact
activation is represented by `active` and by zero local blocks in `C`. This is
useful because the kinematic row layout remains fixed while the active set
changes.

The present file constructs the contact geometry, kinematics, penalty matrix and
state containers. Solver-side active-set elimination and the nonlinear
stick-slip iteration are intentionally handled separately.
"""
function contact(
    U::VectorField;
    master::String,
    slave::String,
    displacement::Problem=U.model,
    LagrangeMultiplierField::Union{Nothing,Problem}=nothing,
    cn=1.0,
    ct=0.0,
    μ=0.0,
    activation_tol::Real=0.0,
    step::Int=U.nsteps,
    normal_sign::Real=1.0,
    self_contact::Bool=(slave == master),
    self_exclusion_layers::Int=1,
    aabb_padding::Real=0.05,
    leaf_size::Int=8,
    projection_tol::Real=1e-10,
    projection_maxiter::Int=40
    )

    options = (
        activation_tol=Float64(activation_tol),
        normal_sign=Float64(normal_sign),
        self_contact=self_contact,
        self_exclusion_layers=self_exclusion_layers,
        aabb_padding=Float64(aabb_padding),
        leaf_size=leaf_size,
        projection_tol=Float64(projection_tol),
        projection_maxiter=projection_maxiter
    )

    data = _contact_build_data(
        U,
        displacement,
        LagrangeMultiplierField,
        slave,
        master,
        cn,
        ct,
        μ;
        step=step,
        options...
    )

    nrows = size(data.G, 1)
    state = Vector{UInt8}(undef, length(data.active))
    @inbounds for i in eachindex(state)
        state[i] = data.active[i] ? CONTACT_STICK : CONTACT_OPEN
    end

    return Contact(
        master,
        slave,
        displacement,
        LagrangeMultiplierField,
        U,
        data.slave_nodes,
        data.master_element_tags,
        data.master_local_coordinates,
        data.master_points,
        data.gap,
        data.G,
        data.C,
        data.n,
        data.t1,
        data.t2,
        data.active,
        state,
        zeros(Float64, nrows),
        zeros(Float64, nrows),
        cn,
        ct,
        μ,
        data.cn_values,
        data.ct_values,
        data.μ_values,
        data.multiplier_dofs,
        options
    )
end

"""
    contact(displacement::Problem; kwargs...) -> Contact

Construct contact in the undeformed configuration (`U = 0`).
"""
function contact(
    displacement::Problem;
    kwargs...
    )

    U = _contact_zero_displacement(displacement)
    return contact(U; displacement=displacement, kwargs...)
end

# Positional compatibility helpers.
contact(U::VectorField, slave::String, master::String; kwargs...) =
    contact(U; slave=slave, master=master, kwargs...)

contact(displacement::Problem, slave::String, master::String; kwargs...) =
    contact(displacement; slave=slave, master=master, kwargs...)

"""
    updateContact!(c::Contact, U::VectorField; step=U.nsteps) -> Contact

Recompute the current contact geometry, active set, `G` and `C` from a new
displacement field. The slave-node ordering remains tied to the physical group.
Existing traction and slip history are preserved when the slave-node set is
unchanged; traction is cleared on newly open points.
"""
function updateContact!(
    c::Contact,
    U::VectorField;
    step::Int=U.nsteps
    )

    data = _contact_build_data(
        U,
        c.displacement,
        c.multiplier,
        c.slave,
        c.master,
        c.cn,
        c.ct,
        c.μ;
        step=step,
        c.options...
    )

    old_nodes = c.slave_nodes
    old_traction = c.traction
    old_slip = c.slip
    old_state = c.state
    pdim = c.displacement.pdim

    c.U = U
    c.slave_nodes = data.slave_nodes
    c.master_element_tags = data.master_element_tags
    c.master_local_coordinates = data.master_local_coordinates
    c.master_points = data.master_points
    c.gap = data.gap
    c.G = data.G
    c.C = data.C
    c.n = data.n
    c.t1 = data.t1
    c.t2 = data.t2
    c.active = data.active
    c.cn_values = data.cn_values
    c.ct_values = data.ct_values
    c.μ_values = data.μ_values
    c.multiplier_dofs = data.multiplier_dofs

    nrows = size(c.G, 1)
    c.traction = zeros(Float64, nrows)
    c.slip = zeros(Float64, nrows)
    c.state = fill(CONTACT_OPEN, length(c.slave_nodes))

    if old_nodes == c.slave_nodes && length(old_traction) == nrows
        c.traction .= old_traction
        c.slip .= old_slip

        @inbounds for i in eachindex(c.slave_nodes)
            rows = (i - 1) * pdim + 1:i * pdim
            if c.active[i]
                # Preserve a previous converged stick/slip state. If the point
                # has just closed, initialize it as stick.
                c.state[i] = old_state[i] == CONTACT_OPEN ? CONTACT_STICK : old_state[i]
            else
                c.state[i] = CONTACT_OPEN
                c.traction[rows] .= 0.0
            end
        end
    else
        @inbounds for i in eachindex(c.slave_nodes)
            c.state[i] = c.active[i] ? CONTACT_STICK : CONTACT_OPEN
        end
    end

    return c
end

"""
    contactPenaltyMatrix(c::Contact) -> SystemMatrix

Return the penalty contact contribution

    Kc = G' * C * G

as a `SystemMatrix` associated with the displacement field.
"""
function contactPenaltyMatrix(c::Contact)
    return SystemMatrix(sparse(c.G' * c.C * c.G), c.displacement)
end

# -----------------------------------------------------------------------------
# Main construction
# -----------------------------------------------------------------------------

function _contact_build_data(
    U::VectorField,
    displacement::Problem,
    multiplier::Union{Nothing,Problem},
    slave::String,
    master::String,
    cn,
    ct,
    μ;
    step::Int,
    activation_tol::Float64,
    normal_sign::Float64,
    self_contact::Bool,
    self_exclusion_layers::Int,
    aabb_padding::Float64,
    leaf_size::Int,
    projection_tol::Float64,
    projection_maxiter::Int
    )

    _contact_check_models(U, displacement, multiplier)

    self_exclusion_layers >= 0 ||
        error("contact: self_exclusion_layers must be non-negative.")

    aabb_padding >= 0 ||
        error("contact: aabb_padding must be non-negative.")

    leaf_size >= 1 ||
        error("contact: leaf_size must be at least one.")

    projection_tol > 0 ||
        error("contact: projection_tol must be positive.")

    projection_maxiter >= 1 ||
        error("contact: projection_maxiter must be at least one.")

    isfinite(normal_sign) && abs(abs(normal_sign) - 1.0) <= 10 * eps(Float64) ||
        error("contact: normal_sign must be either +1 or -1.")

    isfinite(activation_tol) ||
        error("contact: activation_tol must be finite.")

    gmsh.model.setCurrent(displacement.name)

    nodecoords = _contact_deformed_coordinates(displacement, U; step=step)

    slave_elements, slave_dim =
        _contact_group_elements(displacement, slave, nodecoords; aabb_padding=0.0)

    master_elements, master_dim =
        _contact_group_elements(displacement, master, nodecoords; aabb_padding=aabb_padding)

    slave_dim in (1, 2) ||
        error("contact: slave physical group '$slave' must be a curve or surface.")

    master_dim in (1, 2) ||
        error("contact: master physical group '$master' must be a curve or surface.")

    isempty(slave_elements) &&
        error("contact: no finite elements were found in slave group '$slave'.")

    isempty(master_elements) &&
        error("contact: no finite elements were found in master group '$master'.")

    _contact_check_master_element_types(master_elements)

    tree = _contact_build_aabb_tree(
        master_elements,
        collect(eachindex(master_elements));
        leaf_size=leaf_size
    )

    slave_nodes = sort!(unique!(vcat((e.node_tags for e in slave_elements)...)))

    if self_contact
        excluded = _contact_self_exclusion_sets(
            slave_nodes,
            master_elements,
            self_exclusion_layers
        )
    else
        excluded = Dict{Int,Set{Int}}()
    end

    projections = Dict{Int,_ContactProjection}()

    for node in slave_nodes
        xs = @view nodecoords[:, node]
        ex = self_contact ? get(excluded, node, Set{Int}()) : nothing

        p = _contact_nearest_projection(
            tree,
            master_elements,
            xs;
            excluded=ex,
            normal_sign=normal_sign,
            model_dim=displacement.dim,
            projection_tol=projection_tol,
            projection_maxiter=projection_maxiter
        )

        p === nothing &&
            error(
                "contact: no admissible master element was found for " *
                "slave node $node. In self-contact, try reducing " *
                "self_exclusion_layers."
            )

        projections[node] = p
    end

    gap = _contact_gap_field(displacement, slave_elements, projections)
    normalVec = _contact_vector_field(displacement, slave_elements, projections, :normal)
    tangent1 = _contact_vector_field(displacement, slave_elements, projections, :tangent1)
    tangent2 = displacement.pdim == 3 ?
        _contact_vector_field(displacement, slave_elements, projections, :tangent2) : nothing

    G = _contact_matrix(
        displacement,
        slave_nodes,
        projections,
        master_elements
    )

    gap_values = [projections[node].gap for node in slave_nodes]
    active = BitVector(g <= activation_tol for g in gap_values)

    cn_values = _contact_parameter_values(cn, slave_nodes, nodecoords, displacement; step=step, name="cn")
    ct_values = _contact_parameter_values(ct, slave_nodes, nodecoords, displacement; step=step, name="ct")
    μ_values  = _contact_parameter_values(μ,  slave_nodes, nodecoords, displacement; step=step, name="μ")

    any(x -> x < 0.0, cn_values) && error("contact: cn must be non-negative.")
    any(x -> x < 0.0, ct_values) && error("contact: ct must be non-negative.")
    any(x -> x < 0.0, μ_values)  && error("contact: μ must be non-negative.")

    C = _contact_stiffness_matrix(
        displacement.pdim,
        active,
        cn_values,
        ct_values
    )

    master_element_tags = Vector{Int}(undef, length(slave_nodes))
    master_local_coordinates = Vector{Vector{Float64}}(undef, length(slave_nodes))
    master_points = Matrix{Float64}(undef, 3, length(slave_nodes))

    @inbounds for (i, node) in enumerate(slave_nodes)
        p = projections[node]
        master_element_tags[i] = p.element_tag
        master_local_coordinates[i] = copy(p.ξ)
        master_points[:, i] .= p.x
    end

    multiplier_dofs = _contact_multiplier_dofs(multiplier, displacement, slave_nodes)

    return (
        slave_nodes=slave_nodes,
        master_element_tags=master_element_tags,
        master_local_coordinates=master_local_coordinates,
        master_points=master_points,
        gap=gap,
        G=G,
        C=C,
        n=normalVec,
        t1=tangent1,
        t2=tangent2,
        active=active,
        cn_values=cn_values,
        ct_values=ct_values,
        μ_values=μ_values,
        multiplier_dofs=multiplier_dofs
    )
end

# -----------------------------------------------------------------------------
# Model and displacement helpers
# -----------------------------------------------------------------------------

function _contact_check_models(
    U::VectorField,
    displacement::Problem,
    multiplier::Union{Nothing,Problem}
    )

    displacement.pdim in (2, 3) ||
        error(
            "contact: displacement must be a 2D or 3D vector field; " *
            "got pdim=$(displacement.pdim)."
        )

    U.model.name == displacement.name ||
        error("contact: U and displacement must use the same Gmsh model.")

    U.model.non == displacement.non ||
        error("contact: U and displacement must use the same mesh nodes.")

    U.model.pdim == displacement.pdim ||
        error(
            "contact: U and displacement must have the same number of " *
            "components per node."
        )

    if multiplier !== nothing
        multiplier.name == displacement.name ||
            error(
                "contact: LagrangeMultiplierField and displacement must use " *
                "the same Gmsh model."
            )

        multiplier.non == displacement.non ||
            error(
                "contact: LagrangeMultiplierField and displacement must use " *
                "the same mesh nodes."
            )

        multiplier.pdim == displacement.pdim ||
            error(
                "contact: the frictional Lagrange multiplier field must have " *
                "$(displacement.pdim) components per node."
            )
    end

    return nothing
end

function _contact_zero_displacement(problem::Problem)
    problem.pdim in (2, 3) ||
        error("contact: zero displacement can only be created for pdim=2 or 3.")

    type = problem.pdim == 2 ? :v2D : :v3D
    return VectorField(
        Matrix{Float64}[],
        zeros(Float64, ndofs(problem), 1),
        [0.0],
        Int[],
        1,
        type,
        problem
    )
end

"""
    _contact_node_coordinates(problem) -> Matrix{Float64}

Return reference mesh-node coordinates in a `3 × problem.non` matrix indexed
by Gmsh node tag.
"""
function _contact_node_coordinates(problem::Problem)
    gmsh.model.setCurrent(problem.name)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()

    length(node_tags) == problem.non ||
        error(
            "contact: Problem node count ($(problem.non)) differs from the " *
            "current Gmsh mesh ($(length(node_tags)))."
        )

    X = zeros(Float64, 3, problem.non)

    @inbounds for (i, tag0) in enumerate(node_tags)
        tag = Int(tag0)
        1 <= tag <= problem.non ||
            error(
                "contact: non-contiguous Gmsh node tags detected. " *
                "Create the Problem after node renumbering."
            )

        X[1, tag] = coords[3i - 2]
        X[2, tag] = coords[3i - 1]
        X[3, tag] = coords[3i]
    end

    return X
end

"""
    _contact_deformed_coordinates(problem, U; step) -> Matrix{Float64}

Return current coordinates `x = X + U` for contact search and projection.
"""
function _contact_deformed_coordinates(
    problem::Problem,
    U::VectorField;
    step::Int
    )

    1 <= step <= U.nsteps ||
        error("contact: displacement step $step is outside 1:$(U.nsteps).")

    Un = isNodal(U) ? U : elementsToNodes(U)

    size(Un.a, 1) == ndofs(problem) ||
        error(
            "contact: nodal displacement contains $(size(Un.a, 1)) values " *
            "per step, expected $(ndofs(problem))."
        )

    size(Un.a, 2) >= step ||
        error("contact: displacement field does not contain step $step.")

    X = _contact_node_coordinates(problem)
    pdim = problem.pdim

    @inbounds for node in 1:problem.non
        base = (node - 1) * pdim
        for c in 1:pdim
            X[c, node] += Un.a[base + c, step]
        end
    end

    return X
end

function _contact_group_elements(
    problem::Problem,
    phName::String,
    nodecoords::Matrix{Float64};
    aabb_padding::Real
    )

    gmsh.model.setCurrent(problem.name)
    dim_tags = gmsh.model.getEntitiesForPhysicalName(phName)

    isempty(dim_tags) &&
        error("contact: physical group '$phName' was not found.")

    dims = unique(Int(dt[1]) for dt in dim_tags)
    length(dims) == 1 ||
        error(
            "contact: all entities of physical group '$phName' must have " *
            "the same dimension."
        )

    edim = first(dims)
    edim in (1, 2) ||
        error(
            "contact: physical group '$phName' must contain curve (1D) " *
            "or surface (2D) elements."
        )

    elements = _ContactElement[]
    seen = Set{Int}()

    for (_, entity_tag0) in dim_tags
        entity_tag = Int(entity_tag0)
        elem_types, elem_tags, elem_node_tags =
            gmsh.model.mesh.getElements(edim, entity_tag)

        for it in eachindex(elem_types)
            etype = Int(elem_types[it])
            name, elem_dim, order, num_nodes, local_node_coord, _ =
                gmsh.model.mesh.getElementProperties(etype)

            Int(elem_dim) == edim ||
                error("contact: inconsistent Gmsh element dimension for '$name'.")

            local_nodes = _contact_local_node_matrix(
                local_node_coord,
                edim,
                Int(num_nodes)
            )

            tags = elem_tags[it]
            conn = elem_node_tags[it]

            @inbounds for j in eachindex(tags)
                elem_tag = Int(tags[j])
                elem_tag in seen && continue
                push!(seen, elem_tag)

                first_idx = (j - 1) * Int(num_nodes) + 1
                last_idx = j * Int(num_nodes)
                nodes = Int.(conn[first_idx:last_idx])

                X = Matrix{Float64}(undef, 3, Int(num_nodes))
                for a in 1:Int(num_nodes)
                    X[:, a] .= @view nodecoords[:, nodes[a]]
                end

                bmin = vec(minimum(X, dims=2))
                bmax = vec(maximum(X, dims=2))

                diag = norm(bmax - bmin)
                pad = Float64(aabb_padding) * max(diag, eps(Float64))
                bmin .-= pad
                bmax .+= pad

                push!(
                    elements,
                    _ContactElement(
                        elem_tag,
                        etype,
                        edim,
                        String(name),
                        Int(order),
                        nodes,
                        X,
                        local_nodes,
                        bmin,
                        bmax
                    )
                )
            end
        end
    end

    return elements, edim
end

function _contact_local_node_matrix(localCoord, dim::Int, num_nodes::Int)
    data = Float64.(localCoord)

    if length(data) == dim * num_nodes
        return reshape(data, dim, num_nodes)
    elseif length(data) == 3 * num_nodes
        return reshape(data, 3, num_nodes)[1:dim, :]
    else
        error(
            "contact: unexpected number of local node coordinates " *
            "($(length(data))) for a $dim-D, $num_nodes-node element."
        )
    end
end

function _contact_check_master_element_types(elements::Vector{_ContactElement})
    for e in elements
        if e.dim == 1
            occursin("Line", e.name) ||
                error(
                    "contact: unsupported master curve element '$(e.name)'. " *
                    "Only Gmsh Lagrange line elements are currently supported."
                )
        elseif e.dim == 2
            (occursin("Triangle", e.name) || occursin("Quadrilateral", e.name)) ||
                error(
                    "contact: unsupported master surface element '$(e.name)'. " *
                    "Only Gmsh Lagrange triangle and quadrilateral elements " *
                    "are currently supported."
                )
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# AABB tree
# -----------------------------------------------------------------------------

"""
    _contact_build_aabb_tree(elements, indices; leaf_size=8)

Build a binary AABB tree over master elements.
"""
function _contact_build_aabb_tree(
    elements::Vector{_ContactElement},
    indices::Vector{Int};
    leaf_size::Int=8
    )

    isempty(indices) && error("contact: cannot build an AABB tree from no elements.")

    bmin = fill(Inf, 3)
    bmax = fill(-Inf, 3)

    @inbounds for i in indices
        bmin .= min.(bmin, elements[i].bmin)
        bmax .= max.(bmax, elements[i].bmax)
    end

    if length(indices) <= leaf_size
        return _ContactAABBNode(bmin, bmax, nothing, nothing, copy(indices))
    end

    span = bmax - bmin
    axis = argmax(span)

    sorted_indices = sort(
        indices;
        by=i -> 0.5 * (elements[i].bmin[axis] + elements[i].bmax[axis])
    )

    mid = length(sorted_indices) ÷ 2
    left_indices = sorted_indices[1:mid]
    right_indices = sorted_indices[mid+1:end]

    if isempty(left_indices) || isempty(right_indices)
        return _ContactAABBNode(bmin, bmax, nothing, nothing, copy(indices))
    end

    left = _contact_build_aabb_tree(
        elements,
        collect(left_indices);
        leaf_size=leaf_size
    )

    right = _contact_build_aabb_tree(
        elements,
        collect(right_indices);
        leaf_size=leaf_size
    )

    return _ContactAABBNode(bmin, bmax, left, right, Int[])
end

@inline function _contact_point_aabb_distance2(x, bmin, bmax)
    d2 = 0.0
    @inbounds for k in 1:3
        if x[k] < bmin[k]
            d = bmin[k] - x[k]
            d2 += d * d
        elseif x[k] > bmax[k]
            d = x[k] - bmax[k]
            d2 += d * d
        end
    end
    return d2
end

# -----------------------------------------------------------------------------
# Closest-point search
# -----------------------------------------------------------------------------

function _contact_nearest_projection(
    tree::_ContactAABBNode,
    elements::Vector{_ContactElement},
    xs;
    excluded::Union{Nothing,Set{Int}},
    normal_sign::Float64,
    model_dim::Int,
    projection_tol::Float64,
    projection_maxiter::Int
    )

    best = Ref{Union{Nothing,_ContactProjection}}(nothing)
    best_d2 = Ref(Inf)

    _contact_search_aabb!(
        tree,
        elements,
        xs,
        excluded,
        best,
        best_d2,
        normal_sign,
        model_dim,
        projection_tol,
        projection_maxiter
    )

    return best[]
end

function _contact_search_aabb!(
    node::_ContactAABBNode,
    elements::Vector{_ContactElement},
    xs,
    excluded::Union{Nothing,Set{Int}},
    best::Base.RefValue{Union{Nothing,_ContactProjection}},
    best_d2::Base.RefValue{Float64},
    normal_sign::Float64,
    model_dim::Int,
    projection_tol::Float64,
    projection_maxiter::Int
    )

    _contact_point_aabb_distance2(xs, node.bmin, node.bmax) > best_d2[] &&
        return nothing

    if node.left === nothing && node.right === nothing
        @inbounds for idx in node.elements
            excluded !== nothing && idx in excluded && continue

            elem = elements[idx]
            p = _contact_project_element(
                elem,
                xs;
                normal_sign=normal_sign,
                model_dim=model_dim,
                projection_tol=projection_tol,
                projection_maxiter=projection_maxiter,
                element_index=idx
            )

            if p.distance2 < best_d2[]
                best_d2[] = p.distance2
                best[] = p
            end
        end
        return nothing
    end

    left = node.left
    right = node.right

    if left === nothing
        _contact_search_aabb!(
            right, elements, xs, excluded, best, best_d2,
            normal_sign, model_dim, projection_tol, projection_maxiter
        )
        return nothing
    elseif right === nothing
        _contact_search_aabb!(
            left, elements, xs, excluded, best, best_d2,
            normal_sign, model_dim, projection_tol, projection_maxiter
        )
        return nothing
    end

    dl = _contact_point_aabb_distance2(xs, left.bmin, left.bmax)
    dr = _contact_point_aabb_distance2(xs, right.bmin, right.bmax)

    first_node, second_node = dl <= dr ? (left, right) : (right, left)

    _contact_search_aabb!(
        first_node, elements, xs, excluded, best, best_d2,
        normal_sign, model_dim, projection_tol, projection_maxiter
    )

    _contact_search_aabb!(
        second_node, elements, xs, excluded, best, best_d2,
        normal_sign, model_dim, projection_tol, projection_maxiter
    )

    return nothing
end

# -----------------------------------------------------------------------------
# Parametric element projection
# -----------------------------------------------------------------------------

"""
    _contact_project_element(element, xs; ...) -> _ContactProjection

Compute the closest point of `xs` on one master element using projected
Gauss-Newton iterations in the reference element. The local contact basis is
constructed at the converged master point.
"""
function _contact_project_element(
    element::_ContactElement,
    xs;
    normal_sign::Float64,
    model_dim::Int,
    projection_tol::Float64,
    projection_maxiter::Int,
    element_index::Int
    )

    starts = _contact_projection_starts(element)

    best_d2 = Inf
    best_ξ = Vector{Float64}()
    best_x = zeros(3)
    best_N = Float64[]
    best_J = zeros(3, element.dim)

    for ξ0 in starts
        ξ, x, N, J, d2 = _contact_project_from_start(
            element,
            xs,
            ξ0;
            tol=projection_tol,
            maxiter=projection_maxiter
        )

        if d2 < best_d2
            best_d2 = d2
            best_ξ = ξ
            best_x = x
            best_N = N
            best_J = J
        end
    end

    isfinite(best_d2) ||
        error(
            "contact: closest-point projection failed on master element " *
            "$(element.tag) ($(element.name))."
        )

    separation = Vector{Float64}(xs .- best_x)
    normal = _contact_normal(
        element,
        best_J,
        separation,
        model_dim,
        normal_sign
    )

    tangent1, tangent2 = _contact_tangent_basis(
        element,
        best_J,
        normal,
        model_dim
    )

    gap = dot(separation, normal)

    return _ContactProjection(
        element_index,
        element.tag,
        best_ξ,
        best_x,
        best_N,
        normal,
        tangent1,
        tangent2,
        gap,
        best_d2
    )
end

function _contact_projection_starts(element::_ContactElement)
    starts = Vector{Vector{Float64}}()

    if element.dim == 1
        push!(starts, [0.0])
    elseif occursin("Triangle", element.name)
        push!(starts, [1 / 3, 1 / 3])
        push!(starts, [0.5, 0.0])
        push!(starts, [0.5, 0.5])
        push!(starts, [0.0, 0.5])
    elseif occursin("Quadrilateral", element.name)
        push!(starts, [0.0, 0.0])
        push!(starts, [-1.0, 0.0])
        push!(starts, [1.0, 0.0])
        push!(starts, [0.0, -1.0])
        push!(starts, [0.0, 1.0])
    end

    for a in axes(element.local_nodes, 2)
        push!(starts, Vector{Float64}(element.local_nodes[:, a]))
    end

    return starts
end

function _contact_project_from_start(
    element::_ContactElement,
    xs,
    ξ0::Vector{Float64};
    tol::Float64,
    maxiter::Int
    )

    ξ = _contact_project_reference(element, ξ0)
    x, N, J = _contact_geometry(element, ξ)
    r = x - xs
    f = 0.5 * dot(r, r)

    scale = max(norm(element.bmax - element.bmin), 1.0)
    gtol = tol * scale^2
    ξtol = tol

    for _ in 1:maxiter
        g = J' * r
        H = J' * J

        norm(g) <= gtol && break

        # Small diagonal regularization keeps the local solve stable near
        # degenerate parameter directions without changing the converged point.
        reg = max(opnorm(H, Inf), 1.0) * 1e-14
        Hreg = H + reg * I

        δ = try
            -(Hreg \ g)
        catch
            break
        end

        all(isfinite, δ) || break

        accepted = false
        α = 1.0

        for _ in 1:12
            ξtrial = _contact_project_reference(element, ξ + α * δ)

            if norm(ξtrial - ξ) <= ξtol
                ξ = ξtrial
                accepted = true
                break
            end

            xtrial, Ntrial, Jtrial = _contact_geometry(element, ξtrial)
            rtrial = xtrial - xs
            ftrial = 0.5 * dot(rtrial, rtrial)

            if ftrial <= f
                ξ = ξtrial
                x = xtrial
                N = Ntrial
                J = Jtrial
                r = rtrial
                f = ftrial
                accepted = true
                break
            end

            α *= 0.5
        end

        accepted || break

        # Refresh geometry if the accepted step only changed the projected
        # reference coordinate by less than ξtol.
        x, N, J = _contact_geometry(element, ξ)
        r = x - xs
        f = 0.5 * dot(r, r)

        norm(α * δ) <= ξtol && break
    end

    return ξ, x, N, J, 2f
end

function _contact_geometry(element::_ContactElement, ξ::Vector{Float64})
    lc = zeros(Float64, 3)
    lc[1:element.dim] .= ξ

    _, basis, _ =
        gmsh.model.mesh.getBasisFunctions(element.etype, lc, "Lagrange")

    N = Float64.(basis)
    length(N) == length(element.node_tags) ||
        error(
            "contact: unexpected number of Lagrange basis functions on " *
            "element $(element.tag)."
        )

    _, grad, _ =
        gmsh.model.mesh.getBasisFunctions(element.etype, lc, "GradLagrange")

    dN3 = reshape(Float64.(grad), 3, length(element.node_tags))
    dN = @view dN3[1:element.dim, :]

    x = element.coords * N
    J = element.coords * transpose(dN)

    return Vector{Float64}(x), N, Matrix{Float64}(J)
end

function _contact_project_reference(
    element::_ContactElement,
    ξ::AbstractVector
    )

    if element.dim == 1
        return [clamp(Float64(ξ[1]), -1.0, 1.0)]
    elseif occursin("Quadrilateral", element.name)
        return [
            clamp(Float64(ξ[1]), -1.0, 1.0),
            clamp(Float64(ξ[2]), -1.0, 1.0)
        ]
    elseif occursin("Triangle", element.name)
        return _contact_project_triangle_reference(Float64(ξ[1]), Float64(ξ[2]))
    else
        error("contact: unsupported reference element '$(element.name)'.")
    end
end

function _contact_project_triangle_reference(u::Float64, v::Float64)
    if u >= 0 && v >= 0 && u + v <= 1
        return [u, v]
    end

    candidates = Vector{Vector{Float64}}(undef, 3)
    candidates[1] = [0.0, clamp(v, 0.0, 1.0)]
    candidates[2] = [clamp(u, 0.0, 1.0), 0.0]

    t = clamp((u - v + 1.0) / 2.0, 0.0, 1.0)
    candidates[3] = [t, 1.0 - t]

    best = candidates[1]
    best_d2 = (u - best[1])^2 + (v - best[2])^2

    for i in 2:3
        c = candidates[i]
        d2 = (u - c[1])^2 + (v - c[2])^2
        if d2 < best_d2
            best = c
            best_d2 = d2
        end
    end

    return best
end

function _contact_normal(
    element::_ContactElement,
    J::Matrix{Float64},
    separation::Vector{Float64},
    model_dim::Int,
    normal_sign::Float64
    )

    if element.dim == 2
        a1 = @view J[:, 1]
        a2 = @view J[:, 2]
        n = cross(a1, a2)
        nrm = norm(n)

        nrm > sqrt(eps(Float64)) ||
            error(
                "contact: degenerate master surface tangent basis on " *
                "element $(element.tag)."
            )

        n ./= nrm
        n .*= normal_sign
        return Vector{Float64}(n)
    end

    t = Vector{Float64}(@view J[:, 1])
    tnorm = norm(t)
    tnorm > sqrt(eps(Float64)) ||
        error(
            "contact: degenerate master curve tangent on element " *
            "$(element.tag)."
        )
    t ./= tnorm

    # Planar mechanics: the curve orientation defines a signed normal.
    if model_dim == 2
        n = [-t[2], t[1], 0.0]
        n .*= normal_sign
        return n
    end

    # Space curve: the normal plane is not unique. The closest-point
    # separation supplies the physically relevant contact direction.
    dnorm = norm(separation)
    if dnorm > sqrt(eps(Float64))
        n = separation / dnorm
        n .*= normal_sign
        return Vector{Float64}(n)
    end

    # Coincident points: choose a numerically stable vector in the normal plane.
    if abs(t[1]) <= abs(t[2]) && abs(t[1]) <= abs(t[3])
        axis = [1.0, 0.0, 0.0]
    elseif abs(t[2]) <= abs(t[3])
        axis = [0.0, 1.0, 0.0]
    else
        axis = [0.0, 0.0, 1.0]
    end

    n = cross(t, axis)
    n ./= norm(n)
    n .*= normal_sign
    return Vector{Float64}(n)
end

# -----------------------------------------------------------------------------
# Self-contact topology exclusion
# -----------------------------------------------------------------------------

function _contact_self_exclusion_sets(
    slave_nodes::Vector{Int},
    master_elements::Vector{_ContactElement},
    layers::Int
    )

    node_to_elements = Dict{Int,Vector{Int}}()

    for (idx, e) in enumerate(master_elements)
        for node in e.node_tags
            push!(get!(node_to_elements, node, Int[]), idx)
        end
    end

    result = Dict{Int,Set{Int}}()

    for slave_node in slave_nodes
        excluded = Set{Int}(get(node_to_elements, slave_node, Int[]))
        frontier = copy(excluded)

        for _ in 1:layers
            isempty(frontier) && break
            next_frontier = Set{Int}()

            for elem_idx in frontier
                e = master_elements[elem_idx]
                for node in e.node_tags
                    for neighbour in get(node_to_elements, node, Int[])
                        if !(neighbour in excluded)
                            push!(excluded, neighbour)
                            push!(next_frontier, neighbour)
                        end
                    end
                end
            end

            frontier = next_frontier
        end

        result[slave_node] = excluded
    end

    return result
end


"""
    _contact_tangent_basis(element, J, n, model_dim) -> t1, t2

Construct an orthonormal tangential basis associated with the contact normal.
For 2D contact only `t1` is used and `t2 === nothing`. In 3D the tangential
response is isotropic, so the particular in-plane orientation is immaterial for
`C`, while the basis is still stored for slip tracking and post-processing.
"""
function _contact_tangent_basis(
    element::_ContactElement,
    J::Matrix{Float64},
    n::Vector{Float64},
    model_dim::Int
    )

    if model_dim == 2
        t = [-n[2], n[1], 0.0]
        tnorm = norm(t)
        tnorm > sqrt(eps(Float64)) ||
            error("contact: failed to construct a planar contact tangent.")
        t ./= tnorm

        # Keep the tangent consistent with the master parameter direction.
        if element.dim >= 1
            a = Vector{Float64}(@view J[:, 1])
            if norm(a) > sqrt(eps(Float64)) && dot(t, a) < 0
                t .*= -1.0
            end
        end

        return t, nothing
    end

    # 3D surface: use the first covariant direction and orthogonalize it.
    if element.dim == 2
        t1 = Vector{Float64}(@view J[:, 1])
        t1 .-= dot(t1, n) .* n
        if norm(t1) <= sqrt(eps(Float64))
            t1 = Vector{Float64}(@view J[:, 2])
            t1 .-= dot(t1, n) .* n
        end
    else
        # 3D curve: one tangent is supplied by the curve geometry.
        t1 = Vector{Float64}(@view J[:, 1])
        t1 .-= dot(t1, n) .* n
    end

    t1norm = norm(t1)
    t1norm > sqrt(eps(Float64)) ||
        error(
            "contact: failed to construct the first tangential direction on " *
            "master element $(element.tag)."
        )
    t1 ./= t1norm

    t2 = cross(n, t1)
    t2norm = norm(t2)
    t2norm > sqrt(eps(Float64)) ||
        error(
            "contact: failed to construct the second tangential direction on " *
            "master element $(element.tag)."
        )
    t2 ./= t2norm

    # Re-orthogonalize t1 against the final n-t2 pair.
    t1 = cross(t2, n)
    t1 ./= norm(t1)

    return Vector{Float64}(t1), Vector{Float64}(t2)
end

# -----------------------------------------------------------------------------
# LLFEM output objects and reduced contact matrices
# -----------------------------------------------------------------------------

function _contact_gap_field(
    field_model::Problem,
    slave_elements::Vector{_ContactElement},
    projections::Dict{Int,_ContactProjection}
    )

    A = Vector{Matrix{Float64}}(undef, length(slave_elements))
    num_elem = Vector{Int}(undef, length(slave_elements))

    @inbounds for (i, e) in enumerate(slave_elements)
        values = Matrix{Float64}(undef, length(e.node_tags), 1)

        for (a, node) in enumerate(e.node_tags)
            values[a, 1] = projections[node].gap
        end

        A[i] = values
        num_elem[i] = e.tag
    end

    return ScalarField(
        A,
        [;;],
        [0.0],
        num_elem,
        1,
        :scalar,
        field_model
    )
end

function _contact_vector_field(
    model::Problem,
    slave_elements::Vector{_ContactElement},
    projections::Dict{Int,_ContactProjection},
    component::Symbol
    )

    component in (:normal, :tangent1, :tangent2) ||
        error("contact: unknown contact vector component '$component'.")

    A = Vector{Matrix{Float64}}(undef, length(slave_elements))
    num_elem = Vector{Int}(undef, length(slave_elements))

    @inbounds for (i, e) in enumerate(slave_elements)
        values = Matrix{Float64}(undef, 3 * length(e.node_tags), 1)

        for (a, node) in enumerate(e.node_tags)
            p = projections[node]
            v = if component === :normal
                p.normal
            elseif component === :tangent1
                p.tangent1
            else
                p.tangent2 === nothing &&
                    error("contact: tangent2 is not defined for this contact dimension.")
                p.tangent2
            end

            values[3a - 2, 1] = v[1]
            values[3a - 1, 1] = v[2]
            values[3a, 1] = v[3]
        end

        A[i] = values
        num_elem[i] = e.tag
    end

    return VectorField(
        A,
        [;;],
        [0.0],
        num_elem,
        1,
        :v3D,
        model
    )
end

"""
    _contact_matrix(model, slave_nodes, projections, master_elements)

Build the reduced local contact kinematic operator. For each slave node the row
ordering is `(n,t)` in 2D and `(n,t1,t2)` in 3D.
"""
function _contact_matrix(
    model::Problem,
    slave_nodes::Vector{Int},
    projections::Dict{Int,_ContactProjection},
    master_elements::Vector{_ContactElement}
    )

    pdim = model.pdim
    nrows = length(slave_nodes) * pdim
    ncols = ndofs(model)

    Iidx = Int[]
    Jidx = Int[]
    Vval = Float64[]

    estimated = 0
    for node in slave_nodes
        p = projections[node]
        estimated += pdim * pdim *
            (1 + length(master_elements[p.element_index].node_tags))
    end
    sizehint!(Iidx, estimated)
    sizehint!(Jidx, estimated)
    sizehint!(Vval, estimated)

    @inbounds for (i, slave_node) in enumerate(slave_nodes)
        p = projections[slave_node]
        e = master_elements[p.element_index]

        basis = if pdim == 2
            (p.normal, p.tangent1)
        else
            p.tangent2 === nothing &&
                error("contact: missing second tangent in 3D contact.")
            (p.normal, p.tangent1, p.tangent2)
        end

        for α in 1:pdim
            row = (i - 1) * pdim + α
            q = basis[α]

            for c in 1:pdim
                val = q[c]
                iszero(val) && continue
                push!(Iidx, row)
                push!(Jidx, (slave_node - 1) * pdim + c)
                push!(Vval, val)
            end

            for (a, master_node) in enumerate(e.node_tags)
                Na = p.N[a]
                iszero(Na) && continue

                for c in 1:pdim
                    val = -Na * q[c]
                    iszero(val) && continue
                    push!(Iidx, row)
                    push!(Jidx, (master_node - 1) * pdim + c)
                    push!(Vval, val)
                end
            end
        end
    end

    G = sparse(Iidx, Jidx, Vval, nrows, ncols)
    dropzeros!(G)
    return G
end

function _contact_stiffness_matrix(
    pdim::Int,
    active::BitVector,
    cn::Vector{Float64},
    ct::Vector{Float64}
    )

    nc = length(active)
    length(cn) == nc || error("contact: cn size mismatch.")
    length(ct) == nc || error("contact: ct size mismatch.")

    nrows = nc * pdim
    Iidx = Int[]
    Jidx = Int[]
    Vval = Float64[]
    sizehint!(Iidx, nrows)
    sizehint!(Jidx, nrows)
    sizehint!(Vval, nrows)

    @inbounds for i in 1:nc
        active[i] || continue

        row_n = (i - 1) * pdim + 1
        if !iszero(cn[i])
            push!(Iidx, row_n)
            push!(Jidx, row_n)
            push!(Vval, cn[i])
        end

        for α in 2:pdim
            iszero(ct[i]) && continue
            row_t = (i - 1) * pdim + α
            push!(Iidx, row_t)
            push!(Jidx, row_t)
            push!(Vval, ct[i])
        end
    end

    C = sparse(Iidx, Jidx, Vval, nrows, nrows)
    dropzeros!(C)
    return C
end

function _contact_multiplier_dofs(
    multiplier::Union{Nothing,Problem},
    displacement::Problem,
    slave_nodes::Vector{Int}
    )

    multiplier === nothing && return Int[]

    pdim = displacement.pdim
    dofs = Vector{Int}(undef, length(slave_nodes) * pdim)

    @inbounds for (i, node) in enumerate(slave_nodes)
        for α in 1:pdim
            dofs[(i - 1) * pdim + α] = (node - 1) * pdim + α
        end
    end

    return dofs
end

function _contact_parameter_values(
    parameter,
    slave_nodes::Vector{Int},
    nodecoords::Matrix{Float64},
    model::Problem;
    step::Int,
    name::String
    )

    values = Vector{Float64}(undef, length(slave_nodes))

    if parameter isa Number
        value = Float64(parameter)
        isfinite(value) || error("contact: $name must be finite.")
        fill!(values, value)
        return values
    end

    if parameter isa ScalarField
        parameter.model.name == model.name ||
            error("contact: $name ScalarField must use the same Gmsh model.")

        S = isNodal(parameter) ? parameter : elementsToNodes(parameter)
        sstep = S.nsteps == 1 ? 1 : step
        1 <= sstep <= S.nsteps ||
            error("contact: $name ScalarField does not contain step $step.")

        @inbounds for (i, node) in enumerate(slave_nodes)
            values[i] = S.a[node, sstep]
            isfinite(values[i]) || error("contact: $name contains a non-finite value.")
        end
        return values
    end

    if parameter isa Function
        @inbounds for (i, node) in enumerate(slave_nodes)
            x = nodecoords[1, node]
            y = nodecoords[2, node]
            z = nodecoords[3, node]
            values[i] = Float64(parameter(x, y, z))
            isfinite(values[i]) || error("contact: $name function returned a non-finite value.")
        end
        return values
    end

    error(
        "contact: $name must be a Number, ScalarField, or function f(x,y,z); " *
        "got $(typeof(parameter))."
    )
end
