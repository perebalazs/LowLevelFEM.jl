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

const _CONTACT_REF_LINE = UInt8(1)
const _CONTACT_REF_TRI  = UInt8(2)
const _CONTACT_REF_QUAD = UInt8(3)

"""
Cached polynomial representation of one Gmsh Lagrange basis.

The cache is built once per Gmsh element type. Gmsh is only used while the
cache is constructed; all evaluations inside the closest-point iteration are
performed locally from the cached polynomial coefficients.
"""
struct _ContactBasisCache
    etype::Int
    dim::Int
    name::String
    order::Int
    num_nodes::Int
    refkind::UInt8
    exp_u::Vector{Int}
    exp_v::Vector{Int}
    coeff::Matrix{Float64}                 # monomial coefficients × shape functions
    max_u::Int
    max_v::Int
    starts::Vector{NTuple{2,Float64}}
    use_gmsh::Bool
end

const _CONTACT_BASIS_CACHE = Dict{Int,_ContactBasisCache}()
const _CONTACT_BASIS_CACHE_LOCK = ReentrantLock()

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
    basis::_ContactBasisCache
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

mutable struct _ContactProjectionWorkspace
    N::Vector{Float64}
    dN::Matrix{Float64}
    Ntrial::Vector{Float64}
    upow::Vector{Float64}
    vpow::Vector{Float64}
    x::Vector{Float64}
    J::Matrix{Float64}
    xtrial::Vector{Float64}
    localcoord::Vector{Float64}
end

mutable struct _ContactSearchResult
    element_index::Int
    u::Float64
    v::Float64
    distance2::Float64
end

# -----------------------------------------------------------------------------
# Public contact object
# -----------------------------------------------------------------------------

"""
    Contact

Container for one slave-master contact pair.

`U` is the displacement `Problem` defining the global displacement DoFs, while
`displacement` is the current `VectorField` used to construct the deformed
geometry `x = X + displacement`.

The contact kinematics are stored in a reduced local contact space. Every slave
contact node contributes `U.pdim` rows to `G`: `(normal,tangent)` in 2D and
`(normal,tangent1,tangent2)` in 3D. Therefore, for `nc` slave contact nodes,

    size(G) == (nc * U.pdim, ndofs(U))

and the penalty contribution has the direct mathematical form

    Kc = G' * C * G

where `C` is block diagonal in the local contact basis. Isotropic tangential
behaviour is represented by the same `ct` in both tangential directions in 3D.
Inactive contact points have zero local blocks in `C`.

If a Lagrange multiplier field is supplied, `multiplier_dofs` maps the reduced
rows of `G` to the corresponding multiplier DoFs on the slave nodes. This
mapping is stored for the later `solveField` implementation.
"""
mutable struct Contact
    master::String
    slave::String

    U::Problem
    multiplier::Union{Nothing,Problem}
    displacement::VectorField

    slave_nodes::Vector{Int}
    master_element_tags::Vector{Int}
    master_local_coordinates::Vector{Vector{Float64}}
    master_points::Matrix{Float64}

    gap::ScalarField
    gap_values::Vector{Float64}
    G::SparseMatrixCSC{Float64,Int}
    C::SparseMatrixCSC{Float64,Int}

    n::VectorField
    t1::VectorField
    t2::Union{Nothing,VectorField}

    active::BitVector
    state::Vector{UInt8}

    # Local contact-space history vectors reserved for the Coulomb update.
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
    contact(U::Problem; master, slave, displacement=zero_displacement,
            LagrangeMultiplierField=nothing, cn=1.0, ct=0.0, μ=0.0,
            activation_tol=0.0, step=displacement.nsteps, kwargs...) -> Contact

Construct one node-to-manifold contact pair in the current configuration.

`U` is the displacement problem and `displacement` is the current displacement
field. The geometry used by the contact search is

    x = X + displacement

on both the slave and master sides. Hence the closest-point projection, normal,
tangential basis, signed gap and kinematic matrix `G` are all evaluated in the
current configuration.

# Arguments

- `U::Problem`: displacement problem associated with the columns of `G`.
- `master::String`: master physical group.
- `slave::String`: slave physical group. Contact quantities are discretized on
  the slave nodes.
- `displacement::VectorField`: current displacement field. If omitted, a zero
  nodal displacement field is used.
- `LagrangeMultiplierField::Union{Nothing,Problem}=nothing`: optional vector
  multiplier problem. In 2D it must have two components and in 3D three.
- `cn`: normal penalty stiffness. It may be a number, `ScalarField`, or
  `f(x,y,z)` function.
- `ct`: isotropic tangential penalty stiffness. It may be a number,
  `ScalarField`, or `f(x,y,z)` function. Set `ct=0` for frictionless contact.
- `μ`: isotropic Coulomb friction coefficient reserved for the stick-slip
  update. It may be a number, `ScalarField`, or `f(x,y,z)` function.
- `activation_tol::Real=0.0`: a point is active when `gap <= activation_tol`.
- `step::Int=displacement.nsteps`: displacement step used for the current
  geometry.

# Geometry options

- `normal_sign::Real=1.0`: multiply the master-side normal by `+1` or `-1`.
  The signed gap is `(x_slave - x_master) ⋅ n`.
- `self_contact::Bool=(slave == master)`: enable local-topology exclusion.
- `self_exclusion_layers::Int=1`: number of additional node-connected master
  element layers excluded in self-contact.
- `aabb_padding::Real=0.05`: relative AABB padding.
- `leaf_size::Int=2`: maximum number of elements in an AABB leaf.
- `projection_tol::Real=1e-10`: closest-point solver tolerance.
- `projection_maxiter::Int=40`: maximum projected Gauss-Newton iterations.

# Performance

Gmsh basis functions are sampled only while a cache for a new master element
type is built. Closest-point iterations use a local polynomial evaluator and a
reused workspace, so no Gmsh calls occur in the hot projection loops.
"""
function contact(
    U::Problem;
    master::String,
    slave::String,
    displacement::VectorField=_contact_zero_displacement(U),
    LagrangeMultiplierField::Union{Nothing,Problem}=nothing,
    cn=1.0,
    ct=0.0,
    μ=0.0,
    activation_tol::Real=0.0,
    step::Int=displacement.nsteps,
    normal_sign::Real=1.0,
    self_contact::Bool=(slave == master),
    self_exclusion_layers::Int=1,
    aabb_padding::Real=0.05,
    leaf_size::Int=2,
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
        U,
        LagrangeMultiplierField,
        displacement,
        data.slave_nodes,
        data.master_element_tags,
        data.master_local_coordinates,
        data.master_points,
        data.gap,
        data.gap_values,
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
    contact(displacement::VectorField; master, slave, U=displacement.model, kwargs...)

Compatibility/convenience method using a displacement field as the first
argument. The canonical internal naming remains `U::Problem` and
`displacement::VectorField`.
"""
function contact(
    displacement::VectorField;
    master::String,
    slave::String,
    U::Problem=displacement.model,
    kwargs...
    )

    return contact(
        U;
        master=master,
        slave=slave,
        displacement=displacement,
        kwargs...
    )
end

# Positional compatibility helpers.
contact(displacement::VectorField, slave::String, master::String; kwargs...) =
    contact(displacement; slave=slave, master=master, kwargs...)

contact(U::Problem, slave::String, master::String; kwargs...) =
    contact(U; slave=slave, master=master, kwargs...)

"""
    updateContact!(c::Contact, displacement::VectorField; step=displacement.nsteps)

Recompute the current contact geometry, active set, `G` and `C` from a new
current displacement field. The previous master element and local coordinates
are used as a warm start and as an initial upper bound for the exact AABB
branch-and-bound search. The search therefore remains global and may switch to
any closer master element. Existing traction and slip history arrays are reused
when the slave-node set is unchanged; traction is cleared on newly open points.
"""
function updateContact!(
    c::Contact,
    displacement::VectorField;
    step::Int=displacement.nsteps
    )

    # Keep the previous geometric association as a warm start.  The previous
    # master element and local coordinates provide a cheap current upper bound
    # for the exact AABB search; the global search is still performed, so the
    # result is not restricted to the previous element or its neighbourhood.
    old_nodes = c.slave_nodes
    old_master_element_tags = c.master_element_tags
    old_master_local_coordinates = c.master_local_coordinates

    data = _contact_build_data(
        c.U,
        displacement,
        c.multiplier,
        c.slave,
        c.master,
        c.cn,
        c.ct,
        c.μ;
        step=step,
        previous_slave_nodes=old_nodes,
        previous_master_element_tags=old_master_element_tags,
        previous_master_local_coordinates=old_master_local_coordinates,
        c.options...
    )

    pdim = c.U.pdim
    nrows = size(data.G, 1)
    preserve_history = old_nodes == data.slave_nodes && length(c.traction) == nrows

    c.displacement = displacement
    c.slave_nodes = data.slave_nodes
    c.master_element_tags = data.master_element_tags
    c.master_local_coordinates = data.master_local_coordinates
    c.master_points = data.master_points
    c.gap = data.gap
    c.gap_values = data.gap_values
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

    if preserve_history
        # Reuse the existing history arrays.  This avoids three allocations on
        # every contact update and, more importantly, keeps the history in the
        # same slave-node ordering.
        @inbounds for i in eachindex(c.slave_nodes)
            firstrow = (i - 1) * pdim + 1
            if c.active[i]
                c.state[i] = c.state[i] == CONTACT_OPEN ? CONTACT_STICK : c.state[i]
            else
                c.state[i] = CONTACT_OPEN
                for j in 0:pdim-1
                    c.traction[firstrow + j] = 0.0
                end
            end
        end
    else
        c.traction = zeros(Float64, nrows)
        c.slip = zeros(Float64, nrows)
        c.state = Vector{UInt8}(undef, length(c.slave_nodes))
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

as a `SystemMatrix` associated with the displacement problem `c.U`.
"""
function contactPenaltyMatrix(c::Contact)
    return SystemMatrix(sparse(c.G' * c.C * c.G), c.U)
end

# -----------------------------------------------------------------------------
# Main construction
# -----------------------------------------------------------------------------

function _contact_build_data(
    U::Problem,
    displacement::VectorField,
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
    projection_maxiter::Int,
    previous_slave_nodes::Union{Nothing,Vector{Int}}=nothing,
    previous_master_element_tags::Union{Nothing,Vector{Int}}=nothing,
    previous_master_local_coordinates::Union{Nothing,Vector{Vector{Float64}}}=nothing
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

    gmsh.model.setCurrent(U.name)

    nodecoords = _contact_deformed_coordinates(U, displacement; step=step)

    slave_elements, slave_dim =
        _contact_group_elements(U, slave, nodecoords; aabb_padding=0.0)

    master_elements, master_dim =
        _contact_group_elements(U, master, nodecoords; aabb_padding=aabb_padding)

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

    # One workspace is reused for every slave point and every candidate master
    # element. The current search is serial; a future threaded implementation
    # should allocate one workspace per thread.
    workspace = _contact_projection_workspace(master_elements)
    projections = Dict{Int,_ContactProjection}()
    sizehint!(projections, length(slave_nodes))
    empty_exclusion = Set{Int}()

    #gap_values = zeros(Float64, length(slave_nodes))

    # Prepare optional warm-start lookup tables.  In the normal update path the
    # slave-node ordering is unchanged, so no node dictionary is needed.
    have_previous =
        previous_slave_nodes !== nothing &&
        previous_master_element_tags !== nothing &&
        previous_master_local_coordinates !== nothing &&
        length(previous_slave_nodes) == length(previous_master_element_tags) &&
        length(previous_slave_nodes) == length(previous_master_local_coordinates)

    same_slave_order = have_previous && previous_slave_nodes == slave_nodes
    previous_node_index = nothing
    if have_previous && !same_slave_order
        previous_node_index = Dict{Int,Int}()
        sizehint!(previous_node_index, length(previous_slave_nodes))
        @inbounds for (i, node) in enumerate(previous_slave_nodes)
            previous_node_index[node] = i
        end
    end

    master_tag_to_index = nothing
    if have_previous
        master_tag_to_index = Dict{Int,Int}()
        sizehint!(master_tag_to_index, length(master_elements))
        @inbounds for (i, element) in enumerate(master_elements)
            master_tag_to_index[element.tag] = i
        end
    end

    @inbounds for (i, node) in enumerate(slave_nodes)
        xs = @view nodecoords[:, node]
        ex = self_contact ? get(excluded, node, empty_exclusion) : nothing

        previous_element_index = 0
        previous_local_coordinate = nothing

        if have_previous
            old_i = same_slave_order ? i : get(previous_node_index, node, 0)
            if old_i != 0
                old_tag = previous_master_element_tags[old_i]
                previous_element_index = get(master_tag_to_index, old_tag, 0)
                if previous_element_index != 0 &&
                   (ex === nothing || !(previous_element_index in ex))
                    previous_local_coordinate = previous_master_local_coordinates[old_i]
                else
                    previous_element_index = 0
                end
            end
        end

        p = _contact_nearest_projection(
            tree,
            master_elements,
            xs,
            workspace;
            excluded=ex,
            normal_sign=normal_sign,
            model_dim=U.dim,
            projection_tol=projection_tol,
            projection_maxiter=projection_maxiter,
            previous_element_index=previous_element_index,
            previous_local_coordinate=previous_local_coordinate
        )

        p === nothing &&
            error(
                "contact: no admissible master element was found for " *
                "slave node $node. In self-contact, try reducing " *
                "self_exclusion_layers."
            )

        projections[node] = p
    end

    gap = _contact_gap_field(U, slave_elements, projections)
    normalVec = _contact_vector_field(U, slave_elements, projections, :normal)
    tangent1 = _contact_vector_field(U, slave_elements, projections, :tangent1)
    tangent2 = U.pdim == 3 ?
        _contact_vector_field(U, slave_elements, projections, :tangent2) : nothing

    G = _contact_matrix(
        U,
        slave_nodes,
        projections,
        master_elements
    )

    gap_values = [projections[node].gap for node in slave_nodes]
    active = BitVector(g <= activation_tol for g in gap_values)

    cn_values = _contact_parameter_values(cn, slave_nodes, nodecoords, U; step=step, name="cn")
    ct_values = _contact_parameter_values(ct, slave_nodes, nodecoords, U; step=step, name="ct")
    μ_values  = _contact_parameter_values(μ,  slave_nodes, nodecoords, U; step=step, name="μ")

    any(x -> x < 0.0, cn_values) && error("contact: cn must be non-negative.")
    any(x -> x < 0.0, ct_values) && error("contact: ct must be non-negative.")
    any(x -> x < 0.0, μ_values)  && error("contact: μ must be non-negative.")

    C = _contact_stiffness_matrix(
        U.pdim,
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

    multiplier_dofs = _contact_multiplier_dofs(multiplier, U, slave_nodes)

    return (
        slave_nodes=slave_nodes,
        master_element_tags=master_element_tags,
        master_local_coordinates=master_local_coordinates,
        master_points=master_points,
        gap=gap,
        gap_values = gap_values,
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
    U::Problem,
    displacement::VectorField,
    multiplier::Union{Nothing,Problem}
    )

    U.pdim in (2, 3) ||
        error(
            "contact: U must be a 2D or 3D displacement problem; " *
            "got pdim=$(U.pdim)."
        )

    displacement.model.name == U.name ||
        error("contact: U and displacement must use the same Gmsh model.")

    displacement.model.non == U.non ||
        error("contact: U and displacement must use the same mesh nodes.")

    displacement.model.pdim == U.pdim ||
        error(
            "contact: U and displacement must have the same number of " *
            "components per node."
        )

    if multiplier !== nothing
        multiplier.name == U.name ||
            error(
                "contact: LagrangeMultiplierField and U must use the same " *
                "Gmsh model."
            )

        multiplier.non == U.non ||
            error(
                "contact: LagrangeMultiplierField and U must use the same " *
                "mesh nodes."
            )

        multiplier.pdim == U.pdim ||
            error(
                "contact: the frictional Lagrange multiplier field must have " *
                "$(U.pdim) components per node."
            )
    end

    return nothing
end

function _contact_zero_displacement(U::Problem)
    U.pdim in (2, 3) ||
        error("contact: zero displacement can only be created for pdim=2 or 3.")

    type = U.pdim == 2 ? :v2D : :v3D
    return VectorField(
        Matrix{Float64}[],
        zeros(Float64, ndofs(U), 1),
        [0.0],
        Int[],
        1,
        type,
        U
    )
end

"""
    _contact_node_coordinates(U) -> Matrix{Float64}

Return reference mesh-node coordinates in a `3 × U.non` matrix indexed by Gmsh
node tag.
"""
function _contact_node_coordinates(U::Problem)
    gmsh.model.setCurrent(U.name)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()

    length(node_tags) == U.non ||
        error(
            "contact: Problem node count ($(U.non)) differs from the " *
            "current Gmsh mesh ($(length(node_tags)))."
        )

    X = zeros(Float64, 3, U.non)

    @inbounds for (i, tag0) in enumerate(node_tags)
        tag = Int(tag0)
        1 <= tag <= U.non ||
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
    _contact_deformed_coordinates(U, displacement; step) -> Matrix{Float64}

Return current coordinates `x = X + displacement` for contact search and
projection.
"""
function _contact_deformed_coordinates(
    U::Problem,
    displacement::VectorField;
    step::Int
    )

    1 <= step <= displacement.nsteps ||
        error(
            "contact: displacement step $step is outside " *
            "1:$(displacement.nsteps)."
        )

    un = isNodal(displacement) ? displacement : elementsToNodes(displacement)

    size(un.a, 1) == ndofs(U) ||
        error(
            "contact: nodal displacement contains $(size(un.a, 1)) values " *
            "per step, expected $(ndofs(U))."
        )

    size(un.a, 2) >= step ||
        error("contact: displacement field does not contain step $step.")

    X = _contact_node_coordinates(U)
    pdim = U.pdim

    @inbounds for node in 1:U.non
        base = (node - 1) * pdim
        for c in 1:pdim
            X[c, node] += un.a[base + c, step]
        end
    end

    return X
end

function _contact_group_elements(
    U::Problem,
    phName::String,
    nodecoords::Matrix{Float64};
    aabb_padding::Real
    )

    gmsh.model.setCurrent(U.name)
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
            name0, elem_dim, order0, num_nodes0, local_node_coord, _ =
                gmsh.model.mesh.getElementProperties(etype)

            name = String(name0)
            order = Int(order0)
            num_nodes = Int(num_nodes0)

            Int(elem_dim) == edim ||
                error("contact: inconsistent Gmsh element dimension for '$name'.")

            local_nodes = _contact_local_node_matrix(
                local_node_coord,
                edim,
                num_nodes
            )

            # This is the only place where basis information is obtained from
            # Gmsh. The cache is shared by all elements of the same Gmsh type.
            basis = _contact_get_basis_cache(
                etype,
                edim,
                name,
                order,
                num_nodes,
                local_nodes
            )

            tags = elem_tags[it]
            conn = elem_node_tags[it]

            @inbounds for j in eachindex(tags)
                elem_tag = Int(tags[j])
                elem_tag in seen && continue
                push!(seen, elem_tag)

                first_idx = (j - 1) * num_nodes + 1
                last_idx = j * num_nodes
                nodes = Int.(conn[first_idx:last_idx])

                X = Matrix{Float64}(undef, 3, num_nodes)
                for a in 1:num_nodes
                    node = nodes[a]
                    X[1, a] = nodecoords[1, node]
                    X[2, a] = nodecoords[2, node]
                    X[3, a] = nodecoords[3, node]
                end

                bmin = [Inf, Inf, Inf]
                bmax = [-Inf, -Inf, -Inf]
                for a in 1:num_nodes
                    for k in 1:3
                        x = X[k, a]
                        x < bmin[k] && (bmin[k] = x)
                        x > bmax[k] && (bmax[k] = x)
                    end
                end

                dx = bmax[1] - bmin[1]
                dy = bmax[2] - bmin[2]
                dz = bmax[3] - bmin[3]
                diag = sqrt(dx * dx + dy * dy + dz * dz)
                pad = Float64(aabb_padding) * max(diag, eps(Float64))
                for k in 1:3
                    bmin[k] -= pad
                    bmax[k] += pad
                end

                push!(
                    elements,
                    _ContactElement(
                        elem_tag,
                        etype,
                        edim,
                        name,
                        order,
                        nodes,
                        X,
                        local_nodes,
                        bmin,
                        bmax,
                        basis
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

function _contact_reference_kind(name::String, dim::Int)
    if dim == 1 && occursin("Line", name)
        return _CONTACT_REF_LINE
    elseif dim == 2 && occursin("Triangle", name)
        return _CONTACT_REF_TRI
    elseif dim == 2 && occursin("Quadrilateral", name)
        return _CONTACT_REF_QUAD
    end

    error(
        "contact: unsupported master element '$name'. " *
        "Lagrange line, triangle and quadrilateral elements are supported."
    )
end

function _contact_check_master_element_types(elements::Vector{_ContactElement})
    for e in elements
        _contact_reference_kind(e.name, e.dim)
    end
    return nothing
end

"""
    _contact_get_basis_cache(...)

Return a cached polynomial evaluator for a Gmsh Lagrange element type. The
Gmsh basis is sampled in one batched call while the cache is built; subsequent
contact projections do not call `gmsh.model.mesh.getBasisFunctions`.
"""
function _contact_get_basis_cache(
    etype::Int,
    dim::Int,
    name::String,
    order::Int,
    num_nodes::Int,
    local_nodes::Matrix{Float64}
    )

    lock(_CONTACT_BASIS_CACHE_LOCK)
    try
        if haskey(_CONTACT_BASIS_CACHE, etype)
            return _CONTACT_BASIS_CACHE[etype]
        end

        cache = _contact_build_basis_cache(
            etype,
            dim,
            name,
            order,
            num_nodes,
            local_nodes
        )
        _CONTACT_BASIS_CACHE[etype] = cache
        return cache
    finally
        unlock(_CONTACT_BASIS_CACHE_LOCK)
    end
end

function _contact_build_basis_cache(
    etype::Int,
    dim::Int,
    name::String,
    order::Int,
    num_nodes::Int,
    local_nodes::Matrix{Float64}
    )

    refkind = _contact_reference_kind(name, dim)
    starts = _contact_projection_starts(local_nodes, refkind, dim)

    # Standard Gmsh Lagrange bases are polynomial. For complete line/triangle/
    # quadrilateral elements the nominal order is sufficient. The small
    # adaptive loop also covers incomplete/serendipity variants without
    # putting Gmsh calls into the projection loop.
    max_degree = max(order + 3, 2 * order + 1)
    best_residual = Inf

    for degree in order:max_degree
        exp_u, exp_v = _contact_basis_exponents(refkind, degree)
        sample_u, sample_v = _contact_basis_sample_points(refkind, degree)
        ns = length(sample_u)
        nm = length(exp_u)

        ns >= nm || continue

        local_coord = Vector{Float64}(undef, 3 * ns)
        @inbounds for q in 1:ns
            local_coord[3q - 2] = sample_u[q]
            local_coord[3q - 1] = sample_v[q]
            local_coord[3q] = 0.0
        end

        num_components, basis_raw, _ =
            gmsh.model.mesh.getBasisFunctions(etype, local_coord, "Lagrange")

        Int(num_components) == 1 ||
            error(
                "contact: expected scalar Lagrange basis for element type " *
                "$etype, got $num_components components."
            )

        basis_data = Float64.(basis_raw)
        length(basis_data) == ns * num_nodes ||
            error(
                "contact: unexpected number of basis values while caching " *
                "element '$name' (type $etype)."
            )

        # Gmsh stores all basis values of point 1 first, then point 2, ...
        B = Matrix(transpose(reshape(basis_data, num_nodes, ns)))
        Φ = Matrix{Float64}(undef, ns, nm)
        _contact_fill_monomial_matrix!(Φ, sample_u, sample_v, exp_u, exp_v)

        coeff = Φ \ B
        residual = norm(Φ * coeff - B, Inf) / max(norm(B, Inf), 1.0)
        best_residual = min(best_residual, residual)

        if isfinite(residual) && residual <= 1e-9
            return _ContactBasisCache(
                etype,
                dim,
                name,
                order,
                num_nodes,
                refkind,
                exp_u,
                exp_v,
                coeff,
                maximum(exp_u),
                maximum(exp_v),
                starts,
                false
            )
        end
    end

    @warn(
        "contact: polynomial basis cache verification failed for '$name' " *
        "(Gmsh type $etype, order $order; best relative residual = " *
        "$(best_residual)). Falling back to direct Gmsh basis evaluation " *
        "for this element type. Generality is preserved, but projection " *
        "will be slower for this type."
    )

    return _ContactBasisCache(
        etype,
        dim,
        name,
        order,
        num_nodes,
        refkind,
        Int[],
        Int[],
        zeros(Float64, 0, num_nodes),
        0,
        0,
        starts,
        true
    )
end

function _contact_basis_exponents(refkind::UInt8, degree::Int)
    exp_u = Int[]
    exp_v = Int[]

    if refkind == _CONTACT_REF_LINE
        sizehint!(exp_u, degree + 1)
        sizehint!(exp_v, degree + 1)
        for a in 0:degree
            push!(exp_u, a)
            push!(exp_v, 0)
        end
    elseif refkind == _CONTACT_REF_TRI
        n = (degree + 1) * (degree + 2) ÷ 2
        sizehint!(exp_u, n)
        sizehint!(exp_v, n)
        for total in 0:degree
            for a in 0:total
                push!(exp_u, a)
                push!(exp_v, total - a)
            end
        end
    elseif refkind == _CONTACT_REF_QUAD
        n = (degree + 1)^2
        sizehint!(exp_u, n)
        sizehint!(exp_v, n)
        for b in 0:degree
            for a in 0:degree
                push!(exp_u, a)
                push!(exp_v, b)
            end
        end
    else
        error("contact: unknown reference-element kind.")
    end

    return exp_u, exp_v
end

function _contact_basis_sample_points(refkind::UInt8, degree::Int)
    # Use an overdetermined set. This makes the cache construction robust and
    # simultaneously verifies that the selected polynomial space reproduces
    # the Gmsh basis.
    q = max(2 * degree + 3, 5)
    u = Float64[]
    v = Float64[]

    if refkind == _CONTACT_REF_LINE
        sizehint!(u, q)
        sizehint!(v, q)
        for i in 0:q-1
            s = -1.0 + 2.0 * i / (q - 1)
            push!(u, s)
            push!(v, 0.0)
        end
    elseif refkind == _CONTACT_REF_QUAD
        sizehint!(u, q * q)
        sizehint!(v, q * q)
        for j in 0:q-1
            η = -1.0 + 2.0 * j / (q - 1)
            for i in 0:q-1
                ξ = -1.0 + 2.0 * i / (q - 1)
                push!(u, ξ)
                push!(v, η)
            end
        end
    elseif refkind == _CONTACT_REF_TRI
        n = q * (q + 1) ÷ 2
        sizehint!(u, n)
        sizehint!(v, n)
        m = q - 1
        for j in 0:m
            for i in 0:m-j
                push!(u, i / m)
                push!(v, j / m)
            end
        end
    else
        error("contact: unknown reference-element kind.")
    end

    return u, v
end

function _contact_fill_monomial_matrix!(
    Φ::Matrix{Float64},
    u::Vector{Float64},
    v::Vector{Float64},
    exp_u::Vector{Int},
    exp_v::Vector{Int}
    )

    @inbounds for q in eachindex(u)
        uq = u[q]
        vq = v[q]
        for k in eachindex(exp_u)
            Φ[q, k] = uq^exp_u[k] * vq^exp_v[k]
        end
    end
    return Φ
end

function _contact_projection_starts(
    local_nodes::Matrix{Float64},
    refkind::UInt8,
    dim::Int
    )

    starts = NTuple{2,Float64}[]

    if refkind == _CONTACT_REF_LINE
        push!(starts, (0.0, 0.0))
    elseif refkind == _CONTACT_REF_TRI
        push!(starts, (1 / 3, 1 / 3))
        push!(starts, (0.5, 0.0))
        push!(starts, (0.5, 0.5))
        push!(starts, (0.0, 0.5))
    elseif refkind == _CONTACT_REF_QUAD
        push!(starts, (0.0, 0.0))
        push!(starts, (-1.0, 0.0))
        push!(starts, (1.0, 0.0))
        push!(starts, (0.0, -1.0))
        push!(starts, (0.0, 1.0))
    end

    for a in axes(local_nodes, 2)
        u = Float64(local_nodes[1, a])
        v = dim == 2 ? Float64(local_nodes[2, a]) : 0.0
        push!(starts, (u, v))
    end

    unique!(starts)
    return starts
end

function _contact_projection_workspace(elements::Vector{_ContactElement})
    max_nodes = maximum(e.basis.num_nodes for e in elements)
    max_degree = maximum(max(e.basis.max_u, e.basis.max_v) for e in elements)

    return _ContactProjectionWorkspace(
        zeros(Float64, max_nodes),
        zeros(Float64, 2, max_nodes),
        zeros(Float64, max_nodes),
        ones(Float64, max_degree + 1),
        ones(Float64, max_degree + 1),
        zeros(Float64, 3),
        zeros(Float64, 3, 2),
        zeros(Float64, 3),
        zeros(Float64, 3)
    )
end

# -----------------------------------------------------------------------------
# AABB tree
# -----------------------------------------------------------------------------

"""
    _contact_build_aabb_tree(elements, indices; leaf_size=2)

Build a binary AABB tree over master elements.
"""
function _contact_build_aabb_tree(
    elements::Vector{_ContactElement},
    indices::Vector{Int};
    leaf_size::Int=2
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
    xs,
    workspace::_ContactProjectionWorkspace;
    excluded::Union{Nothing,Set{Int}},
    normal_sign::Float64,
    model_dim::Int,
    projection_tol::Float64,
    projection_maxiter::Int,
    previous_element_index::Int=0,
    previous_local_coordinate::Union{Nothing,AbstractVector}=nothing
)

    best = _ContactSearchResult(0, 0.0, 0.0, Inf)

    # Store the previous projection separately for diagnostics.
    previous_u = 0.0
    previous_v = 0.0
    previous_d2 = Inf

    # Warm start from the previous contact association. This does not replace
    # the global search: it only supplies a finite initial upper bound for the
    # branch-and-bound traversal. Hence a different master element is still
    # selected whenever it gives a smaller distance.
    if previous_element_index != 0 &&
       previous_local_coordinate !== nothing

        element = elements[previous_element_index]

        u0 = Float64(previous_local_coordinate[1])
        v0 = element.dim == 1 ||
             length(previous_local_coordinate) < 2 ?
             0.0 : Float64(previous_local_coordinate[2])

        u, v, d2 = _contact_project_from_start!(
            element,
            xs,
            u0,
            v0,
            workspace;
            tol=projection_tol,
            maxiter=projection_maxiter
        )

        if isfinite(d2)
            best.element_index = previous_element_index
            best.u = u
            best.v = v
            best.distance2 = d2

            # Save the warm-start result before the global search can overwrite it.
            previous_u = u
            previous_v = v
            previous_d2 = d2
        end
    end

    _contact_search_aabb!(
        tree,
        elements,
        xs,
        excluded,
        best,
        workspace,
        projection_tol,
        projection_maxiter
    )

    # Temporary diagnostic for master-element switching.
    if previous_element_index != 0 &&
       isfinite(previous_d2) &&
       best.element_index != 0 &&
       best.element_index != previous_element_index

        d_old = sqrt(previous_d2)
        d_new = sqrt(best.distance2)

        println(
            "master switch: ",
            previous_element_index, " -> ", best.element_index,
            ", d_old = ", d_old,
            ", d_new = ", d_new,
            ", Δd = ", d_old - d_new,
            ", ξ_old = (", previous_u, ", ", previous_v, ")",
            ", ξ_new = (", best.u, ", ", best.v, ")"
        )
    end

    best.element_index == 0 && return nothing

    return _contact_finalize_projection(
        elements[best.element_index],
        xs,
        best.u,
        best.v,
        best.distance2,
        best.element_index,
        workspace;
        normal_sign=normal_sign,
        model_dim=model_dim
    )
end

function _contact_search_aabb!(
    node::_ContactAABBNode,
    elements::Vector{_ContactElement},
    xs,
    excluded::Union{Nothing,Set{Int}},
    best::_ContactSearchResult,
    workspace::_ContactProjectionWorkspace,
    projection_tol::Float64,
    projection_maxiter::Int
    )

    _contact_point_aabb_distance2(xs, node.bmin, node.bmax) > best.distance2 &&
        return nothing

    if node.left === nothing && node.right === nothing
        @inbounds for idx in node.elements
            excluded !== nothing && idx in excluded && continue

            element = elements[idx]

            # Exact lower bound for this individual element AABB.
            # If it is already farther away than the current best projection,
            # this element cannot improve the solution.
            _contact_point_aabb_distance2(
                xs,
                element.bmin,
                element.bmax
            ) > best.distance2 && continue

            d2, u, v = _contact_project_element_distance!(
                element,
                xs,
                workspace;
                tol=projection_tol,
                maxiter=projection_maxiter
            )

            if d2 < best.distance2
                best.element_index = idx
                best.u = u
                best.v = v
                best.distance2 = d2
            end
        end

        return nothing
    end

    left = node.left
    right = node.right

    if left === nothing
        _contact_search_aabb!(
            right, elements, xs, excluded, best, workspace,
            projection_tol, projection_maxiter
        )
        return nothing
    elseif right === nothing
        _contact_search_aabb!(
            left, elements, xs, excluded, best, workspace,
            projection_tol, projection_maxiter
        )
        return nothing
    end

    dl = _contact_point_aabb_distance2(xs, left.bmin, left.bmax)
    dr = _contact_point_aabb_distance2(xs, right.bmin, right.bmax)

    if dl <= dr
        _contact_search_aabb!(
            left, elements, xs, excluded, best, workspace,
            projection_tol, projection_maxiter
        )
        _contact_search_aabb!(
            right, elements, xs, excluded, best, workspace,
            projection_tol, projection_maxiter
        )
    else
        _contact_search_aabb!(
            right, elements, xs, excluded, best, workspace,
            projection_tol, projection_maxiter
        )
        _contact_search_aabb!(
            left, elements, xs, excluded, best, workspace,
            projection_tol, projection_maxiter
        )
    end

    return nothing
end

# -----------------------------------------------------------------------------
# Cached parametric element projection
# -----------------------------------------------------------------------------

"""
    _contact_project_element_distance!(element, xs, workspace; tol, maxiter)

Return the squared closest distance and the corresponding reference coordinate.
The routine is allocation-free in the iterative hot path. It is fully generic
with respect to the polynomial order because the element-specific Gmsh
Lagrange basis is represented by `_ContactBasisCache`.
"""
function _contact_project_element_distance!(
    element::_ContactElement,
    xs,
    workspace::_ContactProjectionWorkspace;
    tol::Float64,
    maxiter::Int
    )

    best_d2 = Inf
    best_u = 0.0
    best_v = 0.0

    @inbounds for start in element.basis.starts
        u, v, d2 = _contact_project_from_start!(
            element,
            xs,
            start[1],
            start[2],
            workspace;
            tol=tol,
            maxiter=maxiter
        )

        if d2 < best_d2
            best_d2 = d2
            best_u = u
            best_v = v
        end
    end

    return best_d2, best_u, best_v
end

function _contact_project_from_start!(
    element::_ContactElement,
    xs,
    u0::Float64,
    v0::Float64,
    workspace::_ContactProjectionWorkspace;
    tol::Float64,
    maxiter::Int
    )

    u, v = _contact_project_reference_uv(element.basis.refkind, u0, v0)
    _contact_geometry!(workspace, element, u, v)

    r1 = workspace.x[1] - xs[1]
    r2 = workspace.x[2] - xs[2]
    r3 = workspace.x[3] - xs[3]
    d2 = r1 * r1 + r2 * r2 + r3 * r3

    dx = element.bmax[1] - element.bmin[1]
    dy = element.bmax[2] - element.bmin[2]
    dz = element.bmax[3] - element.bmin[3]
    scale = max(sqrt(dx * dx + dy * dy + dz * dz), 1.0)
    gtol = tol * scale^2
    gtol2 = gtol * gtol
    ξtol2 = tol * tol

    for _ in 1:maxiter
        j11 = workspace.J[1, 1]
        j21 = workspace.J[2, 1]
        j31 = workspace.J[3, 1]

        g1 = j11 * r1 + j21 * r2 + j31 * r3
        h11 = j11 * j11 + j21 * j21 + j31 * j31

        δu = 0.0
        δv = 0.0

        if element.dim == 1
            g1 * g1 <= gtol2 && break

            reg = max(h11, 1.0) * 1e-14
            denom = h11 + reg
            abs(denom) > eps(Float64) || break
            δu = -g1 / denom
        else
            j12 = workspace.J[1, 2]
            j22 = workspace.J[2, 2]
            j32 = workspace.J[3, 2]

            g2 = j12 * r1 + j22 * r2 + j32 * r3
            g1 * g1 + g2 * g2 <= gtol2 && break

            h12 = j11 * j12 + j21 * j22 + j31 * j32
            h22 = j12 * j12 + j22 * j22 + j32 * j32
            hnorm = max(abs(h11) + abs(h12), abs(h12) + abs(h22), 1.0)
            reg = hnorm * 1e-14

            a = h11 + reg
            b = h12
            d = h22 + reg
            detH = a * d - b * b

            abs(detH) > eps(Float64) * hnorm * hnorm || break

            δu = (-d * g1 + b * g2) / detH
            δv = ( b * g1 - a * g2) / detH
        end

        isfinite(δu) && isfinite(δv) || break

        accepted = false
        α = 1.0

        for _ in 1:12
            utrial, vtrial = _contact_project_reference_uv(
                element.basis.refkind,
                u + α * δu,
                v + α * δv
            )

            du = utrial - u
            dv = vtrial - v
            step2 = du * du + dv * dv

            if step2 <= ξtol2
                u = utrial
                v = vtrial
                accepted = true
                break
            end

            _contact_position_trial!(workspace, element, utrial, vtrial)
            rt1 = workspace.xtrial[1] - xs[1]
            rt2 = workspace.xtrial[2] - xs[2]
            rt3 = workspace.xtrial[3] - xs[3]
            d2trial = rt1 * rt1 + rt2 * rt2 + rt3 * rt3

            if d2trial <= d2
                u = utrial
                v = vtrial
                d2 = d2trial
                accepted = true
                break
            end

            α *= 0.5
        end

        accepted || break

        # Refresh N, x and J only once after the accepted line-search step.
        _contact_geometry!(workspace, element, u, v)
        r1 = workspace.x[1] - xs[1]
        r2 = workspace.x[2] - xs[2]
        r3 = workspace.x[3] - xs[3]
        d2 = r1 * r1 + r2 * r2 + r3 * r3

        α * α * (δu * δu + δv * δv) <= ξtol2 && break
    end

    return u, v, d2
end

function _contact_finalize_projection(
    element::_ContactElement,
    xs,
    u::Float64,
    v::Float64,
    distance2::Float64,
    element_index::Int,
    workspace::_ContactProjectionWorkspace;
    normal_sign::Float64,
    model_dim::Int
    )

    _contact_geometry!(workspace, element, u, v)

    nn = element.basis.num_nodes
    N = copy(@view workspace.N[1:nn])
    x = copy(workspace.x)
    J = copy(workspace.J)

    separation = [
        xs[1] - x[1],
        xs[2] - x[2],
        xs[3] - x[3]
    ]

    normal = _contact_normal(
        element,
        J,
        separation,
        model_dim,
        normal_sign
    )

    tangent1, tangent2 = _contact_tangent_basis(
        element,
        J,
        normal,
        model_dim
    )

    gap = separation[1] * normal[1] +
          separation[2] * normal[2] +
          separation[3] * normal[3]

    ξ = element.dim == 1 ? [u] : [u, v]

    return _ContactProjection(
        element_index,
        element.tag,
        ξ,
        x,
        N,
        normal,
        tangent1,
        tangent2,
        gap,
        distance2
    )
end

function _contact_geometry!(
    workspace::_ContactProjectionWorkspace,
    element::_ContactElement,
    u::Float64,
    v::Float64
    )

    _contact_eval_basis!(workspace, element.basis, u, v, true, false)

    nn = element.basis.num_nodes
    x1 = 0.0
    x2 = 0.0
    x3 = 0.0
    j11 = 0.0
    j21 = 0.0
    j31 = 0.0
    j12 = 0.0
    j22 = 0.0
    j32 = 0.0

    @inbounds for a in 1:nn
        Na = workspace.N[a]
        dNu = workspace.dN[1, a]

        X1 = element.coords[1, a]
        X2 = element.coords[2, a]
        X3 = element.coords[3, a]

        x1 += X1 * Na
        x2 += X2 * Na
        x3 += X3 * Na

        j11 += X1 * dNu
        j21 += X2 * dNu
        j31 += X3 * dNu

        if element.dim == 2
            dNv = workspace.dN[2, a]
            j12 += X1 * dNv
            j22 += X2 * dNv
            j32 += X3 * dNv
        end
    end

    workspace.x[1] = x1
    workspace.x[2] = x2
    workspace.x[3] = x3

    workspace.J[1, 1] = j11
    workspace.J[2, 1] = j21
    workspace.J[3, 1] = j31
    workspace.J[1, 2] = j12
    workspace.J[2, 2] = j22
    workspace.J[3, 2] = j32

    return nothing
end

function _contact_position_trial!(
    workspace::_ContactProjectionWorkspace,
    element::_ContactElement,
    u::Float64,
    v::Float64
    )

    _contact_eval_basis!(workspace, element.basis, u, v, false, true)

    nn = element.basis.num_nodes
    x1 = 0.0
    x2 = 0.0
    x3 = 0.0

    @inbounds for a in 1:nn
        Na = workspace.Ntrial[a]
        x1 += element.coords[1, a] * Na
        x2 += element.coords[2, a] * Na
        x3 += element.coords[3, a] * Na
    end

    workspace.xtrial[1] = x1
    workspace.xtrial[2] = x2
    workspace.xtrial[3] = x3
    return nothing
end

function _contact_eval_basis!(
    workspace::_ContactProjectionWorkspace,
    cache::_ContactBasisCache,
    u::Float64,
    v::Float64,
    derivatives::Bool,
    trial::Bool
    )

    if cache.use_gmsh
        _contact_eval_basis_gmsh!(workspace, cache, u, v, derivatives, trial)
        return nothing
    end

    N = trial ? workspace.Ntrial : workspace.N
    nn = cache.num_nodes

    @inbounds for a in 1:nn
        N[a] = 0.0
        if derivatives
            workspace.dN[1, a] = 0.0
            workspace.dN[2, a] = 0.0
        end
    end

    workspace.upow[1] = 1.0
    @inbounds for a in 1:cache.max_u
        workspace.upow[a + 1] = workspace.upow[a] * u
    end

    workspace.vpow[1] = 1.0
    @inbounds for b in 1:cache.max_v
        workspace.vpow[b + 1] = workspace.vpow[b] * v
    end

    coeff = cache.coeff
    exp_u = cache.exp_u
    exp_v = cache.exp_v

    if derivatives
        @inbounds for k in eachindex(exp_u)
            aexp = exp_u[k]
            bexp = exp_v[k]
            φ = workspace.upow[aexp + 1] * workspace.vpow[bexp + 1]
            dφdu = aexp == 0 ? 0.0 :
                aexp * workspace.upow[aexp] * workspace.vpow[bexp + 1]
            dφdv = bexp == 0 ? 0.0 :
                bexp * workspace.upow[aexp + 1] * workspace.vpow[bexp]

            for a in 1:nn
                c = coeff[k, a]
                N[a] += φ * c
                workspace.dN[1, a] += dφdu * c
                workspace.dN[2, a] += dφdv * c
            end
        end
    else
        @inbounds for k in eachindex(exp_u)
            φ = workspace.upow[exp_u[k] + 1] * workspace.vpow[exp_v[k] + 1]
            for a in 1:nn
                N[a] += φ * coeff[k, a]
            end
        end
    end

    return nothing
end

function _contact_eval_basis_gmsh!(
    workspace::_ContactProjectionWorkspace,
    cache::_ContactBasisCache,
    u::Float64,
    v::Float64,
    derivatives::Bool,
    trial::Bool
    )

    workspace.localcoord[1] = u
    workspace.localcoord[2] = v
    workspace.localcoord[3] = 0.0

    N = trial ? workspace.Ntrial : workspace.N
    _, basis, _ = gmsh.model.mesh.getBasisFunctions(
        cache.etype,
        workspace.localcoord,
        "Lagrange"
    )

    length(basis) == cache.num_nodes ||
        error(
            "contact: unexpected number of Lagrange basis functions for " *
            "element type $(cache.etype)."
        )

    @inbounds for a in 1:cache.num_nodes
        N[a] = Float64(basis[a])
    end

    if derivatives
        _, grad, _ = gmsh.model.mesh.getBasisFunctions(
            cache.etype,
            workspace.localcoord,
            "GradLagrange"
        )

        length(grad) == 3 * cache.num_nodes ||
            error(
                "contact: unexpected number of Lagrange basis gradients for " *
                "element type $(cache.etype)."
            )

        @inbounds for a in 1:cache.num_nodes
            workspace.dN[1, a] = Float64(grad[3a - 2])
            workspace.dN[2, a] = cache.dim == 2 ? Float64(grad[3a - 1]) : 0.0
        end
    end

    return nothing
end

@inline function _contact_project_reference_uv(
    refkind::UInt8,
    u::Float64,
    v::Float64
    )

    if refkind == _CONTACT_REF_LINE
        return clamp(u, -1.0, 1.0), 0.0
    elseif refkind == _CONTACT_REF_QUAD
        return clamp(u, -1.0, 1.0), clamp(v, -1.0, 1.0)
    elseif refkind == _CONTACT_REF_TRI
        return _contact_project_triangle_reference_uv(u, v)
    end

    error("contact: unknown reference-element kind.")
end

@inline function _contact_project_triangle_reference_uv(u::Float64, v::Float64)
    if u >= 0.0 && v >= 0.0 && u + v <= 1.0
        return u, v
    end

    # Edge u = 0.
    u1 = 0.0
    v1 = clamp(v, 0.0, 1.0)
    d1 = (u - u1)^2 + (v - v1)^2

    # Edge v = 0.
    u2 = clamp(u, 0.0, 1.0)
    v2 = 0.0
    d2 = (u - u2)^2 + (v - v2)^2

    # Edge u + v = 1.
    t = clamp((u - v + 1.0) / 2.0, 0.0, 1.0)
    u3 = t
    v3 = 1.0 - t
    d3 = (u - u3)^2 + (v - v3)^2

    if d1 <= d2 && d1 <= d3
        return u1, v1
    elseif d2 <= d3
        return u2, v2
    else
        return u3, v3
    end
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
    U::Problem,
    slave_nodes::Vector{Int}
    )

    multiplier === nothing && return Int[]

    pdim = U.pdim
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
