###############################################################################
#                                                                             #
#                              Contact                                        #
#                                                                             #
###############################################################################

export Contact, ContactSet, ContactGap, ContactStiffness, contact, updateContact!
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
    coeff::Matrix{Float64}
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
    coords::Matrix{Float64}
    local_nodes::Matrix{Float64}
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

struct _ContactReferenceVertex
    node::Int
    local_index::Int
    ξ::NTuple{2,Float64}
end

struct _ContactReferenceEdge
    key::Tuple{Int,Int}
    node_a::Int
    node_b::Int
    ξa::NTuple{2,Float64}
    ξb::NTuple{2,Float64}
end

struct _ContactMasterTopology
    element_vertices::Vector{Vector{_ContactReferenceVertex}}
    element_edges::Vector{Vector{_ContactReferenceEdge}}
    vertex_elements::Dict{Int,Vector{Tuple{Int,_ContactReferenceVertex}}}
    edge_elements::Dict{Tuple{Int,Int},Vector{Tuple{Int,_ContactReferenceEdge}}}
end



# -----------------------------------------------------------------------------
# Gauss-point warm-start state
# -----------------------------------------------------------------------------

mutable struct _ContactGPWarmStart
    master_element_tag::Int
    local_coordinate::Vector{Float64}
end

_ContactGPWarmStart() = _ContactGPWarmStart(0, [0.0, 0.0])

# -----------------------------------------------------------------------------
# Contact geometry object
# -----------------------------------------------------------------------------

"""
    Contact

Geometry/search state of one slave-master contact pair.

`Contact` intentionally stores geometry and search data only. Contact weak-form
operators are created by [`ContactGap`](@ref), while penalty or multiplier
matrices are assembled by the ordinary `∫` syntax.

The current configuration is `x = X + displacement`. `updateContact!` refreshes
closest-point projections, normals, tangents, nodal gap postprocessing fields,
and all search caches needed by `ContactGap` integration.
"""
mutable struct Contact
    master::String
    slave::String
    U::Problem
    displacement::VectorField
    step::Int

    slave_nodes::Vector{Int}
    projections::Vector{_ContactProjection}
    master_element_tags::Vector{Int}
    master_local_coordinates::Vector{Vector{Float64}}
    master_points::Matrix{Float64}

    gap::ScalarField
    gap_values::Vector{Float64}
    n::VectorField
    t1::VectorField
    t2::Union{Nothing,VectorField}
    active::BitVector

    nodecoords::Matrix{Float64}
    slave_elements::Vector{_ContactElement}
    master_elements::Vector{_ContactElement}
    tree::_ContactAABBNode
    topology::_ContactMasterTopology
    excluded::Dict{Int,Set{Int}}
    master_tag_to_index::Dict{Int,Int}

    # Warm-start data for Gauss-point closest-point projections. The key is
    # `(gauss_rule, slave_element_tag)`. The cache survives `updateContact!`
    # calls because master element tags and slave element tags are stable on a
    # fixed mesh; tags are remapped to current master-element indices before
    # every search.
    gp_warm_start::Dict{Any,Vector{_ContactGPWarmStart}}

    options::NamedTuple
end

function Base.show(io::IO, c::Contact)
    nc = length(c.slave_nodes)
    na = count(c.active)
    print(
        io,
        "Contact(\"$(c.slave)\" -> \"$(c.master)\", " *
        "$nc candidate nodes, $na active)"
    )
end

"""
    ContactSet(contacts...)

Lightweight container for several independent slave-master contact pairs.
Matrices assembled from individual contact pairs can be added directly, e.g.

    Kc = ∫(ContactGap(C1) ⋅ D1 ⋅ ContactGap(C1)) +
         ∫(ContactGap(C2) ⋅ D2 ⋅ ContactGap(C2))
"""
struct ContactSet
    contacts::Vector{Contact}

    function ContactSet(contacts::Vector{Contact})
        isempty(contacts) && error("ContactSet: at least one contact pair is required.")
        U = contacts[1].U
        seen = Set{Tuple{String,String}}()
        for (i, c) in enumerate(contacts)
            c.U === U || error("ContactSet: contact pair $i uses a different displacement Problem.")
            pair = (c.slave, c.master)
            pair in seen && error("ContactSet: duplicate contact pair slave='$(c.slave)', master='$(c.master)'.")
            push!(seen, pair)
        end
        new(contacts)
    end
end

ContactSet(c::Contact, cs::Contact...) = ContactSet(Contact[c, cs...])
ContactSet(cs::AbstractVector{<:Contact}) = ContactSet(Contact[cs...])
Base.length(cs::ContactSet) = length(cs.contacts)
Base.getindex(cs::ContactSet, i::Int) = cs.contacts[i]
Base.iterate(cs::ContactSet, state...) = iterate(cs.contacts, state...)
Base.firstindex(cs::ContactSet) = firstindex(cs.contacts)
Base.lastindex(cs::ContactSet) = lastindex(cs.contacts)
Base.eltype(::Type{ContactSet}) = Contact

function Base.show(io::IO, cs::ContactSet)
    nc = sum(length(c.slave_nodes) for c in cs.contacts)
    na = sum(count(c.active) for c in cs.contacts)
    print(io, "ContactSet($(length(cs)) pairs, $nc candidate nodes, $na active)")
end


# -----------------------------------------------------------------------------
# Local contact constitutive coefficient
# -----------------------------------------------------------------------------

struct _ContactStiffnessMatrix{Tn,Tt} <: AbstractMatrix{Any}
    dim::Int
    cn::Tn
    ct::Tt
end

Base.size(D::_ContactStiffnessMatrix) = (D.dim, D.dim)
Base.IndexStyle(::Type{<:_ContactStiffnessMatrix}) = IndexCartesian()

function Base.getindex(D::_ContactStiffnessMatrix, i::Int, j::Int)
    1 <= i <= D.dim || throw(BoundsError(D, (i,j)))
    1 <= j <= D.dim || throw(BoundsError(D, (i,j)))
    i == j || return 0.0
    return i == 1 ? D.cn : D.ct
end

"""
    ContactStiffness(C::Contact, cn; ct=0.0)

Return the local contact constitutive matrix used between two full
`ContactGap(C; components=:all)` operators.

In 2D it represents

    [cn   0
      0  ct]

and in 3D

    [cn   0   0
      0  ct   0
      0   0  ct]

`cn` and `ct` may be numbers or nodal `ScalarField`s. The object behaves as an
`AbstractMatrix`, so it can be inserted directly into the standard LLFEM
matrix-chain syntax without exposing a Julia matrix literal in user code.
"""
function ContactStiffness(c::Contact, cn; ct=0.0)
    return _ContactStiffnessMatrix(c.U.pdim, cn, ct)
end

# -----------------------------------------------------------------------------
# Contact weak-form operator
# -----------------------------------------------------------------------------

"""
    ContactGap(C::Contact; components=:normal, active=:current)

Create a contact kinematic operator for the LLFEM weak-form DSL.

At each slave-side integration point the operator is assembled directly from
the slave interpolation, the closest-point master interpolation and the current
local contact basis. No nodal contact matrix is interpolated.

For `components=:normal`, the operator output is the scalar normal relative
position. For `components=:all`, the output ordering is `(n,t)` in 2D and
`(n,t1,t2)` in 3D.

`active=:current` integrates only Gauss points whose current normal gap satisfies
`gap <= C.options.activation_tol`. `active=:all` integrates the entire slave
candidate manifold.

Examples
--------
Frictionless penalty:

    Gn = ContactGap(C; components=:normal)
    Kc = ∫(Gn ⋅ cn ⋅ Gn)

Normal and tangential penalty:

    G = ContactGap(C; components=:all)
    D = Diagonal([cn, ct, ct])
    Kc = ∫(G ⋅ D ⋅ G)

Mixed Lagrange-multiplier coupling:

    Gn = ContactGap(C; components=:normal)
    B = ∫(Λ ⋅ Gn)

The returned applied operator is also callable on a displacement `VectorField`
for postprocessing:

    gap = ContactGap(C; components=:normal)(u)
    d   = ContactGap(C; components=:all)(u)

The field evaluation adds the reference nodal coordinates internally, i.e. it
evaluates the relative current position corresponding to `X + u` using the
current frozen contact projection and local basis.
"""
struct ContactGapOp <: AbstractOp
    contact::Contact
    components::Symbol
    active::Symbol
end

function ContactGap(
    c::Contact;
    components::Symbol=:normal,
    active::Symbol=:current
    )

    components in (:all, :normal) ||
        error("ContactGap: components must be :all or :normal.")
    active in (:current, :all) ||
        error("ContactGap: active must be :current or :all.")

    return OpApplied(c.U, ContactGapOp(c, components, active))
end

op_outdim(op::ContactGapOp, P::Problem) =
    op.components === :normal ? 1 : P.pdim

# Contact CSC patterns contain structural entries that can be numerically zero
# for a particular normal orientation or active set. Preserve them so the
# returned matrix can be reused through `csc_matrix=Kc.A`.
_preserve_csc_pattern(::ContactGapOp) = true

# -----------------------------------------------------------------------------
# Contact construction and update
# -----------------------------------------------------------------------------

"""
    contact(U::Problem; master, slave, displacement=zero_displacement,
            activation_tol=0.0, step=displacement.nsteps, kwargs...) -> Contact

Construct a slave-master contact geometry in the current configuration.

The object contains only geometry, closest-point and search data. Penalty and
Lagrange-multiplier operators are built with `ContactGap(C)` and assembled with
`∫`.
"""
function contact(
    U::Problem;
    master::String,
    slave::String,
    displacement::VectorField=_contact_zero_displacement(U),
    activation_tol::Real=0.0,
    step::Int=displacement.nsteps,
    normal_sign::Real=1.0,
    self_contact::Bool=(slave == master),
    self_exclusion_layers::Int=1,
    aabb_padding::Real=0.05,
    leaf_size::Int=2,
    projection_tol::Real=1e-10,
    projection_maxiter::Int=40,
    topology_tol::Real=1e-3,
    topology_angle::Real=45.0,
    LagrangeMultiplierField=nothing
    )

    # Kept only as a source-compatible keyword. Multiplier fields now enter the
    # weak form directly, e.g. `∫(Λ ⋅ ContactGap(C; components=:normal))`.
    options = (
        activation_tol=Float64(activation_tol),
        normal_sign=Float64(normal_sign),
        self_contact=self_contact,
        self_exclusion_layers=self_exclusion_layers,
        aabb_padding=Float64(aabb_padding),
        leaf_size=leaf_size,
        projection_tol=Float64(projection_tol),
        projection_maxiter=projection_maxiter,
        topology_tol=Float64(topology_tol),
        topology_angle=Float64(topology_angle)
    )

    data = _contact_build_geometry(
        U,
        displacement,
        slave,
        master;
        step=step,
        options...
    )

    return Contact(
        master,
        slave,
        U,
        displacement,
        step,
        data.slave_nodes,
        data.projections,
        data.master_element_tags,
        data.master_local_coordinates,
        data.master_points,
        data.gap,
        data.gap_values,
        data.n,
        data.t1,
        data.t2,
        data.active,
        data.nodecoords,
        data.slave_elements,
        data.master_elements,
        data.tree,
        data.topology,
        data.excluded,
        data.master_tag_to_index,
        Dict{Any,Vector{_ContactGPWarmStart}}(),
        options
    )
end

function contact(
    displacement::VectorField;
    master::String,
    slave::String,
    U::Problem=displacement.model,
    kwargs...
    )
    return contact(U; master=master, slave=slave, displacement=displacement, kwargs...)
end

contact(displacement::VectorField, slave::String, master::String; kwargs...) =
    contact(displacement; slave=slave, master=master, kwargs...)

contact(U::Problem, slave::String, master::String; kwargs...) =
    contact(U; slave=slave, master=master, kwargs...)

"""
    updateContact!(C::Contact, displacement::VectorField; step=displacement.nsteps)

Refresh current closest-point geometry and search caches while preserving the
`Contact` object identity used by existing `ContactGap(C)` operators.
"""
function updateContact!(
    c::Contact,
    displacement::VectorField;
    step::Int=displacement.nsteps
    )

    data = _contact_build_geometry(
        c.U,
        displacement,
        c.slave,
        c.master;
        step=step,
        previous_slave_nodes=c.slave_nodes,
        previous_master_element_tags=c.master_element_tags,
        previous_master_local_coordinates=c.master_local_coordinates,
        c.options...
    )

    c.displacement = displacement
    c.step = step
    c.slave_nodes = data.slave_nodes
    c.projections = data.projections
    c.master_element_tags = data.master_element_tags
    c.master_local_coordinates = data.master_local_coordinates
    c.master_points = data.master_points
    c.gap = data.gap
    c.gap_values = data.gap_values
    c.n = data.n
    c.t1 = data.t1
    c.t2 = data.t2
    c.active = data.active
    c.nodecoords = data.nodecoords
    c.slave_elements = data.slave_elements
    c.master_elements = data.master_elements
    c.tree = data.tree
    c.topology = data.topology
    c.excluded = data.excluded
    c.master_tag_to_index = data.master_tag_to_index

    return c
end

function updateContact!(
    cs::ContactSet,
    displacement::VectorField;
    step::Int=displacement.nsteps
    )
    for c in cs.contacts
        updateContact!(c, displacement; step=step)
    end
    return cs
end

# -----------------------------------------------------------------------------
# Geometry builder used by contact/updateContact!
# -----------------------------------------------------------------------------

function _contact_build_geometry(
    U::Problem,
    displacement::VectorField,
    slave::String,
    master::String;
    step::Int,
    activation_tol::Float64,
    normal_sign::Float64,
    self_contact::Bool,
    self_exclusion_layers::Int,
    aabb_padding::Float64,
    leaf_size::Int,
    projection_tol::Float64,
    projection_maxiter::Int,
    topology_tol::Float64,
    topology_angle::Float64,
    previous_slave_nodes::Union{Nothing,Vector{Int}}=nothing,
    previous_master_element_tags::Union{Nothing,Vector{Int}}=nothing,
    previous_master_local_coordinates::Union{Nothing,Vector{Vector{Float64}}}=nothing
    )

    _contact_check_models(U, displacement, nothing)

    self_exclusion_layers >= 0 || error("contact: self_exclusion_layers must be non-negative.")
    aabb_padding >= 0 || error("contact: aabb_padding must be non-negative.")
    leaf_size >= 1 || error("contact: leaf_size must be at least one.")
    projection_tol > 0 || error("contact: projection_tol must be positive.")
    projection_maxiter >= 1 || error("contact: projection_maxiter must be at least one.")
    topology_tol >= 0 || error("contact: topology_tol must be non-negative.")
    isfinite(topology_angle) && 0.0 <= topology_angle <= 180.0 ||
        error("contact: topology_angle must be between 0 and 180 degrees.")
    isfinite(normal_sign) && abs(abs(normal_sign) - 1.0) <= 10 * eps(Float64) ||
        error("contact: normal_sign must be either +1 or -1.")
    isfinite(activation_tol) || error("contact: activation_tol must be finite.")

    gmsh.model.setCurrent(U.name)
    nodecoords = _contact_deformed_coordinates(U, displacement; step=step)

    slave_elements, slave_dim =
        _contact_group_elements(U, slave, nodecoords; aabb_padding=0.0)
    master_elements, master_dim =
        _contact_group_elements(U, master, nodecoords; aabb_padding=aabb_padding)

    slave_dim in (1, 2) || error("contact: slave physical group '$slave' must be a curve or surface.")
    master_dim in (1, 2) || error("contact: master physical group '$master' must be a curve or surface.")
    isempty(slave_elements) && error("contact: no finite elements were found in slave group '$slave'.")
    isempty(master_elements) && error("contact: no finite elements were found in master group '$master'.")

    _contact_check_master_element_types(master_elements)
    topology = _contact_master_topology(master_elements)
    tree = _contact_build_aabb_tree(
        master_elements,
        collect(eachindex(master_elements));
        leaf_size=leaf_size
    )

    slave_nodes = sort!(unique!(vcat((e.node_tags for e in slave_elements)...)))

    excluded = self_contact ?
        _contact_self_exclusion_sets(slave_nodes, master_elements, self_exclusion_layers) :
        Dict{Int,Set{Int}}()

    workspace = _contact_projection_workspace(vcat(slave_elements, master_elements))
    projection_dict = Dict{Int,_ContactProjection}()
    sizehint!(projection_dict, length(slave_nodes))
    empty_exclusion = Set{Int}()

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

    master_tag_to_index = Dict{Int,Int}()
    sizehint!(master_tag_to_index, length(master_elements))
    @inbounds for (i, element) in enumerate(master_elements)
        master_tag_to_index[element.tag] = i
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
            topology=topology,
            topology_tol=topology_tol,
            topology_angle=topology_angle,
            previous_element_index=previous_element_index,
            previous_local_coordinate=previous_local_coordinate
        )

        p === nothing && error(
            "contact: no admissible master element was found for slave node $node. " *
            "In self-contact, try reducing self_exclusion_layers."
        )

        projection_dict[node] = p
    end

    gap = _contact_gap_field(U, slave_elements, projection_dict)
    normalVec = _contact_vector_field(U, slave_elements, projection_dict, :normal)
    tangent1 = _contact_vector_field(U, slave_elements, projection_dict, :tangent1)
    tangent2 = U.pdim == 3 ?
        _contact_vector_field(U, slave_elements, projection_dict, :tangent2) : nothing

    projections = Vector{_ContactProjection}(undef, length(slave_nodes))
    gap_values = Vector{Float64}(undef, length(slave_nodes))
    master_element_tags = Vector{Int}(undef, length(slave_nodes))
    master_local_coordinates = Vector{Vector{Float64}}(undef, length(slave_nodes))
    master_points = Matrix{Float64}(undef, 3, length(slave_nodes))

    @inbounds for (i, node) in enumerate(slave_nodes)
        p = projection_dict[node]
        projections[i] = p
        gap_values[i] = p.gap
        master_element_tags[i] = p.element_tag
        master_local_coordinates[i] = copy(p.ξ)
        master_points[:, i] .= p.x
    end

    active = BitVector(g <= activation_tol for g in gap_values)

    return (
        slave_nodes=slave_nodes,
        projections=projections,
        master_element_tags=master_element_tags,
        master_local_coordinates=master_local_coordinates,
        master_points=master_points,
        gap=gap,
        gap_values=gap_values,
        n=normalVec,
        t1=tangent1,
        t2=tangent2,
        active=active,
        nodecoords=nodecoords,
        slave_elements=slave_elements,
        master_elements=master_elements,
        tree=tree,
        topology=topology,
        excluded=excluded,
        master_tag_to_index=master_tag_to_index
    )
end
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

"""
    _contact_refresh_element_geometry!(
        elements,
        nodecoords;
        aabb_padding
    )

Update the current nodal coordinates and bounding boxes of contact elements
without rebuilding their connectivity, basis cache, or reference topology.
"""
function _contact_refresh_element_geometry!(
    elements::Vector{_ContactElement},
    nodecoords::Matrix{Float64};
    aabb_padding::Real
    )

    padding = Float64(aabb_padding)

    @inbounds for element in elements

        X = element.coords
        nodes = element.node_tags
        nn = length(nodes)

        bmin = element.bmin
        bmax = element.bmax

        bmin[1] = Inf
        bmin[2] = Inf
        bmin[3] = Inf

        bmax[1] = -Inf
        bmax[2] = -Inf
        bmax[3] = -Inf

        for a in 1:nn
            node = nodes[a]

            x = nodecoords[1, node]
            y = nodecoords[2, node]
            z = nodecoords[3, node]

            X[1, a] = x
            X[2, a] = y
            X[3, a] = z

            x < bmin[1] && (bmin[1] = x)
            y < bmin[2] && (bmin[2] = y)
            z < bmin[3] && (bmin[3] = z)

            x > bmax[1] && (bmax[1] = x)
            y > bmax[2] && (bmax[2] = y)
            z > bmax[3] && (bmax[3] = z)
        end

        dx = bmax[1] - bmin[1]
        dy = bmax[2] - bmin[2]
        dz = bmax[3] - bmin[3]

        diag = sqrt(dx^2 + dy^2 + dz^2)
        pad = padding * max(diag, eps(Float64))

        bmin[1] -= pad
        bmin[2] -= pad
        bmin[3] -= pad

        bmax[1] += pad
        bmax[2] += pad
        bmax[3] += pad
    end

    return nothing
end

"""
    _contact_update_geometry!(
        C,
        displacement;
        step=displacement.nsteps
    )

Update only the deformed contact geometry required by Gauss-point contact
integration.

Unlike `updateContact!`, this function does not recompute nodal closest-point
projections, nodal gap values, normals, tangents, or the nodal active set.
Gauss-point warm-start data are preserved.
"""
function _contact_update_geometry!(
    c::Contact,
    displacement::VectorField;
    step::Int=displacement.nsteps
    )

    _contact_check_models(c.U, displacement, nothing)

    # x = X + u
    nodecoords =
        _contact_deformed_coordinates(
            c.U,
            displacement;
            step=step
        )

    # Update element coordinates in place.
    _contact_refresh_element_geometry!(
        c.slave_elements,
        nodecoords;
        aabb_padding=0.0
    )

    _contact_refresh_element_geometry!(
        c.master_elements,
        nodecoords;
        aabb_padding=c.options.aabb_padding
    )

    # Only the spatial search tree depends on the current coordinates.
    # Connectivity/topology/exclusion sets remain unchanged.
    c.tree =
        _contact_build_aabb_tree(
            c.master_elements,
            collect(eachindex(c.master_elements));
            leaf_size=c.options.leaf_size
        )

    c.nodecoords = nodecoords
    c.displacement = displacement
    c.step = step

    return c
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
# Master topology and boundary stabilization
# -----------------------------------------------------------------------------

function _contact_reference_corners(refkind::UInt8)
    if refkind == _CONTACT_REF_LINE
        return [(-1.0, 0.0), (1.0, 0.0)]
    elseif refkind == _CONTACT_REF_TRI
        return [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    elseif refkind == _CONTACT_REF_QUAD
        return [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
    end
    error("contact: unknown reference-element kind.")
end

function _contact_master_topology(elements::Vector{_ContactElement})
    ne = length(elements)
    element_vertices = Vector{Vector{_ContactReferenceVertex}}(undef, ne)
    element_edges = Vector{Vector{_ContactReferenceEdge}}(undef, ne)

    vertex_elements =
        Dict{Int,Vector{Tuple{Int,_ContactReferenceVertex}}}()
    edge_elements =
        Dict{Tuple{Int,Int},Vector{Tuple{Int,_ContactReferenceEdge}}}()

    @inbounds for (ei, e) in enumerate(elements)
        corners = _contact_reference_corners(e.basis.refkind)
        vertices = _ContactReferenceVertex[]

        for ξ in corners
            best_a = 0
            best_d2 = Inf

            for a in axes(e.local_nodes, 2)
                u = Float64(e.local_nodes[1, a])
                v = e.dim == 2 ? Float64(e.local_nodes[2, a]) : 0.0
                d2 = (u - ξ[1])^2 + (v - ξ[2])^2

                if d2 < best_d2
                    best_d2 = d2
                    best_a = a
                end
            end

            best_a != 0 && best_d2 <= 1e-12 ||
                error(
                    "contact: failed to identify a corner node of master " *
                    "element $(e.tag)."
                )

            vertex = _ContactReferenceVertex(
                e.node_tags[best_a],
                best_a,
                ξ
            )
            push!(vertices, vertex)
            push!(
                get!(
                    vertex_elements,
                    vertex.node,
                    Tuple{Int,_ContactReferenceVertex}[]
                ),
                (ei, vertex)
            )
        end

        element_vertices[ei] = vertices

        edges = _ContactReferenceEdge[]
        if e.dim == 2
            pairs = if e.basis.refkind == _CONTACT_REF_TRI
                ((1, 2), (2, 3), (3, 1))
            elseif e.basis.refkind == _CONTACT_REF_QUAD
                ((1, 2), (2, 3), (3, 4), (4, 1))
            else
                ()
            end

            for (ia, ib) in pairs
                va = vertices[ia]
                vb = vertices[ib]
                key = va.node < vb.node ?
                    (va.node, vb.node) : (vb.node, va.node)

                edge = _ContactReferenceEdge(
                    key,
                    va.node,
                    vb.node,
                    va.ξ,
                    vb.ξ
                )
                push!(edges, edge)
                push!(
                    get!(
                        edge_elements,
                        key,
                        Tuple{Int,_ContactReferenceEdge}[]
                    ),
                    (ei, edge)
                )
            end
        end

        element_edges[ei] = edges
    end

    return _ContactMasterTopology(
        element_vertices,
        element_edges,
        vertex_elements,
        edge_elements
    )
end

@inline function _contact_reference_segment_projection(
    u::Float64,
    v::Float64,
    a::NTuple{2,Float64},
    b::NTuple{2,Float64}
    )

    du = b[1] - a[1]
    dv = b[2] - a[2]
    denom = du * du + dv * dv
    denom > 0.0 || return 0.0, (u - a[1])^2 + (v - a[2])^2

    t = clamp(
        ((u - a[1]) * du + (v - a[2]) * dv) / denom,
        0.0,
        1.0
    )
    ur = a[1] + t * du
    vr = a[2] + t * dv
    d2 = (u - ur)^2 + (v - vr)^2
    return t, d2
end

function _contact_reference_entity(
    projection::_ContactProjection,
    topology::_ContactMasterTopology,
    tol::Float64
    )

    ei = projection.element_index
    u = projection.ξ[1]
    v = length(projection.ξ) == 1 ? 0.0 : projection.ξ[2]
    tol2 = tol * tol

    # Vertices take precedence over edges.
    best_vertex = nothing
    best_d2 = Inf
    for vertex in topology.element_vertices[ei]
        d2 = (u - vertex.ξ[1])^2 + (v - vertex.ξ[2])^2
        if d2 < best_d2
            best_d2 = d2
            best_vertex = vertex
        end
    end

    if best_vertex !== nothing &&
       best_d2 <= tol2 &&
       length(get(topology.vertex_elements, best_vertex.node, Tuple{Int,_ContactReferenceVertex}[])) >= 2
        return :vertex, best_vertex, 0.0
    end

    best_edge = nothing
    best_t = 0.0
    best_d2 = Inf

    for edge in topology.element_edges[ei]
        t, d2 = _contact_reference_segment_projection(
            u, v, edge.ξa, edge.ξb
        )
        if d2 < best_d2
            best_d2 = d2
            best_t = t
            best_edge = edge
        end
    end

    if best_edge !== nothing &&
       best_d2 <= tol2 &&
       length(get(topology.edge_elements, best_edge.key, Tuple{Int,_ContactReferenceEdge}[])) >= 2
        return :edge, best_edge, best_t
    end

    return :interior, nothing, 0.0
end

@inline function _contact_jacobian_weight(
    element::_ContactElement,
    J::Matrix{Float64}
    )
    if element.dim == 2
        return norm(cross(@view(J[:, 1]), @view(J[:, 2])))
    end
    return norm(@view J[:, 1])
end

function _contact_average_vertex_normal!(
    topology::_ContactMasterTopology,
    elements::Vector{_ContactElement},
    vertex::_ContactReferenceVertex,
    xs,
    workspace::_ContactProjectionWorkspace;
    normal_sign::Float64,
    model_dim::Int,
    angle_cos::Float64
    )

    incident = get(
        topology.vertex_elements,
        vertex.node,
        Tuple{Int,_ContactReferenceVertex}[]
    )
    length(incident) >= 2 || return nothing

    nsum = zeros(Float64, 3)
    nref = nothing

    for (ei, local_vertex) in incident
        e = elements[ei]

        # Boundary stabilization is currently defined for planar curve contact
        # and 3D surface contact. Other dimensional combinations keep the
        # original manifold projection.
        if !((model_dim == 2 && e.dim == 1) ||
             (model_dim == 3 && e.dim == 2))
            return nothing
        end

        u = local_vertex.ξ[1]
        v = local_vertex.ξ[2]
        _contact_geometry!(workspace, e, u, v)

        separation = [
            xs[1] - workspace.x[1],
            xs[2] - workspace.x[2],
            xs[3] - workspace.x[3]
        ]

        ni = _contact_normal(
            e,
            workspace.J,
            separation,
            model_dim,
            normal_sign
        )

        if nref === nothing
            nref = copy(ni)
        else
            d = dot(ni, nref)
            if d < 0.0
                ni .*= -1.0
                d = -d
            end
            d >= angle_cos || return nothing
        end

        w = _contact_jacobian_weight(e, workspace.J)
        isfinite(w) && w > sqrt(eps(Float64)) || continue
        nsum .+= w .* ni
    end

    nrm = norm(nsum)
    nrm > sqrt(eps(Float64)) || return nothing
    nsum ./= nrm

    # Keep the sign consistent with the first incident normal.
    if nref !== nothing && dot(nsum, nref) < 0.0
        nsum .*= -1.0
    end

    return nsum
end

function _contact_average_edge_normal!(
    topology::_ContactMasterTopology,
    elements::Vector{_ContactElement},
    edge::_ContactReferenceEdge,
    canonical_t::Float64,
    xs,
    workspace::_ContactProjectionWorkspace;
    normal_sign::Float64,
    model_dim::Int,
    angle_cos::Float64
    )

    model_dim == 3 || return nothing

    incident = get(
        topology.edge_elements,
        edge.key,
        Tuple{Int,_ContactReferenceEdge}[]
    )
    length(incident) >= 2 || return nothing

    nsum = zeros(Float64, 3)
    nref = nothing

    for (ei, local_edge) in incident
        e = elements[ei]
        e.dim == 2 || return nothing

        t = local_edge.node_a == edge.key[1] ?
            canonical_t : 1.0 - canonical_t

        u = local_edge.ξa[1] +
            t * (local_edge.ξb[1] - local_edge.ξa[1])
        v = local_edge.ξa[2] +
            t * (local_edge.ξb[2] - local_edge.ξa[2])

        _contact_geometry!(workspace, e, u, v)

        separation = [
            xs[1] - workspace.x[1],
            xs[2] - workspace.x[2],
            xs[3] - workspace.x[3]
        ]

        ni = _contact_normal(
            e,
            workspace.J,
            separation,
            model_dim,
            normal_sign
        )

        if nref === nothing
            nref = copy(ni)
        else
            d = dot(ni, nref)
            if d < 0.0
                ni .*= -1.0
                d = -d
            end
            d >= angle_cos || return nothing
        end

        w = _contact_jacobian_weight(e, workspace.J)
        isfinite(w) && w > sqrt(eps(Float64)) || continue
        nsum .+= w .* ni
    end

    nrm = norm(nsum)
    nrm > sqrt(eps(Float64)) || return nothing
    nsum ./= nrm

    if nref !== nothing && dot(nsum, nref) < 0.0
        nsum .*= -1.0
    end

    return nsum
end

function _contact_stable_tangent_basis(
    n::Vector{Float64},
    model_dim::Int
    )
    if model_dim == 2
        t1 = [-n[2], n[1], 0.0]
        t1 ./= norm(t1)
        return t1, nothing
    end

    # Deterministic global-axis construction. This avoids inheriting a
    # tangential direction from whichever adjacent face happened to win the
    # closest-point search.
    axis = if abs(n[1]) <= abs(n[2]) && abs(n[1]) <= abs(n[3])
        [1.0, 0.0, 0.0]
    elseif abs(n[2]) <= abs(n[3])
        [0.0, 1.0, 0.0]
    else
        [0.0, 0.0, 1.0]
    end

    t1 = cross(axis, n)
    t1 ./= norm(t1)
    t2 = cross(n, t1)
    t2 ./= norm(t2)
    return Vector{Float64}(t1), Vector{Float64}(t2)
end

function _contact_project_edge!(
    element::_ContactElement,
    edge::_ContactReferenceEdge,
    xs,
    t0::Float64,
    workspace::_ContactProjectionWorkspace;
    tol::Float64,
    maxiter::Int=30
    )

    t = clamp(t0, 0.0, 1.0)
    du = edge.ξb[1] - edge.ξa[1]
    dv = edge.ξb[2] - edge.ξa[2]

    u = edge.ξa[1] + t * du
    v = edge.ξa[2] + t * dv
    _contact_geometry!(workspace, element, u, v)

    r = [
        workspace.x[1] - xs[1],
        workspace.x[2] - xs[2],
        workspace.x[3] - xs[3]
    ]
    d2 = dot(r, r)

    for _ in 1:maxiter
        dxdt = [
            workspace.J[1, 1] * du + workspace.J[1, 2] * dv,
            workspace.J[2, 1] * du + workspace.J[2, 2] * dv,
            workspace.J[3, 1] * du + workspace.J[3, 2] * dv
        ]

        h = dot(dxdt, dxdt)
        h > eps(Float64) || break

        δt = -dot(dxdt, r) / h
        abs(δt) <= tol && break

        accepted = false
        α = 1.0

        for _ in 1:12
            ttrial = clamp(t + α * δt, 0.0, 1.0)
            abs(ttrial - t) <= tol && begin
                t = ttrial
                accepted = true
                break
            end

            utrial = edge.ξa[1] + ttrial * du
            vtrial = edge.ξa[2] + ttrial * dv
            _contact_position_trial!(
                workspace,
                element,
                utrial,
                vtrial
            )

            rt = [
                workspace.xtrial[1] - xs[1],
                workspace.xtrial[2] - xs[2],
                workspace.xtrial[3] - xs[3]
            ]
            d2trial = dot(rt, rt)

            if d2trial <= d2
                t = ttrial
                d2 = d2trial
                accepted = true
                break
            end

            α *= 0.5
        end

        accepted || break

        u = edge.ξa[1] + t * du
        v = edge.ξa[2] + t * dv
        _contact_geometry!(workspace, element, u, v)

        r[1] = workspace.x[1] - xs[1]
        r[2] = workspace.x[2] - xs[2]
        r[3] = workspace.x[3] - xs[3]
        d2 = dot(r, r)
    end

    u = edge.ξa[1] + t * du
    v = edge.ξa[2] + t * dv
    return t, u, v, d2
end

function _contact_vertex_projection(
    original::_ContactProjection,
    elements::Vector{_ContactElement},
    topology::_ContactMasterTopology,
    vertex::_ContactReferenceVertex,
    xs,
    workspace::_ContactProjectionWorkspace;
    normal_sign::Float64,
    model_dim::Int,
    angle_cos::Float64
    )

    normal = _contact_average_vertex_normal!(
        topology,
        elements,
        vertex,
        xs,
        workspace;
        normal_sign=normal_sign,
        model_dim=model_dim,
        angle_cos=angle_cos
    )
    normal === nothing && return original

    e = elements[original.element_index]
    N = zeros(Float64, e.basis.num_nodes)
    N[vertex.local_index] = 1.0

    x = Vector{Float64}(e.coords[:, vertex.local_index])
    separation = [
        xs[1] - x[1],
        xs[2] - x[2],
        xs[3] - x[3]
    ]

    t1, t2 = _contact_stable_tangent_basis(normal, model_dim)
    gap = dot(separation, normal)
    d2 = dot(separation, separation)

    ξ = e.dim == 1 ?
        [vertex.ξ[1]] :
        [vertex.ξ[1], vertex.ξ[2]]

    return _ContactProjection(
        original.element_index,
        e.tag,
        ξ,
        x,
        N,
        normal,
        t1,
        t2,
        gap,
        d2
    )
end

function _contact_edge_projection(
    original::_ContactProjection,
    elements::Vector{_ContactElement},
    topology::_ContactMasterTopology,
    edge::_ContactReferenceEdge,
    t0::Float64,
    xs,
    workspace::_ContactProjectionWorkspace;
    projection_tol::Float64,
    normal_sign::Float64,
    model_dim::Int,
    angle_cos::Float64
    )

    e = elements[original.element_index]
    e.dim == 2 && model_dim == 3 || return original

    t, u, v, d2 = _contact_project_edge!(
        e,
        edge,
        xs,
        t0,
        workspace;
        tol=max(projection_tol, 1e-12)
    )

    # Convert to a canonical edge parameter independent of which incident face
    # supplied the projection.
    canonical_t =
        edge.node_a == edge.key[1] ? t : 1.0 - t

    normal = _contact_average_edge_normal!(
        topology,
        elements,
        edge,
        canonical_t,
        xs,
        workspace;
        normal_sign=normal_sign,
        model_dim=model_dim,
        angle_cos=angle_cos
    )
    normal === nothing && return original

    # Refresh the selected element after the averaging loop reused workspace.
    _contact_geometry!(workspace, e, u, v)

    N = copy(@view workspace.N[1:e.basis.num_nodes])
    x = copy(workspace.x)
    J = copy(workspace.J)

    du = edge.ξb[1] - edge.ξa[1]
    dv = edge.ξb[2] - edge.ξa[2]
    tangent = [
        J[1, 1] * du + J[1, 2] * dv,
        J[2, 1] * du + J[2, 2] * dv,
        J[3, 1] * du + J[3, 2] * dv
    ]

    if edge.node_a != edge.key[1]
        tangent .*= -1.0
    end

    tangent .-= dot(tangent, normal) .* normal
    if norm(tangent) <= sqrt(eps(Float64))
        t1, t2 = _contact_stable_tangent_basis(normal, model_dim)
    else
        tangent ./= norm(tangent)
        t1 = Vector{Float64}(tangent)
        t2 = cross(normal, t1)
        t2 ./= norm(t2)
        t2 = Vector{Float64}(t2)
    end

    separation = [
        xs[1] - x[1],
        xs[2] - x[2],
        xs[3] - x[3]
    ]
    gap = dot(separation, normal)

    return _ContactProjection(
        original.element_index,
        e.tag,
        [u, v],
        x,
        N,
        normal,
        t1,
        t2,
        gap,
        d2
    )
end

function _contact_stabilize_topological_projection(
    projection::_ContactProjection,
    elements::Vector{_ContactElement},
    topology::_ContactMasterTopology,
    xs,
    workspace::_ContactProjectionWorkspace;
    tol::Float64,
    projection_tol::Float64,
    normal_sign::Float64,
    model_dim::Int,
    angle_cos::Float64
    )

    kind, entity, t = _contact_reference_entity(
        projection,
        topology,
        tol
    )

    if kind === :vertex
        return _contact_vertex_projection(
            projection,
            elements,
            topology,
            entity,
            xs,
            workspace;
            normal_sign=normal_sign,
            model_dim=model_dim,
            angle_cos=angle_cos
        )
    elseif kind === :edge
        return _contact_edge_projection(
            projection,
            elements,
            topology,
            entity,
            t,
            xs,
            workspace;
            projection_tol=projection_tol,
            normal_sign=normal_sign,
            model_dim=model_dim,
            angle_cos=angle_cos
        )
    end

    return projection
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
    topology::_ContactMasterTopology,
    topology_tol::Float64,
    topology_angle::Float64,
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

    #=
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
    =#

    best.element_index == 0 && return nothing

    projection = _contact_finalize_projection(
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

    topology_tol <= 0 && return projection

    return _contact_stabilize_topological_projection(
        projection,
        elements,
        topology,
        xs,
        workspace;
        tol=topology_tol,
        projection_tol=projection_tol,
        normal_sign=normal_sign,
        model_dim=model_dim,
        angle_cos=cosd(topology_angle)
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
`C`, while the basis is still available for tangential kinematics and post-processing.
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

# -----------------------------------------------------------------------------
# ContactGap Gauss-point assembly
# -----------------------------------------------------------------------------

@inline function _contact_gauss_order(order::Int, gauss)
    if gauss === :reduced
        return max(1, 2order - 1)
    elseif gauss === :full
        return 2order + 1
    elseif gauss isa Int
        return max(1, 2order + 1 + gauss)
    else
        error("ContactGap integration: gauss must be :full, :reduced, or an integer offset.")
    end
end

function _contact_quadrature(element::_ContactElement, gauss)
    gorder = _contact_gauss_order(element.order, gauss)
    intPoints, intWeights =
        gmsh.model.mesh.getIntegrationPoints(
            element.etype,
            "Gauss" * string(gorder)
        )

    nip = length(intWeights)

    _, fun, _ = gmsh.model.mesh.getBasisFunctions(
        element.etype,
        intPoints,
        "Lagrange"
    )
    N = reshape(fun, :, nip)

    _, dfun, _ = gmsh.model.mesh.getBasisFunctions(
        element.etype,
        intPoints,
        "GradLagrange"
    )
    dN = reshape(dfun, :, nip)

    return (
        points=Float64.(intPoints),
        weights=Float64.(intWeights),
        N=N,
        dN=dN,
        nip=nip
    )
end

function _contact_slave_geometry_at_gp(
    element::_ContactElement,
    qdata,
    q::Int
    )

    nn = length(element.node_tags)
    N = @view qdata.N[1:nn, q]

    xs = zeros(Float64, 3)
    J = zeros(Float64, 3, 2)

    @inbounds for a in 1:nn
        Na = N[a]
        xs[1] += Na * element.coords[1, a]
        xs[2] += Na * element.coords[2, a]
        xs[3] += Na * element.coords[3, a]

        dNu = qdata.dN[3a - 2, q]
        J[1, 1] += dNu * element.coords[1, a]
        J[2, 1] += dNu * element.coords[2, a]
        J[3, 1] += dNu * element.coords[3, a]

        if element.dim == 2
            dNv = qdata.dN[3a - 1, q]
            J[1, 2] += dNv * element.coords[1, a]
            J[2, 2] += dNv * element.coords[2, a]
            J[3, 2] += dNv * element.coords[3, a]
        end
    end

    measure = if element.dim == 1
        sqrt(J[1,1]^2 + J[2,1]^2 + J[3,1]^2)
    else
        c1 = J[2,1] * J[3,2] - J[3,1] * J[2,2]
        c2 = J[3,1] * J[1,2] - J[1,1] * J[3,2]
        c3 = J[1,1] * J[2,2] - J[2,1] * J[1,2]
        sqrt(c1^2 + c2^2 + c3^2)
    end

    measure > sqrt(eps(Float64)) ||
        error(
            "ContactGap integration: degenerate slave element $(element.tag) " *
            "at integration point $q."
        )

    return xs, N, measure
end

function _contact_element_exclusion(c::Contact, element::_ContactElement)
    c.options.self_contact || return nothing

    ex = Set{Int}()
    for node in element.node_tags
        union!(ex, get(c.excluded, node, Set{Int}()))
    end
    return ex
end

@inline function _contact_gp_is_active(op::ContactGapOp, p::_ContactProjection)
    return op.active === :all || p.gap <= op.contact.options.activation_tol
end

function _contact_local_frame(
    p::_ContactProjection,
    pdim::Int,
    components::Symbol
    )

    if components === :normal
        Q = Matrix{Float64}(undef, 1, pdim)
        @inbounds for j in 1:pdim
            Q[1, j] = p.normal[j]
        end
        return Q
    end

    Q = Matrix{Float64}(undef, pdim, pdim)
    @inbounds for j in 1:pdim
        Q[1, j] = p.normal[j]
        Q[2, j] = p.tangent1[j]
    end

    if pdim == 3
        p.tangent2 === nothing &&
            error("ContactGap integration: missing second tangent in 3D contact.")
        @inbounds for j in 1:3
            Q[3, j] = p.tangent2[j]
        end
    end

    return Q
end

function _contact_gap_B(
    c::Contact,
    slave_element::_ContactElement,
    Ns,
    p::_ContactProjection,
    components::Symbol
    )

    pdim = c.U.pdim
    Q = _contact_local_frame(p, pdim, components)
    ncomp = size(Q, 1)

    master_element = c.master_elements[p.element_index]
    Nm = p.N

    ns = length(slave_element.node_tags)
    nm = length(master_element.node_tags)

    B = zeros(Float64, ncomp, pdim * (ns + nm))
    dofs = Vector{Int}(undef, pdim * (ns + nm))

    @inbounds for a in 1:ns
        node = slave_element.node_tags[a]
        Na = Ns[a]
        base = (a - 1) * pdim
        for j in 1:pdim
            dofs[base + j] = (node - 1) * pdim + j
            for k in 1:ncomp
                B[k, base + j] = Na * Q[k, j]
            end
        end
    end

    @inbounds for a in 1:nm
        node = master_element.node_tags[a]
        Na = Nm[a]
        base = pdim * ns + (a - 1) * pdim
        for j in 1:pdim
            dofs[base + j] = (node - 1) * pdim + j
            for k in 1:ncomp
                B[k, base + j] = -Na * Q[k, j]
            end
        end
    end

    return B, dofs
end

function _contact_id_B(P::Problem, slave_element::_ContactElement, Ns)
    pdim = P.pdim
    ns = length(slave_element.node_tags)

    B = zeros(Float64, pdim, pdim * ns)
    dofs = Vector{Int}(undef, pdim * ns)

    @inbounds for a in 1:ns
        node = slave_element.node_tags[a]
        Na = Ns[a]
        base = (a - 1) * pdim
        for j in 1:pdim
            B[j, base + j] = Na
            dofs[base + j] = (node - 1) * pdim + j
        end
    end

    return B, dofs
end

function _contact_scalar_at_gp(
    f::ScalarField,
    element::_ContactElement,
    N,
    step::Int
    )

    isNodal(f) ||
        error("ContactGap integration: ScalarField coefficients must be nodal.")

    s = f.nsteps == 1 ? 1 : step
    1 <= s <= f.nsteps ||
        error(
            "ContactGap integration: ScalarField coefficient does not contain step $step."
        )

    value = 0.0
    @inbounds for a in eachindex(element.node_tags)
        node = element.node_tags[a]
        value += N[a] * f.a[node, s]
    end
    return value
end

_contact_coefficient_factor(x::Number, element, N, step) = Float64(x)
_contact_coefficient_factor(x::ScalarField, element, N, step) =
    _contact_scalar_at_gp(x, element, N, step)

function _contact_coefficient_factor(A::AbstractMatrix, element, N, step)
    M = Matrix{Float64}(undef, size(A,1), size(A,2))
    @inbounds for j in axes(A,2), i in axes(A,1)
        x = A[i,j]
        if x isa Number
            M[i,j] = Float64(x)
        elseif x isa ScalarField
            M[i,j] = _contact_scalar_at_gp(x, element, N, step)
        else
            error(
                "ContactGap integration: coefficient matrix entries must be " *
                "Number or ScalarField, got $(typeof(x))."
            )
        end
    end
    return M
end

function _contact_coefficient_at_gp(coefficient, element, N, step)
    if coefficient isa Number || coefficient isa ScalarField || coefficient isa AbstractMatrix
        return _contact_coefficient_factor(coefficient, element, N, step)
    elseif coefficient isa AbstractVector
        isempty(coefficient) && error("ContactGap integration: empty coefficient chain.")
        value = nothing
        for factor in coefficient
            fv = _contact_coefficient_factor(factor, element, N, step)
            value = value === nothing ? fv : value * fv
        end
        return value
    else
        error(
            "ContactGap integration: unsupported coefficient type $(typeof(coefficient))."
        )
    end
end

function _contact_weight_at_gp(weight, element, N, step)
    weight === nothing && return 1.0
    value = _contact_coefficient_factor(weight, element, N, step)
    value isa Number || error("ContactGap integration: weight must evaluate to a scalar.")
    return Float64(value)
end

function _contact_local_bilinear(Bs, Cgp, Bu, scale::Float64)
    if Cgp isa Number
        size(Bs,1) == size(Bu,1) ||
            error(
                "ContactGap integration: scalar coefficient requires equal " *
                "operator dimensions $(size(Bs,1)) and $(size(Bu,1))."
            )
        return scale * Float64(Cgp) .* (transpose(Bs) * Bu)
    end

    size(Cgp,1) == size(Bs,1) ||
        error(
            "ContactGap integration: coefficient has $(size(Cgp,1)) rows, " *
            "expected $(size(Bs,1))."
        )
    size(Cgp,2) == size(Bu,1) ||
        error(
            "ContactGap integration: coefficient has $(size(Cgp,2)) columns, " *
            "expected $(size(Bu,1))."
        )

    return scale .* (transpose(Bs) * (Cgp * Bu))
end

function _contact_scatter_block!(I, J, V, Ke, rows, cols)
    @inbounds for j in eachindex(cols)
        cj = cols[j]
        for i in eachindex(rows)
            v = Ke[i,j]
            iszero(v) && continue
            push!(I, rows[i])
            push!(J, cj)
            push!(V, v)
        end
    end
    return nothing
end


@inline _contact_gauss_cache_key(gauss) =
    gauss isa Symbol ? gauss : Int(gauss)

function _contact_gp_warm_slots!(
    c::Contact,
    slave_element::_ContactElement,
    gauss,
    nip::Int
    )

    key = (_contact_gauss_cache_key(gauss), slave_element.tag)
    slots = get(c.gp_warm_start, key, nothing)

    if slots === nothing || length(slots) != nip
        slots = [_ContactGPWarmStart() for _ in 1:nip]
        c.gp_warm_start[key] = slots
    end

    return slots
end

function _contact_projection_at_gp(
    c::Contact,
    slave_element::_ContactElement,
    xs,
    workspace::_ContactProjectionWorkspace,
    warm_start::Union{Nothing,_ContactGPWarmStart}=nothing
    )

    ex = _contact_element_exclusion(c, slave_element)

    previous_element_index = 0
    previous_local_coordinate = nothing

    if warm_start !== nothing && warm_start.master_element_tag != 0
        previous_element_index =
            get(c.master_tag_to_index, warm_start.master_element_tag, 0)

        if previous_element_index != 0 &&
           (ex === nothing || !(previous_element_index in ex))
            previous_local_coordinate = warm_start.local_coordinate
        else
            previous_element_index = 0
        end
    end

    p = _contact_nearest_projection(
        c.tree,
        c.master_elements,
        xs,
        workspace;
        excluded=ex,
        normal_sign=c.options.normal_sign,
        model_dim=c.U.dim,
        projection_tol=c.options.projection_tol,
        projection_maxiter=c.options.projection_maxiter,
        topology=c.topology,
        topology_tol=c.options.topology_tol,
        topology_angle=c.options.topology_angle,
        previous_element_index=previous_element_index,
        previous_local_coordinate=previous_local_coordinate
    )

    p === nothing &&
        error(
            "ContactGap integration: no admissible master projection was found " *
            "for slave element $(slave_element.tag)."
        )

    if warm_start !== nothing
        warm_start.master_element_tag = p.element_tag
        warm_start.local_coordinate[1] = p.ξ[1]
        warm_start.local_coordinate[2] =
            length(p.ξ) >= 2 ? p.ξ[2] : 0.0
    end

    return p
end



# -----------------------------------------------------------------------------
# Optimized Gauss-point projection and CSC assembly helpers
# -----------------------------------------------------------------------------


@inline function _contact_resolve_element_chunk_size(
    nel::Int,
    num_threads::Int,
    element_chunk_size
    )

    chunk =
        element_chunk_size === :auto ?
        max(1, min(4096, cld(max(nel, 1), num_threads))) :
        element_chunk_size isa Integer ? Int(element_chunk_size) :
        error(
            "ContactGap integration: element_chunk_size must be " *
            ":auto or a positive integer."
        )

    chunk > 0 ||
        error("ContactGap integration: element_chunk_size must be positive.")

    return chunk
end


function _contact_prepare_gp_projection_data(
    c::Contact,
    gauss,
    threads,
    element_chunk_size=:auto
    )

    nel = length(c.slave_elements)
    num_threads = resolve_num_threads(threads)
    chunk = _contact_resolve_element_chunk_size(
        nel,
        num_threads,
        element_chunk_size
    )

    # Gmsh quadrature and basis data depend only on the element type and rule.
    qcache = Dict{Int,Any}()
    qdata_by_element = Vector{Any}(undef, nel)
    warm_by_element = Vector{Any}(undef, nel)

    @inbounds for e in 1:nel
        element = c.slave_elements[e]
        qdata = get!(qcache, element.etype) do
            _contact_quadrature(element, gauss)
        end
        qdata_by_element[e] = qdata
        warm_by_element[e] =
            _contact_gp_warm_slots!(c, element, gauss, qdata.nip)
    end

    projections =
        [Vector{_ContactProjection}(undef, qdata_by_element[e].nip)
         for e in 1:nel]
    measures =
        [Vector{Float64}(undef, qdata_by_element[e].nip)
         for e in 1:nel]

    all_elements = vcat(c.slave_elements, c.master_elements)
    workspaces = [
        _contact_projection_workspace(all_elements)
        for _ in 1:num_threads
    ]

    _run_workers(num_threads) do worker
        ws = workspaces[worker]
        stride = chunk * num_threads

        for chunk_first in
            (1 + (worker - 1) * chunk):stride:nel

            chunk_last = min(nel, chunk_first + chunk - 1)

            @inbounds for e in chunk_first:chunk_last
                element = c.slave_elements[e]
                qdata = qdata_by_element[e]
                warm_slots = warm_by_element[e]

                for q in 1:qdata.nip
                    xs, _, measure =
                        _contact_slave_geometry_at_gp(element, qdata, q)

                    projections[e][q] =
                        _contact_projection_at_gp(
                            c,
                            element,
                            xs,
                            ws,
                            warm_slots[q]
                        )
                    measures[e][q] = measure
                end
            end
        end
    end

    return qdata_by_element, projections, measures
end

function _contact_unique_master_indices(projections)
    ids = Int[]
    for p in projections
        push!(ids, p.element_index)
    end
    sort!(ids)
    unique!(ids)
    return ids
end

function _contact_finalize_node_csc_pattern(
    row_nodes_by_col_node::Vector{Vector{Int}},
    Ps::Problem,
    Pu::Problem
    )

    nrows = ndofs(Ps)
    ncols = ndofs(Pu)

    nrow_nodes, rem_s = divrem(nrows, Ps.pdim)
    ncol_nodes, rem_u = divrem(ncols, Pu.pdim)

    rem_s == 0 ||
        error("ContactGap CSC pattern: invalid test-space DoF count.")
    rem_u == 0 ||
        error("ContactGap CSC pattern: invalid trial-space DoF count.")
    length(row_nodes_by_col_node) == ncol_nodes ||
        error("ContactGap CSC pattern: invalid column-node adjacency size.")

    total_node_nonzeros = 0

    for rows in row_nodes_by_col_node
        sort!(rows)

        if !isempty(rows)
            write_pos = 1
            previous = rows[1]

            @inbounds for read_pos in 2:length(rows)
                current = rows[read_pos]

                if current != previous
                    write_pos += 1
                    rows[write_pos] = current
                    previous = current
                end
            end

            resize!(rows, write_pos)
        end

        total_node_nonzeros += length(rows)
    end

    total_nonzeros =
        Base.Checked.checked_mul(
            Base.Checked.checked_mul(
                total_node_nonzeros,
                Ps.pdim
            ),
            Pu.pdim
        )

    colptr = Vector{Int}(undef, ncols + 1)
    rowval = Vector{Int}(undef, total_nonzeros)

    position = 1

    @inbounds for col_node in 1:ncol_nodes
        rows = row_nodes_by_col_node[col_node]

        for comp_u in 1:Pu.pdim
            col = (col_node - 1) * Pu.pdim + comp_u
            colptr[col] = position

            for row_node in rows
                first_row = (row_node - 1) * Ps.pdim

                for comp_s in 1:Ps.pdim
                    rowval[position] = first_row + comp_s
                    position += 1
                end
            end
        end
    end

    colptr[ncols + 1] = position

    return SparseMatrixCSC(
        nrows,
        ncols,
        colptr,
        rowval,
        zeros(Float64, total_nonzeros)
    )
end

function _contact_build_contact_csc_pattern(
    c::Contact,
    projections_by_element,
    Ps::Problem,
    Pu::Problem
    )

    row_nodes_by_col_node = [Int[] for _ in 1:Pu.non]

    @inbounds for e in eachindex(c.slave_elements)
        slave_element = c.slave_elements[e]
        slave_nodes = slave_element.node_tags

        for master_index in
            _contact_unique_master_indices(projections_by_element[e])

            master_nodes = c.master_elements[master_index].node_tags
            nodes = vcat(slave_nodes, master_nodes)

            for col_node in nodes
                append!(row_nodes_by_col_node[col_node], nodes)
            end
        end
    end

    return _contact_finalize_node_csc_pattern(
        row_nodes_by_col_node,
        Ps,
        Pu
    )
end

function _contact_build_mixed_csc_pattern(
    c::Contact,
    projections_by_element,
    Ps::Problem,
    Pu::Problem
    )

    row_nodes_by_col_node = [Int[] for _ in 1:Pu.non]

    @inbounds for e in eachindex(c.slave_elements)
        slave_element = c.slave_elements[e]
        slave_nodes = slave_element.node_tags

        for master_index in
            _contact_unique_master_indices(projections_by_element[e])

            master_nodes = c.master_elements[master_index].node_tags
            trial_nodes = vcat(slave_nodes, master_nodes)

            for col_node in trial_nodes
                append!(row_nodes_by_col_node[col_node], slave_nodes)
            end
        end
    end

    return _contact_finalize_node_csc_pattern(
        row_nodes_by_col_node,
        Ps,
        Pu
    )
end

@inline function _contact_find_csc_position(
    rowval::Vector{Int},
    first::Int,
    last::Int,
    row::Int
    )

    lo = first
    hi = last

    @inbounds while lo <= hi
        mid = lo + ((hi - lo) >>> 1)
        current = rowval[mid]

        if current < row
            lo = mid + 1
        elseif current > row
            hi = mid - 1
        else
            return mid
        end
    end

    return 0
end

function _contact_scatter_csc!(
    nzval::Vector{Float64},
    colptr::Vector{Int},
    rowval::Vector{Int},
    Ke,
    rows,
    cols
    )

    @inbounds for j in eachindex(cols)
        col = cols[j]
        first = colptr[col]
        last = colptr[col + 1] - 1

        for i in eachindex(rows)
            value = Ke[i, j]
            iszero(value) && continue

            p = _contact_find_csc_position(
                rowval,
                first,
                last,
                rows[i]
            )

            p != 0 || error(
                "ContactGap CSC pattern no longer covers the current " *
                "slave-master projection. Reassemble once without " *
                "`csc_matrix`, then reuse the new pattern."
            )

            nzval[p] += value
        end
    end

    return nothing
end

function _contact_prepare_csc_buffers(
    K,
    pattern_builder,
    nrows::Int,
    ncols::Int,
    num_threads::Int
    )

    Kcsc =
        K === nothing ?
        pattern_builder() :
        K

    Kcsc isa SparseMatrixCSC{Float64,Int} ||
        error(
            "ContactGap CSC assembly requires " *
            "`SparseMatrixCSC{Float64,Int}` for csc_matrix."
        )

    size(Kcsc) == (nrows, ncols) ||
        error(
            "ContactGap CSC matrix has size $(size(Kcsc)); " *
            "expected ($(nrows), $(ncols))."
        )

    nzval_buffers = Vector{Vector{Float64}}(undef, num_threads)
    nzval_buffers[1] = Kcsc.nzval

    @inbounds for worker in 2:num_threads
        nzval_buffers[worker] = zeros(Float64, length(Kcsc.nzval))
    end

    return Kcsc, nzval_buffers
end


function _contact_contact_csc_worker!(
    c::Contact,
    op_u::ContactGapOp,
    op_s::ContactGapOp,
    coefficient,
    weight,
    qdata_by_element,
    projections_by_element,
    measures_by_element,
    nzval::Vector{Float64},
    colptr::Vector{Int},
    rowval::Vector{Int},
    chunk::Int,
    num_threads::Int,
    worker::Int
    )

    nel = length(c.slave_elements)
    stride = chunk * num_threads

    for chunk_first in
        (1 + (worker - 1) * chunk):stride:nel

        chunk_last = min(nel, chunk_first + chunk - 1)

        @inbounds for e in chunk_first:chunk_last
            slave_element = c.slave_elements[e]
            qdata = qdata_by_element[e]

            for q in 1:qdata.nip
                p = projections_by_element[e][q]

                _contact_gp_is_active(op_u, p) || continue

                nn = length(slave_element.node_tags)
                Ns = @view qdata.N[1:nn, q]

                Bu, cols =
                    _contact_gap_B(
                        c,
                        slave_element,
                        Ns,
                        p,
                        op_u.components
                    )
                Bs, rows =
                    _contact_gap_B(
                        c,
                        slave_element,
                        Ns,
                        p,
                        op_s.components
                    )

                Cgp =
                    _contact_coefficient_at_gp(
                        coefficient,
                        slave_element,
                        Ns,
                        c.step
                    )
                wcoef =
                    _contact_weight_at_gp(
                        weight,
                        slave_element,
                        Ns,
                        c.step
                    )
                scale =
                    measures_by_element[e][q] *
                    qdata.weights[q] *
                    wcoef

                Ke = _contact_local_bilinear(
                    Bs,
                    Cgp,
                    Bu,
                    scale
                )

                _contact_scatter_csc!(
                    nzval,
                    colptr,
                    rowval,
                    Ke,
                    rows,
                    cols
                )
            end
        end
    end

    return nothing
end


function _contact_mixed_csc_worker!(
    c::Contact,
    op_u::ContactGapOp,
    Ps::Problem,
    coefficient,
    weight,
    qdata_by_element,
    projections_by_element,
    measures_by_element,
    nzval::Vector{Float64},
    colptr::Vector{Int},
    rowval::Vector{Int},
    chunk::Int,
    num_threads::Int,
    worker::Int
    )

    nel = length(c.slave_elements)
    stride = chunk * num_threads

    for chunk_first in
        (1 + (worker - 1) * chunk):stride:nel

        chunk_last = min(nel, chunk_first + chunk - 1)

        @inbounds for e in chunk_first:chunk_last
            slave_element = c.slave_elements[e]
            qdata = qdata_by_element[e]

            for q in 1:qdata.nip
                p = projections_by_element[e][q]

                _contact_gp_is_active(op_u, p) || continue

                nn = length(slave_element.node_tags)
                Ns = @view qdata.N[1:nn, q]

                Bu, cols =
                    _contact_gap_B(
                        c,
                        slave_element,
                        Ns,
                        p,
                        op_u.components
                    )
                Bs, rows = _contact_id_B(Ps, slave_element, Ns)

                Cgp =
                    _contact_coefficient_at_gp(
                        coefficient,
                        slave_element,
                        Ns,
                        c.step
                    )
                wcoef =
                    _contact_weight_at_gp(
                        weight,
                        slave_element,
                        Ns,
                        c.step
                    )
                scale =
                    measures_by_element[e][q] *
                    qdata.weights[q] *
                    wcoef

                Ke = _contact_local_bilinear(
                    Bs,
                    Cgp,
                    Bu,
                    scale
                )

                _contact_scatter_csc!(
                    nzval,
                    colptr,
                    rowval,
                    Ke,
                    rows,
                    cols
                )
            end
        end
    end

    return nothing
end

function _contact_check_operator_domain(domain)
    domain === nothing ||
        error(
            "ContactGap already defines its integration manifold through " *
            "Contact.slave. Omit Ω/Γ from the ∫ call."
        )
    return nothing
end

function _contact_check_same_pair(a::ContactGapOp, b::ContactGapOp)
    a.contact === b.contact ||
        error("ContactGap bilinear form requires both operators to use the same Contact object.")
    a.active == b.active ||
        error("ContactGap bilinear form requires identical active policies on both sides.")
    return a.contact
end

"""
Specialized contact-contact assembler used automatically by expressions such as

    ∫(ContactGap(C) ⋅ D ⋅ ContactGap(C); updateFrom=u)

The contact operator is evaluated directly at slave Gauss points and includes
both slave and projected master interpolation in the local matrix.

With `assembly=:csc` the matrix is assembled directly into CSC storage with
worker-local value buffers and parallel reduction. Gauss-point closest-point
projections are warm-started from the preceding assembly. `updateFrom=u`
performs exactly one `updateContact!(C, u)` before the Gauss-point projection
and assembly pass.
"""

function assemble_operator(
    Pu::Problem,
    op_u::ContactGapOp,
    Ps::Problem,
    op_s::ContactGapOp;
    coefficient=1.0,
    weight=nothing,
    domain=nothing,
    gauss=:full,
    assembly::Symbol=:csc,
    threads=:auto,
    K=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto,
    updateFrom=nothing,
    kwargs...
    )

    _contact_check_operator_domain(domain)
    c = _contact_check_same_pair(op_s, op_u)

    Pu === c.U || error("ContactGap trial Problem must be Contact.U.")
    Ps === c.U || error("ContactGap test Problem must be Contact.U.")

    assembly in (:csc, :matrix, :ijv, :triplets) ||
        error("ContactGap integration: unsupported assembly mode $assembly.")

    updateFrom === nothing || begin
        updateFrom isa VectorField ||
            error("ContactGap integration: updateFrom must be a VectorField.")

        _contact_update_geometry!(c, updateFrom)
    end

    gmsh.model.setCurrent(c.U.name)

    if assembly === :csc
        num_threads = resolve_num_threads(threads)
        old_blas_threads = LinearAlgebra.BLAS.get_num_threads()

        try
            num_threads > 1 && LinearAlgebra.BLAS.set_num_threads(1)

            qdata_by_element, projections_by_element, measures_by_element =
                _contact_prepare_gp_projection_data(
                    c,
                    gauss,
                    threads,
                    element_chunk_size
                )

            Kcsc, nzval_buffers =
                _contact_prepare_csc_buffers(
                    K,
                    () -> _contact_build_contact_csc_pattern(
                        c,
                        projections_by_element,
                        Ps,
                        Pu
                    ),
                    ndofs(Ps),
                    ndofs(Pu),
                    num_threads
                )

            colptr = Kcsc.colptr
            rowval = Kcsc.rowval
            chunk = _contact_resolve_element_chunk_size(
                length(c.slave_elements),
                num_threads,
                element_chunk_size
            )

            _run_workers(num_threads) do worker
                _contact_contact_csc_worker!(
                    c,
                    op_u,
                    op_s,
                    coefficient,
                    weight,
                    qdata_by_element,
                    projections_by_element,
                    measures_by_element,
                    nzval_buffers[worker],
                    colptr,
                    rowval,
                    chunk,
                    num_threads,
                    worker
                )
            end

            reduce_csc_buffers!(
                Kcsc.nzval,
                nzval_buffers,
                num_threads
            )

            return SystemMatrix(Kcsc, Pu, Ps)
        finally
            num_threads > 1 &&
                LinearAlgebra.BLAS.set_num_threads(old_blas_threads)
        end
    end

    # Legacy triplet path retained for validation/debugging.
    I = Int[]
    J = Int[]
    V = Float64[]

    qdata_by_element, projections_by_element, measures_by_element =
        _contact_prepare_gp_projection_data(
            c,
            gauss,
            1,
            element_chunk_size
        )

    @inbounds for e in eachindex(c.slave_elements)
        slave_element = c.slave_elements[e]
        qdata = qdata_by_element[e]

        for q in 1:qdata.nip
            p = projections_by_element[e][q]
            _contact_gp_is_active(op_u, p) || continue

            nn = length(slave_element.node_tags)
            Ns = @view qdata.N[1:nn, q]

            Bu, cols =
                _contact_gap_B(
                    c,
                    slave_element,
                    Ns,
                    p,
                    op_u.components
                )
            Bs, rows =
                _contact_gap_B(
                    c,
                    slave_element,
                    Ns,
                    p,
                    op_s.components
                )

            Cgp =
                _contact_coefficient_at_gp(
                    coefficient,
                    slave_element,
                    Ns,
                    c.step
                )
            wcoef =
                _contact_weight_at_gp(
                    weight,
                    slave_element,
                    Ns,
                    c.step
                )
            scale =
                measures_by_element[e][q] *
                qdata.weights[q] *
                wcoef

            Ke = _contact_local_bilinear(
                Bs,
                Cgp,
                Bu,
                scale
            )

            _contact_scatter_block!(
                I,
                J,
                V,
                Ke,
                rows,
                cols
            )
        end
    end

    A = sparse(I, J, V, ndofs(Ps), ndofs(Pu))
    dropzeros!(A)
    return SystemMatrix(A, Pu, Ps)
end


"""
Specialized mixed assembler for a standard identity test field and a contact-gap
trial operator. It is used automatically by

    B = ∫(Λ ⋅ ContactGap(C; components=:normal))

and returns a rectangular `SystemMatrix` mapping the displacement field to the
multiplier test field.
"""

function assemble_operator(
    Pu::Problem,
    op_u::ContactGapOp,
    Ps::Problem,
    op_s::IdOp;
    coefficient=1.0,
    weight=nothing,
    domain=nothing,
    gauss=:full,
    assembly::Symbol=:csc,
    threads=:auto,
    K=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto,
    updateFrom=nothing,
    kwargs...
    )

    _contact_check_operator_domain(domain)
    c = op_u.contact

    Pu === c.U || error("ContactGap trial Problem must be Contact.U.")
    Ps.name == c.U.name ||
        error("Multiplier/test Problem must use the same Gmsh model as Contact.U.")
    Ps.non == c.U.non ||
        error("Multiplier/test Problem must use the same mesh nodes as Contact.U.")

    assembly in (:csc, :matrix, :ijv, :triplets) ||
        error("ContactGap integration: unsupported assembly mode $assembly.")

    updateFrom === nothing || begin
        updateFrom isa VectorField ||
            error("ContactGap integration: updateFrom must be a VectorField.")

        _contact_update_geometry!(c, updateFrom)
    end

    gmsh.model.setCurrent(c.U.name)

    if assembly === :csc
        num_threads = resolve_num_threads(threads)
        old_blas_threads = LinearAlgebra.BLAS.get_num_threads()

        try
            num_threads > 1 && LinearAlgebra.BLAS.set_num_threads(1)

            qdata_by_element, projections_by_element, measures_by_element =
                _contact_prepare_gp_projection_data(
                    c,
                    gauss,
                    threads,
                    element_chunk_size
                )

            Bcsc, nzval_buffers =
                _contact_prepare_csc_buffers(
                    K,
                    () -> _contact_build_mixed_csc_pattern(
                        c,
                        projections_by_element,
                        Ps,
                        Pu
                    ),
                    ndofs(Ps),
                    ndofs(Pu),
                    num_threads
                )

            colptr = Bcsc.colptr
            rowval = Bcsc.rowval
            chunk = _contact_resolve_element_chunk_size(
                length(c.slave_elements),
                num_threads,
                element_chunk_size
            )

            _run_workers(num_threads) do worker
                _contact_mixed_csc_worker!(
                    c,
                    op_u,
                    Ps,
                    coefficient,
                    weight,
                    qdata_by_element,
                    projections_by_element,
                    measures_by_element,
                    nzval_buffers[worker],
                    colptr,
                    rowval,
                    chunk,
                    num_threads,
                    worker
                )
            end

            reduce_csc_buffers!(
                Bcsc.nzval,
                nzval_buffers,
                num_threads
            )

            return SystemMatrix(Bcsc, Pu, Ps)
        finally
            num_threads > 1 &&
                LinearAlgebra.BLAS.set_num_threads(old_blas_threads)
        end
    end

    # Legacy triplet path retained for validation/debugging.
    I = Int[]
    J = Int[]
    V = Float64[]

    qdata_by_element, projections_by_element, measures_by_element =
        _contact_prepare_gp_projection_data(
            c,
            gauss,
            1,
            element_chunk_size
        )

    @inbounds for e in eachindex(c.slave_elements)
        slave_element = c.slave_elements[e]
        qdata = qdata_by_element[e]

        for q in 1:qdata.nip
            p = projections_by_element[e][q]
            _contact_gp_is_active(op_u, p) || continue

            nn = length(slave_element.node_tags)
            Ns = @view qdata.N[1:nn, q]

            Bu, cols =
                _contact_gap_B(
                    c,
                    slave_element,
                    Ns,
                    p,
                    op_u.components
                )
            Bs, rows = _contact_id_B(Ps, slave_element, Ns)

            Cgp =
                _contact_coefficient_at_gp(
                    coefficient,
                    slave_element,
                    Ns,
                    c.step
                )
            wcoef =
                _contact_weight_at_gp(
                    weight,
                    slave_element,
                    Ns,
                    c.step
                )
            scale =
                measures_by_element[e][q] *
                qdata.weights[q] *
                wcoef

            Ke = _contact_local_bilinear(
                Bs,
                Cgp,
                Bu,
                scale
            )

            _contact_scatter_block!(
                I,
                J,
                V,
                Ke,
                rows,
                cols
            )
        end
    end

    A = sparse(I, J, V, ndofs(Ps), ndofs(Pu))
    dropzeros!(A)
    return SystemMatrix(A, Pu, Ps)
end

function _contact_transpose_coefficient(coefficient)
    if coefficient isa Number || coefficient isa ScalarField
        return coefficient
    elseif coefficient isa AbstractMatrix
        return transpose(coefficient)
    elseif coefficient isa AbstractVector
        out = Any[]
        sizehint!(out, length(coefficient))
        for factor in Iterators.reverse(coefficient)
            if factor isa Number || factor isa ScalarField
                push!(out, factor)
            elseif factor isa AbstractMatrix
                push!(out, transpose(factor))
            else
                error("ContactGap integration: unsupported coefficient factor $(typeof(factor)).")
            end
        end
        return out
    end
    error("ContactGap integration: unsupported coefficient type $(typeof(coefficient)).")
end

"""
Transpose mixed contact assembler. This makes expressions with ContactGap on the
test side work consistently with the standard LLFEM DSL.
"""
function assemble_operator(
    Pu::Problem,
    op_u::IdOp,
    Ps::Problem,
    op_s::ContactGapOp;
    coefficient=1.0,
    weight=nothing,
    domain=nothing,
    gauss=:full,
    assembly::Symbol=:csc,
    threads=:auto,
    kwargs...
    )

    B = assemble_operator(
        Ps,
        op_s,
        Pu,
        op_u;
        coefficient=_contact_transpose_coefficient(coefficient),
        weight=weight,
        domain=domain,
        gauss=gauss,
        assembly=assembly,
        threads=threads,
        kwargs...
    )

    A = sparse(transpose(B.A))
    return SystemMatrix(A, Pu, Ps)
end

# -----------------------------------------------------------------------------
# ContactGap field evaluation / postprocessing
# -----------------------------------------------------------------------------

function _contact_absolute_coordinates_from_field(
    c::Contact,
    u::VectorField;
    step::Int,
    absolute::Bool
    )

    un = isNodal(u) ? u : elementsToNodes(u)
    size(un.a, 1) == ndofs(c.U) ||
        error("ContactGap field evaluation: incompatible VectorField size.")

    s = un.nsteps == 1 ? 1 : step
    1 <= s <= un.nsteps ||
        error("ContactGap field evaluation: field does not contain step $step.")

    X = _contact_node_coordinates(c.U)
    pdim = c.U.pdim

    @inbounds for node in 1:c.U.non
        base = (node - 1) * pdim
        for j in 1:pdim
            if absolute
                X[j, node] = un.a[base + j, s]
            else
                X[j, node] += un.a[base + j, s]
            end
        end
    end

    return X
end

function _contact_gap_gauss_field_values(
    op::ContactGapOp,
    u::VectorField;
    gauss=:full,
    step::Int=u.nsteps,
    absolute::Bool=false
    )

    absolute &&
        error(
            "ContactGap Gauss-point evaluation currently expects a displacement field. " *
            "Use absolute=false and call updateContact!(C, u) first."
        )

    c = op.contact

    # Keep the same explicit update semantics as the nodal ContactGap evaluator.
    # The Gauss-point path uses the current Contact geometry/search caches, i.e.
    # the same geometry used by ∫(ContactGap(C) ⋅ ... ⋅ ContactGap(C)).
    if !(c.displacement === u && c.step == step)
        @warn(
            "ContactGap(C, u; gauss=...) uses the current Contact geometry. " *
            "Call updateContact!(C, u; step=step) first when u differs from the " *
            "state stored in C."
        )
    end

    gmsh.model.setCurrent(c.U.name)

    ncomp = op.components === :normal ? 1 : c.U.pdim
    pdim = c.U.pdim

    # ------------------------------------------------------------------
    # Global L2 projection on the slave Lagrange space
    #
    #     M g_h = b,
    #
    #     M = ∫_Γ Nᵀ N dΓ,
    #     b = ∫_Γ Nᵀ g_q dΓ.
    #
    # Only slave-surface nodes are included in the reduced projection
    # system. Shared nodes therefore receive one common coefficient and the
    # reconstructed field is C0-continuous across slave element boundaries.
    # ------------------------------------------------------------------

    slave_nodes = c.slave_nodes
    nsurf = length(slave_nodes)

    node_to_local = Dict{Int,Int}(
        node => i for (i, node) in enumerate(slave_nodes)
    )

    I = Int[]
    J = Int[]
    V = Float64[]
    rhs = zeros(Float64, nsurf, ncomp)

    # Reserve roughly one dense local mass block per slave element. This is
    # only a hint; mixed element types/orders remain supported.
    if !isempty(c.slave_elements)
        nnmax = maximum(length(e.node_tags) for e in c.slave_elements)
        sizehint!(I, length(c.slave_elements) * nnmax^2)
        sizehint!(J, length(c.slave_elements) * nnmax^2)
        sizehint!(V, length(c.slave_elements) * nnmax^2)
    end

    workspace =
        _contact_projection_workspace(vcat(c.slave_elements, c.master_elements))
    qcache = Dict{Int,Any}()

    @inbounds for slave_element in c.slave_elements
        qdata = get!(qcache, slave_element.etype) do
            _contact_quadrature(slave_element, gauss)
        end

        nn = length(slave_element.node_tags)
        local_ids = Vector{Int}(undef, nn)
        for a in 1:nn
            node = slave_element.node_tags[a]
            local_ids[a] = get(node_to_local, node, 0)
            local_ids[a] != 0 ||
                error(
                    "ContactGap Gauss-point projection: slave element " *
                    "$(slave_element.tag) contains node $node which is not " *
                    "present in Contact.slave_nodes."
                )
        end

        Me = zeros(Float64, nn, nn)
        be = zeros(Float64, nn, ncomp)

        for q in 1:qdata.nip
            xs, Ns, measure =
                _contact_slave_geometry_at_gp(slave_element, qdata, q)
            p = _contact_projection_at_gp(c, slave_element, xs, workspace)

            wq = measure * qdata.weights[q]

            gq = zeros(Float64, ncomp)
            gq[1] = p.gap

            if op.components === :all
                dx1 = xs[1] - p.x[1]
                dx2 = xs[2] - p.x[2]
                dx3 = xs[3] - p.x[3]

                gq[2] =
                    dx1 * p.tangent1[1] +
                    dx2 * p.tangent1[2] +
                    dx3 * p.tangent1[3]

                if pdim == 3
                    p.tangent2 === nothing &&
                        error(
                            "ContactGap Gauss-point evaluation: missing second " *
                            "tangent in 3D contact."
                        )
                    gq[3] =
                        dx1 * p.tangent2[1] +
                        dx2 * p.tangent2[2] +
                        dx3 * p.tangent2[3]
                end
            end

            # Local mass matrix and right-hand side.
            for a in 1:nn
                Na = Ns[a]

                for k in 1:ncomp
                    be[a, k] += wq * Na * gq[k]
                end

                for b in 1:nn
                    Me[a, b] += wq * Na * Ns[b]
                end
            end
        end

        # Scatter the local projection system to the reduced global
        # slave-surface system.
        for a in 1:nn
            ia = local_ids[a]

            for k in 1:ncomp
                rhs[ia, k] += be[a, k]
            end

            for b in 1:nn
                push!(I, ia)
                push!(J, local_ids[b])
                push!(V, Me[a, b])
            end
        end
    end

    M = sparse(I, J, V, nsurf, nsurf)

    # The consistent slave-surface mass matrix is symmetric positive definite
    # when the chosen quadrature sufficiently resolves the Lagrange space.
    # A failed factorization is therefore a useful indication that a denser
    # Gauss rule is required (e.g. for high-order elements with :reduced).
    coeff = try
        F = cholesky(Symmetric(M))
        F \ rhs
    catch err
        error(
            "ContactGap Gauss-point L2 projection failed. The slave-surface " *
            "projection matrix is singular or not positive definite for " *
            "gauss=$gauss. Use a denser quadrature rule (for example gauss=0, " *
            "gauss=2, or larger). Original error: $(sprint(showerror, err))"
        )
    end

    if op.components === :normal
        values = zeros(Float64, c.U.non, 1)
        for (i, node) in enumerate(slave_nodes)
            values[node, 1] = coeff[i, 1]
        end

        return ScalarField(
            Matrix{Float64}[],
            values,
            [0.0],
            Int[],
            1,
            :scalar,
            c.U
        )
    end

    values = zeros(Float64, ndofs(c.U), 1)
    for (i, node) in enumerate(slave_nodes)
        base = (node - 1) * pdim
        for k in 1:pdim
            values[base + k, 1] = coeff[i, k]
        end
    end

    type = pdim == 2 ? :v2D : :v3D
    return VectorField(
        Matrix{Float64}[],
        values,
        [0.0],
        Int[],
        1,
        type,
        c.U
    )
end

function _contact_gap_field_values(
    op::ContactGapOp,
    u::VectorField;
    step::Int=u.nsteps,
    absolute::Bool=false
    )

    c = op.contact
    X = _contact_absolute_coordinates_from_field(c, u; step=step, absolute=absolute)
    pdim = c.U.pdim

    if op.components === :normal
        values = zeros(Float64, c.U.non, 1)
    else
        values = zeros(Float64, ndofs(c.U), 1)
    end

    @inbounds for (i, slave_node) in enumerate(c.slave_nodes)
        p = c.projections[i]
        master_element = c.master_elements[p.element_index]

        xm1 = 0.0
        xm2 = 0.0
        xm3 = 0.0
        for a in eachindex(master_element.node_tags)
            node = master_element.node_tags[a]
            Na = p.N[a]
            xm1 += Na * X[1, node]
            xm2 += Na * X[2, node]
            xm3 += Na * X[3, node]
        end

        dx1 = X[1, slave_node] - xm1
        dx2 = X[2, slave_node] - xm2
        dx3 = X[3, slave_node] - xm3

        dn = dx1 * p.normal[1] + dx2 * p.normal[2] + dx3 * p.normal[3]

        if op.components === :normal
            values[slave_node, 1] = dn
        else
            base = (slave_node - 1) * pdim
            values[base + 1, 1] = dn
            values[base + 2, 1] =
                dx1 * p.tangent1[1] + dx2 * p.tangent1[2] + dx3 * p.tangent1[3]

            if pdim == 3
                p.tangent2 === nothing &&
                    error("ContactGap field evaluation: missing second tangent in 3D contact.")
                values[base + 3, 1] =
                    dx1 * p.tangent2[1] + dx2 * p.tangent2[2] + dx3 * p.tangent2[3]
            end
        end
    end

    if op.components === :normal
        return ScalarField(
            Matrix{Float64}[],
            values,
            [0.0],
            Int[],
            1,
            :scalar,
            c.U
        )
    end

    type = pdim == 2 ? :v2D : :v3D
    return VectorField(
        Matrix{Float64}[],
        values,
        [0.0],
        Int[],
        1,
        type,
        c.U
    )
end

"""
    (G::OpApplied)(u::VectorField; gauss=nothing, step=u.nsteps, absolute=false)

Evaluate a `ContactGap(C)` operator on a displacement field for postprocessing.

With `gauss=nothing` (default), the gap is evaluated at slave nodes and returned
as a nodal `ScalarField` for `components=:normal`, or a nodal `VectorField` for
`components=:all`.

With `gauss=:full`, `gauss=:reduced`, or an integer Gauss-order offset, the gap
is evaluated at the same slave Gauss points used by contact integration and
globally L2-projected onto the continuous slave-side Lagrange space. The result
is a nodal field whose coefficients solve

    ∫_Γ Nᵀ N dΓ * g_h = ∫_Γ Nᵀ g_q dΓ.

Increasing the Gauss order improves the numerical projection of the generally
non-polynomial closest-point gap without changing the interpolation order.

For ordinary nodal evaluation `u` is a displacement field and the reference
coordinates are added internally. Set `absolute=true` only when `u` already
stores absolute nodal positions. Gauss-point evaluation currently requires
`absolute=false`.

Call `updateContact!(C, u)` first when the closest-point geometry itself should
be refreshed to the supplied displacement.
"""
function (G::OpApplied)(
    u::VectorField;
    gauss=nothing,
    step::Int=u.nsteps,
    absolute::Bool=false
    )

    G.op isa ContactGapOp ||
        error("Only ContactGap applied operators are callable on VectorField objects.")

    if gauss === nothing
        return _contact_gap_field_values(
            G.op,
            u;
            step=step,
            absolute=absolute
        )
    end

    return _contact_gap_gauss_field_values(
        G.op,
        u;
        gauss=gauss,
        step=step,
        absolute=absolute
    )
end

"""
    ContactGap(C::Contact, u::VectorField; components=:normal, active=:current,
               gauss=nothing, step=u.nsteps, absolute=false)

Convenience form equivalent to creating `ContactGap(C; ...)` and immediately
evaluating it on `u` for postprocessing.

Use `gauss=nothing` for the original nodal field, or specify `gauss=:full`,
`gauss=:reduced`, or an integer offset (for example `gauss=2`) to evaluate the
gap at slave Gauss points and globally L2-project it onto the continuous
slave-side Lagrange space.
"""
function ContactGap(
    c::Contact,
    u::VectorField;
    components::Symbol=:normal,
    active::Symbol=:current,
    gauss=nothing,
    step::Int=u.nsteps,
    absolute::Bool=false
    )

    G = ContactGap(c; components=components, active=active)
    return G(u; gauss=gauss, step=step, absolute=absolute)
end
