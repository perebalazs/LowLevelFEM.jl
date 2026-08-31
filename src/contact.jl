###############################################################################
#                                                                             #
#                              Contact                                        #
#                                                                             #
###############################################################################

export contact

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
    coords::Matrix{Float64}          # 3 × number of element nodes
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
    gap::Float64
    distance2::Float64
end

# -----------------------------------------------------------------------------
# Public interface
# -----------------------------------------------------------------------------

"""
    contact(model::Problem, slave::String, master::String; kwargs...)
        -> gap::ScalarField, G::SystemMatrix, normalVec::VectorField

    contact(model::Problem, test_model::Problem,
            slave::String, master::String; kwargs...)
        -> gap::ScalarField, G::SystemMatrix, normalVec::VectorField

Construct the geometric quantities required by a frictionless node-to-manifold
contact formulation.

For every node of the `slave` physical group, the closest point on the
`master` curve or surface is determined. The function returns:

- `gap`: signed normal gap as an element-wise `ScalarField` on the slave side,
- `G`: sparse contact kinematic matrix as a `SystemMatrix`,
- `normalVec`: contact normal as an element-wise `VectorField` on the slave side.

The contact matrix satisfies, for each slave node `s`,

    (G * u)_s = n_s ⋅ (u_s - Σ_a N_a u_a)

where `N_a` are the master-element shape functions evaluated at the closest
point and `n_s` is the contact normal.

With one supplied field, `G.model === G.test_model === model`. This form is
intended for penalty formulations, e.g.

    gap, G, n = contact(Pu, "slave", "master")
    Kc = c_n * G' * G
    Ktot = K + Kc

`G` itself is rectangular: it has one row per global mesh node and one column
per displacement degree of freedom. Only `G' * G` is square.

With two supplied fields, the second field is interpreted as the scalar
Lagrange-multiplier field:

    gap, G, n = contact(Pu, Pp, "slave", "master")

and the metadata are

    G.model      === Pu
    G.test_model === Pp

so `G` and `G'` can be inserted directly into a multifield block system.

# Supported geometry

The master physical group may contain Gmsh Lagrange line, triangle or
quadrilateral elements of arbitrary order. The slave physical group may be a
curve or a surface. Consequently the same implementation covers curve-curve,
curve-surface, surface-curve and surface-surface mappings. A planar 2D contact
problem is the special case of curve-curve contact in the xy plane.

For a 3D master curve no unique geometric normal exists. In that case the
contact normal is chosen in the direction from the closest master point to the
slave point. If the two points coincide, a stable vector perpendicular to the
curve tangent is used.

# Self-contact

If `slave == master`, self-contact mode is enabled automatically. Master
elements belonging to the local topological neighbourhood of the current slave
node are excluded from the closest-point search. This is a one-pass
node-to-manifold self-contact mapping; reciprocal self-contact constraints are
not explicitly deduplicated.

# Keyword Arguments

- `normal_sign::Real=1.0`: multiply all contact normals by this value. Use `-1`
  to reverse the master-side orientation.
- `self_contact::Bool=(slave == master)`: enable local-topology exclusion.
- `self_exclusion_layers::Int=1`: number of additional node-connected master
  element layers excluded in self-contact. The elements containing the slave
  node are always excluded.
- `aabb_padding::Real=0.05`: relative padding of master-element AABBs. This is
  useful for curved higher-order elements whose geometry may extend slightly
  outside the nodal bounding box.
- `leaf_size::Int=8`: maximum number of elements stored in an AABB-tree leaf.
- `projection_tol::Real=1e-10`: relative tolerance of the closest-point solver.
- `projection_maxiter::Int=40`: maximum projected Gauss-Newton iterations for
  one starting point.

# Sign convention

For surfaces and planar curves,

    gap = (x_slave - x_master) ⋅ n

so a positive gap denotes separation when the master normal points toward the
slave side. Reverse the convention with `normal_sign=-1`.
"""
function contact(
    model::Problem,
    slave::String,
    master::String;
    kwargs...
    )

    return _contact_impl(
        model,
        model,
        slave,
        master;
        lagrange=false,
        kwargs...
    )
end

function contact(
    model::Problem,
    test_model::Problem,
    slave::String,
    master::String;
    kwargs...
    )

    return _contact_impl(
        model,
        test_model,
        slave,
        master;
        lagrange=true,
        kwargs...
    )
end

# -----------------------------------------------------------------------------
# Main implementation
# -----------------------------------------------------------------------------

function _contact_impl(
    model::Problem,
    test_model::Problem,
    slave::String,
    master::String;
    lagrange::Bool,
    normal_sign::Real=1.0,
    self_contact::Bool=(slave == master),
    self_exclusion_layers::Int=1,
    aabb_padding::Real=0.05,
    leaf_size::Int=8,
    projection_tol::Real=1e-10,
    projection_maxiter::Int=40
    )

    _contact_check_models(model, test_model, lagrange)

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

    isfinite(normal_sign) && abs(abs(Float64(normal_sign)) - 1.0) <= 10 * eps(Float64) ||
        error("contact: normal_sign must be either +1 or -1.")

    gmsh.model.setCurrent(model.name)

    nodecoords = _contact_node_coordinates(model)

    slave_elements, slave_dim =
        _contact_group_elements(model, slave, nodecoords; aabb_padding=0.0)

    master_elements, master_dim =
        _contact_group_elements(model, master, nodecoords; aabb_padding=aabb_padding)

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
            normal_sign=Float64(normal_sign),
            model_dim=model.dim,
            projection_tol=Float64(projection_tol),
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

    gap = _contact_gap_field(
        lagrange ? test_model : model,
        slave_elements,
        projections
    )

    normalVec = _contact_normal_field(
        model,
        slave_elements,
        projections
    )

    Gmat = _contact_matrix(
        model,
        test_model,
        slave_nodes,
        projections,
        master_elements;
        lagrange=lagrange
    )

    if lagrange
        G = SystemMatrix(Gmat, model, test_model)
    else
        G = SystemMatrix(Gmat, model, model)
    end

    return gap, G, normalVec
end

# -----------------------------------------------------------------------------
# Model and mesh helpers
# -----------------------------------------------------------------------------

"""
    _contact_check_models(model, test_model, lagrange)

Validate field dimensions and mesh compatibility required by the contact
operator.
"""
function _contact_check_models(
    model::Problem,
    test_model::Problem,
    lagrange::Bool
    )

    model.pdim in (2, 3) ||
        error(
            "contact: model must be a 2D or 3D vector displacement field; " *
            "got pdim=$(model.pdim)."
        )

    model.name == test_model.name ||
        error("contact: model and test_model must use the same Gmsh model.")

    model.non == test_model.non ||
        error("contact: model and test_model must use the same mesh nodes.")

    if lagrange
        test_model.pdim == 1 ||
            error(
                "contact: the Lagrange-multiplier test_model must be scalar " *
                "(pdim == 1)."
            )
    end

    return nothing
end

"""
    _contact_node_coordinates(problem) -> Matrix{Float64}

Return global mesh-node coordinates in a `3 × problem.non` matrix indexed by
Gmsh node tag.
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
    _contact_group_elements(problem, phName, nodecoords; aabb_padding)
        -> elements, dimension

Collect curve or surface elements belonging to a physical group.
"""
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

function _contact_local_node_matrix(local, dim::Int, num_nodes::Int)
    data = Float64.(local)

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
Gauss-Newton iterations in the reference element.
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

    gap = dot(separation, normal)

    return _ContactProjection(
        element_index,
        element.tag,
        best_ξ,
        best_x,
        best_N,
        normal,
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

# -----------------------------------------------------------------------------
# Contact normal
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# LLFEM output objects
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

function _contact_normal_field(
    model::Problem,
    slave_elements::Vector{_ContactElement},
    projections::Dict{Int,_ContactProjection}
    )

    A = Vector{Matrix{Float64}}(undef, length(slave_elements))
    num_elem = Vector{Int}(undef, length(slave_elements))

    @inbounds for (i, e) in enumerate(slave_elements)
        values = Matrix{Float64}(undef, 3 * length(e.node_tags), 1)

        for (a, node) in enumerate(e.node_tags)
            n = projections[node].normal
            values[3a - 2, 1] = n[1]
            values[3a - 1, 1] = n[2]
            values[3a, 1] = n[3]
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

function _contact_matrix(
    model::Problem,
    test_model::Problem,
    slave_nodes::Vector{Int},
    projections::Dict{Int,_ContactProjection},
    master_elements::Vector{_ContactElement};
    lagrange::Bool
    )

    nrows = lagrange ? ndofs(test_model) : model.non
    ncols = ndofs(model)

    lagrange && nrows != model.non &&
        error(
            "contact: the Lagrange multiplier field must provide exactly one " *
            "degree of freedom per mesh node."
        )

    pdim = model.pdim

    Iidx = Int[]
    Jidx = Int[]
    Vval = Float64[]

    # Each slave-node row contains pdim slave coefficients and pdim coefficients
    # for every node of the corresponding master element.
    estimated = 0
    for node in slave_nodes
        p = projections[node]
        estimated += pdim * (1 + length(master_elements[p.element_index].node_tags))
    end
    sizehint!(Iidx, estimated)
    sizehint!(Jidx, estimated)
    sizehint!(Vval, estimated)

    @inbounds for slave_node in slave_nodes
        p = projections[slave_node]
        e = master_elements[p.element_index]
        row = slave_node
        n = p.normal

        for c in 1:pdim
            val = n[c]
            iszero(val) && continue
            push!(Iidx, row)
            push!(Jidx, (slave_node - 1) * pdim + c)
            push!(Vval, val)
        end

        for (a, master_node) in enumerate(e.node_tags)
            Na = p.N[a]
            iszero(Na) && continue

            for c in 1:pdim
                val = -Na * n[c]
                iszero(val) && continue
                push!(Iidx, row)
                push!(Jidx, (master_node - 1) * pdim + c)
                push!(Vval, val)
            end
        end
    end

    G = sparse(Iidx, Jidx, Vval, nrows, ncols)
    dropzeros!(G)
    return G
end
