export reductionMatrices

"""
    reductionMatrices(P::Problem)

Construct sparse transformation matrices for reducing the polynomial order
of a C0 Lagrange finite element field from order `p` to `p - 1`.

The returned matrices satisfy

    u_full = T * u_reduced
    u_reduced = R * u_full

for fields representable in the reduced space.

The transformation is constructed from the Gmsh Lagrange basis functions.
All element types belonging to the problem must have the same polynomial
order. An error is thrown for first-order meshes.

The matrices are expanded automatically according to `P.pdim`.
"""
function reductionMatrices(P::Problem)

    gmsh.model.setCurrent(P.name)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    function element_family(name::String)
        if occursin("Line", name)
            return "Line"
        elseif occursin("Triangle", name)
            return "Triangle"
        elseif occursin("Quadrilateral", name) ||
               occursin("Quadrangle", name)
            return "Quadrangle"
        elseif occursin("Tetrahedron", name)
            return "Tetrahedron"
        elseif occursin("Hexahedron", name)
            return "Hexahedron"
        elseif occursin("Prism", name)
            return "Prism"
        elseif occursin("Pyramid", name)
            return "Pyramid"
        else
            error(
                "reductionMatrices: unsupported Gmsh element family \"$name\"."
            )
        end
    end

    function local_to_3d(localCoord, dim, n)
        ξ = zeros(Float64, 3 * n)

        @inbounds for a in 1:n
            o3 = 3 * (a - 1)
            od = dim * (a - 1)

            for d in 1:dim
                ξ[o3+d] = localCoord[od+d]
            end
        end

        return ξ
    end

    # ------------------------------------------------------------------
    # Node coordinates
    #
    # LLFEM already uses node tags directly for global DOF numbering,
    # therefore the same assumption is used here: node tags are in
    # 1:P.non.
    # ------------------------------------------------------------------

    nodeTags, coord, _ =
        gmsh.model.mesh.getNodes(-1, -1, false, false)

    nodeCoord = Matrix{Float64}(undef, 3, P.non)

    @inbounds for i in eachindex(nodeTags)

        node = Int(nodeTags[i])

        1 <= node <= P.non ||
            error(
                "reductionMatrices: node tag $node is incompatible with P.non=$(P.non)."
            )

        nodeCoord[1, node] = coord[3*i-2]
        nodeCoord[2, node] = coord[3*i-1]
        nodeCoord[3, node] = coord[3*i]
    end

    # ------------------------------------------------------------------
    # Element-type data
    #
    # Parallel concretely typed arrays are used instead of Dict{Int,Any}.
    # Each element stores only the integer index of its element type.
    # ------------------------------------------------------------------

    type_id = Dict{Int,Int}()

    nHighs = Int[]
    nLows = Int[]
    nPrimarys = Int[]
    nPrimaryLows = Int[]

    Tes = Matrix{Float64}[]
    Res = Matrix{Float64}[]

    # Reusable workspaces for physical node coordinates
    Xhigh_work = Matrix{Float64}[]
    Xlow_work = Matrix{Float64}[]

    # ------------------------------------------------------------------
    # Elements belonging to the Problem
    # ------------------------------------------------------------------

    elem_type_id = Int[]
    elem_nodes = Vector{UInt64}[]

    seen_elements = Set{UInt64}()

    p_global = 0

    xmin = Inf
    ymin = Inf
    zmin = Inf

    xmax = -Inf
    ymax = -Inf
    zmax = -Inf

    for mat in P.material

        dimTags =
            gmsh.model.getEntitiesForPhysicalName(mat.phName)

        for (edim, etag) in dimTags

            edim == P.dim || continue

            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for it in eachindex(elemTypes)

                et = Int(elemTypes[it])

                # ------------------------------------------------------
                # Initialize element-type cache once
                # ------------------------------------------------------

                tid = get(type_id, et, 0)

                if tid == 0

                    name, dim, p, nHigh, ξHigh0, nPrimary =
                        gmsh.model.mesh.getElementProperties(et)

                    p > 1 ||
                        error(
                            "reductionMatrices: polynomial order must be greater than one."
                        )

                    if p_global == 0
                        p_global = p
                    elseif p != p_global
                        error(
                            "reductionMatrices: the mesh must have homogeneous polynomial order; " *
                            "found both order $p_global and order $p."
                        )
                    end

                    q = p - 1

                    family = element_family(name)

                    etLow =
                        gmsh.model.mesh.getElementType(
                            family,
                            q,
                            false
                        )

                    _, dimLow, qCheck,
                    nLow, ξLow0, nPrimaryLow =
                        gmsh.model.mesh.getElementProperties(etLow)

                    dimLow == dim ||
                        error(
                            "reductionMatrices: incompatible reduced element dimension."
                        )

                    qCheck == q ||
                        error(
                            "reductionMatrices: could not construct order-$q element."
                        )

                    # --------------------------------------------------
                    # Local prolongation matrix
                    #
                    # u_high = Te * u_low
                    # --------------------------------------------------

                    ξHigh =
                        local_to_3d(
                            ξHigh0,
                            dim,
                            nHigh
                        )

                    _, funT, _ =
                        gmsh.model.mesh.getBasisFunctions(
                            etLow,
                            ξHigh,
                            "Lagrange"
                        )

                    Te =
                        Matrix(
                            transpose(
                                reshape(
                                    funT,
                                    nLow,
                                    nHigh
                                )
                            )
                        )

                    # --------------------------------------------------
                    # Local restriction matrix
                    #
                    # u_low = Re * u_high
                    # --------------------------------------------------

                    ξLow =
                        local_to_3d(
                            ξLow0,
                            dim,
                            nLow
                        )

                    _, funR, _ =
                        gmsh.model.mesh.getBasisFunctions(
                            et,
                            ξLow,
                            "Lagrange"
                        )

                    Re =
                        Matrix(
                            transpose(
                                reshape(
                                    funR,
                                    nHigh,
                                    nLow
                                )
                            )
                        )

                    push!(nHighs, nHigh)
                    push!(nLows, nLow)
                    push!(nPrimarys, nPrimary)
                    push!(nPrimaryLows, nPrimaryLow)

                    push!(Tes, Te)
                    push!(Res, Re)

                    # Store coordinates as node × xyz:
                    #
                    # Xlow = Re * Xhigh
                    #
                    push!(
                        Xhigh_work,
                        Matrix{Float64}(undef, nHigh, 3)
                    )

                    push!(
                        Xlow_work,
                        Matrix{Float64}(undef, nLow, 3)
                    )

                    tid = length(nHighs)
                    type_id[et] = tid

                else

                    # Homogeneous-order check also for already known types
                    _, _, p, _, _, _ =
                        gmsh.model.mesh.getElementProperties(et)

                    p == p_global ||
                        error(
                            "reductionMatrices: the mesh must have homogeneous polynomial order."
                        )
                end

                # ------------------------------------------------------
                # Collect connectivity
                # ------------------------------------------------------

                nHigh = nHighs[tid]

                tags = elemTags[it]
                conn = elemNodeTags[it]

                @inbounds for e in eachindex(tags)

                    elemTag = tags[e]

                    elemTag in seen_elements && continue
                    push!(seen_elements, elemTag)

                    o = (e - 1) * nHigh

                    nodes = Vector{UInt64}(undef, nHigh)

                    copyto!(
                        nodes,
                        1,
                        conn,
                        o + 1,
                        nHigh
                    )

                    push!(elem_type_id, tid)
                    push!(elem_nodes, nodes)

                    # Bounding box at no additional mesh traversal cost
                    for nodeTag in nodes

                        node = Int(nodeTag)

                        x = nodeCoord[1, node]
                        y = nodeCoord[2, node]
                        z = nodeCoord[3, node]

                        xmin = min(xmin, x)
                        ymin = min(ymin, y)
                        zmin = min(zmin, z)

                        xmax = max(xmax, x)
                        ymax = max(ymax, y)
                        zmax = max(zmax, z)
                    end
                end
            end
        end
    end

    isempty(elem_nodes) &&
        error(
            "reductionMatrices: no domain elements found."
        )

    p_global > 1 ||
        error(
            "reductionMatrices: polynomial order must be greater than one."
        )

    # ------------------------------------------------------------------
    # Coordinate comparison tolerance
    # ------------------------------------------------------------------

    L =
        max(
            xmax - xmin,
            ymax - ymin,
            zmax - zmin
        )

    L > 0 ||
        error(
            "reductionMatrices: degenerate mesh geometry."
        )

    tol2 = (1e-10 * L)^2

    # ------------------------------------------------------------------
    # Global reduced-node numbering
    # ------------------------------------------------------------------

    nElem = length(elem_nodes)

    reduced_conn =
        Vector{Vector{Int}}(undef, nElem)

    # Global physical coordinates of reduced nodes.
    # A coordinate is stored only once per global reduced node.
    xr = Float64[]
    yr = Float64[]
    zr = Float64[]

    # Original primary node -> reduced node.
    # Direct indexing is possible because LLFEM already relies on
    # contiguous Gmsh node tags.
    vertex_to_reduced =
        zeros(Int, P.non)

    # Primary node -> already processed elements touching that node.
    #
    # Vectors are allocated lazily only for primary nodes actually used.
    node_to_elements =
        Vector{Union{Nothing,Vector{Int}}}(undef, P.non)

    fill!(node_to_elements, nothing)

    # Workspace for neighbour detection without allocating a Dict
    # for every element.
    neighbour_hits = zeros(Int, nElem)

    touched = Int[]
    sizehint!(touched, 32)

    nReduced = 0

    for ie in 1:nElem

        tid = elem_type_id[ie]

        nodes = elem_nodes[ie]

        nHigh = nHighs[tid]
        nLow = nLows[tid]
        nPrimary = nPrimarys[tid]
        nPrimaryLow = nPrimaryLows[tid]

        Re = Res[tid]

        Xhigh = Xhigh_work[tid]
        Xlow = Xlow_work[tid]

        # --------------------------------------------------------------
        # Physical coordinates of the high-order element nodes
        # --------------------------------------------------------------

        @inbounds for a in 1:nHigh

            node = Int(nodes[a])

            Xhigh[a, 1] = nodeCoord[1, node]
            Xhigh[a, 2] = nodeCoord[2, node]
            Xhigh[a, 3] = nodeCoord[3, node]
        end

        # --------------------------------------------------------------
        # Physical coordinates of reduced nodes
        #
        # Xlow = Re * Xhigh
        #
        # No temporary matrix is allocated.
        # --------------------------------------------------------------

        mul!(
            Xlow,
            Re,
            Xhigh
        )

        rconn = Vector{Int}(undef, nLow)

        # --------------------------------------------------------------
        # Find previously processed topological neighbours.
        #
        # A shared non-vertex C0 node can only lie on a common edge/face.
        # Such elements share at least two primary nodes.
        # --------------------------------------------------------------

        empty!(touched)

        @inbounds for a in 1:nPrimary

            node = Int(nodes[a])

            lst = node_to_elements[node]

            lst === nothing && continue

            for je in lst

                if neighbour_hits[je] == 0
                    push!(touched, je)
                end

                neighbour_hits[je] += 1
            end
        end

        # --------------------------------------------------------------
        # Number reduced nodes
        # --------------------------------------------------------------

        @inbounds for a in 1:nLow

            # ----------------------------------------------------------
            # Vertex nodes can be identified exactly by original
            # Gmsh node tag.
            # ----------------------------------------------------------

            if a <= nPrimaryLow

                node = Int(nodes[a])

                rnode = vertex_to_reduced[node]

                if rnode == 0

                    nReduced += 1
                    rnode = nReduced

                    vertex_to_reduced[node] = rnode

                    push!(xr, Xlow[a, 1])
                    push!(yr, Xlow[a, 2])
                    push!(zr, Xlow[a, 3])
                end

                rconn[a] = rnode

                continue
            end

            # ----------------------------------------------------------
            # Edge/face/internal reduced node
            # ----------------------------------------------------------

            xa = Xlow[a, 1]
            ya = Xlow[a, 2]
            za = Xlow[a, 3]

            found = 0

            for je in touched

                # Sharing only one primary node means vertex contact.
                neighbour_hits[je] >= 2 || continue

                old_conn = reduced_conn[je]

                for rnode in old_conn

                    dx = xa - xr[rnode]
                    dy = ya - yr[rnode]
                    dz = za - zr[rnode]

                    if dx * dx +
                       dy * dy +
                       dz * dz <= tol2

                        found = rnode
                        break
                    end
                end

                found != 0 && break
            end

            if found == 0

                nReduced += 1
                found = nReduced

                push!(xr, xa)
                push!(yr, ya)
                push!(zr, za)
            end

            rconn[a] = found
        end

        reduced_conn[ie] = rconn

        # --------------------------------------------------------------
        # Reset neighbour workspace
        # --------------------------------------------------------------

        @inbounds for je in touched
            neighbour_hits[je] = 0
        end

        # --------------------------------------------------------------
        # Register current element at its primary nodes
        # --------------------------------------------------------------

        @inbounds for a in 1:nPrimary

            node = Int(nodes[a])

            lst = node_to_elements[node]

            if lst === nothing

                new_list = Int[ie]
                node_to_elements[node] = new_list

            else

                push!(lst, ie)
            end
        end
    end

    # ------------------------------------------------------------------
    # Assemble scalar Tn
    # ------------------------------------------------------------------

    IT = Int[]
    JT = Int[]
    VT = Float64[]

    max_nLow = maximum(nLows)

    sizehint!(IT, P.non * max_nLow)
    sizehint!(JT, P.non * max_nLow)
    sizehint!(VT, P.non * max_nLow)

    full_seen = falses(P.non)

    zero_tol = 100 * eps(Float64)

    @inbounds for ie in 1:nElem

        tid = elem_type_id[ie]

        nodes = elem_nodes[ie]
        rconn = reduced_conn[ie]

        Te = Tes[tid]

        nHigh = nHighs[tid]
        nLow = nLows[tid]

        for a in 1:nHigh

            node = Int(nodes[a])

            full_seen[node] && continue
            full_seen[node] = true

            for b in 1:nLow

                v = Te[a, b]

                abs(v) <= zero_tol && continue

                push!(IT, node)
                push!(JT, rconn[b])
                push!(VT, v)
            end
        end
    end

    Tn =
        sparse(
            IT,
            JT,
            VT,
            P.non,
            nReduced
        )

    # ------------------------------------------------------------------
    # Assemble scalar Rn
    #
    # For a shared reduced node, the interpolation row from the first
    # processed element containing that node is used.
    # ------------------------------------------------------------------

    IR = Int[]
    JR = Int[]
    VR = Float64[]

    max_nHigh = maximum(nHighs)

    sizehint!(IR, nReduced * max_nHigh)
    sizehint!(JR, nReduced * max_nHigh)
    sizehint!(VR, nReduced * max_nHigh)

    reduced_seen = falses(nReduced)

    @inbounds for ie in 1:nElem

        tid = elem_type_id[ie]

        nodes = elem_nodes[ie]
        rconn = reduced_conn[ie]

        Re = Res[tid]

        nHigh = nHighs[tid]
        nLow = nLows[tid]

        for a in 1:nLow

            rnode = rconn[a]

            reduced_seen[rnode] && continue
            reduced_seen[rnode] = true

            for b in 1:nHigh

                v = Re[a, b]

                abs(v) <= zero_tol && continue

                push!(IR, rnode)
                push!(JR, Int(nodes[b]))
                push!(VR, v)
            end
        end
    end

    Rn =
        sparse(
            IR,
            JR,
            VR,
            nReduced,
            P.non
        )

    # ------------------------------------------------------------------
    # Expand according to P.pdim
    # ------------------------------------------------------------------

    d = P.pdim

    if d == 1
        return Tn, Rn
    end

    Id =
        spdiagm(
            0 => ones(Float64, d)
        )

    T = kron(Tn, Id)
    R = kron(Rn, Id)

    return T, R
end