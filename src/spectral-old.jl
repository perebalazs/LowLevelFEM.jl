###############################################################################
#                                                                             #
#                              spectral.jl                                    #
#                                                                             #
#  Experimental spectral-element basis transformations.                      #
#                                                                             #
#  Current scope:                                                            #
#      - 1D Line elements                                                     #
#      - 2D complete tensor-product Quadrangle elements                       #
#      - Legendre-Gauss-Lobatto nodal basis                                   #
#      - homogeneous polynomial order                                         #
#                                                                             #
#  The ordinary LowLevelFEM/Gmsh Lagrange space remains the global storage    #
#  space. Only the solver-side spectral space is compact.                     #
#                                                                             #
###############################################################################

using LinearAlgebra
using SparseArrays

export spectralTransformation


# =============================================================================
# Legendre-Gauss-Lobatto helpers
# =============================================================================

"""
    _sem_legendre_pair(p, x)

Return `(Pp, Ppm1)`, where `Pp = P_p(x)` and `Ppm1 = P_{p-1}(x)`.
"""
function _sem_legendre_pair(
    p::Integer,
    x::Real
)
    p >= 0 ||
        throw(
            ArgumentError(
                "Polynomial order must be non-negative."
            )
        )

    p == 0 &&
        return (1.0, 0.0)

    P0 = 1.0
    P1 = Float64(x)

    p == 1 &&
        return (P1, P0)

    for k in 2:Int(p)
        Pk =
            ((2k - 1) * x * P1 - (k - 1) * P0) / k

        P0, P1 =
            P1, Pk
    end

    return P1, P0
end


"""
    _sem_gll(p; tol=1e-14, maxiter=100)

Compute the `p+1` Legendre-Gauss-Lobatto nodes and quadrature weights on
`[-1, 1]`.
"""
function _sem_gll(
    p::Integer;
    tol::Real=1e-14,
    maxiter::Integer=100
)
    p >= 1 ||
        throw(
            ArgumentError(
                "Polynomial order must be at least one."
            )
        )

    p0 = Int(p)

    ξ =
        zeros(Float64, p0 + 1)

    ξ[1] = -1.0
    ξ[end] = 1.0

    for k in 1:p0-1
        x =
            -cos(pi * k / p0)

        converged = false

        for _ in 1:Int(maxiter)
            Pp, Ppm1 =
                _sem_legendre_pair(
                    p0,
                    x
                )

            dPp =
                p0 *
                (Ppm1 - x * Pp) /
                (1 - x^2)

            ddPp =
                (
                    2x * dPp -
                    p0 * (p0 + 1) * Pp
                ) /
                (1 - x^2)

            Δx =
                dPp / ddPp

            x -= Δx

            if abs(Δx) <=
               tol * max(1.0, abs(x))
                converged = true
                break
            end
        end

        converged ||
            error(
                "_sem_gll: Newton iteration did not converge " *
                "for polynomial order $p0."
            )

        ξ[k + 1] = x
    end

    sort!(ξ)

    w =
        similar(ξ)

    for i in eachindex(ξ)
        Pp, _ =
            _sem_legendre_pair(
                p0,
                ξ[i]
            )

        w[i] =
            2 /
            (
                p0 *
                (p0 + 1) *
                Pp^2
            )
    end

    return ξ, w
end


"""
    _sem_barycentric_weights(nodes)

Compute barycentric interpolation weights for distinct interpolation nodes.
"""
function _sem_barycentric_weights(
    nodes::AbstractVector
)
    n =
        length(nodes)

    w =
        ones(Float64, n)

    for j in 1:n
        for k in 1:n
            j == k &&
                continue

            w[j] /=
                nodes[j] - nodes[k]
        end
    end

    return w
end


"""
    _sem_lagrange_matrix(nodes, points; atol=1e-13)

Return the Lagrange interpolation matrix that maps nodal values at `nodes`
to values at `points`.
"""
function _sem_lagrange_matrix(
    nodes::AbstractVector,
    points::AbstractVector;
    atol::Real=1e-13
)
    xnodes =
        Float64.(nodes)

    xpoints =
        Float64.(points)

    n =
        length(xnodes)

    m =
        length(xpoints)

    bw =
        _sem_barycentric_weights(
            xnodes
        )

    N =
        zeros(Float64, m, n)

    for q in 1:m
        x =
            xpoints[q]

        hit =
            findfirst(
                j ->
                    isapprox(
                        x,
                        xnodes[j];
                        atol=atol,
                        rtol=0.0
                    ),
                1:n
            )

        if hit !== nothing
            N[q, hit] = 1.0
        else
            a =
                bw ./ (x .- xnodes)

            N[q, :] .=
                a ./ sum(a)
        end
    end

    return N
end


# =============================================================================
# Sparse assembly helpers
# =============================================================================

"""
    _sem_insert_unique!(entries, i, j, value; atol=1e-12)

Insert one sparse matrix entry. If an entry has already been generated by a
neighbouring element, verify that both values agree instead of summing them.
"""
function _sem_insert_unique!(
    entries::Dict{Tuple{Int,Int},Float64},
    i::Integer,
    j::Integer,
    value::Real;
    atol::Real=1e-12
)
    v =
        Float64(value)

    abs(v) <=
        100eps(Float64) &&
        return nothing

    key =
        (Int(i), Int(j))

    if haskey(entries, key)
        isapprox(
            entries[key],
            v;
            atol=atol,
            rtol=atol
        ) ||
            error(
                "spectralTransformation: inconsistent transformation " *
                "entry ($(key[1]), $(key[2])): " *
                "$(entries[key]) vs. $v."
            )
    else
        entries[key] = v
    end

    return nothing
end


"""
    _sem_sparse_from_dict(entries, m, n)

Construct a sparse matrix from dictionary entries.
"""
function _sem_sparse_from_dict(
    entries::Dict{Tuple{Int,Int},Float64},
    m::Integer,
    n::Integer
)
    I = Int[]
    J = Int[]
    V = Float64[]

    sizehint!(I, length(entries))
    sizehint!(J, length(entries))
    sizehint!(V, length(entries))

    for ((i, j), value) in entries
        push!(I, i)
        push!(J, j)
        push!(V, value)
    end

    return sparse(
        I,
        J,
        V,
        Int(m),
        Int(n)
    )
end


# =============================================================================
# Current LowLevelFEM global numbering adapter
# =============================================================================

"""
    _sem_check_global_numbering(P, nodeTags)

Verify the global Gmsh numbering convention currently used by LowLevelFEM.

At present LowLevelFEM uses Gmsh node tags directly in the global algebraic
numbering. Therefore the full model must use contiguous node tags `1:P.non`.
If Gmsh leaves gaps in the numbering, call `renumberNodes!` before creating
the `Problem`.

This check is intentionally isolated so that a future compact global node
numbering can be introduced without changing the spectral-element logic.
"""
function _sem_check_global_numbering(
    P::Problem,
    nodeTags
)
    tags =
        sort!(
            Int.(collect(nodeTags))
        )

    expected =
        collect(1:P.non)

    tags == expected ||
        error(
            "spectralTransformation: the current LowLevelFEM global " *
            "numbering requires contiguous Gmsh node tags 1:P.non. " *
            "Call renumberNodes! before creating the Problem."
        )

    return nothing
end


"""
    _sem_global_node(P, nodeTag)

Map a Gmsh node tag to the current LowLevelFEM global node index.

Today this is the identity map after `_sem_check_global_numbering`. Keeping
the mapping behind a helper makes the SEM implementation ready for a future
core numbering refactor.
"""
@inline function _sem_global_node(
    P::Problem,
    nodeTag::Integer
)
    return Int(nodeTag)
end


# =============================================================================
# 1D element collection and spectral node numbering
# =============================================================================

"""
    _sem_collect_1d_elements(P)

Collect all one-dimensional domain elements belonging to `P` and cache the
local FEM-to-GLL basis transformations for each Gmsh element type.
"""
function _sem_collect_1d_elements(
    P::Problem
)
    gmsh.model.setCurrent(P.name)

    P.dim == 1 ||
        error(
            "spectralTransformation: only one-dimensional Problems " *
            "are currently supported."
        )

    nodeTags, coord, _ =
        gmsh.model.mesh.getNodes(
            -1,
            -1,
            false,
            false
        )

    _sem_check_global_numbering(
        P,
        nodeTags
    )

    nodeCoord =
        Dict{
            UInt64,
            NTuple{3,Float64}
        }()

    sizehint!(
        nodeCoord,
        length(nodeTags)
    )

    @inbounds for i in eachindex(nodeTags)
        nodeCoord[nodeTags[i]] = (
            coord[3i - 2],
            coord[3i - 1],
            coord[3i]
        )
    end

    elements =
        NamedTuple[]

    seen_elements =
        Set{UInt64}()

    orders =
        Set{Int}()

    type_cache =
        Dict{Int,Any}()

    for mat in P.material
        dimTags =
            gmsh.model.getEntitiesForPhysicalName(
                mat.phName
            )

        for (edim, etag) in dimTags
            edim == P.dim ||
                continue

            elemTypes,
            elemTags,
            elemNodeTags =
                gmsh.model.mesh.getElements(
                    edim,
                    etag
                )

            for it in eachindex(elemTypes)
                et =
                    Int(elemTypes[it])

                name,
                dim,
                p0,
                nFem0,
                ξFem0,
                nPrimary0 =
                    gmsh.model.mesh.getElementProperties(
                        et
                    )

                dim == 1 ||
                    error(
                        "spectralTransformation: element \"$name\" " *
                        "is not one-dimensional."
                    )

                occursin("Line", name) ||
                    error(
                        "spectralTransformation: unsupported 1D element " *
                        "family \"$name\"."
                    )

                p =
                    Int(p0)

                nFem =
                    Int(nFem0)

                nPrimary =
                    Int(nPrimary0)

                push!(
                    orders,
                    p
                )

                if !haskey(type_cache, et)
                    ξFem =
                        Float64.(ξFem0)

                    length(ξFem) == nFem ||
                        error(
                            "spectralTransformation: invalid local node " *
                            "coordinates for element \"$name\"."
                        )

                    ξGLL, wGLL =
                        _sem_gll(p)

                    nSem =
                        length(ξGLL)

                    nSem == nFem ||
                        error(
                            "spectralTransformation: FEM and spectral " *
                            "spaces must have the same local dimension."
                        )

                    Te =
                        _sem_lagrange_matrix(
                            ξGLL,
                            ξFem
                        )

                    Re =
                        _sem_lagrange_matrix(
                            ξFem,
                            ξGLL
                        )

                    left_fem =
                        argmin(ξFem)

                    right_fem =
                        argmax(ξFem)

                    isapprox(
                        ξFem[left_fem],
                        -1.0;
                        atol=1e-12,
                        rtol=0.0
                    ) ||
                        error(
                            "spectralTransformation: could not identify " *
                            "the left endpoint of element \"$name\"."
                        )

                    isapprox(
                        ξFem[right_fem],
                        1.0;
                        atol=1e-12,
                        rtol=0.0
                    ) ||
                        error(
                            "spectralTransformation: could not identify " *
                            "the right endpoint of element \"$name\"."
                        )

                    type_cache[et] = (
                        name=name,
                        p=p,
                        nFem=nFem,
                        nSem=nSem,
                        nPrimary=nPrimary,
                        ξFem=ξFem,
                        ξGLL=ξGLL,
                        wGLL=wGLL,
                        Te=Te,
                        Re=Re,
                        left_fem=left_fem,
                        right_fem=right_fem
                    )
                end

                cache =
                    type_cache[et]

                tags =
                    elemTags[it]

                conn =
                    elemNodeTags[it]

                @inbounds for e in eachindex(tags)
                    elemTag =
                        UInt64(tags[e])

                    elemTag in seen_elements &&
                        continue

                    push!(
                        seen_elements,
                        elemTag
                    )

                    offset =
                        (e - 1) * cache.nFem

                    nodes =
                        UInt64.(
                            conn[
                                offset + 1 : offset + cache.nFem
                            ]
                        )

                    push!(
                        elements,
                        (
                            tag=elemTag,
                            et=et,
                            nodes=nodes
                        )
                    )
                end
            end
        end
    end

    isempty(elements) &&
        error(
            "spectralTransformation: no 1D domain elements found."
        )

    length(orders) == 1 ||
        error(
            "spectralTransformation: homogeneous polynomial order " *
            "is currently required; found $(sort!(collect(orders)))."
        )

    return (
        elements=elements,
        type_cache=type_cache,
        nodeCoord=nodeCoord,
        order=first(orders)
    )
end


"""
    _sem_number_spectral_nodes_1d(elements, type_cache, nodeCoord)

Construct compact global spectral-node numbering for a conforming 1D mesh.

Element endpoints are shared through their original Gmsh endpoint tags.
Interior GLL nodes are element-local. The returned physical coordinates are
computed from the original isoparametric Gmsh geometry.
"""
function _sem_number_spectral_nodes_1d(
    elements,
    type_cache,
    nodeCoord
)
    nElem =
        length(elements)

    spectral_conn =
        Vector{Vector{Int}}(
            undef,
            nElem
        )

    spectral_coordinates =
        NTuple{3,Float64}[]

    endpoint_to_spectral =
        Dict{UInt64,Int}()

    nSpectral =
        0

    for ie in 1:nElem
        elem =
            elements[ie]

        cache =
            type_cache[elem.et]

        Xfem =
            zeros(
                Float64,
                cache.nFem,
                3
            )

        for a in 1:cache.nFem
            x, y, z =
                nodeCoord[
                    elem.nodes[a]
                ]

            Xfem[a, 1] = x
            Xfem[a, 2] = y
            Xfem[a, 3] = z
        end

        Xgll =
            cache.Re * Xfem

        sconn =
            zeros(
                Int,
                cache.nSem
            )

        for a in 1:cache.nSem
            if a == 1 || a == cache.nSem
                fem_local =
                    a == 1 ?
                    cache.left_fem :
                    cache.right_fem

                endpoint_tag =
                    elem.nodes[fem_local]

                spectral_node =
                    get(
                        endpoint_to_spectral,
                        endpoint_tag,
                        0
                    )

                if spectral_node == 0
                    nSpectral += 1
                    spectral_node = nSpectral

                    endpoint_to_spectral[
                        endpoint_tag
                    ] = spectral_node

                    push!(
                        spectral_coordinates,
                        (
                            Xgll[a, 1],
                            Xgll[a, 2],
                            Xgll[a, 3]
                        )
                    )
                end

                sconn[a] =
                    spectral_node
            else
                nSpectral += 1

                sconn[a] =
                    nSpectral

                push!(
                    spectral_coordinates,
                    (
                        Xgll[a, 1],
                        Xgll[a, 2],
                        Xgll[a, 3]
                    )
                )
            end
        end

        spectral_conn[ie] =
            sconn
    end

    Xspectral =
        zeros(
            Float64,
            nSpectral,
            3
        )

    for i in 1:nSpectral
        Xspectral[i, :] .=
            spectral_coordinates[i]
    end

    return (
        connectivity=spectral_conn,
        coordinates=Xspectral,
        count=nSpectral
    )
end



# =============================================================================
# 2D quadrilateral element collection and spectral node numbering
# =============================================================================

"""
    _sem_local_to_3d(localCoord, dim, n)

Convert Gmsh reference coordinates returned by `getElementProperties` to the
`(u,v,w)` triplet layout expected by `getBasisFunctions`.
"""
function _sem_local_to_3d(
    localCoord,
    dim::Integer,
    n::Integer
)
    ξ =
        zeros(
            Float64,
            3 * Int(n)
        )

    @inbounds for a in 1:Int(n)
        for d in 1:Int(dim)
            ξ[3 * (a - 1) + d] =
                localCoord[
                    Int(dim) * (a - 1) + d
                ]
        end
    end

    return ξ
end


"""
    _sem_find_local_node_2d(ξ, η, ξ0, η0, name; atol=1e-12)

Return the local FEM node index at the requested reference coordinate.
"""
function _sem_find_local_node_2d(
    ξ,
    η,
    ξ0::Real,
    η0::Real,
    name::AbstractString;
    atol::Real=1e-12
)
    for a in eachindex(ξ)
        if isapprox(
               ξ[a],
               ξ0;
               atol=atol,
               rtol=0.0
           ) &&
           isapprox(
               η[a],
               η0;
               atol=atol,
               rtol=0.0
           )
            return a
        end
    end

    error(
        "spectralTransformation: could not identify reference point " *
        "($ξ0, $η0) in element \"$name\"."
    )
end


"""
    _sem_edge_key(a, b)

Return an orientation-independent key for a topological edge.
"""
@inline function _sem_edge_key(
    a::UInt64,
    b::UInt64
)
    return a < b ?
           (a, b) :
           (b, a)
end


"""
    _sem_distance2(a, b)

Squared Euclidean distance between two 3D coordinate tuples.
"""
@inline function _sem_distance2(
    a::NTuple{3,Float64},
    b::NTuple{3,Float64}
)
    dx = a[1] - b[1]
    dy = a[2] - b[2]
    dz = a[3] - b[3]

    return dx * dx +
           dy * dy +
           dz * dz
end


"""
    _sem_collect_2d_quad_elements(P)

Collect all two-dimensional quadrilateral domain elements belonging to `P`
and cache the local FEM-to-GLL basis transformations.

Only complete tensor-product Gmsh quadrangles are currently supported.
Triangular spectral elements are deliberately rejected with a
`not yet implemented` error.
"""
function _sem_collect_2d_quad_elements(
    P::Problem
)
    gmsh.model.setCurrent(P.name)

    P.dim == 2 ||
        error(
            "spectralTransformation: _sem_collect_2d_quad_elements " *
            "requires a two-dimensional Problem."
        )

    nodeTags, coord, _ =
        gmsh.model.mesh.getNodes(
            -1,
            -1,
            false,
            false
        )

    _sem_check_global_numbering(
        P,
        nodeTags
    )

    nodeCoord =
        Dict{
            UInt64,
            NTuple{3,Float64}
        }()

    sizehint!(
        nodeCoord,
        length(nodeTags)
    )

    @inbounds for i in eachindex(nodeTags)
        nodeCoord[nodeTags[i]] = (
            coord[3i - 2],
            coord[3i - 1],
            coord[3i]
        )
    end

    elements =
        NamedTuple[]

    seen_elements =
        Set{UInt64}()

    orders =
        Set{Int}()

    type_cache =
        Dict{Int,Any}()

    for mat in P.material
        dimTags =
            gmsh.model.getEntitiesForPhysicalName(
                mat.phName
            )

        for (edim, etag) in dimTags
            edim == P.dim ||
                continue

            elemTypes,
            elemTags,
            elemNodeTags =
                gmsh.model.mesh.getElements(
                    edim,
                    etag
                )

            for it in eachindex(elemTypes)
                et =
                    Int(elemTypes[it])

                name,
                dim,
                p0,
                nFem0,
                localFem0,
                nPrimary0 =
                    gmsh.model.mesh.getElementProperties(
                        et
                    )

                dim == 2 ||
                    error(
                        "spectralTransformation: element \"$name\" " *
                        "is not two-dimensional."
                    )

                if occursin(
                    "Triangle",
                    name
                )
                    error(
                        "spectralTransformation: triangular spectral " *
                        "elements are not yet implemented."
                    )
                end

                (
                    occursin(
                        "Quadrangle",
                        name
                    ) ||
                    occursin(
                        "Quadrilateral",
                        name
                    )
                ) ||
                    error(
                        "spectralTransformation: unsupported 2D element " *
                        "family \"$name\". Only quadrilateral elements " *
                        "are currently supported."
                    )

                p =
                    Int(p0)

                nFem =
                    Int(nFem0)

                nPrimary =
                    Int(nPrimary0)

                push!(
                    orders,
                    p
                )

                if !haskey(type_cache, et)
                    n1 =
                        p + 1

                    nSem =
                        n1^2

                    nFem == nSem ||
                        error(
                            "spectralTransformation: 2D spectral elements " *
                            "currently require complete tensor-product " *
                            "quadrangles with (p+1)^2 nodes. Element " *
                            "\"$name\" has $nFem nodes for p=$p."
                        )

                    nPrimary == 4 ||
                        error(
                            "spectralTransformation: quadrilateral element " *
                            "\"$name\" must have four primary vertices."
                        )

                    length(localFem0) ==
                    2 * nFem ||
                        error(
                            "spectralTransformation: invalid local node " *
                            "coordinates for element \"$name\"."
                        )

                    ξFem =
                        Vector{Float64}(
                            undef,
                            nFem
                        )

                    ηFem =
                        Vector{Float64}(
                            undef,
                            nFem
                        )

                    @inbounds for a in 1:nFem
                        ξFem[a] =
                            localFem0[
                                2a - 1
                            ]

                        ηFem[a] =
                            localFem0[
                                2a
                            ]
                    end

                    ξGLL, wGLL =
                        _sem_gll(p)

                    # --------------------------------------------------
                    # Local GLL ordering:
                    #
                    #   b = (j - 1) * (p + 1) + i
                    #
                    # with ξ varying fastest.
                    # --------------------------------------------------

                    localGLL =
                        zeros(
                            Float64,
                            nSem,
                            2
                        )

                    b = 0

                    for j in 1:n1
                        for i in 1:n1
                            b += 1

                            localGLL[b, 1] =
                                ξGLL[i]

                            localGLL[b, 2] =
                                ξGLL[j]
                        end
                    end

                    # --------------------------------------------------
                    # Local prolongation:
                    #
                    #     u_FEM = Te * u_GLL
                    #
                    # Evaluate the tensor-product GLL basis at the
                    # original Gmsh FEM nodes.
                    # --------------------------------------------------

                    Lξ =
                        _sem_lagrange_matrix(
                            ξGLL,
                            ξFem
                        )

                    Lη =
                        _sem_lagrange_matrix(
                            ξGLL,
                            ηFem
                        )

                    Te =
                        zeros(
                            Float64,
                            nFem,
                            nSem
                        )

                    @inbounds for a in 1:nFem
                        b = 0

                        for j in 1:n1
                            for i in 1:n1
                                b += 1

                                Te[a, b] =
                                    Lξ[a, i] *
                                    Lη[a, j]
                            end
                        end
                    end

                    # --------------------------------------------------
                    # Local restriction:
                    #
                    #     u_GLL = Re * u_FEM
                    #
                    # Evaluate the original Gmsh Lagrange basis at the
                    # tensor-product GLL points. This keeps the
                    # implementation independent of Gmsh local node
                    # ordering.
                    # --------------------------------------------------

                    localGLL3 =
                        zeros(
                            Float64,
                            3 * nSem
                        )

                    @inbounds for a in 1:nSem
                        localGLL3[
                            3a - 2
                        ] =
                            localGLL[a, 1]

                        localGLL3[
                            3a - 1
                        ] =
                            localGLL[a, 2]
                    end

                    ncomp,
                    funR,
                    norient =
                        gmsh.model.mesh.getBasisFunctions(
                            et,
                            localGLL3,
                            "Lagrange"
                        )

                    ncomp == 1 ||
                        error(
                            "spectralTransformation: expected scalar " *
                            "Lagrange basis functions for element \"$name\"."
                        )

                    Re =
                        Matrix(
                            transpose(
                                reshape(
                                    funR,
                                    nFem,
                                    nSem
                                )
                            )
                        )

                    # Local consistency check before global assembly.
                    local_inverse_error =
                        norm(
                            Re * Te -
                            Matrix{Float64}(
                                I,
                                nSem,
                                nSem
                            ),
                            Inf
                        )

                    local_inverse_error <= 1e-10 ||
                        error(
                            "spectralTransformation: local 2D FEM/GLL " *
                            "transformations are inconsistent for element " *
                            "\"$name\"; ||Re*Te-I||∞ = " *
                            "$local_inverse_error."
                        )

                    corner_fem = (
                        mm=
                            _sem_find_local_node_2d(
                                ξFem,
                                ηFem,
                                -1.0,
                                -1.0,
                                name
                            ),
                        pm=
                            _sem_find_local_node_2d(
                                ξFem,
                                ηFem,
                                1.0,
                                -1.0,
                                name
                            ),
                        pp=
                            _sem_find_local_node_2d(
                                ξFem,
                                ηFem,
                                1.0,
                                1.0,
                                name
                            ),
                        mp=
                            _sem_find_local_node_2d(
                                ξFem,
                                ηFem,
                                -1.0,
                                1.0,
                                name
                            )
                    )

                    type_cache[et] = (
                        name=name,
                        p=p,
                        n1=n1,
                        nFem=nFem,
                        nSem=nSem,
                        nPrimary=nPrimary,
                        ξFem=ξFem,
                        ηFem=ηFem,
                        ξGLL=ξGLL,
                        wGLL=wGLL,
                        localGLL=localGLL,
                        Te=Te,
                        Re=Re,
                        corner_fem=corner_fem
                    )
                end

                cache =
                    type_cache[et]

                tags =
                    elemTags[it]

                conn =
                    elemNodeTags[it]

                @inbounds for e in eachindex(tags)
                    elemTag =
                        UInt64(tags[e])

                    elemTag in seen_elements &&
                        continue

                    push!(
                        seen_elements,
                        elemTag
                    )

                    offset =
                        (e - 1) *
                        cache.nFem

                    nodes =
                        UInt64.(
                            conn[
                                offset + 1 : offset + cache.nFem
                            ]
                        )

                    push!(
                        elements,
                        (
                            tag=elemTag,
                            et=et,
                            nodes=nodes
                        )
                    )
                end
            end
        end
    end

    isempty(elements) &&
        error(
            "spectralTransformation: no 2D quadrilateral domain " *
            "elements found."
        )

    length(orders) == 1 ||
        error(
            "spectralTransformation: homogeneous polynomial order " *
            "is currently required; found $(sort!(collect(orders)))."
        )

    return (
        elements=elements,
        type_cache=type_cache,
        nodeCoord=nodeCoord,
        order=first(orders)
    )
end


"""
    _sem_number_spectral_nodes_2d(elements, type_cache, nodeCoord)

Construct compact global spectral-node numbering for a conforming
quadrilateral mesh.

Vertices are shared by their original Gmsh corner tags. Spectral nodes on a
shared edge are identified by the orientation-independent pair of its corner
tags and by physical position. Element-interior GLL nodes remain local to the
element.

The physical GLL coordinates are obtained by isoparametric interpolation of
the original Gmsh geometry.
"""
function _sem_number_spectral_nodes_2d(
    elements,
    type_cache,
    nodeCoord
)
    nElem =
        length(elements)

    spectral_conn =
        Vector{Vector{Int}}(
            undef,
            nElem
        )

    spectral_coordinates =
        NTuple{3,Float64}[]

    vertex_to_spectral =
        Dict{
            UInt64,
            Int
        }()

    edge_to_spectral =
        Dict{
            Tuple{UInt64,UInt64},
            Vector{Int}
        }()

    # Characteristic length for coordinate matching on shared edges.
    used_nodes =
        Set{UInt64}()

    for elem in elements
        union!(
            used_nodes,
            elem.nodes
        )
    end

    xmin = Inf
    ymin = Inf
    zmin = Inf
    xmax = -Inf
    ymax = -Inf
    zmax = -Inf

    for node in used_nodes
        x, y, z =
            nodeCoord[node]

        xmin = min(xmin, x)
        ymin = min(ymin, y)
        zmin = min(zmin, z)

        xmax = max(xmax, x)
        ymax = max(ymax, y)
        zmax = max(zmax, z)
    end

    L =
        max(
            xmax - xmin,
            ymax - ymin,
            zmax - zmin
        )

    L > 0 ||
        error(
            "spectralTransformation: degenerate 2D mesh geometry."
        )

    tol =
        1e-10 * L

    tol2 =
        tol^2

    nSpectral =
        0

    for ie in 1:nElem
        elem =
            elements[ie]

        cache =
            type_cache[elem.et]

        Xfem =
            zeros(
                Float64,
                cache.nFem,
                3
            )

        for a in 1:cache.nFem
            x, y, z =
                nodeCoord[
                    elem.nodes[a]
                ]

            Xfem[a, 1] = x
            Xfem[a, 2] = y
            Xfem[a, 3] = z
        end

        Xgll =
            cache.Re *
            Xfem

        sconn =
            zeros(
                Int,
                cache.nSem
            )

        corner_tags = (
            mm=
                elem.nodes[
                    cache.corner_fem.mm
                ],
            pm=
                elem.nodes[
                    cache.corner_fem.pm
                ],
            pp=
                elem.nodes[
                    cache.corner_fem.pp
                ],
            mp=
                elem.nodes[
                    cache.corner_fem.mp
                ]
        )

        for j in 1:cache.n1
            for i in 1:cache.n1
                a =
                    (j - 1) *
                    cache.n1 +
                    i

                xgll = (
                    Xgll[a, 1],
                    Xgll[a, 2],
                    Xgll[a, 3]
                )

                on_left =
                    i == 1

                on_right =
                    i == cache.n1

                on_bottom =
                    j == 1

                on_top =
                    j == cache.n1

                spectral_node =
                    0

                # --------------------------------------------------
                # Corner
                # --------------------------------------------------

                if (on_left || on_right) &&
                   (on_bottom || on_top)

                    corner_tag =
                        if on_left &&
                           on_bottom
                            corner_tags.mm
                        elseif on_right &&
                               on_bottom
                            corner_tags.pm
                        elseif on_right &&
                               on_top
                            corner_tags.pp
                        else
                            corner_tags.mp
                        end

                    spectral_node =
                        get(
                            vertex_to_spectral,
                            corner_tag,
                            0
                        )

                    if spectral_node == 0
                        nSpectral += 1

                        spectral_node =
                            nSpectral

                        vertex_to_spectral[
                            corner_tag
                        ] =
                            spectral_node

                        push!(
                            spectral_coordinates,
                            xgll
                        )
                    end

                # --------------------------------------------------
                # Edge interior
                # --------------------------------------------------

                elseif on_left ||
                       on_right ||
                       on_bottom ||
                       on_top

                    edge_key =
                        if on_bottom
                            _sem_edge_key(
                                corner_tags.mm,
                                corner_tags.pm
                            )
                        elseif on_right
                            _sem_edge_key(
                                corner_tags.pm,
                                corner_tags.pp
                            )
                        elseif on_top
                            _sem_edge_key(
                                corner_tags.mp,
                                corner_tags.pp
                            )
                        else
                            _sem_edge_key(
                                corner_tags.mm,
                                corner_tags.mp
                            )
                        end

                    candidates =
                        get!(
                            edge_to_spectral,
                            edge_key,
                            Int[]
                        )

                    for candidate in candidates
                        if _sem_distance2(
                               spectral_coordinates[
                                   candidate
                               ],
                               xgll
                           ) <= tol2

                            spectral_node =
                                candidate

                            break
                        end
                    end

                    if spectral_node == 0
                        nSpectral += 1

                        spectral_node =
                            nSpectral

                        push!(
                            spectral_coordinates,
                            xgll
                        )

                        push!(
                            candidates,
                            spectral_node
                        )
                    end

                # --------------------------------------------------
                # Element interior
                # --------------------------------------------------

                else
                    nSpectral += 1

                    spectral_node =
                        nSpectral

                    push!(
                        spectral_coordinates,
                        xgll
                    )
                end

                sconn[a] =
                    spectral_node
            end
        end

        spectral_conn[ie] =
            sconn
    end

    Xspectral =
        zeros(
            Float64,
            nSpectral,
            3
        )

    for i in 1:nSpectral
        Xspectral[i, :] .=
            spectral_coordinates[i]
    end

    return (
        connectivity=spectral_conn,
        coordinates=Xspectral,
        count=nSpectral
    )
end


# =============================================================================
# Public transformation
# =============================================================================

"""
    spectralTransformation(P::Problem; details=false, check=false)

Construct a Legendre-Gauss-Lobatto spectral basis transformation for the
finite-element field `P`.

Supported element families are currently 1D `Line` elements and complete
tensor-product 2D `Quadrangle` elements.

The ordinary LowLevelFEM/Gmsh nodal representation remains the global storage
space. The spectral space is compact and contains only the DOFs belonging to
`P`.

The returned matrices satisfy

    u_global = T * u_spectral
    u_spectral = R * u_global

for fields represented by the spectral space.

Consequently, an already assembled LowLevelFEM system can be transformed as

    K_spectral = T' * K * T
    f_spectral = T' * f

and the solution is reconstructed in the ordinary global LowLevelFEM
representation by

    u_global = T * u_spectral

# Current limitations

- 1D `Line` and complete tensor-product 2D `Quadrangle` elements are supported,
- triangular spectral elements are not yet implemented,
- all domain elements must have the same polynomial order,
- the current LowLevelFEM global numbering requires contiguous Gmsh node tags
  `1:P.non`; use `renumberNodes!` if necessary,
- spectral basis combined with reduced-order interpolation is not intended,
- MPC compatibility is not yet validated.

# Keyword arguments

- `details=false`: return only `(T, R)`.
- `details=true`: return a named tuple with the transformation matrices and
  diagnostic spectral mesh information.
- `check=false`: if `true`, verify `R*T ≈ I` and that the represented global
  rows coincide with `allDoFs(P)`.

The implementation deliberately isolates the current Gmsh-tag-to-global-node
mapping so that a future LowLevelFEM compact-numbering refactor does not
require redesigning the spectral transformation.
"""
function spectralTransformation(
    P::Problem;
    details::Bool=false,
    check::Bool=false
)
    data,
    numbered =
        if P.dim == 1
            data1 =
                _sem_collect_1d_elements(P)

            numbered1 =
                _sem_number_spectral_nodes_1d(
                    data1.elements,
                    data1.type_cache,
                    data1.nodeCoord
                )

            data1,
            numbered1

        elseif P.dim == 2
            data2 =
                _sem_collect_2d_quad_elements(P)

            numbered2 =
                _sem_number_spectral_nodes_2d(
                    data2.elements,
                    data2.type_cache,
                    data2.nodeCoord
                )

            data2,
            numbered2

        else
            error(
                "spectralTransformation: dimension $(P.dim) is not yet " *
                "implemented. Supported dimensions are 1 and 2."
            )
        end

    elements =
        data.elements

    type_cache =
        data.type_cache

    spectral_conn =
        numbered.connectivity

    nSpectral =
        numbered.count

    d =
        P.pdim

    nFullDofs =
        ndofs(P)

    nSpectralDofs =
        nSpectral * d

    # ------------------------------------------------------------------
    # Prolongation
    #
    #     u_global = T * u_spectral
    #
    # Rows use the existing global LowLevelFEM DOF numbering.
    # Columns use compact spectral DOF numbering.
    # ------------------------------------------------------------------

    Tentries =
        Dict{
            Tuple{Int,Int},
            Float64
        }()

    for ie in eachindex(elements)
        elem =
            elements[ie]

        cache =
            type_cache[elem.et]

        sconn =
            spectral_conn[ie]

        for a in 1:cache.nFem
            global_node =
                _sem_global_node(
                    P,
                    elem.nodes[a]
                )

            for comp in 1:d
                row =
                    nodeToDof(
                        P,
                        global_node,
                        comp
                    )

                for b in 1:cache.nSem
                    col =
                        (sconn[b] - 1) *
                        d +
                        comp

                    _sem_insert_unique!(
                        Tentries,
                        row,
                        col,
                        cache.Te[a, b]
                    )
                end
            end
        end
    end

    T =
        _sem_sparse_from_dict(
            Tentries,
            nFullDofs,
            nSpectralDofs
        )

    dropzeros!(T)

    # ------------------------------------------------------------------
    # Restriction
    #
    #     u_spectral = R * u_global
    #
    # A shared spectral endpoint is taken from the first processed
    # element containing it.
    # ------------------------------------------------------------------

    Rentries =
        Dict{
            Tuple{Int,Int},
            Float64
        }()

    spectral_seen =
        falses(
            nSpectral
        )

    for ie in eachindex(elements)
        elem =
            elements[ie]

        cache =
            type_cache[elem.et]

        sconn =
            spectral_conn[ie]

        for a in 1:cache.nSem
            spectral_node =
                sconn[a]

            spectral_seen[
                spectral_node
            ] &&
                continue

            spectral_seen[
                spectral_node
            ] = true

            for comp in 1:d
                row =
                    (spectral_node - 1) *
                    d +
                    comp

                for b in 1:cache.nFem
                    global_node =
                        _sem_global_node(
                            P,
                            elem.nodes[b]
                        )

                    col =
                        nodeToDof(
                            P,
                            global_node,
                            comp
                        )

                    _sem_insert_unique!(
                        Rentries,
                        row,
                        col,
                        cache.Re[a, b]
                    )
                end
            end
        end
    end

    R =
        _sem_sparse_from_dict(
            Rentries,
            nSpectralDofs,
            nFullDofs
        )

    dropzeros!(R)

    # ------------------------------------------------------------------
    # Active-space consistency
    # ------------------------------------------------------------------

    active =
        sort!(
            unique(
                Int.(allDoFs(P))
            )
        )

    rowsT, _, _ =
        findnz(T)

    represented =
        sort!(
            unique(rowsT)
        )

    represented == active ||
        error(
            "spectralTransformation: represented global DOFs do not " *
            "match allDoFs(P)."
        )

    inverse_error =
        nothing

    if check
        Isem =
            spdiagm(
                0 =>
                ones(
                    Float64,
                    nSpectralDofs
                )
            )

        inverse_error =
            norm(
                R * T - Isem,
                Inf
            )

        inverse_error <= 1e-10 ||
            error(
                "spectralTransformation: R*T differs from identity; " *
                "||R*T-I||∞ = $inverse_error."
            )
    end

    if !details
        return T, R
    end

    first_cache =
        type_cache[
            elements[1].et
        ]

    return (
        T=T,
        R=R,
        dimension=P.dim,
        order=data.order,
        ξ_gll=first_cache.ξGLL,
        w_gll=first_cache.wGLL,
        local_gll_coordinates=
            P.dim == 1 ?
            reshape(
                copy(first_cache.ξGLL),
                :,
                1
            ) :
            copy(first_cache.localGLL),
        spectral_coordinates=
            numbered.coordinates,
        spectral_connectivity=
            spectral_conn,
        element_tags=
            [elem.tag for elem in elements],
        active_dofs=
            active,
        global_dofs=
            nFullDofs,
        spectral_dofs=
            nSpectralDofs,
        inverse_error=
            inverse_error
    )
end
