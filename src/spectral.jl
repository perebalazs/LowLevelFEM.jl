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
#      - experimental GLL quadrature for the bilinear form P ⋅ P              #
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
# Direct CSC transformation assembly helpers
# =============================================================================

const _SEM_ENTRY_TOL = 100eps(Float64)

@inline _sem_keep_entry(v::Real) = abs(v) > _SEM_ENTRY_TOL

"""
    _sem_finalize_column_pattern!(columns)

Sort and deduplicate row indices collected for each CSC column.
"""
function _sem_finalize_column_pattern!(
    columns::Vector{Vector{Int}}
)
    for rows in columns
        sort!(rows)
        unique!(rows)
    end

    return columns
end

"""
    _sem_csc_from_column_pattern(m, n, columns)

Allocate an empty `SparseMatrixCSC{Float64,Int}` with an already known
column-wise sparsity pattern. No I/J/V triplet representation is created.
"""
function _sem_csc_from_column_pattern(
    m::Integer,
    n::Integer,
    columns::Vector{Vector{Int}}
)
    length(columns) == n ||
        error("Internal spectral CSC error: invalid number of columns.")

    _sem_finalize_column_pattern!(columns)

    colptr = Vector{Int}(undef, Int(n) + 1)
    colptr[1] = 1

    total_nnz = 0
    @inbounds for j in 1:Int(n)
        total_nnz += length(columns[j])
        colptr[j + 1] = total_nnz + 1
    end

    rowval = Vector{Int}(undef, total_nnz)
    nzval = zeros(Float64, total_nnz)

    pos = 1
    @inbounds for j in 1:Int(n)
        rows = columns[j]
        if !isempty(rows)
            copyto!(rowval, pos, rows, 1, length(rows))
            pos += length(rows)
        end
    end

    return SparseMatrixCSC(
        Int(m),
        Int(n),
        colptr,
        rowval,
        nzval
    )
end

@inline function _sem_csc_position(
    rowval::Vector{Int},
    lo::Int,
    hi::Int,
    row::Int
)
    while lo <= hi
        mid = (lo + hi) >>> 1
        r = rowval[mid]

        if r < row
            lo = mid + 1
        elseif r > row
            hi = mid - 1
        else
            return mid
        end
    end

    return 0
end

"""
    _sem_set_unique!(A, row, col, value; atol=1e-12)

Write a value into a preallocated CSC pattern. If a neighbouring element has
already generated the same global entry, verify that both values agree rather
than summing them.
"""
function _sem_set_unique!(
    A::SparseMatrixCSC{Float64,Int},
    row::Integer,
    col::Integer,
    value::Real;
    atol::Real=1e-12
)
    v = Float64(value)
    _sem_keep_entry(v) || return nothing

    j = Int(col)
    i = Int(row)
    lo = A.colptr[j]
    hi = A.colptr[j + 1] - 1
    pos = _sem_csc_position(A.rowval, lo, hi, i)

    pos != 0 ||
        error(
            "Internal spectral CSC error: entry ($i, $j) is missing " *
            "from the precomputed pattern."
        )

    old = A.nzval[pos]

    if iszero(old)
        A.nzval[pos] = v
    else
        isapprox(old, v; atol=atol, rtol=atol) ||
            error(
                "spectralTransformation: inconsistent transformation " *
                "entry ($i, $j): $old vs. $v."
            )
    end

    return nothing
end

"""
    _sem_nodal_prolongation(P, data, numbered)

Build the scalar nodal FEM-from-GLL prolongation directly in CSC format.
"""
function _sem_nodal_prolongation(
    P::Problem,
    data,
    numbered
)
    elements = data.elements
    type_cache = data.type_cache
    spectral_conn = numbered.connectivity
    nSpectral = numbered.count

    columns = [Int[] for _ in 1:nSpectral]

    for ie in eachindex(elements)
        elem = elements[ie]
        cache = type_cache[elem.et]
        sconn = spectral_conn[ie]

        @inbounds for b in 1:cache.nSem
            col = sconn[b]
            rows = columns[col]

            for a in 1:cache.nFem
                v = cache.Te[a, b]
                _sem_keep_entry(v) || continue

                push!(
                    rows,
                    _sem_global_node(P, elem.nodes[a])
                )
            end
        end
    end

    T = _sem_csc_from_column_pattern(
        P.non,
        nSpectral,
        columns
    )

    for ie in eachindex(elements)
        elem = elements[ie]
        cache = type_cache[elem.et]
        sconn = spectral_conn[ie]

        @inbounds for b in 1:cache.nSem
            col = sconn[b]

            for a in 1:cache.nFem
                v = cache.Te[a, b]
                _sem_keep_entry(v) || continue

                row = _sem_global_node(P, elem.nodes[a])
                _sem_set_unique!(T, row, col, v)
            end
        end
    end

    return T
end

"""
    _sem_spectral_owners(data, numbered)

For each shared spectral node, store the first element/local-node pair that
represents it. This reproduces the restriction convention used by the original
implementation without a global entry dictionary.
"""
function _sem_spectral_owners(
    data,
    numbered
)
    nSpectral = numbered.count
    owner_element = zeros(Int, nSpectral)
    owner_local = zeros(Int, nSpectral)

    for ie in eachindex(data.elements)
        cache = data.type_cache[data.elements[ie].et]
        sconn = numbered.connectivity[ie]

        @inbounds for a in 1:cache.nSem
            s = sconn[a]

            if owner_element[s] == 0
                owner_element[s] = ie
                owner_local[s] = a
            end
        end
    end

    all(x -> !iszero(x), owner_element) ||
        error("Internal spectral error: unowned spectral node detected.")

    return owner_element, owner_local
end

"""
    _sem_nodal_restriction(P, data, numbered)

Build the scalar nodal GLL-from-FEM restriction directly in CSC format.
"""
function _sem_nodal_restriction(
    P::Problem,
    data,
    numbered
)
    nSpectral = numbered.count
    owner_element, owner_local =
        _sem_spectral_owners(data, numbered)

    # R has full FEM nodes as columns and compact spectral nodes as rows.
    columns = [Int[] for _ in 1:P.non]

    @inbounds for spectral_node in 1:nSpectral
        ie = owner_element[spectral_node]
        a = owner_local[spectral_node]
        elem = data.elements[ie]
        cache = data.type_cache[elem.et]

        for b in 1:cache.nFem
            v = cache.Re[a, b]
            _sem_keep_entry(v) || continue

            col = _sem_global_node(P, elem.nodes[b])
            push!(columns[col], spectral_node)
        end
    end

    R = _sem_csc_from_column_pattern(
        nSpectral,
        P.non,
        columns
    )

    @inbounds for spectral_node in 1:nSpectral
        ie = owner_element[spectral_node]
        a = owner_local[spectral_node]
        elem = data.elements[ie]
        cache = data.type_cache[elem.et]

        for b in 1:cache.nFem
            v = cache.Re[a, b]
            _sem_keep_entry(v) || continue

            col = _sem_global_node(P, elem.nodes[b])
            _sem_set_unique!(R, spectral_node, col, v)
        end
    end

    return R
end

"""
    _sem_expand_components(A, pdim)

Expand a scalar nodal transformation to node-major field DOF ordering without
forming a Kronecker product or an intermediate triplet representation.
"""
function _sem_expand_components(
    A::SparseMatrixCSC{Float64,Int},
    pdim::Integer
)
    d = Int(pdim)
    d >= 1 || error("Internal spectral error: pdim must be positive.")
    d == 1 && return A

    m, n = size(A)
    nnzA = nnz(A)

    colptr = Vector{Int}(undef, n * d + 1)
    rowval = Vector{Int}(undef, nnzA * d)
    nzval = Vector{Float64}(undef, nnzA * d)

    pos = 1
    colptr[1] = 1

    @inbounds for j in 1:n
        p1 = A.colptr[j]
        p2 = A.colptr[j + 1] - 1

        for comp in 1:d
            jfull = (j - 1) * d + comp

            for p in p1:p2
                rowval[pos] = (A.rowval[p] - 1) * d + comp
                nzval[pos] = A.nzval[p]
                pos += 1
            end

            colptr[jfull + 1] = pos
        end
    end

    pos == nnzA * d + 1 ||
        error("Internal spectral CSC error while expanding field components.")

    return SparseMatrixCSC(
        m * d,
        n * d,
        colptr,
        rowval,
        nzval
    )
end

"""
    _sem_transformation_matrices(P, data, numbered)

Build the full field prolongation and restriction matrices directly in CSC
format. The scalar nodal pattern is assembled once and then expanded over field
components.
"""
function _sem_transformation_matrices(
    P::Problem,
    data,
    numbered
)
    Tnode = _sem_nodal_prolongation(P, data, numbered)
    Rnode = _sem_nodal_restriction(P, data, numbered)

    T = _sem_expand_components(Tnode, P.pdim)
    R = _sem_expand_components(Rnode, P.pdim)

    return T, R
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
# Cached spectral approximation data
# =============================================================================

"""
    SpectralData

Basis-specific metadata cached by a spectral [`Problem`](@ref). The large
solver transformations themselves are stored once in `P.approximation.T` and
`P.approximation.R`; this object keeps the GLL topology and element data needed
by quadrature and diagnostics.
"""
struct SpectralData
    dimension::Int
    order::Int
    ξ_gll::Vector{Float64}
    w_gll::Vector{Float64}
    local_gll_coordinates::Matrix{Float64}
    spectral_coordinates::Matrix{Float64}
    spectral_connectivity::Vector{Vector{Int}}
    element_tags::Vector{UInt64}
    active_dofs::Vector{Int}
    global_dofs::Int
    spectral_dofs::Int

    # Compact internal data retained for GLL quadrature. The temporary
    # transformation build cache (including local Te/Re matrices and node
    # dictionaries) is deliberately not retained.
    elements::Vector{NamedTuple}
    integration_cache::Dict{Int,Any}
    node_coordinates::Matrix{Float64}
end

function _sem_prepare_spectral_mesh(P::Problem)
    if P.dim == 1
        data = _sem_collect_1d_elements(P)
        numbered = _sem_number_spectral_nodes_1d(
            data.elements,
            data.type_cache,
            data.nodeCoord
        )
        return data, numbered

    elseif P.dim == 2
        data = _sem_collect_2d_quad_elements(P)
        numbered = _sem_number_spectral_nodes_2d(
            data.elements,
            data.type_cache,
            data.nodeCoord
        )
        return data, numbered

    else
        error(
            "spectralTransformation: dimension $(P.dim) is not yet " *
            "implemented. Supported dimensions are 1 and 2."
        )
    end
end

function _sem_active_rows(T::SparseMatrixCSC{Float64,Int})
    represented = falses(size(T, 1))

    @inbounds for row in T.rowval
        represented[row] = true
    end

    return findall(represented)
end

function _sem_inverse_error(T, R)
    n = size(T, 2)
    Isem = spdiagm(0 => ones(Float64, n))
    return norm(R * T - Isem, Inf)
end

function _build_spectral_approximation(
    P::Problem;
    check::Bool=false
)
    P.reducedOrder &&
        error(
            "spectralTransformation: basis=:spectral combined with " *
            "reducedOrder=true is not yet implemented."
        )

    data, numbered = _sem_prepare_spectral_mesh(P)
    T, R = _sem_transformation_matrices(P, data, numbered)

    active = sort!(unique(Int.(allDoFs(P))))
    represented = _sem_active_rows(T)

    represented == active ||
        error(
            "spectralTransformation: represented global DOFs do not " *
            "match allDoFs(P)."
        )

    inverse_error = check ? _sem_inverse_error(T, R) : nothing

    inverse_error === nothing || inverse_error <= 1e-10 ||
        error(
            "spectralTransformation: R*T differs from identity; " *
            "||R*T-I||∞ = $inverse_error."
        )

    first_cache = data.type_cache[data.elements[1].et]

    local_gll =
        P.dim == 1 ?
        reshape(copy(first_cache.ξGLL), :, 1) :
        copy(first_cache.localGLL)

    integration_cache = Dict{Int,Any}()

    for (et, cache) in data.type_cache
        integration_cache[et] = (
            nFem=cache.nFem,
            nSem=cache.nSem,
            n1=P.dim == 2 ? cache.n1 : cache.nSem,
            ξGLL=copy(cache.ξGLL),
            wGLL=copy(cache.wGLL),
            localGLL=P.dim == 2 ?
                copy(cache.localGLL) :
                zeros(Float64, 0, 0)
        )
    end

    node_coordinates = zeros(Float64, P.non, 3)

    for (tag, xyz) in data.nodeCoord
        i = Int(tag)
        node_coordinates[i, 1] = xyz[1]
        node_coordinates[i, 2] = xyz[2]
        node_coordinates[i, 3] = xyz[3]
    end

    spectral_data = SpectralData(
        P.dim,
        data.order,
        copy(first_cache.ξGLL),
        copy(first_cache.wGLL),
        local_gll,
        numbered.coordinates,
        numbered.connectivity,
        UInt64[elem.tag for elem in data.elements],
        active,
        ndofs(P),
        numbered.count * P.pdim,
        data.elements,
        integration_cache,
        node_coordinates
    )

    return T, R, spectral_data, inverse_error
end

"""
    _initialize_spectral_approximation!(P)

Build the spectral transformation and metadata once and store them in the
mutable approximation cache owned by `P`.
"""
function _initialize_spectral_approximation!(P::Problem)
    P.basis === :spectral ||
        error("Internal error: spectral cache requested for a non-spectral Problem.")

    T, R, metadata, _ = _build_spectral_approximation(P)

    P.approximation.T = T
    P.approximation.R = R
    P.approximation.metadata = metadata

    return P.approximation
end

function _spectral_data(P::Problem)
    data = _ensure_approximation_data!(P)

    data.metadata isa SpectralData ||
        error("Spectral metadata is not available for this Problem.")

    return data.metadata
end

function _spectral_details(
    P::Problem,
    data::SpectralData;
    inverse_error=nothing
)
    return (
        T=P.approximation.T,
        R=P.approximation.R,
        dimension=data.dimension,
        order=data.order,
        ξ_gll=data.ξ_gll,
        w_gll=data.w_gll,
        local_gll_coordinates=data.local_gll_coordinates,
        spectral_coordinates=data.spectral_coordinates,
        spectral_connectivity=data.spectral_connectivity,
        element_tags=data.element_tags,
        active_dofs=data.active_dofs,
        global_dofs=data.global_dofs,
        spectral_dofs=data.spectral_dofs,
        inverse_error=inverse_error
    )
end

# =============================================================================
# Public transformation
# =============================================================================

"""
    spectralTransformation(P::Problem; details=false, check=false)

Return a Legendre-Gauss-Lobatto spectral basis transformation for `P`.

For `basis=:spectral`, the transformation and spectral mesh metadata are built
when the `Problem` is initialized and reused here. The ordinary
LowLevelFEM/Gmsh Lagrange representation remains the global storage space:

    u_global = T * u_spectral
    u_spectral = R * u_global

and an assembled system is projected as

    K_spectral = T' * K * T
    f_spectral = T' * f

The transformation matrices are assembled directly in CSC format; no global
`Dict{Tuple{Int,Int},...}` or I/J/V sparse construction is used.

# Current limitations

- 1D `Line` and complete tensor-product 2D `Quadrangle` elements are supported,
- triangular spectral elements are not yet implemented,
- all domain elements must have the same polynomial order,
- `basis=:spectral, reducedOrder=true` is reserved by the API but not yet
  implemented,
- spectral basis combined with MPCs is not yet implemented by the solver.

# Keyword arguments

- `details=false`: return only `(T, R)`.
- `details=true`: return a named tuple with cached transformation and spectral
  mesh information.
- `check=false`: if `true`, verify `R*T ≈ I`.

Calling `spectralTransformation` on a non-spectral `Problem` remains supported
as an explicit diagnostic utility; in that case the result is built on demand
and is not cached in the `Problem`.
"""
function spectralTransformation(
    P::Problem;
    details::Bool=false,
    check::Bool=false
)
    if P.basis === :spectral
        cache = _ensure_approximation_data!(P)
        metadata = _spectral_data(P)

        inverse_error = check ?
            _sem_inverse_error(cache.T, cache.R) :
            nothing

        inverse_error === nothing || inverse_error <= 1e-10 ||
            error(
                "spectralTransformation: R*T differs from identity; " *
                "||R*T-I||∞ = $inverse_error."
            )

        return details ?
            _spectral_details(
                P,
                metadata;
                inverse_error=inverse_error
            ) :
            (cache.T, cache.R)
    end

    T, R, metadata, inverse_error =
        _build_spectral_approximation(P; check=check)

    if !details
        return T, R
    end

    return (
        T=T,
        R=R,
        dimension=metadata.dimension,
        order=metadata.order,
        ξ_gll=metadata.ξ_gll,
        w_gll=metadata.w_gll,
        local_gll_coordinates=metadata.local_gll_coordinates,
        spectral_coordinates=metadata.spectral_coordinates,
        spectral_connectivity=metadata.spectral_connectivity,
        element_tags=metadata.element_tags,
        active_dofs=metadata.active_dofs,
        global_dofs=metadata.global_dofs,
        spectral_dofs=metadata.spectral_dofs,
        inverse_error=inverse_error
    )
end


# =============================================================================
# Experimental GLL quadrature
# =============================================================================

"""
    _sem_restriction_matrix(P, data, numbered)

Construct the full field restriction directly in CSC format from already
prepared spectral mesh data. This helper is retained for diagnostic/on-demand
paths; cached spectral Problems normally reuse `P.approximation.R` instead.
"""
function _sem_restriction_matrix(
    P::Problem,
    data,
    numbered
)
    Rnode = _sem_nodal_restriction(P, data, numbered)
    return _sem_expand_components(Rnode, P.pdim)
end


"""
    _sem_gll_reference_points_3d(P, cache)

Return the local GLL points in the three-coordinate layout expected by Gmsh.
"""
function _sem_gll_reference_points_3d(
    P::Problem,
    cache
)
    points = zeros(Float64, 3 * cache.nSem)

    if P.dim == 1
        @inbounds for q in 1:cache.nSem
            points[3q - 2] = cache.ξGLL[q]
        end

    elseif P.dim == 2
        @inbounds for q in 1:cache.nSem
            points[3q - 2] = cache.localGLL[q, 1]
            points[3q - 1] = cache.localGLL[q, 2]
        end

    else
        error(
            "GLL quadrature: dimension $(P.dim) is not yet implemented. " *
            "Supported dimensions are 1 and 2."
        )
    end

    return points
end


"""
    _sem_geometry_gradients_at_gll(P, element_type, cache)

Evaluate the original Gmsh geometry basis gradients at the GLL points. The
returned matrix has one column per GLL point and uses the same `3*numNodes`
row layout as the ordinary LowLevelFEM element kernel.
"""
function _sem_geometry_gradients_at_gll(
    P::Problem,
    element_type::Integer,
    cache
)
    points =
        _sem_gll_reference_points_3d(
            P,
            cache
        )

    ncomp, dfun, _ =
        gmsh.model.mesh.getBasisFunctions(
            Int(element_type),
            points,
            "GradLagrange"
        )

    ncomp == 3 ||
        error(
            "GLL quadrature: expected three reference-gradient components " *
            "from Gmsh, got $ncomp."
        )

    return reshape(
        Float64.(dfun),
        3 * cache.nFem,
        cache.nSem
    )
end


"""
    _sem_geometry_measure_at_gll(P, elem, cache, nodeCoord, grad, q)

Return the physical integration measure at GLL point `q` using the original
high-order Gmsh isoparametric geometry.
"""
@inline function _sem_geometry_measure_at_gll(
    P::Problem,
    elem,
    cache,
    nodeCoord,
    grad,
    q::Integer
)
    if P.dim == 1
        dx_dξ = 0.0

        @inbounds for a in 1:cache.nFem
            x = nodeCoord[Int(elem.nodes[a]), 1]
            dNa_dξ = grad[3a - 2, q]
            dx_dξ += x * dNa_dξ
        end

        measure = abs(dx_dξ)

    elseif P.dim == 2
        dx_dξ = 0.0
        dx_dη = 0.0
        dy_dξ = 0.0
        dy_dη = 0.0

        @inbounds for a in 1:cache.nFem
            node = Int(elem.nodes[a])
            x = nodeCoord[node, 1]
            y = nodeCoord[node, 2]

            dNa_dξ = grad[3a - 2, q]
            dNa_dη = grad[3a - 1, q]

            dx_dξ += x * dNa_dξ
            dx_dη += x * dNa_dη
            dy_dξ += y * dNa_dξ
            dy_dη += y * dNa_dη
        end

        measure =
            abs(
                dx_dξ * dy_dη -
                dx_dη * dy_dξ
            )

    else
        error(
            "GLL quadrature: dimension $(P.dim) is not yet implemented."
        )
    end

    measure > 0.0 ||
        error(
            "GLL quadrature: singular geometry in element $(elem.tag) " *
            "at local GLL point $q."
        )

    return measure
end


"""
    _sem_gll_quadrature_weight(P, cache, q)

Return the reference-domain tensor-product GLL quadrature weight for local
spectral node `q`.
"""
@inline function _sem_gll_quadrature_weight(
    P::Problem,
    cache,
    q::Integer
)
    if P.dim == 1
        return cache.wGLL[q]

    elseif P.dim == 2
        i = mod(q - 1, cache.n1) + 1
        j = div(q - 1, cache.n1) + 1
        return cache.wGLL[i] * cache.wGLL[j]

    else
        error(
            "GLL quadrature: dimension $(P.dim) is not yet implemented."
        )
    end
end


"""
    _sem_domain_element_tags(P, domain)

Return the element tags selected by a volume domain. `nothing` means that all
domain elements of `P` are used. Boundary GLL integration is not yet
implemented.
"""
function _sem_domain_element_tags(
    P::Problem,
    domain
)
    domain === nothing &&
        return nothing

    hasproperty(domain, :kind) &&
    hasproperty(domain, :name) ||
        error(
            "GLL quadrature: unsupported domain specification $(typeof(domain))."
        )

    domain.kind === :Ω ||
        error(
            "GLL quadrature is currently implemented only for volume/domain " *
            "integrals (Ω), not boundary integrals (Γ)."
        )

    gmsh.model.setCurrent(P.name)

    dimTags =
        gmsh.model.getEntitiesForPhysicalName(
            domain.name
        )

    isempty(dimTags) &&
        error(
            "GLL quadrature: physical group \"$(domain.name)\" not found."
        )

    selected = Set{UInt64}()

    for (edim, etag) in dimTags
        edim == P.dim ||
            error(
                "GLL quadrature: Ω=\"$(domain.name)\" has dimension $edim, " *
                "but problem.dim=$(P.dim)."
            )

        _, elemTags, _ =
            gmsh.model.mesh.getElements(
                edim,
                etag
            )

        for tags in elemTags
            for tag in tags
                push!(selected, UInt64(tag))
            end
        end
    end

    isempty(selected) &&
        error(
            "GLL quadrature: no elements found in Ω=\"$(domain.name)\"."
        )

    return selected
end


"""
    _spectral_gll_mass_matrix(P; domain=nothing)

Assemble the unit-coefficient mass matrix of `P` with GLL quadrature in the
spectral nodal basis and immediately pull it back to the ordinary LowLevelFEM
global representation.

The spectral-space matrix is diagonal. The returned `SystemMatrix` generally
is not diagonal because it is represented in the ordinary Gmsh/Lagrange global
storage basis:

    M_global = R' * M_GLL * R

The existing solver-side spectral transformation therefore recovers
`M_GLL` (up to floating-point roundoff) without changing the `solveField` API.
"""
function _spectral_gll_mass_matrix(
    P::Problem;
    domain=nothing
)
    approximation = _ensure_approximation_data!(P)
    spectral = _spectral_data(P)

    selected =
        _sem_domain_element_tags(
            P,
            domain
        )

    if selected !== nothing
        available =
            Set(
                elem.tag for elem in spectral.elements
            )

        missing = setdiff(selected, available)

        isempty(missing) ||
            error(
                "GLL quadrature: Ω=\"$(domain.name)\" contains elements that " *
                "are not part of the spectral Problem."
            )
    end

    d = P.pdim
    nSpectralDofs = spectral.spectral_dofs
    diagonal = zeros(Float64, nSpectralDofs)

    grad_cache = Dict{Int,Matrix{Float64}}()
    integrated_elements = 0

    for ie in eachindex(spectral.elements)
        elem = spectral.elements[ie]

        selected !== nothing &&
        !(elem.tag in selected) &&
            continue

        integrated_elements += 1

        cache = spectral.integration_cache[elem.et]
        sconn = spectral.spectral_connectivity[ie]

        grad =
            get!(grad_cache, elem.et) do
                _sem_geometry_gradients_at_gll(
                    P,
                    elem.et,
                    cache
                )
            end

        @inbounds for q in 1:cache.nSem
            measure =
                _sem_geometry_measure_at_gll(
                    P,
                    elem,
                    cache,
                    spectral.node_coordinates,
                    grad,
                    q
                )

            w =
                _sem_gll_quadrature_weight(
                    P,
                    cache,
                    q
                )

            contribution = measure * w
            spectral_node = sconn[q]
            first_dof = (spectral_node - 1) * d

            for comp in 1:d
                diagonal[first_dof + comp] += contribution
            end
        end
    end

    integrated_elements > 0 ||
        error(
            "GLL quadrature: no spectral elements were selected for integration."
        )

    R = approximation.R

    Mhat =
        spdiagm(
            0 => diagonal
        )

    Mglobal =
        sparse(
            transpose(R) *
            Mhat *
            R
        )

    dropzeros!(Mglobal)

    return SystemMatrix(
        Mglobal,
        P,
        P
    )
end


"""
    spectralIntegral(term::BilinearTerm; domain=nothing, weight=nothing,
                     threads=:auto)

Assemble a bilinear form with Gauss-Lobatto-Legendre quadrature.

Current scope is deliberately narrow:

- both test and trial fields must use `basis=:spectral`,
- test and trial must be the same `Problem` instance,
- only `Id(P) ⋅ Id(P)` is implemented,
- only the unit coefficient is implemented,
- additional `weight` factors are not yet implemented,
- only volume/domain integration is implemented,
- supported elements are the same 1D lines and complete 2D tensor-product
  quadrangles supported by `spectralTransformation`.

The GLL matrix is assembled in the compact spectral space and then immediately
returned in the ordinary LowLevelFEM global representation so that existing
solvers and boundary-condition handling remain unchanged.
"""
function spectralIntegral(
    term::BilinearTerm;
    domain=nothing,
    weight=nothing,
    threads=:auto
)
    Ps = term.a.P
    Pu = term.b.P

    Ps.basis === :spectral ||
        error(
            "GLL quadrature requires basis=:spectral on the test field."
        )

    Pu.basis === :spectral ||
        error(
            "GLL quadrature requires basis=:spectral on the trial field."
        )

    Ps === Pu ||
        error(
            "GLL quadrature currently supports only P ⋅ P with the same " *
            "Problem instance on the test and trial sides."
        )

    term.a.op isa IdOp ||
        error(
            "GLL quadrature currently supports only P ⋅ P " *
            "(identity operator on the test field)."
        )

    term.b.op isa IdOp ||
        error(
            "GLL quadrature currently supports only P ⋅ P " *
            "(identity operator on the trial field)."
        )

    term.coef isa Number && isone(term.coef) ||
        error(
            "GLL quadrature currently supports only the unit-coefficient " *
            "bilinear form P ⋅ P."
        )

    weight === nothing ||
        error(
            "GLL quadrature with an additional weight is not yet implemented."
        )

    # Accepted for API compatibility. The first implementation is serial; the
    # expensive global transformation/assembly path can be optimized later.
    _ = threads

    return _spectral_gll_mass_matrix(
        Pu;
        domain=domain
    )
end
