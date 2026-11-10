mutable struct OperatorWorkspace
    Jac::Matrix{Float64}
    jacDet::Vector{Float64}
    invJac::Matrix{Float64}
    ∂h::Matrix{Float64}
    Bu::Matrix{Float64}
    Bs::Matrix{Float64}
    Ke::Matrix{Float64}
    tmp::Matrix{Float64}

    patBu::MatrixPattern
    patBsT::MatrixPattern
    patTmp::MatrixPattern
    patCgp::Union{Nothing,MatrixPattern}

    Cgp::Union{Nothing,Matrix{Float64}}
end

@inline function _eval_coefficient_at_gp!(
    W::Matrix{Float64},
    Cprep::AbstractMatrix,
    elem,
    hgp
    )
    @inbounds for I in CartesianIndices(Cprep)
        cij = Cprep[I]

        if cij isa Number
            W[I] = Float64(cij)
        elseif cij isa Dict{Int,Matrix{Float64}}
            W[I] = _coeff_at_gp(cij, elem, hgp)
        else
            W[I] = _coeff_at_gp(cij, elem, hgp)
        end
    end

    return W
end

"""
    OperatorWorkspace(
        dim,
        num_nodes,
        num_integration_points,
        out_u,
        out_s,
        ndofs_u,
        ndofs_s
    )

Allocate reusable element-level work arrays for operator assembly.
"""
function OperatorWorkspace(
    dim::Int,
    num_nodes::Int,
    num_integration_points::Int,
    out_u::Int,
    out_s::Int,
    ndofs_u::Int,
    ndofs_s::Int
    )
    Bu = zeros(out_u, ndofs_u)
    Bs = zeros(out_s, ndofs_s)
    tmp = zeros(ndofs_s, out_u)

    patBu  = MatrixPattern(out_u, ndofs_u)
    patBsT = MatrixPattern(ndofs_s, out_s)
    patTmp = MatrixPattern(ndofs_s, out_u)

    sizehint!(patBu.rows, length(Bu))
    sizehint!(patBu.cols, length(Bu))

    sizehint!(patBsT.rows, length(Bs))
    sizehint!(patBsT.cols, length(Bs))

    sizehint!(patTmp.rows, length(tmp))
    sizehint!(patTmp.cols, length(tmp))

    return OperatorWorkspace(
        zeros(3, 3 * num_integration_points),
        zeros(num_integration_points),
        zeros(3, 3 * num_integration_points),
        zeros(dim, num_nodes * num_integration_points),
        Bu,
        Bs,
        zeros(ndofs_s, ndofs_u),
        tmp,
        patBu,
        patBsT,
        patTmp,
        nothing,
        nothing
    )
end

"""
    dense_node_coordinates(problem)

Return a `3 × number_of_nodes` coordinate table indexed by node tag.

The common Gmsh case, in which node tags are already ordered as `1:n`, returns
a reshaped view of the coordinate vector without copying it. A single compact
reordering copy is used otherwise. LowLevelFEM currently requires dense node
tags because global degrees of freedom are derived directly from node tags.
"""
function dense_node_coordinates(problem::Problem)
    gmsh.model.setCurrent(problem.name)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()

    number_of_nodes, remainder = divrem(ndofs(problem), problem.pdim)
    remainder == 0 || error("dense_node_coordinates: invalid number of degrees of freedom.")
    length(node_tags) == number_of_nodes || error(
        "dense_node_coordinates: the mesh contains $(length(node_tags)) nodes, " *
        "but the problem expects $number_of_nodes."
    )

    ordered = true
    @inbounds for i in eachindex(node_tags)
        if Int(node_tags[i]) != i
            ordered = false
            break
        end
    end

    if ordered
        return reshape(coordinates, 3, number_of_nodes)
    end

    dense = Matrix{Float64}(undef, 3, number_of_nodes)
    @inbounds for source in eachindex(node_tags)
        node = Int(node_tags[source])
        1 <= node <= number_of_nodes || error(
            "dense_node_coordinates: node tag $node is outside 1:$number_of_nodes."
        )
        source3 = 3 * (source - 1)
        dense[1, node] = coordinates[source3 + 1]
        dense[2, node] = coordinates[source3 + 2]
        dense[3, node] = coordinates[source3 + 3]
    end

    return dense
end

@inline function jacobian_smatrix(Jac::AbstractMatrix{<:Real}, first_col::Int)
    return SMatrix{3,3,Float64,9}((
        Jac[1, first_col],     Jac[2, first_col],     Jac[3, first_col],
        Jac[1, first_col + 1], Jac[2, first_col + 1], Jac[3, first_col + 1],
        Jac[1, first_col + 2], Jac[2, first_col + 2], Jac[3, first_col + 2]
    ))
end

"""
    compute_element_geometry!(ws, coordinates, connectivity, first_node,
                              num_nodes, num_integration_points, dim, edim, grad_h)

Compute Jacobians and integration measures for one element directly from its
nodal coordinates. Only the reusable arrays in `ws` are written; no
element-sized heap allocation is required.
"""
function compute_element_geometry!(
    ws::OperatorWorkspace,
    coordinates::AbstractMatrix{<:Real},
    connectivity::AbstractVector{<:Integer},
    first_node::Integer,
    numNodes::Integer,
    numIntPoints::Integer,
    dim::Integer,
    edim::Integer,
    ∇h
)
    Jac = ws.Jac
    jacDet = ws.jacDet
    fill!(Jac, 0.0)

    @inbounds for k in 1:numIntPoints
        first_col = 3k - 2

        for a in 1:numNodes
            node = Int(connectivity[first_node + a])
            grad_row = 3a - 2

            for α in 1:edim
                dNa = ∇h[grad_row + α - 1, k]
                col = first_col + α - 1
                Jac[1, col] += coordinates[1, node] * dNa
                Jac[2, col] += coordinates[2, node] * dNa
                Jac[3, col] += coordinates[3, node] * dNa
            end
        end

        if edim == dim
            if dim == 1
                jacDet[k] = abs(Jac[1, first_col])
                Jac[2, first_col + 1] = 1.0
                Jac[3, first_col + 2] = 1.0
            elseif dim == 2
                jacDet[k] = abs(
                    Jac[1, first_col] * Jac[2, first_col + 1] -
                    Jac[1, first_col + 1] * Jac[2, first_col]
                )
                Jac[3, first_col + 2] = 1.0
            else
                Jk = jacobian_smatrix(Jac, first_col)
                jacDet[k] = abs(det(Jk))
            end
        elseif edim == 1
            measure2 = 0.0
            for i in 1:dim
                measure2 += Jac[i, first_col]^2
            end
            jacDet[k] = sqrt(measure2)
        elseif edim == 2
            g11 = 0.0
            g12 = 0.0
            g22 = 0.0
            for i in 1:dim
                a1 = Jac[i, first_col]
                a2 = Jac[i, first_col + 1]
                g11 += a1 * a1
                g12 += a1 * a2
                g22 += a2 * a2
            end
            jacDet[k] = sqrt(max(0.0, g11 * g22 - g12 * g12))
        else
            error(
                "compute_element_geometry!: unsupported entity dimension " *
                "$edim for problem dimension $dim."
            )
        end

        jacDet[k] > 0.0 || error(
            "compute_element_geometry!: singular geometry at integration point $k."
        )
    end

    return nothing
end

"""
    assemble_element_matrix!(ws, elem, gidx, ...)

Compute the local element matrix from the batched Gmsh geometry and store it
in `ws.Ke`. The Jacobian data are accessed through views; no element-level
geometry copy is made.
"""
function assemble_element_matrix!(
    ws::OperatorWorkspace,
    elem::Integer,
    gidx::Integer,
    jac_all,
    det_all,
    nip::Integer,
    numIntPoints::Integer,
    dim::Integer,
    edim::Integer,
    numNodes::Integer,
    ∇h,
    h,
    intWeights,
    Cprep,
    Wprep,
    op_u::AbstractOp,
    Pu::Problem,
    op_s::AbstractOp,
    Ps::Problem
)
    nip == numIntPoints || error(
        "assemble_element_matrix!: inconsistent number of integration points."
    )

    jac_slice =
        @view jac_all[gidx * 9 * nip + 1 : (gidx + 1) * 9 * nip]
    jacDet =
        @view det_all[gidx * nip + 1 : (gidx + 1) * nip]
    Jac = reshape(jac_slice, 3, :)

    return assemble_element_matrix_with_geometry!(
        ws,
        Jac,
        jacDet,
        elem,
        numIntPoints,
        dim,
        edim,
        numNodes,
        ∇h,
        h,
        intWeights,
        Cprep,
        Wprep,
        op_u,
        Pu,
        op_s,
        Ps
    )
end

"""
    assemble_element_matrix_from_geometry!(ws, elem, ...)

Compute the local element matrix using geometry already stored in `ws.Jac`
and `ws.jacDet`. This is shared by the legacy batched-Gmsh path and the
low-memory on-the-fly geometry path.
"""
function assemble_element_matrix_from_geometry!(
    ws::OperatorWorkspace,
    elem::Integer,
    numIntPoints::Integer,
    dim::Integer,
    edim::Integer,
    numNodes::Integer,
    ∇h,
    h,
    intWeights,
    Cprep,
    Wprep,
    op_u::AbstractOp,
    Pu::Problem,
    op_s::AbstractOp,
    Ps::Problem
)
    return assemble_element_matrix_with_geometry!(
        ws,
        ws.Jac,
        ws.jacDet,
        elem,
        numIntPoints,
        dim,
        edim,
        numNodes,
        ∇h,
        h,
        intWeights,
        Cprep,
        Wprep,
        op_u,
        Pu,
        op_s,
        Ps
    )
end

"""
    assemble_element_matrix_with_geometry!(ws, Jac, jacDet, elem, ...)

Compute the local element matrix using explicitly supplied Jacobian and
integration-measure arrays. This keeps the element kernel shared while allowing
the IJV path to use zero-copy Gmsh views and the CSC path to use workspace-owned
geometry.
"""
function assemble_element_matrix_with_geometry!(
    ws::OperatorWorkspace,
    Jac::AbstractMatrix{<:Real},
    jacDet::AbstractVector{<:Real},
    elem::Integer,
    numIntPoints::Integer,
    dim::Integer,
    edim::Integer,
    numNodes::Integer,
    ∇h,
    h,
    intWeights,
    Cprep,
    Wprep,
    op_u::AbstractOp,
    Pu::Problem,
    op_s::AbstractOp,
    Ps::Problem
    )
    invJac = ws.invJac
    ∂h     = ws.∂h
    Bu     = ws.Bu
    Bs     = ws.Bs
    Ke     = ws.Ke
    tmp    = ws.tmp

    @inbounds for k in 1:numIntPoints
        if edim == dim
            Jk = jacobian_smatrix(Jac, 3k - 2)
            invJac[1:3, 3k-2:3k] .= transpose(inv3x3_fast(Jk))
        end
    end

    # Physical / tangential gradients
    fill!(∂h, 0.0)

    @inbounds for k in 1:numIntPoints
        if edim == dim
            @views begin
                invJk = invJac[1:dim, 3k-2:3k-(3-dim)]

                for a in 1:numNodes
                    gha = ∇h[a*3-2:a*3-(3-dim), k]
                    dha = ∂h[1:dim, (k-1)*numNodes+a]
                    mul!(dha, invJk, gha)
                end
            end
        else
            Jk = @view Jac[1:dim, 3k-2:3k]
            Jtan = Matrix(Jk[:, 1:edim])
            Gtan = Jtan * inv(Jtan' * Jtan)

            for a in 1:numNodes
                gha = ∇h[a*3-2:a*3-3+edim, k]
                ∂h[1:dim, (k-1)*numNodes+a] .= Gtan * gha
            end
        end
    end

    fill!(Ke, 0.0)

    @inbounds for k in 1:numIntPoints
        if Cprep isa AbstractVector
            Cgp = _eval_coefficient_chain_at_gp(
                Cprep, elem, view(h, :, k)
            )

        elseif Cprep isa AbstractMatrix &&
               !(eltype(Cprep) <: Number)

            W = ws.Cgp

            if W === nothing || size(W) != size(Cprep)
                W = Matrix{Float64}(undef, size(Cprep)...)
                ws.Cgp = W
            end

            Cgp = _eval_coefficient_at_gp!(
                W, Cprep, elem, view(h, :, k)
            )

        else
            Cgp = _eval_coefficient_at_gp(
                Cprep, elem, view(h, :, k)
            )
        end

        #=
        if Cprep isa AbstractVector
            Cgp = _eval_coefficient_chain_at_gp(
                Cprep, elem, view(h, :, k)
            )
        else
            Cgp = _eval_coefficient_at_gp(
                Cprep, elem, view(h, :, k)
            )
        end
        =#

        wcoef =
            Wprep === nothing ?
            1.0 :
            _eval_coefficient_at_gp(Wprep, elem, view(h, :, k))

        if Cgp isa Number
            w = jacDet[k] * intWeights[k] * Cgp * wcoef

            if op_u isa SurfaceSymGradOp
                t1, t2 = surface_basis_from_J(Jac, k)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t1, t2)
            elseif op_u isa AxialGradOp
                t = @view Jac[1:Pu.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
            elseif op_u isa TangentialGradOp
                t = Jac[1:Pu.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
            else
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes)
            end

            if op_s isa SurfaceSymGradOp
                t1, t2 = surface_basis_from_J(Jac, k)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t1, t2)
            elseif op_s isa AxialGradOp
                t = @view Jac[1:Ps.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
            elseif op_s isa TangentialGradOp
                t = Jac[1:Ps.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
            else
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes)
            end

            detect_pattern!(ws.patBsT, transpose(Bs))
            detect_pattern!(ws.patBu, Bu)

            mul_opt!(
                Ke,
                transpose(Bs),
                Bu;
                patA = ws.patBsT,
                patB = ws.patBu,
                alpha = w,
                beta = 1.0
            )

        elseif Cgp isa AbstractMatrix
            w = jacDet[k] * intWeights[k] * wcoef

            if op_u isa SurfaceSymGradOp
                t1, t2 = surface_basis_from_J(Jac, k)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t1, t2)
            elseif op_u isa AxialGradOp
                t = @view Jac[1:Pu.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
            elseif op_u isa TangentialGradOp
                t = Jac[1:Pu.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
            else
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes)
            end

            if op_s isa SurfaceSymGradOp
                t1, t2 = surface_basis_from_J(Jac, k)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t1, t2)
            elseif op_s isa AxialGradOp
                t = @view Jac[1:Ps.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
            elseif op_s isa TangentialGradOp
                t = Jac[1:Ps.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
            else
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes)
            end

            if ws.patCgp === nothing
                ws.patCgp =
                    MatrixPattern(size(Cgp, 1), size(Cgp, 2))

                sizehint!(ws.patCgp.rows, length(Cgp))
                sizehint!(ws.patCgp.cols, length(Cgp))
            end

            detect_pattern!(ws.patBsT, transpose(Bs))
            detect_pattern!(ws.patCgp, Cgp)

            mul_opt!(
                tmp,
                transpose(Bs),
                Cgp;
                patA = ws.patBsT,
                patB = ws.patCgp,
                alpha = 1.0,
                beta = 0.0
            )

            detect_pattern!(ws.patTmp, tmp)
            detect_pattern!(ws.patBu, Bu)

            mul_opt!(
                Ke,
                tmp,
                Bu;
                patA = ws.patTmp,
                patB = ws.patBu,
                alpha = w,
                beta = 1.0
            )
        else
            error("assemble_element_matrix!: unsupported coefficient type")
        end
    end

    return nothing
end

"""
    build_csc_pattern(Pu::Problem, Ps::Problem; domain=nothing, Ω=nothing, Γ=nothing)

Build the exact structural CSC pattern required by the selected mesh domain.

The pattern is constructed from node connectivity without allocating full
degree-of-freedom triplet arrays. The returned matrix has zero-valued `nzval`
entries and sorted row indices in every column.

The construction is independent of element shape and interpolation order: it
uses the connectivity returned by Gmsh and expands every connected node pair to
the test/trial field components. Consequently, it supports all Gmsh Lagrange
element types that are otherwise supported by the selected LowLevelFEM
operators and element kernel.

The returned pattern can be supplied as `csc_matrix` to a bilinear or compound
integral. Serial assembly writes directly to its `nzval` array. Parallel
assembly uses the same array as the first worker's accumulator and allocates
one additional private `nzval` buffer per remaining worker, followed by an
allocation-free reduction.

Before reusing the pattern for a new independent assembly, reset its numerical
values with `fill!(K.nzval, 0.0)`. Do not call `dropzeros!(K)`, because that
would remove structural entries required by subsequent direct CSC scatter.

The caller is responsible for reusing a pattern only with the same mesh,
domain, and compatible test/trial spaces. Matrix dimensions alone do not prove
structural compatibility.
"""
function build_csc_pattern(
    Pu::Problem,
    Ps::Problem;
    domain=nothing,
    Ω=nothing,
    Γ=nothing
)
    @assert Pu.name == Ps.name "Both problems must refer to the same gmsh model/mesh."
    @assert Pu.dim == Ps.dim "Both problems must have the same spatial dimension."

    gmsh.model.setCurrent(Pu.name)

    if Ω !== nothing || Γ !== nothing
        domain === nothing || error(
            "build_csc_pattern: use either domain=... or Ω/Γ, not both."
        )
        domain = _domain_spec(; Ω=Ω, Γ=Γ)
    end

    nrows = ndofs(Ps)
    ncols = ndofs(Pu)

    nrow_nodes, rem_s = divrem(nrows, Ps.pdim)
    ncol_nodes, rem_u = divrem(ncols, Pu.pdim)

    rem_s == 0 || error("build_csc_pattern: invalid number of test-space degrees of freedom.")
    rem_u == 0 || error("build_csc_pattern: invalid number of trial-space degrees of freedom.")

    # Collect node-level adjacency first. Expanding to components only after
    # duplicate removal keeps temporary storage much smaller than DoF-level I/J.
    row_nodes_by_col_node = [Int[] for _ in 1:ncol_nodes]

    phNames = String[]
    domkind = nothing

    if domain === nothing
        phNames = [mat.phName for mat in Pu.material]
    elseif domain isa DomainSpec
        phNames = [domain.name]
        domkind = domain.kind
    elseif domain isa AbstractString
        phNames = [String(domain)]
    else
        error(
            "build_csc_pattern: unsupported domain type $(typeof(domain)). " *
            "Expected nothing, DomainSpec or String."
        )
    end

    for phName in phNames
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        isempty(dimTags) && error("build_csc_pattern: physical group \"$phName\" not found.")

        for (edim, etag) in dimTags
            if domkind === :Ω
                edim == Pu.dim ||
                    error("Ω=\"$phName\" has dim=$edim but problem.dim=$(Pu.dim)")
            elseif domkind === :Γ
                edim < Pu.dim ||
                    error("Γ=\"$phName\" has dim=$edim but expected < $(Pu.dim)")
            end

            elemTypes = gmsh.model.mesh.getElementTypes(edim, etag)

            for et in elemTypes
                _, _, _, numNodes::Int64, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                _, connectivity =
                    gmsh.model.mesh.getElementsByType(et, etag)
                nel = length(connectivity) ÷ numNodes

                @inbounds for e in 1:nel
                    first_node = (e - 1) * numNodes

                    for b in 1:numNodes
                        col_node = Int(connectivity[first_node + b])

                        1 <= col_node <= ncol_nodes || error(
                            "build_csc_pattern: node tag $col_node is outside " *
                            "the trial-space node range 1:$ncol_nodes."
                        )

                        rows = row_nodes_by_col_node[col_node]

                        for a in 1:numNodes
                            row_node = Int(connectivity[first_node + a])

                            1 <= row_node <= nrow_nodes || error(
                                "build_csc_pattern: node tag $row_node is outside " *
                                "the test-space node range 1:$nrow_nodes."
                            )

                            push!(rows, row_node)
                        end
                    end
                end
            end
        end
    end

    # Sort and remove duplicates in place. Unlike unique!, this pass needs no
    # auxiliary hash set because the vectors are already sorted.
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
            sizehint!(rows, write_pos)
        end

        total_node_nonzeros += length(rows)
    end

    total_nonzeros = Base.Checked.checked_mul(
        Base.Checked.checked_mul(total_node_nonzeros, Ps.pdim),
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
    @assert position == total_nonzeros + 1

    return SparseMatrixCSC(
        nrows,
        ncols,
        colptr,
        rowval,
        zeros(Float64, total_nonzeros)
    )
end

# The CSC scatter kernel addresses each nodal component block by offsets.
# Validate caller-owned patterns before using that optimization: dropping only
# some structural zeros can otherwise shift writes to unrelated matrix entries.
function _validate_csc_component_pattern(K::SparseMatrixCSC, pdim_u::Int, pdim_s::Int)
    pdim_u == 1 && pdim_s == 1 && return nothing
    invalid() = throw(ArgumentError(
        "csc_matrix has an incomplete component-block pattern. Rebuild it with " *
        "build_csc_pattern(...); keep structural zeros and avoid matrix " *
        "operations that remove them (use weight=... for scalar factors)."
    ))
    rows, columns = K.rowval, K.colptr
    @inbounds for first_col in 1:pdim_u:size(K, 2)
        first = columns[first_col]
        after = columns[first_col + 1]
        count = after - first
        count % pdim_s == 0 || invalid()
        if pdim_s > 1
            for position in first:pdim_s:after-1
                row = rows[position]
                (row - 1) % pdim_s == 0 || invalid()
                for component in 2:pdim_s
                    rows[position + component - 1] == row + component - 1 || invalid()
                end
            end
        end
        for component in 2:pdim_u
            col = first_col + component - 1
            start = columns[col]
            columns[col + 1] - start == count || invalid()
            for offset in 0:count-1
                rows[start + offset] == rows[first + offset] || invalid()
            end
        end
    end
    return nothing
end

"""
    find_csc_position(rowval, first, last, row)

Return the `nzval` position of `row` in a sorted CSC column range.

The search is allocation-free. An error indicates that the prebuilt sparsity
pattern and the element connectivity used for assembly are inconsistent.
"""
@inline function find_csc_position(
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

    error("find_csc_position: row $row is missing from the CSC pattern.")
end

"""
    resolve_num_threads(threads) -> Int

Resolve the requested number of assembly workers.

- `threads=:auto` uses all Julia threads in the default thread pool.
- `threads=1` selects serial assembly.
- A positive integer selects exactly that many worker tasks, up to
  `Threads.nthreads()`.
"""
function resolve_num_threads(threads)
    available = Threads.nthreads()

    if threads === :auto
        return available
    elseif threads isa Bool
        throw(ArgumentError(
            "`threads` must be `:auto` or a positive integer; " *
            "use `threads=1` instead of `false`."
        ))
    elseif threads isa Integer
        threads >= 1 || throw(ArgumentError(
            "`threads` must be `:auto` or a positive integer."
        ))
        threads <= available || throw(ArgumentError(
            "Requested threads=$threads, but Julia has only " *
            "$available thread(s) in the default thread pool."
        ))
        return Int(threads)
    else
        throw(ArgumentError(
            "`threads` must be `:auto` or a positive integer; " *
            "got $(repr(threads))."
        ))
    end
end

@inline function _worker_range(nitems::Int, nworkers::Int, worker::Int)
    base, extra = divrem(nitems, nworkers)
    count = base + (worker <= extra)
    first = (worker - 1) * base + min(worker - 1, extra) + 1
    return first:(first + count - 1)
end

function _run_workers(f, nworkers::Int)
    if nworkers == 1
        f(1)
        return nothing
    end

    @sync begin
        for worker in 1:nworkers
            let worker = worker
                Threads.@spawn f(worker)
            end
        end
    end
    return nothing
end

@inline function reduce_csc_worker!(
    target::Vector{Float64},
    buffers::Vector{Vector{Float64}},
    nworkers::Int,
    worker::Int
)
    @inbounds for p in _worker_range(length(target), nworkers, worker)
        value = target[p]
        for source in 2:nworkers
            value += buffers[source][p]
        end
        target[p] = value
    end

    return nothing
end

"""
    reduce_csc_buffers!(target, buffers, nworkers)

Accumulate private worker CSC value buffers into `target`. With one worker,
`target` already contains the complete assembly, so the reduction is skipped.
"""
function reduce_csc_buffers!(
    target::Vector{Float64},
    buffers::Vector{Vector{Float64}},
    nworkers::Int
)
    nworkers == 1 && return nothing

    _run_workers(nworkers) do worker
        reduce_csc_worker!(target, buffers, nworkers, worker)
    end

    return nothing
end

"""
    assemble_csc_worker!(ws, nzval, element_tags, connectivity, coordinates,
                         colptr, rowval, chunk, num_threads, worker,
                         num_nodes, num_integration_points, dim, edim,
                         grad_h, h, integration_weights, coefficient, weight,
                         op_u, Pu, op_s, Ps)

Assemble one worker's element subset directly into its private CSC value
buffer. All data used by the element loop are explicit arguments so Julia can
specialize the complete geometry, element-matrix, and scatter kernel without
capturing variables from `assemble_operator`.
"""
function assemble_csc_worker!(
    ws::OperatorWorkspace,
    nzval::Vector{Float64},
    element_tags::AbstractVector{<:Integer},
    connectivity::AbstractVector{<:Integer},
    coordinates::AbstractMatrix{<:Real},
    colptr::Vector{Int},
    rowval::Vector{Int},
    chunk::Int,
    num_threads::Int,
    worker::Int,
    num_nodes::Integer,
    num_integration_points::Integer,
    dim::Integer,
    edim::Integer,
    grad_h,
    h,
    integration_weights,
    coefficient,
    weight,
    op_u::OU,
    Pu::Problem,
    op_s::OS,
    Ps::Problem
) where {OU<:AbstractOp,OS<:AbstractOp}
    pdim_u = Pu.pdim
    pdim_s = Ps.pdim
    stride = chunk * num_threads

    for chunk_first in (1 + (worker - 1) * chunk):stride:length(element_tags)
        chunk_last = min(length(element_tags), chunk_first + chunk - 1)

        @inbounds for e in chunk_first:chunk_last
            elem = element_tags[e]
            first_node = (e - 1) * num_nodes

            compute_element_geometry!(
                ws,
                coordinates,
                connectivity,
                first_node,
                num_nodes,
                num_integration_points,
                dim,
                edim,
                grad_h
            )

            assemble_element_matrix_from_geometry!(
                ws,
                elem,
                num_integration_points,
                dim,
                edim,
                num_nodes,
                grad_h,
                h,
                integration_weights,
                coefficient,
                weight,
                op_u,
                Pu,
                op_s,
                Ps
            )

            Ke = ws.Ke

            # Columns belonging to the components of one trial node have
            # identical row patterns. Search once per node pair, then address
            # the complete component block by fixed offsets in nzval.
            for node_b in 1:num_nodes
                Jb_node = Int(connectivity[first_node + node_b])
                first_col = (Jb_node - 1) * pdim_u + 1
                first = colptr[first_col]
                last = colptr[first_col + 1] - 1
                column_nnz = last - first + 1
                b_loc0 = (node_b - 1) * pdim_u

                for node_a in 1:num_nodes
                    Ia_node = Int(connectivity[first_node + node_a])
                    first_row = (Ia_node - 1) * pdim_s + 1
                    block_position = find_csc_position(
                        rowval,
                        first,
                        last,
                        first_row
                    )
                    a_loc0 = (node_a - 1) * pdim_s

                    for comp_b in 1:pdim_u
                        p = block_position + (comp_b - 1) * column_nnz
                        b_loc = b_loc0 + comp_b

                        for comp_a in 1:pdim_s
                            a_loc = a_loc0 + comp_a
                            nzval[p + comp_a - 1] += Ke[a_loc, b_loc]
                        end
                    end
                end
            end
        end
    end

    return nothing
end

"""
    assemble_operator(Pu::Problem, 
                      op_u::AbstractOp, 
                      Ps::Problem, 
                      op_s::AbstractOp; 
                      coefficient::Union{Number,ScalarField}=1.0, 
                      weight=nothing, 
                      domain=nothing,
                      gauss=:full)

Assemble the sparse matrix corresponding to the bilinear form

    ∫ (Op_s v) ⋅ (Op_u u) dΩ

where

- `u` belongs to trial space `Pu`
- `v` belongs to test space `Ps`

# Arguments

`Pu`

Trial field problem.

`Ps`

Test field problem.

`op_u`

Operator applied to the trial field.

`op_s`

Operator applied to the test field.

# Keyword arguments

`coefficient`

Scalar or tensor coefficient multiplying the integrand.

Supported types

    Number
    ScalarField
    Matrix{Number}
    Matrix{ScalarField}
    coefficient chain containing matrices and scalar factors
        (Number or ScalarField)

The vector form represents a matrix chain

    A ⋅ B ⋅ C

which is evaluated at Gauss points as

    Cgp = A(x_gp) * B(x_gp) * C(x_gp)

`weight`

Optional scalar quadrature weight multiplying the integrand at each Gauss
point. It may be a `Number` or a `ScalarField`.

This is useful for geometry-related weights such as the axisymmetric factor
`2πr`:

    ∫(B' ⋅ D ⋅ B; Ω="body", weight=2π*r)

The weight is not treated as a matrix coefficient.

`domain`

Optional domain specification (`Ω` or `Γ`).

`gauss`

Number of integration points.
    `:full`       2order+1
    `:reduced`    2order-1
    Int           2order+1 + Int

`assembly=:csc`

Use direct CSC assembly. The CSC pattern is created first, element Jacobians
are computed on demand from nodal coordinates, and element matrices are
accumulated into one `nzval` array per worker. The first worker writes directly
to the final matrix; the remaining worker-local arrays are reduced into it.
With one worker, no private value buffer or reduction pass is required.

The same CSC path is used for scalar, matrix, and matrix-chain coefficients.
Matrices may contain `Number` and `ScalarField` entries; their values are
evaluated at Gauss points before the element matrix is scattered.

`K`

Optional prebuilt `SparseMatrixCSC{Float64,Int}` pattern. Its dimensions and
structural entries must match the current mesh, domain, and test/trial spaces.
Assembly adds to its existing `nzval` values; use
`fill!(K.nzval, 0.0)` before a new independent assembly.

`node_coordinates`

Optional reusable dense `3 × number_of_nodes` coordinate array for internal
compound assembly.

`threads`

Number of assembly workers. `:auto` uses all Julia threads, while `1` selects
the serial low-memory path.

`element_chunk_size`

Element scheduling block size for `assembly=:csc`. `:auto` uses at most 4096
elements and adapts to the element and worker counts. Chunks are distributed
cyclically among the requested workers.

`_nzval_buffers`

Internal reusable worker buffers used by compound CSC assembly. The first
buffer must be identical to `K.nzval`. This keyword is not part of the public
API.

# Returns

`SystemMatrix`

Sparse matrix representing the discretized bilinear form.
"""
function assemble_operator(
    Pu::Problem, op_u::AbstractOp,
    Ps::Problem, op_s::AbstractOp;
    coefficient::Union{Number,ScalarField,AbstractMatrix,AbstractVector}=1.0,
    weight=nothing,
    domain=nothing,
    gauss=:full,
    assembly::Symbol=:csc,
    threads=:auto,
    I=nothing,
    J=nothing,
    V=nothing,
    K=nothing,
    node_coordinates=nothing,
    _nzval_buffers=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto)

    gauss === :gll && error(
        "assemble_operator: gauss=:gll must be requested through ∫ for a " *
        "supported spectral bilinear form."
    )

    direct_csc = assembly === :csc
    num_threads = resolve_num_threads(threads)

    old_blas_threads = LinearAlgebra.BLAS.get_num_threads()
    if old_blas_threads != 1
        LinearAlgebra.BLAS.set_num_threads(1)
    end

    try
        assembly ∈ (:matrix, :triplets, :add, :csc) ||
        error("""
        assemble_operator: invalid assembly mode $assembly.
   
        Valid modes are:
            :matrix
            :triplets
            :add
            :csc
        """)
        @assert Pu.name == Ps.name "Both problems must refer to the same gmsh model/mesh."
        @assert Pu.dim == Ps.dim "Both problems must have the same spatial dimension."
        gmsh.model.setCurrent(Pu.name)
   
        # operator output dimensions
        out_u = op_outdim(op_u, Pu)
        out_s = op_outdim(op_s, Ps)
        
        if coefficient isa Number || coefficient isa ScalarField
            @assert out_u == out_s
        
        elseif coefficient isa AbstractMatrix
            if size(coefficient,1) != out_s || size(coefficient,2) != out_u
                error("""
                    Coefficient matrix size mismatch.
   
                    Test operator:
                        $(typeof(op_s))
                    output dimension = $out_s
   
                    Trial operator:
                        $(typeof(op_u))
                    output dimension = $out_u
   
                    Expected coefficient size:
                        ($out_s, $out_u)
   
                    Got:
                        $(size(coefficient))
   
                    Coefficient type:
                        $(typeof(coefficient))
                    """)
            end
        
        elseif coefficient isa AbstractVector

            isempty(coefficient) &&
                error("Coefficient chain must not be empty.")

            # Every factor must be either a scalar coefficient or a matrix.
            for C in coefficient
                (
                    C isa Number ||
                    C isa ScalarField ||
                    C isa AbstractMatrix
                ) ||
                error(
                    "Coefficient chain entries must be Number, ScalarField " *
                    "or AbstractMatrix, got $(typeof(C)).")
            end

            # Scalar factors do not change dimensions.
            matrices = [C for C in coefficient if C isa AbstractMatrix]

            if isempty(matrices)

                # Pure scalar chain:
                #   B_s' * c1 * c2 * ... * B_u
                out_s == out_u ||
                    error(
                        "Scalar coefficient chain requires equal test and trial " *
                        "operator dimensions: $out_s != $out_u."
                    )

            else

                A1 = first(matrices)
                An = last(matrices)

                size(A1, 1) == out_s ||
                    error(
                        "Coefficient chain size mismatch: first matrix has " *
                        "$(size(A1,1)) rows, expected $out_s."
                    )

                size(An, 2) == out_u ||
                    error(
                        "Coefficient chain size mismatch: last matrix has " *
                        "$(size(An,2)) columns, expected $out_u."
                    )

                # Scalar factors between matrices do not affect compatibility.
                for i in 1:length(matrices)-1
                    Ai = matrices[i]
                    Aj = matrices[i+1]

                    size(Ai, 2) == size(Aj, 1) ||
                        error(
                            "Coefficient chain matrix size mismatch: " *
                            "$(size(Ai)) cannot be multiplied by $(size(Aj))."
                        )
                end
            end
        end
#=
        elseif coefficient isa AbstractVector
            if length(coefficient) == 1 &&
               (coefficient[1] isa Number || coefficient[1] isa ScalarField)
        
                @assert out_u == out_s
        
            else
                A1 = coefficient[1]
                An = coefficient[end]
        
                @assert A1 isa AbstractMatrix
                @assert An isa AbstractMatrix
        
                if size(A1,1) != out_s || size(An,2) != out_u
                    error("""
                    Coefficient chain size mismatch.
   
                        Test operator:
                            $(typeof(op_s))
                        output dimension = $out_s
   
                        Trial operator:
                            $(typeof(op_u))
                        output dimension = $out_u
   
                        Expected:
                            first matrix rows = $out_s
                            last matrix cols  = $out_u
   
                        Got:
                            size(A1) = $(size(A1))
                            size(An) = $(size(An))
   
                        Coefficient chain types:
                            $(typeof.(coefficient))
                        """)
                end
        
                for i in 1:length(coefficient)-1
                    Ai = coefficient[i]
                    Aj = coefficient[i+1]
        
                    @assert Ai isa AbstractMatrix
                    @assert Aj isa AbstractMatrix
                    @assert size(Ai,2) == size(Aj,1)
                end
            end
        end
        =#
   
        # ------------------------------------------------------------
        # weight dimension check
        # ------------------------------------------------------------
        if weight !== nothing
            # optional: csak Number vagy ScalarField legyen
            (weight isa Number || weight isa ScalarField) ||
            error("weight must be Number or ScalarField")
        end
   
        # ------------------------------------------------------------
        # prepare coefficient
        # ------------------------------------------------------------
   
        # ------------------------------------------------------------
        # Assembly storage
        # ------------------------------------------------------------
   
        build_pattern = assembly !== :add && !direct_csc

        Kcsc = nothing
        coordinates = nothing
        nzval_buffers = nothing

        if direct_csc
            I === nothing ||
                error("assemble_operator: I must be nothing for assembly=:csc.")
            J === nothing ||
                error("assemble_operator: J must be nothing for assembly=:csc.")
            V === nothing ||
                error("assemble_operator: V must be nothing for assembly=:csc.")

            if K === nothing
                Kcsc = build_csc_pattern(Pu, Ps; domain=domain)
            else
                K isa SparseMatrixCSC{Float64,Int} || error(
                    "assemble_operator: assembly=:csc requires " *
                    "K::SparseMatrixCSC{Float64,Int}; got $(typeof(K))."
                )

                size(K) == (ndofs(Ps), ndofs(Pu)) || error(
                    "assemble_operator: supplied CSC matrix has size $(size(K)); " *
                    "expected $((ndofs(Ps), ndofs(Pu)))."
                )

                Kcsc = K
                _nzval_buffers === nothing &&
                    _validate_csc_component_pattern(Kcsc, Pu.pdim, Ps.pdim)
            end

            coordinates = node_coordinates === nothing ?
                dense_node_coordinates(Pu) : node_coordinates

            size(coordinates, 1) == 3 || error(
                "assemble_operator: node_coordinates must have three rows."
            )
            if _nzval_buffers === nothing
                nzval_buffers = Vector{Vector{Float64}}(undef, num_threads)
                nzval_buffers[1] = Kcsc.nzval
                for worker in 2:num_threads
                    nzval_buffers[worker] = zeros(Float64, length(Kcsc.nzval))
                end
            else
                _nzval_buffers isa Vector{Vector{Float64}} || error(
                    "assemble_operator: invalid internal CSC accumulator storage."
                )
                length(_nzval_buffers) == num_threads || error(
                    "assemble_operator: CSC accumulator count does not match threads."
                )
                _nzval_buffers[1] === Kcsc.nzval || error(
                    "assemble_operator: the first CSC accumulator must be K.nzval."
                )
                all(length(buffer) == length(Kcsc.nzval) for buffer in _nzval_buffers) ||
                    error("assemble_operator: CSC accumulator length mismatch.")
                nzval_buffers = _nzval_buffers
                for worker in 2:num_threads
                    fill!(nzval_buffers[worker], 0.0)
                end
            end
        elseif build_pattern
            I === nothing ||
                error("assemble_operator: I must be nothing for assembly=$assembly.")
   
            J === nothing ||
                error("assemble_operator: J must be nothing for assembly=$assembly.")
   
            V === nothing ||
                error("assemble_operator: V must be nothing for assembly=$assembly.")
   
            lengthOfIJV =
                LowLevelFEM.estimateLengthOfIJV(Pu) *
                max(1, Ps.pdim) *
                max(1, Pu.pdim)
   
            I = Vector{Int}(undef, lengthOfIJV)
            J = Vector{Int}(undef, lengthOfIJV)
            V = Vector{Float64}(undef, lengthOfIJV)
   
        else
            I === nothing ||
                error("assemble_operator: I must be nothing for assembly=:add.")
   
            J === nothing ||
                error("assemble_operator: J must be nothing for assembly=:add.")
   
            V isa Vector{Float64} ||
                error(
                    "assemble_operator: assembly=:add requires " *
                    "V::Vector{Float64}; got $(typeof(V))."
                )
        end
   
        pos = 1
   
        dim = Pu.dim
   
        # ------------------------------------------------------------
        # Domain/material selection FIX
        # ------------------------------------------------------------
        phNames = String[]
        domkind = nothing
   
        if domain === nothing
            # original behavior: loop over all materials (each has phName)
            phNames = [mat.phName for mat in Pu.material]
        elseif domain isa DomainSpec
            phNames = [domain.name]
            domkind = domain.kind   # :Ω or :Γ
        elseif domain isa AbstractString
            # (optional) accept plain String domain too
            phNames = [String(domain)]
            domkind = nothing
        else
            error("assemble_operator: unsupported domain type $(typeof(domain)). Expected nothing, DomainSpec or String.")
        end
   
        # loop physical groups
        for phName in phNames
            dom_here =
                domain === nothing ? DomainSpec(:Ω, phName) :
                domain isa DomainSpec ? domain :
                DomainSpec(:Ω, String(phName))
   
            Cprep = _prepare_coefficient(coefficient, dom_here)
            
            # Unwrap a coefficient chain containing a single factor.
            if Cprep isa AbstractVector && length(Cprep) == 1
                Cprep = only(Cprep)
            end
            
            Wprep = weight === nothing ? nothing : _prepare_coefficient(weight, dom_here)
   
            dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
   
            isempty(dimTags) && error("assemble_operator: physical group \"$phName\" not found.")
   
            for (edim, etag) in dimTags
   
                # optional dimension checks (matches your DSL design)
                if domkind === :Ω
                    edim == Pu.dim || error("Ω=\"$phName\" has dim=$edim but problem.dim=$(Pu.dim)")
                elseif domkind === :Γ
                    edim < Pu.dim || error("Γ=\"$phName\" has dim=$edim but expected < $(Pu.dim)")
                end
   
                if direct_csc
                    elemTypes = gmsh.model.mesh.getElementTypes(edim, etag)
                    elemTags = nothing
                    elemNodeTags = nothing
                else
                    elemTypes, elemTags, elemNodeTags =
                        gmsh.model.mesh.getElements(edim, etag)
                end
   
                for itype in eachindex(elemTypes)
                    et = elemTypes[itype]
                    _, _, order, numNodes::Int64, _, _ = gmsh.model.mesh.getElementProperties(et)

                    if direct_csc
                        element_tags, connectivity =
                            gmsh.model.mesh.getElementsByType(et, etag)
                    else
                        element_tags = elemTags[itype]
                        connectivity = elemNodeTags[itype]
                    end
   
                    gorder = gauss == :reduced ? max(1, 2order - 1) :
                        gauss == :full    ? (2order + 1) :
                        gauss isa Int     ? max(1, 2order+1 + gauss) :
                        (2order + 1)
   
                    intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(gorder))
                    #intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                    numIntPoints = length(intWeights)
   
                    _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                    h = reshape(fun, :, numIntPoints)
   
                    _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                    ∇h = reshape(dfun, :, numIntPoints)
   
                    if !direct_csc
                        # Legacy IJV path: keep the original batched Gmsh
                        # geometry retrieval until that path is retired.
                        elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)
                        jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)

                        numElementsGlobal = length(elemTags_global)
                        nip = length(det_all) ÷ numElementsGlobal
                        @assert nip == numIntPoints

                        TagType = eltype(elemTags_global)
                        idxmap = Dict{TagType,Int}()
                        sizehint!(idxmap, numElementsGlobal)

                        @inbounds for (gi, gtag) in enumerate(elemTags_global)
                            idxmap[gtag] = gi - 1
                        end
                    end
   
                    # buffers
                    nel = length(element_tags)
   
                    ndofs_u_loc = Pu.pdim * numNodes
                    ndofs_s_loc = Ps.pdim * numNodes
   
                    # The CSC path uses the flat Gmsh connectivity directly.
                    # The legacy IJV path retains its old dense table for now.
                    nnet = direct_csc ? nothing : zeros(Int, nel, numNodes)
                    if !direct_csc
                        @inbounds for e in 1:nel
                            for a in 1:numNodes
                                nnet[e, a] = connectivity[(e-1)*numNodes+a]
                            end
                        end
                    end
   
                    entries_per_element = ndofs_s_loc * ndofs_u_loc
                    block_start = pos
                    block_entries = nel * entries_per_element
   
                    if !direct_csc && build_pattern && block_start + block_entries - 1 > length(V)
                        newlen = max(
                            block_start + block_entries - 1,
                            Int(ceil(1.5 * length(V))) + 1000
                        )
   
                        resize!(I, newlen)
                        resize!(J, newlen)
                        resize!(V, newlen)
                    elseif !direct_csc && !build_pattern && block_start + block_entries - 1 > length(V)
                        error("""
                        assemble_operator: assembly pattern mismatch.
                    
                        The current element block requires more entries than the existing
                        V vector contains.
                        """)
                    end
   
                    # `Threads.threadid()` is a global thread identifier, not
                    # an index local to the default thread pool. In notebook
                    # environments the task that enters `@threads` may run on
                    # an `:interactive` thread whose ID is larger than
                    # `Threads.nthreads()`. The static IJV path indexes by the
                    # actual thread ID, so provide a workspace for every
                    # possible global thread ID. Logical-worker paths continue
                    # to need only `num_threads` workspaces.
                    workspace_count =
                        !direct_csc && threads === :auto && num_threads > 1 ?
                        Threads.maxthreadid() : num_threads

                    workspaces = [
                        OperatorWorkspace(
                            dim,
                            numNodes,
                            numIntPoints,
                            out_u,
                            out_s,
                            ndofs_u_loc,
                            ndofs_s_loc
                        )
                        for _ in 1:workspace_count
                    ]
                    
                    if direct_csc
                        colptr = Kcsc.colptr
                        rowval = Kcsc.rowval

                        chunk = element_chunk_size === :auto ?
                            max(1, min(4096, cld(nel, num_threads))) :
                            element_chunk_size isa Integer ? Int(element_chunk_size) :
                            error(
                                "assemble_operator: element_chunk_size must be " *
                                ":auto or a positive integer."
                            )
                        chunk > 0 || error(
                            "assemble_operator: element_chunk_size must be positive."
                        )

                        _run_workers(num_threads) do worker
                            assemble_csc_worker!(
                                workspaces[worker],
                                nzval_buffers[worker],
                                element_tags,
                                connectivity,
                                coordinates,
                                colptr,
                                rowval,
                                chunk,
                                num_threads,
                                worker,
                                numNodes,
                                numIntPoints,
                                dim,
                                edim,
                                ∇h,
                                h,
                                intWeights,
                                Cprep,
                                Wprep,
                                op_u,
                                Pu,
                                op_s,
                                Ps
                            )
                        end
                    else
                        Vthread = V::Vector{Float64}
                        Ithread = build_pattern ? (I::Vector{Int}) : Int[]
                        Jthread = build_pattern ? (J::Vector{Int}) : Int[]
   
                        let Ithread = Ithread,
                        Jthread = Jthread,
                        Vthread = Vthread,
                        build_pattern = build_pattern
   
                        if threads === :auto && num_threads > 1
                            # Restore the original v1.13.28 scheduling for the
                            # default multithreaded IJV path. Static scheduling
                            # keeps each thread on a fixed element subset and
                            # associates its workspace with the actual thread.
                            Threads.@threads :static for e in 1:nel
                                ws = workspaces[Threads.threadid()]
                                elem = elemTags[itype][e]
                                gidx = idxmap[elem]

                                assemble_element_matrix!(
                                    ws,
                                    elem,
                                    gidx,
                                    jac_all,
                                    det_all,
                                    nip,
                                    numIntPoints,
                                    dim,
                                    edim,
                                    numNodes,
                                    ∇h,
                                    h,
                                    intWeights,
                                    Cprep,
                                    Wprep,
                                    op_u,
                                    Pu,
                                    op_s,
                                    Ps
                                )

                                Ke = ws.Ke

                                element_start =
                                    block_start + (e - 1) * entries_per_element

                                # scatter Ke(s,u) -> global IJV
                                @inbounds for a_loc in 1:ndofs_s_loc
                                    node_a = div(a_loc - 1, Ps.pdim) + 1
                                    comp_a = mod(a_loc - 1, Ps.pdim) + 1
                                    Ia_node = nnet[e, node_a]
                                    Ia = (Ia_node - 1) * Ps.pdim + comp_a

                                    row_start =
                                        element_start + (a_loc - 1) * ndofs_u_loc

                                    @inbounds for b_loc in 1:ndofs_u_loc
                                        node_b = div(b_loc - 1, Pu.pdim) + 1
                                        comp_b = mod(b_loc - 1, Pu.pdim) + 1
                                        Jb_node = nnet[e, node_b]
                                        Jb = (Jb_node - 1) * Pu.pdim + comp_b

                                        p = row_start + b_loc - 1

                                        if build_pattern
                                            Ithread[p] = Ia
                                            Jthread[p] = Jb
                                            Vthread[p] = Ke[a_loc, b_loc]
                                        else
                                            Vthread[p] += Ke[a_loc, b_loc]
                                        end
                                    end
                                end
                            end
                        else
                            # Keep the common function barrier for serial and
                            # explicitly limited worker counts. With one worker,
                            # _run_workers calls the closure directly and creates
                            # no task.
                            _run_workers(num_threads) do worker
                                ws = workspaces[worker]
                                for e in _worker_range(nel, num_threads, worker)
                                    elem = elemTags[itype][e]
                                    gidx = idxmap[elem]

                                    assemble_element_matrix!(
                                        ws,
                                        elem,
                                        gidx,
                                        jac_all,
                                        det_all,
                                        nip,
                                        numIntPoints,
                                        dim,
                                        edim,
                                        numNodes,
                                        ∇h,
                                        h,
                                        intWeights,
                                        Cprep,
                                        Wprep,
                                        op_u,
                                        Pu,
                                        op_s,
                                        Ps
                                    )

                                    Ke = ws.Ke

                                    element_start =
                                        block_start + (e - 1) * entries_per_element

                                    # scatter Ke(s,u) -> global IJV
                                    @inbounds for a_loc in 1:ndofs_s_loc
                                        node_a = div(a_loc - 1, Ps.pdim) + 1
                                        comp_a = mod(a_loc - 1, Ps.pdim) + 1
                                        Ia_node = nnet[e, node_a]
                                        Ia = (Ia_node - 1) * Ps.pdim + comp_a

                                        row_start =
                                            element_start + (a_loc - 1) * ndofs_u_loc

                                        @inbounds for b_loc in 1:ndofs_u_loc
                                            node_b = div(b_loc - 1, Pu.pdim) + 1
                                            comp_b = mod(b_loc - 1, Pu.pdim) + 1
                                            Jb_node = nnet[e, node_b]
                                            Jb = (Jb_node - 1) * Pu.pdim + comp_b

                                            p = row_start + b_loc - 1

                                            if build_pattern
                                                Ithread[p] = Ia
                                                Jthread[p] = Jb
                                                Vthread[p] = Ke[a_loc, b_loc]
                                            else
                                                Vthread[p] += Ke[a_loc, b_loc]
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        end
                    end
                    pos += block_entries
                end
            end
        end
   
        nentries = pos - 1
        nrows = ndofs(Ps)
        ncols = ndofs(Pu)
        
        @assert nrows > 0
        @assert ncols > 0

        if direct_csc
            Kcsc_typed = Kcsc::SparseMatrixCSC{Float64,Int}
            nzval_buffers_typed =
                nzval_buffers::Vector{Vector{Float64}}

            reduce_csc_buffers!(
                Kcsc_typed.nzval,
                nzval_buffers_typed,
                num_threads
            )

            return SystemMatrix(Kcsc_typed, Pu, Ps)
        end
        
        if assembly === :add
            nentries == length(V) ||
                error("""
                    assemble_operator: assembly pattern mismatch.
        
                The existing V vector contains:
                    $(length(V)) entries
        
                The current operator generated:
                    $nentries entries
        
                All compound terms must traverse exactly the same
                elements and local matrix entries in the same order.
                """)
        
            return V
        end
        
        resize!(I, nentries)
        resize!(J, nentries)
        resize!(V, nentries)
        
        if nentries > 0
            @assert maximum(I) <= nrows
            @assert maximum(J) <= ncols
            @assert minimum(I) >= 1
            @assert minimum(J) >= 1
        end
        
        if assembly === :triplets
            return I, J, V, nrows, ncols
        end
        
        K = sparse(I, J, V, nrows, ncols)
        dropzeros!(K)
   
        return SystemMatrix(K, Pu, Ps)
    finally
        if old_blas_threads != 1
            LinearAlgebra.BLAS.set_num_threads(old_blas_threads)
        end
    end
end

