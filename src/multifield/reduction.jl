struct SymmetricSystemMatrix
    parent::SystemMatrix
    uplo::Symbol
end

function LinearAlgebra.Symmetric(K::SystemMatrix, uplo::Symbol=:U)
    uplo in (:U, :L) ||
        throw(ArgumentError("uplo must be :U or :L"))
    return SymmetricSystemMatrix(K, uplo)
end

"""
    _canonicalize_mpc_rep!(rep)

Resolve chained MPC relations so that every constrained DOF points directly
to its final master DOF.

Cycles are rejected.
"""
function _canonicalize_mpc_rep!(rep::Vector{Int})

    n = length(rep)

    for i in 1:n
        j = i
        path = Int[]

        while rep[j] != j

            j in path &&
                error("MPC: cyclic master-slave relation detected.")

            push!(path, j)
            j = rep[j]
        end

        root = j

        for k in path
            rep[k] = root
        end

        rep[i] = root
    end

    return rep
end

"""
    _mpc_master_dofs(P, mpcs)

Return all DOFs belonging to MPC master nodes of field `P`.

A master physical group may contain either a single node or multiple nodes,
as in periodic constraints. All components are returned, including components
explicitly disabled in the MPC. This preserves unconstrained master DOFs as
genuine system DOFs.
"""
function _mpc_master_dofs(
    P::Problem,
    mpcs::Vector{MPC}
    )

    gmsh.model.setCurrent(P.name)

    dofs = Int[]

    for c in mpcs

        if c.problem !== nothing && c.problem !== P
            continue
        end
        #c.problem === P || continue

        tag =
            getTagForPhysicalName(c.master)

        nodes, _ =
            gmsh.model.mesh.getNodesForPhysicalGroup(
                -1,
                tag
            )

        isempty(nodes) &&
            error(
                "MPC: master physical group " *
                "'$(c.master)' contains no nodes."
            )

        for node0 in nodes

            node = Int(node0)

            for comp in 1:P.pdim

                push!(
                    dofs,
                    P.pdim * (node - 1) + comp
                )
            end
        end
    end

    sort!(unique!(dofs))

    return dofs
end

#=
"""
    _mpc_master_dofs(P, mpcs)

Return all DOFs belonging to MPC master nodes of field `P`.

All components are returned, including components explicitly disabled in
the MPC. This preserves unconstrained master DOFs as genuine system DOFs.
"""
function _mpc_master_dofs(
    P::Problem,
    mpcs::Vector{MPC}
    )

    gmsh.model.setCurrent(P.name)

    dofs = Int[]

    for c in mpcs

        c.problem === P || continue

        tag = getTagForPhysicalName(c.master)

        nodes, _ =
            gmsh.model.mesh.getNodesForPhysicalGroup(-1, tag)

        length(nodes) == 1 ||
            error(
                "MPC: master physical group '$(c.master)' " *
                "must contain exactly one node."
            )

        node = Int(nodes[1])

        for comp in 1:P.pdim
            push!(
                dofs,
                P.pdim * (node - 1) + comp
            )
        end
    end

    sort!(unique!(dofs))

    return dofs
end
=#

"""
    _represented_rows(T)

Return a Boolean vector marking full-space DOFs represented by `T`.
"""
function _represented_rows(T::SparseMatrixCSC)

    represented = falses(size(T, 1))

    rows = rowvals(T)
    vals = nonzeros(T)

    for j in axes(T, 2)
        for k in nzrange(T, j)

            iszero(vals[k]) && continue

            represented[rows[k]] = true
        end
    end

    return represented
end

"""
    _initialize_reduced_approximation!(P)

Build and cache the intrinsic `p -> p-1` Lagrange transformation of a
`reducedOrder=true` field. Active full-space DOFs that are not represented by
`reductionMatrices(P)` (for example remote reference points) are appended once
as identity coordinates, so solve-specific transformations only need to deal
with MPC-dependent additions.
"""
function _initialize_reduced_approximation!(P::Problem)
    P.basis === :lagrange ||
        error(
            "Internal error: reduced-order cache currently expects " *
            "basis=:lagrange."
        )

    P.reducedOrder ||
        error(
            "Internal error: reduced-order cache requested for " *
            "reducedOrder=false."
        )

    T, R = reductionMatrices(P)

    represented = _represented_rows(T)
    wanted = sort!(unique(Int.(allDoFs(P))))

    extra = Int[
        i for i in wanted
        if !represented[i]
    ]

    if !isempty(extra)
        nfull = ndofs(P)
        nextra = length(extra)

        Te = sparse(
            extra,
            1:nextra,
            ones(Float64, nextra),
            nfull,
            nextra
        )

        Re = sparse(
            1:nextra,
            extra,
            ones(Float64, nextra),
            nextra,
            nfull
        )

        T = hcat(T, Te)
        R = vcat(R, Re)
    end

    P.approximation.T = T
    P.approximation.R = R
    P.approximation.metadata = nothing

    return P.approximation
end

"""
    _base_reduction(P, mpcs)

Construct the field-level base transformation before applying MPCs.

For full-order fields, active physical DOFs and MPC master DOFs are retained.

For reduced-order fields, the polynomial-order reduction is used and
additional active DOFs not represented by that reduction (e.g. remote
reference points) are appended as identity DOFs.
"""
function _base_reduction(
    P::Problem,
    mpcs::Vector{MPC}
    )

    if P.basis === :spectral
        P.reducedOrder &&
            error(
                "Approximation: basis=:spectral combined with " *
                "reducedOrder=true is not yet implemented."
            )

        isempty(mpcs) ||
            error(
                "Approximation: spectral basis combined with MPCs " *
                "is not yet implemented."
            )

        return _approximation_matrices(P)
    end

    master_dofs = _mpc_master_dofs(P, mpcs)

    if !P.reducedOrder

        active =
            sort!(
                unique(
                    vcat(
                        allDoFs(P),
                        master_dofs
                    )
                )
            )

        nfull = ndofs(P)
        nred = length(active)

        T = sparse(
            active,
            1:nred,
            ones(Float64, nred),
            nfull,
            nred
        )

        R = sparse(
            1:nred,
            active,
            ones(Float64, nred),
            nred,
            nfull
        )

        return T, R
    end

    # Polynomial-order reduction is intrinsic to the Problem and cached.
    T, R = _approximation_matrices(P)

    # Remote points and other lower-dimensional active DOFs are not
    # necessarily represented by reductionMatrices().
    represented = _represented_rows(T)

    wanted =
        sort!(
            unique(
                vcat(
                    allDoFs(P),
                    master_dofs
                )
            )
        )

    extra =
        [
            i for i in wanted
            if !represented[i]
        ]

    isempty(extra) && return T, R

    nfull = ndofs(P)
    nextra = length(extra)

    Te = sparse(
        extra,
        1:nextra,
        ones(Float64, nextra),
        nfull,
        nextra
    )

    Re = sparse(
        1:nextra,
        extra,
        ones(Float64, nextra),
        nextra,
        nfull
    )

    return hcat(T, Te), vcat(R, Re)
end

"""
    _apply_mpc_to_reduction(T, R, P, mpcs)

Apply equality MPC constraints to an existing field transformation.

The constraints are imposed in the original full nodal space,

    u_slave = u_master

while the unknowns may already belong to a reduced-order basis.

Returns the modified prolongation and restriction matrices.
"""
function _apply_mpc_to_reduction(
    T::SparseMatrixCSC{Float64,Int},
    R::SparseMatrixCSC{Float64,Int},
    P::Problem,
    mpcs::Vector{MPC}
    )

    isempty(mpcs) && return T, R

    rep = mpcRepresentativeMap(P, mpcs)

    _canonicalize_mpc_rep!(rep)

    # ------------------------------------------------------
    # Fast path for full-order fields
    # ------------------------------------------------------

    if !P.reducedOrder

        # DOFs represented by the active material space
        represented = _represented_rows(T)
        active = findall(represented)

        # Final MPC representatives of all active DOFs
        roots =
            sort!(
                unique(
                    rep[active]
                )
            )

        # Every representative must itself belong to the
        # active space. Remote masters are already added by
        # _base_reduction().
        for root in roots
            represented[root] ||
                error(
                    "MPC: representative DOF $root is outside " *
                    "the active field space."
                )
        end

        # full representative DOF -> reduced DOF
        root_to_reduced =
            zeros(
                Int,
                length(rep)
            )

        for (j, root) in enumerate(roots)
            root_to_reduced[root] = j
        end

        # Every active full-space DOF gets exactly one
        # reduced coordinate.
        J =
            Vector{Int}(
                undef,
                length(active)
            )

        @inbounds for k in eachindex(active)

            dof = active[k]

            J[k] =
                root_to_reduced[
                    rep[dof]
                ]

            J[k] != 0 ||
                error(
                    "MPC: failed to map DOF $dof " *
                    "to the reduced field space."
                )
        end

        nfull = size(T, 1)
        nred = length(roots)

        # Prolongation:
        #
        #     u = Tnew * ur
        #
        Tnew =
            sparse(
                active,
                J,
                ones(Float64, length(active)),
                nfull,
                nred
            )

        # Restriction selects the representative/master value.
        Rnew =
            sparse(
                1:nred,
                roots,
                ones(Float64, nred),
                nred,
                nfull
            )

        return Tnew, Rnew
    end

    represented = _represented_rows(T)

    constraints =
        [
            (i, rep[i])
            for i in eachindex(rep)
            if rep[i] != i &&
            represented[i]
        ]
    
    for (slave, master) in constraints

        represented[master] ||
            error(
                "MPC: master DOF $master of slave DOF $slave " *
                "is not represented in the active field space."
            )
    end

    #constraints =
    #    [
    #        (i, rep[i])
    #        for i in eachindex(rep)
    #        if rep[i] != i
    #    ]

    isempty(constraints) && return T, R

    # Final root masters are protected from elimination whenever possible.
    master_full =
        unique(last(c) for c in constraints)

    for (slave, master) in constraints

        # Current constraint in the current reduced coordinates:
        #
        #     (T_slave - T_master) * q = 0
        #
        c =
            vec(
                Array(
                    T[slave, :] - T[master, :]
                )
            )

        scale = maximum(abs, c; init=0.0)
        tol = 100eps(Float64) * max(1.0, scale)

        nz =
            findall(x -> abs(x) > tol, c)

        # Constraint already follows from previous constraints.
        isempty(nz) && continue

        # ------------------------------------------------------
        # Protect current master coordinates whenever possible.
        # ------------------------------------------------------

        protected = falses(size(T, 2))

        for mdof in master_full

            row = vec(Array(T[mdof, :]))

            for j in eachindex(row)
                if abs(row[j]) > tol
                    protected[j] = true
                end
            end
        end

        candidates =
            [
                j for j in nz
                if !protected[j]
            ]

        # This can occur for chained constraints.
        isempty(candidates) &&
            (candidates = nz)

        pivot =
            candidates[
                argmax(abs.(c[candidates]))
            ]

        cp = c[pivot]

        abs(cp) > tol ||
            error("MPC: failed to select a stable elimination pivot.")

        nold = size(T, 2)

        keep =
            [
                j for j in 1:nold
                if j != pivot
            ]

        nnew = nold - 1

        I = Int[]
        J = Int[]
        V = Float64[]

        sizehint!(I, 2nnew)
        sizehint!(J, 2nnew)
        sizehint!(V, 2nnew)

        for (newj, oldj) in enumerate(keep)

            # Identity part
            push!(I, oldj)
            push!(J, newj)
            push!(V, 1.0)

            # Eliminated pivot:
            #
            # q_p = -Σ c_j/c_p q_j
            val = -c[oldj] / cp

            if abs(val) > tol
                push!(I, pivot)
                push!(J, newj)
                push!(V, val)
            end
        end

        E = sparse(
            I,
            J,
            V,
            nold,
            nnew
        )

        # Natural restriction corresponding to E.
        S = sparse(
            1:nnew,
            keep,
            ones(Float64, nnew),
            nnew,
            nold
        )

        T = T * E
        R = S * R
    end

    return T, R
end

"""
    _singlefield_mpc_transformation(P, mpcs)

Construct the active field transformation and apply equality MPCs.
"""
function _singlefield_mpc_transformation(
    P::Problem,
    mpcs::Vector{MPC}
    )

    T, R =
        _base_reduction(
            P,
            mpcs
        )

    T, R =
        _apply_mpc_to_reduction(
            T,
            R,
            P,
            mpcs
        )

    return T, R
end

function _singlefield_mpc_bc_data(
    P::Problem,
    mpcs::Vector{MPC},
    fixed::Vector{Int},
    uD::AbstractVector
    )

    rep =
        mpcRepresentativeMap(
            P,
            mpcs
        )

    _canonicalize_mpc_rep!(rep)

    fixed_new = Int[]
    uD_new = zeros(Float64, length(uD))

    assigned =
        Dict{Int,Float64}()

    for dof in fixed

        master = rep[dof]
        value = uD[dof]

        if haskey(assigned, master)

            isapprox(
                assigned[master],
                value
            ) ||
                error(
                    "MPC: incompatible prescribed values on tied DOFs."
                )

        else

            assigned[master] = value
            push!(fixed_new, master)
            uD_new[master] = value
        end
    end

    sort!(unique!(fixed_new))

    return fixed_new, uD_new
end

"""
    _multifield_mpc_bc_data(K, mpcs, fixed, xD)

Move prescribed values on tied slave DOFs to their final master DOFs.

Conflicting prescribed values belonging to the same MPC DOF are rejected.
"""
function _multifield_mpc_bc_data(
    K::SystemMatrix,
    mpcs::Vector{MPC},
    fixed::Vector{Int},
    xD::AbstractVector
    )

    problems = K.problems
    offsets = K.offsets

    ndof = size(K.A, 1)

    global_rep = collect(1:ndof)

    for (i, P) in enumerate(problems)

        local_mpcs =
            [
                c for c in mpcs
                if c.problem === P
            ]

        isempty(local_mpcs) && continue

        rep =
            mpcRepresentativeMap(
                P,
                local_mpcs
            )

        _canonicalize_mpc_rep!(rep)

        off = offsets[i]

        for j in eachindex(rep)
            global_rep[off + j] =
                off + rep[j]
        end
    end

    fixed_new = Int[]
    xD_new = zeros(Float64, ndof)

    assigned = Dict{Int,Float64}()

    for dof in fixed

        master = global_rep[dof]
        value = xD[dof]

        if haskey(assigned, master)

            isapprox(
                assigned[master],
                value
            ) ||
                error(
                    "MPC: incompatible prescribed values on tied DOFs."
                )

        else

            assigned[master] = value
            push!(fixed_new, master)
            xD_new[master] = value
        end
    end

    sort!(unique!(fixed_new))

    return fixed_new, xD_new
end

function _remove_inactive_mpc_dofs(
    Kr::SparseMatrixCSC,
    Br::AbstractMatrix,
    free::Vector{Int},
    fixed::Vector{Int},
    protected::Vector{Int},
    prunable::Vector{Int}
    )

    n = size(Kr, 1)

    active = falses(n)

    rows = rowvals(Kr)
    vals = nonzeros(Kr)

    for j in axes(Kr, 2)
        for k in nzrange(Kr, j)

            iszero(vals[k]) && continue

            active[j] = true
            active[rows[k]] = true
        end
    end

    for i in 1:n
        if any(!iszero, @view Br[i, :])
            active[i] = true
        end
    end

    active[fixed] .= true
    active[protected] .= true

    can_prune = falses(n)
    can_prune[prunable] .= true

    return [
        i for i in free
        if !can_prune[i] || active[i]
    ]
end



###############################################################
