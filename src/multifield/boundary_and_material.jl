"""
    multifield_bc_data(K::SystemMatrix,
                       bc::Vector{BoundaryCondition};
                       nsteps::Int=1)

Compute global constrained and free DOFs together with the prescribed
Dirichlet values for a multifield block system.

Returns
- `free  :: Vector{Int}`
- `fixed :: Vector{Int}`
- `xD    :: Matrix{Float64}` of size `(ndof, nsteps)`

Notes
- This function requires `K` to be a block `SystemMatrix`, i.e.
  `K.problems !== nothing` and `K.offsets !== nothing`.
- In multifield problems every `BoundaryCondition` must explicitly
  refer to a `problem`.
- Prescribed values are assembled from `applyBoundaryConditions(problem, [bc]; steps=nsteps)`.
"""
function multifield_bc_data(
    K::SystemMatrix,
    bc::Vector{BoundaryCondition};
    nsteps::Int=1)

    K.problems === nothing &&
        error("multifield_bc_data: K is not a block SystemMatrix (K.problems === nothing).")

    K.offsets === nothing &&
        error("multifield_bc_data: K is not a block SystemMatrix (K.offsets === nothing).")

    problems = K.problems
    offsets = K.offsets

    length(problems) == length(offsets) ||
        error("multifield_bc_data: inconsistent metadata in K (length(problems) != length(offsets)).")

    ndof = size(K.A, 1)
    size(K.A, 1) == size(K.A, 2) ||
        error("multifield_bc_data: K must be square.")

    fixed = Int[]
    xD = zeros(ndof, nsteps)

    for bci in bc
        bci.problem === nothing &&
            error("multifield_bc_data: in multifield systems every BoundaryCondition must have an explicit `problem`.")

        idx = findfirst(P -> P === bci.problem, problems)
        idx === nothing &&
            error("multifield_bc_data: BC refers to a Problem that is not present in K.problems.")

        P = problems[idx]
        off = offsets[idx]

        local_fixed = constrainedDoFs(P, [bci])
        isempty(local_fixed) && continue

        global_fixed = off .+ local_fixed
        append!(fixed, global_fixed)

        field_bc = applyBoundaryConditions(P, [bci]; steps=nsteps)

        # safety checks
        size(field_bc.a, 1) == ndofs(P) ||
            error("multifield_bc_data: applyBoundaryConditions returned incompatible field size for problem $(P.field).")

        if field_bc.nsteps == 1
            @inbounds for it in 1:nsteps
                xD[global_fixed, it] .= field_bc.a[local_fixed, 1]
            end
        else
            field_bc.nsteps == nsteps ||
                error("multifield_bc_data: BC field nsteps = $(field_bc.nsteps), expected 1 or $nsteps.")
            xD[global_fixed, :] .= field_bc.a[local_fixed, 1:nsteps]
        end
    end

    fixed = unique(fixed)
    sort!(fixed)

    #free = setdiff(collect(1:ndof), fixed)
    active = Int[]

    for (i, P) in enumerate(problems)
        off = offsets[i]
        append!(active, off .+ allDoFs(P))
    end

    sort!(active)
    unique!(active)

    free = setdiff(active, fixed)

    return free, fixed, xD
end

"""
    multifield_constrainedDoFs(K::SystemMatrix,
                               bc::Vector{BoundaryCondition})

Return global constrained DOFs for a multifield block system.
"""
function multifield_constrainedDoFs(
    K::SystemMatrix,
    bc::Vector{BoundaryCondition}
)
    _, fixed, _ = multifield_bc_data(K, bc; nsteps=1)
    return fixed
end


"""
    multifield_freeDoFs(K::SystemMatrix,
                        bc::Vector{BoundaryCondition})

Return global free DOFs for a multifield block system.
"""
function multifield_freeDoFs(
    K::SystemMatrix,
    bc::Vector{BoundaryCondition}
)
    free, _, _ = multifield_bc_data(K, bc; nsteps=1)
    return free
end

"""
    constrainedDoFs(K::SystemMatrix, 
                    support::Vector{BoundaryCondition})

Return global constrained DOFs for single- or multi-field systems.
"""
function constrainedDoFs(
    K::SystemMatrix,
    support::Vector{BoundaryCondition}
    )
    if K.problems === nothing
        return constrainedDoFs(K.model, support)
    else
        _, fixed, _ = multifield_bc_data(K, support; nsteps=1)
        return fixed
    end
end

"""
    freeDoFs(K::SystemMatrix, 
              support::Vector{BoundaryCondition})

Return global free DOFs for single- or multi-field systems.
"""
function freeDoFs(
    K::SystemMatrix,
    support::Vector{BoundaryCondition}
    )
    if K.problems === nothing
        return freeDoFs(K.model, support)
    else
        free, _, _ = multifield_bc_data(K, support; nsteps=1)
        return free
    end
end

function active_identity(P::Problem)

    active = allDoFs(P)

    n = ndofs(P)

    I = sparse(
        active,
        1:length(active),
        ones(Float64, length(active)),
        n,
        length(active)
    )

    R = sparse(
        1:length(active),
        active,
        ones(Float64, length(active)),
        length(active),
        n
    )

    return I, R
end

"""
    consistentToLumped(M::SystemMatrix) -> SystemMatrix

Converts a **consistent mass (or capacity) matrix** into a **lumped (diagonal) matrix**
using the row-sum technique.

Each diagonal entry is computed as
```

M_ii = ∑_j M_ij

```

The resulting matrix is diagonal and preserves the total mass associated with each DOF.

---

### Arguments
- `M::SystemMatrix`  
  Consistent mass (or capacity) matrix.

---

### Returns
- `Ml::SystemMatrix`  
  Lumped (diagonal) matrix with the same structure (`model`, `test_model`,
  `problems`, `offsets`) as the input.

---

### Notes
- This is the standard **row-sum lumping** technique.
- The method is commonly used for **explicit time integration schemes**
  (e.g. central difference), where a diagonal mass matrix is required.
- Off-diagonal entries are discarded, and their contributions are accumulated
  onto the diagonal.
- The operation preserves the total mass but modifies the spectral properties
  of the system.

---

### Performance
The implementation avoids dense row access and works efficiently with sparse matrices.

"""
function consistentToLumped(M::SystemMatrix)
    n = size(M.A, 1)

    # Row sums (efficient for sparse matrices)
    d = vec(sum(M.A, dims=2))

    # Build diagonal sparse matrix
    A = spdiagm(d)

    return SystemMatrix(A, M.model, M.test_model, M.problems, M.offsets)
end

"""
    rigidRotationMap(mpc_u::MPC, mpc_φ::MPC)

Construct the rigid-body rotational kinematic map associated with a remote
point constraint.

The function creates the displacement contribution generated by rotations
about a single master node. The displacement and rotation MPCs must refer to
the same master and slave physical groups.

The returned matrix represents the relation

```math
u_{rot} = R \\varphi,
```

where `R` is the rigid-rotation map.

The full displacement field of the slave region can therefore be written as

# Arguments

- `mpc_u::MPC`:  
  MPC associated with the translational displacement field.

- `mpc_φ::MPC`:  
  MPC associated with the rotational field.

Both MPCs must specify their fields and must use identical `master` and  
`slave` physical groups.

# Returns

- `SystemMatrix`:  
  Kinematic transformation matrix mapping the rotation field to the  
  displacement field.
  
  The trial field is the rotation field and the test field is the  
  displacement field.

# Supported Configurations

- 2D:  
  displacement `VectorField` with two components and scalar rotation field.

- 3D:  
  displacement `VectorField` with three components and rotation  
  `VectorField` with three components.

The master physical group must contain exactly one node.

# Example

```julia
mpc_u = MPC(master="remote", slave="right", field=U)
mpc_φ = MPC(master="remote", slave="right", field=Φ)

R = rigidRotationMap(mpc_u, mpc_φ)

Kuφ = Ku * R
Kφ  = R' * Ku * R

K = SystemMatrix([
    Ku    Kuφ
    Kuφ'  Kφ
])

u, φ = solveField(K, F; support=[bc], mpc=[mpc_u, mpc_φ])

u_phys = u + R * φ
```

# Notes

The component selections stored in the MPCs are respected when constructing  
the map.  
"""
function rigidRotationMap(mpc_u::MPC, mpc_φ::MPC)

    # ----------------------------------------------------------
    # Basic compatibility checks
    # ----------------------------------------------------------

    mpc_u.problem !== nothing ||
        error(
            "rigidRotationMap: displacement field must be specified in the MPC."
        )

    mpc_φ.problem !== nothing ||
        error(
            "rigidRotationMap: rotation field must be specified in the MPC."
        )

    mpc_u.master == mpc_φ.master ||
        error(
            "rigidRotationMap: displacement and rotation MPCs must have the same master."
        )

    mpc_u.slave == mpc_φ.slave ||
        error(
            "rigidRotationMap: displacement and rotation MPCs must have the same slave."
        )

    U = mpc_u.problem
    Φ = mpc_φ.problem

    U.name == Φ.name ||
        error(
            "rigidRotationMap: displacement and rotation fields must belong to the same Gmsh model."
        )

    U.dim == Φ.dim ||
        error(
            "rigidRotationMap: displacement and rotation fields must have the same spatial dimension."
        )

    is2D = U.dim == 2 && U.pdim == 2 && Φ.pdim == 1
    is3D = U.dim == 3 && U.pdim == 3 && Φ.pdim == 3

    (is2D || is3D) ||
        error(
            "rigidRotationMap: expected either a 2D VectorField/ScalarField " *
            "or a 3D VectorField/VectorField pair."
        )

    gmsh.model.setCurrent(U.name)

    # ----------------------------------------------------------
    # Active MPC components
    # ----------------------------------------------------------

    function active_components(mpc::MPC, problem::Problem)

        active = trues(problem.pdim)

        comp_map = _bc_component_map(problem)
        prefix = String(problem.field)

        for (sym, value) in mpc.components

            s = String(sym)

            startswith(s, prefix) ||
                error(
                    "rigidRotationMap: unknown MPC component '$sym' " *
                    "for field $(problem.field)."
                )

            # Unicode-safe, important e.g. for fieldName=:φ
            comp = String(chopprefix(s, prefix))

            haskey(comp_map, comp) ||
                error(
                    "rigidRotationMap: invalid component '$comp' " *
                    "for field $(problem.field)."
                )

            active[comp_map[comp]] = value
        end

        return active
    end

    active_u = active_components(mpc_u, U)
    active_φ = active_components(mpc_φ, Φ)

    # ----------------------------------------------------------
    # Master and slave nodes
    # ----------------------------------------------------------

    masterTag = getTagForPhysicalName(mpc_u.master)
    slaveTag  = getTagForPhysicalName(mpc_u.slave)

    masterNodes, masterCoord =
        gmsh.model.mesh.getNodesForPhysicalGroup(-1, masterTag)

    slaveNodes, slaveCoord =
        gmsh.model.mesh.getNodesForPhysicalGroup(-1, slaveTag)

    length(masterNodes) == 1 ||
        error(
            "rigidRotationMap: master physical group '$(mpc_u.master)' " *
            "must contain exactly one node."
        )

    masterNode = Int(masterNodes[1])

    xm = masterCoord[1]
    ym = masterCoord[2]
    zm = masterCoord[3]

    # ----------------------------------------------------------
    # Sparse rigid-rotation map
    # ----------------------------------------------------------

    I = Int[]
    J = Int[]
    V = Float64[]

    if is2D

        φ_active = active_φ[1]

        if φ_active
            for (k, node) in enumerate(slaveNodes)

                slaveNode = Int(node)

                slaveNode == masterNode && continue

                x = slaveCoord[3k - 2]
                y = slaveCoord[3k - 1]

                dx = x - xm
                dy = y - ym

                φdof = slaveNode

                # ux = -dy * φ
                if active_u[1] && !iszero(dy)

                    uxdof =
                        U.pdim * (slaveNode - 1) + 1

                    push!(I, uxdof)
                    push!(J, φdof)
                    push!(V, -dy)
                end

                # uy = +dx * φ
                if active_u[2] && !iszero(dx)

                    uydof =
                        U.pdim * (slaveNode - 1) + 2

                    push!(I, uydof)
                    push!(J, φdof)
                    push!(V, dx)
                end
            end
        end

    else

        for (k, node) in enumerate(slaveNodes)

            slaveNode = Int(node)

            slaveNode == masterNode && continue

            x = slaveCoord[3k - 2]
            y = slaveCoord[3k - 1]
            z = slaveCoord[3k]

            dx = x - xm
            dy = y - ym
            dz = z - zm

            # Local map:
            #
            #       φx    φy    φz
            #
            # ux     0    dz   -dy
            # uy   -dz     0    dx
            # uz    dy   -dx     0

            Rd = (
                ( 0.0,  dz,  -dy),
                (-dz,  0.0,  dx),
                ( dy, -dx,  0.0)
            )

            for i in 1:3

                active_u[i] || continue

                udof =
                    U.pdim * (slaveNode - 1) + i

                for j in 1:3

                    active_φ[j] || continue

                    value = Rd[i][j]

                    iszero(value) && continue

                    φdof =
                        Φ.pdim * (slaveNode - 1) + j

                    push!(I, udof)
                    push!(J, φdof)
                    push!(V, value)
                end
            end
        end
    end

    RA = sparse(
        I,
        J,
        V,
        ndofs(U),
        ndofs(Φ)
    )

    # R : Φ -> U
    return SystemMatrix(
        RA,
        Φ,
        U
    )
end

"""
    collapseMPC(field, mpc)

Collapse generalized nodal quantities from MPC slave degrees of freedom onto
their corresponding master degrees of freedom.

The operation applies the transpose of the MPC kinematic transformation to
the supplied nodal field. It is primarily intended for generalized forces,
reactions and moments associated with constrained degrees of freedom.

For an MPC transformation

```math
q = T q_r,
```

the corresponding generalized force transformation is

`collapseMPC` performs this force-side reduction and returns the result in the  
original field representation, with the collapsed values assigned to the  
master degrees of freedom.

# Arguments

- `field::Union{ScalarField,VectorField,TensorField}`:  
  Nodal field containing generalized forces or other quantities to be  
  collapsed.

- `mpc::MPC`:  
  Multi-point constraint defining the master-slave relation.

# Returns

A field of the same type as `field`, containing the collapsed generalized  
quantities.

# Example

For a rotational field in a coupled beam formulation:

```julia
mφ = Kuφ' * u + Kφ * φ
M = collapseMPC(mφ, mpc_φ)
```

The resulting field `M` contains the resultant generalized moments transferred  
from the constrained slave degrees of freedom to their MPC master degrees of  
freedom.

# Notes

This function is intended for nodal generalized quantities. It should not be  
used as an elementwise integration operation.  
"""
function collapseMPC(
    field::Union{ScalarField,VectorField,TensorField},
    mpc::MPC
    )

    problem = field.model

    rep = mpcRepresentativeMap(problem, [mpc])
    T, p, _ = mpcTransformation(rep)

    f = elementsToNodes(field)

    fred = T' * f.a

    out = copy(f)
    fill!(out.a, 0.0)

    out.a[p, :] .= fred

    return out
end

"""
    ConstitutiveMatrix(type::Symbol)
    ConstitutiveMatrix(type::Symbol, mat::Material)
    ConstitutiveMatrix(type::Symbol, E, ν)

Return the isotropic elastic constitutive matrix for the selected model.

Supported models:
- `:PlaneStress`
- `:PlaneStrain`
- `:Axisymmetric`
- `:Solid`

Engineering shear strain convention is used.

# Alias

`D` is a shorthand alias for `ConstitutiveMatrix` and accepts the same
arguments.

```julia
D(:Solid, mat)
D(:PlaneStress, E, ν)
```

"""
function ConstitutiveMatrix(
    type::Symbol,
    E = 2.1e5,
    ν = 0.3
)
    if type === :PlaneStress

        c = E / (1 - ν^2)

        return [
            c      c*ν    0
            c*ν    c      0
            0      0      c*(1-ν)/2
        ]

    elseif type === :PlaneStrain

        λ = E * ν / ((1 + ν) * (1 - 2ν))
        μ = E / (2 * (1 + ν))

        return [
            λ + 2μ   λ         0
            λ        λ + 2μ    0
            0        0         μ
        ]

    elseif type === :Axisymmetric

        λ = E * ν / ((1 + ν) * (1 - 2ν))
        μ = E / (2 * (1 + ν))

        return [
            λ + 2μ   λ         λ         0
            λ        λ + 2μ    λ         0
            λ        λ         λ + 2μ    0
            0        0         0         μ
        ]

    elseif type === :Solid

        λ = E * ν / ((1 + ν) * (1 - 2ν))
        μ = E / (2 * (1 + ν))

        return [
            λ + 2μ   λ         λ         0   0   0
            λ        λ + 2μ    λ         0   0   0
            λ        λ         λ + 2μ    0   0   0
            0        0         0         μ   0   0
            0        0         0         0   μ   0
            0        0         0         0   0   μ
        ]

    else
        error(
            "ConstitutiveMatrix: unsupported type '$type'. " *
            "Supported types are :PlaneStress, :PlaneStrain, " *
            ":Axisymmetric and :Solid."
        )
    end
end

ConstitutiveMatrix(type::Symbol, mat::Material) =
    ConstitutiveMatrix(type, mat.E, mat.ν)

"""
    D(args...)

Shorthand alias for [`ConstitutiveMatrix`](@ref).

Accepts the same arguments and returns the same constitutive matrix.

# Examples

```julia
D(:Solid, mat)
D(:PlaneStress, 210e3, 0.3)
```

"""  
D(args...) = ConstitutiveMatrix(args...)

#=
"""
    mpcReducedBCData(full_to_reduced, fixed, uD::AbstractMatrix)

Map a time history of prescribed full-space values to the MPC-reduced space.

The set of constrained reduced DOFs must remain identical for all time steps.
"""
function mpcReducedBCData(
    full_to_reduced,
    fixed,
    uD::AbstractMatrix
    )

    nsteps = size(uD, 2)

    free_r, fixed_r, uD1 =
        mpcReducedBCData(
            full_to_reduced,
            fixed,
            @view(uD[:, 1])
        )

    UDr = zeros(Float64, length(uD1), nsteps)
    UDr[:, 1] .= uD1

    for i in 2:nsteps

        free_i, fixed_i, uDi =
            mpcReducedBCData(
                full_to_reduced,
                fixed,
                @view(uD[:, i])
            )

        free_i == free_r ||
            error(
                "MPC: reduced free DOFs change between time steps."
            )

        fixed_i == fixed_r ||
            error(
                "MPC: reduced constrained DOFs change between time steps."
            )

        UDr[:, i] .= uDi
    end

    return free_r, fixed_r, UDr
end
=#

"""
    _multifield_mpc_bc_data(K, mpcs, fixed, xD::AbstractMatrix)

Map a time history of prescribed multifield values through MPC relations.
"""
function _multifield_mpc_bc_data(
    K::SystemMatrix,
    mpcs::Vector{MPC},
    fixed::Vector{Int},
    xD::AbstractMatrix
    )

    nsteps = size(xD, 2)
    ndof = size(xD, 1)

    fixed_new, xD1 =
        _multifield_mpc_bc_data(
            K,
            mpcs,
            fixed,
            @view(xD[:, 1])
        )

    XDnew = zeros(Float64, ndof, nsteps)
    XDnew[:, 1] .= xD1

    for i in 2:nsteps

        fixed_i, xDi =
            _multifield_mpc_bc_data(
                K,
                mpcs,
                fixed,
                @view(xD[:, i])
            )

        fixed_i == fixed_new ||
            error(
                "MPC: constrained multifield DOFs change between time steps."
            )

        XDnew[:, i] .= xDi
    end

    return fixed_new, XDnew
end

"""
    _multifield_mpc_transformation(K, mpcs)

Construct the global multifield transformation including reduced-order
interpolation and equality MPC constraints.

Returns

    T, R, prunable, protected
"""
function _multifield_mpc_transformation(
    K::SystemMatrix,
    mpcs::Vector{MPC}
    )

    K.problems === nothing &&
        error(
            "_multifield_mpc_transformation: K must be a block SystemMatrix."
        )

    problems = K.problems
    offsets = K.offsets

    for c in mpcs

        c.problem === nothing &&
            error(
                "MPC: in a multifield system every MPC must have an explicit field."
            )

        any(P -> P === c.problem, problems) ||
            error(
                "MPC: constraint refers to a field that is not present in the system."
            )
    end

    Tblocks = SparseMatrixCSC{Float64,Int}[]
    Rblocks = SparseMatrixCSC{Float64,Int}[]

    prunable = Int[]
    master_full_global = Int[]

    final_offset = 0

    for (i, P) in enumerate(problems)

        local_mpcs =
            [
                c for c in mpcs
                if c.problem === P
            ]

        Tp, Rp =
            _base_reduction(
                P,
                local_mpcs
            )

        Tp, Rp =
            _apply_mpc_to_reduction(
                Tp,
                Rp,
                P,
                local_mpcs
            )

        push!(Tblocks, Tp)
        push!(Rblocks, Rp)

        if !isempty(local_mpcs)

            append!(
                prunable,
                final_offset .+
                collect(1:size(Tp, 2))
            )

            local_master =
                _mpc_master_dofs(
                    P,
                    local_mpcs
                )

            append!(
                master_full_global,
                offsets[i] .+ local_master
            )
        end

        final_offset += size(Tp, 2)
    end

    T = blockdiag(Tblocks...)
    R = blockdiag(Rblocks...)

    protected = Int[]

    if !isempty(master_full_global)

        RR = R[:, master_full_global]

        I, _, V = findnz(RR)

        for k in eachindex(I)
            iszero(V[k]) && continue
            push!(protected, I[k])
        end

        sort!(unique!(protected))
    end

    return T, R, prunable, protected
end

function _remove_inactive_mpc_dofs(
    matrices::Tuple,
    rhs,
    free::Vector{Int},
    fixed::Vector{Int},
    protected::Vector{Int},
    prunable::Vector{Int}
    )

    n = size(first(matrices), 1)
    active = falses(n)

    for A in matrices

        rows = rowvals(A)
        vals = nonzeros(A)

        for j in axes(A, 2)
            for k in nzrange(A, j)

                iszero(vals[k]) && continue

                active[j] = true
                active[rows[k]] = true
            end
        end
    end

    if rhs !== nothing

        for i in 1:n
            if any(!iszero, @view rhs[i, :])
                active[i] = true
            end
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

