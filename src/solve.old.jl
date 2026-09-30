###############################################################################
#                                                                             #
#                              solve.jl                                       #
#                                                                             #
#  Linear finite-element solution pipeline.                                  #
#                                                                             #
#      prepare  ->  solve_linear_system  ->  reconstruct                      #
#                                                                             #
#  Supports single-field and multifield systems, reduced-order interpolation, #
#  MPCs, local coordinate systems, multiple right-hand sides, and direct or   #
#  iterative linear solvers.                                                 #
#                                                                             #
###############################################################################

using LinearAlgebra
using SparseArrays
using IterativeSolvers

export Global, solveField

const _SolveFieldType = Union{ScalarField,VectorField,TensorField}

"""
    GlobalField{F}

Internal wrapper marking right-hand-side data as expressed in the global
Cartesian coordinate system. Construct it through [`Global`](@ref).
"""
struct GlobalField{F}
    field::F
end

"""
    SolveRHS{F}

Internal symbolic sum of right-hand-side data expressed partly in local solver
coordinates and partly in global Cartesian coordinates.
"""
struct SolveRHS{F}
    local_part::Union{Nothing,F}
    global_part::Union{Nothing,F}
end

"""
    NodalCoordinateSystem(e1, e2=nothing)

Local coordinate-system definition based on direction `VectorField`s.

`e1` defines the first local basis direction. In 3D, `e2` supplies a second
reference direction and is orthogonalized against `e1`. In 2D, `e2` is
optional; the second basis direction is obtained by a 90-degree rotation of
`e1`, while an optional `e2` can select its sign.

The support of the direction fields determines the nodes on which the local
coordinate system is active.
"""
struct NodalCoordinateSystem
    e1::VectorField
    e2::Union{VectorField,Nothing}
end

"""
    Global(f)

Mark right-hand-side data as expressed in the global Cartesian coordinate
system.

When `coordSys` is supplied to [`solveField`](@ref), ordinary unwrapped loads
are interpreted directly in the local solver components. `Global(f)` is the
explicit exception: its Cartesian components are transformed to the local
solver basis before the linear system is solved.

`Global` accepts scalar, vector and tensor fields as well as multifield
`SystemVector`s. It can also be combined with local load terms, for example

```julia
f = f_local + Global(f_gravity)
```
"""
Global(f::_SolveFieldType) = GlobalField(f)
Global(F::SystemVector) = GlobalField(F)

# VectorField-based constructors for local nodal coordinate systems.
CoordinateSystem(e1::VectorField) = NodalCoordinateSystem(e1, nothing)
CoordinateSystem(e1::VectorField, e2::VectorField) = NodalCoordinateSystem(e1, e2)


# =============================================================================
# RHS expression handling
# =============================================================================

function _check_rhs_compatibility(
    a::SystemVector,
    b::SystemVector
)
    a.problems == b.problems ||
        error("solveField: multifield RHS Problem ordering mismatch.")

    a.offsets == b.offsets ||
        error("solveField: multifield RHS offsets are incompatible.")

    size(a.a) == size(b.a) ||
        error("solveField: multifield RHS vectors have incompatible sizes.")

    return nothing
end

function _merge_rhs_fields(
    a::SystemVector,
    b::SystemVector
)
    _check_rhs_compatibility(a, b)

    return SystemVector(
        a.a + b.a,
        a.model,
        a.problems,
        a.offsets
    )
end

function _check_rhs_compatibility(a, b)
    typeof(a) === typeof(b) ||
        error("solveField: local and global load terms must have the same field type.")

    a.model === b.model ||
        error("solveField: local and global load terms must belong to the same Problem.")

    size(a.a) == size(b.a) ||
        error("solveField: local and global load terms have incompatible sizes.")

    return nothing
end

_merge_rhs_fields(::Nothing, b) = b
_merge_rhs_fields(a, ::Nothing) = a

function _merge_rhs_fields(a, b)
    _check_rhs_compatibility(a, b)
    return a + b
end

function Base.:+(a::_SolveFieldType, b::GlobalField)
    _check_rhs_compatibility(a, b.field)
    return SolveRHS(a, b.field)
end

Base.:+(a::GlobalField, b::_SolveFieldType) = b + a

function Base.:+(a::GlobalField, b::GlobalField)
    _check_rhs_compatibility(a.field, b.field)
    return SolveRHS(nothing, a.field + b.field)
end

function Base.:+(
    a::SystemVector,
    b::GlobalField{<:SystemVector}
)
    _check_rhs_compatibility(a, b.field)
    return SolveRHS(a, b.field)
end

Base.:+(
    a::GlobalField{<:SystemVector},
    b::SystemVector
) = b + a

function Base.:+(
    a::GlobalField{<:SystemVector},
    b::GlobalField{<:SystemVector}
)
    return SolveRHS(
        nothing,
        _merge_rhs_fields(a.field, b.field)
    )
end

function Base.:+(
    a::SolveRHS{SystemVector},
    b::SystemVector
)
    return SolveRHS(
        _merge_rhs_fields(a.local_part, b),
        a.global_part
    )
end

Base.:+(
    a::SystemVector,
    b::SolveRHS{SystemVector}
) = b + a

function Base.:+(a::SolveRHS, b::_SolveFieldType)
    return SolveRHS(
        _merge_rhs_fields(a.local_part, b),
        a.global_part
    )
end

Base.:+(a::_SolveFieldType, b::SolveRHS) = b + a

function Base.:+(a::SolveRHS, b::GlobalField)
    return SolveRHS(
        a.local_part,
        _merge_rhs_fields(a.global_part, b.field)
    )
end

Base.:+(a::GlobalField, b::SolveRHS) = b + a

function Base.:+(a::SolveRHS, b::SolveRHS)
    return SolveRHS(
        _merge_rhs_fields(a.local_part, b.local_part),
        _merge_rhs_fields(a.global_part, b.global_part)
    )
end

_rhs_parts(f::_SolveFieldType) = (f, nothing)
_rhs_parts(F::SystemVector) = (F, nothing)
_rhs_parts(f::GlobalField) = (nothing, f.field)
_rhs_parts(f::SolveRHS) = (f.local_part, f.global_part)

function _rhs_template(rhs)
    local_part, global_part = _rhs_parts(rhs)

    if local_part !== nothing
        return local_part
    elseif global_part !== nothing
        return global_part
    else
        error("solveField: empty right-hand side.")
    end
end


# =============================================================================
# Linear solver backend
# =============================================================================

"""
    solve_linear_system(A, B; solver=:backslash, preconditioner=nothing,
                        reltol=sqrt(eps()), maxiter=size(A, 1))

Solve `A * X = B` for one or multiple right-hand sides stored columnwise.

Supported solver symbols are `:backslash`, `:auto_symmetric`, `:lu`,
`:lu_no_ordering`, `:cholesky`, `:qr`, `:cg`, and `:gmres`.

`:auto_symmetric` first attempts a Cholesky factorization. If the symmetric
matrix is not positive definite, it automatically falls back to sparse LU.
"""
function solve_linear_system(
    A,
    B;
    solver=:backslash,
    preconditioner=nothing,
    reltol=sqrt(eps()),
    maxiter=size(A, 1)
)
    if solver === :backslash
        return A \ B

    elseif solver === :auto_symmetric
        try
            F = cholesky(A)
            return F \ B
        catch err
            if err isa PosDefException
                @info "Symmetric system is not positive definite; using LU factorization."
                F = lu(sparse(A))
                return F \ B
            end
            rethrow()
        end

    elseif solver === :lu
        A0 = A isa Union{Symmetric,Hermitian} ? sparse(A) : A
        F = lu(A0)
        return F \ B

    elseif solver === :lu_no_ordering
        A0 = A isa Union{Symmetric,Hermitian} ? sparse(A) : A
        F = lu(A0, q=nothing)
        return F \ B

    elseif solver === :cholesky
        F = cholesky(A)
        return F \ B

    elseif solver === :qr
        F = qr(A)
        return F \ B

    elseif solver === :cg
        X = similar(B)

        for j in axes(B, 2)
            if preconditioner === nothing
                X[:, j] = cg(
                    A,
                    @view(B[:, j]);
                    reltol=reltol,
                    maxiter=maxiter
                )
            else
                X[:, j] = cg(
                    A,
                    @view(B[:, j]);
                    Pl=preconditioner,
                    reltol=reltol,
                    maxiter=maxiter
                )
            end
        end

        return X

    elseif solver === :gmres
        X = similar(B)

        for j in axes(B, 2)
            if preconditioner === nothing
                X[:, j] = gmres(
                    A,
                    @view(B[:, j]);
                    reltol=reltol,
                    maxiter=maxiter
                )
            else
                X[:, j] = gmres(
                    A,
                    @view(B[:, j]);
                    Pl=preconditioner,
                    reltol=reltol,
                    maxiter=maxiter
                )
            end
        end

        return X

    else
        error("Unknown linear solver: $solver.")
    end
end

"""
    default_solver(K)

Return the default linear solver associated with the matrix wrapper.
"""
default_solver(::SystemMatrix) = :backslash
default_solver(::SymmetricSystemMatrix) = :auto_symmetric


# =============================================================================
# Kinematic reduction and boundary-condition helpers
# =============================================================================

"""
    field_transformation(P, mpc=MPC[])

Return prolongation and restriction matrices containing the full kinematic
transformation of a single field, including reduced-order interpolation and
multi-point constraints.
"""
function field_transformation(
    P::Problem,
    mpc::Vector{MPC}=MPC[]
)
    return _singlefield_mpc_transformation(P, mpc)
end

"""
    singlefield_mpc_bc_data_matrix(P, mpc, fixed, xD)

Transfer prescribed values from MPC slave DOFs to their final master DOFs for
one or multiple right-hand sides.
"""
function singlefield_mpc_bc_data_matrix(
    P::Problem,
    mpc::Vector{MPC},
    fixed::Vector{Int},
    xD::AbstractMatrix
)
    isempty(mpc) && return fixed, xD

    nsteps = size(xD, 2)
    xD_eff = zeros(eltype(xD), size(xD))
    fixed_eff = Int[]

    for j in 1:nsteps
        fixed_j, xD_j = _singlefield_mpc_bc_data(
            P,
            mpc,
            fixed,
            @view xD[:, j]
        )

        if j == 1
            fixed_eff = fixed_j
        elseif fixed_j != fixed_eff
            error("MPC: inconsistent constrained DOFs between right-hand sides.")
        end

        xD_eff[:, j] .= xD_j
    end

    return fixed_eff, xD_eff
end

"""
    reduced_bc_data_matrix(R, fixed, xD)

Map full-space Dirichlet data with one or multiple right-hand sides to the
reduced space.
"""
function reduced_bc_data_matrix(
    R::SparseMatrixCSC,
    fixed::AbstractVector{<:Integer},
    xD::AbstractMatrix
)
    nr = size(R, 1)

    if isempty(fixed)
        return (
            collect(1:nr),
            Int[],
            zeros(eltype(xD), nr, size(xD, 2))
        )
    end

    free_r, fixed_r, _ = reduced_bc_data(
        R,
        fixed,
        @view xD[:, 1]
    )

    xD_r = R * xD

    return free_r, fixed_r, xD_r
end

parent_system(K::SystemMatrix) = K
parent_system(K::SymmetricSystemMatrix) = K.parent

reduced_matrix_wrapper(::SystemMatrix, A) = A
reduced_matrix_wrapper(K::SymmetricSystemMatrix, A) = Symmetric(A, K.uplo)


# =============================================================================
# VectorField-based coordinate systems
# =============================================================================

function _coordinate_field_support_nodes(v::VectorField)
    P = v.model
    gmsh.model.setCurrent(P.name)

    if isElementwise(v)
        nodes = Int[]

        for elem in v.numElem
            _, nodeTags, _, _ = gmsh.model.mesh.getElement(elem)
            append!(nodes, Int.(nodeTags))
        end

        sort!(unique!(nodes))
        return nodes
    end

    ncomp = v.type == :v2D ? 2 : 3
    nodes = Int[]
    tol = 100 * eps(Float64)

    for node in 1:P.non
        i0 = (node - 1) * ncomp
        value = @view v.a[i0 + 1:i0 + ncomp, 1]
        norm(value) > tol && push!(nodes, node)
    end

    return nodes
end

function _coordinate_nodal_field(v::VectorField)
    v.nsteps == 1 ||
        error("CoordinateSystem: direction fields must contain exactly one step.")

    return isNodal(v) ? v : elementsToNodes(v)
end

function _coordinate_direction(
    vn::VectorField,
    node::Int,
    dim::Int
)
    isNodal(vn) ||
        error("CoordinateSystem: internal direction field must be nodal.")

    ncomp = vn.type == :v2D ? 2 : 3

    ncomp >= dim ||
        error("CoordinateSystem: direction field has too few components.")

    i0 = (node - 1) * ncomp
    return Vector{Float64}(vn.a[i0 + 1:i0 + dim, 1])
end

"""
    _coordinate_transformation(P, cs)

Build a nodal `Transformation` from a VectorField-based coordinate-system
definition. The returned matrix maps local nodal components to global nodal
components.
"""
function _coordinate_transformation(
    P::Problem,
    cs::NodalCoordinateSystem
)
    P.pdim == P.dim ||
        error("solveField: coordSys requires a vector field with pdim == dim.")

    P.dim in (2, 3) ||
        error("solveField: coordSys currently supports only 2D and 3D vector fields.")

    cs.e1.model === P ||
        error("CoordinateSystem: e1 belongs to a different Problem.")

    if cs.e2 !== nothing
        cs.e2.model === P ||
            error("CoordinateSystem: e2 belongs to a different Problem.")
    end

    nodes1 = _coordinate_field_support_nodes(cs.e1)
    isempty(nodes1) && error("CoordinateSystem: e1 has empty support.")

    e1n = _coordinate_nodal_field(cs.e1)
    e2n = nothing

    if cs.e2 !== nothing
        nodes2 = _coordinate_field_support_nodes(cs.e2)
        nodes1 == nodes2 ||
            error("CoordinateSystem: e1 and e2 must have identical nodal support.")
        e2n = _coordinate_nodal_field(cs.e2)
    elseif P.dim == 3
        error("CoordinateSystem: a second direction field is required in 3D.")
    end

    nd = ndofs(P)
    active_dofs = falses(nd)

    for node in nodes1
        i0 = (node - 1) * P.dim
        active_dofs[i0 + 1:i0 + P.dim] .= true
    end

    I = Int[]
    J = Int[]
    V = Float64[]

    sizehint!(I, nd + length(nodes1) * P.dim^2)
    sizehint!(J, nd + length(nodes1) * P.dim^2)
    sizehint!(V, nd + length(nodes1) * P.dim^2)

    for dof in 1:nd
        if !active_dofs[dof]
            push!(I, dof)
            push!(J, dof)
            push!(V, 1.0)
        end
    end

    tol = sqrt(eps(Float64))

    for node in nodes1
        e1 = _coordinate_direction(e1n, node, P.dim)
        n1 = norm(e1)
        n1 > tol ||
            error("CoordinateSystem: zero first direction at node $node.")
        e1 ./= n1

        if P.dim == 2
            e2 = [-e1[2], e1[1]]

            if cs.e2 !== nothing
                e2ref = _coordinate_direction(e2n, node, P.dim)
                if norm(e2ref) > tol && dot(e2, e2ref) < 0
                    e2 .*= -1
                end
            end

            Qnode = hcat(e1, e2)
        else
            e2raw = _coordinate_direction(e2n, node, P.dim)
            e2 = e2raw - dot(e1, e2raw) * e1
            n2 = norm(e2)
            n2 > tol ||
                error("CoordinateSystem: first and second directions are parallel at node $node.")
            e2 ./= n2

            e3 = cross(e1, e2)
            n3 = norm(e3)
            n3 > tol ||
                error("CoordinateSystem: degenerate local basis at node $node.")
            e3 ./= n3

            Qnode = hcat(e1, e2, e3)
        end

        i0 = (node - 1) * P.dim

        for j in 1:P.dim
            for i in 1:P.dim
                value = Qnode[i, j]
                if value != 0.0
                    push!(I, i0 + i)
                    push!(J, i0 + j)
                    push!(V, value)
                end
            end
        end
    end

    Q = sparse(I, J, V, nd, nd)
    dropzeros!(Q)

    return Transformation(Q, P.non, P.dim), nodes1
end

"""
    _build_coordinate_transformation(P, coordSys)

Build the combined nodal coordinate transformation. Coordinate systems must
have disjoint nodal support.
"""
function _build_coordinate_transformation(
    P::Problem,
    coordSys
)
    isempty(coordSys) && return nothing

    P.pdim == P.dim ||
        error("solveField: coordSys requires a vector field with pdim == dim.")

    Q = Transformation(
        spdiagm(0 => ones(Float64, ndofs(P))),
        P.non,
        P.dim
    )

    used_nodes = Set{Int}()

    for cs in coordSys
        cs isa NodalCoordinateSystem ||
            error(
                "solveField: coordSys entries must be created from VectorField " *
                "directions, e.g. CoordinateSystem(e1, e2)."
            )

        Qi, nodes = _coordinate_transformation(P, cs)
        overlap = intersect(used_nodes, Set(nodes))

        isempty(overlap) ||
            error(
                "solveField: overlapping coordinate systems are not supported yet; " *
                "$(length(overlap)) node(s) are assigned more than once."
            )

        union!(used_nodes, nodes)
        Q = Q * Qi
    end

    return Q
end

"""
    rhs_in_solver_coordinates(rhs, Q)

Return a template field and the RHS matrix expressed in solver coordinates.
Ordinary RHS terms are local. Terms wrapped in `Global(...)` are transformed
from global Cartesian components by `Q'`.
"""
function rhs_in_solver_coordinates(
    rhs,
    Q
)
    local_part, global_part = _rhs_parts(rhs)
    template = _rhs_template(rhs)

    if local_part !== nothing && global_part !== nothing
        _check_rhs_compatibility(local_part, global_part)
    end

    B = zeros(eltype(template.a), size(template.a))

    if local_part !== nothing
        B .+= local_part.a
    end

    if global_part !== nothing
        if Q === nothing
            B .+= global_part.a
        else
            B .+= Q.T' * global_part.a
        end
    end

    return template, B
end


# =============================================================================
# Single-field system preparation and reconstruction
# =============================================================================

"""
    _prepare_singlefield_system(K, rhs, support; mpc=MPC[], coordSys=[])

Prepare a single-field finite-element system for the common solve pipeline.

With a coordinate transformation `Q` and a kinematic transformation `T`,

    u_global = Q * T * u_reduced

and the reduced system is formed as

    K_reduced = T' * Q' * K * Q * T.

`coordSys` and MPC may be combined. In that case the MPC identifies local
components and the coordinate transformation gives those components their
global physical directions. Combining `coordSys` with reduced-order
interpolation on the same field is not yet implemented.
"""
function _prepare_singlefield_system(
    K0,
    rhs,
    support;
    mpc::Vector{MPC}=MPC[],
    coordSys=NodalCoordinateSystem[]
)
    K = parent_system(K0)
    P = K.model

    K.problems === nothing ||
        error("prepare_singlefield_system: K must be a single-field SystemMatrix.")

    P === nothing &&
        error("prepare_singlefield_system: K has no associated Problem.")

    size(K.A, 1) == size(K.A, 2) ||
        error("prepare_singlefield_system: system matrix must be square.")

    template = _rhs_template(rhs)

    template.model === P ||
        error("prepare_singlefield_system: RHS and system matrix belong to different Problems.")

    size(template.a, 1) == size(K.A, 1) ||
        error("prepare_singlefield_system: incompatible matrix and RHS sizes.")

    Q = _build_coordinate_transformation(P, coordSys)

    if Q !== nothing
        P.reducedOrder &&
            error("solveField: CoordinateSystem + reducedOrder=true on the same field is not yet implemented.")

        template isa VectorField ||
            error("solveField: coordSys currently supports VectorField problems only.")
    end

    rhs_template, F0 = rhs_in_solver_coordinates(rhs, Q)
    nsteps = size(F0, 2)

    A0 = Q === nothing ? K.A : Q.T' * K.A * Q.T

    fixed = constrainedDoFs(P, support)
    xDfield = applyBoundaryConditions(P, support; steps=nsteps)
    xD = xDfield.a

    T, R = field_transformation(P, mpc)

    fixed_eff, xD_eff = singlefield_mpc_bc_data_matrix(
        P,
        mpc,
        fixed,
        xD
    )

    free_r, fixed_r, xD_r = reduced_bc_data_matrix(
        R,
        fixed_eff,
        xD_eff
    )

    Kr = T' * A0 * T
    Br = T' * F0

    B = copy(Br[free_r, :])

    if !isempty(fixed_r)
        B .-= Kr[free_r, fixed_r] * xD_r[fixed_r, :]
    end

    A = reduced_matrix_wrapper(K0, Kr[free_r, free_r])

    return (
        A=A,
        B=B,
        T=T,
        Q=Q,
        template=rhs_template,
        xD_r=xD_r,
        free_r=free_r,
        fixed_r=fixed_r,
        nr=size(T, 2)
    )
end

"""
    prepare_singlefield_system(K, rhs, support; mpc=MPC[], coordSys=[])

Prepare a single-field system without solving it.
"""
function prepare_singlefield_system(
    K::Union{SystemMatrix,SymmetricSystemMatrix},
    rhs::Union{_SolveFieldType,GlobalField,SolveRHS},
    support::Vector{BoundaryCondition};
    mpc::Vector{MPC}=MPC[],
    coordSys=NodalCoordinateSystem[]
)
    return _prepare_singlefield_system(
        K,
        rhs,
        support;
        mpc=mpc,
        coordSys=coordSys
    )
end

"""
    reconstruct_singlefield_solution(Xfree, prepared)

Reconstruct the full nodal finite-element field from the reduced free-DOF
solution and return it in global Cartesian components.
"""
function reconstruct_singlefield_solution(
    Xfree,
    prepared
)
    nsteps = size(Xfree, 2)
    xr = zeros(eltype(Xfree), prepared.nr, nsteps)

    xr[prepared.free_r, :] .= Xfree

    if !isempty(prepared.fixed_r)
        xr[prepared.fixed_r, :] .= prepared.xD_r[prepared.fixed_r, :]
    end

    x_local = prepared.T * xr
    x_global = prepared.Q === nothing ? x_local : prepared.Q.T * x_local

    u = copy(prepared.template)
    u.a .= x_global

    return u
end



# =============================================================================
# Multifield coordinate-system helpers
# =============================================================================

"""
    _build_multifield_coordinate_transformation(K, coordSys)

Build a block-diagonal full-space coordinate transformation for a multifield
system.

MPCs, boundary conditions, and unwrapped loads are interpreted in these local
component bases. `Global(...)` is the explicit exception for right-hand-side
data supplied in global Cartesian components.

Each `NodalCoordinateSystem` identifies its field through `e1.model`. Fields
without a local coordinate system receive an identity block. A coordinate
system may coexist with reduced-order interpolation on another field, but not
yet on the same field.
"""
function _build_multifield_coordinate_transformation(
    K::SystemMatrix,
    coordSys
)
    isempty(coordSys) && return nothing

    K.problems === nothing &&
        error(
            "_build_multifield_coordinate_transformation: " *
            "K must be a multifield SystemMatrix."
        )

    problems = K.problems

    for cs in coordSys
        cs isa NodalCoordinateSystem ||
            error(
                "solveField: coordSys entries must be created from VectorField " *
                "directions, e.g. CoordinateSystem(e1, e2)."
            )

        P = cs.e1.model

        any(Q -> Q === P, problems) ||
            error(
                "solveField: coordinate system refers to a field " *
                "that is not present in the multifield system."
            )

        if cs.e2 !== nothing
            cs.e2.model === P ||
                error(
                    "CoordinateSystem: e1 and e2 must belong to the same Problem."
                )
        end
    end

    Qblocks = SparseMatrixCSC{Float64,Int}[]

    for P in problems
        local_coordSys = [
            cs for cs in coordSys
            if cs.e1.model === P
        ]

        if isempty(local_coordSys)
            push!(
                Qblocks,
                spdiagm(0 => ones(Float64, ndofs(P)))
            )
            continue
        end

        P.pdim == P.dim ||
            error(
                "solveField: coordSys can only be attached to vector fields " *
                "with pdim == dim."
            )

        P.reducedOrder &&
            error(
                "solveField: CoordinateSystem + reducedOrder=true on the same field " *
                "is not yet implemented."
            )

        Qp = _build_coordinate_transformation(
            P,
            local_coordSys
        )

        push!(Qblocks, Qp.T)
    end

    Q = blockdiag(Qblocks...)
    dropzeros!(Q)

    return Q
end

"""
    rhs_in_multifield_solver_coordinates(rhs, Q)

Return the multifield `SystemVector` template and its matrix expressed in the
solver coordinate basis.

An unwrapped `SystemVector` is interpreted in local solver coordinates.
`Global(F)` marks a `SystemVector` as global Cartesian data and applies `Q'`.
A `SolveRHS{SystemVector}` may contain both local and global contributions.
"""
function rhs_in_multifield_solver_coordinates(
    rhs,
    Q
)
    local_part, global_part = _rhs_parts(rhs)
    template = _rhs_template(rhs)

    template isa SystemVector ||
        error(
            "rhs_in_multifield_solver_coordinates: RHS must be a SystemVector."
        )

    if local_part !== nothing && global_part !== nothing
        _check_rhs_compatibility(local_part, global_part)
    end

    B = zeros(Float64, size(template.a))

    if local_part !== nothing
        B .+= local_part.a
    end

    if global_part !== nothing
        if Q === nothing
            B .+= global_part.a
        else
            B .+= Q' * global_part.a
        end
    end

    return template, B
end


# =============================================================================
# Multifield system preparation and reconstruction
# =============================================================================

"""
    prepare_multifield_system(K, rhs, support; mpc=[], coordSys=[])

Prepare a block finite-element system for the common solve pipeline.

Each field receives its own kinematic transformation `Tᵢ` and restriction
`Rᵢ`. Coordinate-system transformations are assembled separately in the
original full space,

    Q = blockdiag(Q₁, Q₂, ...),

and the physical and reduced unknowns are related by

    x_global = Q * T * x_reduced.

Consequently,

    K_reduced = T' * Q' * K * Q * T.

Fields without a coordinate transformation receive an identity `Qᵢ`.
A coordinate system may coexist with reduced-order interpolation on a
different field and may also be combined with field-specific MPCs.

MPC relations are imposed on the local solver components. Therefore, if
master and slave DOFs use different local coordinate bases, the corresponding
global relation contains the relative rotation between those bases. This is
intentional and can represent rotated periodic or sector-symmetry relations.

The caller is responsible for assigning coordinate systems that give the
intended physical meaning to tied DOFs. Attaching `coordSys` and
`reducedOrder=true` to the same field is not yet implemented.
"""
function prepare_multifield_system(
    K0::Union{SystemMatrix,SymmetricSystemMatrix},
    rhs,
    support::Vector{BoundaryCondition};
    mpc::Vector{MPC}=MPC[],
    coordSys=NodalCoordinateSystem[]
)
    K = parent_system(K0)

    K.problems === nothing &&
        error("prepare_multifield_system: K must be a multifield SystemMatrix.")

    template = _rhs_template(rhs)

    template isa SystemVector ||
        error("prepare_multifield_system: RHS must be a SystemVector.")

    template.problems === nothing &&
        error("prepare_multifield_system: RHS must be a multifield SystemVector.")

    K.problems == template.problems ||
        error(
            "prepare_multifield_system: Problem ordering mismatch between K and RHS."
        )

    size(K.A, 1) == size(K.A, 2) ||
        error("prepare_multifield_system: system matrix must be square.")

    size(template.a, 1) == size(K.A, 1) ||
        error("prepare_multifield_system: incompatible matrix and RHS sizes.")

    problems = K.problems
    offsets = K.offsets

    # ----------------------------------------------------------
    # 1) Coordinate basis in the original full multifield space
    # ----------------------------------------------------------

    Q = _build_multifield_coordinate_transformation(
        K,
        coordSys
    )

    rhs_template, F0 =
        rhs_in_multifield_solver_coordinates(
            rhs,
            Q
        )

    nsteps = size(F0, 2)

    A0 =
        Q === nothing ?
        K.A :
        Q' * K.A * Q

    # ----------------------------------------------------------
    # 2) Full-space Dirichlet data
    #
    # With coordSys present the prescribed component values are
    # interpreted directly in the local solver basis.
    # ----------------------------------------------------------

    _, fixed, xD = multifield_bc_data(
        K,
        support;
        nsteps=nsteps
    )

    # ----------------------------------------------------------
    # 3) Transfer prescribed values through MPC relations
    # ----------------------------------------------------------

    fixed_eff, xD_eff = _multifield_mpc_bc_data(
        K,
        mpc,
        fixed,
        xD
    )

    # ----------------------------------------------------------
    # 4) Global field-wise reduction + MPC transformation
    # ----------------------------------------------------------

    T, R, prunable, protected =
        _multifield_mpc_transformation(
            K,
            mpc
        )

    # ----------------------------------------------------------
    # 5) Dirichlet data in the transformed space
    # ----------------------------------------------------------

    free_r, fixed_r, xD_r = reduced_bc_data_matrix(
        R,
        fixed_eff,
        xD_eff
    )

    # ----------------------------------------------------------
    # 6) Galerkin projection
    # ----------------------------------------------------------

    Kr = T' * A0 * T
    Br = T' * F0

    # ----------------------------------------------------------
    # 7) Remove intentional inactive MPC auxiliary DOFs
    # ----------------------------------------------------------

    free_r = _remove_inactive_mpc_dofs(
        Kr,
        Br,
        free_r,
        fixed_r,
        protected,
        prunable
    )

    # ----------------------------------------------------------
    # 8) Reduced free-DOF system
    # ----------------------------------------------------------

    B = copy(Br[free_r, :])

    if !isempty(fixed_r)
        B .-= Kr[free_r, fixed_r] * xD_r[fixed_r, :]
    end

    A = reduced_matrix_wrapper(
        K0,
        Kr[free_r, free_r]
    )

    return (
        A=A,
        B=B,
        T=T,
        R=R,
        Q=Q,
        rhs_template=rhs_template,
        xD_r=xD_r,
        free_r=free_r,
        fixed_r=fixed_r,
        nr=size(T, 2),
        problems=problems,
        offsets=offsets,
        nsteps=nsteps,
        prunable=prunable,
        protected=protected
    )
end

"""
    reconstruct_multifield_solution(Xfree, prepared)

Reconstruct all physical fields from the reduced multifield solution.
The return value is a tuple ordered exactly as `prepared.problems`.
"""
function reconstruct_multifield_solution(
    Xfree,
    prepared
)
    nsteps = size(Xfree, 2)
    xr = zeros(eltype(Xfree), prepared.nr, nsteps)

    xr[prepared.free_r, :] .= Xfree

    if !isempty(prepared.fixed_r)
        xr[prepared.fixed_r, :] .=
            prepared.xD_r[prepared.fixed_r, :]
    end

    x_local = prepared.T * xr
    x = prepared.Q === nothing ? x_local : prepared.Q * x_local

    results = Vector{Any}(undef, length(prepared.problems))
    t = collect(0.0:(nsteps - 1))

    for (i, P) in enumerate(prepared.problems)
        offset = prepared.offsets[i]
        nloc = ndofs(P)
        xloc = Matrix(@view x[offset + 1:offset + nloc, :])

        if P.pdim == 1
            results[i] = ScalarField(
                Matrix{Float64}[],
                xloc,
                t,
                Int[],
                nsteps,
                :scalar,
                P
            )
        elseif P.pdim == 2 || P.pdim == 3
            type = P.pdim == 2 ? :v2D : :v3D
            results[i] = VectorField(
                Matrix{Float64}[],
                xloc,
                t,
                Int[],
                nsteps,
                type,
                P
            )
        elseif P.pdim == 9
            results[i] = TensorField(
                Matrix{Float64}[],
                xloc,
                t,
                Int[],
                nsteps,
                :tensor,
                P
            )
        else
            error(
                "reconstruct_multifield_solution: unsupported pdim $(P.pdim)."
            )
        end
    end

    return tuple(results...)
end

function _resolve_solver(
    K,
    solver,
    iterative,
    ordering
)
    if iterative !== nothing
        iterative isa Bool ||
            error("solveField: `iterative` must be `true`, `false`, or omitted.")

        solver === :auto ||
            error("solveField: use either `solver` or the legacy `iterative` keyword, not both.")

        iterative && return :cg
    end

    if ordering !== nothing
        ordering isa Bool ||
            error("solveField: `ordering` must be `true`, `false`, or omitted.")

        solver === :auto ||
            error("solveField: use either `solver` or the legacy `ordering` keyword, not both.")

        ordering || return :lu_no_ordering
    end

    return solver === :auto ? default_solver(K) : solver
end

function _solve_prepared_system(
    K,
    prepared;
    solver=:auto,
    iterative=nothing,
    ordering=nothing,
    preconditioner=nothing,
    reltol=sqrt(eps()),
    maxiter=size(K isa SymmetricSystemMatrix ? K.parent.A : K.A, 1)
)
    solver0 = _resolve_solver(
        K,
        solver,
        iterative,
        ordering
    )

    return solve_linear_system(
        prepared.A,
        prepared.B;
        solver=solver0,
        preconditioner=preconditioner,
        reltol=reltol,
        maxiter=maxiter
    )
end


# =============================================================================
# Public solver API
# =============================================================================

"""
    solveField(K, rhs;
               support=BoundaryCondition[],
               mpc=MPC[],
               coordSys=NodalCoordinateSystem[],
               solver=:auto,
               preconditioner=nothing,
               reltol=sqrt(eps()),
               maxiter=...)

Solve a linear single-field finite-element system.

The solution pipeline is

```text
prepare -> solve_linear_system -> reconstruct
```

and uses the transformation

```math
u_\\mathrm{global} = Q T u_\\mathrm{reduced},
```

where `Q` is the nodal coordinate-basis transformation and `T` contains
kinematic transformations such as reduced-order interpolation and MPCs.

# Right-hand side

`rhs` may be a `ScalarField`, `VectorField`, `TensorField`, `Global(...)`
wrapper, or a symbolic sum of local and global load terms. Multiple
right-hand sides are supported columnwise.

With `coordSys` present, ordinary loads and prescribed component values are
interpreted in the local solver basis. `Global(...)` marks right-hand-side
data that are supplied in global Cartesian components.

# Coordinate systems and MPCs

MPCs identify local solver components. Therefore master and slave nodes may
deliberately use different local bases. In global coordinates the relation is

```math
u_{g,s} = Q_s Q_m^T u_{g,m},
```

which can represent rotated periodic or sector-symmetry constraints.

`CoordinateSystem + reducedOrder=true` on the same field is not yet
implemented.

# Linear solvers

`solver` may be `:auto`, `:backslash`, `:lu`, `:lu_no_ordering`,
`:cholesky`, `:qr`, `:cg`, or `:gmres`.

For an ordinary `SystemMatrix`, `:auto` uses Julia's backslash solver. For a
`SymmetricSystemMatrix`, `:auto` first attempts Cholesky and falls back to LU
when the matrix is symmetric but not positive definite.

The legacy keywords `iterative` and `ordering` are still accepted for backward
compatibility. `iterative=true` selects conjugate gradient; `ordering=false`
selects LU without column reordering. New code should prefer `solver`.
"""
function solveField(
    K::Union{SystemMatrix,SymmetricSystemMatrix},
    rhs::Union{_SolveFieldType,GlobalField,SolveRHS};
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[],
    coordSys=NodalCoordinateSystem[],
    solver=:auto,
    iterative=nothing,
    ordering=nothing,
    preconditioner=nothing,
    reltol=sqrt(eps()),
    maxiter=size(K isa SymmetricSystemMatrix ? K.parent.A : K.A, 1)
)
    prepared = prepare_singlefield_system(
        K,
        rhs,
        support;
        mpc=mpc,
        coordSys=coordSys
    )

    Xfree = _solve_prepared_system(
        K,
        prepared;
        solver=solver,
        iterative=iterative,
        ordering=ordering,
        preconditioner=preconditioner,
        reltol=reltol,
        maxiter=maxiter
    )

    return reconstruct_singlefield_solution(
        Xfree,
        prepared
    )
end


"""
    solveField(K, F::SystemVector; kwargs...)

Solve a linear multifield finite-element system.

Full-order and reduced-order fields may be mixed in the same block system.
Each field may have its own MPCs and local coordinate systems. A coordinate
system may coexist with `reducedOrder=true` on another field, but
`CoordinateSystem + reducedOrder=true` on the same field is not yet
implemented.

For a multifield system,

```math
x_\\mathrm{global} = Q T x_\\mathrm{reduced},
```

with block-diagonal field-wise coordinate transformation `Q`. MPC relations
are imposed on local solver components, so deliberately different master and
slave bases naturally produce rotated periodic relations.

`F` may be a `SystemVector`, `Global(F)`, or a symbolic sum containing local
and global `SystemVector` contributions. Multiple right-hand sides are
supported.

The solver keywords and backward-compatibility behavior are the same as for
the single-field method.
"""
function solveField(
    K::Union{SystemMatrix,SymmetricSystemMatrix},
    rhs::Union{
        SystemVector,
        GlobalField{<:SystemVector},
        SolveRHS{SystemVector}
    };
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[],
    coordSys=NodalCoordinateSystem[],
    solver=:auto,
    iterative=nothing,
    ordering=nothing,
    preconditioner=nothing,
    reltol=sqrt(eps()),
    maxiter=size(K isa SymmetricSystemMatrix ? K.parent.A : K.A, 1)
)
    prepared = prepare_multifield_system(
        K,
        rhs,
        support;
        mpc=mpc,
        coordSys=coordSys
    )

    Xfree = _solve_prepared_system(
        K,
        prepared;
        solver=solver,
        iterative=iterative,
        ordering=ordering,
        preconditioner=preconditioner,
        reltol=reltol,
        maxiter=maxiter
    )

    return reconstruct_multifield_solution(
        Xfree,
        prepared
    )
end
