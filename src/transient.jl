# Transient solvers: prepare a common constrained space, integrate, reconstruct.
export CDM, HHT, FDM, CDMaccuracyAnalysis, HHTaccuracyAnalysis, FDMaccuracyAnalysis

const _TransientMatrix = Union{SystemMatrix,SymmetricSystemMatrix}
const _TransientData = Union{_SolveFieldType,SystemVector,GlobalField,SolveRHS,
                             Tuple,AbstractVector{<:_SolveFieldType}}

function check_multifield_system_compatibility(K::SystemMatrix, C::SystemMatrix)
    K.problems !== nothing && C.problems !== nothing ||
        error("Transient/eigen system: both matrices must be block SystemMatrices.")
    K.offsets !== nothing && C.offsets !== nothing ||
        error("Transient/eigen system: missing block offsets.")
    length(K.problems) == length(C.problems) &&
        all(a === b for (a, b) in zip(K.problems, C.problems)) &&
        K.offsets == C.offsets || error("Transient/eigen system: field ordering differs.")
    size(K.A) == size(C.A) && size(K.A, 1) == size(K.A, 2) ||
        error("Transient/eigen system: incompatible matrix dimensions.")
    return nothing
end

function _transient_check_matrix(K, A)
    size(K.A) == size(A.A) && size(A.A, 1) == size(A.A, 2) ||
        throw(DimensionMismatch("Transient: matrices must have matching square dimensions."))
    if K.problems === nothing
        A.problems === nothing && K.model !== nothing && K.model === A.model &&
            K.test_model === K.model && A.test_model === A.model ||
            error("Transient: matrices must refer to the same test/trial Problem.")
        size(K.A, 1) == ndofs(K.model) ||
            throw(DimensionMismatch("Transient: matrix dimension does not match the Problem."))
    else
        check_multifield_system_compatibility(K, A)
    end
    all(isfinite, nonzeros(A.A)) || error("Transient: non-finite matrix entries.")
end

function _transient_time(n, dt, t0)
    n >= 2 || throw(ArgumentError("Transient: n must be at least 2 (including the initial sample)."))
    isfinite(dt) && dt > 0 && isfinite(t0) ||
        throw(ArgumentError("Transient: require finite Δt > 0 and finite t0."))
    t = Float64[t0 + (j - 1) * dt for j in 1:n]
    all(isfinite, t) || throw(ArgumentError("Transient: non-finite time grid."))
    return t
end

function _transient_field(P, a, t)
    kind = P.pdim == 1 ? :scalar : P.pdim == 2 ? :v2D :
           P.pdim == 3 ? :v3D : P.pdim == 9 ? :tensor :
           error("Transient: unsupported component count $(P.pdim).")
    FT = P.pdim == 1 ? ScalarField : P.pdim in (2, 3) ? VectorField : TensorField
    return FT(Matrix{Float64}[], Matrix{Float64}(a), collect(t), Int[], length(t), kind, P)
end

function split_multifield_solution(X::AbstractMatrix, problems::Vector{Problem},
                                  offsets::Vector{Int}, t::AbstractVector)
    results = map(eachindex(problems)) do i
        P = problems[i]
        _transient_field(P, @view(X[offsets[i]+1:offsets[i]+ndofs(P), :]), t)
    end
    return length(results) == 1 ? only(results) : tuple(results...)
end

# Input values are sampled columnwise. A single column is constant in time.
# Initial values, like ordinary loads and BC components, use the solver basis;
# Global(initial) explicitly requests Cartesian input.
function _transient_values(data, K, Q, n)
    if data isa GlobalField
        A = _transient_values(data.field, K, nothing, n)
        return Q === nothing ? A : Q' * A
    elseif data isa SolveRHS
        a = data.local_part === nothing ? nothing :
            _transient_values(data.local_part, K, nothing, n)
        b = data.global_part === nothing ? nothing :
            _transient_values(Global(data.global_part), K, Q, n)
        return a === nothing ? b : b === nothing ? a : a .+ b
    elseif data isa Tuple || data isa AbstractVector{<:_SolveFieldType}
        return _transient_values(SystemVector([elementsToNodes(f) for f in data]), K, Q, n)
    end
    data isa _SolveFieldType && (data = elementsToNodes(data))
    if K.problems === nothing
        data isa Union{_SolveFieldType,SystemVector} && data.model === K.model ||
            error("Transient: input belongs to a different Problem.")
        data isa SystemVector && data.problems !== nothing &&
            error("Transient: expected single-field data.")
    else
        data isa SystemVector || error("Transient: block data must be a SystemVector or field tuple.")
        data.problems !== nothing && length(data.problems) == length(K.problems) &&
            all(a === b for (a, b) in zip(data.problems, K.problems)) &&
            data.offsets == K.offsets || error("Transient: input block ordering differs from the matrices.")
    end
    A = data.a
    size(A, 1) == size(K.A, 1) || throw(DimensionMismatch("Transient: incompatible input DOF count."))
    size(A, 2) in (1, n) || throw(DimensionMismatch("Transient: inputs must have 1 or n columns."))
    all(isfinite, A) || error("Transient: non-finite input values.")
    return A
end

# Sample only constrained rows. Physical-group lookup and component resolution
# happen once per BC, never inside the time loop. Prefer a spatial signature
# when a varargs/optional-argument function supports both signatures.
function _transient_support(P, support, t)
    entries = []
    fixed = Int[]
    gmsh.model.setCurrent(P.name)
    components = _bc_component_map(P)
    prefix = String(P.field)
    for bc in support
        bc.problem !== nothing && bc.problem !== P && continue
        tag = getTagForPhysicalName(bc.phName)
        nodes, coordinates = gmsh.model.mesh.getNodesForPhysicalGroup(-1, tag)
        for (key, value) in bc.values
            name = String(key)
            startswith(name, prefix) || error("Transient: unknown BC key '$key'.")
            suffix = String(chopprefix(name, prefix))
            haskey(components, suffix) || error("Transient: invalid BC component '$suffix'.")
            dofs = (Int.(nodes) .- 1) .* P.pdim .+ components[suffix]
            timed = false
            if value isa Function
                spatial = applicable(value, 0.0, 0.0, 0.0)
                timed = !spatial && applicable(value, 0.0, 0.0, 0.0, first(t))
                spatial || timed || throw(ArgumentError("Transient: BC functions must accept (x,y,z) or (x,y,z,t)."))
            elseif value isa ScalarField
                value = elementsToNodes(value)
                size(value.a, 2) in (1, length(t)) ||
                    throw(DimensionMismatch("Transient: BC ScalarField must have 1 or n columns."))
            elseif !(value isa Number)
                throw(ArgumentError("Transient: BC components must be numbers, ScalarFields or functions."))
            end
            push!(entries, (; dofs, nodes, coordinates, value, timed))
            append!(fixed, dofs)
        end
    end
    sort!(unique!(fixed))
    row = Dict(dof => i for (i, dof) in enumerate(fixed))
    D = zeros(length(fixed), length(t))
    for entry in entries
        (; dofs, nodes, coordinates, value, timed) = entry
        for (i, dof) in enumerate(dofs)
            r = row[dof]
            if value isa Function
                x, y, z = coordinates[3i-2], coordinates[3i-1], coordinates[3i]
                if timed
                    for (j, tj) in enumerate(t)
                        D[r, j] = value(x, y, z, tj)
                    end
                else
                    D[r, :] .= value(x, y, z)
                end
            elseif value isa ScalarField
                D[r, :] .= @view value.a[nodes[i], :]
            else
                D[r, :] .= value
            end
        end
    end
    return fixed, D
end

# Transfer compact prescribed rows through representative maps once, then
# restrict only the constrained reduced rows. No ndof-by-n temporary is needed.
function _transient_boundary_data(problems, offsets, mpc, R, support, t)
    histories = Matrix{Float64}[]
    masters = Int[]
    for (i, P) in enumerate(problems)
        local_fixed, values = _transient_support(P, support, t)
        push!(histories, values)
        if isempty(mpc)
            append!(masters, offsets[i] .+ local_fixed)
        else
            rep = mpcRepresentativeMap(P, mpc)
            _canonicalize_mpc_rep!(rep)
            append!(masters, offsets[i] .+ rep[local_fixed])
        end
    end
    D = vcat(histories...)
    all(isfinite, D) || error("Transient: non-finite prescribed values.")
    fixed_eff = sort!(unique(masters))
    effective = zeros(length(fixed_eff), length(t))
    rows = Dict(dof => i for (i, dof) in enumerate(fixed_eff))
    assigned = falses(length(fixed_eff))
    for (i, master) in enumerate(masters)
        r = rows[master]
        if assigned[r]
            all(isapprox.(effective[r, :], D[i, :])) ||
                error("MPC: incompatible prescribed values on tied DOFs.")
        else
            effective[r, :] .= @view D[i, :]
            assigned[r] = true
        end
    end
    free, fixed_r, _ = reduced_bc_data(R, fixed_eff, spzeros(size(R, 2)))
    Dr = R[fixed_r, fixed_eff] * effective
    return free, fixed_r, Dr
end

"""
    prepare_transient_system(K, M, C, rhs, U0, V0, n, Δt; ...)

Common coordinate/approximation/MPC preparation for all transient methods.
`M` and `V0` may be `nothing` for a first-order system. Matrices and input
fields are never modified. Dirichlet values are transferred through the same
MPC and approximation helpers as `solveField` and `solveEigen`.
`U0` and `V0` contain one initial column. All prescribed primary-variable
histories come from `support`; the integrator determines their derivatives.
"""
function prepare_transient_system(K0, M0, C0, rhs, U0, V0, n, dt;
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[], coordSys=CoordinateSystem[], t0=0.0, _modifier_data=false)
    t = _transient_time(n, dt, t0)
    K = parent_system(K0)
    M = M0 === nothing ? nothing : parent_system(M0)
    C = C0 === nothing ? nothing : parent_system(C0)
    for A in (K, M, C)
        A === nothing || _transient_check_matrix(K, A)
    end
    problems = K.problems === nothing ? [K.model] : K.problems
    offsets = K.problems === nothing ? [0] : K.offsets
    for bc in support
        if K.problems !== nothing && bc.problem === nothing
            error("Transient: multifield supports must specify field=Problem.")
        end
        bc.problem === nothing || any(P -> P === bc.problem, problems) ||
            error("Transient: support refers to a field outside the system.")
    end
    # Follow the existing solveField/solveEigen limitation explicitly.
    for P in problems
        if _has_approximation(P) && any(cs -> cs.e1.model === P, coordSys)
            error("Transient: coordSys + reduced/spectral basis on the same field is not yet implemented.")
        end
    end
    if K.problems === nothing
        basis = _build_coordinate_transformation(K.model, coordSys)
        Q = basis === nothing ? nothing : basis.T
        T, R = field_transformation(K.model, mpc)
    else
        Q = _build_multifield_coordinate_transformation(K, coordSys)
        T, R, prunable, protected = _multifield_mpc_transformation(K, mpc)
    end
    L = Q === nothing ? T : Q*T
    project(A) = _project_system_matrix(A, L)
    Kr = project(system_matrix_data(K0))
    Mr = M === nothing ? spzeros(size(Kr)...) : project(system_matrix_data(M0))
    Cr = C === nothing ? spzeros(size(Kr)...) : project(system_matrix_data(C0))
    F0 = _transient_values(rhs, K, Q, n)
    F = T' * F0
    initial = _transient_values(U0, K, Q, 1)
    u0 = Vector(R * @view(initial[:, 1]))
    v0 = V0 === nothing ? zeros(size(T, 2)) :
        Vector(R * @view(_transient_values(V0, K, Q, 1)[:, 1]))
    free, fixed_r, Dr = _transient_boundary_data(problems, offsets, mpc, R, support, t)
    if K.problems !== nothing && !isempty(mpc)
        free = _remove_inactive_mpc_dofs((Kr, Mr, Cr), F, free, fixed_r, protected, prunable)
    end
    return (; K=Kr, M=Mr, C=Cr, F=Matrix(F[free, :]), u0=u0[free], v0=v0[free],
            D=Dr, free, fixed=fixed_r, T, Q, problems, offsets, t,
            layout=K, L=_modifier_data ? L : nothing,
            Fsource=_modifier_data ? F0 : nothing,
            inactive=_modifier_data ? setdiff(1:size(T,2),union(free,fixed_r)) : Int[])
end

"""
    reconstruct_transient_solution(U, D, prepared)

Insert free values `U` and constrained values `D`, prolong through `T`, rotate
through `Q`, and return ordinary nodal Fields with physical time labels.
Multifield output follows the matrix block order.
"""
function reconstruct_transient_solution(U, D, p)
    X = zeros(size(p.T, 2), size(U, 2))
    X[p.free, :] .= U
    X[p.fixed, :] .= D
    Xfull = p.T * X
    p.Q === nothing || (Xfull = p.Q * Xfull)
    return split_multifield_solution(Xfull, p.problems, p.offsets, p.t)
end

# Second-order endpoint differences; two samples define a linear history.
function _transient_derivatives(D, dt)
    V, A = zeros(size(D)), zeros(size(D))
    n = size(D, 2)
    if n == 2
        V[:, 1] .= V[:, 2] .= (D[:, 2] - D[:, 1]) / dt
    else
        V[:, 1] .= (-3D[:, 1] + 4D[:, 2] - D[:, 3]) / (2dt)
        V[:, n] .= (3D[:, n] - 4D[:, n-1] + D[:, n-2]) / (2dt)
        for j in 2:n-1
            V[:, j] .= (D[:, j+1] - D[:, j-1]) / (2dt)
        end
        for j in 2:n-1
            A[:, j] .= (D[:, j+1] - 2D[:, j] + D[:, j-1]) / dt^2
        end
        if n >= 4
            A[:, 1] .= (2D[:, 1] - 5D[:, 2] + 4D[:, 3] - D[:, 4]) / dt^2
            A[:, n] .= (2D[:, n] - 5D[:, n-1] + 4D[:, n-2] - D[:, n-3]) / dt^2
        else
            A[:, 1] .= A[:, n] .= A[:, 2]
        end
    end
    return V, A
end

function _transient_factor(A; mass=false)
    if isdiag(A)
        # Sparse diagonal storage would make Diagonal solves return sparse
        # state vectors, which cannot be passed to the cached UMFPACK solve.
        d = Vector(diag(A))
        all(x -> isfinite(x) && (mass ? x > 0 : !iszero(x)), d) ||
            error("Transient: $(mass ? "mass/capacity" : "effective") matrix has invalid diagonal entries.")
        return Diagonal(d)
    end
    return lu(A)
end

# Preserve symmetry only when every contributing input matrix explicitly
# declares it. In particular, symmetric K/M do not imply symmetric damping.
function _prepare_transient_linear_solver(A, matrices; solver=:auto,
    preconditioner=:auto, solveroptions=(;), preconditioneroptions=(;))
    symmetric = all(K -> K === nothing || K isa SymmetricSystemMatrix, matrices)
    matrix = symmetric ? Symmetric(A) : A
    backend = solver === :auto ? (symmetric ? :auto_symmetric : :backslash) : solver
    return prepare_linear_solver(matrix; solver=backend, preconditioner,
                                 solveroptions, preconditioneroptions)
end

_transient_column(F, j) = @view F[:, size(F, 2) == 1 ? 1 : j]

"""
    HHT(K, M, C, f, U0, V0, n, Δt; kwargs...)
    HHT(K, M, f, U0, V0, n, Δt; damping=nothing, kwargs...)

Integrate a second-order system with the **HHT-α / Newmark method**:

```math
M\\ddot{u} + C\\dot{u} + Ku = f.
```

# Arguments

| Argument | Meaning |
|:--|:--|
| `K`, `M`, `C` | Field, second-derivative and first-derivative operators; single- or multifield `SystemMatrix`s. |
| `f` | Constant or sampled RHS Field, `SystemVector`, field tuple or local/global sum. |
| `U0`, `V0` | Initial primary field and first derivative, one column each, on free DOFs. |
| `n` | Number of stored samples, **including the initial sample**; `n ≥ 2`. |
| `Δt` | Constant finite positive time increment. |

Omit positional `C` for zero physical damping. Only in that calling form,
`damping=C` supplies the same operator by keyword (default `nothing`).
The positional-support alias `HHT(K,M,[C,]f,bc,U0,V0,n,Δt)` is retained.

# HHT parameters

| Keyword | Default | Meaning |
|:--|:--|:--|
| `α` | `0.0` | Algorithmic damping, in `[-1/3,0]`; distinct from physical `C`. |
| `δ` | `nothing` | Optional override of the relation `γ=0.5+δ`. |
| `γ` | `0.5-α` | Velocity parameter when `δ` is omitted; explicit values override it. |
| `β` | `(0.5+γ)^2/4` | Displacement parameter; positive, and explicitly overridable. |

The standard choice is `γ=0.5-α`, `β=(1-α)^2/4`; `α=0` gives average-acceleration
Newmark. Explicit `δ=0` recovers the historical independent parameter family.
Custom `γ ≥ 0`, `β > 0` are accepted but need their own stability assessment.

# Data, constraints and coordinates

| Keyword | Default | Meaning |
|:--|:--|:--|
| `support` | `BoundaryCondition[]` | Prescribed components and their complete time histories. |
| `mpc` | `MPC[]` | Multi-point constraints, including multifield remote couplings. |
| `coordSys` | `CoordinateSystem[]` | Local bases for loads, initial values, supports and MPC components. |
| `t0` | `0.0` | Time of the initial sample. |

Loads have one constant column or exactly `n` sampled columns. Initial fields
have **one column**. Input time labels are not resampled. `support` is the only
source of prescribed histories: use numbers, spatial `(x,y,z)` functions,
time-dependent `(x,y,z,t)` functions or sampled `ScalarField` components.
Supported values override initial data on constrained DOFs; no external
`applyBoundaryConditions!` call is needed.

Ordinary loads and initial fields use local solver components. Use `Global(...)`
for Cartesian input or combine local/global loads. Output is always global.
For multiple fields, use a `SystemVector` or tuple in matrix block order and
identify `field=Problem` in every BC/MPC. Reduced-order and SEM spaces follow
the shared solver transformations.

Prescribed derivatives come from the support history, using centered interior
and one-sided endpoint differences. With two samples acceleration is zero;
with three it is the same second difference at all samples; endpoint acceleration
is second-order for at least four samples. Only free unknowns use Newmark
recurrence. Prescribed inertia, damping and stiffness enter the effective load
with the same HHT endpoint weighting as applied loads.

# Linear solver options

| Keyword | Default | Meaning |
|:--|:--|:--|
| `solver` | `:auto` | Select a direct or iterative solve; choices below. |
| `preconditioner` | `:auto` | `:ich` for CG, `:ilu` for GMRES; also `:none` or a supplied object. |
| `solveroptions` | `(;)` | Native iterative options, e.g. `(reltol=1e-9, maxiter=1000)`. |
| `preconditioneroptions` | `(;)` | Options for a symbolic preconditioner, e.g. `(memory=2,)` for `:ich`. |

- Direct choices: `:auto`, `:backslash`, `:lu`, `:lu_no_ordering`,
  `:cholesky`, `:ldlt`, `:qr`.
- Iterative choices: `:cg` (positive-definite symmetric systems) and `:gmres`.
  Load `IterativeSolvers`; `:ich` additionally needs `Preconditioners`,
  and `:ilu` needs `IncompleteLU`.
- Wrap every contributing operator in `Symmetric(...)` to declare the effective
  matrix symmetric. Then `:auto` tries Cholesky with LU fallback; otherwise it
  uses the general sparse factorization. Symmetry alone does not imply positive definiteness.
- Nonempty option tuples and actual preconditioners are for iterative solves.
  `log`, `Pl`, `Pr` and `initially_zero` are managed internally. Constructor
  options cannot accompany `:none` or a preconstructed preconditioner.

For unchanged coefficients, the effective factorization/preconditioner is reused.
See [`solveField`](@ref) for the shared backend contract.

The initial acceleration uses a separate mass solve with the same solver
settings; diagonal mass uses direct division.

# Operator updates and nonlinear iteration

The generic names are `A2=M`, `A1=C`, `A0=K`.

| Keyword | Default | Updated quantity |
|:--|:--|:--|
| `modifyA2` | `nothing` | Second-derivative operator, e.g. mass. |
| `modifyA1` | `nothing` | First-derivative operator, e.g. damping. |
| `modifyA0` | `nothing` | Primary-field operator, e.g. stiffness or tangent. |
| `modifyRHS` | `nothing` | RHS at the evaluation time. |
| `nonlinearoptions` | `(;)` | Nonlinear iteration controls listed below. |

HHT reevaluates callbacks within each step until discrete equilibrium converges.

The callback signature is `callback(operator, state)`. Matrices are full/global
`SystemMatrix`s (optionally `Symmetric`); RHS results are one-column global Fields,
field tuples or `SystemVector`s. `state` includes trial fields and derivatives,
Problems/material data, time, step, iteration, previous accepted snapshot and
current operators. Callback order is **A2, A0, A1, RHS**.
See the [complete callback guide](@ref transient-callbacks).

The operator argument is a private copy. Return it to keep in-place edits;
return `nothing` to keep the accepted matrix or the original prescribed RHS
sample. State data are read-only. Mesh, field ordering, supports, MPCs and bases
stay fixed; inactive MPC DOFs cannot be activated.

| `nonlinearoptions` member | Default | Meaning |
|:--|:--|:--|
| `maxiter` | `30` | Maximum nonlinear updates; positive integer. |
| `reltol` | `1e-9` | Relative residual tolerance; nonnegative. |
| `abstol` | `1e-11` | Absolute residual tolerance, in equilibrium-equation units; nonnegative. |
| `relaxation` | `1.0` | Fraction of each correction, in `(0,1]`. |

For Newton iteration, return the tangent `J(u)` as A0 and
`external - internal(u) + J(u)*u` as RHS. State-dependent A2/A1 also require
consistent linearization for Newton convergence. Picard updates are permitted.
Nonconvergence throws a step/time/residual error. Without callbacks, leave
`nonlinearoptions=(;)` to use the linear path.

# Returns

Return `(u,v)` as global nodal Fields, with `n` columns at
`t[j] = t0 + (j-1)*Δt`. For multiple fields, return tuples in block order:
`(u1,u2), (v1,v2) = HHT(...)`.

# Choosing the time step

For linear positive-definite mass, nonnegative stiffness and physical damping,
the standard HHT parameter family is unconditionally stable. This is a
**stability property, not a time-resolution guarantee**. Start with about
20–50 samples per shortest period of interest, resolve rapid loads/support
motion, then halve `Δt` and compare the quantities of interest.
Use [`HHTaccuracyAnalysis`](@ref) to inspect damping and phase effects;
see the [time-step guide](@ref transient-timestep) for a practical workflow.

# Example

```julia
bc = [BoundaryCondition("base"; field=U,
                        ux=(x,y,z,t)->0.001*sin(20t), uy=0.0)]
u, v = HHT(K, M, C, f, u0, v0, n, Δt; support=bc, α=-0.1)
```

# Limitations

Mass must be nonsingular on free reduced coordinates; arbitrary mass-singular
DAEs are not supported. Combining coordinates with reduced/SEM interpolation
on the same field, or SEM with MPC on the same field, follows the shared
backend's current restrictions. See also [`CDM`](@ref), [`FDM`](@ref).
"""
function HHT(K::_TransientMatrix, M::_TransientMatrix, C::_TransientMatrix,
             f::_TransientData, U0, V0, n::Int, dt::Real;
             α=0.0, δ=nothing, γ=0.5+(δ === nothing ? -α : δ), β=0.25*(0.5+γ)^2, kwargs...)
    return _HHT(K, M, C, f, U0, V0, n, dt; α, γ, β, kwargs...)
end
function HHT(K::_TransientMatrix, M::_TransientMatrix, f::_TransientData, U0, V0, n::Int, dt::Real;
             damping=nothing, α=0.0, δ=nothing, γ=0.5+(δ === nothing ? -α : δ), β=0.25*(0.5+γ)^2, kwargs...)
    return _HHT(K, M, damping, f, U0, V0, n, dt; α, γ, β, kwargs...)
end
function _HHT(K, M, C, f, U0, V0, n, dt; α, γ, β,
    solver=:auto, preconditioner=:auto, solveroptions=(;), preconditioneroptions=(;),
    modifyA2=nothing, modifyA1=nothing, modifyA0=nothing, modifyRHS=nothing,
    nonlinearoptions=(;), kwargs...)
    all(isfinite, (α, γ, β)) && -1/3 <= α <= 0 && γ >= 0 && β > 0 ||
        throw(ArgumentError("HHT: require -1/3 ≤ α ≤ 0, γ ≥ 0 and β > 0, all finite."))
    if _has_transient_modifiers(modifyA2,modifyA1,modifyA0,modifyRHS)
        return _modified_hht(K,M,C,f,U0,V0,n,dt; α,γ,β,solver,preconditioner,
            solveroptions,preconditioneroptions,modifyA2,modifyA1,modifyA0,modifyRHS,
            nonlinearoptions,kwargs...)
    end
    nonlinearoptions isa NamedTuple && isempty(nonlinearoptions) ||
        throw(ArgumentError("nonlinearoptions require a modify callback and must be a NamedTuple."))
    p = prepare_transient_system(K, M, C, f, U0, V0, n, dt; kwargs...)
    free, fixed = p.free, p.fixed
    Kff, Mff, Cff = p.K[free, free], p.M[free, free], p.C[free, free]
    Kfc, Mfc, Cfc = p.K[free, fixed], p.M[free, fixed], p.C[free, fixed]
    Vd, Ad = _transient_derivatives(p.D, dt)
    U, V = zeros(length(free), n), zeros(length(free), n)
    U[:, 1], V[:, 1] = p.u0, p.v0
    if !isempty(free)
        effective = _prepare_transient_linear_solver(
            Mff + (1+α)*γ*dt*Cff + (1+α)*β*dt^2*Kff, (K, M, C);
            solver, preconditioner, solveroptions, preconditioneroptions)
        boundary_force(j) = Kfc * view(p.D, :, j) + Cfc * view(Vd, :, j) + Mfc * view(Ad, :, j)
        oldforce = Kff*p.u0 + Cff*p.v0 + boundary_force(1) - _transient_column(p.F, 1)
        rhs0 = -oldforce
        if isdiag(Mff)
            a = _transient_factor(Mff; mass=true) \ rhs0
        else
            mass = _prepare_transient_linear_solver(Mff, (M,);
                solver, preconditioner, solveroptions, preconditioneroptions)
            a = solve_prepared_linear_system(mass, rhs0)
        end
        for j in 2:n
            u, v = @view(U[:, j-1]), @view(V[:, j-1])
            up = u + dt*v + dt^2*(0.5-β)*a
            vp = v + dt*(1-γ)*a
            force = Kff*up + Cff*vp + boundary_force(j) - _transient_column(p.F, j)
            rhs = -(1+α)*force
            iszero(α) || (rhs .+= α .* oldforce)
            a = solve_prepared_linear_system(effective, rhs)
            U[:, j] .= up + β*dt^2*a
            V[:, j] .= vp + γ*dt*a
            if !iszero(α)
                oldforce = force + β*dt^2*(Kff*a) + γ*dt*(Cff*a)
            end
        end
    end
    return reconstruct_transient_solution(U, p.D, p), reconstruct_transient_solution(V, Vd, p)
end

function _transient_cdm_options(dampingScheme, velocityScheme)
    dampingScheme in (:lagged, :central) ||
        throw(ArgumentError("CDM: dampingScheme must be :lagged or :central."))
    velocityScheme in (:backward, :central) ||
        throw(ArgumentError("CDM: velocityScheme must be :backward or :central."))
end

"""
    CDM(K, M, C, f, U0, V0, n, Δt; kwargs...)
    CDM(K, M, f, U0, V0, n, Δt; damping=nothing, kwargs...)

Integrate a second-order system by **central differences**:

```math
M\\ddot{u} + C\\dot{u} + Ku = f.
```

# Arguments

| Argument | Meaning |
|:--|:--|
| `K`, `M`, `C` | Primary-field, second-derivative and first-derivative operators. |
| `f` | Constant or sampled RHS; Field, `SystemVector`, tuple or local/global sum. |
| `U0`, `V0` | Initial primary field and first derivative; one column each. |
| `n`, `Δt` | Stored sample count (`n ≥ 2`, including the initial sample) and positive constant increment. |

Single-field and multifield `SystemMatrix`s are accepted. Omit `C` for zero
damping; only that calling form accepts `damping=C` (default `nothing`).
The positional-support alias `CDM(K,M,[C,]f,bc,U0,V0,n,Δt)` is retained.

# CDM scheme options

| Keyword | Default | Meaning |
|:--|:--|:--|
| `dampingScheme` | `:lagged` | Backward damping velocity; `:central` uses centered damping. |
| `velocityScheme` | `:backward` | Historical half-step free-velocity estimates; `:central` returns sample-time estimates. |

- `:lagged` damping uses `C*(u[j]-u[j-1])/Δt` and the free mass matrix.
- `:central` damping uses `C*(u[j+1]-u[j-1])/(2Δt)` and the effective matrix
  `M + Δt*C/2`; it is second-order in the damping term.
- Central output velocities use centered interior and one-sided final
  differences; the initial free value remains `V0`.

Diagonal systems use componentwise division; otherwise a factorization is
cached. A consistent or projected mass is supported and is **not automatically
lumped**. CDM has no `solver` or preconditioner keywords.

# Data, constraints and coordinates

| Keyword | Default | Meaning |
|:--|:--|:--|
| `support` | `BoundaryCondition[]` | Prescribed components and their complete time histories. |
| `mpc` | `MPC[]` | Multi-point constraints, including multifield remote couplings. |
| `coordSys` | `CoordinateSystem[]` | Local bases for loads, initial values, supports and MPC components. |
| `t0` | `0.0` | Time of the initial sample. |

Loads have one constant column or exactly `n` sampled columns. Initial fields
have **one column**. Input time labels are not resampled. `support` is the only
source of prescribed histories: use numbers, spatial `(x,y,z)` functions,
time-dependent `(x,y,z,t)` functions or sampled `ScalarField` components.
Supported values override initial data on constrained DOFs; no external
`applyBoundaryConditions!` call is needed.

Ordinary loads and initial fields use local solver components. Use `Global(...)`
for Cartesian input or combine local/global loads. Output is always global.
For multiple fields, use a `SystemVector` or tuple in matrix block order and
identify `field=Problem` in every BC/MPC. Reduced-order and SEM spaces follow
the shared solver transformations.

Constrained derivatives always follow the prescribed history by finite
differences, independently of `velocityScheme`. Two samples define linear
motion with zero acceleration; four samples permit second-order acceleration
at both endpoints. The unsupported coordinate/reduced/SEM combinations are
the same as for [`HHT`](@ref); mass must be nonsingular on free coordinates.

# Operator updates

| Keyword | Default | Updated quantity |
|:--|:--|:--|
| `modifyA2` | `nothing` | A2=M: second-derivative operator. |
| `modifyA1` | `nothing` | A1=C: first-derivative operator. |
| `modifyA0` | `nothing` | A0=K: primary-field operator. |
| `modifyRHS` | `nothing` | RHS at the known current state. |
| `nonlinearoptions` | `(;)` | `maxiter=30`, `reltol=1e-9`, `abstol=1e-11`, `relaxation=1.0`. |

Callbacks use `callback(operator,state)` in full global coordinates, with the
same state, return types, private-copy and read-only contracts as [`HHT`](@ref).
Returning `nothing` preserves the accepted operator or prescribed RHS sample.
Displacement updates remain explicit. A consistency iteration handles callback
dependence on acceleration or on the centered damping velocity; it does not
introduce an implicit displacement Newton solve. Nonlinear options control
that residual iteration and require at least one callback. `maxiter` is positive,
tolerances are nonnegative, and `relaxation` lies in `(0,1]`.

# Returns

Return `(u,v)` as global nodal Fields at `t0 + (j-1)*Δt`, `j=1:n`.
Multifield outputs are tuples in block order, for example
`(u1,u2), (v1,v2) = CDM(...)`. Backward free velocities are half-step estimates
stored under these same time labels; choose `velocityScheme=:central` for
sample-time derivatives.

# Choosing the time step

For an undamped linear system, let `ωmax` be the largest free-coordinate
frequency from `K*ϕ = ω²*M*ϕ`. Require **Δt < 2/ωmax** and leave a safety margin.
For classical modal damping with `dampingScheme=:lagged`, a conservative limit is

```math
\\Delta t_{\\mathrm{crit}}
= \\frac{2}{\\omega_{\\max}}
  \\left(\\sqrt{1+\\xi_{\\max}^2}-\\xi_{\\max}\\right),
```

where `ξmax` bounds modal damping ratios over the retained spectrum. Centered
damping keeps the undamped modal limit `2/ωmax` for nonnegative classical
damping. These formulas do not establish stability for arbitrary nonsymmetric
or nonproportionally damped systems. Resolve forcing/support motion as well,
perform step refinement, and reassess the limit when callbacks change operators.
Use [`CDMaccuracyAnalysis`](@ref) with the same `dampingScheme`; the
[time-step guide](@ref transient-timestep) gives an eigenvalue-based example.

# Example

```julia
u, v = CDM(K, M, C, f, u0, v0, n, Δt;
           support=bc, dampingScheme=:central, velocityScheme=:central)
```
"""
function CDM(K::_TransientMatrix, M::_TransientMatrix, C::_TransientMatrix,
             f::_TransientData, U0, V0, n::Int, dt::Real; kwargs...)
    return _CDM(K, M, C, f, U0, V0, n, dt; kwargs...)
end
function CDM(K::_TransientMatrix, M::_TransientMatrix, f::_TransientData, U0, V0, n::Int, dt::Real;
             damping=nothing, kwargs...)
    return _CDM(K, M, damping, f, U0, V0, n, dt; kwargs...)
end
function _CDM(K, M, C, f, U0, V0, n, dt;
    modifyA2=nothing, modifyA1=nothing, modifyA0=nothing, modifyRHS=nothing,
    nonlinearoptions=(;), dampingScheme=:lagged, velocityScheme=:backward, kwargs...)
    _transient_cdm_options(dampingScheme, velocityScheme)
    if _has_transient_modifiers(modifyA2,modifyA1,modifyA0,modifyRHS)
        return _modified_cdm(K,M,C,f,U0,V0,n,dt; modifyA2,modifyA1,modifyA0,modifyRHS,
                             nonlinearoptions,dampingScheme,velocityScheme,kwargs...)
    end
    nonlinearoptions isa NamedTuple && isempty(nonlinearoptions) ||
        throw(ArgumentError("nonlinearoptions require a modify callback and must be a NamedTuple."))
    p = prepare_transient_system(K, M, C, f, U0, V0, n, dt; kwargs...)
    free, fixed = p.free, p.fixed
    Kff, Mff, Cff = p.K[free, free], p.M[free, free], p.C[free, free]
    Kfc, Mfc, Cfc = p.K[free, fixed], p.M[free, fixed], p.C[free, fixed]
    Vd, Ad = _transient_derivatives(p.D, dt)
    U, V = zeros(length(free), n), zeros(length(free), n)
    U[:, 1], V[:, 1] = p.u0, p.v0
    if !isempty(free)
        mass = _transient_factor(Mff; mass=true)
        a0 = mass \ (_transient_column(p.F, 1) - Kff*p.u0 - Cff*p.v0 -
                     Kfc*p.D[:, 1] - Cfc*Vd[:, 1] - Mfc*Ad[:, 1])
        previous = p.u0 - dt*p.v0 + 0.5dt^2*a0
        centered = dampingScheme === :central
        effective = centered ? _transient_factor(Mff + dt/2*Cff) : mass
        for j in 2:n
            u = @view U[:, j-1]
            rhs = _transient_column(p.F, j-1) - Kff*u -
                  Cff*((u-previous)/dt) - Kfc*p.D[:, j-1] -
                  Cfc*Vd[:, j-1] - Mfc*Ad[:, j-1]
            U[:, j] .= 2u - previous + dt^2*(effective \ rhs)
            V[:, j] .= (U[:, j] - u) / dt
            previous = copy(u)
        end
    end
    if velocityScheme === :central
        V, _ = _transient_derivatives(U, dt)
        V[:, 1] .= p.v0
    end
    return reconstruct_transient_solution(U, p.D, p), reconstruct_transient_solution(V, Vd, p)
end

"""
    FDM(K, C, q, X0, n, Δt; kwargs...)

Integrate a first-order system with the **theta method**:

```math
C\\dot{x} + Kx = q.
```

# Arguments

| Argument | Meaning |
|:--|:--|
| `K`, `C` | Primary-field and first-derivative/capacity operators; single- or multifield matrices. |
| `q` | Constant or sampled RHS Field, `SystemVector`, tuple or local/global sum. |
| `X0` | Initial primary field, one column, on free coordinates. |
| `n`, `Δt` | Sample count (`n ≥ 2`, including the initial sample) and positive constant increment. |

The positional-support alias `FDM(K,C,q,bc,X0,n,Δt)` is retained. Here `C` is
capacity/first-derivative data, whereas it denotes damping in HHT/CDM.

# Theta parameter

| Keyword | Default | Meaning |
|:--|:--|:--|
| `ϑ` | `0.5` | Endpoint weight in `[0,1]`; choices below. |

| `ϑ` | Scheme | Time accuracy for smooth solutions |
|:--|:--|:--|
| `0.0` | Explicit Euler | First order; conditionally stable. |
| `0.5` | Crank–Nicolson | Second order. |
| `1.0` | Backward Euler | First order; strong decay of stiff modes. |

# Data, constraints and coordinates

| Keyword | Default | Meaning |
|:--|:--|:--|
| `support` | `BoundaryCondition[]` | Prescribed components and their complete time histories. |
| `mpc` | `MPC[]` | Multi-point constraints, including multifield remote couplings. |
| `coordSys` | `CoordinateSystem[]` | Local bases for loads, initial values, supports and MPC components. |
| `t0` | `0.0` | Time of the initial sample. |

Loads have one constant column or exactly `n` sampled columns. Initial fields
have **one column**. Input time labels are not resampled. `support` is the only
source of prescribed histories: use numbers, spatial `(x,y,z)` functions,
time-dependent `(x,y,z,t)` functions or sampled `ScalarField` components.
Supported values override initial data on constrained DOFs; no external
`applyBoundaryConditions!` call is needed.

Ordinary loads and initial fields use local solver components. Use `Global(...)`
for Cartesian input or combine local/global loads. Output is always global.
For multiple fields, use a `SystemVector` or tuple in matrix block order and
identify `field=Problem` in every BC/MPC. Reduced-order and SEM spaces follow
the shared solver transformations.

Prescribed capacity coupling includes the complete increment
`-C_fc*(d[j+1]-d[j])/Δt`; all histories are supplied through `support`.

# Linear solver options

| Keyword | Default | Meaning |
|:--|:--|:--|
| `solver` | `:auto` | Select a direct or iterative solve; choices below. |
| `preconditioner` | `:auto` | `:ich` for CG, `:ilu` for GMRES; also `:none` or a supplied object. |
| `solveroptions` | `(;)` | Native iterative options, e.g. `(reltol=1e-9, maxiter=1000)`. |
| `preconditioneroptions` | `(;)` | Options for a symbolic preconditioner, e.g. `(memory=2,)` for `:ich`. |

- Direct choices: `:auto`, `:backslash`, `:lu`, `:lu_no_ordering`,
  `:cholesky`, `:ldlt`, `:qr`.
- Iterative choices: `:cg` (positive-definite symmetric systems) and `:gmres`.
  Load `IterativeSolvers`; `:ich` additionally needs `Preconditioners`,
  and `:ilu` needs `IncompleteLU`.
- Wrap every contributing operator in `Symmetric(...)` to declare the effective
  matrix symmetric. Then `:auto` tries Cholesky with LU fallback; otherwise it
  uses the general sparse factorization. Symmetry alone does not imply positive definiteness.
- Nonempty option tuples and actual preconditioners are for iterative solves.
  `log`, `Pl`, `Pr` and `initially_zero` are managed internally. Constructor
  options cannot accompany `:none` or a preconstructed preconditioner.

For unchanged coefficients, the effective factorization/preconditioner is reused.
See [`solveField`](@ref) for the shared backend contract.

For constant operators the effective free matrix is `C_ff + ϑ*Δt*K_ff`.
Wrap both operators in `Symmetric(...)` for the symmetric automatic path.

# Operator updates and nonlinear iteration

| Keyword | Default | Updated quantity |
|:--|:--|:--|
| `modifyA1` | `nothing` | A1=C: first-derivative/capacity operator. |
| `modifyA0` | `nothing` | A0=K: primary-field operator or tangent. |
| `modifyRHS` | `nothing` | RHS at the trial endpoint. |
| `modifyA2` | `nothing` | Reserved: any value other than `nothing` is rejected by FDM. |
| `nonlinearoptions` | `(;)` | Nonlinear iteration controls listed below. |

FDM reevaluates callbacks within each step. It interpolates capacity as
`A1θ=(1-ϑ)*A1_old+ϑ*A1_new` and uses each endpoint's own `A0*x-rhs`.
There is no second-derivative operator.

The callback signature is `callback(operator, state)`. Matrices are full/global
`SystemMatrix`s (optionally `Symmetric`); RHS results are one-column global Fields,
field tuples or `SystemVector`s. `state` includes trial fields and derivatives,
Problems/material data, time, step, iteration, previous accepted snapshot and
current operators. Callback order is **A2, A0, A1, RHS**.
See the [complete callback guide](@ref transient-callbacks).

The operator argument is a private copy. Return it to keep in-place edits;
return `nothing` to keep the accepted matrix or the original prescribed RHS
sample. State data are read-only. Mesh, field ordering, supports, MPCs and bases
stay fixed; inactive MPC DOFs cannot be activated.

| `nonlinearoptions` member | Default | Meaning |
|:--|:--|:--|
| `maxiter` | `30` | Maximum nonlinear updates; positive integer. |
| `reltol` | `1e-9` | Relative residual tolerance; nonnegative. |
| `abstol` | `1e-11` | Absolute residual tolerance, in equilibrium-equation units; nonnegative. |
| `relaxation` | `1.0` | Fraction of each correction, in `(0,1]`. |

For Newton iteration, return the tangent `J(u)` as A0 and
`external - internal(u) + J(u)*u` as RHS. State-dependent A2/A1 also require
consistent linearization for Newton convergence. Picard updates are permitted.
Nonconvergence throws a step/time/residual error. Without callbacks, leave
`nonlinearoptions=(;)` to use the linear path.

At FDM initialization both derivative entries in `state` are `nothing`;
a first derivative is available during step iterations. The same coordinate,
reduced-order and SEM restrictions as in [`HHT`](@ref) apply.

# Returns

Return a global nodal Field `x` with `n` columns at `t0 + (j-1)*Δt`.
A multifield result is a tuple in matrix block order: `x1,x2 = FDM(...)`.

# Choosing the time step

For a linear diffusion-type pencil `K*ϕ=λ*C*ϕ` with real nonnegative eigenvalues:

- `ϑ ≥ 0.5`: unconditionally stable, but `Δt` still controls accuracy.
- `ϑ < 0.5`: require `Δt < 2/((1-2ϑ)*λmax)` with a safety margin.

For complex/nonsymmetric spectra, assess the full stability region instead.
Crank–Nicolson may alternate stiff modes when `λ*Δt` is large; backward Euler
damps them more strongly. Resolve the fastest decay/forcing of interest and
halve `Δt` to check accuracy. See [`FDMaccuracyAnalysis`](@ref).

# Example

```julia
x = FDM(K, C, q, x0, n, Δt; support=bc, ϑ=0.5)
x = FDM(Symmetric(K), Symmetric(C), q, x0, n, Δt;
        support=bc, ϑ=1.0, solver=:cholesky)
```
"""
function FDM(K::_TransientMatrix, C::_TransientMatrix, q, X0, n::Int, dt::Real;
             ϑ=0.5, solver=:auto, preconditioner=:auto,
             solveroptions=(;), preconditioneroptions=(;),
             modifyA2=nothing, modifyA1=nothing, modifyA0=nothing, modifyRHS=nothing,
             nonlinearoptions=(;), kwargs...)
    isfinite(ϑ) && 0 <= ϑ <= 1 || throw(ArgumentError("FDM: require 0 ≤ ϑ ≤ 1."))
    modifyA2 === nothing || throw(ArgumentError("FDM has no A2 operator; use modifyA1 for capacity."))
    if _has_transient_modifiers(modifyA2,modifyA1,modifyA0,modifyRHS)
        return _modified_fdm(K,C,q,X0,n,dt; ϑ,solver,preconditioner,solveroptions,
            preconditioneroptions,modifyA1,modifyA0,modifyRHS,nonlinearoptions,kwargs...)
    end
    nonlinearoptions isa NamedTuple && isempty(nonlinearoptions) ||
        throw(ArgumentError("nonlinearoptions require a modify callback and must be a NamedTuple."))
    p = prepare_transient_system(K, nothing, C, q, X0, nothing, n, dt; kwargs...)
    free, fixed = p.free, p.fixed
    Kff, Cff = p.K[free, free], p.C[free, free]
    Kfc, Cfc = p.K[free, fixed], p.C[free, fixed]
    U = zeros(length(free), n)
    U[:, 1] .= p.u0
    if !isempty(free)
        effective = _prepare_transient_linear_solver(Cff + ϑ*dt*Kff, (K, C);
            solver, preconditioner, solveroptions, preconditioneroptions)
        previous = Cff - (1-ϑ)*dt*Kff
        for j in 2:n
            rhs = previous*U[:, j-1] +
                  dt*((1-ϑ)*_transient_column(p.F, j-1) + ϑ*_transient_column(p.F, j)) -
                  dt*Kfc*((1-ϑ)*p.D[:, j-1] + ϑ*p.D[:, j]) - Cfc*(p.D[:, j]-p.D[:, j-1])
            U[:, j] .= solve_prepared_linear_system(effective, rhs)
        end
    end
    return reconstruct_transient_solution(U, p.D, p)
end

# Historical positional BC lists route through exactly the same pipeline.
for method in (:HHT, :CDM)
    @eval begin
        $method(K::_TransientMatrix, M::_TransientMatrix, f::_TransientData, bc::Vector{BoundaryCondition},
                U0, V0, n::Int, dt::Real; kwargs...) =
            $method(K, M, f, U0, V0, n, dt; support=bc, kwargs...)
        $method(K::_TransientMatrix, M::_TransientMatrix, C::_TransientMatrix, f::_TransientData,
                bc::Vector{BoundaryCondition}, U0, V0, n::Int, dt::Real; kwargs...) =
            $method(K, M, C, f, U0, V0, n, dt; support=bc, kwargs...)
    end
end
FDM(K::_TransientMatrix, C::_TransientMatrix, q, bc::Vector{BoundaryCondition},
    X0, n::Int, dt::Real; kwargs...) = FDM(K, C, q, X0, n, dt; support=bc, kwargs...)
FDM_old(args...; kwargs...) = FDM(args...; kwargs...)

"""
    CDMaccuracyAnalysis(ωmin, ωmax, Δt, type; kwargs...)

Evaluate a **scalar modal amplification model** of [`CDM`](@ref).
Frequencies are angular frequencies from `K*ϕ=ω²*M*ϕ`.

# Arguments and output

- `ωmin`, `ωmax`: positive frequency interval in radians/time.
- `Δt`: positive candidate time increment.
- `type`: analysis listed below.
- Return `(x,y)` vectors of length `n`, for plotting without a GUI dependency.

| `type` | `y` | `x` |
|:--|:--|:--|
| `:SR` | Spectral radius of the amplification matrix. | Natural log of `ω*Δt/(2π)`. |
| `:PDR` | Physical damping contribution from the β series. | `ω*Δt/(2π)`. |
| `:ADR` | Algorithmic damping indicator. | `ω*Δt/(2π)`. |
| `:PE` | Relative period-error indicator. | `ω*Δt/(2π)`. |

# Keyword arguments

| Keyword | Default | Meaning |
|:--|:--|:--|
| `n` | `100` | Number of frequency samples. |
| `dampingScheme` | `:lagged` | Match the integrator: `:lagged` or `:central`. |
| `α` | `0.0` | Legacy mass-proportional damping parameter; not HHT algorithmic α. |
| `ξ` | `0.01` | Scalar/vector of β-series damping contributions at `ωmax`. |
| `β` | Derived from `ξ` | Explicit Caughey-series coefficients; overrides the derived values. |
| `show_β` | `false` | Print β coefficients. |
| `show_ξ` | `false` | Print the corresponding β-series contributions at `ωmax`. |

The derived coefficients are `β[i]=2ξ[i]/ωmax^(2i-1)`. The amplification
matrix uses the existing convention
`ξeff(ω)=α/ω + sum(β[i]*ω^(2i-1)/2)`. Thus, for a physical mass term `a*M`,
use `α=a/2`. The `:PDR` curve and `show_ξ` report only the β-series contribution;
they do not include the mass term. Explicit β permits stiffness-proportional
or higher-order Caughey damping models.

# Use and limitations

Use `:SR` to inspect a candidate stable step; use `:PE` and `:ADR` to judge
numerical distortion. Select the same damping scheme as the integration run.
This is a modal/proportional model, not a stability proof for an arbitrary
coupled or nonproportionally damped matrix system.

```julia
x, radius = CDMaccuracyAnalysis(ωmin, ωmax, Δt, :SR;
                              ξ=0.02, dampingScheme=:central)
```

See Serfőző and Pere, *An effective reduction method with Caughey damping for
spurious oscillations in dynamic problems*, Meccanica (2025),
<https://doi.org/10.1007/s11012-025-02036-9>.
"""
function CDMaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, ξ=0.01, β=[2ξ[i]/(ωₘₐₓ)^(2i-1) for i in 1:length(ξ)], show_β=false, show_ξ=false, dampingScheme=:lagged)
    _transient_cdm_options(dampingScheme, :backward)
    if show_β == true
        println("β = $β")
    end
    if show_ξ == true
        ξ = [β[i] / 2 * ωₘₐₓ^((2i-1)) for i in 1:length(β)]
        println("ξ = $ξ")
    end
    #Tₘᵢₙ /= √(1-ξₘₐₓ^2)
    x = zeros(n)
    y = zeros(n)
    ω = range(ωₘᵢₙ, length=n, stop=ωₘₐₓ)
    for i ∈ 1:n
        ξ = α/ω[i]
        for j in 1:length(β)
            ξ += β[j] / 2 * ω[i]^(2j-1)
        end
        Ω = Δt * ω[i]
        A = dampingScheme === :central ?
            [(2-Ω^2)/(1+ξ*Ω) (ξ*Ω-1)/(1+ξ*Ω); 1 0] :
            [2-2ξ*Ω-Ω^2 2ξ*Ω-1; 1 0]

        eig = eigen(A)
        ρ, idx = findmax(abs, eig.values)
        λ = eig.values[idx]
        σ = real(λ)
        ε = imag(λ)
        if type == :SR
            x[i] = log((ω[i] / 2π) * Δt)
            y[i] = ρ
        elseif type == :ADR
            x[i] = (ω[i] / 2π) * Δt
            Ω0 = √((log(ρ))^2 + (atan(ε,σ))^2 / 4)
            y[i] = -log(ρ) / 2Ω0
        elseif type == :PDR
            x[i] = (ω[i] / 2π) * Δt
            for j in 1:length(β)
                y[i] += β[j] / 2 * (2π * x[i] / Δt) ^ (2j-1)
            end
        elseif type == :PE
            x[i] = (ω[i] / 2π) * Δt
            Ω0 = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = 1 - Ω0/(Δt*ω[i])
        else
            str1 = "CDMaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\n:SR - spectral radius\n"
            str4 = ":PDR - physical damping ratio\n"
            str5 = ":ADR - algorithmic damping ratio\n"
            str6 = ":PE - period error\n"
            str7 = "\nFor details see Serfőző, D., Pere, B.: A method to accurately define arbitrary\n"
            str8 = "algorithmic damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).\n"
            str9 = "https://doi.org/10.1007/s00419-023-02454-9\n"
            error(str1*str2*str3*str4*str5*str6*str7*str8*str9)
        end
    end
    return x, y
end

"""
    HHTaccuracyAnalysis(ωmin, ωmax, Δt, type; kwargs...)

Evaluate the **HHT modal amplification matrix**, including physical damping.
Use it to compare candidate time increments and HHT parameter choices.

# Arguments and output

- `ωmin`, `ωmax`: positive angular-frequency interval from `K*ϕ=ω²*M*ϕ`.
- `Δt`: finite positive candidate increment.
- `type`: analysis listed below.
- Return `(x,y)` vectors of length `n`.

| `type` | `y` | `x` |
|:--|:--|:--|
| `:SR` | Amplification spectral radius. | Natural log of `ω*Δt/(2π)`. |
| `:PDR` | Total physical modal damping ratio. | `ω*Δt/(2π)`. |
| `:ADR` | Algorithmic damping indicator. | `ω*Δt/(2π)`. |
| `:PE` | Relative period-error indicator. | `ω*Δt/(2π)`. |

# Keyword arguments

| Keyword | Default | Meaning |
|:--|:--|:--|
| `n` | `100` | Number of samples; at least 2. |
| `α` | `0.0` | HHT algorithmic damping in `[-1/3,0]`. |
| `δ` | `nothing` | Optional `γ=0.5+δ` relation. |
| `γ` | `0.5-α` | Standard default, or `0.5+δ` when δ is supplied. |
| `β` | `(0.5+γ)^2/4` | Displacement parameter; explicit values override the default. |
| `ξ` | `0.0` | Physical modal ratio, a number or function `ξ(ω)`. |
| `massDamping` | `0.0` | Coefficient of M in Rayleigh damping. |
| `stiffnessDamping` | `0.0` | Coefficient of K in Rayleigh damping. |

The physical modal ratio is
`ξ(ω) + massDamping/(2ω) + stiffnessDamping*ω/2`.
These terms affect the amplification matrix; `:PDR` reports their sum.
HHT defaults and explicit parameter overrides follow [`HHT`](@ref).

# Example and limitations

```julia
x, radius = HHTaccuracyAnalysis(ωmin, ωmax, Δt, :SR;
                               α=-0.1, ξ=0.02, stiffnessDamping=1e-4)
```

Use the same coefficients as the transient calculation. This scalar analysis
assumes modal/proportional damping; an arbitrary nonproportional matrix cannot
be reconstructed from the undamped frequency interval alone. Stability of the
standard linear scheme does not replace a time-step refinement check.

See Belytschko and Hughes, *Computational Methods for Transient Analysis*
(North-Holland, 1983), and Serfőző and Pere, *A method to accurately define
arbitrary algorithmic damping character as viscous damping*, Arch. Appl. Mech.
93 (2023), <https://doi.org/10.1007/s00419-023-02454-9>.
"""
function HHTaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, δ=nothing,
    γ=0.5 + (δ === nothing ? -α : δ), β=0.25 * (0.5 + γ)^2, ξ=0.0, massDamping=0.0, stiffnessDamping=0.0)
    _transient_accuracy_check(ωₘᵢₙ, ωₘₐₓ, Δt, n)
    all(isfinite, (α, γ, β, massDamping, stiffnessDamping)) &&
        -1/3 <= α <= 0 && γ >= 0 && β > 0 ||
        throw(ArgumentError("HHTaccuracyAnalysis: invalid HHT/damping parameters."))
    type in (:SR, :ADR, :PE, :PDR) || throw(ArgumentError("HHTaccuracyAnalysis: use :SR, :ADR, :PE or :PDR."))
    x = zeros(n)
    y = similar(x)
    invT = range(ωₘᵢₙ/2π, length=n, stop=ωₘₐₓ/2π)
    for i ∈ 1:n
        ω = 2π * invT[i]
        ratio = ξ isa Function ? ξ(ω) : ξ
        isfinite(ratio) || throw(ArgumentError("HHTaccuracyAnalysis: non-finite damping ratio."))
        c = 2ratio*ω + massDamping + stiffnessDamping*ω^2
        A1 = [1 0 -Δt^2*β
        0 1 -Δt*γ
        (1+α)*ω^2 (1+α)*c 1]
        A2 = [1 Δt Δt^2*(0.5-β)
        0 1 Δt*(1-γ)
        α*ω^2 α*c 0]

        A = A1 \ A2

        eig = eigen(A)
        ρ, idx = findmax(abs, eig.values)
        λ = eig.values[idx]
        σ = real(λ)
        ε = imag(λ)
        if type == :SR
            x[i] = log(invT[i] * Δt)
            y[i] = ρ
        elseif type == :PDR
            x[i] = invT[i] * Δt
            y[i] = c / (2ω)
        elseif type == :ADR
            x[i] = invT[i] * Δt
            Ω = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = -log(ρ) / 2Ω
            #y[i] = -log(ρ) / atan(ε, σ)
        elseif type == :PE
            x[i] = invT[i] * Δt
            Ω = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = 1 - Ω/(2π*Δt*invT[i])
        else
            str1 = "HHTaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\n:SR - spectral radius\n"
            str5 = ":ADR - algorithmic damping ratio\n"
            str6 = ":PE - period error\n"
            str7 = "\nFor details see Serfőző, D., Pere, B.: A method to accurately define arbitrary\n"
            str8 = "algorithmic damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).\n"
            str9 = "https://doi.org/10.1007/s00419-023-02454-9\n"
            error(str1*str2*str3*str5*str6*str7*str8*str9)
        end
    end
    return x, y
end


"""
    FDMaccuracyAnalysis(λmin, λmax, Δt; type=:SR, n=100, ϑ=0.5)

Evaluate the **theta-method amplification factor** for [`FDM`](@ref).

# Arguments

| Argument | Meaning |
|:--|:--|
| `λmin`, `λmax` | Nonnegative real interval from `K*ϕ=λ*C*ϕ`. |
| `Δt` | Finite positive candidate time increment. |

# Keyword arguments

| Keyword | Default | Meaning |
|:--|:--|:--|
| `type` | `:SR` | Spectral-radius analysis; the only supported choice. |
| `n` | `100` | Number of eigenvalue samples; at least 2. |
| `ϑ` | `0.5` | Theta weight in `[0,1]`, matching the integration run. |

# Returns and interpretation

Return `(x,radius)` vectors, with `x=λ*Δt` and
`radius=abs((1-(1-ϑ)*x)/(1+ϑ*x))`. Values above one indicate modal growth.
For `ϑ<0.5`, the nonnegative-real stability boundary is
`λmax*Δt < 2/(1-2ϑ)`. Complex spectra require a full stability-region check.

```julia
x, radius = FDMaccuracyAnalysis(0.0, λmax, Δt; ϑ=0.5)
```
"""
function FDMaccuracyAnalysis(λₘᵢₙ, λₘₐₓ, Δt; type=:SR, n=100, ϑ=0.5)
    _transient_accuracy_check(λₘᵢₙ, λₘₐₓ, Δt, n; positive=false)
    isfinite(ϑ) && 0 <= ϑ <= 1 ||
        throw(ArgumentError("FDMaccuracyAnalysis: require 0 ≤ ϑ ≤ 1."))
    type === :SR || throw(ArgumentError("FDMaccuracyAnalysis: only :SR is supported."))
    x = collect(range(λₘᵢₙ, λₘₐₓ; length=n)) .* Δt
    radius = abs.((1 .- (1-ϑ).*x) ./ (1 .+ ϑ.*x))
    return x, radius
end


function _transient_accuracy_check(lo, hi, dt, n; positive=true)
    isfinite(lo) && isfinite(hi) && lo <= hi && (positive ? lo > 0 : lo >= 0) ||
        throw(ArgumentError("Accuracy analysis: invalid spectral interval."))
    isfinite(dt) && dt > 0 && n >= 1 ||
        throw(ArgumentError("Accuracy analysis: require finite Δt > 0 and n ≥ 1."))
end
