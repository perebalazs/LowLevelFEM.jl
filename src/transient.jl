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
    HHT(K, M, [C,] f, U0, V0, n, Δt; support=[], mpc=[], coordSys=[],
        α=0.0, δ=nothing, γ=0.5+(δ === nothing ? -α : δ), β=0.25*(0.5+γ)^2, t0=0.0,
        solver=:auto, preconditioner=:auto, solveroptions=(;), preconditioneroptions=(;),
        modifyA2=nothing, modifyA1=nothing, modifyA0=nothing, modifyRHS=nothing,
        nonlinearoptions=(;))

Integrate `M*ü + C*u̇ + K*u = f` with HHT/Newmark. Omitting `C` means zero
physical damping; `α` controls algorithmic damping. The default `δ=nothing`
selects the standard HHT choice `γ=0.5-α`, `β=(1-α)^2/4`. Explicit `δ`, `γ`
and `β` override that choice; `δ=0` recovers the historical parameter family. The equilibrium equation is
`M*a[j+1] + (1+α)*(C*v[j+1]+K*u[j+1]-f[j+1]) -
α*(C*v[j]+K*u[j]-f[j]) = 0` in free coordinates. Prescribed inertia,
damping and stiffness are first transferred to the effective RHS and receive
the same endpoint weights as applied loads.

All transient solvers use `Q*T`: local coordinates, intrinsic reducedOrder/SEM
approximation, then MPC constraints, following `solveField`. The existing
restriction on combining coordinates and an approximation on the same field
also applies here. Loads and initial fields are in local solver components;
`Global(...)` accepts Cartesian data. Results are always global nodal Fields.

`f` accepts a Field, a multifield `SystemVector`/tuple, `Global(...)`, or a
local/global sum. Its columns are constant (1 column) or sampled at
`t0 + (j-1)*Δt` (`n` columns); input Field time labels are not resampled.
`U0` and `V0` contain one initial column. The complete boundary history is
assembled internally from `support`, including time-dependent ScalarField
components and `(x,y,z,t)` functions. No external applyBoundaryConditions!
call is needed. Boundary data override initial values on supported DOFs;
`V0` supplies initial free velocities only.

Constrained velocities and accelerations always follow the prescribed primary
history through second-order finite differences at every sample; prescribed
DOFs are not advanced with Newmark recurrence. At least three samples are needed to estimate
nonzero initial support acceleration. To prescribe a velocity history,
integrate it and put the resulting scalar components into BoundaryConditions.

Return `(u, v)`. For multiple fields each of `u` and `v` is a tuple in matrix
block order, e.g. `(u1,u2), (v1,v2) = HHT(...)`. Initial values on constrained
DOFs are replaced by prescribed histories. The factorization is reused.
The positional `..., f, support, U0, V0, n, Δt` form is also supported.
HHT/CDM require a nonsingular mass on free reduced coordinates; arbitrary
mass-singular differential-algebraic systems are outside this interface.
SEM plus MPC on the same field follows the shared backend's current rejection.

The generic equation is `A2*ü + A1*u̇ + A0*u = rhs`, with the positional
aliases `A0=K`, `A2=M`, `A1=C`. The optional `modifyA2`, `modifyA1`,
`modifyA0`, `modifyRHS` callbacks default to `nothing`. Each is called as
`callback(operator, state)` in the full global field space. See the transient
callback documentation in `docs/src/reference/solution.md` for the state and
linearization contract. Matrix callbacks return a SystemMatrix (optionally
wrapped in Symmetric); modifyRHS returns a one-column global Field/SystemVector.
Returning `nothing` keeps the previous accepted matrix, or the prescribed RHS
sample at the evaluation time. Callback inputs are private copies; return the
copy to retain in-place edits. `state` data must be treated as read-only.

HHT reevaluates callbacks inside each time step until the discrete equilibrium
residual converges. `nonlinearoptions` accepts `maxiter=30`, `reltol=1e-9`,
`abstol=1e-11`, and `relaxation=1.0`. Nonconvergence throws an error with step,
time and residual. For Newton iteration, return a consistent tangent in A0 and
`rhs = external - internal + A0*u` in modifyRHS. Otherwise the same loop can
perform a Picard iteration. State-dependent A2/A1 need a consistent user
linearization for Newton convergence. Fixed layout, supports, MPCs, coordinates
and approximation transformations are retained; inactive MPC DOFs cannot be
activated by updates. Equal matrices reuse the cached linear solver.

Linear solvers follow `solveField`: `:auto`, `:backslash`, `:lu`,
`:lu_no_ordering`, `:cholesky`, `:ldlt`, `:qr`, `:cg` and `:gmres`.
`preconditioner`, `solveroptions` and `preconditioneroptions` have the same
meaning and optional-package requirements as in `solveField`. The effective
matrix factorization or preconditioner is prepared once and reused at every
step. For a nondiagonal mass, the initial acceleration uses a separate solve
with the same settings; a diagonal mass uses direct division.
Wrap all contributing matrices in `Symmetric(...)` to declare the effective
system symmetric, including `C` when supplied. Then `:auto` tries Cholesky
with LU fallback, and `:cg` with `preconditioner=:auto` uses incomplete
Cholesky. General systems use sparse backslash by default, while `:gmres`
with automatic preconditioning uses ILU. CG requires positive definiteness;
a symmetric wrapper alone does not guarantee it.
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
    CDM(K, M, [C,] f, U0, V0, n, Δt; support=[], mpc=[], coordSys=[],
        dampingScheme=:lagged, velocityScheme=:backward,
        modifyA2=nothing, modifyA1=nothing, modifyA0=nothing, modifyRHS=nothing,
        nonlinearoptions=(;))

Central differences with `dampingScheme=:lagged` (historical backward damping)
or `:central`, which uses `C*(u[j+1]-u[j-1])/(2Δt)` and factors `M+Δt*C/2`.
The central choice is second-order in damping and avoids the additional lagged
stability restriction. Accepts the same
load, initial-value, coordinate, approximation, MPC and prescribed-history
keywords as [`HHT`](@ref). Returns `(u,v)` with multifield tuple grouping as in
HHT. A diagonal mass uses componentwise division. A consistent or projected
nondiagonal mass uses a cached factorization, without automatic lumping.
The scheme is conditionally stable; choose Δt for the reduced dynamic system.
Constrained derivatives are estimated from the U0/support displacement history
with second-order finite differences. `velocityScheme=:backward` preserves
historical free velocities at half-step locations (the first sample is V0).
Choose `velocityScheme=:central` for second-order velocities at the returned
sample times, using a one-sided final stencil and preserving the initial V0.
Operator callbacks use the same global-space contract as HHT. Displacement
updates remain explicit: callbacks are evaluated at the known current field.
An acceleration consistency iteration handles any callback dependence on the
current acceleration; no implicit displacement Newton solve is introduced.
Updated coefficients must satisfy the CDM time-step stability requirement.
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
    FDM(K, C, q, X0, n, Δt; ϑ=0.5, support=[], mpc=[], coordSys=[],
        t0=0.0, solver=:auto, preconditioner=:auto, solveroptions=(;), preconditioneroptions=(;),
        modifyA1=nothing, modifyA0=nothing, modifyRHS=nothing, nonlinearoptions=(;))

Theta integration of `C*ẋ + K*x = q`, sharing the constrained preparation and
sampled-data conventions of [`HHT`](@ref). `X0` has one initial column and
all boundary histories are supplied by `support`. Includes the complete prescribed
increment `-C_fc*(d[j+1]-d[j])/Δt` and theta-weighted loads and stiffness terms.
Return one global nodal Field, or a tuple in matrix block order. `ϑ=0` uses
explicit Euler, `ϑ=0.5` Crank–Nicolson and `ϑ=1` backward Euler. The effective
matrix is factorized once; diagonal effective matrices use direct division.
The positional `..., q, support, X0, n, Δt` form is retained.
Linear solver and preconditioner options follow `solveField` and HHT. The
factorization or preconditioner of `C_ff + ϑ*Δt*K_ff` is built once and reused.
Wrap both `K` and `C` in `Symmetric(...)` for the symmetric automatic path.
Generic operators are A0=K and A1=C (the first-derivative/capacity matrix).
modifyA2 is rejected for this first-order method. With modify callbacks, the
time-step loop checks the nonlinear theta residual. A1 is interpolated as
`(1-ϑ)*A1_old + ϑ*A1_new`; the old and new A0*u-rhs terms use their own
endpoint operators. Callback and nonlinear options follow HHT.
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
    CDMaccuracyAnalysis(ωₘᵢₙ::Float64,
                        ωₘₐₓ::Float64,
                        Δt::Float64,
                        type::Symbol;
                        n=100,
                        α=0.0,
                        ξ=0.01,
                        β=[2ξ[i]/(ωₘₐₓ)^(2i-1) for i in 1:length(ξ)],
                        show_β=false,
                        show_ξ=false)

Gives some functions (graphs) for accuracy analysis of the CDM method.
`ωₘᵢₙ` and `ωₘₐₓ` are the square root of smallest and largest eigenvalues of the
**Kϕ**=ω²**Mϕ** eigenvalue problem, `Δt` is the time step size. `type` is one of the
following values:
- `:SR`: spectral radius
- `:PDR`: physical damping ratio
- `:ADR`: algorithmic damping ratio
- `:PE`: period error
For details see [^3]. `dampingScheme=:lagged` matches historical CDM;
`:central` matches its centered-damping variant.
`n` is the number of points in the graph. The damping matrix is assembled in the
following ways: **C**=α**M**+β**K** or **C**=α**M**+β₁**K**+β₂**KM⁻¹K**+β₃**KM⁻¹KM⁻¹K**+⋅⋅⋅.
The latter corresponds to the damping characteristic characterized by a power series consisting of powers
of the natural frequencies with odd exponents. ξᵢ (`ξ` in the argument list) are the values ​​of the
individual members of the series corresponding to the ωₘₐₓ value. βᵢ (`β` in the argument list) are the
coefficients of the series. (see [^4]) Either `ξ` or `β` must be specified. `ξ` or `β` are scalars or
vectors. If `show_β` or `show_ξ` is `true`, the corresponding `β` or `ξ` values will be
sent to the output.
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

[^4]: Serfőző, D., Pere, B.: *An effective reduction method with Caughey damping for
    spurious oscillations in dynamic problems*, Meccanica, <https://doi.org/10.1007/s11012-025-02036-9>

Return: `xy`

Types:
- `ωₘᵢₙ`: Float64
- `ωₘₐₓ`: Float64
- `Δt`: Float64
- `n`: Int64
- `α`: Float64
- `β`: Float64 of Vector{Float64}
- `ξ`: Float64 of Vector{Float64}
- `show_β`: Boolean
- `show_ξ`: Boolean
- `xy`: Tuple{Vector{Float64},Vector{Float64}}
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
    HHTaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, δ=nothing, γ=0.5 + (δ === nothing ? -α : δ), β=0.25 * (0.5 + γ)^2)

Gives some functions (graphs) for accuracy analysis of the HHT-α method described by [`HHT`](@ref).
`ωₘᵢₙ` and `ωₘₐₓ` are the square root of smallest and largest eigenvalues of the
**Kϕ**=ω²**Mϕ** eigenvalue problem, `Δt` is the time step size. `type` is one of the
following values:
- `:SR`: spectral radius
- `:ADR`: algorithmic damping ratio
- `:PE`: period error
For details see [^2] and [^3].

Physical damping is included in the amplification matrix. `ξ` is a modal damping
ratio (a number or function of angular frequency). For Rayleigh damping
`C = massDamping*M + stiffnessDamping*K`, supply the two coefficients; their
contribution is added to `ξ`. `:PDR` returns the resulting physical damping ratio.
A general nonproportional damping matrix cannot be characterized from the
undamped frequency interval alone; the time integrator accepts that matrix.

`n` is the number of points in the graph. For the meaning of `α`, `β` and `γ`
see [`HHT`](@ref). If `δ` is given, γ=0.5+δ and β=0.25⋅(0.5+γ)².
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

[^2]: Belytschko, Ted, and Thomas JR, Hughes: *Computational methods for
    transient analysis*, North-Holland, (1983).

[^3]: Serfőző, D., Pere, B.: *A method to accurately define arbitrary algorithmic
    damping character as viscous damping*. Arch Appl Mech 93, 3581–3595 (2023).
    <https://doi.org/10.1007/s00419-023-02454-9>

Return: `xy`

Types:
- `ωₘᵢₙ`: Float64
- `ωₘₐₓ`: Float64
- `Δt`: Float64
- `n`: Int64
- `α`: Float64
- `β`: Float64
- `γ`: Float64
- `δ`: Float64
- `xy`: Tuple{Vector{Float64},Vector{Float64}}
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

Return `(x, radius)` for the theta method on `C*ẋ + K*x = 0`, where
`x = λ*Δt` and `λ` ranges over the generalized spectrum `K*ϕ = λ*C*ϕ`.
The amplification factor is `(1-(1-ϑ)*x)/(1+ϑ*x)`. Only `type=:SR` is
supported. See [`FDM`](@ref) for the time-integration method.
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
