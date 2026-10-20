###############################################################################
# eigen.jl — generalized eigenvalue problems (single-field and multifield)
#
# Preparation:   coordinates Q -> approximation/MPC T -> homogeneous BCs
# Algebra:       A*x = λ*B*x
# Reconstruction: phi_global = Q*T*phi_reduced
#
# No physical interpretation is attached to λ. In particular, negative and
# complex eigenvalues/eigenvectors must be preserved.
###############################################################################

export solveEigen, solveEigenProblem
export solveEigenModes, solveEigenFields, solveBucklingModes
export solveModalAnalysis, solveBuckling
export smallestEigenValue, largestEigenValue
export smallestPeriodTime, largestPeriodTime

# Existing symmetric wrapper retains the owning SystemMatrix metadata.
_eigen_parent(K::SystemMatrix) = K
_eigen_parent(K::SymmetricSystemMatrix) = K.parent

_eigen_symmetric(K) = K isa SymmetricSystemMatrix

"""
    prepare_eigen_system(A, B; support=BoundaryCondition[], mpc=MPC[], coordSys=[])

Prepare `A*phi = λ*B*phi` in a common constrained solver space.

Coordinate systems follow `solveField`: the same `coordSys` objects, the same
coordinate builders, and the same order of transformations (`Q*T`). Dirichlet
constraints prescribe zero *modal perturbations*, regardless of their static
values. For mixed-field systems, BCs and MPCs must identify their `Problem`.

Return a named tuple with reduced matrices `A`, `B`, transformation `T`,
coordinate matrix `Q`, the active reduced DOFs, and the full problem dimension.
"""
function prepare_eigen_system(
    A0::Union{SystemMatrix,SymmetricSystemMatrix},
    B0::Union{SystemMatrix,SymmetricSystemMatrix};
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[],
    coordSys=[]
)
    A = _eigen_parent(A0)
    B = _eigen_parent(B0)

    size(A.A) == size(B.A) ||
        error("solveEigen: incompatible matrix sizes $(size(A.A)) and $(size(B.A)).")
    size(A.A, 1) == size(A.A, 2) ||
        error("solveEigen: matrices must be square.")
    size(A.A, 1) > 0 || error("solveEigen: the system is empty.")

    if A.problems === nothing
        B.problems === nothing ||
            error("solveEigen: both matrices must be single-field or both multifield.")
        A.model !== nothing || error("solveEigen: missing Problem metadata in A.")
        A.model === B.model ||
            error("solveEigen: A and B must refer to the same Problem object.")
        A.test_model === A.model && B.test_model === B.model ||
            error("solveEigen: a single-field eigenproblem requires square test/trial field metadata.")
    else
        B.problems !== nothing ||
            error("solveEigen: both matrices must be single-field or both multifield.")
        check_multifield_system_compatibility(A, B)
    end

    nfull = size(A.A, 1)
    symA = _eigen_symmetric(A0)
    symB = _eigen_symmetric(B0)

    # Fast path: full-order Lagrange, no MPC, no local coordinates. Do not
    # allocate any global identity or project large sparse matrices needlessly.
    problems = A.problems === nothing ? [A.model] : A.problems
    is_plain = isempty(mpc) && isempty(coordSys) &&
               all(P -> !_has_approximation(P), problems)

    if is_plain
        free = A.problems === nothing ?
            freeDoFs(A.model, support) :
            first(multifield_bc_data(A, support; nsteps=1))

        isempty(free) && error("solveEigen: no free degrees of freedom.")
        return (
            A=A.A[free, free], B=B.A[free, free], T=nothing,
            Q=nothing, free=free, nfull=nfull,
            symmetricA=symA, symmetricB=symB,
            uploA=symA ? A0.uplo : :U,
            uploB=symB ? B0.uplo : :U
        )
    end

    # Local coordinate systems are defined in solve.jl. Never maintain a
    # second implementation of their nodal directions or multifield mapping.
    Q = nothing
    if !isempty(coordSys)
        if A.problems === nothing
            isdefined(@__MODULE__, :_build_coordinate_transformation) ||
                error("solveEigen: coordinate builders from solve.jl are unavailable.")
            P = A.model
            _has_approximation(P) &&
                error("solveEigen: coordSys + reduced/spectral basis on the same field is not yet implemented.")
            QQ = _build_coordinate_transformation(P, coordSys)
            Q = QQ === nothing ? nothing : QQ.T
        else
            isdefined(@__MODULE__, :_build_multifield_coordinate_transformation) ||
                error("solveEigen: multifield coordinate builders from solve.jl are unavailable.")
            # The solveField builder validates that coordinate systems belong
            # to fields in A and that unsupported combinations are rejected.
            Q = _build_multifield_coordinate_transformation(A, coordSys)
            for P in problems
                if _has_approximation(P) &&
                   any(cs -> cs.e1.model === P, coordSys)
                    error("solveEigen: coordSys + reduced/spectral basis on the same field is not yet implemented.")
                end
            end
        end
    end

    # Constraints are imposed on LOCAL field components, *after* Q.
    T = nothing
    R = nothing
    free = Int[]
    fixed_r = Int[]
    prunable = Int[]
    protected = Int[]

    if A.problems === nothing
        P = A.model
        T, R = _singlefield_mpc_transformation(P, mpc)
        fixed = constrainedDoFs(P, support)
        fixed_eff, uD_eff = _singlefield_mpc_bc_data(
            P, mpc, fixed, zeros(Float64, ndofs(P))
        )
        free, fixed_r, _ = reduced_bc_data(R, fixed_eff, uD_eff)
    else
        T, R, prunable, protected = _multifield_mpc_transformation(A, mpc)
        _, fixed, _ = multifield_bc_data(A, support; nsteps=1)
        fixed_eff, xD_eff = _multifield_mpc_bc_data(
            A, mpc, fixed, zeros(Float64, nfull)
        )
        free, fixed_r, _ = reduced_bc_data(R, fixed_eff, xD_eff)
    end

    Acoord = Q === nothing ? A.A : Q' * A.A * Q
    Bcoord = Q === nothing ? B.A : Q' * B.A * Q

    # Avoid allocating two large nested sparse products when possible.
    Ar = T' * (Acoord * T)
    Br = T' * (Bcoord * T)

    if A.problems !== nothing && !isempty(mpc)
        free = _remove_inactive_mpc_dofs(
            (Ar, Br), nothing, free, fixed_r, protected, prunable
        )
    end

    isempty(free) && error("solveEigen: no free degrees of freedom after reduction.")

    return (
        A=Ar[free, free], B=Br[free, free], T=T,
        Q=Q, free=free, nfull=nfull,
        symmetricA=symA, symmetricB=symB,
        uploA=symA ? A0.uplo : :U,
        uploB=symB ? B0.uplo : :U
    )
end

"""
    solve_eigen_reduced(prepared; n=6, shift=0.0, which=:LM,
                        solver=:auto, solveroptions=(;))

Solve the prepared generalized eigenproblem. `:direct` returns the selected
`n` eigenpairs of the dense generalized problem; `:arpack` uses sparse
shift-invert or ordinary Arnoldi/Lanczos iterations.

When both inputs were marked `Symmetric`, symmetric algorithms are attempted
only if the reduced mass/metric matrix is positive definite. Indefinite and
nonsymmetric pairs use the general solver instead.
"""
function solve_eigen_reduced(
    prepared;
    n::Int=6,
    shift=0.0,
    which::Symbol=:LM,
    solver::Symbol=:auto,
    solveroptions=(;)
    )
    A, B = prepared.A, prepared.B
    N = size(A, 1)
    1 <= n <= N || error("solveEigen: require 1 <= n <= $N, got n=$n.")
    solver in (:auto, :direct, :arpack) ||
        error("solveEigen: solver must be :auto, :direct, or :arpack.")
    solveroptions isa NamedTuple ||
        error("solveEigen: solveroptions must be a NamedTuple.")

    # ARPACK cannot reliably return an almost-complete eigenspectrum;
    # choose a dense solve automatically for small/reduced problems.
    chosen = solver === :auto ? (N <= 32 || n >= N - 1 ? :direct : :arpack) : solver

    if chosen === :arpack && n >= N - 1
        error("solveEigen: ARPACK requires fewer eigenpairs than N-1; use solver=:direct.")
    end

    # A symmetric pencil requires a positive-definite B to use Hermitian
    # generalized solvers. A symmetric but indefinite B (e.g. buckling)
    # must NOT be sent blindly to Lanczos/Cholesky.
    symmetric = prepared.symmetricA && prepared.symmetricB &&
                (prepared.uploA == prepared.uploB)
    if symmetric
        try
            cholesky(Symmetric(B, prepared.uploB))
        catch
            symmetric = false
        end
    end

    AA = symmetric ? Symmetric(A, prepared.uploA) : A
    BB = symmetric ? Symmetric(B, prepared.uploB) : B

    if chosen === :direct
        # One dense solve obtains the full spectrum. Keep complex eigenvalues
        # intact and select the requested part using the same target logic.
        vals, vecs = if symmetric
            eigen(Symmetric(Matrix(A), prepared.uploA),
                  Symmetric(Matrix(B), prepared.uploB))
        else
            eigen(Matrix(A), Matrix(B))
        end
        ix = _eigen_selected_indices(vals, min(n, length(vals)), shift, which)
        return vals[ix], vecs[:, ix]
    end

    options = merge((
        ; nev=n, 
        which=which, 
        maxiter=10000, 
        explicittransform=:none), 
        solveroptions)
    if shift !== nothing
        options = merge(options, (; sigma=shift))
    end
    result = Arpack.eigs(AA, BB; options...)
    λ, ϕ, nconv = result[1], result[2], result[3]
    nconv >= n || error("solveEigen: ARPACK converged only $nconv of $n eigenpairs.")
    return λ, ϕ
end

function _eigen_selected_indices(λ, n::Int, shift, which::Symbol)
    if shift !== nothing
        # Shift-invert: :LM means eigenvalues closest to the supplied shift.
        distances = abs.(λ .- shift)
        ordering = which in (:SM, :SR, :SI) ?
            sortperm(distances; rev=true) : sortperm(distances)
    elseif which === :SM
        ordering = sortperm(abs.(λ))
    elseif which === :LM
        ordering = sortperm(abs.(λ); rev=true)
    elseif which === :SR
        ordering = sortperm(real.(λ))
    elseif which === :LR
        ordering = sortperm(real.(λ); rev=true)
    elseif which === :SI
        ordering = sortperm(imag.(λ))
    elseif which === :LI
        ordering = sortperm(imag.(λ); rev=true)
    else
        error("solveEigen: unknown which=$which.")
    end
    return ordering[1:n]
end

"""
    reconstruct_eigenvectors(ϕr, prepared; normalize=:unit)

Reconstruct eigenvectors in the original full, Cartesian DOF space. Zero
Dirichlet perturbations are inserted before restoring approximation/MPC and
local-coordinate transformations. `:mass` means normalization with respect
to the supplied `B` matrix, and requires a positive real modal metric.
"""
function reconstruct_eigenvectors(ϕr, prepared; normalize::Symbol=:unit, Bfull=nothing)
    normalize in (:none, :unit, :mass, :B) ||
        error("solveEigen: normalize must be :none, :unit, :mass, or :B.")

    ϕ = zeros(eltype(ϕr), prepared.nfull, size(ϕr, 2))
    if prepared.T === nothing
        ϕ[prepared.free, :] .= ϕr
    else
        X = zeros(eltype(ϕr), size(prepared.T, 2), size(ϕr, 2))
        X[prepared.free, :] .= ϕr
        ϕ .= prepared.T * X
    end
    if prepared.Q !== nothing
        ϕ = prepared.Q * ϕ
    end

    for j in axes(ϕ, 2)
        v = @view ϕ[:, j]
        if normalize === :unit
            d = norm(v)
            d > 0 || error("solveEigen: eigenvector $j has zero norm.")
            v ./= d
        elseif normalize in (:mass, :B)
            Bfull === nothing && error("solveEigen: normalize=:mass requires B.")
            z = dot(v, Bfull * v)
            (real(z) > 0 && abs(imag(z)) <= 1e-8 * max(1.0, abs(real(z)))) ||
                error("solveEigen: B-normalization requires positive real phi'*B*phi, got $z.")
            v ./= sqrt(real(z))
        end
    end
    return ϕ
end

"""
    solveEigen(A, B; n=6, shift=0.0, which=:LM,
               support=BoundaryCondition[], mpc=MPC[], coordSys=[],
               solver=:auto, solveroptions=(;), normalize=:unit)

Solve `A*phi = λ*B*phi` for single-field or multifield `SystemMatrix`
arguments. The result is `(λ, ϕ)` with concrete numerical element types
(`Vector{Float64}` / `Vector{ComplexF64}` and matching matrices as needed).
Neither negative eigenvalues nor imaginary parts are discarded.

The `coordSys` keyword is identical to that of `solveField`, relying on its
coordinate builders. Spectral basis + MPC remains not implemented.
"""
function solveEigen(
    A::Union{SystemMatrix,SymmetricSystemMatrix},
    B::Union{SystemMatrix,SymmetricSystemMatrix};
    n::Int=6,
    shift=0.0,
    which::Symbol=:LM,
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[],
    coordSys=[],
    solver::Symbol=:auto,
    solveroptions=(;),
    normalize::Symbol=:unit
)
    prepared = prepare_eigen_system(A, B;
        support=support, mpc=mpc, coordSys=coordSys)
    λ, ϕr = solve_eigen_reduced(prepared;
        n=n, shift=shift, which=which, solver=solver,
        solveroptions=solveroptions)
    Bfull = normalize in (:mass, :B) ? _eigen_parent(B).A : nothing
    ϕ = reconstruct_eigenvectors(ϕr, prepared;
        normalize=normalize, Bfull=Bfull)
    return λ, ϕ
end

"""
    solveEigenProblem(A, B; directSolver=false, kwargs...)

Backward-compatible alias of [`solveEigen`](@ref). Unlike the historical
implementation the new routine respects `n` also with the direct solver.
"""
function solveEigenProblem(
    A::Union{SystemMatrix,SymmetricSystemMatrix},
    B::Union{SystemMatrix,SymmetricSystemMatrix};
    directSolver::Bool=false,
    kwargs...
)
    return solveEigen(A, B;
        solver=directSolver ? :direct : :arpack,
        kwargs...)
end

"""
    reduced_system_matrices(A, B, support[, mpc])

Compatibility helper returning `(Ar, Br, T, free)` for legacy callers.
Unlike `prepare_eigen_system`, this helper does not accept `coordSys`.
"""
function reduced_system_matrices(
    A::SystemMatrix, B::SystemMatrix,
    support::Vector{BoundaryCondition}, mpc::Vector{MPC}=MPC[]
)
    p = prepare_eigen_system(A, B; support=support, mpc=mpc)
    return p.A, p.B, p.T, p.free
end

# -----------------------------------------------------------------------------
# Legacy mechanical adapters. These are deliberately separate from solveEigen.
# -----------------------------------------------------------------------------

function _eigen_real_modes(λ, ϕ; context="modal analysis")
    # A mode's negligible imaginary rounding noise is not physical, but a
    # genuinely complex mode must never be silently flattened for Gmsh.
    imλ = isempty(λ) ? 0.0 : maximum(abs, imag.(λ))
    imϕ = isempty(ϕ) ? 0.0 : maximum(abs, imag.(ϕ))
    scλ = max(1.0, isempty(λ) ? 0.0 : maximum(abs, λ))
    scϕ = max(1.0, isempty(ϕ) ? 0.0 : maximum(abs, ϕ))
    (imλ <= 1e-8*scλ && imϕ <= 1e-8*scϕ) ||
        error("$context: complex eigenpairs cannot be displayed as real mode shapes. Use solveEigen to obtain the full complex result.")
    return real.(λ), real.(ϕ)
end

"""
    solveEigenModes(K, M; n=6, fₘᵢₙ=0.01, support=[], mpc=[],
                    coordSys=[], solver=:auto, directSolver=false)

Compute natural frequencies and real mode shapes for a single mechanical
field and return the existing Gmsh-compatible `Eigen` container.
"""
function solveEigenModes(
    K::SystemMatrix, M::SystemMatrix;
    n::Int=6, fₘᵢₙ=0.01,
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[], coordSys=[],
    solver::Symbol=:auto, solveroptions=(;),
    directSolver::Bool=false,
    normalize::Symbol=:mass
)
    K.problems === nothing ||
        error("solveEigenModes: use solveEigenFields for multifield systems.")
    λ, ϕ = solveEigen(K, M;
        n=n, shift=(2π*fₘᵢₙ)^2, which=:LM,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=directSolver ? :direct : solver,
        solveroptions=solveroptions, normalize=normalize)
    ω2, modes = _eigen_real_modes(λ, ϕ; context="solveEigenModes")
    minimum(ω2) >= -1e-9 * max(1.0, maximum(abs, ω2)) ||
        error("solveEigenModes: negative eigenvalues detected; use solveEigen to inspect instability.")
    freq = sqrt.(max.(ω2, 0.0)) ./ (2π)
    return Eigen(freq, modes, K.model, nothing, nothing)
end

"""
    solveEigenFields(K, M; n=6, fmin=0.0, support=[], mpc=[], coordSys=[])

Return a tuple of `Eigen` containers, one per physical field, for real
multifield vibration modes. Use `solveEigen` directly for complex problems.
"""
function solveEigenFields(
    K::SystemMatrix, M::SystemMatrix;
    n::Int=6, fmin=0.0,
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[], coordSys=[],
    solver::Symbol=:auto, solveroptions=(;),
    directSolver::Bool=false,
    normalize::Symbol=:mass
)
    K.problems !== nothing ||
        error("solveEigenFields: use solveEigenModes for single-field systems.")
    λ, ϕ = solveEigen(K, M;
        n=n, shift=(2π*fmin)^2, which=:LM,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=directSolver ? :direct : solver,
        solveroptions=solveroptions, normalize=normalize)
    ω2, modes = _eigen_real_modes(λ, ϕ; context="solveEigenFields")
    minimum(ω2) >= -1e-9 * max(1.0, maximum(abs, ω2)) ||
        error("solveEigenFields: negative eigenvalues detected; use solveEigen.")
    f = sqrt.(max.(ω2, 0.0)) ./ (2π)
    return splitEigenToEigen(Eigen(f, modes, nothing, K.problems, K.offsets))
end

"""
    solveBucklingModes(K, Kg; n=6, support=[], mpc=[], coordSys=[])

Compute critical load multipliers from `K*phi = -λ*Kg*phi` and return a
Gmsh-compatible `Eigen` container for a real, single-field problem.
"""
function solveBucklingModes(
    K::SystemMatrix, Kg::SystemMatrix;
    n::Int=6,
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[], coordSys=[],
    solver::Symbol=:auto, solveroptions=(;),
    directSolver::Bool=false
)
    K.problems === nothing ||
        error("solveBucklingModes: single-field system expected.")
    Kg.model === K.model ||
        error("solveBucklingModes: K and Kg must refer to the same Problem.")
    B = SystemMatrix(-Kg.A, Kg.model)
    λ, ϕ = solveEigen(K, B;
        n=n, shift=0.0, which=:LM,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=directSolver ? :direct : solver,
        solveroptions=solveroptions, normalize=:unit)
    vals, modes = _eigen_real_modes(λ, ϕ; context="solveBucklingModes")
    return Eigen(vals, modes, K.model, nothing, nothing)
end

"""
    splitEigenToEigen(eig::Eigen)

Split a global multifield `Eigen` container into one `Eigen` object per
Problem while preserving the common spectrum and physical DOF ordering.
"""
function splitEigenToEigen(eig::Eigen)
    eig.problems === nothing &&
        error("splitEigenToEigen: a multifield Eigen container is required.")
    results = Vector{Eigen}(undef, length(eig.problems))
    for (i, P) in enumerate(eig.problems)
        off = eig.offsets[i]
        nr = ndofs(P)
        results[i] = Eigen(
            eig.f, copy(eig.ϕ[off+1:off+nr, :]), P, nothing, nothing
        )
    end
    return tuple(results...)
end

"""
    smallestEigenValue(K, C; support=[], mpc=[], coordSys=[])

Return the generalized eigenvalue closest to zero.
"""
function smallestEigenValue(
    K::SystemMatrix, C::SystemMatrix;
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[], coordSys=[], solver::Symbol=:auto
)
    λ, _ = solveEigen(K, C; n=1, shift=1e-8, which=:LM,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=solver, normalize=:none)
    return only(λ)
end

"""
    largestEigenValue(K, C; support=[], mpc=[], coordSys=[])

Return the generalized eigenvalue of largest magnitude.
"""
function largestEigenValue(
    K::SystemMatrix, C::SystemMatrix;
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[], coordSys=[], solver::Symbol=:auto
)
    λ, _ = solveEigen(K, C; n=1, shift=nothing, which=:LM,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=solver, normalize=:none)
    return only(λ)
end

function _eigen_period(λ)
    abs(imag(λ)) <= 1e-8 * max(1.0, abs(real(λ))) ||
        error("Eigenperiod: a complex eigenvalue has no real undamped period.")
    real(λ) > 0 ||
        error("Eigenperiod: the eigenvalue must be positive, got $λ.")
    return 2π / sqrt(real(λ))
end

"""Return the longest undamped vibration period for positive eigenvalues."""
function largestPeriodTime(K::SystemMatrix, M::SystemMatrix;
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[], coordSys=[])
    return _eigen_period(smallestEigenValue(K, M;
        support=support, mpc=mpc, coordSys=coordSys))
end

"""Return the shortest undamped vibration period for positive eigenvalues."""
function smallestPeriodTime(K::SystemMatrix, M::SystemMatrix;
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[], coordSys=[])
    return _eigen_period(largestEigenValue(K, M;
        support=support, mpc=mpc, coordSys=coordSys))
end

# -----------------------------------------------------------------------------
# Prestress convenience workflows: geometric tangent is calculated separately.
# The present PR only changes how the final eigenproblem is solved. These
# routines remain available, but a general load-enabled solveEigen is deferred.
# -----------------------------------------------------------------------------

"""
    solveModalAnalysis(problem; support=[], load=[], n=6,
                       fₘᵢₙ=1e-5, directSolver=false)

Compute real undamped modes of a structure, optionally around the current
iteratively updated prestressed equilibrium. The eigensolver itself receives
the full matrices and applies the support constraints exactly once.
"""
function solveModalAnalysis(
    problem::Problem;
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    load::Vector{BoundaryCondition}=BoundaryCondition[],
    n::Int=6, fₘᵢₙ=1e-5, directSolver::Bool=false
)
    K = stiffnessMatrix(problem)
    M = massMatrix(problem)
    if !isempty(load)
        f = loadVector(problem, load)
        q = solveDisplacement(K, f; support=support)
        err = Inf
        iteration = 0
        while err > 1e-3 && iteration < 10
            iteration += 1
            q0 = copy(q)
            Kg = nonLinearStiffnessMatrix(q)
            q = solveDisplacement(K + Kg, f; support=support)
            denom = max(sum(abs, q0.a), eps(Float64))
            err = sum(abs, q.a - q0.a) / denom
        end
        err > 1e-3 && @warn "solveModalAnalysis: prestress iteration did not converge (relative change $err)."
        K = K + nonLinearStiffnessMatrix(q)
    end
    return solveEigenModes(K, M;
        n=n, fₘᵢₙ=fₘᵢₙ, support=support, directSolver=directSolver)
end

"""
    solveBuckling(problem; load=[], support=[], n=6)

Compute critical load multipliers after calculating the geometric stiffness
of the prestressed reference configuration. Contact is not included.
"""
function solveBuckling(
    problem::Problem;
    load::Vector{BoundaryCondition}=BoundaryCondition[],
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    n::Int=6
)
    K = stiffnessMatrix(problem)
    f = loadVector(problem, load)
    q = solveDisplacement(K, f; support=support)
    err = Inf
    iteration = 0
    while err > 1e-3 && iteration < 10
        iteration += 1
        q0 = copy(q)
        Kg = nonLinearStiffnessMatrix(q)
        q = solveDisplacement(K + Kg, f; support=support)
        denom = max(sum(abs, q0.a), eps(Float64))
        err = sum(abs, q.a - q0.a) / denom
    end
    err > 1e-3 && @warn "solveBuckling: prestress iteration did not converge (relative change $err)."
    Kg = nonLinearStiffnessMatrix(q)
    return solveBucklingModes(K, Kg; n=n, support=support)
end
