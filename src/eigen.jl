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

export solveEigen, solveEigenProblem, eigenVector, eigenValue, showEigenResults
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

# ARPACK keyword names supported by its generalized eigs(A, B; ...) API.
const _EIGEN_ARPACK_OPTIONS = (
    :nev, :ncv, :which, :tol, :maxiter, :sigma,
    :ritzvec, :explicittransform, :v0, :check
)
const _EIGEN_WHICH_OPTIONS = (:LM, :SM, :LR, :SR, :LI, :SI, :BE)

function _eigen_check_choice(name::Symbol, value, choices)
    value in choices || throw(ArgumentError(
        "solveEigen: invalid $name=$(repr(value)). " *
        "Available values: $(join(string.(choices), ", "))."
    ))
    return value
end

"""
    _eigen_selected_indices(λ, n; sigma=nothing, which=:SR, shift=nothing)

Select from an already computed spectrum. In native ARPACK mode, `which`
refers to the transformed eigenvalues `1/(λ-sigma)` when `sigma` is set.
For simple threshold mode, select the smallest `λ > shift` by real part.
The returned indices are always sorted by `(real(λ), imag(λ))`.
"""
function _eigen_selected_indices(λ, n::Int;
    sigma=nothing, which::Symbol=:SR, shift=nothing)

    candidates = [i for i in eachindex(λ)
                  if isfinite(real(λ[i])) && isfinite(imag(λ[i])) &&
                     (shift === nothing || real(λ[i]) > shift)]

    if shift !== nothing || (sigma === nothing && which === :SR)
        order = sort(candidates; by=i -> (real(λ[i]), imag(λ[i])))
    else
        # ARPACK's `which` targets the *transformed* spectrum when sigma is set.
        transformed = sigma === nothing ? λ : 1 ./ (λ .- sigma)
        if which === :BE
            small = sort(candidates; by=i -> real(transformed[i]))
            nlo = fld(n, 2)
            nhi = n - nlo
            order = vcat(small[1:min(nlo, length(small))],
                         reverse(small[max(1, length(small)-nhi+1):end]))
        else
            score(i) = which === :LM ? abs(transformed[i]) :
                       which === :SM ? -abs(transformed[i]) :
                       which === :LR ? real(transformed[i]) :
                       which === :SR ? -real(transformed[i]) :
                       which === :LI ? imag(transformed[i]) :
                       -imag(transformed[i])  # :SI
            order = sort(candidates; by=score, rev=true)
        end
    end

    length(order) >= n || throw(ArgumentError(
        "solveEigen: only $(length(order)) eigenvalues satisfy the requested " *
        "selection; $n requested." *
        (shift === nothing ? "" : " Try a smaller shift or request fewer modes.")
    ))
    selected = order[1:n]
    return sort(selected; by=i -> (real(λ[i]), imag(λ[i])))
end

"""
    solve_eigen_reduced(prepared; n=6, shift=nothing, sigma=nothing,
                        which=nothing, solver=:auto, solveroptions=(;))

Solve a generalized eigenproblem and return exactly `n` eigenpairs (sorted
by increasing real part, then imaginary part).

Two mutually exclusive search modes are supported:

- Simple mode: `shift=lower_bound` selects the first `n` real eigenvalues
  strictly greater than `lower_bound`. For sparse iterative calculation, this
  is guaranteed only for real symmetric matrices with positive-definite B.
  ARPACK uses one-sided shift-invert (`sigma≈shift`, `which=:LR`).
- Native ARPACK mode: `sigma` and `which` retain their original ARPACK
  meanings. Without `sigma` or `shift`, the default is `which=:SR`, selecting
  eigenvalues with smallest real part.

Native ARPACK options may be passed via `solveroptions`, and are checked
before dispatch. ARPACK's `ritzvec=false` is incompatible with this API.
"""
function _eigen_nearly_symmetric(M)
    eltype(M) <: Real || return false
    scale = opnorm(M, Inf)
    scale == 0 && return true
    return opnorm(M - M', Inf) <= 1e-12 * scale
end

function solve_eigen_reduced(
    prepared;
    n::Int=6,
    shift=nothing,
    sigma=nothing,
    which=nothing,
    solver::Symbol=:auto,
    solveroptions=(;)
)
    A, B = prepared.A, prepared.B
    N = size(A, 1)
    1 <= n <= N || throw(ArgumentError(
        "solveEigen: n must be between 1 and $N (received $n)."
    ))
    _eigen_check_choice(:solver, solver, (:auto, :direct, :arpack))
    solveroptions isa NamedTuple || throw(ArgumentError(
        "solveEigen: solveroptions must be a NamedTuple, e.g. (ncv=40, tol=1e-9)."
    ))
    for key in keys(solveroptions)
        key in _EIGEN_ARPACK_OPTIONS || throw(ArgumentError(
            "solveEigen: unknown ARPACK option :$key. Available options: " *
            join(string.(_EIGEN_ARPACK_OPTIONS), ", ") * "."
        ))
    end
    if :nev in keys(solveroptions)
        solveroptions.nev == n || throw(ArgumentError(
            "solveEigen: conflicting n=$n and solveroptions.nev=$(solveroptions.nev). " *
            "Use n or nev to set the number of eigenpairs."
        ))
    end
    if :sigma in keys(solveroptions)
        (sigma === nothing || isequal(sigma, solveroptions.sigma)) ||
            throw(ArgumentError("solveEigen: conflicting sigma and solveroptions.sigma."))
        sigma = solveroptions.sigma
    end
    if :which in keys(solveroptions)
        (which === nothing || isequal(which, solveroptions.which)) ||
            throw(ArgumentError("solveEigen: conflicting which and solveroptions.which."))
        which = solveroptions.which
    end
    shift === nothing || (shift isa Real && isfinite(shift)) ||
        throw(ArgumentError("solveEigen: shift must be a finite real lower bound or nothing."))
    (sigma === nothing || sigma isa Number) ||
        throw(ArgumentError("solveEigen: sigma must be a number or nothing (ARPACK option)."))
    !(shift !== nothing && sigma !== nothing) || throw(ArgumentError(
        "solveEigen: shift and sigma are different search modes; use only one. " *
        "shift is a lower bound; sigma is the original ARPACK target."
    ))
    !(shift !== nothing && which !== nothing) || throw(ArgumentError(
        "solveEigen: which is a native ARPACK option and cannot accompany shift. " *
        "Use shift alone, or sigma together with which."
    ))
    if which !== nothing
        _eigen_check_choice(:which, which, _EIGEN_WHICH_OPTIONS)
    end
    selection = shift !== nothing ? :threshold : :native
    arpack_which = selection === :threshold ? :LR :
                   which === nothing ? (sigma === nothing ? :SR : :LM) : which

    for (key, choices) in ((:explicittransform, (:auto, :none, :shiftinvert)),
                           (:check, (0, 1, 2)))
        if key in keys(solveroptions)
            _eigen_check_choice(key, getproperty(solveroptions, key), choices)
        end
    end
    if (:ritzvec in keys(solveroptions)) && solveroptions.ritzvec !== true
        throw(ArgumentError(
            "solveEigen: ritzvec must be true because solveEigen returns eigenvectors ϕ."
        ))
    end
    if selection === :threshold && :explicittransform in keys(solveroptions) &&
       solveroptions.explicittransform !== :none
        throw(ArgumentError(
            "solveEigen: shift requires explicittransform=:none; " *
            "valid values in native ARPACK mode are :none, :auto, :shiftinvert."
        ))
    end

    symmetric = eltype(A) <: Real && eltype(B) <: Real &&
                (prepared.symmetricA || _eigen_nearly_symmetric(A)) &&
                (prepared.symmetricB || _eigen_nearly_symmetric(B))
    if symmetric
        try
            cholesky(Symmetric(B, prepared.uploB))
        catch
            symmetric = false
        end
    end
    if arpack_which === :BE && !symmetric
        throw(ArgumentError(
            "solveEigen: which=:BE requires real symmetric A and positive-definite B. " *
            "Available options: $(join(string.(_EIGEN_WHICH_OPTIONS), ", "))."
        ))
    end
    if symmetric && arpack_which in (:LI, :SI)
        throw(ArgumentError(
            "solveEigen: which=$arpack_which cannot be used for real symmetric eigenproblems. " *
            "Choose :LM, :SM, :LR, :SR or :BE."
        ))
    end

    chosen = solver === :auto ?
        (N <= 32 || n >= (symmetric ? N : N-1) ? :direct : :arpack) : solver
    if selection === :threshold && chosen === :arpack && !symmetric
        throw(ArgumentError(
            "solveEigen: iterative shift selection requires real symmetric A and B " *
            "with positive-definite B. Use Symmetric(A), Symmetric(B), " *
            "solver=:direct, or native ARPACK sigma/which for nonsymmetric systems."
        ))
    end
    if chosen === :arpack && n > (symmetric ? N-1 : N-2)
        throw(ArgumentError(
            "solveEigen: ARPACK supports at most $(symmetric ? N-1 : max(0,N-2)) " *
            "eigenpairs at this size; use solver=:direct."
        ))
    end

    AA = symmetric ? Symmetric(A, prepared.uploA) : A
    BB = symmetric ? Symmetric(B, prepared.uploB) : B

    if chosen === :direct
        nontrivial = [key for key in keys(solveroptions)
                      if key ∉ (:sigma, :which, :nev)]
        isempty(nontrivial) || throw(ArgumentError(
            "solveEigen: ARPACK-only options $(join(string.(nontrivial), ", ")) " *
            "cannot be used with solver=:direct. Supported direct options: " *
            "n/nev, shift OR sigma, which, normalize."
        ))
        vals, vecs = if symmetric
            eigen(Symmetric(Matrix(A), prepared.uploA),
                  Symmetric(Matrix(B), prepared.uploB))
        else
            eigen(Matrix(A), Matrix(B))
        end
        ix = _eigen_selected_indices(vals, n;
            shift=shift, sigma=sigma, which=arpack_which)
        return vals[ix], vecs[:, ix]
    end

    # ARPACK requires ritzvec=true: eigenvectors must be returned.
    opts = merge((;
        nev=n, which=arpack_which, maxiter=10000,
        explicittransform=:none, ritzvec=true
    ), solveroptions)
    if selection === :threshold
        # A small offset above the strict bound avoids singular A-shift*B
        # and places eigenvalues exactly at the bound on the *negative*
        # side of the shift-invert spectrum. Near-zero numerical noise in
        # rigid-body modes is then rejected by the one-sided ARPACK search.
        size_scale = opnorm(A, Inf) / max(opnorm(B, Inf), floatmin(Float64))
        gap = 64 * eps(Float64) * max(abs(Float64(shift)), size_scale)
        target = Float64(shift) + max(gap, eps(Float64)*max(1.0, abs(Float64(shift))))
        opts = merge(opts, (; which=:LR, sigma=target, explicittransform=:none))
    else
        opts = merge(opts, (; sigma=sigma, which=arpack_which))
    end

    result = Arpack.eigs(AA, BB; opts...)
    λ, ϕ, nconv = result[1], result[2], result[3]
    nconv >= n || error("solveEigen: ARPACK converged only $nconv of $n eigenpairs.")
    ix = _eigen_selected_indices(λ, n;
        shift=shift, sigma=sigma, which=arpack_which)
    return λ[ix], ϕ[:, ix]
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
    solveEigen(A, B; n=6, nev=nothing, shift=nothing, sigma=nothing,
               which=nothing, support=[], mpc=[], coordSys=[],
               solver=:auto, solveroptions=(;), normalize=:unit,
               output=:field, part=:real, kwargs...)

Solve `A*ϕ = λ*B*ϕ` for single-field or multifield problems.

By default return one standard nodal Field per Problem (as `solveField` does),
with columns representing mode shapes and `field.t` holding the real or imaginary
parts of the eigenvalues. `part=:real`, `:imag`, or `:both` selects components;
`:both` returns the two Fields consecutively for each Problem.

With `output=:raw`, return `(λ, ϕ)` with the complete real/complex generalized
eigenpairs. `part` is ignored in this mode. Eigenpairs are ordered by increasing
real part of λ (then imaginary part).

# Selection modes

- `shift=lo`: Return the first `n` eigenpairs with `real(λ) > lo`. For ARPACK,
  this one-sided selection is currently implemented only for real symmetric
  pencils with positive-definite B. Use `solver=:direct` for other pencils.
- `sigma=s, which=:LM`: Native ARPACK shift-invert selection; modes on both
  sides of `s` may be included. `which` and `sigma` retain ARPACK semantics.
- Neither `shift` nor `sigma`: `which=:SR` by default (smallest real part).
  For large positive elastic systems, `shift=0.0` is usually faster.

Use `n` or its alias `nev` for the number of eigenpairs. The original ARPACK
keywords `ncv`, `tol`, `maxiter`, `v0`, `check`, `explicittransform`, `which`
and `sigma` can be passed directly, or in `solveroptions=(; ...)`.
`ritzvec=false` is unsupported because eigenvectors are always returned.

Coordinates, boundary conditions and MPCs share the `solveField` reduction.
"""
function solveEigen(
    A::Union{SystemMatrix,SymmetricSystemMatrix},
    B::Union{SystemMatrix,SymmetricSystemMatrix};
    n::Union{Nothing,Int}=nothing,
    nev::Union{Nothing,Int}=nothing,
    shift=nothing,
    sigma=nothing,
    which=nothing,
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[],
    coordSys=[],
    solver::Symbol=:auto,
    solveroptions=(;),
    normalize::Symbol=:unit,
    output::Symbol=:field,
    part::Symbol=:real,
    kwargs...
)
    _eigen_check_choice(:normalize, normalize, (:none, :unit, :mass, :B))
    _eigen_check_choice(:output, output, (:field, :raw))
    if output === :field
        _eigen_check_choice(:part, part, (:real, :imag, :both))
    end
    if n !== nothing && nev !== nothing
        throw(ArgumentError("solveEigen: use n or nev, not both (n=$n, nev=$nev)."))
    end
    solveroptions isa NamedTuple || throw(ArgumentError(
        "solveEigen: solveroptions must be a NamedTuple, e.g. (ncv=40, tol=1e-9)."
    ))
    native = (; kwargs...)
    for key in keys(native)
        key in _EIGEN_ARPACK_OPTIONS || throw(ArgumentError(
            "solveEigen: unknown keyword :$key. Available ARPACK keywords: " *
            join(string.(_EIGEN_ARPACK_OPTIONS), ", ") *
            "; LowLevelFEM options: n, nev, shift, sigma, which, support, " *
            "mpc, coordSys, solver, normalize, output, part, solveroptions."
        ))
    end
    for key in keys(native)
        if key in keys(solveroptions) && !isequal(native[key], solveroptions[key])
            throw(ArgumentError("solveEigen: conflicting :$key in keywords and solveroptions."))
        end
    end
    native = merge(solveroptions, native)
    amount = n === nothing ? (nev === nothing ? get(native, :nev, 6) : nev) : n
    amount isa Int || throw(ArgumentError("solveEigen: n/nev must be an integer."))
    if :nev in keys(native) && amount != native.nev
        throw(ArgumentError("solveEigen: conflicting n/nev and solveroptions.nev."))
    end
    prepared = prepare_eigen_system(A, B;
        support=support, mpc=mpc, coordSys=coordSys)
    λ, ϕr = solve_eigen_reduced(prepared;
        n=amount, shift=shift, sigma=sigma, which=which,
        solver=solver, solveroptions=native)
    Bfull = normalize in (:mass, :B) ? _eigen_parent(B).A : nothing
    ϕ = reconstruct_eigenvectors(ϕr, prepared;
        normalize=normalize, Bfull=Bfull)
    output === :raw && return λ, ϕ
    return _eigen_fields(λ, ϕ, _eigen_parent(A); part=part)
end

# =============================================================================
# Field-valued eigenresults (without changing the algebraic eigenvalue solver)
# =============================================================================

const _EigenFieldType = Union{ScalarField,VectorField,TensorField}
const _EIGEN_FIELD_KINDS = (:eigen_real, :eigen_imag)

"""
    _eigen_field(P, modes, eigenvalues, kind)

Wrap a real-valued block of nodal eigenvectors in the standard LowLevelFEM
Field type for `P`. Each column is one mode and each entry of `t` is the
corresponding real-valued eigenvalue component (not a physical time).
"""
function _eigen_field(P::Problem, modes::AbstractMatrix,
                      eigenvalues::AbstractVector, kind::Symbol)
    kind in _EIGEN_FIELD_KINDS || throw(ArgumentError(
        "eigen field: kind must be :eigen_real or :eigen_imag."))
    size(modes, 1) == ndofs(P) || throw(DimensionMismatch(
        "eigen field: $(size(modes,1)) mode DOFs, expected $(ndofs(P)) for field $(P.field)."))
    size(modes, 2) == length(eigenvalues) || throw(DimensionMismatch(
        "eigen field: eigenvector column count must equal eigenvalue count."))

    type = P.pdim == 1 ? :scalar :
           P.pdim == 2 ? :v2D :
           P.pdim == 3 ? :v3D :
           P.pdim == 9 ? :tensor :
           throw(ArgumentError("eigen field: unsupported component count $(P.pdim); expected 1, 2, 3 or 9."))

    T = P.pdim == 1 ? ScalarField :
        P.pdim in (2, 3) ? VectorField : TensorField

    return T(Matrix{Float64}[], Matrix{Float64}(modes),
             Float64.(eigenvalues), Int[], length(eigenvalues),
             type, P, kind)
end

"""
    _eigen_fields(λ, ϕ, K; part=:real)

Convert the complete eigenpairs of a single- or multifield SystemMatrix into
ordinary nodal Fields. `:real` and `:imag` select a component of both eigenvalues
and eigenvectors; `:both` returns paired real and imaginary Fields in problem
order: `(ur, ui, vr, vi, ...)`. Normalization must already have been applied to
the complete complex eigenvectors before this conversion.
"""
function _eigen_fields(λ::AbstractVector, ϕ::AbstractMatrix,
                       K::SystemMatrix; part::Symbol=:real)
    _eigen_check_choice(:part, part, (:real, :imag, :both))
    size(ϕ, 1) == size(K.A, 1) || throw(DimensionMismatch(
        "solveEigen: eigenvector DOFs and SystemMatrix dimension disagree."))
    size(ϕ, 2) == length(λ) || throw(DimensionMismatch(
        "solveEigen: eigenvalue and eigenvector counts disagree."))

    problems = K.problems === nothing ? [K.model] : K.problems
    offsets = K.offsets === nothing ? [0] : K.offsets
    length(problems) == length(offsets) ||
        error("solveEigen: multifield Problem/offset metadata are inconsistent.")

    results = _EigenFieldType[]
    for (i, P) in enumerate(problems)
        rows = (offsets[i] + 1):(offsets[i] + ndofs(P))
        if part in (:real, :both)
            push!(results, _eigen_field(P, real.(@view(ϕ[rows, :])),
                                        real.(λ), :eigen_real))
        end
        if part in (:imag, :both)
            push!(results, _eigen_field(P, imag.(@view(ϕ[rows, :])),
                                        imag.(λ), :eigen_imag))
        end
    end
    return length(results) == 1 ? only(results) : tuple(results...)
end

"""
    eigenValue(u::Union{ScalarField,VectorField,TensorField}; n=1)

Return the stored real or imaginary part of the `n`-th eigenvalue, as `Float64`.
Only Fields marked `:eigen_real` or `:eigen_imag` are accepted. For complex
eigenvalues, use `part=:both` to retain both parts separately.
"""
function eigenValue(u::_EigenFieldType; n::Integer=1)
    u.kind in _EIGEN_FIELD_KINDS || throw(ArgumentError(
        "eigenValue: this is an ordinary Field (kind=$(u.kind)), not an eigenmode Field."))
    1 <= n <= u.nsteps || throw(ArgumentError(
        "eigenValue: n must be between 1 and $(u.nsteps); got $n."))
    return u.t[n]
end

"""
    eigenVector(u::Union{ScalarField,VectorField,TensorField}; n=1)

Extract one mode as a single-step Field of the same type and representation.
The selected eigenvalue component remains in `t[1]`, and `kind` is preserved.
"""
function eigenVector(u::_EigenFieldType; n::Integer=1)
    u.kind in _EIGEN_FIELD_KINDS || throw(ArgumentError(
        "eigenVector: expected an eigenmode Field, got kind=$(u.kind)."))
    1 <= n <= u.nsteps || throw(ArgumentError(
        "eigenVector: n must be between 1 and $(u.nsteps); got $n."))
    A = isElementwise(u) ?
        Matrix{Float64}[copy(@view(Ae[:, n:n])) for Ae in u.A] :
        Matrix{Float64}[]
    a = isNodal(u) ? copy(@view(u.a[:, n:n])) : Matrix{Float64}(undef, 0, 0)
    return typeof(u)(A, a, [u.t[n]], copy(u.numElem),
                     1, u.type, u.model, u.kind)
end

"""
    showEigenResults(u; name=:eigen, visible=false, ff=2, factor=0)

Display the modes of an eigenmode Field in a Gmsh view. Each stored column
is a mode; its (real-valued) eigenvalue component is shown via the Field's
`t` data. `ff=1` and `ff=2` display mode shapes as displacements with eigenvalue labels;
`ff=0` selects an ordinary vector/field view with time labels.
The default `ff=2` uses eigenvalue labeling and displacement visualization. Returns the Gmsh view tag.
"""
function showEigenResults(u::_EigenFieldType;
                          name=:eigen, visible::Bool=false,
                          ff::Int=2, factor=0)
    u.kind in _EIGEN_FIELD_KINDS || throw(ArgumentError(
        "showEigenResults: expected an eigenmode Field; got kind=$(u.kind)."))
    ff in (0, 1, 2) || throw(ArgumentError(
        "showEigenResults: ff must be 0 (ordinary), 1 (modal), or 2 (buckling)."))
    comp = u isa ScalarField ? :scalar :
           u isa VectorField ? :vector : :tensor
    tag = showDoFResults(u, comp;
        name=string(name), visible=visible, ff=ff, factor=factor)
    # Show the eigenvalue label also when a single mode is displayed.
    if ff == 0
        gmsh.view.option.setNumber(tag, "ShowTime", 1)
    end
    return tag
end

"""
    showEigenResults(fields::Tuple; name=:eigen, visible=false, ff=2, factor=0)

Create one Gmsh view for each Field in a multifield eigenresult tuple.
Returns a tuple of view tags in the same order as the input Fields.
"""
function showEigenResults(fields::Tuple;
                          name=:eigen, visible::Bool=false,
                          ff::Int=2, factor=0)
    all(u -> u isa _EigenFieldType, fields) ||
        throw(ArgumentError("showEigenResults: tuple must contain only Field objects."))
    return tuple((showEigenResults(u;
        name="$(name) ($(u.model.field), $(u.kind), $i)",
        visible=visible, ff=ff, factor=factor)
        for (i, u) in enumerate(fields))...)
end

"""
    showModalResults(u::Field; name=:modal, visible=false, ff=1, factor=0)

Render a real eigenmode Field as an undamped vibration mode. The temporary
Gmsh view displays frequencies `sqrt(λ)/(2π)` rather than `λ`; the input Field
and its original eigenvalue labels are not modified.
"""
function showModalResults(u::_EigenFieldType;
                          name=:modal, visible::Bool=false,
                          ff::Int=1, factor=0)
    u.kind === :eigen_real || throw(ArgumentError(
        "showModalResults: use a :eigen_real Field containing real ω² eigenvalues."))
    scale = max(1.0, maximum(abs, u.t))
    minimum(u.t) >= -1e-10 * scale || throw(ArgumentError(
        "showModalResults: negative eigenvalues cannot be converted to real frequencies."))
    freq = sqrt.(max.(u.t, 0.0)) ./ (2π)
    viewfield = typeof(u)(u.A, u.a, freq, u.numElem, u.nsteps,
                          u.type, u.model, u.kind)
    return showEigenResults(viewfield;
        name=name, visible=visible, ff=ff, factor=factor)
end

"""
    showBucklingResults(u::Field; name=:buckling, visible=false, ff=2, factor=0)

Display buckling eigenvectors and their corresponding eigenvalue components.
"""
function showBucklingResults(u::_EigenFieldType;
                            name=:buckling, visible::Bool=false,
                            ff::Int=2, factor=0)
    return showEigenResults(u;
        name=name, visible=visible, ff=ff, factor=factor)
end

"""
    solveEigenProblem(A, B; shift=0.0, directSolver=false, kwargs...)

Backward-compatible alias. The historical `shift` means the ARPACK numerical
shift (`sigma`), *not* the `solveEigen` lower spectral bound.
"""
function solveEigenProblem(
    A::Union{SystemMatrix,SymmetricSystemMatrix},
    B::Union{SystemMatrix,SymmetricSystemMatrix};
    shift=0.0,
    directSolver::Bool=false,
    which::Symbol=:LM,
    kwargs...
)
    return solveEigen(A, B;
        sigma=shift, which=which,
        solver=directSolver ? :direct : :arpack,
        output=:raw,
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
        n=n, shift=(2π*fₘᵢₙ)^2,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=directSolver ? :direct : solver,
        solveroptions=solveroptions, normalize=normalize, output=:raw)
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
        n=n, shift=(2π*fmin)^2,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=directSolver ? :direct : solver,
        solveroptions=solveroptions, normalize=normalize, output=:raw)
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
        n=n, sigma=0.0, which=:LM,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=directSolver ? :direct : solver,
        solveroptions=solveroptions, normalize=:unit, output=:raw)
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

# -----------------------------------------------------------------------------
# Legacy scalar eigenvalue and vibration-period helpers.
# The full (possibly signed/complex) spectrum is available through
# solveEigen(...; output=:raw). These wrappers keep the historical Float64
# return value abs(real(λ)) used by older transient analysis notebooks.
# -----------------------------------------------------------------------------

"""
    smallestEigenValue(K, C; support=BoundaryCondition[], mpc=MPC[],
                       coordSys=[], solver=:auto, sigma=1e-8, which=:LR,
                       solveroptions=(;), kwargs...)

Return the smallest generalized eigenvalue near zero for `K*ϕ = λ*C*ϕ`.

For compatibility with the original implementation, the result is a `Float64`
containing `abs(real(λ))`. By default ARPACK uses a small numerical shift
(`sigma=1e-8`, `which=:LR`); if that solve fails numerically, it retries with
`which=:SM` and no numerical shift. Boundary conditions and MPCs use the same
reduction as `solveEigen`, for single- and multifield systems.

This legacy wrapper discards the sign and imaginary part of an eigenvalue.
Use `solveEigen(K, C; output=:raw)` to retain the original eigenvalue.
"""
function smallestEigenValue(
    K::Union{SystemMatrix,SymmetricSystemMatrix},
    C::Union{SystemMatrix,SymmetricSystemMatrix};
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[],
    coordSys=[],
    solver::Symbol=:auto,
    sigma::Number=1e-8,
    which::Symbol=:LR,
    solveroptions=(;),
    kwargs...
)
    # Use the raw result: solveEigen now returns Field objects by default.
    # Forward additional supported ARPACK keywords (e.g. tol, ncv, maxiter).
    value = try
        λ, _ = solveEigen(K, C;
            nev=1, sigma=sigma, which=which,
            support=support, mpc=mpc, coordSys=coordSys,
            solver=solver, solveroptions=solveroptions,
            normalize=:none, output=:raw, kwargs...)
        only(λ)
    catch err
        # Keep the old ARPACK fallback, but do not conceal invalid user options.
        if solver === :direct ||
           err isa ArgumentError || err isa DimensionMismatch ||
           err isa MethodError || which !== :LR || sigma != 1e-8
            rethrow()
        end
        @warn "smallestEigenValue: near-zero shift-invert failed; retrying with which=:SM and no sigma." exception=(err, catch_backtrace())
        λ, _ = solveEigen(K, C;
            nev=1, which=:SM,
            support=support, mpc=mpc, coordSys=coordSys,
            solver=solver, solveroptions=solveroptions,
            normalize=:none, output=:raw, kwargs...)
        only(λ)
    end
    return Float64(abs(real(value)))
end

"""
    largestEigenValue(K, C; support=BoundaryCondition[], mpc=MPC[],
                      coordSys=[], solver=:auto, which=:LM,
                      solveroptions=(;), kwargs...)

Return the generalized eigenvalue of largest magnitude for `K*ϕ = λ*C*ϕ`.

The default selection is the original unshifted ARPACK `which=:LM`. Return
`abs(real(λ))::Float64`, preserving the legacy API. This is intended for real
positive spectra such as undamped mechanical vibration problems; use the raw
`solveEigen` result for instability or complex spectral analysis.
"""
function largestEigenValue(
    K::Union{SystemMatrix,SymmetricSystemMatrix},
    C::Union{SystemMatrix,SymmetricSystemMatrix};
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    mpc::Vector{MPC}=MPC[],
    coordSys=[],
    solver::Symbol=:auto,
    which::Symbol=:LM,
    solveroptions=(;),
    kwargs...
)
    λ, _ = solveEigen(K, C;
        nev=1, which=which,
        support=support, mpc=mpc, coordSys=coordSys,
        solver=solver, solveroptions=solveroptions,
        normalize=:none, output=:raw, kwargs...)
    return Float64(abs(real(only(λ))))
end

"""
    largestPeriodTime(K, M; kwargs...)

Return the longest vibration period `2π/sqrt(λmin)` in a real mechanical
system. Supports the same boundary-condition, MPC and solver keywords as
`smallestEigenValue`. A zero eigenvalue yields `Inf`, as in the legacy formula.
"""
function largestPeriodTime(
    K::Union{SystemMatrix,SymmetricSystemMatrix},
    M::Union{SystemMatrix,SymmetricSystemMatrix};
    kwargs...
)
    λmin = smallestEigenValue(K, M; kwargs...)
    return λmin == 0 ? Inf : 2π / sqrt(λmin)
end

"""
    smallestPeriodTime(K, M; kwargs...)

Return the shortest vibration period `2π/sqrt(λmax)` in a real mechanical
system. Supports the same boundary-condition, MPC and solver keywords as
`largestEigenValue`.
"""
function smallestPeriodTime(
    K::Union{SystemMatrix,SymmetricSystemMatrix},
    M::Union{SystemMatrix,SymmetricSystemMatrix};
    kwargs...
)
    λmax = largestEigenValue(K, M; kwargs...)
    return λmax == 0 ? Inf : 2π / sqrt(λmax)
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
