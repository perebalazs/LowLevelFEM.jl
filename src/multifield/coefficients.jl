struct DomainSpec
    kind::Symbol     # :Ω vagy :Γ
    name::String
end

Base.show(io::IO, d::DomainSpec) =
    print(io, "$(d.kind)=\"$(d.name)\"")

#=
function _covers_domain(coef::ScalarField, domain::DomainSpec)
    dimTags = gmsh.model.getEntitiesForPhysicalName(domain.name)
    elems = Int[]

    for (edim, etag) in dimTags
        elemTypes, elemTags, _ = gmsh.model.mesh.getElements(edim, etag)
        for tags in elemTags
            append!(elems, Int.(tags))
        end
    end

    return all(e -> e in coef.numElem, elems)
end
=#

function _covers_domain(coef::ScalarField, domain::DomainSpec)
    selected = Set(coef.numElem)

    for (edim, etag) in gmsh.model.getEntitiesForPhysicalName(domain.name)
        _, elementTags, _ = gmsh.model.mesh.getElements(edim, etag)

        for tags in elementTags
            for tag in tags
                tag in selected || return false
            end
        end
    end

    return true
end

@inline function _build_elemwise_coeff_dict(coef::ScalarField, domain::DomainSpec)
    if isElementwise(coef) && _covers_domain(coef, domain)
        p = coef
    else
        p = elementsToNodes(coef)
        p = nodesToElements(p, onPhysicalGroup=domain.name)
    end

    return Dict(zip(p.numElem, p.A))
end

#=
@inline function _build_elemwise_coeff_dict(coef::ScalarField, domain::DomainSpec)
    p = elementsToNodes(coef)
    p = nodesToElements(p, onPhysicalGroup=domain.name)
    return Dict(zip(p.numElem, p.A))  # elemTag => coeff nodal vector(s)
end
=#

"""
    _prepare_coefficient(C, domain)

Prepare coefficient for assembly.

Returns one of:

Float64
Dict(elem => nodal array)
Matrix{Any} with entries Float64 or Dict(elem => nodal array)
"""
function _prepare_coefficient(C, domain)

    # scalar constant
    if C isa Number
        return Float64(C)
    #end

    # scalar field
    elseif C isa ScalarField
        return _build_elemwise_coeff_dict(C, domain)
    #end

    # tensor coefficient
    elseif C isa AbstractMatrix

        # Constant numeric tensors are reused at every Gauss point. Keeping a
        # concrete Float64 matrix here avoids one matrix allocation per point.
        if eltype(C) <: Number || all(x -> x isa Number, C)
            return Float64.(C)
        end

        W = Matrix{Any}(undef, size(C)...)

        for I in CartesianIndices(C)

            cij = C[I]

            if cij isa Number
                W[I] = Float64(cij)

            elseif cij isa ScalarField
                W[I] = _build_elemwise_coeff_dict(cij, domain)

            else
                error("Matrix coefficient entries must be Number or ScalarField, got $(typeof(cij))")
            end

        end

        return W
    #end

    elseif C isa AbstractVector
        mats = [_prepare_coefficient(M, domain) for M in C]
        return mats

    #elseif C isa VectorField
    #    if C.type == :v2D
    #        nc = 2
    #    elseif C.type == :v3D
    #        nc = 3
    #    else
    #        error("Unknown VectorField type $(C.type)")
    #    end
    #    return [_prepare_coefficient(C[i]) for i in 1:nc]
    #
    #elseif C isa TensorField
    #    return [_prepare_coefficient(C[i]) for i in 1:9]

    end

    error("Unsupported coefficient type $(typeof(C))")
end

#@inline function _coeff_at_gp(pa::Dict{<:Integer,<:AbstractMatrix}, elem::Integer, hcol::AbstractVector)
#    return dot(view(pa[elem], :, 1), hcol)
#end

@inline function _coeff_at_gp(
    pa::Dict{K, V},
    elem::Integer,
    hcol::AbstractVector) where {K<:Integer,V<:AbstractMatrix}

    vals = get(pa, Int(elem), nothing)

    vals === nothing && return 0.0

    s = 0.0

    @inbounds for i in eachindex(hcol)
        s += vals[i, 1] * hcol[i]
    end

    return s

    #return dot(view(vals, :, 1), hcol)
end

@inline function _coeff_at_gp_old(
    pa::Dict{<:Integer,<:AbstractMatrix},
    elem::Integer,
    hcol::AbstractVector)

    return dot(view(pa[Int(elem)], :, 1), hcol)
end

"""
    _eval_tensor_at_gp(C, pa, elem, hgp)

Evaluate scalar or tensor coefficient at Gauss point.

Supports:
    Number
    ScalarField
    Matrix{Number}
    Matrix{ScalarField}
"""
function _eval_tensor_at_gp(C, pa, elem, hgp)

    # scalar constant
    if C isa Number
        return C
    end

    # scalar field
    if C isa ScalarField
        return _coeff_at_gp(pa, elem, hgp)
    end

    # tensor
    if C isa AbstractMatrix
        m, n = size(C)
        W = zeros(Float64, m, n)

        for i in 1:m, j in 1:n
            cij = C[i, j]

            if cij isa Number
                W[i, j] = cij
            else
                # ScalarField
                W[i, j] = _coeff_at_gp(pa, elem, hgp)
            end

        end

        return W
    end

    error("Unsupported coefficient type $(typeof(C))")

end

"""
    _eval_coefficient_at_gp(Cprep, elem, hgp)

Evaluate prepared coefficient at Gauss point.

Supports:
Float64
Dict(elem => nodal array)
Matrix{Any} with entries Float64 or Dict(elem => nodal array)
"""
function _eval_coefficient_at_gp(Cprep, elem, hgp)

    # scalar constant
    if Cprep isa Number
        return Cprep
    end

    # scalar field
    if Cprep isa Dict
        return _coeff_at_gp(Cprep, elem, hgp)
    end

    # tensor coefficient
    if Cprep isa AbstractMatrix

        eltype(Cprep) <: Number && return Cprep

        m, n = size(Cprep)
        W = Matrix{Float64}(undef, m, n)

        @inbounds for I in CartesianIndices(Cprep)
            cij = Cprep[I]

            if cij isa Number
                W[I] = Float64(cij)

            elseif cij isa Dict{Int,Matrix{Float64}}
                # Fast path: concrete coefficient dictionary
                W[I] = _coeff_at_gp(cij, elem, hgp)

            elseif cij isa AbstractDict
                # Generic fallback
                W[I] = _coeff_at_gp(cij, elem, hgp)

            else
                error("Unsupported prepared coefficient type $(typeof(cij))")
            end
        end
        #=
        for I in CartesianIndices(Cprep)

            cij = Cprep[I]

            if cij isa Number
                W[I] = cij
            else
                W[I] = _coeff_at_gp(cij, elem, hgp)
            end

        end
        =#

        return W
    end

    # vector coefficient
    if Cprep isa AbstractVector

        n = length(Cprep)
        W = Vector{Float64}(undef, n)

        for i in 1:n
        
            cij = Cprep[i]

            if cij isa Number
                W[i] = cij
            else
                W[i] = _coeff_at_gp(cij, elem, hgp)
            end

        end

        return W
    end

    error("Unsupported prepared coefficient type $(typeof(Cprep))")
end


# ---------------------------------------------------------------------------
# Small local matrix multiplication with optional structural sparsity
# ---------------------------------------------------------------------------
#
# This layer is intentionally hidden below the weak-form DSL.  The user can keep
# writing paper-like expressions such as
#
#     ∫( SymGrad(Pu) ⋅ C ⋅ SymGrad(Pu) )
#
# while the local Gauss-point products try to avoid multiplications by known
# zero entries.  If a factor is not sparse enough, the code falls back to dense
# BLAS-backed mul!.

mutable struct MatrixPattern
    rows::Vector{Int}
    cols::Vector{Int}
    nnz::Int
    m::Int
    n::Int
end

MatrixPattern(m::Integer, n::Integer) = MatrixPattern(Int[], Int[], 0, Int(m), Int(n))

@inline density(p::MatrixPattern) = p.nnz / (p.m * p.n)

"""
    detect_pattern(M; tol=0.0)

Detect the nonzero pattern of a small local matrix.

The default `tol=0.0` is deliberately strict: entries that are exactly zero are
treated as structural zeros, but accidentally small nonzero entries are kept.
For finite element operator matrices this is usually the safe choice.
"""
function detect_pattern(M::AbstractMatrix; tol::Real=0.0)
    m, n = size(M)
    rows = Int[]
    cols = Int[]
    sizehint!(rows, length(M))
    sizehint!(cols, length(M))

    @inbounds for j in 1:n
        for i in 1:m
            mij = M[i, j]
            if tol == 0
                if !iszero(mij)
                    push!(rows, i)
                    push!(cols, j)
                end
            else
                if abs(mij) > tol
                    push!(rows, i)
                    push!(cols, j)
                end
            end
        end
    end

    return MatrixPattern(rows, cols, length(rows), m, n)
end

"""
    detect_pattern!(pattern, M; tol=0.0)

Update a reusable matrix pattern from the current numerical values of `M`.
"""
function detect_pattern!(
    pattern::MatrixPattern,
    M::AbstractMatrix;
    tol::Real = 0.0
)
    m, n = size(M)

    pattern.m == m && pattern.n == n ||
        error(
            "detect_pattern!: matrix size changed from " *
            "($(pattern.m), $(pattern.n)) to ($m, $n)."
        )

    empty!(pattern.rows)
    empty!(pattern.cols)

    @inbounds for j in 1:n
        for i in 1:m
            mij = M[i, j]

            if tol == 0
                if !iszero(mij)
                    push!(pattern.rows, i)
                    push!(pattern.cols, j)
                end
            elseif abs(mij) > tol
                push!(pattern.rows, i)
                push!(pattern.cols, j)
            end
        end
    end

    pattern.nnz = length(pattern.rows)

    return pattern
end

@inline function sparse_enough(
    p::MatrixPattern;
    density_limit::Real = 0.40,
    min_saved_entries::Integer = 2
)
    total = p.m * p.n
    return (total - p.nnz) >= min_saved_entries && density(p) <= density_limit
end

@inline function _scale_or_zero!(R, beta)
    if beta == 0
        fill!(R, 0.0)
    elseif beta != 1
        R .*= beta
    end
    return R
end

"""
    mul_dense_sparse!(R, A, B, patB, alpha=1.0, beta=0.0)

Compute `R = alpha*A*B + beta*R`, using only the nonzero entries of `B`.
"""
function mul_dense_sparse!(
    R::AbstractMatrix,
    A::AbstractMatrix,
    B::AbstractMatrix,
    patB::MatrixPattern,
    alpha::Number = 1.0,
    beta::Number = 0.0
)
    _scale_or_zero!(R, beta)

    @inbounds for q in 1:patB.nnz
        k = patB.rows[q]
        j = patB.cols[q]
        bkj = B[k, j]
        iszero(bkj) && continue

        αb = alpha * bkj
        for i in axes(A, 1)
            R[i, j] += A[i, k] * αb
        end
    end

    return R
end

"""
    mul_sparse_dense!(R, A, B, patA, alpha=1.0, beta=0.0)

Compute `R = alpha*A*B + beta*R`, using only the nonzero entries of `A`.
"""
function mul_sparse_dense!(
    R::AbstractMatrix,
    A::AbstractMatrix,
    B::AbstractMatrix,
    patA::MatrixPattern,
    alpha::Number = 1.0,
    beta::Number = 0.0
)
    _scale_or_zero!(R, beta)

    @inbounds for q in 1:patA.nnz
        i = patA.rows[q]
        k = patA.cols[q]
        aik = A[i, k]
        iszero(aik) && continue

        αa = alpha * aik
        for j in axes(B, 2)
            R[i, j] += αa * B[k, j]
        end
    end

    return R
end

"""
    mul_sparse_sparse!(R, A, B, patA, patB, alpha=1.0, beta=0.0)

Compute `R = alpha*A*B + beta*R`, using the nonzero entries of both factors.

This is fully general and does not assume symmetry.  The matrices are expected
to be small local FE matrices, so a simple pattern loop is often adequate.
"""
function mul_sparse_sparse!(
    R::AbstractMatrix,
    A::AbstractMatrix,
    B::AbstractMatrix,
    patA::MatrixPattern,
    patB::MatrixPattern,
    alpha::Number = 1.0,
    beta::Number = 0.0
)
    _scale_or_zero!(R, beta)

    @inbounds for qa in 1:patA.nnz
        i = patA.rows[qa]
        k = patA.cols[qa]
        aik = A[i, k]
        iszero(aik) && continue

        αa = alpha * aik

        for qb in 1:patB.nnz
            kb = patB.rows[qb]
            kb == k || continue

            j = patB.cols[qb]
            bkj = B[kb, j]
            iszero(bkj) && continue

            R[i, j] += αa * bkj
        end
    end

    return R
end

"""
    mul_opt!(R, A, B; patA=nothing, patB=nothing, alpha=1.0, beta=0.0, ...)

General local matrix product with sparse-pattern acceleration and dense fallback.

No symmetry is assumed.  If neither side is sparse enough according to the
detected/provided pattern, this falls back to `mul!(R, A, B, alpha, beta)`.
"""
function mul_opt!(
    R::AbstractMatrix,
    A::AbstractMatrix,
    B::AbstractMatrix;
    patA::Union{MatrixPattern,Nothing}=nothing,
    patB::Union{MatrixPattern,Nothing}=nothing,
    alpha::Number = 1.0,
    beta::Number = 0.0,
    density_limit::Real = 0.40,
    tol::Real = 0.0
)
    patA === nothing && (patA = detect_pattern(A; tol=tol))
    patB === nothing && (patB = detect_pattern(B; tol=tol))

    sparseA = sparse_enough(patA; density_limit=density_limit)
    sparseB = sparse_enough(patB; density_limit=density_limit)

    if sparseA && sparseB
        return mul_sparse_sparse!(R, A, B, patA, patB, alpha, beta)
    elseif sparseA
        return mul_sparse_dense!(R, A, B, patA, alpha, beta)
    elseif sparseB
        return mul_dense_sparse!(R, A, B, patB, alpha, beta)
    else
        return mul!(R, A, B, alpha, beta)
    end
end

# Convenience allocation wrapper used in coefficient chains.
function mul_opt(
    A::AbstractMatrix,
    B::AbstractMatrix;
    patA::Union{MatrixPattern,Nothing}=nothing,
    patB::Union{MatrixPattern,Nothing}=nothing,
    density_limit::Real = 0.40,
    tol::Real = 0.0
)
    R = Matrix{promote_type(eltype(A), eltype(B))}(undef, size(A, 1), size(B, 2))
    return mul_opt!(R, A, B;
        patA=patA,
        patB=patB,
        alpha=1.0,
        beta=0.0,
        density_limit=density_limit,
        tol=tol
    )
end

# Backward-compatible helper names kept for experiments/notebooks.
function sparse_like_mul(A, B)
    patA = detect_pattern(A)
    patB = detect_pattern(B)
    C = zeros(promote_type(eltype(A), eltype(B)), size(A,1), size(B,2))
    return mul_sparse_sparse!(C, A, B, patA, patB)
end

function _small_matmul(A, B)
    if A isa Number || B isa Number
        return A * B
    end
    return mul_opt(A, B)
end

function _eval_coefficient_chain_at_gp(Cprep::AbstractVector, elem, hgp)
    Cgp = _eval_coefficient_at_gp(Cprep[1], elem, hgp)

    @inbounds for i in 2:length(Cprep)
        Mi = _eval_coefficient_at_gp(Cprep[i], elem, hgp)

        if Cgp isa Number || Mi isa Number
            Cgp = Cgp * Mi
        else
            Cgp = mul_opt(Cgp, Mi)
        end
    end

    return Cgp
end

const PreparedChainFactor = Union{Float64,Matrix{Float64},Matrix{Any},Dict{Int,Matrix{Float64}}}

# Prepared bilinear chains follow the same AbstractVector dispatch as raw chains.
struct CoefficientChainPlan <: AbstractVector{Any}
    factors::Vector{PreparedChainFactor}
end
Base.size(p::CoefficientChainPlan) = size(p.factors)
Base.getindex(p::CoefficientChainPlan, i::Int) = p.factors[i]
Base.IndexStyle(::Type{CoefficientChainPlan}) = IndexLinear()

struct CoefficientChainWorkspace
    plan::CoefficientChainPlan
    evaluated::Vector{Matrix{Float64}}
    products::Vector{Matrix{Float64}}
end

_is_constant_chain_factor(x) = x isa Number || (x isa AbstractMatrix && eltype(x) <: Number)
_normalize_chain_factor(x::AbstractMatrix{<:Number}) = Matrix{Float64}(x)
_normalize_chain_factor(x) = x

function _prepare_coefficient_chain(mats::AbstractVector)
    # Keep the old path for nonfinite literal data, including its propagation
    # of NaN/Inf through sparse-pattern-aware multiplication.
    for factor in mats
        if factor isa Number
            isfinite(factor) || return mats
        elseif factor isa AbstractMatrix{<:Number}
            all(isfinite,factor) || return mats
        end
    end
    factors = PreparedChainFactor[]
    constant_prefix = true
    for factor in mats
        factor = _normalize_chain_factor(factor)
        constant_prefix &= _is_constant_chain_factor(factor)
        # Only fold the leading literal prefix. Folding A*B in c(x)*A*B
        # reassociates floating-point products and can overflow/underflow even
        # when the original left-to-right evaluation remains finite.
        if constant_prefix && !isempty(factors)
            factors[end] = _small_matmul(last(factors), factor)
        else
            push!(factors, factor)
        end
    end
    length(factors) == 1 && return only(factors)
    any(x -> x isa AbstractMatrix, factors) || return mats
    return CoefficientChainPlan(factors)
end

function CoefficientChainWorkspace(plan::CoefficientChainPlan)
    evaluated = Matrix{Float64}[]
    products = Matrix{Float64}[]
    rows, cols = 0, 0
    for factor in plan.factors
        push!(evaluated, factor isa Matrix{Any} ? zeros(size(factor)) : zeros(0,0))
        if factor isa AbstractMatrix
            if rows == 0
                rows, cols = size(factor)
            else
                cols = size(factor, 2)
            end
        end
        push!(products, zeros(rows, cols))
    end
    return CoefficientChainWorkspace(plan, evaluated, products)
end

@inline _eval_coefficient_chain_at_gp!(ws, plan::AbstractVector, elem, hgp) =
    _eval_coefficient_chain_at_gp(plan, elem, hgp)

function _eval_coefficient_chain_at_gp!(ws, plan::CoefficientChainPlan, elem, hgp)
    # Owned by one assembly worker; no result or evaluated field buffer is shared.
    scratch = ws.chain_workspace
    if scratch === nothing || scratch.plan !== plan
        scratch = CoefficientChainWorkspace(plan)
        ws.chain_workspace = scratch
    end
    scalar = 1.0
    product = nothing
    @inbounds for i in eachindex(plan.factors)
        factor = plan.factors[i]
        if factor isa AbstractMatrix
            value = factor isa Matrix{Any} ?
                _eval_coefficient_at_gp!(scratch.evaluated[i], factor, elem, hgp) : factor
            out = scratch.products[i]
            if product === nothing
                @. out = scalar * value
            else
                mul!(out, product, value)
            end
            product = out
        else
            value = _eval_coefficient_at_gp(factor, elem, hgp)
            if product === nothing
                scalar *= value
            else
                out = scratch.products[i]
                @. out = product * value
                product = out
            end
        end
    end
    return product
end

@inline function surface_basis_from_J(Jac, k)
    a1 = SVector(Jac[1,3k-2], Jac[2,3k-2], Jac[3,3k-2])
    a2 = SVector(Jac[1,3k-1], Jac[2,3k-1], Jac[3,3k-1])

    t1 = a1 / norm(a1)

    a2o = a2 - dot(a2, t1) * t1
    n2 = norm(a2o)

    if n2 < 1e-12
        error("surface_basis_from_J: degenerate surface basis.")
    end

    t2 = a2o / n2

    return t1, t2
end

