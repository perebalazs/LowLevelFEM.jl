###############################################################
# Compound operators ##########################################
###############################################################

"""
    ChainSumTerm(chain, gauss)

Internal storage for one term of a compound operator expression.

`chain` is usually an `OpApplied` or a `MatrixChain`.
`gauss` stores the preferred integration rule for this term.
"""
struct ChainSumTerm
    chain::Any
    gauss::Any
end

"""
    ChainSum(terms)

Symbolic sum of operator chains.

Example:

```julia
B = Grad(Pu) ⋅ A + Pu ⋅ G
```

stores two terms without assembling anything.
"""
struct ChainSum
    terms::Vector{ChainSumTerm}
end

"""
TransposedChain(chain)

Symbolic transpose marker for compound expressions.

This does not assemble or numerically transpose anything immediately.
It is interpreted during expansion of the bilinear form.
"""
struct TransposedChain
    chain::Any
end

"""
CompoundChain(left, mats)

Temporary object representing

```
left ⋅ M1 ⋅ M2 ⋯
```

before the right-hand side is supplied.
"""
struct CompoundChain
    left::Any
    mats::Vector{Any}
end

"""
CompoundBilinear(left, mats, right)

Temporary object representing

```
left ⋅ M1 ⋅ M2 ⋯ ⋅ right
```

where `left` and/or `right` may be `ChainSum`.
"""
struct CompoundBilinear
    left::Any
    mats::Vector{Any}
    right::Any
end

#=
struct CompoundBilinearWeight
    CB::CompoundBilinear
    weight::Any
end

*(A::CompoundBilinear, B) = CompoundBilinearWeight(A, B)

struct LinearTermWeight
    CB::LinearTerm
    weight::Any
end

*(A::LinearTerm, B) = LinearTermWeight(A, B)
=#

struct WeightedTerm
    term::Any
    weight::Any
end

promote_term(t::WeightedTerm) = WeakTerm(t.weight, t.term)

Base.:+(a::WeightedTerm, b::WeightedTerm) =
    WeakSum(promote_term(a), promote_term(b))

Base.:+(a::WeightedTerm, b::BilinearTerm) =
    WeakSum(promote_term(a), promote_term(b))

Base.:+(a::BilinearTerm, b::WeightedTerm) =
    WeakSum(promote_term(a), promote_term(b))

Base.:-(t::WeightedTerm) =
    WeakTerm(-t.weight, t.term)

Base.:-(a::WeightedTerm, b::WeightedTerm) = a + (-b)
Base.:-(a::WeightedTerm, b::BilinearTerm) = a + (-b)
Base.:-(a::BilinearTerm, b::WeightedTerm) = a + (-b)

Base.:*(t::Union{BilinearTerm,LinearTerm,CompoundBilinear}, w::Union{Number,ScalarField}) =
    WeightedTerm(t, w)

Base.:*(w::Union{Number,ScalarField}, t::Union{BilinearTerm,LinearTerm,CompoundBilinear}) =
    WeightedTerm(t, w)

function ∫(wt::WeightedTerm;
    Ω=nothing, Γ=nothing, weight=nothing,
    gauss=:full,
    threads=:auto,
    assembly::Symbol=:csc,
    csc_matrix=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto)

    weight === nothing || error("Weight specified both as `* weight` and keyword `weight=...`.")

    if wt.term isa LinearTerm
        assembly in (:csc, :ijv) || error(
            "Unsupported matrix assembly mode $assembly. Expected :csc or :ijv."
        )
        csc_matrix === nothing || error(
            "csc_matrix is only available for bilinear forms."
        )

        return ∫(
            wt.term;
            Ω=Ω,
            Γ=Γ,
            weight=wt.weight,
            gauss=gauss,
            threads=threads
        )
    end

    return ∫(
        wt.term;
        Ω=Ω,
        Γ=Γ,
        weight=wt.weight,
        gauss=gauss,
        threads=threads,
        assembly=assembly,
        csc_matrix=csc_matrix,
        element_chunk_size=element_chunk_size
    )
end

withgauss(chain, gauss) = ChainSumTerm(chain, gauss)

"""
    full(chain)
    reduced(chain)

Attach a preferred quadrature rule to one term of a compound operator sum.

Example:

    B = full(A ⋅ Grad(Pu)) + reduced(G ⋅ Pu)
    K = ∫(B' ⋅ D ⋅ B; Ω="body", gauss=:auto)

With `gauss=:auto`, equal-rule products keep their rule; mixed products
default to `:full`.
"""
full(chain) = withgauss(chain, :full)
"""
    reduced(chain)

Mark one compound-operator term for reduced quadrature. Combine with
[`full`](@ref) and assemble the expanded bilinear expression using
`gauss=:auto` to retain per-term rules; mixed-rule products use full
quadrature. This selects integration, not `Problem(reducedOrder=true)`.
"""
reduced(chain) = withgauss(chain, :reduced)

_as_sumterm(t::ChainSumTerm) = t
_as_sumterm(x) = ChainSumTerm(x, :full)

_chain_sum(xs...) = ChainSum([_as_sumterm(x) for x in xs])

_is_chain_object(x) =
    x isa OpApplied ||
    x isa MatrixChain ||
    x isa ChainSum ||
    x isa ChainSumTerm ||
    x isa TransposedChain

function _coeff_mats(c)
    if c isa AbstractMatrix
        _check_coeff_matrix(c)
        return Any[c]
    elseif c isa TensorField
        return Any[tensorfield_to_matrix(c)]
    elseif c isa Number || c isa ScalarField
        return Any[c]
    else
        error("Compound operator: invalid coefficient type $(typeof(c)).")
    end
end

function _op_and_mats(x)
    if x isa OpApplied
        return x, Any[]
    elseif x isa MatrixChain
        return x.a, copy(x.mats)
    else
        error("Compound operator term must be OpApplied or MatrixChain, got $(typeof(x)).")
    end
end

struct _SideTerm
    chain::Any
    transposed::Bool
    gauss::Any
end

function _side_terms(x; transposed=false)
    if x isa TransposedChain
        return _side_terms(x.chain; transposed=(!transposed))
    elseif x isa ChainSum
        return [_SideTerm(t.chain, transposed, t.gauss) for t in x.terms]

    elseif x isa ChainSumTerm
        return [_SideTerm(x.chain, transposed, x.gauss)]

    elseif x isa OpApplied || x isa MatrixChain
        return [_SideTerm(x, transposed, :full)]

    else
        error("Compound operator: invalid side type $(typeof(x)).")
    end
end

function _maybe_transpose_mats(mats, transposed::Bool)
    if transposed
        return Any[adjoint(M) for M in reverse(mats)]
    else
        return copy(mats)
    end
end

function _select_gauss(left::_SideTerm, right::_SideTerm, gauss)
    gauss !== :auto && return gauss
    if left.gauss == right.gauss
        return left.gauss
    end

    return :full
end

function _compress_mats(mats)
    isempty(mats) && return 1.0

    C = mats[1]
    for i in 2:length(mats)
        C = C * mats[i]
    end

    return C
end

function _make_bilinear_term(left::_SideTerm, middle::Vector{Any}, right::_SideTerm)
    left_op, left_mats = _op_and_mats(left.chain)
    right_op, right_mats = _op_and_mats(right.chain)
    mats = Any[
        #_maybe_transpose_mats(left_mats, left.transposed)...,
        left_mats...,
        middle...,
        #_maybe_transpose_mats(right_mats, right.transposed)...
        adjoint.(reverse(right_mats))...
    ]

    coef = isempty(mats) ? 1.0 : mats

    #coef = _compress_mats(mats)

    return BilinearTerm(left_op, coef, right_op)
end

#function _coeff_mats(c)
#    if c isa AbstractMatrix
#        return Any[c]
#    elseif c isa Number || c isa ScalarField || c isa TensorField
#        return Any[c]
#    else
#        error("Compound operator: invalid coefficient type $(typeof(c)).")
#    end
#end

# ---------------------------------------------------------------------------

# Addition: build ChainSum

# ---------------------------------------------------------------------------

+(a::ChainSum, b::ChainSum) =
    ChainSum([a.terms...; b.terms...])

+(a::ChainSum, b::Union{OpApplied,MatrixChain,ChainSumTerm,TransposedChain}) =
    ChainSum([a.terms...; _as_sumterm(b)])

+(a::Union{OpApplied,MatrixChain,ChainSumTerm,TransposedChain}, b::ChainSum) =
    ChainSum([_as_sumterm(a); b.terms...])

+(a::Union{OpApplied,MatrixChain,ChainSumTerm,TransposedChain},
    b::Union{OpApplied,MatrixChain,ChainSumTerm,TransposedChain}) =
    _chain_sum(a, b)

# ---------------------------------------------------------------------------

# Symbolic transpose

# ---------------------------------------------------------------------------

adjoint(x::Union{ChainSum,ChainSumTerm,MatrixChain,OpApplied}) =
    TransposedChain(x)

transpose(x::Union{ChainSum,ChainSumTerm,MatrixChain,OpApplied}) =
    TransposedChain(x)

# ---------------------------------------------------------------------------

# Dot products involving ChainSum / TransposedChain

# ---------------------------------------------------------------------------

function ⋅(left::Union{ChainSum,TransposedChain,ChainSumTerm}, c)
    return CompoundChain(left, _coeff_mats(c))
end

function ⋅(cc::CompoundChain, c)
    append!(cc.mats, _coeff_mats(c))
    return cc
end

function ⋅(cc::CompoundChain, right::Union{ChainSum,TransposedChain,ChainSumTerm,OpApplied,MatrixChain})
    return CompoundBilinear(cc.left, cc.mats, right)
end

function ⋅(left::Union{ChainSum,TransposedChain,ChainSumTerm},
    right::Union{ChainSum,TransposedChain,ChainSumTerm,OpApplied,MatrixChain})
    return CompoundBilinear(left, Any[], right)
end

function ⋅(left::Union{OpApplied,MatrixChain},
    right::Union{ChainSum,TransposedChain,ChainSumTerm})
    return CompoundBilinear(left, Any[], right)
end

function ⋅(A::AbstractMatrix, op::OpApplied)
    return op ⋅ adjoint(A)
end

function ⋅(A::AbstractMatrix, ch::MatrixChain)
    return ch ⋅ adjoint(A)
end

function ⋅(A::Union{Number,ScalarField}, op::OpApplied)
    return op ⋅ A
end

# ---------------------------------------------------------------------------

# Compound block routing / assembly

# ---------------------------------------------------------------------------

"""
    _side_problem(term::_SideTerm)

Return the `Problem` associated with one compound-side term.
"""
@inline function _side_problem(term::_SideTerm)
    op, _ = _op_and_mats(term.chain)
    return op.P
end

"""
    _group_side_terms(terms)

Group compound-side terms by `Problem` while preserving the order of first
occurrence.

Returns

    problems, groups

where `problems[i]` owns all terms stored in `groups[i]`.
"""
function _group_side_terms(terms::AbstractVector{_SideTerm})
    problems = Problem[]
    groups = Vector{Vector{_SideTerm}}()

    for term in terms
        P = _side_problem(term)
        idx = findfirst(Q -> Q === P, problems)

        if idx === nothing
            push!(problems, P)
            push!(groups, _SideTerm[term])
        else
            push!(groups[idx], term)
        end
    end

    return problems, groups
end

"""
    _assemble_compound_block(left_terms, right_terms, middle; ...)

Assemble one test/trial block of a compound bilinear form.

All terms in `left_terms` must belong to one test `Problem`, and all terms in
`right_terms` must belong to one trial `Problem`. Within a block, the existing
compound assembly strategy is retained: direct CSC assembly builds one
structural pattern and reuses it, together with the worker value buffers and
node-coordinate table, for every expanded term pair.

The optional `csc_matrix` is therefore a block-local pattern. It is used by the
single-block compound path; multifield compound assembly creates a separate
pattern for every block.
"""
function _assemble_compound_block(
    left_terms::AbstractVector{_SideTerm},
    right_terms::AbstractVector{_SideTerm},
    middle::Vector{Any};
    dom=nothing,
    weight=nothing,
    gauss=:auto,
    threads=:auto,
    assembly::Symbol=:csc,
    csc_matrix=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto
)
    isempty(left_terms) && error("Cannot assemble an empty compound test block.")
    isempty(right_terms) && error("Cannot assemble an empty compound trial block.")

    first_bt = _make_bilinear_term(left_terms[1], middle, right_terms[1])
    Pu = first_bt.b.P
    Ps = first_bt.a.P

    if assembly === :csc
        num_threads = resolve_num_threads(threads)

        if csc_matrix === nothing
            Kcsc = build_csc_pattern(Pu, Ps; domain=dom)
        else
            csc_matrix isa SparseMatrixCSC{Float64,Int} || error(
                "csc_matrix must be SparseMatrixCSC{Float64,Int}; " *
                "got $(typeof(csc_matrix))."
            )
            size(csc_matrix) == (ndofs(Ps), ndofs(Pu)) || error(
                "csc_matrix has size $(size(csc_matrix)); expected " *
                "$((ndofs(Ps), ndofs(Pu)))."
            )
            Kcsc = csc_matrix
            _validate_csc_component_pattern(Kcsc, Pu.pdim, Ps.pdim)
        end

        coordinates = dense_node_coordinates(Pu)

        nzval_buffers = Vector{Vector{Float64}}(undef, num_threads)
        nzval_buffers[1] = Kcsc.nzval
        for worker in 2:num_threads
            nzval_buffers[worker] = zeros(Float64, length(Kcsc.nzval))
        end

        for LL in left_terms
            for R in right_terms
                bt = _make_bilinear_term(LL, middle, R)

                bt.a.P === Ps || error(
                    "Compound block contains more than one test Problem."
                )
                bt.b.P === Pu || error(
                    "Compound block contains more than one trial Problem."
                )

                g = _select_gauss(LL, R, gauss)

                Kterm = assemble_operator(
                    bt.b.P,
                    bt.b.op,
                    bt.a.P,
                    bt.a.op;
                    coefficient = bt.coef,
                    domain = dom,
                    weight = weight,
                    gauss = g,
                    assembly = :csc,
                    threads = num_threads,
                    K = Kcsc,
                    node_coordinates = coordinates,
                    _nzval_buffers = nzval_buffers,
                    element_chunk_size = element_chunk_size
                )

                Kcsc = Kterm.A
            end
        end

        # Preserve a complete externally supplied reusable pattern. Internally
        # created block patterns are finalized in the usual way.
        csc_matrix === nothing && dropzeros!(Kcsc)

        return SystemMatrix(Kcsc, Pu, Ps)
    end

    # ------------------------------------------------------------------
    # Legacy IJV path, local to one test/trial block
    # ------------------------------------------------------------------

    first = true

    I = nothing
    J = nothing
    V = nothing
    nr = 0
    nc = 0

    for LL in left_terms
        for R in right_terms
            bt = _make_bilinear_term(LL, middle, R)

            bt.a.P === Ps || error(
                "Compound block contains more than one test Problem."
            )
            bt.b.P === Pu || error(
                "Compound block contains more than one trial Problem."
            )

            g = _select_gauss(LL, R, gauss)

            if first
                I, J, V, nr, nc = assemble_operator(
                    bt.b.P,
                    bt.b.op,
                    bt.a.P,
                    bt.a.op;
                    coefficient = bt.coef,
                    domain = dom,
                    weight = weight,
                    gauss = g,
                    assembly = :triplets,
                    threads = threads
                )

                first = false
            else
                assemble_operator(
                    bt.b.P,
                    bt.b.op,
                    bt.a.P,
                    bt.a.op;
                    coefficient = bt.coef,
                    domain = dom,
                    weight = weight,
                    gauss = g,
                    assembly = :add,
                    V = V,
                    threads = threads
                )
            end
        end
    end

    K = sparse(I, J, V, nr, nc)
    dropzeros!(K)

    return SystemMatrix(K, Pu, Ps)
end

# ---------------------------------------------------------------------------

# Integral wrapper

# ---------------------------------------------------------------------------

"""
    ∫(expr::CompoundBilinear; Ω=nothing, Γ=nothing, weight=nothing,
      gauss=:auto, threads=:auto, assembly=:csc,
      csc_matrix=nothing, element_chunk_size=:auto)

Expand and assemble a compound bilinear operator expression.

## Single-block compound forms

If all test-side terms belong to one `Problem` and all trial-side terms belong
to one `Problem`, the original compound assembly path is used. For
`assembly=:csc`, one structural pattern is built for the block and reused for
all expanded left/right term pairs, together with the final `nzval` array,
worker-local value buffers, and dense node-coordinate array.

A supplied `csc_matrix` is supported for this single-block case. Its complete
structural pattern is preserved; reset `csc_matrix.nzval` before a new
independent assembly, but not between internally expanded terms.

## Multifield compound forms

If either side contains terms associated with multiple `Problem`s, terms are
first grouped by `Problem` in order of first occurrence. The test-side groups
define block rows and the trial-side groups define block columns. Each block is
then assembled independently from the Cartesian product of the terms belonging
to that test/trial pair, and the completed blocks are combined with
`SystemMatrix(blocks)`.

For a square multifield system, test and trial fields must appear in the same
order. The user is responsible for writing the compound expression with a
consistent field order. For example,

```julia
B = Au ⋅ Grad(Pu) + Aφ ⋅ Grad(Pφ) + Gφ ⋅ Id(Pφ)
K = ∫(B' ⋅ D ⋅ B; Γ="beam")
```

is grouped as `[Pu, Pφ]` on both sides and automatically produces the four
blocks `Kuu`, `Kuφ`, `Kφu`, and `Kφφ`.

With direct CSC assembly, every multifield block builds its own structural
pattern because blocks may have different dimensions. Compound terms inside a
single block still reuse that block's pattern and worker buffers exactly as in
the original homogeneous compound path. Consequently, `csc_matrix` is not
accepted for a multifield compound form; external pattern reuse remains a
single-block option.

Matrix chains and matrices containing `Number` and `ScalarField` entries use
the same coefficient evaluation and element kernel as ordinary bilinear terms.
Each compound term may select `full(...)` or `reduced(...)` quadrature; with
`gauss=:auto`, mixed-rule products use full integration.

Use `assembly=:ijv` to select the legacy triplet path. Multifield IJV assembly
is likewise performed independently block by block.
"""
function ∫(expr::CompoundBilinear;
    Ω=nothing, Γ=nothing, weight=nothing,
    gauss=:auto,
    threads=:auto,
    assembly::Symbol=:csc,
    csc_matrix=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto)

    gauss === :gll && error(
        "GLL quadrature for compound or multifield bilinear forms is not yet implemented. " *
        "Use a single spectral term of the form P ⋅ P."
    )

    assembly ∈ (:ijv, :csc) || error(
        "Unsupported matrix assembly mode $assembly. " *
        "Expected :csc or :ijv."
    )

    assembly === :ijv && csc_matrix !== nothing && error(
        "csc_matrix is only valid with assembly=:csc."
    )

    left_terms = _side_terms(expr.left)
    right_terms = _side_terms(expr.right)

    isempty(left_terms) && error("Cannot assemble an empty compound bilinear form.")
    isempty(right_terms) && error("Cannot assemble an empty compound bilinear form.")

    left_problems, left_groups = _group_side_terms(left_terms)
    right_problems, right_groups = _group_side_terms(right_terms)

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    # ------------------------------------------------------------------
    # Original homogeneous / single-block path
    # ------------------------------------------------------------------

    if length(left_problems) == 1 && length(right_problems) == 1
        return _assemble_compound_block(
            left_terms,
            right_terms,
            expr.mats;
            dom = dom,
            weight = weight,
            gauss = gauss,
            threads = threads,
            assembly = assembly,
            csc_matrix = csc_matrix,
            element_chunk_size = element_chunk_size
        )
    end

    # ------------------------------------------------------------------
    # Multifield block routing
    # ------------------------------------------------------------------

    csc_matrix === nothing || error(
        "csc_matrix cannot represent a multifield compound system because " *
        "different blocks may have different dimensions. Each block builds " *
        "and reuses its own CSC pattern."
    )

    length(left_problems) == length(right_problems) || error(
        "Compound multifield system is not square in field sense: " *
        "$(length(left_problems)) test fields and " *
        "$(length(right_problems)) trial fields were found."
    )

    for i in eachindex(left_problems)
        left_problems[i] === right_problems[i] || error(
            "Compound multifield field ordering mismatch at position $i. " *
            "Test and trial fields must appear in the same order."
        )
    end

    nfields = length(left_problems)
    blocks = Matrix{SystemMatrix}(undef, nfields, nfields)

    for i in 1:nfields
        for j in 1:nfields
            blocks[i, j] = _assemble_compound_block(
                left_groups[i],
                right_groups[j],
                expr.mats;
                dom = dom,
                weight = weight,
                gauss = gauss,
                threads = threads,
                assembly = assembly,
                csc_matrix = nothing,
                element_chunk_size = element_chunk_size
            )
        end
    end

    return SystemMatrix(blocks)
end

#∫(expr::CompoundBilinearWeight; Ω=nothing, Γ=nothing, gauss=:auto) = ∫(expr.CB; Ω=Ω, Γ=Γ, weight=expr.weight, gauss=gauss)
#∫(expr::LinearTermWeight; Ω=nothing, Γ=nothing, gauss=:auto) = ∫(expr.CB; Ω=Ω, Γ=Γ, weight=expr.weight, gauss=gauss)

########################################################################

