# Weak form DSL layer for LowLevelFEM
###############################################################

# ------------------------------------------------------------
# Applied operator
# ------------------------------------------------------------

struct OpApplied
    P
    op
end

"""
    Grad(P)

Create the weak-form gradient of a `Problem`/`Field` descriptor. A scalar
field produces `P.dim` rows. A vector field produces `P.dim*P.pdim` rows in
component-major order: derivatives of component 1, then component 2, etc.
For 2D vectors the order is `[u₁,x, u₁,y, u₂,x, u₂,y]`.

This constructs a symbolic operator; it does not differentiate a solved field.
For field values use [`∇`](@ref). Example: `K = ∫(Grad(P) ⋅ Grad(P))`.
"""
Grad(P) = OpApplied(P, GradOp())

"""
    Div(P)

Construct the scalar divergence operator for a vector unknown with
`P.pdim == P.dim`. The output is `sum(∂uᵢ/∂xᵢ)`.
Example: `B = ∫(Q ⋅ Div(U))` couples scalar test field `Q` to vector trial
field `U`. Use `∇(u; nabla=:div)` to differentiate a computed field.
"""
Div(P) = OpApplied(P, DivOp())

"""
    Curl(P)

Construct the weak-form curl of a vector unknown (`P.pdim == P.dim`).
In 2D it has one component `∂u₂/∂x - ∂u₁/∂y`; in 3D it has the three
Cartesian curl components. Example: `K = ∫(Curl(P) ⋅ Curl(P))`.
For a computed vector field use `∇(u; nabla=:curl)`.
"""
Curl(P) = OpApplied(P, CurlOp())

"""
    SymGrad(P)

Construct the engineering-strain operator for a vector unknown with
`P.pdim == P.dim`. The component order is `[xx, yy, xy]` in 2D and
`[xx, yy, zz, xy, yz, zx]` in 3D. For `P.dim == P.pdim == 1`, the output is the single axial derivative.
Shear entries are engineering strains:
`γij = ∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ`, twice the corresponding tensor strain.

Use a constitutive matrix with the same order, for example
`K = ∫(SymGrad(U) ⋅ D(:Solid, mat) ⋅ SymGrad(U))`.
`ε` is an alias. A solved-field tensor strain is instead
`(u ∘ ∇ + ∇ ∘ u)/2`, whose off-diagonal components are tensor strains.
"""
SymGrad(P) = OpApplied(P, SymGradOp())

"""
    Id(P)

Create a weak-form DSL identity operator applied to `P`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object for identity mapping.

# Example
```julia
M = ∫(Id(Pu) ⋅ Id(Pu); Ω="solid")
```
"""
Id(P) = OpApplied(P, IdOp())

"""
    TensorDiv(P)

Construct the divergence of a 2D or 3D tensor unknown. `P.pdim` must be 9:
full 3×3 tensors are stored column-major as
`[xx, yx, zx, xy, yy, zy, xz, yz, zz]`. The output has `P.dim` components,
`(div T)ᵢ = sum(∂Tᵢⱼ/∂xⱼ, j=1:P.dim)`. In 2D only the in-plane block
and x/y derivatives contribute; no symmetry is assumed.
Example: `K = ∫(TensorDiv(P) ⋅ TensorDiv(P))`.
"""
TensorDiv(P) = OpApplied(P, TensorDivOp())

"""
    Adv(P)

Reserved convenience constructor for scalar advection in the weak-form DSL.
The current implementation calls the internal `AdvOp` without its required
velocity argument, so `Adv(P)` currently raises a `MethodError`.

Use an explicit coefficient chain for working scalar advection, for example
`A = ∫(P ⋅ [bx by] ⋅ Grad(P))` in 2D. A `1×2` coefficient row contracts
the gradient with the prescribed velocity. No solver fix is made by this
interface documentation.
"""
Adv(P) = OpApplied(P, AdvOp())

function _check_scalarfield(sf::ScalarField)

    if !isNodal(sf)
        error(
            "ScalarField must be nodal for integration.\n" *
            "Field: $(sf)\n" *
            "Only nodal fields are supported because Gauss values\n" *
            "are obtained via Lagrange interpolation."
        )
    end

end

"""
    AxialGrad(P)

Construct the axial derivative/strain operator along the local element tangent
`t`. For a scalar unknown this is `t ⋅ grad(s)`; for a vector unknown it is
`t ⋅ grad(u) ⋅ t`, a scalar axial strain. The direction is obtained from the
line-element geometry, not supplied as a separate user axis.

Typical use: `K = ∫(AxialGrad(U) ⋅ (E*A) ⋅ AxialGrad(U); Ω="beam")`.
[`TangentialGrad`](@ref) retains vector components of the directional derivative.
"""
AxialGrad(P) = OpApplied(P, AxialGradOp())

"""
    SurfaceGrad(P)

Surface gradient operator for fields defined on embedded surfaces.

Computes tangential derivatives on a 2D manifold embedded in 3D using
local tangent bases evaluated at Gauss points.

For scalar fields, returns the tangential surface gradient.
For vector fields, returns the full tangential displacement gradient.

The operator is orientation independent and invariant with respect to the
global coordinate system.

Typical applications include:

- membrane mechanics,
- surface PDEs,
- Laplace–Beltrami operators,
- surface diffusion.

See also: `SurfaceDiv`, `SurfaceSymGrad`
"""
SurfaceGrad(P) = OpApplied(P, SurfaceGradOp())

"""
    SurfaceDiv(P)

Surface divergence operator.

Computes the tangential divergence of vector fields defined on embedded
surfaces.

The operator uses local tangent bases evaluated at Gauss points.

See also: `SurfaceGrad`, `SurfaceSymGrad`
"""
SurfaceDiv(P) = OpApplied(P, SurfaceDivOp())

"""
    SurfaceSymGrad(P)

Symmetric tangential gradient operator for membrane mechanics.

Computes the membrane strain vector on a surface embedded in 3D.

The operator evaluates strains in a local tangent coordinate system
constructed at each Gauss point.

For vector fields, the operator returns:

    [ε11, ε22, γ12]

where `1` and `2` denote the local tangent directions.

Only in-surface membrane strains are included.
Bending strains are not part of this operator.

Typical usage:

```julia
K = ∫(SurfaceSymGrad(Pu) ⋅ C ⋅ SurfaceSymGrad(Pu))
```

See also: SurfaceGrad, SurfaceDiv
"""
SurfaceSymGrad(P) = OpApplied(P, SurfaceSymGradOp())

function _check_scalarfields(expr)

    if expr isa WeakTerm

        c = expr.term.coef

        if c isa ScalarField
            _check_scalarfield(c)

        elseif c isa AbstractMatrix
            for x in c
                if x isa ScalarField
                    _check_scalarfield(x)
                end
            end
        end

    elseif expr isa WeakSum

        _check_scalarfields(expr.a)
        _check_scalarfields(expr.b)

    end

end

# ------------------------------------------------------------
# Bilinear term
# ------------------------------------------------------------
"""
    BilinearTerm(a, coef, b)

Internal representation of a bilinear weak-form term.

Represents

    ∫ (a(v))' * coef * b(u)

where `a` is the test-side operator and `b` is the trial-side operator.
`coef` may be a scalar, a matrix, or a matrix chain containing `Number`
and `ScalarField` entries.
"""
struct BilinearTerm
    a::OpApplied
    coef
    b::OpApplied
    gauss::Any
    function BilinearTerm(a::OpApplied, coef, b::OpApplied)
        return new(a, coef, b, :full)
    end
end

"""
    LinearTerm(chain::MatrixChain)

Represents a linear weak-form term where the test-field operator
appears once and the other factors are known coefficient fields
or matrices.

Examples
--------
    f ⋅ Pu
    p ⋅ n ⋅ Pu
    σT ⋅ SymGrad(Pu)
    Pu ⋅ A ⋅ g

Used internally by the `∫` assembler.
"""
struct LinearTerm
    chain
end

# check

@inline function _check_coeff_matrix(C)

    for x in C
        if !(x isa Number || x isa ScalarField)
            error("Matrix coefficient entries must be Number or ScalarField, got $(typeof(x))")
        end
    end

end

# ------------------------------------------------------------
# Operator combination
# ------------------------------------------------------------

function _matvec_sf(A, v::AbstractVector)

    m, n = size(A)
    @assert n == length(v)

    w = Vector{Any}(undef, m)

    for i in 1:m
        s = nothing
        for j in 1:n
            term = A[i,j] * v[j]
            s = s === nothing ? term : s + term
        end
        w[i] = s
    end

    return w
end

_to_components(v::Number) = [v]

_to_components(v::AbstractVector) = v

_to_components(v::AbstractMatrix) = v

_to_components(v::ScalarField) = [v]

_to_components(v::VectorField) = [v[i] for i in 1:(v.type == :v2D ? 2 : v.type == :v3D ? 3 : error("_to_components: wrong vector type: $(v.type)"))]

function _to_components(T::TensorField)
    if T.model.dim == 2
        return [T[1], T[2], T[4], T[5]]   # 2×2 → 4 komponens
    else
        return [T[i] for i in 1:9]        # 3×3 → 9 komponens
    end
end

"""
    tensorfield_to_matrix(F::TensorField)

Convert a TensorField to a matrix of ScalarField components.
"""
function tensorfield_to_matrix(F::TensorField)

    dim = F.model.pdim

    if dim == 2
        return [
            F[1,1]  F[1,2]
            F[2,1]  F[2,2]
        ]
    elseif dim == 3
        return [
            F[1,1]  F[1,2]  F[1,3]
            F[2,1]  F[2,2]  F[2,3]
            F[3,1]  F[3,2]  F[3,3]
        ]
    end
end

# Collapse operator chain A⋅B⋅C⋅v → A(B(C(v)))
function collapse_chain(mats, v)

    w = _to_components(v)

    for i in length(mats):-1:1
        w = _matvec_sf(mats[i], w)
    end

    return w
end

function chain_dims(mats)
    M = mats[1]
    for i in 2:length(mats)
        M = M * mats[i]
    end
    return size(M)
end

"""
    MatrixChain

Internal representation of chained tensor coefficients.

Constructed automatically by expressions such as

    Grad(Pu) ⋅ A ⋅ B ⋅ C ⋅ Grad(Pu)
    Pu ⋅ A ⋅ B ⋅ g

The matrices are stored in `mats` and multiplied during
assembly at Gauss points.

The chain is later collapsed during assembly.
Users normally never construct this type directly.
"""
struct MatrixChain
    a::OpApplied
    mats::Vector{Any}
end

###################################################################
# Weak-form dot operator
###################################################################

"""
Weak-form inner product operator used in the DSL.

General pattern

    P1 ⋅ M1 ⋅ M2 ⋅ ... ⋅ Mn ⋅ P2

where

    P1,P2 : OpApplied
    Mi    : matrices or scalar coefficients

If the chain ends with an operator → BilinearTerm

    Grad(Pu) ⋅ C ⋅ Grad(Pu)

If the chain ends with a field → LinearTerm

    Pu ⋅ g
"""
###################################################################
# Operator – Operator  → bilinear
###################################################################

⋅(a::OpApplied, b::OpApplied) =
    BilinearTerm(a, 1.0, b)


###################################################################
# Operator – Matrix  → start matrix chain
###################################################################

function ⋅(a::OpApplied, C::AbstractMatrix)
    _check_coeff_matrix(C)
    return MatrixChain(a, Any[C])
end

⋅(a::OpApplied, F::TensorField) =
    MatrixChain(a, Any[tensorfield_to_matrix(F)])

⋅(mc::MatrixChain, F::TensorField) =
    MatrixChain(mc.a, Any[mc.mats..., tensorfield_to_matrix(F)])

###################################################################
# Operator – scalar coefficient
###################################################################

⋅(a::OpApplied, c::Union{Number,ScalarField}) =
    MatrixChain(a, Any[c])


###################################################################
# Continue matrix chain
###################################################################

function ⋅(mc::MatrixChain, C::AbstractMatrix)
    _check_coeff_matrix(C)
    push!(mc.mats, C)
    return mc
end


###################################################################
# Continue scalar coefficient chain
###################################################################

function ⋅(mc::MatrixChain, c::Union{Number,ScalarField})
    push!(mc.mats, c)
    return mc
end


###################################################################
# Chain closes with operator → bilinear
###################################################################

⋅(mc::MatrixChain, b::OpApplied) =
    BilinearTerm(mc.a, mc.mats, b)


###################################################################
# Chain closes with field → linear
###################################################################

⋅(mc::MatrixChain, g::Union{
        Number,
        ScalarField,
        VectorField,
        TensorField,
        AbstractVector
    }) =
    LinearTerm(MatrixChain(mc.a, Any[mc.mats..., g]))

⋅(a::OpApplied, g::Union{
        Number,
        ScalarField,
        VectorField,
        TensorField,
        AbstractVector
    }) =
    LinearTerm(MatrixChain(a, [g])) #maybe Any[g]

###################################################################
# DSL sugar
###################################################################

⋅(P::Problem, x) = Id(P) ⋅ x
⋅(x, P::Problem) = x ⋅ Id(P)
⋅(P1::Problem, P2::Problem) = Id(P1) ⋅ Id(P2)

# ------------------------------------------------------------
# Weak expression tree
# ------------------------------------------------------------

"""
Abstract type representing a weak-form expression.

Expressions are built using operator combinations and
arithmetic operations and are later assembled using `∫`.

Example

    expr = Grad(Pu) ⋅ Grad(Pu)
    K = ∫(expr; Ω="domain")
"""
abstract type WeakExpr end

"""
Single bilinear term in a weak-form expression.

Represents

    coef * (Op_s v ⋅ Op_u u)

Terms are automatically constructed when combining
operators using the DSL.
"""
struct WeakTerm{T} <: WeakExpr
    coef::Number
    term::T
end

"""
Sum of two weak-form expressions.

Used internally to represent expressions like

    a + b

inside weak forms.
"""
struct WeakSum <: WeakExpr
    a::WeakExpr
    b::WeakExpr
end


# ------------------------------------------------------------
# Expression building
# ------------------------------------------------------------

import Base: +, -, *, adjoint, transpose

+(a::WeakExpr, b::WeakExpr) = WeakSum(a, b)

+(a::WeakTerm, b::WeakTerm) = WeakSum(a, b)

###############################################################
# Automatic promotion to WeakTerm
###############################################################

promote_term(t::BilinearTerm) = WeakTerm(1.0, t)
promote_term(t::LinearTerm)   = WeakTerm(1.0, t)
promote_term(t::WeakTerm) = t



###############################################################
# Scalar multiplication
###############################################################

#=
*(c::Union{Number,ScalarField}, t::BilinearTerm) =
    BilinearTerm(t.a, c, t.b)

*(t::BilinearTerm, c::Union{Number,ScalarField}) =
    BilinearTerm(t.a, c, t.b)
=#

*(op::OpApplied, c::Union{Number,ScalarField}) =
    MatrixChain(op, Any[c])

*(c::Union{Number,ScalarField}, op::OpApplied) =
    MatrixChain(op, Any[c])

###############################################################
# Addition
###############################################################

+(a::BilinearTerm, b::BilinearTerm) =
    WeakSum(promote_term(a), promote_term(b))

+(a::WeakTerm, b::BilinearTerm) =
    WeakSum(a, promote_term(b))

+(a::BilinearTerm, b::WeakTerm) =
    WeakSum(promote_term(a), b)


###############################################################
# Subtraction
###############################################################

-(a::WeakExpr, b::WeakExpr) = a + (-b)

-(t::WeakTerm) = WeakTerm(-t.coef, t.term)

-(t::BilinearTerm) = WeakTerm(-1.0, t)

-(a::WeakTerm, b::BilinearTerm) =
    WeakSum(a, -promote_term(b))

-(a::BilinearTerm, b::WeakTerm) =
    WeakSum(promote_term(a), -b)

-(a::BilinearTerm, b::BilinearTerm) =
    WeakSum(promote_term(a), -promote_term(b))


###################################################
# Domains
###################################################

function _domain_spec(; Ω=nothing, Γ=nothing)

    if Ω === nothing && Γ === nothing
        return nothing
    elseif Ω !== nothing && Γ === nothing
        return DomainSpec(:Ω, String(Ω))
    elseif Γ !== nothing && Ω === nothing
        return DomainSpec(:Γ, String(Γ))
    else
        error("Use either Ω= or Γ= in ∫, not both.")
    end

end

function _check_domain_dim(problem, dom::DomainSpec)

    dimTags = gmsh.model.getEntitiesForPhysicalName(dom.name)

    isempty(dimTags) &&
        error("Physical group \"$(dom.name)\" not found.")

    for (dim, tag) in dimTags

        if dom.kind === :Ω
            dim == problem.dim ||
                error("Ω=\"$(dom.name)\" dim=$dim but problem.dim=$(problem.dim)")
        else
            dim < problem.dim ||
                error("Γ=\"$(dom.name)\" dim=$dim but expected < $(problem.dim)")
        end

    end

end

# ------------------------------------------------------------
# Internal assembly
# ------------------------------------------------------------

function _assemble(
    term::WeakTerm{BilinearTerm}, dom, weight=nothing, gauss=:full,
    threads=:auto, assembly::Symbol=:csc)

    Pu = term.term.a.P

    if dom !== nothing
        gmsh.model.setCurrent(Pu.name)
        _check_domain_dim(Pu, dom)
    end

    coef = term.term.coef

    K = assemble_operator(
        term.term.a.P,
        term.term.a.op,
        term.term.b.P,
        term.term.b.op,
        coefficient=coef, #isa AbstractMatrix ? 1.0 : coef,
        weight=weight, #coef isa AbstractMatrix ? coef : nothing,
        domain=dom,
        gauss=gauss,
        assembly = assembly === :csc ? :csc : :matrix,
        threads=threads
    )

    assembly === :csc && dropzeros!(K.A)
    return K

end

function _assemble(
    expr::WeakSum, dom, weight, gauss=:full,
    threads=:auto, assembly::Symbol=:csc)

    _assemble(expr.a, dom, weight, gauss, threads, assembly) +
        _assemble(expr.b, dom, weight, gauss, threads, assembly)

end

function _assemble(
    term::WeakTerm{LinearTerm}, dom, weight=nothing, gauss=:full,
    threads=:auto, assembly::Symbol=:csc)

    t = term.term
    chain = t.chain

    if chain isa MatrixChain
        op   = chain.a
        mats = chain.mats
        g    = mats[end]
    elseif chain isa Tuple && length(chain) == 2
        op, g = chain
    else
        op = chain
        g  = chain
    end

    P = op.P

    if dom !== nothing
        gmsh.model.setCurrent(P.name)
        _check_domain_dim(P, dom)
    end

    assemble_linear(
        op.P,
        op.op,
        g;
        weight = weight,
        domain = dom,
        gauss=gauss,
        threads=threads
    )

end

# ------------------------------------------------------------
# Integral
# ------------------------------------------------------------

"""
    ∫(expr::WeakExpr; Ω=nothing, Γ=nothing, weight=nothing, gauss=:full,
      threads=:auto, assembly=:csc)

Assemble a finite element operator from a weak-form expression.

Examples
--------

Diffusion

    K = ∫( Grad(Pu) ⋅ Grad(Pu); Ω="domain")

Elasticity

    K = ∫( SymGrad(Pu) ⋅ C ⋅ SymGrad(Pu); Ω="solid")

Tensor chain

    K = ∫( Grad(Pu) ⋅ F' ⋅ S ⋅ F ⋅ Grad(Pu); Ω="solid")

Elastic support

    K = ∫( Id(Pu) ⋅ [kx 0; 0 ky] ⋅ Id(Pu); Γ="boundary")

or

    K = ∫( Pu ⋅ [kx 0; 0 ky] ⋅ Pu; Γ="boundary")

Mixed formulation

    A = ∫( Div(Pu) ⋅ Pp )
    B = ∫( Pp ⋅ Div(Pu) )

Linear form

    f = ∫(Pu ⋅ g)

With operator chain

    f = ∫(Pu ⋅ A ⋅ g)

With coefficient

    f = ∫(PT ⋅ PT * h, Γ="right")

The bilinear- and linear-form assembly uses all available Julia threads by
default. Set `threads` to a positive integer to control the worker count.

Serial assembly can be selected for individual integrals with

    threads=1

For example

    K = ∫(SymGrad(Pu) ⋅ C ⋅ SymGrad(Pu);
          Ω="solid",
          threads=1)

Direct CSC matrix assembly is the default. The legacy IJV path remains
available explicitly with `assembly=:ijv`.

# Arguments

`expr`

Weak-form expression composed of operators and coefficients.

`gauss`

`:full` or `:reduced` integration, or an `Int` offset to the Gauss rule
order (not a direct point count). For polynomial order `p`, the standard
paths request `Gauss(2p+1)` for full and `Gauss(max(1,2p-1))` for reduced;
an integer `k` requests `Gauss(max(1,2p+1+k))`. `:gll` requests experimental
Gauss-Lobatto-Legendre quadrature and is currently supported only for a single
spectral bilinear term of the form `P ⋅ P`.

# Keyword arguments

`Ω`

Volume physical group name.

`Γ`

Boundary physical group name.

# Returns

A bilinear form returns a `SystemMatrix`; a linear form returns a nodal
`ScalarField`, `VectorField`, or `TensorField` of assembled load entries.
Mixed-field sums may return block `SystemMatrix`s or `SystemVector`s.
These linear-form entries are equivalent nodal loads, not point samples of the
integrand.

`Ω` and `Γ` are mutually exclusive physical-group selectors. With neither,
the field's material groups define the integration domain. `weight` is an
additional numeric or nodal scalar-field factor. Field-valued coefficients
must satisfy the nodal-storage requirements checked by the selected path;
convert explicitly with `elementsToNodes` when needed.

Single bilinear and compound forms additionally accept `csc_matrix`,
`element_chunk_size`, and, for a single bilinear term, `updateFrom`; these
keywords are not available for every weak sum or linear-form method.
`csc_matrix` accumulates into its existing values. `gauss=:auto` is used by
compound forms to honor `full(...)`/`reduced(...)` annotations.

This assembly API differs from `∫(P, group, f)` (a scalar spatial integral)
and `∫(field)` (a time integral). Use `integrate(P, group, f)` when a scalar
spatial result is intended.
"""
function ∫(
    expr::WeakExpr; Ω=nothing, Γ=nothing, weight=nothing, gauss=:full,
    threads=:auto, assembly::Symbol=:csc)

    gauss === :gll && error(
        "GLL quadrature currently supports only a single bilinear term of " *
        "the form P ⋅ P. Weak sums and other WeakExpr forms are not yet implemented."
    )

    _check_scalarfields(expr)

    assembly in (:csc, :ijv) || error(
        "Unsupported matrix assembly mode $assembly. Expected :csc or :ijv."
    )

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    return _assemble(expr, dom, weight, gauss, threads, assembly)

end

"""
    ∫(t::BilinearTerm; Ω=nothing, Γ=nothing, weight=nothing,
      gauss=:full, threads=:auto, assembly=:csc,
      csc_matrix=nothing, element_chunk_size=:auto, updateFrom=nothing)

Assemble one bilinear term.

Direct CSC assembly is the default. If `csc_matrix` is omitted, the structural
pattern is built once, filled, and finalized by removing numerical zeros. If a
prebuilt `SparseMatrixCSC{Float64,Int}` is supplied, its complete pattern is
preserved so that it can be reused. Assembly adds to its current `nzval`
contents; reset them with `fill!(csc_matrix.nzval, 0.0)` before starting an
independent assembly.

Use `assembly=:ijv` for the legacy triplet-based path. `csc_matrix` and
`element_chunk_size` apply only to direct CSC assembly.

`updateFrom` is forwarded only when supplied. It is intended for specialized
operators such as `ContactGap`, which can update their geometry from the given
field immediately before assembly.
"""
function ∫(t::BilinearTerm;
    Ω=nothing, Γ=nothing, weight=nothing,
    gauss=:full,
    threads=:auto,
    assembly::Symbol=:csc,
    csc_matrix=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto,
    updateFrom=nothing)

    assembly ∈ (:ijv, :csc) || error(
        "Unsupported matrix assembly mode $assembly. " *
        "Expected :csc or :ijv."
    )

    assembly === :ijv && csc_matrix !== nothing && error(
        "csc_matrix is only valid with assembly=:csc."
    )

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    Pu = t.a.P

    if dom !== nothing
        gmsh.model.setCurrent(Pu.name)
        _check_domain_dim(Pu, dom)
    end

    if gauss === :gll
        assembly === :csc || error(
            "GLL quadrature currently supports only assembly=:csc."
        )
        csc_matrix === nothing || error(
            "GLL quadrature does not yet support csc_matrix reuse."
        )
        updateFrom === nothing || error(
            "GLL quadrature does not support updateFrom."
        )
        element_chunk_size === :auto || error(
            "GLL quadrature does not yet support element_chunk_size."
        )
        t.a.P.basis === :spectral || error(
            "GLL quadrature requires basis=:spectral on the test field."
        )
        t.b.P.basis === :spectral || error(
            "GLL quadrature requires basis=:spectral on the trial field."
        )
        isdefined(@__MODULE__, :spectralIntegral) || error(
            "GLL quadrature requires spectral.jl to be loaded."
        )

        return spectralIntegral(
            t;
            domain=dom,
            weight=weight,
            threads=threads
        )
    end

    K = if updateFrom === nothing
        assemble_operator(
            t.b.P,
            t.b.op,
            t.a.P,
            t.a.op;
            coefficient = t.coef,
            domain = dom,
            weight = weight,
            gauss = gauss,
            assembly = assembly === :csc ? :csc : :matrix,
            threads=threads,
            K=csc_matrix,
            element_chunk_size=element_chunk_size
        )
    else
        assemble_operator(
            t.b.P,
            t.b.op,
            t.a.P,
            t.a.op;
            coefficient = t.coef,
            domain = dom,
            weight = weight,
            gauss = gauss,
            assembly = assembly === :csc ? :csc : :matrix,
            threads=threads,
            K=csc_matrix,
            element_chunk_size=element_chunk_size,
            updateFrom=updateFrom
        )
    end

    preserve_pattern =
        _preserve_csc_pattern(t.a.op) ||
        _preserve_csc_pattern(t.b.op)

    assembly === :csc &&
        csc_matrix === nothing &&
        !preserve_pattern &&
        dropzeros!(K.A)

    return K

end

"""
    ∫(a::OpApplied, b::OpApplied; kwargs...)

Assemble the bilinear form defined by two applied operators with an identity
coefficient. The CSC-related keywords and pattern-reuse rules are identical to
those of `∫(::BilinearTerm)`.
"""
function ∫(a::OpApplied, b::OpApplied;
    Ω=nothing, Γ=nothing, weight=nothing,
    gauss=:full,
    threads=:auto,
    assembly::Symbol=:csc,
    csc_matrix=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto,
    updateFrom=nothing)

    if gauss === :gll
        return ∫(
            BilinearTerm(a, 1.0, b);
            Ω=Ω,
            Γ=Γ,
            weight=weight,
            gauss=:gll,
            threads=threads,
            assembly=assembly,
            csc_matrix=csc_matrix,
            element_chunk_size=element_chunk_size,
            updateFrom=updateFrom
        )
    end

    assembly ∈ (:ijv, :csc) || error(
        "Unsupported matrix assembly mode $assembly. " *
        "Expected :csc or :ijv."
    )

    assembly === :ijv && csc_matrix !== nothing && error(
        "csc_matrix is only valid with assembly=:csc."
    )

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    Pu = a.P

    if dom !== nothing
        gmsh.model.setCurrent(Pu.name)
        _check_domain_dim(Pu, dom)
    end

    K = if updateFrom === nothing
        assemble_operator(
            b.P,
            b.op,
            a.P,
            a.op;
            coefficient = 1.0,
            domain = dom,
            weight = weight,
            gauss = gauss,
            assembly = assembly === :csc ? :csc : :matrix,
            threads=threads,
            K=csc_matrix,
            element_chunk_size=element_chunk_size
        )
    else
        assemble_operator(
            b.P,
            b.op,
            a.P,
            a.op;
            coefficient = 1.0,
            domain = dom,
            weight = weight,
            gauss = gauss,
            assembly = assembly === :csc ? :csc : :matrix,
            threads=threads,
            K=csc_matrix,
            element_chunk_size=element_chunk_size,
            updateFrom=updateFrom
        )
    end

    preserve_pattern =
        _preserve_csc_pattern(a.op) ||
        _preserve_csc_pattern(b.op)

    assembly === :csc &&
        csc_matrix === nothing &&
        !preserve_pattern &&
        dropzeros!(K.A)

    return K

end

function ∫(t::LinearTerm; Ω=nothing, Γ=nothing, weight=nothing, gauss=:full, threads=:auto)

    gauss === :gll && error(
        "GLL quadrature for linear forms is not yet implemented."
    )

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    mc = t.chain
    a  = mc.a

    mats = mc.mats

    rhs = mats[end]
    coeffs = mats[1:end-1]

    rhs = collapse_chain(coeffs, rhs)

    P  = a.P
    op = a.op

    if dom !== nothing
        gmsh.model.setCurrent(P.name)
        _check_domain_dim(P, dom)
    end

    return assemble_linear(
        P,
        op,
        rhs;
        domain = dom,
        weight = weight,
        gauss = gauss,
        threads=threads
    )

end

function ∫(mc::MatrixChain; Ω=nothing, Γ=nothing, weight=nothing, gauss=:full, threads=:auto)

    return ∫(LinearTerm(mc); Ω=Ω, Γ=Γ, weight=weight, gauss=gauss, threads=threads)

end

"""
    ∫Ω(name, expr; kwargs...)

Assemble `expr` on physical group `name`, equivalent to `∫(expr; Ω=name)`.
A bilinear expression returns a `SystemMatrix`; a linear expression returns
an assembled nodal load field. Additional keywords are forwarded to [`∫`](@ref),
including quadrature, weights, thread count, and CSC reuse where supported
by the selected expression method.
"""
∫Ω(name, expr; kwargs...) = ∫(expr; Ω=name, kwargs...)

"""
    ∫Γ(name, expr; kwargs...)

Assemble `expr` on physical group `name`, equivalent to `∫(expr; Γ=name)`.
A bilinear expression returns a `SystemMatrix`; a linear expression returns
an assembled nodal load field. Additional keywords are forwarded to [`∫`](@ref),
including quadrature, weights, thread count, and CSC reuse where supported
by the selected expression method.
"""
∫Γ(name, expr; kwargs...) = ∫(expr; Γ=name, kwargs...)

const ε = SymGrad

