"""
LowLevelFEM Multifield Operator Kernel

This module implements the low-level operator infrastructure used to
assemble finite element matrices from weak form expressions.

The file contains three main layers:

1. Operator kernel

   Defines differential operators (`GradOp`, `DivOp`, `SymGradOp`, …)
   and constructs element operator matrices through `build_B!`.

2. Matrix assembly

   The function

       assemble_operator

   performs element integration and builds the global sparse matrix.

3. Weak-form DSL

   A lightweight domain-specific language allows the construction
   of weak forms using mathematical notation.

Examples

Scalar diffusion

    ∫( Grad(Pu) ⋅ Grad(Pu); Ω="domain")

Linear elasticity

    ∫( SymGrad(Pu) ⋅ C ⋅ SymGrad(Pu); Ω="solid")

Matrix chain coefficients

    ∫( Grad(Pu) ⋅ A ⋅ B ⋅ C ⋅ Grad(Pu); Ω="solid")

where `A,B,C` may be matrices containing `Number` or `ScalarField`
entries.

Compound operator sums

    B = A ⋅ Grad(Pu) + G ⋅ Pu
    K = ∫(B' ⋅ D ⋅ B; Ω="body", weight=2π*r)
or
    K = ∫(B' ⋅ D ⋅ B * (2π*r); Ω="body")

The sum is expanded into ordinary bilinear terms and assembled using the
standard bilinear assembler. Each term may carry its own quadrature rule via
`full(...)` and `reduced(...)`.

The implementation is designed to be

- transparent
- extensible
- suitable for multifield problems
- compatible with Gmsh meshes.

See also:
`SystemMatrix`, `SystemVector`, `solveField`.
"""

export Grad
export Div
export Curl
export SymGrad
export Id
export TensorDiv
export Adv
export AxialGrad
export TangentialGrad
export SurfaceGrad
export SurfaceDiv
export SurfaceSymGrad
export full, reduced
export build_csc_pattern

export ∫
export ∫Ω
export ∫Γ

export solveField

export ε

export ⋅

export solveEigenFields
export consistentToLumped


"""
Abstract base type for finite element operators.

Operators act on fields defined by a `Problem` and are used to build
weak-form expressions.

Examples include

    IdOp       identity operator
    GradOp     gradient
    DivOp      divergence
    CurlOp     curl
    SymGradOp  symmetric gradient

Operators are applied to a `Problem` via helper constructors:

    Grad(P)
    Div(P)
    SymGrad(P)

and used inside weak forms such as

    ∫( Grad(Pu) ⋅ Grad(Pu) )
"""
abstract type AbstractOp end

"""
    IdOp()

Identity operator.

Represents the operator

    u

in a weak form expression.

Example

    ∫( Id(Pu) ⋅ Id(Pu) )

which corresponds to

    ∫ u · v dΩ
"""
struct IdOp <: AbstractOp end     # u

"""
    GradOp()

Gradient operator.

For scalar fields

    grad(u)

For vector fields

    ∇u

The output dimension depends on the field type:

Scalar field

    grad(u) → dim

Vector field

    ∇u → dim × pdim
"""
struct GradOp <: AbstractOp end     # ∇u

"""
    DivOp()

Divergence operator acting on a vector field.

Represents

    div(u)

The operator assumes that

    pdim == dim

Example

    ∫( Div(Pu) ⋅ Div(Pu) )
"""
struct DivOp <: AbstractOp end     # div u

"""
    CurlOp()

Curl operator acting on a vector field.

2D

    curl(u) = ∂u_y/∂x − ∂u_x/∂y

3D

    curl(u) = ∇ × u
"""
struct CurlOp <: AbstractOp end    # rot u

"""
    SymGradOp()

Symmetric gradient operator used in linear elasticity.

Represents

    ε(u) = (∇u + ∇uᵀ)/2

Engineering strain components are returned.

2D

    [ εxx
      εyy
      γxy ]

3D

    [ εxx
      εyy
      εzz
      γxy
      γxz
      γyz ]
"""
struct SymGradOp <: AbstractOp end

"""
    TensorDivOp()

Divergence of a second-order tensor field.

Represents

    div(σ)

where σ is a tensor field stored with

    pdim = dim^2

The output is a vector with `dim` components.
"""
struct TensorDivOp <: AbstractOp end

"""
    AdvOp(b)

Advection operator for scalar transport problems.

Represents

    b ⋅ ∇u

where `b` is the advection velocity vector.

# Arguments

`b`

Advection velocity stored as

    NTuple{3,Float64}
"""
struct AdvOp <: AbstractOp
    b::NTuple{3,Float64}
end

"""
    AxialGradOp()

Axial (directional) gradient operator.

For vector fields:
    AxialGrad(u) = t ⋅ ∇u

Returns a scalar field representing the derivative along the element axis.
"""
struct AxialGradOp <: AbstractOp end

#AxialGrad(P) = AxialGradOp()

struct TangentialGradOp <: AbstractOp end

"""
    TangentialGrad(P)

Create a weak-form DSL tangential gradient operator applied to `P`.

# Description
The tangential gradient operator computes the projection of the gradient of a vector field
onto the element axis:

    TangentialGrad(u) = (∇u) ⋅ t

where `t` is the tangent vector of the element in the physical domain.

For vector fields `u ∈ ℝᵈ`, the result is a vector field of dimension `d`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object representing `(∇u) ⋅ t`.

Directional derivative operator along the local element tangent direction.

Computes derivatives projected onto the local tangent vector of the element.

Typical applications include:

- beam and rod formulations,
- directional gradients,
- streamline-like operators,
- embedded directional mechanics.

See also: `AxialGrad`

# Example
```julia
Kg = ∫(TangentialGrad(Pu) ⋅ N0 ⋅ TangentialGrad(Pu); Ω="truss")
```

# Notes

- Unlike `AxialGrad`, which produces a scalar strain measure `t ⋅ ∇u ⋅ t`,  
  `TangentialGrad` returns a vector quantity.

- This operator is useful for constructing geometric stiffness matrices  
  (initial stress effects) in truss and structural stability problems.  
"""
TangentialGrad(P) = OpApplied(P, TangentialGradOp())

"""
Surface gradient operator.

Computes tangential derivatives on embedded lower-dimensional entities,
typically a 2D surface embedded in 3D.
"""
struct SurfaceGradOp <: AbstractOp end

"""
Surface divergence operator.

Uses tangential derivatives on embedded lower-dimensional entities.
"""
struct SurfaceDivOp <: AbstractOp end

"""
Surface symmetric gradient operator.

This is the membrane strain operator for a displacement field living on a
surface embedded in the ambient space.
"""
struct SurfaceSymGradOp <: AbstractOp end

"""
    op_outdim(::IdOp, P::Problem)

Return the number of components produced by applying
operator `op` to a field described by `Problem` `P`.

This determines the number of rows of the element
operator matrix `B`.

Example

Scalar gradient in 2D

    op_outdim(GradOp(), P) = 2

Elastic strain in 3D

    op_outdim(SymGradOp(), P) = 6
"""
function op_outdim(::IdOp, P::Problem)
    return P.pdim
end

function op_outdim(::GradOp, P::Problem)
    # Scalar: grad -> dim
    if P.pdim == 1
        return P.dim
    end
    # Vector: grad(u) -> dim×pdim (full gradient, column per component)
    return P.dim * P.pdim
end

function op_outdim(::DivOp, P::Problem)
    # Vector: div -> 1 (assume pdim==dim)
    @assert P.pdim == P.dim "DivOp currently assumes vector field with pdim == dim."
    return 1
end

function op_outdim(::CurlOp, P::Problem)
    @assert P.pdim == P.dim "CurlOp requires vector field with pdim == dim."
    return (P.dim == 2) ? 1 : 3
end

function op_outdim(::SymGradOp, P::Problem)
    @assert P.pdim == P.dim

    if P.dim == 1
        return 1
    elseif P.dim == 2
        return 3
    else
        return 6
    end
end

#function op_outdim(::SymGradOp, P::Problem)
#    @assert P.pdim == P.dim "SymGradOp requires vector field with pdim == dim."
#    return (P.dim == 2) ? 3 : 6  # engineering strain components
#end

function op_outdim(::TensorDivOp, P::Problem)
    @assert P.pdim == P.dim^2 "TensorDivOp requires pdim == dim^2 (full 2nd-order tensor)."
    return P.dim
end

function op_outdim(op::AdvOp, P::Problem)
    @assert P.pdim == 1  # scalar field
    return 1
end

function op_outdim(::AxialGradOp, P::Problem)
    return 1
end

function op_outdim(::TangentialGradOp, P::Problem)
    return P.dim
end

function op_outdim(::SurfaceGradOp, P::Problem)
    return op_outdim(GradOp(), P)
end

function op_outdim(::SurfaceDivOp, P::Problem)
    return op_outdim(DivOp(), P)
end

#function op_outdim(::SurfaceSymGradOp, P::Problem)
#    return op_outdim(SymGradOp(), P)
#end
function op_outdim(::SurfaceSymGradOp, P::Problem)
    @assert P.pdim == 3 "SurfaceSymGradOp currently expects a 3D displacement field."
    return 3
end

"""
    build_B!(B::AbstractMatrix, ::IdOp, P::Problem, k::Int, h, ∂h, numNodes::Int)

Construct the operator matrix `B` at Gauss point `k`.

The matrix maps element degrees of freedom to the operator
value at the Gauss point.

# Arguments

`B`

Operator matrix of size

    (op_outdim(op,P), P.pdim*numNodes)

`op`

Finite element operator.

`P`

Problem describing the field.

`k`

Gauss point index.

`h`

Shape functions evaluated at Gauss points.

`∂h`

Physical gradients of shape functions.

`numNodes`

Number of nodes per element.

# Notes

This function is called inside the element integration
loop of `assemble_operator`.
"""
function build_B!(B::AbstractMatrix, ::IdOp, P::Problem, k::Int, h, ∂h, numNodes::Integer)
    fill!(B, 0.0)
    pdim = P.pdim
    @inbounds for a in 1:numNodes
        Na = h[(k-1)*numNodes+a]
        @inbounds for c in 1:pdim
            row = c
            col = (a - 1) * pdim + c
            B[row, col] = Na
        end
    end
    return B
end

function build_B!(B::AbstractMatrix, ::GradOp, P::Problem, k::Int, h, ∂h, numNodes::Integer)
    fill!(B, 0.0)
    pdim = P.pdim
    dim = P.dim

    if pdim == 1
        # grad(scalar): rows = dim
        @inbounds for a in 1:numNodes
            col = a # scalar: (a-1)*1 + 1
            @inbounds for d in 1:dim
                row = d
                B[row, col] = ∂h[d, (k-1)*numNodes+a]
            end
        end
    else
        # grad(vector): rows = dim*pdim
        # ordering rows: (comp-1)*dim + d   (component-major)
        @inbounds for a in 1:numNodes
            @inbounds for c in 1:pdim
                col = (a - 1) * pdim + c
                @inbounds for d in 1:dim
                    row = (c - 1) * dim + d
                    B[row, col] = ∂h[d, (k-1)*numNodes+a]
                end
            end
        end
    end
    return B
end

function build_B!(B::AbstractMatrix, ::DivOp, P::Problem, k::Int, h, ∂h, numNodes::Integer)
    fill!(B, 0.0)
    pdim = P.pdim
    dim = P.dim
    @assert pdim == dim

    # div(u) = Σ_i ∂u_i/∂x_i
    # For basis (node a, component i): contribution is ∂N_a/∂x_i
    @inbounds for a in 1:numNodes
        @inbounds for i in 1:dim
            col = (a - 1) * pdim + i
            B[1, col] += ∂h[i, (k-1)*numNodes+a]
        end
    end
    return B
end

function build_B!(B::AbstractMatrix, ::CurlOp,
    P::Problem, k::Int, h, ∂h, numNodes::Integer)
    fill!(B, 0.0)
    dim = P.dim
    pdim = P.pdim
    @assert pdim == dim

    if dim == 2
        # curl(u) = ∂uy/∂x - ∂ux/∂y  (scalar)
        @inbounds for a in 1:numNodes
            colx = (a - 1) * pdim + 1
            coly = (a - 1) * pdim + 2
            B[1, colx] = -∂h[2, (k-1)*numNodes+a]   # -∂N/∂y * ux_a
            B[1, coly] = ∂h[1, (k-1)*numNodes+a]   #  ∂N/∂x * uy_a
        end
    else
        # 3D curl:
        # cx = ∂uz/∂y - ∂uy/∂z
        # cy = ∂ux/∂z - ∂uz/∂x
        # cz = ∂uy/∂x - ∂ux/∂y
        @inbounds for a in 1:numNodes
            colx = (a - 1) * pdim + 1
            coly = (a - 1) * pdim + 2
            colz = (a - 1) * pdim + 3

            dNx = ∂h[1, (k-1)*numNodes+a]
            dNy = ∂h[2, (k-1)*numNodes+a]
            dNz = ∂h[3, (k-1)*numNodes+a]

            # cx row = 1
            B[1, coly] = -dNz
            B[1, colz] = dNy

            # cy row = 2
            B[2, colx] = dNz
            B[2, colz] = -dNx

            # cz row = 3
            B[3, colx] = -dNy
            B[3, coly] = dNx
        end
    end

    return B
end

function build_B!(B::AbstractMatrix, ::SymGradOp,
    P::Problem, k::Int, h, ∂h, numNodes::Integer)
    fill!(B, 0.0)
    dim = P.dim
    pdim = P.pdim
    @assert pdim == dim

    if dim == 1
        # ε = du/dx
        @inbounds for a in 1:numNodes
            col = a  # scalar field
            B[1, col] = ∂h[1, (k-1)*numNodes + a]
        end

    elseif dim == 2
        # rows: [εxx, εyy, γxy]
        @inbounds for a in 1:numNodes
            colx = (a - 1) * pdim + 1
            coly = (a - 1) * pdim + 2
            dNx = ∂h[1, (k-1)*numNodes+a]
            dNy = ∂h[2, (k-1)*numNodes+a]

            B[1, colx] = dNx          # εxx
            B[2, coly] = dNy          # εyy
            B[3, colx] = dNy          # γxy = dux/dy + duy/dx
            B[3, coly] = dNx
        end
    else
        # rows: [εxx, εyy, εzz, γxy, γxz, γyz]
        @inbounds for a in 1:numNodes
            colx = (a - 1) * pdim + 1
            coly = (a - 1) * pdim + 2
            colz = (a - 1) * pdim + 3

            dNx = ∂h[1, (k-1)*numNodes+a]
            dNy = ∂h[2, (k-1)*numNodes+a]
            dNz = ∂h[3, (k-1)*numNodes+a]

            B[1, colx] = dNx          # εxx
            B[2, coly] = dNy          # εyy
            B[3, colz] = dNz          # εzz

            B[4, colx] = dNy          # γxy
            B[4, coly] = dNx

            B[5, colx] = dNz          # γxz
            B[5, colz] = dNx

            B[6, coly] = dNz          # γyz
            B[6, colz] = dNy
        end
    end

    return B
end

function build_B!(B::AbstractMatrix, ::TensorDivOp,
    P::Problem, k::Int, h, ∂h, numNodes::Integer)
    fill!(B, 0.0)
    dim = P.dim
    pdim = P.pdim  # dim^2

    # assume tensor component ordering at node:
    # σ_ij mapped to α = (i-1)*dim + j   (i=row, j=col)
    @inbounds for a in 1:numNodes
        for i in 1:dim
            row = i
            for j in 1:dim
                col = (a - 1) * pdim + (i - 1) * dim + j
                B[row, col] = ∂h[j, (k-1)*numNodes+a]  # ∂/∂x_j
            end
        end
    end
    return B
end

function build_B!(B::AbstractMatrix, op::AdvOp,
    P::Problem, k::Int, h, ∂h, numNodes::Integer)

    fill!(B, 0.0)
    b = op.b
    dim = P.dim

    @inbounds for a in 1:numNodes
        val = 0.0
        for d in 1:dim
            val += b[d] * ∂h[d, (k-1)*numNodes+a]
        end
        B[1, a] = val
    end

    return B
end

function get_tangent(P::Problem, k::Int)
    if P.dim == 2
        return (1.0, 0.0)
    elseif P.dim == 3
        return (1.0, 0.0, 0.0)
    else
        return (1.0,)
    end
end

function compute_tangent(Jk, dim)
    t = Jk[1:dim,1]
    return t / norm(t)
end

function build_B!(B::AbstractMatrix, ::AxialGradOp,
                  P::Problem, k::Int, h, ∂h, numNodes::Integer, t)

    fill!(B, 0.0)

    pdim = P.pdim   # 2
    dim  = P.dim    # 2

    @inbounds for a in 1:numNodes
        for c in 1:pdim
            col = (a - 1)*pdim + c

            val = 0.0

            @inbounds for d in 1:dim
                # 🔥 EZ A HIÁNYZÓ RÉSZ
                val += t[c] * t[d] * ∂h[d, (k-1)*numNodes + a]
            end

            B[1, col] = val
        end
    end

    return B
end

function build_B!(B::AbstractMatrix, ::TangentialGradOp,
                  P::Problem, k::Int, h, ∂h, numNodes::Integer, t)

    fill!(B, 0.0)

    pdim = P.pdim   # = dim
    dim  = P.dim

    # rows = dim
    # cols = pdim*numNodes

    @inbounds for a in 1:numNodes
        for c in 1:pdim
            col = (a - 1)*pdim + c

            @inbounds for i in 1:dim
                row = i

                # (∇u)t komponens i-re
                B[row, col] = ∂h[i, (k-1)*numNodes + a] * t[c]
            end
        end
    end

    return B
end

function build_B!(B::AbstractMatrix, ::SurfaceGradOp,
    P::Problem, k::Int, h, ∂h, numNodes::Integer)

    return build_B!(B, GradOp(), P, k, h, ∂h, numNodes)
end

function build_B!(B::AbstractMatrix, ::SurfaceDivOp,
    P::Problem, k::Int, h, ∂h, numNodes::Integer)

    return build_B!(B, DivOp(), P, k, h, ∂h, numNodes)
end

function build_B!(B::AbstractMatrix, ::SurfaceSymGradOp,
    P::Problem, k::Int, h, ∂h, numNodes::Integer)

    fill!(B, 0.0)
    @assert P.pdim == 3 "SurfaceSymGradOp currently expects a 3D displacement field."

    @inbounds for a in 1:numNodes
        colx = (a - 1) * 3 + 1
        coly = (a - 1) * 3 + 2
        # colz = (a - 1) * 3 + 3

        dNx = ∂h[1, (k-1)*numNodes+a]
        dNy = ∂h[2, (k-1)*numNodes+a]

        # rows: [εxx, εyy, γxy]
        B[1, colx] = dNx
        B[2, coly] = dNy
        B[3, colx] = dNy
        B[3, coly] = dNx
    end

    return B
end

function build_B!(B::AbstractMatrix, ::SurfaceSymGradOp,
                  P::Problem, k::Int, h, ∂h, numNodes::Integer, t1, t2)

    fill!(B, 0.0)

    @assert P.dim == 3
    @assert P.pdim == 3

    @inbounds for a in 1:numNodes
        g = SVector(
            ∂h[1, (k-1)*numNodes+a],
            ∂h[2, (k-1)*numNodes+a],
            ∂h[3, (k-1)*numNodes+a]
        )

        dN1 = dot(g, t1)
        dN2 = dot(g, t2)

        colx = (a - 1)*3 + 1
        coly = (a - 1)*3 + 2
        colz = (a - 1)*3 + 3

        # ε11 = ∂(u⋅t1)/∂s1
        B[1, colx] = dN1 * t1[1]
        B[1, coly] = dN1 * t1[2]
        B[1, colz] = dN1 * t1[3]

        # ε22 = ∂(u⋅t2)/∂s2
        B[2, colx] = dN2 * t2[1]
        B[2, coly] = dN2 * t2[2]
        B[2, colz] = dN2 * t2[3]

        # γ12 = ∂(u⋅t1)/∂s2 + ∂(u⋅t2)/∂s1
        B[3, colx] = dN2 * t1[1] + dN1 * t2[1]
        B[3, coly] = dN2 * t1[2] + dN1 * t2[2]
        B[3, colz] = dN2 * t1[3] + dN1 * t2[3]
    end

    return B
end

struct DomainSpec
    kind::Symbol     # :Ω vagy :Γ
    name::String
end

Base.show(io::IO, d::DomainSpec) =
    print(io, "$(d.kind)=\"$(d.name)\"")

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
    pa::Dict{<:Integer,<:AbstractMatrix},
    elem::Integer,
    hcol::AbstractVector)

    vals = get(pa, Int(elem), nothing)

    vals === nothing && return 0.0

    return dot(view(vals, :, 1), hcol)
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

        for I in CartesianIndices(Cprep)

            cij = Cprep[I]

            if cij isa Number
                W[I] = cij
            else
                W[I] = _coeff_at_gp(cij, elem, hgp)
            end

        end

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

mutable struct OperatorWorkspace
    Jac::Matrix{Float64}
    jacDet::Vector{Float64}
    invJac::Matrix{Float64}
    ∂h::Matrix{Float64}
    Bu::Matrix{Float64}
    Bs::Matrix{Float64}
    Ke::Matrix{Float64}
    tmp::Matrix{Float64}

    patBu::MatrixPattern
    patBsT::MatrixPattern
    patTmp::MatrixPattern
    patCgp::Union{Nothing,MatrixPattern}
end

"""
    OperatorWorkspace(
        dim,
        num_nodes,
        num_integration_points,
        out_u,
        out_s,
        ndofs_u,
        ndofs_s
    )

Allocate reusable element-level work arrays for operator assembly.
"""
function OperatorWorkspace(
    dim::Int,
    num_nodes::Int,
    num_integration_points::Int,
    out_u::Int,
    out_s::Int,
    ndofs_u::Int,
    ndofs_s::Int
)
    Bu = zeros(out_u, ndofs_u)
    Bs = zeros(out_s, ndofs_s)
    tmp = zeros(ndofs_s, out_u)

    patBu  = MatrixPattern(out_u, ndofs_u)
    patBsT = MatrixPattern(ndofs_s, out_s)
    patTmp = MatrixPattern(ndofs_s, out_u)

    sizehint!(patBu.rows, length(Bu))
    sizehint!(patBu.cols, length(Bu))

    sizehint!(patBsT.rows, length(Bs))
    sizehint!(patBsT.cols, length(Bs))

    sizehint!(patTmp.rows, length(tmp))
    sizehint!(patTmp.cols, length(tmp))

    return OperatorWorkspace(
        zeros(3, 3 * num_integration_points),
        zeros(num_integration_points),
        zeros(3, 3 * num_integration_points),
        zeros(dim, num_nodes * num_integration_points),
        Bu,
        Bs,
        zeros(ndofs_s, ndofs_u),
        tmp,
        patBu,
        patBsT,
        patTmp,
        nothing
    )
end

"""
    dense_node_coordinates(problem)

Return a `3 × number_of_nodes` coordinate table indexed by node tag.

The common Gmsh case, in which node tags are already ordered as `1:n`, returns
a reshaped view of the coordinate vector without copying it. A single compact
reordering copy is used otherwise. LowLevelFEM currently requires dense node
tags because global degrees of freedom are derived directly from node tags.
"""
function dense_node_coordinates(problem::Problem)
    gmsh.model.setCurrent(problem.name)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()

    number_of_nodes, remainder = divrem(ndofs(problem), problem.pdim)
    remainder == 0 || error("dense_node_coordinates: invalid number of degrees of freedom.")
    length(node_tags) == number_of_nodes || error(
        "dense_node_coordinates: the mesh contains $(length(node_tags)) nodes, " *
        "but the problem expects $number_of_nodes."
    )

    ordered = true
    @inbounds for i in eachindex(node_tags)
        if Int(node_tags[i]) != i
            ordered = false
            break
        end
    end

    if ordered
        return reshape(coordinates, 3, number_of_nodes)
    end

    dense = Matrix{Float64}(undef, 3, number_of_nodes)
    @inbounds for source in eachindex(node_tags)
        node = Int(node_tags[source])
        1 <= node <= number_of_nodes || error(
            "dense_node_coordinates: node tag $node is outside 1:$number_of_nodes."
        )
        source3 = 3 * (source - 1)
        dense[1, node] = coordinates[source3 + 1]
        dense[2, node] = coordinates[source3 + 2]
        dense[3, node] = coordinates[source3 + 3]
    end

    return dense
end

@inline function jacobian_smatrix(Jac::AbstractMatrix{<:Real}, first_col::Int)
    return SMatrix{3,3,Float64,9}((
        Jac[1, first_col],     Jac[2, first_col],     Jac[3, first_col],
        Jac[1, first_col + 1], Jac[2, first_col + 1], Jac[3, first_col + 1],
        Jac[1, first_col + 2], Jac[2, first_col + 2], Jac[3, first_col + 2]
    ))
end

"""
    compute_element_geometry!(ws, coordinates, connectivity, first_node,
                              num_nodes, num_integration_points, dim, edim, grad_h)

Compute Jacobians and integration measures for one element directly from its
nodal coordinates. Only the reusable arrays in `ws` are written; no
element-sized heap allocation is required.
"""
function compute_element_geometry!(
    ws::OperatorWorkspace,
    coordinates::AbstractMatrix{<:Real},
    connectivity::AbstractVector{<:Integer},
    first_node::Integer,
    numNodes::Integer,
    numIntPoints::Integer,
    dim::Integer,
    edim::Integer,
    ∇h
)
    Jac = ws.Jac
    jacDet = ws.jacDet
    fill!(Jac, 0.0)

    @inbounds for k in 1:numIntPoints
        first_col = 3k - 2

        for a in 1:numNodes
            node = Int(connectivity[first_node + a])
            grad_row = 3a - 2

            for α in 1:edim
                dNa = ∇h[grad_row + α - 1, k]
                col = first_col + α - 1
                Jac[1, col] += coordinates[1, node] * dNa
                Jac[2, col] += coordinates[2, node] * dNa
                Jac[3, col] += coordinates[3, node] * dNa
            end
        end

        if edim == dim
            if dim == 1
                jacDet[k] = abs(Jac[1, first_col])
                Jac[2, first_col + 1] = 1.0
                Jac[3, first_col + 2] = 1.0
            elseif dim == 2
                jacDet[k] = abs(
                    Jac[1, first_col] * Jac[2, first_col + 1] -
                    Jac[1, first_col + 1] * Jac[2, first_col]
                )
                Jac[3, first_col + 2] = 1.0
            else
                Jk = jacobian_smatrix(Jac, first_col)
                jacDet[k] = abs(det(Jk))
            end
        elseif edim == 1
            measure2 = 0.0
            for i in 1:dim
                measure2 += Jac[i, first_col]^2
            end
            jacDet[k] = sqrt(measure2)
        elseif edim == 2
            g11 = 0.0
            g12 = 0.0
            g22 = 0.0
            for i in 1:dim
                a1 = Jac[i, first_col]
                a2 = Jac[i, first_col + 1]
                g11 += a1 * a1
                g12 += a1 * a2
                g22 += a2 * a2
            end
            jacDet[k] = sqrt(max(0.0, g11 * g22 - g12 * g12))
        else
            error(
                "compute_element_geometry!: unsupported entity dimension " *
                "$edim for problem dimension $dim."
            )
        end

        jacDet[k] > 0.0 || error(
            "compute_element_geometry!: singular geometry at integration point $k."
        )
    end

    return nothing
end

"""
    assemble_element_matrix!(ws, elem, gidx, ...)

Compute the local element matrix and store it in `ws.Ke`.
"""
function assemble_element_matrix!(
    ws::OperatorWorkspace,
    elem::Integer,
    gidx::Integer,
    jac_all,
    det_all,
    nip::Integer,
    numIntPoints::Integer,
    dim::Integer,
    edim::Integer,
    numNodes::Integer,
    ∇h,
    h,
    intWeights,
    Cprep,
    Wprep,
    op_u::AbstractOp,
    Pu::Problem,
    op_s::AbstractOp,
    Ps::Problem
)
    nip == numIntPoints || error(
        "assemble_element_matrix!: inconsistent number of integration points."
    )

    first_jac = gidx * 9 * nip + 1
    first_det = gidx * nip + 1

    copyto!(ws.Jac, 1, jac_all, first_jac, 9 * nip)
    copyto!(ws.jacDet, 1, det_all, first_det, nip)

    return assemble_element_matrix_from_geometry!(
        ws,
        elem,
        numIntPoints,
        dim,
        edim,
        numNodes,
        ∇h,
        h,
        intWeights,
        Cprep,
        Wprep,
        op_u,
        Pu,
        op_s,
        Ps
    )
end

"""
    assemble_element_matrix_from_geometry!(ws, elem, ...)

Compute the local element matrix using geometry already stored in `ws.Jac`
and `ws.jacDet`. This is shared by the legacy batched-Gmsh path and the
low-memory on-the-fly geometry path.
"""
function assemble_element_matrix_from_geometry!(
    ws::OperatorWorkspace,
    elem::Integer,
    numIntPoints::Integer,
    dim::Integer,
    edim::Integer,
    numNodes::Integer,
    ∇h,
    h,
    intWeights,
    Cprep,
    Wprep,
    op_u::AbstractOp,
    Pu::Problem,
    op_s::AbstractOp,
    Ps::Problem
)
    Jac    = ws.Jac
    jacDet = ws.jacDet
    invJac = ws.invJac
    ∂h     = ws.∂h
    Bu     = ws.Bu
    Bs     = ws.Bs
    Ke     = ws.Ke
    tmp    = ws.tmp

    @inbounds for k in 1:numIntPoints
        if edim == dim
            Jk = jacobian_smatrix(Jac, 3k - 2)
            invJac[1:3, 3k-2:3k] .= transpose(inv3x3_fast(Jk))
        end
    end

    # Physical / tangential gradients
    fill!(∂h, 0.0)

    @inbounds for k in 1:numIntPoints
        if edim == dim
            @views begin
                invJk = invJac[1:dim, 3k-2:3k-(3-dim)]

                for a in 1:numNodes
                    gha = ∇h[a*3-2:a*3-(3-dim), k]
                    dha = ∂h[1:dim, (k-1)*numNodes+a]
                    mul!(dha, invJk, gha)
                end
            end
        else
            Jk = @view Jac[1:dim, 3k-2:3k]
            Jtan = Matrix(Jk[:, 1:edim])
            Gtan = Jtan * inv(Jtan' * Jtan)

            for a in 1:numNodes
                gha = ∇h[a*3-2:a*3-3+edim, k]
                ∂h[1:dim, (k-1)*numNodes+a] .= Gtan * gha
            end
        end
    end

    fill!(Ke, 0.0)

    @inbounds for k in 1:numIntPoints
        if Cprep isa AbstractVector
            Cgp = _eval_coefficient_chain_at_gp(
                Cprep, elem, view(h, :, k)
            )
        else
            Cgp = _eval_coefficient_at_gp(
                Cprep, elem, view(h, :, k)
            )
        end

        wcoef =
            Wprep === nothing ?
            1.0 :
            _eval_coefficient_at_gp(Wprep, elem, view(h, :, k))

        if Cgp isa Number
            w = jacDet[k] * intWeights[k] * Cgp * wcoef

            if op_u isa SurfaceSymGradOp
                t1, t2 = surface_basis_from_J(Jac, k)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t1, t2)
            elseif op_u isa AxialGradOp
                t = @view Jac[1:Pu.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
            elseif op_u isa TangentialGradOp
                t = Jac[1:Pu.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
            else
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes)
            end

            if op_s isa SurfaceSymGradOp
                t1, t2 = surface_basis_from_J(Jac, k)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t1, t2)
            elseif op_s isa AxialGradOp
                t = @view Jac[1:Ps.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
            elseif op_s isa TangentialGradOp
                t = Jac[1:Ps.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
            else
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes)
            end

            detect_pattern!(ws.patBsT, transpose(Bs))
            detect_pattern!(ws.patBu, Bu)

            mul_opt!(
                Ke,
                transpose(Bs),
                Bu;
                patA = ws.patBsT,
                patB = ws.patBu,
                alpha = w,
                beta = 1.0
            )

        elseif Cgp isa AbstractMatrix
            w = jacDet[k] * intWeights[k] * wcoef

            if op_u isa SurfaceSymGradOp
                t1, t2 = surface_basis_from_J(Jac, k)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t1, t2)
            elseif op_u isa AxialGradOp
                t = @view Jac[1:Pu.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
            elseif op_u isa TangentialGradOp
                t = Jac[1:Pu.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
            else
                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes)
            end

            if op_s isa SurfaceSymGradOp
                t1, t2 = surface_basis_from_J(Jac, k)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t1, t2)
            elseif op_s isa AxialGradOp
                t = @view Jac[1:Ps.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
            elseif op_s isa TangentialGradOp
                t = Jac[1:Ps.dim, 3k-2]
                t = t / norm(t)
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
            else
                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes)
            end

            if ws.patCgp === nothing
                ws.patCgp =
                    MatrixPattern(size(Cgp, 1), size(Cgp, 2))

                sizehint!(ws.patCgp.rows, length(Cgp))
                sizehint!(ws.patCgp.cols, length(Cgp))
            end

            detect_pattern!(ws.patBsT, transpose(Bs))
            detect_pattern!(ws.patCgp, Cgp)

            mul_opt!(
                tmp,
                transpose(Bs),
                Cgp;
                patA = ws.patBsT,
                patB = ws.patCgp,
                alpha = 1.0,
                beta = 0.0
            )

            detect_pattern!(ws.patTmp, tmp)
            detect_pattern!(ws.patBu, Bu)

            mul_opt!(
                Ke,
                tmp,
                Bu;
                patA = ws.patTmp,
                patB = ws.patBu,
                alpha = w,
                beta = 1.0
            )
        else
            error("assemble_element_matrix!: unsupported coefficient type")
        end
    end

    return nothing
end

"""
    build_csc_pattern(Pu::Problem, Ps::Problem; domain=nothing, Ω=nothing, Γ=nothing)

Build the exact structural CSC pattern required by the selected mesh domain.

The pattern is constructed from node connectivity without allocating full
degree-of-freedom triplet arrays. The returned matrix has zero-valued `nzval`
entries and sorted row indices in every column.

The returned pattern is used by both serial and parallel direct CSC assembly.
"""
function build_csc_pattern(
    Pu::Problem,
    Ps::Problem;
    domain=nothing,
    Ω=nothing,
    Γ=nothing
)
    @assert Pu.name == Ps.name "Both problems must refer to the same gmsh model/mesh."
    @assert Pu.dim == Ps.dim "Both problems must have the same spatial dimension."

    gmsh.model.setCurrent(Pu.name)

    if Ω !== nothing || Γ !== nothing
        domain === nothing || error(
            "build_csc_pattern: use either domain=... or Ω/Γ, not both."
        )
        domain = _domain_spec(; Ω=Ω, Γ=Γ)
    end

    nrows = ndofs(Ps)
    ncols = ndofs(Pu)

    nrow_nodes, rem_s = divrem(nrows, Ps.pdim)
    ncol_nodes, rem_u = divrem(ncols, Pu.pdim)

    rem_s == 0 || error("build_csc_pattern: invalid number of test-space degrees of freedom.")
    rem_u == 0 || error("build_csc_pattern: invalid number of trial-space degrees of freedom.")

    # Collect node-level adjacency first. Expanding to components only after
    # duplicate removal keeps temporary storage much smaller than DoF-level I/J.
    row_nodes_by_col_node = [Int[] for _ in 1:ncol_nodes]

    phNames = String[]
    domkind = nothing

    if domain === nothing
        phNames = [mat.phName for mat in Pu.material]
    elseif domain isa DomainSpec
        phNames = [domain.name]
        domkind = domain.kind
    elseif domain isa AbstractString
        phNames = [String(domain)]
    else
        error(
            "build_csc_pattern: unsupported domain type $(typeof(domain)). " *
            "Expected nothing, DomainSpec or String."
        )
    end

    for phName in phNames
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        isempty(dimTags) && error("build_csc_pattern: physical group \"$phName\" not found.")

        for (edim, etag) in dimTags
            if domkind === :Ω
                edim == Pu.dim ||
                    error("Ω=\"$phName\" has dim=$edim but problem.dim=$(Pu.dim)")
            elseif domkind === :Γ
                edim < Pu.dim ||
                    error("Γ=\"$phName\" has dim=$edim but expected < $(Pu.dim)")
            end

            elemTypes = gmsh.model.mesh.getElementTypes(edim, etag)

            for et in elemTypes
                _, _, _, numNodes::Int64, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                _, connectivity =
                    gmsh.model.mesh.getElementsByType(et, etag)
                nel = length(connectivity) ÷ numNodes

                @inbounds for e in 1:nel
                    first_node = (e - 1) * numNodes

                    for b in 1:numNodes
                        col_node = Int(connectivity[first_node + b])

                        1 <= col_node <= ncol_nodes || error(
                            "build_csc_pattern: node tag $col_node is outside " *
                            "the trial-space node range 1:$ncol_nodes."
                        )

                        rows = row_nodes_by_col_node[col_node]

                        for a in 1:numNodes
                            row_node = Int(connectivity[first_node + a])

                            1 <= row_node <= nrow_nodes || error(
                                "build_csc_pattern: node tag $row_node is outside " *
                                "the test-space node range 1:$nrow_nodes."
                            )

                            push!(rows, row_node)
                        end
                    end
                end
            end
        end
    end

    # Sort and remove duplicates in place. Unlike unique!, this pass needs no
    # auxiliary hash set because the vectors are already sorted.
    total_node_nonzeros = 0

    for rows in row_nodes_by_col_node
        sort!(rows)

        if !isempty(rows)
            write_pos = 1
            previous = rows[1]

            @inbounds for read_pos in 2:length(rows)
                current = rows[read_pos]

                if current != previous
                    write_pos += 1
                    rows[write_pos] = current
                    previous = current
                end
            end

            resize!(rows, write_pos)
            sizehint!(rows, write_pos; shrink=true)
        end

        total_node_nonzeros += length(rows)
    end

    total_nonzeros = Base.Checked.checked_mul(
        Base.Checked.checked_mul(total_node_nonzeros, Ps.pdim),
        Pu.pdim
    )

    colptr = Vector{Int}(undef, ncols + 1)
    rowval = Vector{Int}(undef, total_nonzeros)

    position = 1

    @inbounds for col_node in 1:ncol_nodes
        rows = row_nodes_by_col_node[col_node]

        for comp_u in 1:Pu.pdim
            col = (col_node - 1) * Pu.pdim + comp_u
            colptr[col] = position

            for row_node in rows
                first_row = (row_node - 1) * Ps.pdim

                for comp_s in 1:Ps.pdim
                    rowval[position] = first_row + comp_s
                    position += 1
                end
            end
        end
    end

    colptr[ncols + 1] = position
    @assert position == total_nonzeros + 1

    return SparseMatrixCSC(
        nrows,
        ncols,
        colptr,
        rowval,
        zeros(Float64, total_nonzeros)
    )
end

"""
    find_csc_position(rowval, first, last, row)

Return the `nzval` position of `row` in a sorted CSC column range.

The search is allocation-free. An error indicates that the prebuilt sparsity
pattern and the element connectivity used for assembly are inconsistent.
"""
@inline function find_csc_position(
    rowval::Vector{Int},
    first::Int,
    last::Int,
    row::Int
)
    lo = first
    hi = last

    @inbounds while lo <= hi
        mid = lo + ((hi - lo) >>> 1)
        current = rowval[mid]

        if current < row
            lo = mid + 1
        elseif current > row
            hi = mid - 1
        else
            return mid
        end
    end

    error("find_csc_position: row $row is missing from the CSC pattern.")
end

"""
    resolve_num_threads(threads) -> Int

Resolve the requested number of assembly workers.

- `threads=:auto` uses all Julia threads in the default thread pool.
- `threads=1` selects serial assembly.
- A positive integer selects exactly that many worker tasks, up to
  `Threads.nthreads()`.
"""
function resolve_num_threads(threads)
    available = Threads.nthreads()

    if threads === :auto
        return available
    elseif threads isa Bool
        throw(ArgumentError(
            "`threads` must be `:auto` or a positive integer; " *
            "use `threads=1` instead of `false`."
        ))
    elseif threads isa Integer
        threads >= 1 || throw(ArgumentError(
            "`threads` must be `:auto` or a positive integer."
        ))
        threads <= available || throw(ArgumentError(
            "Requested threads=$threads, but Julia has only " *
            "$available thread(s) in the default thread pool."
        ))
        return Int(threads)
    else
        throw(ArgumentError(
            "`threads` must be `:auto` or a positive integer; " *
            "got $(repr(threads))."
        ))
    end
end

@inline function _worker_range(nitems::Int, nworkers::Int, worker::Int)
    base, extra = divrem(nitems, nworkers)
    count = base + (worker <= extra)
    first = (worker - 1) * base + min(worker - 1, extra) + 1
    return first:(first + count - 1)
end

function _run_workers(f, nworkers::Int)
    if nworkers == 1
        f(1)
        return nothing
    end

    @sync begin
        for worker in 1:nworkers
            let worker = worker
                Threads.@spawn f(worker)
            end
        end
    end
    return nothing
end

"""
    assemble_operator(Pu::Problem, 
                      op_u::AbstractOp, 
                      Ps::Problem, 
                      op_s::AbstractOp; 
                      coefficient::Union{Number,ScalarField}=1.0, 
                      weight=nothing, 
                      domain=nothing,
                      gauss=:full)

Assemble the sparse matrix corresponding to the bilinear form

    ∫ (Op_s v) ⋅ (Op_u u) dΩ

where

- `u` belongs to trial space `Pu`
- `v` belongs to test space `Ps`

# Arguments

`Pu`

Trial field problem.

`Ps`

Test field problem.

`op_u`

Operator applied to the trial field.

`op_s`

Operator applied to the test field.

# Keyword arguments

`coefficient`

Scalar or tensor coefficient multiplying the integrand.

Supported types

    Number
    ScalarField
    Matrix{Number}
    Matrix{ScalarField}
    Vector{Matrix}

The vector form represents a matrix chain

    A ⋅ B ⋅ C

which is evaluated at Gauss points as

    Cgp = A(x_gp) * B(x_gp) * C(x_gp)

`weight`

Optional scalar quadrature weight multiplying the integrand at each Gauss
point. It may be a `Number` or a `ScalarField`.

This is useful for geometry-related weights such as the axisymmetric factor
`2πr`:

    ∫(B' ⋅ D ⋅ B; Ω="body", weight=2π*r)

The weight is not treated as a matrix coefficient.

`domain`

Optional domain specification (`Ω` or `Γ`).

`gauss`

Number of integration points.
    `:full`       2order+1
    `:reduced`    2order-1
    Int           2order+1 + Int

`assembly=:csc`

Use direct CSC assembly. The CSC pattern is created first, element Jacobians
are computed on demand from nodal coordinates, and element matrices are
accumulated into one `nzval` array per worker. The first worker writes directly
to the final matrix; the remaining worker-local arrays are reduced into it.

`threads`

Number of assembly workers. `:auto` uses all Julia threads, while `1` selects
the serial low-memory path.

`element_chunk_size`

Element scheduling block size for `assembly=:csc`. `:auto` uses at most 4096
elements and adapts to the element and worker counts. Chunks are distributed
cyclically among the requested workers.

# Returns

`SystemMatrix`

Sparse matrix representing the discretized bilinear form.
"""
function assemble_operator(
    Pu::Problem, op_u::AbstractOp,
    Ps::Problem, op_s::AbstractOp;
    coefficient::Union{Number,ScalarField,AbstractMatrix,AbstractVector}=1.0,
    weight=nothing,
    domain=nothing,
    gauss=:full,
    assembly::Symbol=:ijv,
    threads=:auto,
    I=nothing,
    J=nothing,
    V=nothing,
    K=nothing,
    node_coordinates=nothing,
    _nzval_buffers=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto)

    direct_csc = assembly === :csc
    num_threads = resolve_num_threads(threads)
    use_threads = num_threads > 1

    old_blas_threads = LinearAlgebra.BLAS.get_num_threads()
    if use_threads && old_blas_threads != 1
        LinearAlgebra.BLAS.set_num_threads(1)
    end

    try
        assembly ∈ (:matrix, :triplets, :add, :csc) ||
        error("""
        assemble_operator: invalid assembly mode $assembly.
   
        Valid modes are:
            :matrix
            :triplets
            :add
            :csc
        """)
        @assert Pu.name == Ps.name "Both problems must refer to the same gmsh model/mesh."
        @assert Pu.dim == Ps.dim "Both problems must have the same spatial dimension."
        gmsh.model.setCurrent(Pu.name)
   
        # operator output dimensions
        out_u = op_outdim(op_u, Pu)
        out_s = op_outdim(op_s, Ps)
        
        if coefficient isa Number || coefficient isa ScalarField
            @assert out_u == out_s
        
        elseif coefficient isa AbstractMatrix
            if size(coefficient,1) != out_s || size(coefficient,2) != out_u
                error("""
                    Coefficient matrix size mismatch.
   
                    Test operator:
                        $(typeof(op_s))
                    output dimension = $out_s
   
                    Trial operator:
                        $(typeof(op_u))
                    output dimension = $out_u
   
                    Expected coefficient size:
                        ($out_s, $out_u)
   
                    Got:
                        $(size(coefficient))
   
                    Coefficient type:
                        $(typeof(coefficient))
                    """)
            end
        
        elseif coefficient isa AbstractVector
            if length(coefficient) == 1 &&
               (coefficient[1] isa Number || coefficient[1] isa ScalarField)
        
                @assert out_u == out_s
        
            else
                A1 = coefficient[1]
                An = coefficient[end]
        
                @assert A1 isa AbstractMatrix
                @assert An isa AbstractMatrix
        
                if size(A1,1) != out_s || size(An,2) != out_u
                    error("""
                    Coefficient chain size mismatch.
   
                        Test operator:
                            $(typeof(op_s))
                        output dimension = $out_s
   
                        Trial operator:
                            $(typeof(op_u))
                        output dimension = $out_u
   
                        Expected:
                            first matrix rows = $out_s
                            last matrix cols  = $out_u
   
                        Got:
                            size(A1) = $(size(A1))
                            size(An) = $(size(An))
   
                        Coefficient chain types:
                            $(typeof.(coefficient))
                        """)
                end
        
                for i in 1:length(coefficient)-1
                    Ai = coefficient[i]
                    Aj = coefficient[i+1]
        
                    @assert Ai isa AbstractMatrix
                    @assert Aj isa AbstractMatrix
                    @assert size(Ai,2) == size(Aj,1)
                end
            end
        end
   
        # ------------------------------------------------------------
        # weight dimension check
        # ------------------------------------------------------------
        if weight !== nothing
            # optional: csak Number vagy ScalarField legyen
            (weight isa Number || weight isa ScalarField) ||
            error("weight must be Number or ScalarField")
        end
   
        # ------------------------------------------------------------
        # prepare coefficient
        # ------------------------------------------------------------
   
        # ------------------------------------------------------------
        # Assembly storage
        # ------------------------------------------------------------
   
        build_pattern = assembly !== :add && !direct_csc

        Kcsc = nothing
        coordinates = nothing
        nzval_buffers = nothing

        if direct_csc
            I === nothing ||
                error("assemble_operator: I must be nothing for assembly=:csc.")
            J === nothing ||
                error("assemble_operator: J must be nothing for assembly=:csc.")
            V === nothing ||
                error("assemble_operator: V must be nothing for assembly=:csc.")

            if K === nothing
                Kcsc = build_csc_pattern(Pu, Ps; domain=domain)
                # The node-adjacency vectors used while building the pattern
                # are no longer needed. Reclaim them before loading geometry.
                GC.gc(true)
                
            else
                K isa SparseMatrixCSC{Float64,Int} || error(
                    "assemble_operator: assembly=:csc requires " *
                    "K::SparseMatrixCSC{Float64,Int}; got $(typeof(K))."
                )

                size(K) == (ndofs(Ps), ndofs(Pu)) || error(
                    "assemble_operator: supplied CSC matrix has size $(size(K)); " *
                    "expected $((ndofs(Ps), ndofs(Pu)))."
                )

                Kcsc = K
            end

            coordinates = node_coordinates === nothing ?
                dense_node_coordinates(Pu) : node_coordinates

            size(coordinates, 1) == 3 || error(
                "assemble_operator: node_coordinates must have three rows."
            )
            node_coordinates === nothing && GC.gc(true)

            if _nzval_buffers === nothing
                nzval_buffers = Vector{Vector{Float64}}(undef, num_threads)
                nzval_buffers[1] = Kcsc.nzval
                for worker in 2:num_threads
                    nzval_buffers[worker] = zeros(Float64, length(Kcsc.nzval))
                end
            else
                _nzval_buffers isa Vector{Vector{Float64}} || error(
                    "assemble_operator: invalid internal CSC accumulator storage."
                )
                length(_nzval_buffers) == num_threads || error(
                    "assemble_operator: CSC accumulator count does not match threads."
                )
                _nzval_buffers[1] === Kcsc.nzval || error(
                    "assemble_operator: the first CSC accumulator must be K.nzval."
                )
                all(length(buffer) == length(Kcsc.nzval) for buffer in _nzval_buffers) ||
                    error("assemble_operator: CSC accumulator length mismatch.")
                nzval_buffers = _nzval_buffers
                for worker in 2:num_threads
                    fill!(nzval_buffers[worker], 0.0)
                end
            end
        elseif build_pattern
            I === nothing ||
                error("assemble_operator: I must be nothing for assembly=$assembly.")
   
            J === nothing ||
                error("assemble_operator: J must be nothing for assembly=$assembly.")
   
            V === nothing ||
                error("assemble_operator: V must be nothing for assembly=$assembly.")
   
            lengthOfIJV =
                LowLevelFEM.estimateLengthOfIJV(Pu) *
                max(1, Ps.pdim) *
                max(1, Pu.pdim)
   
            I = Vector{Int}(undef, lengthOfIJV)
            J = Vector{Int}(undef, lengthOfIJV)
            V = Vector{Float64}(undef, lengthOfIJV)
   
        else
            I === nothing ||
                error("assemble_operator: I must be nothing for assembly=:add.")
   
            J === nothing ||
                error("assemble_operator: J must be nothing for assembly=:add.")
   
            V isa Vector{Float64} ||
                error(
                    "assemble_operator: assembly=:add requires " *
                    "V::Vector{Float64}; got $(typeof(V))."
                )
        end
   
        pos = 1
   
        dim = Pu.dim
   
        # ------------------------------------------------------------
        # Domain/material selection FIX
        # ------------------------------------------------------------
        phNames = String[]
        domkind = nothing
   
        if domain === nothing
            # original behavior: loop over all materials (each has phName)
            phNames = [mat.phName for mat in Pu.material]
        elseif domain isa DomainSpec
            phNames = [domain.name]
            domkind = domain.kind   # :Ω or :Γ
        elseif domain isa AbstractString
            # (optional) accept plain String domain too
            phNames = [String(domain)]
            domkind = nothing
        else
            error("assemble_operator: unsupported domain type $(typeof(domain)). Expected nothing, DomainSpec or String.")
        end
   
        # loop physical groups
        for phName in phNames
            dom_here =
                domain === nothing ? DomainSpec(:Ω, phName) :
                domain isa DomainSpec ? domain :
                DomainSpec(:Ω, String(phName))
   
            Cprep = _prepare_coefficient(coefficient, dom_here)
            Wprep = weight === nothing ? nothing : _prepare_coefficient(weight, dom_here)
   
            dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
   
            isempty(dimTags) && error("assemble_operator: physical group \"$phName\" not found.")
   
            for (edim, etag) in dimTags
   
                # optional dimension checks (matches your DSL design)
                if domkind === :Ω
                    edim == Pu.dim || error("Ω=\"$phName\" has dim=$edim but problem.dim=$(Pu.dim)")
                elseif domkind === :Γ
                    edim < Pu.dim || error("Γ=\"$phName\" has dim=$edim but expected < $(Pu.dim)")
                end
   
                if direct_csc
                    elemTypes = gmsh.model.mesh.getElementTypes(edim, etag)
                    elemTags = nothing
                    elemNodeTags = nothing
                else
                    elemTypes, elemTags, elemNodeTags =
                        gmsh.model.mesh.getElements(edim, etag)
                end
   
                for itype in eachindex(elemTypes)
                    et = elemTypes[itype]
                    _, _, order, numNodes::Int64, _, _ = gmsh.model.mesh.getElementProperties(et)

                    if direct_csc
                        element_tags, connectivity =
                            gmsh.model.mesh.getElementsByType(et, etag)
                    else
                        element_tags = elemTags[itype]
                        connectivity = elemNodeTags[itype]
                    end
   
                    gorder = gauss == :reduced ? max(1, 2order - 1) :
                        gauss == :full    ? (2order + 1) :
                        gauss isa Int     ? max(1, 2order+1 + gauss) :
                        (2order + 1)
   
                    intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(gorder))
                    #intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                    numIntPoints = length(intWeights)
   
                    _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                    h = reshape(fun, :, numIntPoints)
   
                    _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                    ∇h = reshape(dfun, :, numIntPoints)
   
                    if !direct_csc
                        # Legacy IJV path: keep the original batched Gmsh
                        # geometry retrieval until that path is retired.
                        elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)
                        jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)

                        numElementsGlobal = length(elemTags_global)
                        nip = length(det_all) ÷ numElementsGlobal
                        @assert nip == numIntPoints

                        TagType = eltype(elemTags_global)
                        idxmap = Dict{TagType,Int}()
                        sizehint!(idxmap, numElementsGlobal)

                        @inbounds for (gi, gtag) in enumerate(elemTags_global)
                            idxmap[gtag] = gi - 1
                        end
                    end
   
                    # buffers
                    nel = length(element_tags)
   
                    ndofs_u_loc = Pu.pdim * numNodes
                    ndofs_s_loc = Ps.pdim * numNodes
   
                    # The CSC path uses the flat Gmsh connectivity directly.
                    # The legacy IJV path retains its old dense table for now.
                    nnet = direct_csc ? nothing : zeros(Int, nel, numNodes)
                    if !direct_csc
                        @inbounds for e in 1:nel
                            for a in 1:numNodes
                                nnet[e, a] = connectivity[(e-1)*numNodes+a]
                            end
                        end
                    end
   
                    entries_per_element = ndofs_s_loc * ndofs_u_loc
                    block_start = pos
                    block_entries = nel * entries_per_element
   
                    if !direct_csc && build_pattern && block_start + block_entries - 1 > length(V)
                        newlen = max(
                            block_start + block_entries - 1,
                            Int(ceil(1.5 * length(V))) + 1000
                        )
   
                        resize!(I, newlen)
                        resize!(J, newlen)
                        resize!(V, newlen)
                    elseif !direct_csc && !build_pattern && block_start + block_entries - 1 > length(V)
                        error("""
                        assemble_operator: assembly pattern mismatch.
                    
                        The current element block requires more entries than the existing
                        V vector contains.
                        """)
                    end
   
                    workspaces = [
                        OperatorWorkspace(
                            dim,
                            numNodes,
                            numIntPoints,
                            out_u,
                            out_s,
                            ndofs_u_loc,
                            ndofs_s_loc
                        )
                        for _ in 1:num_threads
                    ]
                    
                    if direct_csc
                        colptr = Kcsc.colptr
                        rowval = Kcsc.rowval

                        chunk = element_chunk_size === :auto ?
                            max(1, min(4096, cld(nel, num_threads))) :
                            element_chunk_size isa Integer ? Int(element_chunk_size) :
                            error(
                                "assemble_operator: element_chunk_size must be " *
                                ":auto or a positive integer."
                            )
                        chunk > 0 || error(
                            "assemble_operator: element_chunk_size must be positive."
                        )

                        _run_workers(num_threads) do worker
                            ws = workspaces[worker]
                            nzval = nzval_buffers[worker]
                            stride = chunk * num_threads

                            for chunk_first in (1 + (worker - 1) * chunk):stride:nel
                                chunk_last = min(nel, chunk_first + chunk - 1)

                                @inbounds for e in chunk_first:chunk_last
                                    elem = element_tags[e]
                                    first_node = (e - 1) * numNodes

                                    compute_element_geometry!(
                                        ws,
                                        coordinates,
                                        connectivity,
                                        first_node,
                                        numNodes,
                                        numIntPoints,
                                        dim,
                                        edim,
                                        ∇h
                                    )

                                    assemble_element_matrix_from_geometry!(
                                        ws,
                                        elem,
                                        numIntPoints,
                                        dim,
                                        edim,
                                        numNodes,
                                        ∇h,
                                        h,
                                        intWeights,
                                        Cprep,
                                        Wprep,
                                        op_u,
                                        Pu,
                                        op_s,
                                        Ps
                                    )

                                    Ke = ws.Ke

                                # Traverse Ke and CSC column-major. Global DoFs
                                # are derived directly from flat connectivity,
                                # avoiding both nnet and gdof work arrays.
                                    for b_loc in 1:ndofs_u_loc
                                        node_b = div(b_loc - 1, Pu.pdim) + 1
                                        comp_b = mod(b_loc - 1, Pu.pdim) + 1
                                        Jb_node = Int(connectivity[first_node + node_b])
                                        Jb = (Jb_node - 1) * Pu.pdim + comp_b
                                        first = colptr[Jb]
                                        last = colptr[Jb + 1] - 1

                                        for a_loc in 1:ndofs_s_loc
                                            node_a = div(a_loc - 1, Ps.pdim) + 1
                                            comp_a = mod(a_loc - 1, Ps.pdim) + 1
                                            Ia_node = Int(connectivity[first_node + node_a])
                                            Ia = (Ia_node - 1) * Ps.pdim + comp_a
                                            p = find_csc_position(rowval, first, last, Ia)
                                            nzval[p] += Ke[a_loc, b_loc]
                                        end
                                    end
                                end
                            end
                        end
                    else
                        Vthread = V::Vector{Float64}
                        Ithread = build_pattern ? (I::Vector{Int}) : Int[]
                        Jthread = build_pattern ? (J::Vector{Int}) : Int[]
   
                        let Ithread = Ithread,
                        Jthread = Jthread,
                        Vthread = Vthread,
                        build_pattern = build_pattern
   
                        # element loop
                        # TODO: duplicated loop -- separate helper function is needed
                        if use_threads
                            _run_workers(num_threads) do worker
                                ws = workspaces[worker]
                                for e in _worker_range(nel, num_threads, worker)
                            
                                elem = elemTags[itype][e]
                                gidx = idxmap[elem]
                    
                                assemble_element_matrix!(
                                    ws,
                                    elem,
                                    gidx,
                                    jac_all,
                                    det_all,
                                    nip,
                                    numIntPoints,
                                    dim,
                                    edim,
                                    numNodes,
                                    ∇h,
                                    h,
                                    intWeights,
                                    Cprep,
                                    Wprep,
                                    op_u,
                                    Pu,
                                    op_s,
                                    Ps
                                )
                    
                                Ke = ws.Ke
                    
                                element_start = block_start + (e - 1) * entries_per_element
                            
                                # scatter Ke(s,u) -> global IJV
                                @inbounds for a_loc in 1:ndofs_s_loc
                                    node_a = div(a_loc - 1, Ps.pdim) + 1
                                    comp_a = mod(a_loc - 1, Ps.pdim) + 1
                                    Ia_node = nnet[e, node_a]
                                    Ia = (Ia_node - 1) * Ps.pdim + comp_a
                            
                                    row_start =
                                        element_start + (a_loc - 1) * ndofs_u_loc
                            
                                    @inbounds for b_loc in 1:ndofs_u_loc
                                        node_b = div(b_loc - 1, Pu.pdim) + 1
                                        comp_b = mod(b_loc - 1, Pu.pdim) + 1
                                        Jb_node = nnet[e, node_b]
                                        Jb = (Jb_node - 1) * Pu.pdim + comp_b
                            
                                        p = row_start + b_loc - 1
                            
                                        if build_pattern
                                            Ithread[p] = Ia
                                            Jthread[p] = Jb
                                            Vthread[p] = Ke[a_loc, b_loc]
                                        else
                                            Vthread[p] += Ke[a_loc, b_loc]
                                        end
                                    end
                                end
                                end
                            end
                        else
                            @inbounds for e in 1:nel
                                ws = workspaces[1]
                            
                                elem = elemTags[itype][e]
                                gidx = idxmap[elem]
                    
                                assemble_element_matrix!(
                                    ws,
                                    elem,
                                    gidx,
                                    jac_all,
                                    det_all,
                                    nip,
                                    numIntPoints,
                                    dim,
                                    edim,
                                    numNodes,
                                    ∇h,
                                    h,
                                    intWeights,
                                    Cprep,
                                    Wprep,
                                    op_u,
                                    Pu,
                                    op_s,
                                    Ps
                                )
                    
                                Ke = ws.Ke
                    
                                element_start = block_start + (e - 1) * entries_per_element
                            
                                # scatter Ke(s,u) -> global IJV
                                @inbounds for a_loc in 1:ndofs_s_loc
                                    node_a = div(a_loc - 1, Ps.pdim) + 1
                                    comp_a = mod(a_loc - 1, Ps.pdim) + 1
                                    Ia_node = nnet[e, node_a]
                                    Ia = (Ia_node - 1) * Ps.pdim + comp_a
                            
                                    row_start =
                                        element_start + (a_loc - 1) * ndofs_u_loc
                            
                                    @inbounds for b_loc in 1:ndofs_u_loc
                                        node_b = div(b_loc - 1, Pu.pdim) + 1
                                        comp_b = mod(b_loc - 1, Pu.pdim) + 1
                                        Jb_node = nnet[e, node_b]
                                        Jb = (Jb_node - 1) * Pu.pdim + comp_b
                            
                                        p = row_start + b_loc - 1
                            
                                        if build_pattern
                                            Ithread[p] = Ia
                                            Jthread[p] = Jb
                                            Vthread[p] = Ke[a_loc, b_loc]
                                        else
                                            Vthread[p] += Ke[a_loc, b_loc]
                                        end
                                    end
                                end
                            end
                        end
                        end
                    end
                    pos += block_entries
                end
            end
        end
   
        nentries = pos - 1
        nrows = ndofs(Ps)
        ncols = ndofs(Pu)
        
        @assert nrows > 0
        @assert ncols > 0

        if direct_csc
            _run_workers(num_threads) do worker
                @inbounds for p in _worker_range(
                    length(Kcsc.nzval), num_threads, worker)

                    value = Kcsc.nzval[p]
                    for source in 2:num_threads
                        value += nzval_buffers[source][p]
                    end
                    Kcsc.nzval[p] = value
                end
            end
            return SystemMatrix(Kcsc, Pu, Ps)
        end
        
        if assembly === :add
            nentries == length(V) ||
                error("""
                    assemble_operator: assembly pattern mismatch.
        
                The existing V vector contains:
                    $(length(V)) entries
        
                The current operator generated:
                    $nentries entries
        
                All compound terms must traverse exactly the same
                elements and local matrix entries in the same order.
                """)
        
            return V
        end
        
        resize!(I, nentries)
        resize!(J, nentries)
        resize!(V, nentries)
        
        if nentries > 0
            @assert maximum(I) <= nrows
            @assert maximum(J) <= ncols
            @assert minimum(I) >= 1
            @assert minimum(J) >= 1
        end
        
        if assembly === :triplets
            return I, J, V, nrows, ncols
        end
        
        K = sparse(I, J, V, nrows, ncols)
        dropzeros!(K)
   
        return SystemMatrix(K, Pu, Ps)
    finally
        if use_threads && old_blas_threads != 1
            LinearAlgebra.BLAS.set_num_threads(old_blas_threads)
        end
    end
end

"""
    matmul_sf(A, B)

Matrix–matrix or matrix–field multiplication used internally
when collapsing operator chains.

Supports combinations of:
- matrices
- scalar/vector/tensor fields
- numeric arrays

Returns the resulting matrix or field.
"""
function matmul_sf(A, B)

    m,k = size(A)
    k2,n = size(B)

    @assert k==k2

    C = Matrix{Any}(undef,m,n)

    for i in 1:m, j in 1:n

        s = nothing

        for p in 1:k

            a = A[i,p]
            b = B[p,j]

            term = a*b

            if s === nothing
                s = term
            else
                s += term
            end

        end

        C[i,j] = s

    end

    return C
end

function _field_pdim(f)
    if f isa ScalarField
        return 1
    elseif f isa VectorField
        if f.type == :v2D
            return 2
        elseif f.type == :v3D
            return 3
        else
            error("Unknown VectorField type $(f.type)")
        end
    elseif f isa TensorField
        return 9
    elseif f isa Number
        return 1
    elseif f isa AbstractVector
        return length(f)
    elseif f isa AbstractMatrix
        return length(f)
    else
        error("Unsupported coefficient type $(typeof(f)) in assemble_linear.")
    end
end

mutable struct LinearWorkspace
    invJac::Matrix{Float64}
    ∂h::Matrix{Float64}
    B::Matrix{Float64}
    fe::Vector{Float64}
end

"""
    LinearWorkspace(dim, num_nodes, num_integration_points, outdim, ndofs_loc)

Allocate reusable element-level work arrays for linear-form assembly.
"""
function LinearWorkspace(
    dim::Int,
    num_nodes::Integer,
    num_integration_points::Integer,
    outdim::Integer,
    ndofs_loc::Integer
)
    return LinearWorkspace(
        zeros(3, 3 * num_integration_points),
        zeros(dim, num_nodes * num_integration_points),
        zeros(outdim, ndofs_loc),
        zeros(ndofs_loc)
    )
end

"""
    assemble_element_vector!(ws, elem, gidx, ...)

Compute the local element vector and store it in `ws.fe`.
"""
function assemble_element_vector!(
    ws::LinearWorkspace,
    elem::Integer,
    gidx::Integer,
    jac_all,
    det_all,
    nip::Integer,
    numIntPoints::Integer,
    dim::Integer,
    numNodes::Integer,
    ∇h,
    h,
    intWeights,
    Gprep,
    Wprep,
    op::AbstractOp,
    P::Problem
)
    invJac = ws.invJac
    ∂h = ws.∂h
    B = ws.B
    fe = ws.fe

    jac_slice =
        @view jac_all[gidx * 9 * nip + 1 : (gidx + 1) * 9 * nip]

    jacDet =
        @view det_all[gidx * nip + 1 : (gidx + 1) * nip]

    Jac = reshape(jac_slice, 3, :)

    @inbounds for k in 1:numIntPoints
        #@views invJac[1:3, 3k-2:3k] .=
        #    inv(Jac[1:3, 3k-2:3k])'
        Jk = SMatrix{3,3,Float64}(Jac[1:3, 3k-2:3k])
        invJac[1:3, 3k-2:3k] .= transpose(inv3x3_fast(Jk))
    end

    fill!(∂h, 0.0)

    @inbounds for k in 1:numIntPoints, a in 1:numNodes
        @views begin
            invJk = invJac[1:dim, 3k-2:3k-(3-dim)]
            gha = ∇h[a*3-2:a*3-(3-dim), k]
            dha = ∂h[1:dim, (k-1)*numNodes+a]
            mul!(dha, invJk, gha)
        end
    end

    fill!(fe, 0.0)

    @inbounds for k in 1:numIntPoints
        wcoef = 1.0

        if Wprep !== nothing
            wcoef = _eval_coefficient_at_gp(
                Wprep, elem, view(h, :, k)
            )
            wcoef isa Number ||
                error("assemble_linear: weight must evaluate to a scalar.")
        end

        w = jacDet[k] * intWeights[k] * wcoef

        build_B!(B, op, P, k, h, ∂h, numNodes)

        g_gp = _eval_coefficient_at_gp(
            Gprep, elem, view(h, :, k)
        )

        if g_gp isa Number
            @assert size(B, 1) == 1
            for i in eachindex(fe)
                fe[i] += B[1, i] * g_gp * w
            end
        else
            mul!(fe, transpose(B), g_gp, w, 1.0)
        end
    end

    return nothing
end

"""
    assemble_linear(P::Problem, op, rhs; weight=nothing, domain, threads=:auto)

Assemble a linear finite element operator of the form

    ∫ v ⋅ op ⋅ rhs

where `v` is the test field associated with problem `P`.

Arguments
---------
- `P::Problem` : finite element problem definition
- `op` : operator or matrix chain acting on the right-hand side
- `rhs` : scalar, vector, tensor field, or numeric vector
- `weight` : optional quadrature weight or coefficient
- `domain` : domain specification (`Ω` or `Γ`)
- `gauss` : Number of integration points.
        `:full`       2order+1
        `:reduced`    2order-1
        Int           2order+1 + Int
- `threads` : number of assembly workers; `:auto` uses all Julia threads and
  `1` selects serial assembly

Returns
-------
Global right-hand side vector.
"""
function assemble_linear(
    P::Problem,
    op::AbstractOp,
    g;
    weight = nothing,
    domain = nothing,
    gauss = :full,
    threads = :auto)

    num_threads = resolve_num_threads(threads)
    use_threads = num_threads > 1

    old_blas_threads = LinearAlgebra.BLAS.get_num_threads()
    if use_threads && old_blas_threads != 1
        LinearAlgebra.BLAS.set_num_threads(1)
    end

    try

    gmsh.model.setCurrent(P.name)

    outdim = op_outdim(op, P)

    # allow matrix coefficient for tensor operators
    if g isa AbstractMatrix
        if length(g) != outdim
            error("matrix coefficient size $(size(g)) incompatible with operator output dimension $outdim")
        end
        g = vec(g)
    end
    if g isa AbstractVector
        if length(g) != outdim
            error("assemble_linear: vector coefficient length ($(length(g))) does not match operator output dimension ($outdim).")
        end
    else
        pdim_g = _field_pdim(g)
        if pdim_g != outdim
            error("assemble_linear: coefficient field dimension ($pdim_g) does not match operator output dimension ($outdim).")
        end
    end

    nd = ndofs(P)
    I = Int[]
    V = Float64[]

    phNames = domain === nothing ?
        [mat.phName for mat in P.material] :
        [domain.name]
    
    for phName in phNames
        Gprep = domain === nothing ?
            _prepare_coefficient(g, DomainSpec(:Ω, phName)) :
            _prepare_coefficient(g, domain)
        Wprep = weight === nothing ? nothing :
            domain === nothing ?
            _prepare_coefficient(weight, DomainSpec(:Ω, phName)) :
            _prepare_coefficient(weight, domain)

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for (edim, etag) in dimTags

            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for itype in eachindex(elemTypes)

                et = elemTypes[itype]

                _, _, order, numNodes, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                gorder = gauss == :reduced ? max(1, 2order - 1) :
                    gauss == :full    ? (2order + 1) :
                    gauss isa Int     ? max(1,2order+1 + gauss) :
                    (2order + 1)

                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(gorder))
                #intPoints, intWeights =
                #    gmsh.model.mesh.getIntegrationPoints(
                #        et, "Gauss" * string(2order + 1)
                #    )

                numIntPoints = length(intWeights)

                _, fun, _ =
                    gmsh.model.mesh.getBasisFunctions(
                        et, intPoints, "Lagrange"
                    )

                h = reshape(fun, :, numIntPoints)

                _, dfun, _ =
                    gmsh.model.mesh.getBasisFunctions(
                        et, intPoints, "GradLagrange"
                    )

                ∇h = reshape(dfun, :, numIntPoints)

                # Get Jacobians for all elements of this type
                elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)

                jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)

                numElementsGlobal = length(elemTags_global)
                nip = length(det_all) ÷ numElementsGlobal

                @assert nip == numIntPoints

                # Map element tag to its position in jac_all and det_all
                idxmap = Dict{eltype(elemTags_global),Int}()
                sizehint!(idxmap, numElementsGlobal)

                @inbounds for (gi, gtag) in enumerate(elemTags_global)
                    idxmap[gtag] = gi - 1
                end

                nel = length(elemTags[itype])
                ndofs_loc = P.pdim * numNodes

                block_start = length(V) + 1
                block_entries = nel * ndofs_loc

                resize!(I, length(I) + block_entries)
                resize!(V, length(V) + block_entries)

                workspaces = [
                    LinearWorkspace(
                        P.dim,
                        numNodes,
                        numIntPoints,
                        outdim,
                        ndofs_loc
                    )
                    for _ in 1:num_threads
                ]

                let Ithread = I, Vthread = V
                    if use_threads
                        _run_workers(num_threads) do worker
                            ws = workspaces[worker]
                            for e in _worker_range(nel, num_threads, worker)
                            elem = elemTags[itype][e]
                            gidx = idxmap[elem]

                            assemble_element_vector!(
                                ws,
                                elem,
                                gidx,
                                jac_all,
                                det_all,
                                nip,
                                numIntPoints,
                                P.dim,
                                numNodes,
                                ∇h,
                                h,
                                intWeights,
                                Gprep,
                                Wprep,
                                op,
                                P
                            )

                            element_start =
                                block_start + (e - 1) * ndofs_loc

                            @inbounds for a_loc in 1:ndofs_loc
                                node = div(a_loc - 1, P.pdim) + 1
                                comp = mod(a_loc - 1, P.pdim) + 1

                                gnode =
                                    elemNodeTags[itype][
                                        (e - 1) * numNodes + node
                                    ]

                                p = element_start + a_loc - 1
                                Ithread[p] =
                                    (gnode - 1) * P.pdim + comp
                                Vthread[p] = ws.fe[a_loc]
                            end
                            end
                        end
                    else
                        @inbounds for e in 1:nel
                            ws = workspaces[1]
                            elem = elemTags[itype][e]
                            gidx = idxmap[elem]

                            assemble_element_vector!(
                                ws,
                                elem,
                                gidx,
                                jac_all,
                                det_all,
                                nip,
                                numIntPoints,
                                P.dim,
                                numNodes,
                                ∇h,
                                h,
                                intWeights,
                                Gprep,
                                Wprep,
                                op,
                                P
                            )

                            element_start =
                                block_start + (e - 1) * ndofs_loc

                            for a_loc in 1:ndofs_loc
                                node = div(a_loc - 1, P.pdim) + 1
                                comp = mod(a_loc - 1, P.pdim) + 1

                                gnode =
                                    elemNodeTags[itype][
                                        (e - 1) * numNodes + node
                                    ]

                                p = element_start + a_loc - 1
                                Ithread[p] =
                                    (gnode - 1) * P.pdim + comp
                                Vthread[p] = ws.fe[a_loc]
                            end
                        end
                    end
                end
            end
        end
    end

    rhs = collect(sparsevec(I, V, nd))

    if P.pdim == 1
        return ScalarField([], reshape(rhs, :, 1), [0], [], 1, :scalar, P)

    elseif P.pdim == 2 || P.pdim == 3
        type = P.pdim == 2 ? :v2D : :v3D
        return VectorField([], reshape(rhs, :, 1), [0], [], 1, type, P)

    elseif P.pdim == 9
        return TensorField([], reshape(rhs, :, 1), [0], [], 1, :tensor, P)

    else
        error("assemble_linear: unsupported pdim $(P.pdim).")
    end
    finally
        if use_threads && old_blas_threads != 1
            LinearAlgebra.BLAS.set_num_threads(old_blas_threads)
        end
    end
end

"""
    compliance9_iso(E, nu; penalty=1e8)

Return a 9×9 compliance-like matrix acting on vec(σ) with ordering
(11,12,13,21,22,23,31,32,33). Symmetric part follows isotropic linear elasticity,
antisymmetric part is penalized with `penalty`.
"""
function compliance9_iso(E, nu; penalty=1e8)
    G = E / (2 * (1 + nu))

    # Voigt compliance (engineering shear):
    # [ε11, ε22, ε33, γ12, γ13, γ23] = S6 * [σ11, σ22, σ33, σ12, σ13, σ23]
    S6 = zeros(6, 6)
    S6[1, 1] = 1 / E
    S6[2, 2] = 1 / E
    S6[3, 3] = 1 / E
    S6[1, 2] = -nu / E
    S6[1, 3] = -nu / E
    S6[2, 1] = -nu / E
    S6[2, 3] = -nu / E
    S6[3, 1] = -nu / E
    S6[3, 2] = -nu / E
    S6[4, 4] = 1 / G
    S6[5, 5] = 1 / G
    S6[6, 6] = 1 / G

    # Map 9 -> 6 (take symmetric + engineering shear)
    # v9 = [σ11 σ12 σ13 σ21 σ22 σ23 σ31 σ32 σ33]
    # v6 = [σ11, σ22, σ33, σ12, σ13, σ23] with σ12 := (σ12+σ21)/2 etc.
    P = zeros(6, 9)
    # normals
    P[1, 1] = 1.0       # σ11
    P[2, 5] = 1.0       # σ22
    P[3, 9] = 1.0       # σ33
    # shear sym averages
    P[4, 2] = 0.5
    P[4, 4] = 0.5   # σ12
    P[5, 3] = 0.5
    P[5, 7] = 0.5   # σ13
    P[6, 6] = 0.5
    P[6, 8] = 0.5   # σ23

    # Sym part compliance in 9-space: Ssym = P' * S6 * P
    Ssym = P' * S6 * P

    # Antisym penalty: penalize (σ12-σ21), (σ13-σ31), (σ23-σ32)
    # Add penalty * Q'Q where Q*v9 = [σ12-σ21, σ13-σ31, σ23-σ32]
    Q = zeros(3, 9)
    Q[1, 2] = 1.0
    Q[1, 4] = -1.0
    Q[2, 3] = 1.0
    Q[2, 7] = -1.0
    Q[3, 6] = 1.0
    Q[3, 8] = -1.0
    Santi = penalty * (Q' * Q)

    return Ssym + Santi
end

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
    solveField(K::SystemMatrix, 
               F::SystemVector; 
               support::Vector{BoundaryCondition}=BoundaryCondition[])

single field version:

    solveField(K::SystemMatrix, 
               f::Union{ScalarField,VectorField,TensorField}; 
               support::Vector{BoundaryCondition}=BoundaryCondition[], 
               iterative=false, 
               reltol::Real = sqrt(eps()), 
               maxiter::Int = K.model.non * K.model.dim, 
               preconditioner = Identity(), 
               ordering=true)
    
Solve the linear system

    K * x = F

for a multifield finite element problem.

Dirichlet boundary conditions are imposed through
`BoundaryCondition` objects.

# Arguments

`K`

Global system matrix (`SystemMatrix`).

`F`

Global RHS vector (`SystemVector` or `ScalarField` or `VectorField`)

# Keyword arguments

`support`

Vector of boundary conditions.

# Returns

Tuple of fields corresponding to the Problems stored in `K`.

Example

    (u, p) = solveField(K, F)
"""
function solveField(
    K::SystemMatrix,
    F::SystemVector;
    support::Vector{BoundaryCondition}=BoundaryCondition[])

    # ----------------------------------------------------------
    # 1) Consistency checks
    # ----------------------------------------------------------
    K.problems === nothing &&
        error("solveField: SystemMatrix is not a block system.")

    F.problems === nothing &&
        error("solveField: SystemVector is not a block vector.")

    K.problems == F.problems ||
        error("solveField: Problem ordering mismatch between K and F.")

    problems = K.problems
    offsets = K.offsets

    A = K.A
    b = F.a

    ndof, nsteps = size(b)

    # ----------------------------------------------------------
    # 2) Collect global constrained DOFs
    # ----------------------------------------------------------
    fixed = Int[]

    for bc in support

        P = bc.problem
        idx = findfirst(q -> q === P, problems)

        idx === nothing &&
            error("solveField: BC refers to Problem not in system.")

        offset = offsets[idx]

        local_dofs = constrainedDoFs(P, [bc])

        append!(fixed, offset .+ local_dofs)
    end

    fixed = unique(fixed)
    sort!(fixed)

    # ----------------------------------------------------------
    # 3) Define free DOFs
    # ----------------------------------------------------------
    all_dofs = collect(1:ndof)
    free = setdiff(all_dofs, fixed)

    # ----------------------------------------------------------
    # 4) Reduced system with non-homogeneous BC
    # ----------------------------------------------------------

    x = zeros(ndof, nsteps)
    xD = zeros(ndof, 1)

    # fill xD
    for bc in support
        P = bc.problem
        idx = findfirst(q -> q === P, problems)
        offset = offsets[idx]

        x_local = applyBoundaryConditions(P, [bc])
        local_dofs = constrainedDoFs(P, [bc])

        xD[offset.+local_dofs, :] .= x_local.a[local_dofs, :]
    end

    fixed = unique(vcat([offsets[findfirst(q -> q === bc.problem, problems)] .+
                         constrainedDoFs(bc.problem, [bc])
                         for bc in support]...))

    free = setdiff(1:ndof, fixed)

    A_ff = A[free, free]
    for i in 1:nsteps
        b_f = b[free] - A[free, fixed] * xD[fixed]

        x[free, i] = A_ff \ b_f
    end
    x[fixed] = xD[fixed]

    # fixed DOFs remain zero (homogeneous)

    # ----------------------------------------------------------
    # 5) Reconstruct fields
    # ----------------------------------------------------------
    results = Vector{Any}(undef, length(problems))

    for (i, P) in enumerate(problems)

        offset = offsets[i]
        nloc = ndofs(P)

        xloc = x[offset+1:offset+nloc]

        if P.pdim == 1
            results[i] = ScalarField([], reshape(xloc, :, nsteps), [0], [], 1, :scalar, P)

        elseif P.pdim == 2 || P.pdim == 3
            type = P.pdim == 2 ? :v2D : :v3D
            results[i] = VectorField([], reshape(xloc, :, nsteps), [0], [], 1, type, P)

        elseif P.pdim == 9
            results[i] = TensorField([], reshape(xloc, :, nsteps), [0], [], 1, :tensor, P)

        else
            error("solveField: unsupported pdim $(P.pdim).")
        end
    end

    if length(results) == 1
        return results[1]
    else
        return tuple(results...)
    end
end

function solveField(
    Ks::SymmetricSystemMatrix,
    f::Union{ScalarField,VectorField,TensorField};
    support::Vector{BoundaryCondition}=BoundaryCondition[],
    iterative=false,
    reltol::Real=sqrt(eps()),
    maxiter::Int=Ks.parent.model.non * Ks.parent.model.dim,
    preconditioner=Identity(),
    ordering=true
)
    K = Ks.parent
    A = Symmetric(K.A, Ks.uplo)

    problem = K.model
    fixed = constrainedDoFs(problem, support)
    free = freeDoFs(problem, support)

    u = copy(f)
    fill!(u.a, 0.0)
    applyBoundaryConditions!(u, support)

    # Fontos: itt is a Symmetric nézetet kell használni.
    f_kin = A[:, fixed] * u.a[fixed, 1]

    A_ff = Symmetric(K.A[free, free], Ks.uplo)
    b_f = f.a[free, 1] - f_kin[free]

    if iterative
        u.a[free] = cg(
            A_ff,
            b_f;
            Pl=preconditioner,
            reltol=reltol,
            maxiter=maxiter
        )
    else
        u.a[free] = A_ff \ b_f
    end

    return u
end

###############################################################
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

Create a weak-form DSL gradient operator applied to `P`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object for `∇`.

# Example
```julia
K = ∫(Grad(Pu) ⋅ Grad(Pu); Ω="solid")
```
"""
Grad(P) = OpApplied(P, GradOp())

"""
    Div(P)

Create a weak-form DSL divergence operator applied to `P`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object for `∇⋅`.

# Example
```julia
A = ∫(Div(Pu) ⋅ Div(Pu); Ω="domain")
```
"""
Div(P) = OpApplied(P, DivOp())

"""
    Curl(P)

Create a weak-form DSL curl operator applied to `P`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object for curl.

# Example
```julia
A = ∫(Curl(Pu) ⋅ Curl(Pu); Ω="domain")
```
"""
Curl(P) = OpApplied(P, CurlOp())

"""
    SymGrad(P)

Create a weak-form DSL symmetric-gradient operator applied to `P`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object for `ε(u)`.

# Example
```julia
K = ∫(SymGrad(Pu) ⋅ C ⋅ SymGrad(Pu); Ω="solid")
```
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

Create a weak-form DSL tensor-divergence operator applied to `P`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object for tensor divergence.

# Example
```julia
A = ∫(TensorDiv(Pσ) ⋅ TensorDiv(Pσ); Ω="solid")
```
"""
TensorDiv(P) = OpApplied(P, TensorDivOp())

"""
    Adv(P)

Create a weak-form DSL advection operator applied to `P`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object for advection terms.

# Example
```julia
A = ∫(Adv(Pu) ⋅ Id(Pu); Ω="domain")
```
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

Create a weak-form DSL axial gradient operator applied to `P`.

# Arguments
- `P`: Field descriptor (`Problem`) used in the weak form.

# Returns
- `OpApplied`: Operator application object representing `t ⋅ ∇u`.

Directional derivative operator along a prescribed axial direction.

Computes derivatives projected onto a predefined axis.

Unlike `TangentialGrad`, the axial direction is prescribed explicitly and
does not depend on the local element geometry.

Typical applications include:

- anisotropic constitutive laws,
- fiber-reinforced materials,
- projected strain operators,
- directional transport problems.

See also: `TangentialGrad`

# Example
```julia
K = ∫(AxialGrad(Pu) ⋅ (E*A) ⋅ AxialGrad(Pu); Ω="truss")
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
    threads=:auto, assembly::Symbol=:ijv)

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
    threads=:auto, assembly::Symbol=:ijv)

    _assemble(expr.a, dom, weight, gauss, threads, assembly) +
        _assemble(expr.b, dom, weight, gauss, threads, assembly)

end

function _assemble(
    term::WeakTerm{LinearTerm}, dom, weight=nothing, gauss=:full,
    threads=:auto, assembly::Symbol=:ijv)

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

`:full` or `:reduced` integration, or `Int` that is a signed number:  the increment of Gauss points relative to `:full`.

# Keyword arguments

`Ω`

Volume physical group name.

`Γ`

Boundary physical group name.

# Returns

`SystemMatrix` or `ScalarField`, `VectorField`, `TensorField`
"""
function ∫(
    expr::WeakExpr; Ω=nothing, Γ=nothing, weight=nothing, gauss=:full,
    threads=:auto, assembly::Symbol=:ijv)

    _check_scalarfields(expr)

    assembly in (:csc, :ijv) || error(
        "Unsupported matrix assembly mode $assembly. Expected :csc or :ijv."
    )

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    return _assemble(expr, dom, weight, gauss, threads, assembly)

end

function ∫(t::BilinearTerm;
    Ω=nothing, Γ=nothing, weight=nothing,
    gauss=:full,
    threads=:auto,
    assembly::Symbol=:ijv,
    csc_matrix=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto)

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

    K = assemble_operator(
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

    assembly === :csc && csc_matrix === nothing && dropzeros!(K.A)
    return K

end

function ∫(a::OpApplied, b::OpApplied;
    Ω=nothing, Γ=nothing, weight=nothing,
    gauss=:full,
    threads=:auto,
    assembly::Symbol=:ijv,
    csc_matrix=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto)

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

    K = assemble_operator(
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

    assembly === :csc && csc_matrix === nothing && dropzeros!(K.A)
    return K

end

function ∫(t::LinearTerm; Ω=nothing, Γ=nothing, weight=nothing, gauss=:full, threads=:auto)

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
    ∫Ω(name, expr)

Convenience wrapper for volume integration on physical group `name`.

# Arguments
- `name`: Gmsh physical group name used as domain `Ω`.
- `expr::WeakExpr`: Weak-form expression to assemble.

# Returns
- `SystemMatrix`: Assembled matrix over the selected volume.

# Example
```julia
K = ∫Ω("solid", Grad(Pu) ⋅ Grad(Pu))
```
"""
∫Ω(name, expr; kwargs...) = ∫(expr; Ω=name, kwargs...)

"""
    ∫Γ(name, expr)

Convenience wrapper for boundary integration on physical group `name`.

# Arguments
- `name`: Gmsh physical group name used as boundary `Γ`.
- `expr::WeakExpr`: Weak-form expression to assemble.

# Returns
- `SystemMatrix`: Assembled matrix over the selected boundary.

# Example
```julia
KΓ = ∫Γ("loaded_boundary", Id(Pu) ⋅ Id(Pu))
```
"""
∫Γ(name, expr; kwargs...) = ∫(expr; Γ=name, kwargs...)

const ε = SymGrad

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
    assembly::Symbol=:ijv,
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

# Integral wrapper

# ---------------------------------------------------------------------------

function ∫(expr::CompoundBilinear;
    Ω=nothing, Γ=nothing, weight=nothing,
    gauss=:auto,
    threads=:auto,
    assembly::Symbol=:ijv,
    csc_matrix=nothing,
    element_chunk_size::Union{Integer,Symbol}=:auto)
    #@show typeof(expr.left)
    #@show typeof(expr.right)

    left_terms = _side_terms(expr.left)
    right_terms = _side_terms(expr.right)

    assembly ∈ (:ijv, :csc) || error(
        "Unsupported matrix assembly mode $assembly. " *
        "Expected :csc or :ijv."
    )

    assembly === :ijv && csc_matrix !== nothing && error(
        "csc_matrix is only valid with assembly=:csc."
    )

    if assembly === :csc
        dom = _domain_spec(; Ω=Ω, Γ=Γ)
        num_threads = resolve_num_threads(threads)

        isempty(left_terms) && error("Cannot assemble an empty compound bilinear form.")
        isempty(right_terms) && error("Cannot assemble an empty compound bilinear form.")

        first_bt = _make_bilinear_term(left_terms[1], expr.mats, right_terms[1])
        Pu = first_bt.b.P
        Ps = first_bt.a.P

        if csc_matrix === nothing
            Kcsc = build_csc_pattern(Pu, Ps; domain=dom)
            GC.gc(false)
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
        end

        coordinates = dense_node_coordinates(Pu)
        GC.gc(false)

        nzval_buffers = Vector{Vector{Float64}}(undef, num_threads)
        nzval_buffers[1] = Kcsc.nzval
        for worker in 2:num_threads
            nzval_buffers[worker] = zeros(Float64, length(Kcsc.nzval))
        end

        for LL in left_terms
            for R in right_terms
                bt = _make_bilinear_term(LL, expr.mats, R)
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

        Kcsc === nothing && error("Cannot assemble an empty compound bilinear form.")
        # Preserve the complete reusable pattern supplied by the caller.
        # A matrix created internally is finalized like the original IJV path.
        csc_matrix === nothing && dropzeros!(Kcsc)

        return SystemMatrix(Kcsc, Pu, Ps)
    end
    #@show left_terms
    #@show right_terms

    #K = nothing

    #for LL in left_terms
    #    for R in right_terms
    #        bt = _make_bilinear_term(LL, expr.mats, R)
    #        #@show bt
    #        g = _select_gauss(LL, R, gauss)

    #        Ki = ∫(bt; Ω=Ω, Γ=Γ, weight=weight, gauss=g)
    #        K = K === nothing ? Ki : K + Ki
    #    end
    #end

    first = true

    I = nothing
    J = nothing
    V = nothing
    nr = 0
    nc = 0

    Pu = nothing
    Ps = nothing

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    for LL in left_terms
        for R in right_terms

            bt = _make_bilinear_term(LL, expr.mats, R)
            g = _select_gauss(LL, R, gauss)

            if first

                Pu = bt.b.P
                Ps = bt.a.P

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
                    threads=threads
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
                    threads=threads
                )

            end
        end
    end

    K = sparse(I,J,V,nr,nc)
    dropzeros!(K)

    return SystemMatrix(K, Pu, Ps)

    #return K
end

#∫(expr::CompoundBilinearWeight; Ω=nothing, Γ=nothing, gauss=:auto) = ∫(expr.CB; Ω=Ω, Γ=Γ, weight=expr.weight, gauss=gauss)
#∫(expr::LinearTermWeight; Ω=nothing, Γ=nothing, gauss=:auto) = ∫(expr.CB; Ω=Ω, Γ=Γ, weight=expr.weight, gauss=gauss)

########################################################################

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

    free = setdiff(collect(1:ndof), fixed)

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
    check_multifield_system_compatibility(K::SystemMatrix,
                                          C::SystemMatrix)

Validate that two block system matrices are compatible for multifield
time integration.
"""
function check_multifield_system_compatibility(
    K::SystemMatrix,
    C::SystemMatrix
)
    K.problems === nothing &&
        error("check_multifield_system_compatibility: K is not a block SystemMatrix.")

    C.problems === nothing &&
        error("check_multifield_system_compatibility: C is not a block SystemMatrix.")

    K.offsets === nothing &&
        error("check_multifield_system_compatibility: K.offsets === nothing.")

    C.offsets === nothing &&
        error("check_multifield_system_compatibility: C.offsets === nothing.")

    K.problems == C.problems ||
        error("check_multifield_system_compatibility: K.problems and C.problems differ.")

    K.offsets == C.offsets ||
        error("check_multifield_system_compatibility: K.offsets and C.offsets differ.")

    size(K.A) == size(C.A) ||
        error("check_multifield_system_compatibility: matrix sizes differ: size(K.A)=$(size(K.A)) vs size(C.A)=$(size(C.A)).")

    size(K.A, 1) == size(K.A, 2) ||
        error("check_multifield_system_compatibility: K must be square.")

    size(C.A, 1) == size(C.A, 2) ||
        error("check_multifield_system_compatibility: C must be square.")

    return nothing
end

"""
    split_multifield_solution(X::AbstractMatrix,
                              problems::Vector{Problem},
                              offsets::Vector{Int},
                              t::AbstractVector)

Split a global multifield solution matrix into field objects and return
them as a tuple in block order.
"""
function split_multifield_solution(
    X::AbstractMatrix,
    problems::Vector{Problem},
    offsets::Vector{Int},
    t::AbstractVector
)
    nsteps = size(X, 2)
    results = Vector{Any}(undef, length(problems))

    for (i, P) in enumerate(problems)
        off = offsets[i]
        nloc = P.non * P.pdim #ndofs(P)
        Xloc = X[off+1:off+nloc, :]

        if P.pdim == 1
            results[i] = ScalarField([], Matrix(Xloc), collect(t), [], nsteps, :scalar, P)
        elseif P.pdim == 2
            results[i] = VectorField([], Matrix(Xloc), collect(t), [], nsteps, :v2D, P)
        elseif P.pdim == 3
            results[i] = VectorField([], Matrix(Xloc), collect(t), [], nsteps, :v3D, P)
        elseif P.pdim == 9
            results[i] = TensorField([], Matrix(Xloc), collect(t), [], nsteps, :tensor, P)
        else
            error("split_multifield_solution: unsupported pdim $(P.pdim).")
        end
    end

    if length(results) == 1
        return results[1]
    else
        return tuple(results...)
    end
end

"""
    FDM(K::SystemMatrix, 
        C::SystemMatrix, 
        q::Union{ScalarField,VectorField,TensorField,SystemVector}, 
        bc::Vector{BoundaryCondition}, 
        X0::Union{ScalarField,VectorField,TensorField,SystemVector}, 
        n::Int, 
        Δt::Float64; 
        ϑ=0.5)

alias:

    FDM(K, C, q, X0, n, Δt; ϑ=0.5, support=bc)

Solves a transient diffusion-type problem (e.g. heat conduction) using the
finite difference method in time (ϑ-method).

The semi-discrete system
```

C * Ẋ(t) + K * X(t) = q(t)

```
is integrated in time using the ϑ-scheme:
- ϑ = 0   : Forward Euler (explicit)
- ϑ = 1/2 : Crank–Nicolson
- ϑ = 1   : Backward Euler (implicit)
- 0 < ϑ < 1 : intermediate schemes

The method supports:
- time-independent or time-dependent load vectors `q`
- time-independent or time-dependent Dirichlet boundary conditions
- scalar, vector and tensor unknowns (`ScalarField`, `VectorField`, `TensorField`)
- multi-field block systems (`SystemVector`)
- consistent treatment of prescribed DOFs via matrix partitioning

Boundary conditions are applied *solver-side*:
- constrained DOFs are prescribed directly at each time step
- free DOFs are solved from the reduced system
- contributions of constrained DOFs enter the RHS through
```

* K_fc X_c - C_fc Ẋ_c

```

If `q.nsteps == 1`, the load is treated as time-independent.  
If `q.nsteps == n`, a true ϑ-weighted load
```

q^{n+ϑ} = (1-ϑ) q^n + ϑ q^{n+1}

```
is used.

For `ϑ = 0` and diagonal `C`, a fully explicit update is used.

---

### Arguments
- `K::SystemMatrix`  
  Diffusion (conductivity / stiffness) matrix.
- `C::SystemMatrix`  
  Capacity (mass / heat capacity) matrix.
- `q`  
  Load / source vector:
  - single-field: `ScalarField`, `VectorField`, `TensorField`
  - multi-field: `SystemVector`
- `bc::Vector{BoundaryCondition}`  
  Dirichlet boundary conditions (possibly time-dependent).
- `X0`  
  Initial condition (same type as `q`).  
  On constrained DOFs this is overridden by `bc`.
- `n::Int`  
  Number of time steps.
- `Δt::Float64`  
  Time step size.
- `ϑ::Float64` (keyword, default = 0.5)  
  Parameter of the ϑ-method.

---

### Returns

Single-field:
- `X::Union{ScalarField,VectorField,TensorField}`  

Multifield:
- tuple of fields corresponding to each problem

Each result contains:
- nodal values (`ndof × nsteps`)
- time vector `t = 0:Δt:(n-1)Δt`

---

### Notes
- Stability depends on `ϑ` and the spectrum of the generalized eigenproblem
  `K x = λ C x`. For `ϑ ≥ 1/2` the method is unconditionally stable.
- The explicit case (`ϑ = 0`) requires a diagonal `C` matrix.
- For multifield systems, all fields must have the same number of time steps.
- The algorithm is algebraic and independent of the physical interpretation
  of the field.
"""
function FDM(
    K::SystemMatrix,
    C::SystemMatrix,
    q::SystemVector,
    bc::Vector{BoundaryCondition},
    X0::SystemVector,
    n::Int,
    Δt::Float64;
    ϑ=0.5)

    all(isfinite, K.A.nzval) ||
        error("FDM: K contains non-finite entries.")

    all(isfinite, C.A.nzval) ||
        error("FDM: C contains non-finite entries.")

    all(isfinite, q.a) ||
        error("FDM: q contains non-finite entries.")

    all(isfinite, X0.a) ||
        error("FDM: X0 contains non-finite entries.")

    @assert size(K.A,1) == size(K.A,2)
    @assert size(C.A,1) == size(C.A,2)
    @assert size(K.A) == size(C.A)
    @assert size(q.a,1) == size(K.A,1)
    @assert size(X0.a,1) == size(K.A,1)
    @assert all(isfinite, K.A.nzval)
    @assert all(isfinite, C.A.nzval)
    @assert all(isfinite, q.a)
    @assert all(isfinite, X0.a)
    # ------------------------------------------------------------------
    # 1) Compatibility checks
    # ------------------------------------------------------------------
    check_multifield_system_compatibility(K, C)

    q.problems === nothing &&
        error("FDM: q must be a block SystemVector.")

    X0.problems === nothing &&
        error("FDM: X0 must be a block SystemVector.")

    K.problems == q.problems ||
        error("FDM: Problem ordering mismatch between K and q.")

    K.problems == X0.problems ||
        error("FDM: Problem ordering mismatch between K and X0.")

    K.offsets == q.offsets ||
        error("FDM: Offset mismatch between K and q.")

    K.offsets == X0.offsets ||
        error("FDM: Offset mismatch between K and X0.")

    ndof = size(K.A, 1)
    size(q.a,1) == ndof ||
        error("FDM: size(q.a,1) = $(size(q.a,1)) does not match ndof = $ndof.")

    size(X0.a,1) == ndof ||
        error("FDM: size(X0.a,1) = $(size(X0.a,1)) does not match ndof = $ndof.")

    n >= 1 || error("FDM: n must be at least 1.")
    Δt > 0 || error("FDM: Δt must be positive.")

    nq = size(q.a, 2)
    nq == 1 || nq == n ||
        error("FDM: q must contain either 1 or n time steps.")

    nX0 = size(X0.a, 2)
    nX0 == 1 || nX0 == n ||
        error("FDM: X0 must contain either 1 or n time steps.")

    0.0 <= ϑ <= 1.0 ||
        error("FDM: ϑ must be between 0 and 1.")

    # ------------------------------------------------------------------
    # 2) BC data
    # ------------------------------------------------------------------
    free, fix, xD = multifield_bc_data(K, bc; nsteps=n)

    # ------------------------------------------------------------------
    # 3) Allocate history
    # ------------------------------------------------------------------
    X = zeros(ndof, n)
    #t = zeros(n)
    t = collect(0:Δt:(n-1)*Δt)

    # initial condition
    X .= X0.a
    if !isempty(fix)
        X[fix, 1] .= xD[fix, 1]
    end

    # ------------------------------------------------------------------
    # 4) Reduced matrices
    # ------------------------------------------------------------------
    K0 = K.A
    C0 = C.A

    Kff = K0[free, free]
    Cff = C0[free, free]

    Kfc = K0[free, fix]
    Cfc = C0[free, fix]
    #if !isempty(fix)
    #    Kfc = K0[free, fix]
    #    Cfc = C0[free, fix]
    #else
    #    Kfc = zeros(length(free), 0)
    #    Cfc = zeros(length(free), 0)
    #end

    # ------------------------------------------------------------------
    # 5) Time stepping
    # ------------------------------------------------------------------
    is_diag_Cff = isdiag(Cff)

    if ϑ == 0 && is_diag_Cff
        d = diag(Cff)

        all(>(0), d) ||
            error("FDM: non-positive diagonal entry detected in Cff on free DOFs.")

        invd = 1.0 ./ d
        has_fix = !isempty(fix)

        for i in 2:n
            qi = size(q.a, 2) == 1 ? 1 : i - 1

            @views begin
                qn = q.a[free, qi]
                xfree_n = X[free, i-1]
            end

            rhs =
                qn -
                Kff * xfree_n

            if has_fix
                @views begin
                    xc_n = xD[fix, i-1]
                    xc_np1 = xD[fix, i]
                end

                rhs .-=
                    Kfc * xc_n +
                    Cfc * ((xc_np1 - xc_n) ./ Δt)
            end

            @views X[free, i] .=
                xfree_n .+ Δt .* invd .* rhs

            if has_fix
                @views X[fix, i] .= xD[fix, i]
            end
        end

    else
        A = Cff + ϑ * Δt * Kff
        B = Cff - (1 - ϑ) * Δt * Kff
        luA = lu(A)

        has_fix = !isempty(fix)

        if has_fix
            Afc = Cfc + ϑ * Δt * Kfc
            Bfc = Cfc - (1 - ϑ) * Δt * Kfc
        end

        x_n = copy(@view X[free, 1])
        x_np1 = similar(x_n)
        rhs = similar(x_n)

        for i in 2:n
            # rhs = B * x_n
            mul!(rhs, B, x_n)
   
            if size(q.a, 2) == 1
                @views rhs .+= Δt .* q.a[free, 1]
            else
                @views rhs .+= Δt .* ((1 - ϑ) .* q.a[free, i-1] .+ ϑ .* q.a[free, i])
            end
   
            if has_fix
                @views begin
                    xc_n = xD[fix, i-1]
                    xc_np1 = xD[fix, i]
   
                    mul!(rhs, Afc, xc_np1, -1.0, 1.0)
                    mul!(rhs, Bfc, xc_n, 1.0, 1.0)
                end
            end
   
            ldiv!(x_np1, luA, rhs)
   
            @views X[free, i] .= x_np1
   
            if has_fix
                @views X[fix, i] .= xD[fix, i]
            end
   
            x_n, x_np1 = x_np1, x_n
        end
    end

    return split_multifield_solution(X, K.problems, K.offsets, t)
end

FDM(K::SystemMatrix, C::SystemMatrix, q::SystemVector, X0::SystemVector, n::Int, Δt::Float64; ϑ=0.5, support=Vector{BoundaryCondition}()) =
    FDM(K, C, q, support, X0, n, Δt, ϑ=ϑ)

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

"""
    smallestEigenValue(K::SystemMatrix, 
                       C::SystemMatrix; 
                       support=Vector{BoundaryCondition}())

Compute the smallest eigenvalue λₘᵢₙ of the generalized eigenproblem

    K * ϕ = λ * C * ϕ

after applying Dirichlet boundary conditions.

Return: `λₘᵢₙ`

Types:
- `K`: SystemMatrix
- `C`: SystemMatrix
- `support`: Vector{BoundaryCondition}
- `λₘᵢₙ`: Float64
"""
function smallestEigenValue(
    K::SystemMatrix,
    C::SystemMatrix;
    support=Vector{BoundaryCondition}()
)
    # ----------------------------------------------------------
    # 1) Compatibility checks
    # ----------------------------------------------------------
    size(K.A) == size(C.A) ||
        error("smallestEigenValue: K and C must have the same size.")

    size(K.A, 1) == size(K.A, 2) ||
        error("smallestEigenValue: K must be square.")

    size(C.A, 1) == size(C.A, 2) ||
        error("smallestEigenValue: C must be square.")

    # ----------------------------------------------------------
    # 2) Free DOFs (single + multifield)
    # ----------------------------------------------------------
    free = nothing
    if K.problems === nothing
        free = freeDoFs(K.model, support)
    else
        free, _, _ = multifield_bc_data(K, support; nsteps=1)
    end

    # ----------------------------------------------------------
    # 3) Reduced matrices
    # ----------------------------------------------------------
    K0 = K.A[free, free]
    C0 = C.A[free, free]

    # ----------------------------------------------------------
    # 4) Eigen solve (shift-invert near zero)
    # ----------------------------------------------------------
    ϕ = nothing
    λ = nothing
    try
        λ, ϕ = Arpack.eigs(
            K0, C0,
            nev=1,
            which=:LR,
            sigma=1e-8,
            maxiter=10000
        )
    catch
        # fallback (ha shift-invert nem konvergál)
        λ, ϕ = Arpack.eigs(
            K0, C0,
            nev=1,
            which=:SM,
            maxiter=10000
        )
    end

    # ----------------------------------------------------------
    # 5) Error estimate
    # ----------------------------------------------------------
    r = K0 * ϕ[:, 1] - λ[1] * C0 * ϕ[:, 1]
    err = norm(r) / norm(K0 * ϕ[:, 1])

    if err > 1e-3
        @warn("smallestEigenValue: relative residual too large: $err")
    end

    # ----------------------------------------------------------
    # 6) Return
    # ----------------------------------------------------------
    λmin = abs(real(λ[1]))

    return λmin
end

"""
    largestEigenValue(K::SystemMatrix, 
                      C::SystemMatrix; 
                      support=Vector{BoundaryCondition}())

Compute the largest eigenvalue λₘₐₓ of the generalized eigenproblem

    K * ϕ = λ * C * ϕ

after applying Dirichlet boundary conditions.

This value is critical for stability analysis of explicit time integration schemes.

Return: `λₘₐₓ`

Types:
- `K`: SystemMatrix
- `C`: SystemMatrix
- `support`: Vector{BoundaryCondition}
- `λₘₐₓ`: Float64
"""
function largestEigenValue(
    K::SystemMatrix,
    C::SystemMatrix;
    support=Vector{BoundaryCondition}()
)
    # ----------------------------------------------------------
    # 1) Compatibility checks
    # ----------------------------------------------------------
    size(K.A) == size(C.A) ||
        error("largestEigenValue: K and C must have the same size.")

    size(K.A, 1) == size(K.A, 2) ||
        error("largestEigenValue: K must be square.")

    size(C.A, 1) == size(C.A, 2) ||
        error("largestEigenValue: C must be square.")

    # ----------------------------------------------------------
    # 2) Free DOFs (single + multifield)
    # ----------------------------------------------------------
    free = nothing
    if K.problems === nothing
        free = freeDoFs(K.model, support)
    else
        # multifield
        free, _, _ = multifield_bc_data(K, support; nsteps=1)
    end

    # ----------------------------------------------------------
    # 3) Reduced matrices
    # ----------------------------------------------------------
    K0 = K.A[free, free]
    C0 = C.A[free, free]

    # ----------------------------------------------------------
    # 4) Eigen solve
    # ----------------------------------------------------------
    λ, ϕ = Arpack.eigs(K0, C0, nev=1, which=:LM)

    # ----------------------------------------------------------
    # 5) Error estimate
    # ----------------------------------------------------------
    r = K0 * ϕ[:, 1] - λ[1] * C0 * ϕ[:, 1]
    err = norm(r) / norm(K0 * ϕ[:, 1])

    if err > 1e-3
        @warn("largestEigenValue: relative residual too large: $err")
    end

    # ----------------------------------------------------------
    # 6) Return
    # ----------------------------------------------------------
    λmax = abs(real(λ[1]))

    return λmax
end

"""
    solveEigenFields(K::SystemMatrix, M::SystemMatrix; n=6, fmin=0.0)

Solve eigenproblem for multifield system and return field-wise Eigen objects.

Usage:
    u_modes, p_modes = solveEigenFields(K, M)
"""
function solveEigenFields(
    K::SystemMatrix,
    M::SystemMatrix;
    n = 6,
    fmin = 0.0
)

    # ----------------------------------------------------------
    # 1) Check multifield
    # ----------------------------------------------------------
    K.problems === nothing &&
        error("solveEigenFields: use solveEigenModes for single-field problems.")

    K.problems == M.problems ||
        error("solveEigenFields: K and M must have same block structure.")

    # ----------------------------------------------------------
    # 2) Eigen solve (GLOBAL!)
    # ----------------------------------------------------------
    ω2min = (2π * fmin)^2

    λ, ϕ = Arpack.eigs(
        K.A,
        M.A,
        nev = n,
        which = :LR,
        sigma = ω2min,
        maxiter = 10000
    )

    f = sqrt.(abs.(real(λ))) ./ (2π)
    ϕ = real(ϕ)

    # ----------------------------------------------------------
    # 3) Pack global Eigen
    # ----------------------------------------------------------
    eig_global = Eigen(f, ϕ, nothing, K.problems, K.offsets)

    # ----------------------------------------------------------
    # 4) Split → Eigen-ek mezőnként
    # ----------------------------------------------------------
    return splitEigenToEigen(eig_global)
end

function splitEigenToEigen(eig::Eigen)

    problems = eig.problems
    offsets  = eig.offsets
    ϕ        = eig.ϕ
    f        = eig.f

    results = Vector{Eigen}(undef, length(problems))

    for (i, P) in enumerate(problems)

        off = offsets[i]
        nd  = P.non * P.pdim

        ϕloc = ϕ[off+1:off+nd, :]

        results[i] = Eigen(f, ϕloc, P, nothing, nothing)
    end

    return tuple(results...)
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
