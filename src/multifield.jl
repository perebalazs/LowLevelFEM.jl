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

export ∫
export ∫Ω
export ∫Γ

export solveField

export ε

export ⋅

export solveEigenFields


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
struct TangentialGradOp <: AbstractOp end

TangentialGrad(P) = OpApplied(P, TangentialGradOp())

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

struct DomainSpec
    kind::Symbol     # :Ω vagy :Γ
    name::String
end

Base.show(io::IO, d::DomainSpec) =
    print(io, "$(d.kind)=\"$(d.name)\"")

@inline function _build_elemwise_coeff_dict(coef::ScalarField, domain::DomainSpec)
    p = elementsToNodes(coef)
    p = nodesToElements(p, onPhysicalGroup=domain.name)
    return Dict(zip(p.numElem, p.A))  # elemTag => coeff nodal vector(s)
end

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

"""
    assemble_operator(Pu::Problem, 
                      op_u::AbstractOp, 
                      Ps::Problem, 
                      op_s::AbstractOp; 
                      coefficient::Union{Number,ScalarField}=1.0, 
                      weight=nothing, 
                      domain=nothing)

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

Optional weighting tensor applied to the operator value. (Deprecated: kept for backward compatibility.)

`domain`

Optional domain specification (`Ω` or `Γ`).

# Returns

`SystemMatrix`

Sparse matrix representing the discretized bilinear form.
"""
function assemble_operator(
    Pu::Problem, op_u::AbstractOp,
    Ps::Problem, op_s::AbstractOp;
    coefficient::Union{Number,ScalarField,AbstractMatrix,AbstractVector}=1.0,
    weight=nothing,
    domain=nothing)

    @assert Pu.name == Ps.name "Both problems must refer to the same gmsh model/mesh."
    @assert Pu.dim == Ps.dim "Both problems must have the same spatial dimension."
    gmsh.model.setCurrent(Pu.name)

    # output dims must match for the dot-product integrand
    #out_u = op_outdim(op_u, Pu)
    #out_s = op_outdim(op_s, Ps)
    #@assert out_u == out_s "Operator output dims mismatch: $out_u vs $out_s."
    
    # operator output dimensions
    out_u = op_outdim(op_u, Pu)
    out_s = op_outdim(op_s, Ps)
    
    if coefficient isa Number || coefficient isa ScalarField
        @assert out_u == out_s
    
    elseif coefficient isa AbstractMatrix
        @assert size(coefficient,1) == out_s
        @assert size(coefficient,2) == out_u
    
    elseif coefficient isa AbstractVector
        if length(coefficient) == 1 &&
           (coefficient[1] isa Number || coefficient[1] isa ScalarField)
    
            @assert out_u == out_s
    
        else
            A1 = coefficient[1]
            An = coefficient[end]
    
            @assert A1 isa AbstractMatrix
            @assert An isa AbstractMatrix
    
            @assert size(A1,1) == out_s
            @assert size(An,2) == out_u
    
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
        w1, w2 = size(weight)
        if w1 != out_u || w2 != out_u
            error("Weight dimension mismatch: expected $(out_u)×$(out_u) but got $(w1)×$(w2)")
        end
    end

    # ------------------------------------------------------------
    # prepare coefficient
    # ------------------------------------------------------------

    # estimate IJV length: crude but OK for notebook prototype
    lengthOfIJV = LowLevelFEM.estimateLengthOfIJV(Pu) * max(1, Ps.pdim) * max(1, Pu.pdim)
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
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

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        isempty(dimTags) && error("assemble_operator: physical group \"$phName\" not found.")

        for (edim, etag) in dimTags

            # optional dimension checks (matches your DSL design)
            if domkind === :Ω
                edim == Pu.dim || error("Ω=\"$phName\" has dim=$edim but problem.dim=$(Pu.dim)")
            elseif domkind === :Γ
                edim < Pu.dim || error("Γ=\"$phName\" has dim=$edim but expected < $(Pu.dim)")
            end

            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

            for itype in eachindex(elemTypes)
                et = elemTypes[itype]
                _, _, order, numNodes::Int64, _, _ = gmsh.model.mesh.getElementProperties(et)

                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)

                _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)

                _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)

                # buffers
                nel = length(elemTags[itype])
                nnet = zeros(Int, nel, numNodes)
                invJac = zeros(3, 3numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)

                ndofs_u_loc = Pu.pdim * numNodes
                ndofs_s_loc = Ps.pdim * numNodes

                Bu = zeros(out_u, ndofs_u_loc)
                Bs = zeros(out_s, ndofs_s_loc)
                Ke = zeros(ndofs_s_loc, ndofs_u_loc)
                tmp = zeros(size(Bs,2), size(Bu,1))

                # connectivity table
                @inbounds for e in 1:nel
                    for a in 1:numNodes
                        nnet[e, a] = elemNodeTags[itype][(e-1)*numNodes+a]
                    end
                end

                tmpBu = similar(Bu)

                # element loop
                @inbounds for e in 1:nel
                    elem = elemTags[itype][e]

                    jac, jacDet, _ = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)

                    @inbounds for k in 1:numIntPoints
                        invJac[1:3, 3k-2:3k] .= inv(Jac[1:3, 3k-2:3k])'
                    end

                    # physical gradients of basis
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numIntPoints, a in 1:numNodes
                        invJk = invJac[1:dim, 3k-2:3k-(3-dim)]
                        gha = ∇h[a*3-2:a*3-(3-dim), k]
                        ∂h[1:dim, (k-1)*numNodes+a] .= invJk * gha
                    end

                    fill!(Ke, 0.0)

                    # integrate
                    @inbounds for k in 1:numIntPoints
                        if Cprep isa AbstractVector
                            mats = [_eval_coefficient_at_gp(M, elem, view(h, :, k)) for M in Cprep]
                            Cgp = mats[1]
                            for i in 2:length(mats)
                                Cgp = matmul_sf(Cgp, mats[i])
                            end
                        else
                            Cgp = _eval_coefficient_at_gp(Cprep, elem, view(h, :, k))
                        end
                        if Cgp isa Number && weight === nothing

                            w = jacDet[k] * intWeights[k] * Cgp

                            if op_u isa AxialGradOp || op_u isa TangentialGradOp
                                invJk = invJac[1:Pu.dim, 3k-2:3k-(3-Pu.dim)]
                                Jk = inv(invJk')
                                t = Jk[:,1]
                                t = t / norm(t)
                                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
                            else
                                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes)
                            end
                            if op_s isa AxialGradOp || op_s isa TangentialGradOp
                                invJk = invJac[1:Ps.dim, 3k-2:3k-(3-Ps.dim)]
                                Jk = inv(invJk')
                                t = Jk[:,1]
                                t = t / norm(t)
                                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
                            else
                                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes)
                            end

                            mul!(Ke, transpose(Bs), Bu, w, 1.0)

                        elseif Cgp isa AbstractMatrix && weight === nothing

                            w = jacDet[k] * intWeights[k]

                            if op_u isa AxialGradOp
                                invJk = invJac[1:Pu.dim, 3k-2:3k-(3-Pu.dim)]
                                Jk = inv(invJk')
                                t = Jk[:,1]
                                t = t / norm(t)
                                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes, t)
                            else
                                build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes)
                            end
                            if op_s isa AxialGradOp
                                invJk = invJac[1:Ps.dim, 3k-2:3k-(3-Ps.dim)]
                                Jk = inv(invJk')
                                t = Jk[:,1]
                                t = t / norm(t)
                                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes, t)
                            else
                                build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes)
                            end

                            mul!(tmp, transpose(Bs), Cgp) # Here Cgp was transposed
                            mul!(Ke, tmp, Bu, w, 1.0)

                        else
                            error("assemble_operator error")
                        end
                    end


                    # scatter Ke(s,u) -> global IJV
                    @inbounds for a_loc in 1:ndofs_s_loc
                        node_a = div(a_loc - 1, Ps.pdim) + 1
                        comp_a = mod(a_loc - 1, Ps.pdim) + 1
                        Ia_node = nnet[e, node_a]
                        Ia = (Ia_node - 1) * Ps.pdim + comp_a
                
                        @inbounds for b_loc in 1:ndofs_u_loc
                            node_b = div(b_loc - 1, Pu.pdim) + 1
                            comp_b = mod(b_loc - 1, Pu.pdim) + 1
                            Jb_node = nnet[e, node_b]
                            Jb = (Jb_node - 1) * Pu.pdim + comp_b
                            
                            #@assert pos >= 1
                            #@assert pos - 1 <= length(I)
                            #@assert all(isfinite, V[1:pos-1])
                            if pos >= length(I)
                                newlen = Int(ceil(1.5*length(I))) + 1000
                                resize!(I, newlen)
                                resize!(J, newlen)
                                resize!(V, newlen)
                            end
                
                            I[pos] = Ia
                            J[pos] = Jb
                            V[pos] = Ke[a_loc, b_loc]
                            pos += 1
                        end
                    end
                end
            end
        end
    end

    resize!(I, pos - 1)
    resize!(J, pos - 1)
    resize!(V, pos - 1)
    K = sparse(I, J, V, ndofs(Ps), ndofs(Pu))
    dropzeros!(K)
    @assert ndofs(Ps) > 0
    @assert ndofs(Pu) > 0
    @assert maximum(I[1:pos-1]) <= ndofs(Ps)
    @assert maximum(J[1:pos-1]) <= ndofs(Pu)
    @assert minimum(I[1:pos-1]) >= 1
    @assert minimum(J[1:pos-1]) >= 1
    return SystemMatrix(K, Pu, Ps)
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

"""
    assemble_linear(P::Problem, op, rhs; weight=nothing, domain)

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

Returns
-------
Global right-hand side vector.
"""
function assemble_linear(
    P::Problem,
    op::AbstractOp,
    g;
    weight = nothing,
    domain = nothing)

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
    rhs = zeros(nd)

    phNames = domain === nothing ?
        [mat.phName for mat in P.material] :
        [domain.name]
    
    for phName in phNames
        Gprep = domain === nothing ?
            _prepare_coefficient(g, DomainSpec(:Ω, phName)) :
            _prepare_coefficient(g, domain)

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for (edim, etag) in dimTags

            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for itype in eachindex(elemTypes)

                et = elemTypes[itype]

                _, _, order, numNodes, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                intPoints, intWeights =
                    gmsh.model.mesh.getIntegrationPoints(
                        et, "Gauss" * string(2order + 1)
                    )

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

                ndofs_loc = P.pdim * numNodes
                B = zeros(outdim, ndofs_loc)
                fe = zeros(ndofs_loc)

                nel = length(elemTags[itype])

                for e in 1:nel
                    fill!(fe, 0.0)

                    elem = elemTags[itype][e]

                    jac, jacDet, _ =
                        gmsh.model.mesh.getJacobian(elem, intPoints)

                    Jac = reshape(jac, 3, :)
                    invJac = zeros(3, 3 * numIntPoints)
                    ∂h = zeros(P.dim, numNodes * numIntPoints)

                    @inbounds for k in 1:numIntPoints
                        invJac[1:3, 3k-2:3k] .= inv(Jac[1:3, 3k-2:3k])'
                    end

                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numIntPoints, a in 1:numNodes
                        invJk = invJac[1:P.dim, 3k-2:3k-(3-P.dim)]
                        gha = ∇h[a*3-2:a*3-(3-P.dim), k]
                        ∂h[1:P.dim, (k-1)*numNodes+a] .= invJk * gha
                    end

                    for k in 1:numIntPoints

                        w = jacDet[k] * intWeights[k]

                        build_B!(B, op, P, k, h, ∂h, numNodes)

                        g_gp = _eval_coefficient_at_gp(
                            Gprep,
                            elem,
                            view(h, :, k)
                        )

                        gvec = g_gp isa Number ? [g_gp] : g_gp

                        fe .+= (transpose(B) * gvec) * w
                    end

                    # scatter local element vector -> global rhs
                    for a_loc in 1:ndofs_loc
                        node = div(a_loc - 1, P.pdim) + 1
                        comp = mod(a_loc - 1, P.pdim) + 1

                        gnode = elemNodeTags[itype][(e-1)*numNodes + node]
                        row = (gnode - 1) * P.pdim + comp

                        rhs[row] += fe[a_loc]
                    end
                end
            end
        end
    end

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

    return tuple(results...)
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

# Example
```julia
K = ∫(AxialGrad(Pu) ⋅ (E*A) ⋅ AxialGrad(Pu); Ω="truss")
"""
AxialGrad(P) = OpApplied(P, AxialGradOp())

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

struct BilinearTerm
    a::OpApplied
    coef
    b::OpApplied
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

import Base: +, -, *

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


*(c::Union{Number,ScalarField}, t::BilinearTerm) =
    BilinearTerm(t.a, c, t.b)

*(t::BilinearTerm, c::Union{Number,ScalarField}) =
    BilinearTerm(t.a, c, t.b)

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

function _assemble(term::WeakTerm{BilinearTerm}, dom, weight=nothing)

    Pu = term.term.a.P

    if dom !== nothing
        gmsh.model.setCurrent(Pu.name)
        _check_domain_dim(Pu, dom)
    end

    coef = term.term.coef

    assemble_operator(
        term.term.a.P,
        term.term.a.op,
        term.term.b.P,
        term.term.b.op,
        coefficient=coef, #isa AbstractMatrix ? 1.0 : coef,
        weight=nothing, #coef isa AbstractMatrix ? coef : nothing,
        domain=dom
    )

end

function _assemble(expr::WeakSum, dom, weight)

    _assemble(expr.a, dom, weight) + _assemble(expr.b, dom, weight)

end

function _assemble(term::WeakTerm{LinearTerm}, dom, weight=nothing)

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
        domain = dom
    )

end

# ------------------------------------------------------------
# Integral
# ------------------------------------------------------------

"""
    ∫(expr::WeakExpr; Ω=nothing, Γ=nothing, weight=nothing)

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

# Arguments

`expr`

Weak-form expression composed of operators and coefficients.

# Keyword arguments

`Ω`

Volume physical group name.

`Γ`

Boundary physical group name.

# Returns

`SystemMatrix` or `ScalarField`, `VectorField`, `TensorField`
"""
function ∫(expr::WeakExpr; Ω=nothing, Γ=nothing, weight=nothing)

    _check_scalarfields(expr)

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    return _assemble(expr, dom, weight)

end

function ∫(t::BilinearTerm; Ω=nothing, Γ=nothing, weight=nothing)

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    Pu = t.a.P

    if dom !== nothing
        gmsh.model.setCurrent(Pu.name)
        _check_domain_dim(Pu, dom)
    end

    return assemble_operator(
        t.b.P,
        t.b.op,
        t.a.P,
        t.a.op;
        coefficient = t.coef,
        domain = dom,
        weight = weight
    )

end

function ∫(a::OpApplied, b::OpApplied; Ω=nothing, Γ=nothing, weight=nothing)

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    Pu = a.P

    if dom !== nothing
        gmsh.model.setCurrent(Pu.name)
        _check_domain_dim(Pu, dom)
    end

    return assemble_operator(
        b.P,
        b.op,
        a.P,
        a.op;
        coefficient = 1.0,
        domain = dom,
        weight = weight
    )

end

function ∫(t::LinearTerm; Ω=nothing, Γ=nothing, weight=nothing)

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
        weight = weight
    )

end

function ∫(mc::MatrixChain; Ω=nothing, Γ=nothing, weight=nothing)

    return ∫(LinearTerm(mc); Ω=Ω, Γ=Γ, weight=weight)

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
∫Ω(name, expr) = ∫(expr; Ω=name)

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
∫Γ(name, expr) = ∫(expr; Γ=name)

const ε = SymGrad

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
    nsteps::Int=1
)
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

    return tuple(results...)
end

"""
    FDM(K::SystemMatrix,
        C::SystemMatrix,
        q::SystemVector,
        bc::Vector{BoundaryCondition},
        X0::SystemVector,
        n::Int,
        Δt::Float64;
        ϑ=0.5)

Alias: FDM(K, C, q, X0, n, Δt; ϑ=0.5, support=Vector{BoundaryCondition}())

Time integration for a multifield first-order system

    C * ẋ + K * x = q

using the theta-method.

Arguments
- `K`: system matrix
- `C`: capacity / mass-like matrix
- `q`: global RHS vector as `SystemVector`
- `bc`: Dirichlet boundary conditions
- `X0`: initial state as `SystemVector`
- `n`: number of stored time steps
- `Δt`: time step size

Keyword arguments
- `ϑ`: theta parameter
    - `0.0` explicit Euler
    - `0.5` Crank-Nicolson
    - `1.0` implicit Euler

Returns
Tuple of fields in the block order of `K.problems`.
"""
function FDM(
    K::SystemMatrix,
    C::SystemMatrix,
    q::SystemVector,
    bc::Vector{BoundaryCondition},
    X0::SystemVector,
    n::Int,
    Δt::Float64;
    ϑ=0.5
)

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

    # ------------------------------------------------------------------
    # 2) BC data
    # ------------------------------------------------------------------
    free, fix, xD = multifield_bc_data(K, bc; nsteps=n)

    # ------------------------------------------------------------------
    # 3) Allocate history
    # ------------------------------------------------------------------
    X = zeros(ndof, n)
    t = zeros(n)

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

    if !isempty(fix)
        Kfc = K0[free, fix]
        Cfc = C0[free, fix]
    else
        Kfc = zeros(length(free), 0)
        Cfc = zeros(length(free), 0)
    end

    # explicit shortcut only if Cff is diagonal and theta = 0
    is_diag_Cff = nnz(Cff) == length(diag(Cff))

    # ------------------------------------------------------------------
    # 5) Time stepping
    # ------------------------------------------------------------------
    if ϑ == 0 && is_diag_Cff
        invCff = spdiagm(1.0 ./ diag(Cff))

        for i in 2:n
            qi = size(q.a, 2) == 1 ? 1 : i
            qn = q.a[:, qi]

            xc_n = isempty(fix) ? zeros(0) : xD[fix, i-1]
            xc_np1 = isempty(fix) ? zeros(0) : xD[fix, i]

            xfree_n = @view X[free, i-1]

            rhs =
                qn[free] -
                Kff * xfree_n -
                Kfc * xc_n -
                Cfc * ((xc_np1 - xc_n) ./ Δt)

            xfree_np1 = xfree_n + Δt .* (invCff * rhs)

            X[free, i] .= xfree_np1
            if !isempty(fix)
                X[fix, i] .= xc_np1
            end

            t[i] = t[i-1] + Δt
        end

    else
        A = Cff + ϑ * Δt * Kff
        B = Cff - (1 - ϑ) * Δt * Kff
        luA = lu(A)

        for i in 2:n
            qi = size(q.a, 2) == 1 ? 1 : i
            mi = qi == 1 ? 0 : 1
            qn = q.a[:, qi-mi]
            qnp1 = q.a[:, qi]
            @views qth = (1 - ϑ) .* qn[free] .+ ϑ .* qnp1[free]

            xc_n = isempty(fix) ? zeros(0) : xD[fix, i-1]
            xc_np1 = isempty(fix) ? zeros(0) : xD[fix, i]

            rhs =
                B * (@view X[free, i-1]) +
                Δt .* qth

            if !isempty(fix)
                rhs .-= (Cfc + ϑ * Δt * Kfc) * xc_np1
                rhs .+= (Cfc - (1 - ϑ) * Δt * Kfc) * xc_n
            end

            X[free, i] .= luA \ rhs
            if !isempty(fix)
                X[fix, i] .= xc_np1
            end

            t[i] = t[i-1] + Δt
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

