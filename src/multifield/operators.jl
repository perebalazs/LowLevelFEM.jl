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

The sum is expanded into ordinary bilinear terms. With the default direct CSC
assembly, one structural pattern and one set of worker value buffers are reused
for all expanded terms belonging to the same test/trial field pair. If a
compound expression contains multiple fields, terms are grouped by `Problem`
and assembled automatically into a multifield block `SystemMatrix`; each block
owns and reuses its own CSC pattern. Each term may carry its own quadrature rule
via `full(...)` and `reduced(...)`.

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

export consistentToLumped
export rigidRotationMap, collapseMPC


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

# Operators may request preservation of structural CSC zeros so a matrix
# returned by `∫` can be reused safely through `csc_matrix=...`.
_preserve_csc_pattern(::AbstractOp) = false

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
      γyz
      γzx ]
"""
struct SymGradOp <: AbstractOp end

"""
    TensorDivOp()

Divergence of a second-order tensor field.

Represents

    div(σ)

where σ is stored as a full 3×3 tensor in column-major order, with

    pdim = 9

The output is a vector with `dim` components. In 2D, only the in-plane
tensor block and the x/y derivatives contribute.
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

Internal scalar axial derivative/strain operator along the line-element tangent.
For scalar fields return `t ⋅ grad(s)`; for vector fields return
`t ⋅ grad(u) ⋅ t`. Construct through [`AxialGrad`](@ref).
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
    @assert P.dim in (2, 3) && P.pdim == 9 "TensorDivOp requires a 2D or 3D tensor field with nine stored components."
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
        # rows: [εxx, εyy, εzz, γxy, γyz, γzx]
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

            B[5, coly] = dNz          # γyz
            B[5, colz] = dNy

            B[6, colx] = dNz          # γzx
            B[6, colz] = dNx
        end
    end

    return B
end

function build_B!(B::AbstractMatrix, ::TensorDivOp,
    P::Problem, k::Int, h, ∂h, numNodes::Integer)
    fill!(B, 0.0)
    dim = P.dim
    pdim = P.pdim  # full 3×3 tensor, also in 2D

    # Column-major tensor storage: σ_ij maps to i + 3*(j - 1).
    @inbounds for a in 1:numNodes
        for i in 1:dim
            row = i
            for j in 1:dim
                col = (a - 1) * pdim + i + 3 * (j - 1)
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

    pdim = P.pdim
    dim = P.dim

    if pdim == 1
        # Scalar field:
        # AxialGrad(φ) = ∇φ ⋅ t
        @inbounds for a in 1:numNodes
            val = 0.0

            for d in 1:dim
                val += t[d] * ∂h[d, (k-1)*numNodes + a]
            end

            B[1, a] = val
        end

    elseif pdim == dim
        # Vector field:
        # AxialGrad(u) = tᵀ ⋅ ∇u ⋅ t
        @inbounds for a in 1:numNodes
            for c in 1:pdim
                col = (a - 1) * pdim + c

                val = 0.0

                for d in 1:dim
                    val += t[c] * t[d] *
                           ∂h[d, (k-1)*numNodes + a]
                end

                B[1, col] = val
            end
        end

    else
        error(
            "AxialGradOp supports scalar fields or vector fields with pdim == dim."
        )
    end

    return B
end
#=
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
=#

function build_B!(B::AbstractMatrix, ::TangentialGradOp,
                  P::Problem, k::Int, h, ∂h, numNodes::Integer, t)

    fill!(B, 0.0)

    pdim = P.pdim
    dim  = P.dim

    if pdim == 1

        # Scalar field:
        # TangentialGrad(φ) = ∇φ ⋅ t
        @inbounds for a in 1:numNodes
            val = 0.0

            for d in 1:dim
                val += t[d] * ∂h[d, (k-1)*numNodes + a]
            end

            B[1, a] = val
        end

    elseif pdim == dim

        # Vector field:
        # TangentialGrad(u) = (∇u) ⋅ t
        @inbounds for a in 1:numNodes

            val = 0.0

            for d in 1:dim
                val += t[d] * ∂h[d, (k-1)*numNodes + a]
            end

            for c in 1:pdim
                col = (a - 1)*pdim + c
                B[c, col] = val
            end
        end

    else

        error(
            "TangentialGradOp supports scalar fields or vector fields with pdim == dim."
        )
    end

    return B
end
#=
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
=#

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

