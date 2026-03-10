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

export ∫
export ∫Ω
export ∫Γ

export solveField

export ε

export ⋅


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

#"""
#    ndofs(problem::Problem)
#
#Return total number of dofs for a single-field problem.
#"""
#ndofs(P::Problem) = P.non * P.pdim

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
    @assert P.pdim == P.dim "SymGradOp requires vector field with pdim == dim."
    return (P.dim == 2) ? 3 : 6  # engineering strain components
end

function op_outdim(::TensorDivOp, P::Problem)
    @assert P.pdim == P.dim^2 "TensorDivOp requires pdim == dim^2 (full 2nd-order tensor)."
    return P.dim
end

function op_outdim(op::AdvOp, P::Problem)
    @assert P.pdim == 1  # scalar field
    return 1
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
function build_B!(B::AbstractMatrix, ::IdOp, P::Problem, k::Int, h, ∂h, numNodes::Int)
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

function build_B!(B::AbstractMatrix, ::GradOp, P::Problem, k::Int, h, ∂h, numNodes::Int)
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

function build_B!(B::AbstractMatrix, ::DivOp, P::Problem, k::Int, h, ∂h, numNodes::Int)
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
    P::Problem, k::Int, h, ∂h, numNodes::Int)
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
    P::Problem, k::Int, h, ∂h, numNodes::Int)
    fill!(B, 0.0)
    dim = P.dim
    pdim = P.pdim
    @assert pdim == dim

    if dim == 2
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
    P::Problem, k::Int, h, ∂h, numNodes::Int)
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
    P::Problem, k::Int, h, ∂h, numNodes::Int)

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

@inline function _build_elemwise_coeff_dict(coef::ScalarField)
    p = nodesToElements(coef)
    return Dict(zip(p.numElem, p.A))  # elemTag => coeff nodal vector(s)
end

"""
    _prepare_coefficient(C)

Prepare coefficient for assembly.

Returns one of:

Float64
Dict(elem => nodal array)
Matrix{Any} with entries Float64 or Dict(elem => nodal array)
"""
function _prepare_coefficient(C)

    # scalar constant
    if C isa Number
        return Float64(C)
    #end

    # scalar field
    elseif C isa ScalarField
        return _build_elemwise_coeff_dict(C)
    #end

    # tensor coefficient
    elseif C isa AbstractMatrix

        W = Matrix{Any}(undef, size(C)...)

        for I in CartesianIndices(C)

            cij = C[I]

            if cij isa Number
                W[I] = Float64(cij)

            elseif cij isa ScalarField
                W[I] = _build_elemwise_coeff_dict(cij)

            else
                error("Matrix coefficient entries must be Number or ScalarField, got $(typeof(cij))")
            end

        end

        return W
    #end

    elseif C isa AbstractVector
        mats = [_prepare_coefficient(M) for M in C]
        return mats
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
    out_u = op_outdim(op_u, Pu)
    out_s = op_outdim(op_s, Ps)
    @assert out_u == out_s "Operator output dims mismatch: $out_u vs $out_s."

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

    Cprep = _prepare_coefficient(coefficient)

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

                            build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes)
                            build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes)

                            mul!(Ke, transpose(Bs), Bu, w, 1.0)

                        elseif Cgp isa AbstractMatrix && weight === nothing

                            w = jacDet[k] * intWeights[k]

                            build_B!(Bu, op_u, Pu, k, h, ∂h, numNodes)
                            build_B!(Bs, op_s, Ps, k, h, ∂h, numNodes)

                            mul!(tmpBu, Cgp, Bu)
                            mul!(Ke, transpose(Bs), tmpBu, w, 1.0)

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
                
                            if pos > length(I)
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
    return SystemMatrix(K, Pu, Ps)
end

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
               f::Union{ScalarField,VectorField}; 
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

    ndof = length(b)

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

    x = zeros(ndof)
    xD = zeros(ndof)

    # fill xD
    for bc in support
        P = bc.problem
        idx = findfirst(q -> q === P, problems)
        offset = offsets[idx]

        x_local = applyBoundaryConditions(P, [bc])
        local_dofs = constrainedDoFs(P, [bc])

        xD[offset.+local_dofs] .= x_local.a[local_dofs, 1]
    end

    fixed = unique(vcat([offsets[findfirst(q -> q === bc.problem, problems)] .+
                         constrainedDoFs(bc.problem, [bc])
                         for bc in support]...))

    free = setdiff(1:ndof, fixed)

    A_ff = A[free, free]
    b_f = b[free] - A[free, fixed] * xD[fixed]

    x[free] = A_ff \ b_f
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
            results[i] = ScalarField([], reshape(xloc, :, 1), [0], [], 1, :scalar, P)

        elseif P.pdim == 2 || P.pdim == 3
            type = P.pdim == 2 ? :v2D : :v3D
            results[i] = VectorField([], reshape(xloc, :, 1), [0], [], 1, type, P)

        elseif P.pdim == 9
            results[i] = TensorField([], reshape(xloc, :, 1), [0], [], 1, :tensor, P)

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

struct DomainSpec
    kind::Symbol     # :Ω vagy :Γ
    name::String
end

Base.show(io::IO, d::DomainSpec) =
    print(io, "$(d.kind)=\"$(d.name)\"")

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

"""
    ⋅(a,b)

Weak-form inner product operator used in the DSL.

Constructs bilinear forms of the type

    (Op_s v) ⋅ C ⋅ (Op_u u)

Supported forms

Operator pair

    Grad(Pu) ⋅ Grad(Pu)

Operator with tensor coefficient

    SymGrad(Pu) ⋅ C ⋅ SymGrad(Pu)

Matrix chain

    Grad(Pu) ⋅ A ⋅ B ⋅ C ⋅ Grad(Pu)

where matrices may contain

    Number
    ScalarField

entries.

The matrix chain is evaluated at Gauss points during assembly.
"""
⋅(a::OpApplied, b::OpApplied) = BilinearTerm(a, 1.0, b)

⋅(a::OpApplied,
   C::AbstractMatrix,
   b::OpApplied) =
    begin
        _check_coeff_matrix(C)
        BilinearTerm(a, C, b)
    end

"""
Internal representation of chained tensor coefficients.

Constructed automatically by expressions such as

    Grad(Pu) ⋅ A ⋅ B ⋅ C ⋅ Grad(Pu)

The matrices are stored in `mats` and multiplied during
assembly at Gauss points.

Users normally never construct this type directly.
"""
struct MatrixChain
    a::OpApplied
    mats::Vector{Any}
end

function ⋅(a::OpApplied,
    C::AbstractMatrix)
    _check_coeff_matrix(C)
    return MatrixChain(a, [C])
end

function ⋅(t::MatrixChain, C::AbstractMatrix)
    _check_coeff_matrix(C)
    push!(t.mats, C)
    return t
end

function ⋅(t::MatrixChain, b::OpApplied)
    return BilinearTerm(t.a, t.mats, b)
end

⋅(a::OpApplied, P::Problem) =
    BilinearTerm(a, 1.0, Id(P))

⋅(P::Problem, b::OpApplied) =
    BilinearTerm(Id(P), 1.0, b)

⋅(a::OpApplied, C::Union{Number,ScalarField}, P::Problem) =
    BilinearTerm(a, C, Id(P))

⋅(P::Problem, C::Union{Number,ScalarField}, b::OpApplied) =
    BilinearTerm(Id(P), C, b)

*(c::Union{Number,ScalarField}, P::Problem) =
    c * Id(P)

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
struct WeakTerm <: WeakExpr
    coef::Number
    term::BilinearTerm
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
promote_term(t::WeakTerm) = t


###############################################################
# Scalar multiplication
###############################################################


*(c::Union{Number,ScalarField}, t::BilinearTerm) =
    BilinearTerm(t.a, c, t.b)

*(t::BilinearTerm, c::Union{Number,ScalarField}) =
    BilinearTerm(t.a, c, t.b)

*(c::Number, t::WeakTerm) = WeakTerm(c * t.coef, t.term)
*(t::WeakTerm, c::Number) = WeakTerm(c * t.coef, t.term)


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

function _assemble(term::WeakTerm, dom, weight=nothing)

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

# ------------------------------------------------------------
# Integral
# ------------------------------------------------------------

"""
    ∫(expr::WeakExpr; Ω=nothing, Γ=nothing, weight=nothing)

Assemble a finite element operator from a weak-form expression.

Examples

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

# Arguments

`expr`

Weak-form expression composed of operators and coefficients.

# Keyword arguments

`Ω`

Volume physical group name.

`Γ`

Boundary physical group name.

# Returns

`SystemMatrix`
"""
function ∫(expr::WeakExpr; Ω=nothing, Γ=nothing, weight=nothing)

    _check_scalarfields(expr)

    dom = _domain_spec(; Ω=Ω, Γ=Γ)

    return _assemble(expr, dom, weight)

end

# ------------------------------------------------------------
# Fallback forms
# ------------------------------------------------------------

#function ∫(t::BilinearTerm; coef=1.0, Ω=nothing, Γ=nothing, weight=nothing)
#
#    dom = _domain_spec(; Ω=Ω, Γ=Γ)
#
#    assemble_operator(
#        t.a.P,
#        t.a.op,
#        t.b.P,
#        t.b.op,
#        coefficient=coef,
#        weight=weight,
#        domain=dom
#    )
#
#end

#function ∫(a::OpApplied, b::OpApplied; coef=1.0, Ω=nothing, Γ=nothing, weight=nothing)
#
#    dom = _domain_spec(; Ω=Ω, Γ=Γ)
#
#    assemble_operator(
#        a.P,
#        a.op,
#        b.P,
#        b.op,
#        coefficient=coef,
#        weight=weight,
#        domain=dom
#    )
#
#end

function ∫(t::BilinearTerm; coef=1.0, Ω=nothing, Γ=nothing, weight=nothing)
    dom = _domain_spec(; Ω=Ω, Γ=Γ)
    Pu = t.a.P
    if dom !== nothing
        gmsh.model.setCurrent(Pu.name)
        _check_domain_dim(Pu, dom)
    end
    return assemble_operator(t.a.P, t.a.op, t.b.P, t.b.op;
        coefficient=t.coef, domain=dom, weight=nothing)
end

function ∫(a::OpApplied, b::OpApplied; coef=1.0, Ω=nothing, Γ=nothing, weight=nothing)
    dom = _domain_spec(; Ω=Ω, Γ=Γ)
    Pu = a.P
    if dom !== nothing
        gmsh.model.setCurrent(Pu.name)
        _check_domain_dim(Pu, dom)
    end
    return assemble_operator(a.P, a.op, b.P, b.op;
        coefficient=coef, domain=dom, weight=nothing)
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