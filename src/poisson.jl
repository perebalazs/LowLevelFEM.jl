#export stiffnessMatrixPoisson, advectionMatrix, massMatrixPoisson
export poissonMatrix, reactionMatrix, advectionMatrix
export ∫∇oN_c_∇oN_dΩ, ∫∇N_c_∇N_dΩ, ∫∇xN_c_∇xN_dΩ, ∫N_c_N_dΩ, ∫N_c_∂N∂x_dΩ, ∫N_c_∂N∂x_dΩ, ∫N_c_∂N∂y_dΩ, ∫N_c_∂N∂z_dΩ
export ∫N_c_dΩ
export gradDivMatrix, symmetricGradientMatrix, curlCurlMatrix
export gradMatrix, navierStokesAdvectionMatrix
export tensorLaplaceMatrix, traceLaplaceMatrix, beltramiMichellMatrix, tensorDivDivMatrix
export sourceVector, loadTensor
export solveField

"""
    stiffnessMatrixPoissonAllInOne(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global stiffness matrix of a Poisson-type problem using the finite
element method.

The matrix corresponds to the weak form of the Poisson operator
```

K_ab = ∫_Ω (∇N_a · ∇N_b) · α(x) dΩ,

```
where `N_a` are the Lagrange shape functions and `α(x)` is a scalar coefficient.

# Arguments
- `problem::Problem`:
  Finite element problem definition, including geometry, dimension, physical
  groups, and discretization.
- `coefficient::Union{Number,ScalarField}`:
  Scalar coefficient `α(x)` in the Poisson equation.
  - If a `Number`, a constant coefficient is used over the entire domain.
  - If a `ScalarField`, the coefficient is given elementwise and interpolated
    to Gauss points.

# Returns
- `SystemMatrix`:
  The assembled global stiffness matrix associated with the Poisson problem.

# Notes
- The function assembles only the left-hand side operator of the Poisson equation.
- The spatial dimension (2D or 3D) is taken from `problem`.
- Boundary conditions and the right-hand side vector are handled separately.

Everything in one big function.
"""

function stiffnessMatrixPoissonAllInOne(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    #lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
    lengthOfIJV = estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)

    if coefficient isa Number
        p = coefficient
    else
        p = nodesToElements(coefficient)
        pa = Dict(zip(p.numElem, p.A))
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dim = problem.dim
        pdim = problem.pdim

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                #H = zeros(pdim * numIntPoints, pdim * numNodes)
                #@inbounds for k in 1:numIntPoints, l in 1:numNodes
                #    val = h[(k-1)*numNodes + l]
                #    for kk in 1:pdim
                #        row = (k-1)*pdim + kk
                #        col = (l-1)*pdim + kk
                #        H[row, col] = val
                #    end
                #end
                @inbounds for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    @inbounds for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    @inbounds for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] .= inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    fill!(K1, 0.0)
                    if p isa ScalarField
                        pa0 = pa[elem][:,1]
                        @inbounds for k in 1:numIntPoints
                            val = dot(pa0, h[:, k])
                            w = jacDet[k] * intWeights[k] * val
                            @inbounds for a in 1:numNodes
                                @inbounds for b in 1:numNodes
                                    s = 0.0
                                    @inbounds for d in 1:dim
                                        s += ∂h[d, (k-1)*numNodes + a] *
                                            ∂h[d, (k-1)*numNodes + b]
                                    end
                                    K1[a, b] += s * w
                                end
                            end
                        end
                    else
                        @inbounds for k in 1:numIntPoints
                            w = jacDet[k] * intWeights[k] * p
                            @inbounds for a in 1:numNodes
                                @inbounds for b in 1:numNodes
                                    s = 0.0
                                    @inbounds for d in 1:dim
                                        s += ∂h[d, (k-1)*numNodes + a] *
                                            ∂h[d, (k-1)*numNodes + b]
                                    end
                                    K1[a, b] += s * w
                                end
                            end
                        end
                    end
                    for a in 1:numNodes
                        Ia = nnet[j, a]
                        for b in 1:numNodes
                            I[pos] = Ia
                            J[pos] = nnet[j, b]
                            V[pos] = K1[a, b]
                            pos += 1
                        end
                    end
                end
            end
        end
    end
    resize!(I, pos-1)
    resize!(J, pos-1)
    resize!(V, pos-1)
    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end













# --- internal helpers ---------------------------------------------------------

@inline function _build_elemwise_coeff_dict(coefficient::ScalarField)
    p = nodesToElements(coefficient)               # elementwise
    return p, Dict(zip(p.numElem, p.A))            # elemTag => nodal coeff vector(s)
end

#@inline function _coeff_at_gp(::Number, ::Any, ::Any, ::Any) # const
#    return nothing
#end

@inline function _coeff_at_gp(pa::Dict{<:Integer, <:AbstractMatrix}, elem::Integer, hcol::AbstractVector)
    # in your original code: pa0 = pa[elem][:,1]; val = dot(pa0, h[:,k])
    return dot(view(pa[elem], :, 1), hcol)
end

@inline function _add_blockdiag!(Ke::Matrix{Float64}, pdim::Int, a::Int, b::Int, val::Float64)
    # local dof ordering: (node-1)*pdim + comp
    @inbounds for kk in 1:pdim
        ia = (a - 1) * pdim + kk
        ib = (b - 1) * pdim + kk
        Ke[ia, ib] += val
    end
    return nothing
end

function _assemble_poissonlike!(
    I::Vector{Int}, J::Vector{Int}, V::Vector{Float64},
    pos0::Int,
    problem::Problem,
    phName::String,
    coefficient::Union{Number,ScalarField},
    kernel!::Function;
    dir::Int = 1,                    # used by convection kernel (x=1,y=2,z=3)
    )
    gmsh.model.setCurrent(problem.name)

    pdim = problem.pdim
    dim  = problem.dim

    # coefficient preparation (same idea as your code)
    pconst = 0.0
    pa = Dict{Int, Any}()
    if coefficient isa Number
        pconst = Float64(coefficient)
    else
        _, pa = _build_elemwise_coeff_dict(coefficient)
    end

    pos = pos0

    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
    for idm in 1:length(dimTags)
        dimTag = dimTags[idm]
        edim = dimTag[1]
        etag = dimTag[2]

        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

        for itype in 1:length(elemTypes)
            et = elemTypes[itype]
            elementName, _, order, numNodes::Int64, localNodeCoord, numPrimaryNodes =
                gmsh.model.mesh.getElementProperties(et)

            # keep your integration order choice: "Gauss"*(2*order+1)
            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
            numIntPoints = length(intWeights)

            # basis and grad basis
            comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h = reshape(fun, :, numIntPoints)  # (numNodes, numIntPoints) effectively

            comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
            ∇h = reshape(dfun, :, numIntPoints)  # (3*numNodes, numIntPoints)

            # element connectivity and geometry buffers
            nnet   = zeros(Int, length(elemTags[itype]), numNodes)
            invJac = zeros(3, 3numIntPoints)
            ∂h     = zeros(dim, numNodes * numIntPoints)
            Ke     = zeros(pdim * numNodes, pdim * numNodes)

            # build connectivity table (same as your code)
            @inbounds for j in 1:length(elemTags[itype])
                for a in 1:numNodes
                    nnet[j, a] = elemNodeTags[itype][(j-1)*numNodes + a]
                end
            end

            # loop elements
            @inbounds for j in 1:length(elemTags[itype])
                elem = elemTags[itype][j]

                jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                Jac = reshape(jac, 3, :)

                # inv(J)' at each GP (same approach)
                @inbounds for k in 1:numIntPoints
                    invJac[1:3, 3*k-2:3*k] .= inv(Jac[1:3, 3*k-2:3*k])'
                end

                # compute physical gradients of basis: ∂h
                fill!(∂h, 0.0)
                @inbounds for k in 1:numIntPoints, a in 1:numNodes
                    # take only needed rows (dim=2 or 3)
                    invJk = invJac[1:dim, 3*k-2:3*k- (3-dim)]
                    # GradLagrange returns 3 components always; slice accordingly
                    gha  = ∇h[a*3-2 : a*3-(3-dim), k]
                    ∂h[1:dim, (k-1)*numNodes + a] .= invJk * gha
                end

                fill!(Ke, 0.0)

                # integrate
                if coefficient isa ScalarField
                    @inbounds for k in 1:numIntPoints
                        valc = _coeff_at_gp(pa, elem, view(h, :, k))
                        w    = jacDet[k] * intWeights[k] * valc
                        kernel!(Ke, w, k, h, ∂h, numNodes, pdim, dim, elem; dir=dir)
                    end
                else
                    @inbounds for k in 1:numIntPoints
                        w = jacDet[k] * intWeights[k] * pconst
                        kernel!(Ke, w, k, h, ∂h, numNodes, pdim, dim, elem; dir=dir)
                    end
                end

                # scatter to global I,J,V (blockdiag already in Ke)
                @inbounds for a in 1:(pdim*numNodes)
                    # local dof -> node
                    na = (div(a-1, pdim) + 1)
                    Ia_node = nnet[j, na]
                    Ia = (Ia_node - 1) * pdim + (mod(a-1, pdim) + 1)

                    @inbounds for b in 1:(pdim*numNodes)
                        nb = (div(b-1, pdim) + 1)
                        Jb_node = nnet[j, nb]
                        Jb = (Jb_node - 1) * pdim + (mod(b-1, pdim) + 1)

                        I[pos] = Ia
                        J[pos] = Jb
                        V[pos] = Ke[a, b]
                        pos += 1
                    end
                end
            end
        end
    end

    return pos
end

# --- kernels ------------------------------------------------------------------

# Diffusion / stiffness: ∫ (∇Na·∇Nb) * w dA
@inline function _kernel_stiffness!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @inbounds for a in 1:numNodes
        @inbounds for b in 1:numNodes
            s = 0.0
            @inbounds for d in 1:dim
                s += ∂h[d, (k-1)*numNodes + a] * ∂h[d, (k-1)*numNodes + b]
            end
            _add_blockdiag!(Ke, pdim, a, b, s * w)
        end
    end
    return nothing
end

# Convection: ∫ Na * (∂Nb/∂x_dir) * w dA   (dir=1 => x)
@inline function _kernel_convection!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @inbounds for a in 1:numNodes
        Na = h[(k-1)*numNodes + a]
        @inbounds for b in 1:numNodes
            dNb = ∂h[dir, (k-1)*numNodes + b]
            _add_blockdiag!(Ke, pdim, a, b, Na * dNb * w)
        end
    end
    return nothing
end

# Mass / reaction: ∫ Na * Nb * w dA
@inline function _kernel_mass!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @inbounds for a in 1:numNodes
        Na = h[(k-1)*numNodes + a]
        @inbounds for b in 1:numNodes
            Nb = h[(k-1)*numNodes + b]
            _add_blockdiag!(Ke, pdim, a, b, Na * Nb * w)
        end
    end
    return nothing
end

# Grad-div: ∫ (div Na_vec) * (div Nb_vec) * w dΩ
#
# IMPORTANT:
#   - This assumes a standard vector unknown u with pdim == dim (e.g., 2D -> (ux,uy), 3D -> (ux,uy,uz)).
#   - Local dof ordering already matches _add_blockdiag! convention (block diagonal by component),
#     but grad-div COUPLES components, so we must assemble the full pdim×pdim block per node pair.
#
@inline function _kernel_graddiv!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @assert pdim == dim "_kernel_graddiv!: requires pdim == dim (vector field components match spatial dim)"

    # ∂h is stored as ∂h[d, (k-1)*numNodes + a] = ∂N_a/∂x_d at GP k
    @inbounds for a in 1:numNodes
        # divergence contribution of "a" basis into each component:
        # div( N_a * e_i ) = ∂N_a/∂x_i   (only i-th component contributes)
        # so for component i, the divergence is dN_a/dx_i
        @inbounds for b in 1:numNodes
            # assemble the pdim×pdim block for node pair (a,b):
            # Ke[(a,i),(b,j)] += (∂N_a/∂x_i) * (∂N_b/∂x_j) * w
            @inbounds for i in 1:dim
                dNa = ∂h[i, (k-1)*numNodes + a]
                ia  = (a - 1) * pdim + i
                @inbounds for j in 1:dim
                    dNb = ∂h[j, (k-1)*numNodes + b]
                    ib  = (b - 1) * pdim + j
                    Ke[ia, ib] += dNa * dNb * w
                end
            end
        end
    end
    return nothing
end

# Symmetric-gradient (strain) energy:
#   ∫ 2μ ε(u):ε(v) dΩ
# Uses engineering shear strains, so the "D" weights are:
#   3D: diag([2,2,2,1,1,1]) * μ
#   2D: diag([2,2,1]) * μ
#
# IMPORTANT:
#   - requires pdim == dim
#   - couples components (not block diagonal)
@inline function _kernel_symgrad!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @assert pdim == dim "_kernel_symgrad!: requires pdim == dim"

    @inbounds for a in 1:numNodes
        # gradients of shape Na at GP k
        ax = ∂h[1, (k-1)*numNodes + a]
        ay = (dim >= 2) ? ∂h[2, (k-1)*numNodes + a] : 0.0
        az = (dim == 3) ? ∂h[3, (k-1)*numNodes + a] : 0.0

        @inbounds for b in 1:numNodes
            bx = ∂h[1, (k-1)*numNodes + b]
            by = (dim >= 2) ? ∂h[2, (k-1)*numNodes + b] : 0.0
            bz = (dim == 3) ? ∂h[3, (k-1)*numNodes + b] : 0.0

            if dim == 1
                # 1D: 2μ (du/dx)(dv/dx)
                Ke[(a-1)*pdim + 1, (b-1)*pdim + 1] += (2.0 * ax * bx) * w

            elseif dim == 2
                ia1 = (a-1)*pdim + 1  # ux
                ia2 = (a-1)*pdim + 2  # uy
                ib1 = (b-1)*pdim + 1
                ib2 = (b-1)*pdim + 2

                # K11: 2*ax*bx + ay*by
                Ke[ia1, ib1] += (2.0*ax*bx + ay*by) * w
                # K22: 2*ay*by + ax*bx
                Ke[ia2, ib2] += (2.0*ay*by + ax*bx) * w
                # Couplings from gamma_xy = dux/dy + duy/dx
                Ke[ia1, ib2] += (ay*bx) * w
                Ke[ia2, ib1] += (ax*by) * w

            else
                ia1 = (a-1)*pdim + 1  # ux
                ia2 = (a-1)*pdim + 2  # uy
                ia3 = (a-1)*pdim + 3  # uz
                ib1 = (b-1)*pdim + 1
                ib2 = (b-1)*pdim + 2
                ib3 = (b-1)*pdim + 3

                # Diagonal blocks (normal + shear contributions)
                Ke[ia1, ib1] += (2.0*ax*bx + ay*by + az*bz) * w
                Ke[ia2, ib2] += (2.0*ay*by + ax*bx + az*bz) * w
                Ke[ia3, ib3] += (2.0*az*bz + ax*bx + ay*by) * w

                # Off-diagonal couplings from engineering shears:
                # gamma_xy = dux/dy + duy/dx
                Ke[ia1, ib2] += (ay*bx) * w
                Ke[ia2, ib1] += (ax*by) * w
                # gamma_xz = dux/dz + duz/dx
                Ke[ia1, ib3] += (az*bx) * w
                Ke[ia3, ib1] += (ax*bz) * w
                # gamma_yz = duy/dz + duz/dy
                Ke[ia2, ib3] += (az*by) * w
                Ke[ia3, ib2] += (ay*bz) * w
            end
        end
    end
    return nothing
end

# Curl-curl operator:
#   ∫ (∇×u) · (∇×v) dΩ
#
# NOTES:
# - Requires pdim == dim (vector field).
# - Uses standard Lagrange H¹ elements (NOT H(curl)-conforming).
# - Suitable for operator studies, stabilization, and educational purposes.
#
@inline function _kernel_curlcurl!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @assert pdim == dim "_kernel_curlcurl!: requires pdim == dim"

    @inbounds for a in 1:numNodes
        ax = ∂h[1, (k-1)*numNodes + a]
        ay = (dim >= 2) ? ∂h[2, (k-1)*numNodes + a] : 0.0
        az = (dim == 3) ? ∂h[3, (k-1)*numNodes + a] : 0.0

        @inbounds for b in 1:numNodes
            bx = ∂h[1, (k-1)*numNodes + b]
            by = (dim >= 2) ? ∂h[2, (k-1)*numNodes + b] : 0.0
            bz = (dim == 3) ? ∂h[3, (k-1)*numNodes + b] : 0.0

            if dim == 2
                # curl u = ∂uy/∂x - ∂ux/∂y   (scalar)
                ia1 = (a-1)*pdim + 1  # ux
                ia2 = (a-1)*pdim + 2  # uy
                ib1 = (b-1)*pdim + 1
                ib2 = (b-1)*pdim + 2

                # (∂uy/∂x - ∂ux/∂y)(∂vy/∂x - ∂vx/∂y)
                Ke[ia1, ib1] += ( ay * by ) * w
                Ke[ia2, ib2] += ( ax * bx ) * w
                Ke[ia1, ib2] += (-ay * bx) * w
                Ke[ia2, ib1] += (-ax * by) * w

            else
                # 3D curl:
                # cx = ∂uz/∂y - ∂uy/∂z
                # cy = ∂ux/∂z - ∂uz/∂x
                # cz = ∂uy/∂x - ∂ux/∂y
                ia1 = (a-1)*pdim + 1
                ia2 = (a-1)*pdim + 2
                ia3 = (a-1)*pdim + 3
                ib1 = (b-1)*pdim + 1
                ib2 = (b-1)*pdim + 2
                ib3 = (b-1)*pdim + 3

                # cx·cx
                Ke[ia2, ib2] += ( az * bz ) * w
                Ke[ia3, ib3] += ( ay * by ) * w
                Ke[ia2, ib3] += (-az * by) * w
                Ke[ia3, ib2] += (-ay * bz) * w

                # cy·cy
                Ke[ia1, ib1] += ( az * bz ) * w
                Ke[ia3, ib3] += ( ax * bx ) * w
                Ke[ia1, ib3] += (-az * bx) * w
                Ke[ia3, ib1] += (-ax * bz) * w

                # cz·cz
                Ke[ia1, ib1] += ( ay * by ) * w
                Ke[ia2, ib2] += ( ax * bx ) * w
                Ke[ia1, ib2] += (-ay * bx) * w
                Ke[ia2, ib1] += (-ax * by) * w
            end
        end
    end
    return nothing
end

@inline function _kernel_tensorlaplace_sym!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @assert pdim == 9 "tensorLaplaceMatrix requires pdim = 9 (TensorField)"

    # Standard tensor Laplace (componentwise)
    @inbounds for a in 1:numNodes, b in 1:numNodes
        s = 0.0
        @inbounds for d in 1:dim
            s += ∂h[d, (k-1)*numNodes + a] *
                 ∂h[d, (k-1)*numNodes + b]
        end
        val = s * w

        ia = (a-1)*pdim
        ib = (b-1)*pdim
        @inbounds for α in 1:9
            Ke[ia+α, ib+α] += val
        end
    end

    # Project operator to symmetric tensor subspace
    @inbounds for a in 1:numNodes, b in 1:numNodes
        ia = (a-1)*pdim
        ib = (b-1)*pdim

        # xy / yx
        s = 0.5 * (Ke[ia+4, ib+2] + Ke[ia+2, ib+4])
        Ke[ia+4, ib+2] = s
        Ke[ia+2, ib+4] = s

        # xz / zx
        s = 0.5 * (Ke[ia+7, ib+3] + Ke[ia+3, ib+7])
        Ke[ia+7, ib+3] = s
        Ke[ia+3, ib+7] = s

        # yz / zy
        s = 0.5 * (Ke[ia+8, ib+6] + Ke[ia+6, ib+8])
        Ke[ia+8, ib+6] = s
        Ke[ia+6, ib+8] = s
    end

    return nothing
end

@inline function _kernel_tracelaplace!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @assert pdim == 9 "traceLaplaceMatrix requires pdim = 9 (TensorField)"

    @inbounds for a in 1:numNodes, b in 1:numNodes
        s = 0.0
        @inbounds for d in 1:dim
            s += ∂h[d, (k-1)*numNodes + a] *
                 ∂h[d, (k-1)*numNodes + b]
        end
        val = s * w

        ia = (a-1)*pdim
        ib = (b-1)*pdim

        # trace indices: σxx, σyy, σzz
        @inbounds for α in (1,5,9), β in (1,5,9)
            Ke[ia+α, ib+β] += val
        end
    end

    return nothing
end

# Tensor div-div:
# ∫ (div σ) · (div τ) dΩ
#
# σ_ij = N_a * e_i ⊗ e_j
# (div σ)_i = ∑_j ∂σ_ij / ∂x_j = ∑_j ∂N_a/∂x_j
#
@inline function _kernel_tensordivdiv!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @assert pdim == dim^2 "_kernel_tensordivdiv!: requires pdim = dim^2"

    @inbounds for a in 1:numNodes
        @inbounds for b in 1:numNodes
            # tensor indices
            @inbounds for i in 1:dim          # divergence component
                @inbounds for j in 1:dim      # σ_ij
                    dNa = ∂h[j, (k-1)*numNodes + a]
                    ia  = (a - 1) * pdim + (i - 1) * dim + j

                    @inbounds for l in 1:dim  # τ_il
                        dNb = ∂h[l, (k-1)*numNodes + b]
                        ib  = (b - 1) * pdim + (i - 1) * dim + l

                        Ke[ia, ib] += dNa * dNb * w
                    end
                end
            end
        end
    end
    return nothing
end

# --- public API ---------------------------------------------------------------

"""
    poissonMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    ∫∇oN_c_∇oN_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global stiffness matrix for a Poisson-type diffusion operator

```

K_ab = ∫_Ω (∇N_a · ∇N_b) α(x) dΩ

```

`coefficient` can be a constant (`Number`) or an elementwise `ScalarField`
(interpolated to Gauss points using the Lagrange basis, as in the original implementation).
"""
function poissonMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    # same spirit as your preallocation: compute an upper estimate from all elements
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    # NOTE: exact sizing is tricky without re-walking physical groups; we keep a safe estimate
    #       and then resize down at the end (same as your pattern).
    estElems = sum(length.(elemTags))
    # worst-case: each element contributes (pdim*numNodes)^2 entries; approximate with (pdim*maxNodes)^2
    maxNodes = 0
    for et in elemTypes
        _, _, _, numNodes::Int64, _, _ = gmsh.model.mesh.getElementProperties(et)
        maxNodes = max(maxNodes, numNodes)
    end
    pdim = problem.pdim
    #lengthOfIJV = estElems * (pdim * maxNodes)^2
    lengthOfIJV = estimateLengthOfIJV(problem)

    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        pos = _assemble_poissonlike!(I, J, V, pos, problem, phName, coefficient, _kernel_stiffness!)
    end

    resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)

    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end
    
∫∇oN_c_∇oN_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0) = poissonMatrix(problem, coefficient=coefficient)

"""
    advectionMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0, dir::Int=1)
    ∫N_c_∂N∂x_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    ∫N_c_∂N∂y_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    ∫N_c_∂N∂z_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global convection (advection/drift) matrix

```

C_ab = ∫_Ω N_a (∂N_b/∂x_dir) β(x) dΩ

```

- `dir=1` corresponds to x, `dir=2` to y, `dir=3` to z (when available).
- `coefficient` can be a constant (`Number`) or an elementwise `ScalarField`,
  interpolated to Gauss points as in the stiffness implementation.
"""
function advectionMatrix(
    problem::Problem;
    coefficient::Union{Number,ScalarField}=1.0,
    dir::Int = 1)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    @assert 1 ≤ dir ≤ problem.dim "advectionMatrix: dir out of bounds"

    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    estElems = sum(length.(elemTags))
    maxNodes = 0
    for et in elemTypes
        _, _, _, numNodes::Int64, _, _ = gmsh.model.mesh.getElementProperties(et)
        maxNodes = max(maxNodes, numNodes)
    end
    pdim = problem.pdim
    #lengthOfIJV = estElems * (pdim * maxNodes)^2
    lengthOfIJV = estimateLengthOfIJV(problem)

    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        pos = _assemble_poissonlike!(I, J, V, pos, problem, phName, coefficient, _kernel_convection!; dir=dir)
    end

    resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)

    dof = problem.pdim * problem.non
    C = sparse(I, J, V, dof, dof)
    dropzeros!(C)
    return SystemMatrix(C, problem)
end

∫N_c_∂N∂x_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0) = advectionMatrix(problem, coefficient=coefficient, dir=1)
∫N_c_∂N∂y_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0) = advectionMatrix(problem, coefficient=coefficient, dir=2)
∫N_c_∂N∂z_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0) = advectionMatrix(problem, coefficient=coefficient, dir=3)

"""
    reactionMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    ∫N_c_N_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global (weighted) mass / reaction matrix

```

M_ab = ∫_Ω N_a N_b c(x) dΩ

```

`coefficient` can be a constant (`Number`) or an elementwise `ScalarField`,
interpolated to Gauss points using the Lagrange basis (same mechanism as stiffness).
"""
function reactionMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    estElems = sum(length.(elemTags))
    maxNodes = 0
    for et in elemTypes
        _, _, _, numNodes::Int64, _, _ = gmsh.model.mesh.getElementProperties(et)
        maxNodes = max(maxNodes, numNodes)
    end
    pdim = problem.pdim
    lengthOfIJV = estElems * (pdim * maxNodes)^2

    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        pos = _assemble_poissonlike!(I, J, V, pos, problem, phName, coefficient, _kernel_mass!)
    end

    resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)

    dof = problem.pdim * problem.non
    M = sparse(I, J, V, dof, dof)
    dropzeros!(M)
    return SystemMatrix(M, problem)
end

∫N_c_N_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0) = reactionMatrix(problem, coefficient=coefficient)

"""
    gradDivMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    ∫∇N_c_∇N_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global matrix corresponding to the grad-div bilinear form

```

G_ab = ∫_Ω (∇·u_h) (∇·v_h) α(x) dΩ

```

This is the weak form associated with the operator `∇(∇·u)` (up to sign conventions),
commonly appearing in linear elasticity and in grad-div stabilization.

Notes
- Requires `problem.pdim == problem.dim` (vector unknown with one component per spatial dimension).
- `coefficient` can be a constant (`Number`) or an elementwise `ScalarField`, interpolated to Gauss points
  using the Lagrange basis (same mechanism as in `poissonMatrix`).
"""
function gradDivMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    @assert problem.pdim == problem.dim "gradDivMatrix: requires problem.pdim == problem.dim"

    lengthOfIJV = estimateLengthOfIJV(problem)

    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        pos = _assemble_poissonlike!(I, J, V, pos, problem, phName, coefficient, _kernel_graddiv!)
    end

    resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)

    dof = problem.pdim * problem.non
    G = sparse(I, J, V, dof, dof)
    dropzeros!(G)
    return SystemMatrix(G, problem)
end

∫∇N_c_∇N_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0) = gradDivMatrix(problem, coefficient=coefficient)

"""
    symmetricGradientMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global matrix for the symmetric-gradient (strain) bilinear form

```

A(u,v) = ∫_Ω 2μ ε(u) : ε(v) dΩ
ε(u)   = 1/2 (∇u + ∇uᵀ)

```

Notes
- Requires `problem.pdim == problem.dim`.
- `coefficient` is typically the shear modulus `μ` (constant or `ScalarField`).
"""
function symmetricGradientMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    @assert problem.pdim == problem.dim "symmetricGradientMatrix: requires problem.pdim == problem.dim"

    lengthOfIJV = estimateLengthOfIJV(problem)

    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        pos = _assemble_poissonlike!(I, J, V, pos, problem, phName, coefficient, _kernel_symgrad!)
    end

    resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)

    dof = problem.pdim * problem.non
    A = sparse(I, J, V, dof, dof)
    dropzeros!(A)
    return SystemMatrix(A, problem)
end

"""
    curlCurlMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    ∫∇xN_c_∇xN_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global curl–curl matrix

```

C(u,v) = ∫_Ω (∇×u) · (∇×v) α(x) dΩ

```

Notes
- Requires `problem.pdim == problem.dim` (vector field).
- Uses standard Lagrange H¹ elements (NOT H(curl)-conforming).
- Intended for operator studies, stabilization terms, and educational use.
- Not suitable as a primary operator for Maxwell-type problems.
"""
function curlCurlMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    @assert problem.pdim == problem.dim "curlCurlMatrix: requires problem.pdim == problem.dim"

    lengthOfIJV = estimateLengthOfIJV(problem)

    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        pos = _assemble_poissonlike!(I, J, V, pos, problem, phName,
                                     coefficient, _kernel_curlcurl!)
    end

    resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)

    dof = problem.pdim * problem.non
    C = sparse(I, J, V, dof, dof)
    dropzeros!(C)
    return SystemMatrix(C, problem)
end
    
∫∇xN_c_∇xN_dΩ(problem::Problem; coefficient::Union{Number,ScalarField}=1.0) = curlCurlMatrix(problem, coefficient=coefficient)

"""
    tensorLaplaceMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the tensor Laplace operator

    ∫_Ω ∇(sym σ) : ∇(sym τ) dΩ

for a nodal TensorField with `pdim = 9`.

Notes
- Operator is restricted to the symmetric tensor subspace.
- `coefficient` may be constant or an elementwise ScalarField.
"""
function tensorLaplaceMatrix(
    problem::Problem;
    coefficient::Union{Number,ScalarField} = 1.0
)
    if problem.type == :dummy
        return nothing
    end
    @assert problem.pdim == 9 "tensorLaplaceMatrix requires pdim = 9"

    lengthOfIJV = estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        ph = problem.material[ipg].phName
        pos = _assemble_poissonlike!(
            I, J, V, pos,
            problem, ph,
            coefficient,
            _kernel_tensorlaplace_sym!
        )
    end

    resize!(I, pos-1)
    resize!(J, pos-1)
    resize!(V, pos-1)

    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)

    return SystemMatrix(K, problem)
end

"""
    traceLaplaceMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the trace Laplace operator

    ∫_Ω ∇tr(σ) · ∇tr(τ) dΩ

for a nodal TensorField with `pdim = 9`.
"""
function traceLaplaceMatrix(
    problem::Problem;
    coefficient::Union{Number,ScalarField} = 1.0
)
    if problem.type == :dummy
        return nothing
    end
    @assert problem.pdim == 9 "traceLaplaceMatrix requires pdim = 9"

    lengthOfIJV = estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        ph = problem.material[ipg].phName
        pos = _assemble_poissonlike!(
            I, J, V, pos,
            problem, ph,
            coefficient,
            _kernel_tracelaplace!
        )
    end

    resize!(I, pos-1)
    resize!(J, pos-1)
    resize!(V, pos-1)

    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)

    return SystemMatrix(K, problem)
end

function beltramiMichellMatrix(
    problem::Problem;
    coeff_laplace::Union{Number,ScalarField} = 1.0,
    coeff_trace::Union{Number,ScalarField} = 1.0
)
    K1 = tensorLaplaceMatrix(problem; coefficient = coeff_laplace)
    K2 = traceLaplaceMatrix(problem; coefficient = coeff_trace)
    return K1 + K2
end

"""
    tensorDivDivMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the tensor div–div matrix

```

D(σ,τ) = ∫_Ω (∇·σ_h) · (∇·τ_h) α(x) dΩ

```

Notes
- Requires `problem.pdim == dim^2` (second-order tensor field).
- Acts as an equilibrium-enforcing operator in stress-based formulations
  (e.g. Beltrami–Michell).
- `coefficient` may be a constant (`Number`) or an elementwise `ScalarField`.
"""
function tensorDivDivMatrix(
    problem::Problem;
    coefficient::Union{Number,ScalarField} = 1.0
)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    @assert problem.pdim == problem.dim^2 "tensorDivDivMatrix: requires pdim = dim^2"

    lengthOfIJV = estimateLengthOfIJV(problem)

    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        pos = _assemble_poissonlike!(
            I, J, V, pos,
            problem, phName,
            coefficient,
            _kernel_tensordivdiv!
        )
    end

    resize!(I, pos-1)
    resize!(J, pos-1)
    resize!(V, pos-1)

    dof = problem.pdim * problem.non
    D = sparse(I, J, V, dof, dof)
    dropzeros!(D)

    return SystemMatrix(D, problem)
end




function estimateLengthOfIJV(problem::Problem)
    gmsh.model.setCurrent(problem.name)

    lengthOfIJV = 0
    pdim = problem.pdim

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for dimTag in dimTags
            edim, etag = dimTag
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for i in 1:length(elemTypes)
                numElems = length(elemTags[i])
                numNodes =
                    div(length(elemNodeTags[i]), numElems)

                nloc = pdim * numNodes
                lengthOfIJV += numElems * nloc^2
            end
        end
    end

    return lengthOfIJV
end

"""
    sourceVector(problem, sources)

Assembles the right-hand-side vector corresponding to a volumetric source term
in Poisson-type and other scalar or vector-valued PDEs expressed in weak form.

Mathematically, this function assembles

    f_a = ∫_Ω N_a f dΩ

where the source term `f` may be given as a constant, a spatial function,
or a `ScalarField`.

This function is an alias of `loadVector` and shares the same implementation.
The difference is purely interpretational: `sourceVector` emphasizes the role
of the right-hand side as a PDE source term rather than a mechanical load.

Axisymmetric or spherically symmetric problems can be handled by including
the appropriate geometric weighting (e.g. `2πr`, `4πr²`) as a coefficient field.

In LowLevelFEM, right-hand-side vectors are assembled independently of the governing equation.
The same numerical machinery can therefore represent mechanical loads, heat sources,
or generic PDE source terms.

Returns a `ScalarField` or `VectorField`, depending on the problem field dimension.
"""
sourceVector = loadVector
∫N_c_dΩ = sourceVector

"""
    loadTensor(problem::Problem;
               source::Union{Matrix{Float64},TensorField} = zeros(3,3))

Assembles an L2 right-hand-side vector for tensor-valued problems (pdim = 9):

    f_{a,α} = ∫_Ω N_a(x) * S_α(x) dΩ

where `S` is either
- a constant 3×3 tensor (`Matrix{Float64}`), or
- a nodal `TensorField`.

Notes
- Intended for Beltrami–Michell and other stress-based formulations.
- No traction, pressure, or surface force interpretation.
"""
function loadTensor(
    problem::Problem;
    source::Union{Matrix{Float64},TensorField} = zeros(3,3)
    )
    @assert problem.pdim == 9 "loadTensor requires pdim = 9 (TensorField)"

    gmsh.model.setCurrent(problem.name)

    pdim = 9
    non  = problem.non
    fp   = zeros(pdim * non)

    # --- prepare source -------------------------------------------------------
    const_source = source isa Matrix
    field_source = source isa TensorField

    if const_source
        @assert size(source) == (3,3)
        # flatten column-wise: (xx,yx,zx,xy,yy,zy,xz,yz,zz)
        Sconst = vec(source)
        if all(iszero, Sconst)
            # trivial RHS
            return TensorField([], reshape(fp, :, 1), [0.0], [], 1, :tensor, problem)
        end
    else
        src_elem = nodesToElements(source)   # elemTag => 9-vector(s)
    end

    # --- integration ----------------------------------------------------------
    # We integrate:
    #   ∫ N_a * S dΩ
    # Only where S is defined / nonzero

    # Decide which elements to loop on
    elem_list = const_source ? gmsh.model.mesh.getElements(problem.dim, -1)[2] :
                               collect(keys(src_elem.A))

    for elem in elem_list
        # element type and properties
        etype, _, _, _ = gmsh.model.mesh.getElement(elem)
        _, _, order, numNodes, _, _ =
            gmsh.model.mesh.getElementProperties(etype)

        intPoints, intWeights =
            gmsh.model.mesh.getIntegrationPoints(
                etype, "Gauss" * string(2order + 1)
            )

        comp, fun, ori =
            gmsh.model.mesh.getBasisFunctions(etype, intPoints, "Lagrange")
        h = reshape(fun, :, length(intWeights))

        nodeTags = gmsh.model.mesh.getElement(elem)[2]
        jac, jacDet, coord =
            gmsh.model.mesh.getJacobian(elem, intPoints)

        for k in 1:length(intWeights)
            w = jacDet[k] * intWeights[k]

            # tensor source at GP
            if const_source
                S = Sconst
            else
                # nodal tensor projected with shape functions
                S = zeros(pdim)
                Se = src_elem.A[elem][:,1]
                for a in 1:numNodes
                    S .+= h[a,k] * Se
                end
            end

            for a in 1:numNodes
                Na = h[a,k]
                base = (nodeTags[a]-1)*pdim
                @inbounds for α in 1:pdim
                    fp[base+α] += Na * S[α] * w
                end
            end
        end
    end

    return TensorField([], reshape(fp, :, 1), [0.0], [], 1, :tensor, problem)
end

"""
    solveField(K, f; support=BoundaryCondition[], iterative=false,
               reltol=sqrt(eps()), maxiter=..., preconditioner=Identity(),
               ordering=true)

Solves a linear static field problem with Dirichlet boundary conditions.

The algebraic system
```

K * u = f

```
is solved for the unknown field `u`, where prescribed values are enforced
*solver-side* via boundary conditions.

The solution is obtained by partitioning the degrees of freedom into
free and constrained sets. On constrained DOFs the values prescribed
by `support` override the solution, while on free DOFs the reduced system
is solved:
```

K_ff * u_f = f_f − K_fc * u_c

```

The function supports both direct and iterative linear solvers and works
uniformly for `ScalarField` and `VectorField` unknowns.

---

### Arguments
- `K::SystemMatrix`  
  Global system matrix.
- `f::Union{ScalarField,VectorField}`  
  Right-hand-side vector.
- `support::Vector{BoundaryCondition}` (keyword, default = `BoundaryCondition[]`)  
  Dirichlet-type boundary conditions.
- `iterative::Bool` (keyword, default = `false`)  
  If `true`, use an iterative solver (conjugate gradient).
- `reltol::Real` (keyword, default = `sqrt(eps())`)  
  Relative tolerance for the iterative solver.
- `maxiter::Int` (keyword, default = `K.model.non * K.model.dim`)  
  Maximum number of iterations for the iterative solver.
- `preconditioner` (keyword, default = `Identity()`)  
  Preconditioner for the iterative solver.
- `ordering::Bool` (keyword, default = `true`)  
  If `false`, use an explicit LU factorization without reordering.

---

### Returns
- `u::Union{ScalarField,VectorField}`  
  Solution field with prescribed values enforced on constrained DOFs.

---

### Notes
- Boundary conditions are applied *inside* the solver; the input field `f`
  is not modified.
- The algorithm itself is agnostic to the physical meaning of the field
  (scalar, vector, tensor), as long as `K` and `f` are dimensionally consistent.
- When `iterative = true`, the system is solved using conjugate gradient on
  the reduced matrix `K_ff`.
"""
function solveField(K::SystemMatrix, f::Union{ScalarField,VectorField}; 
                           support::Vector{BoundaryCondition}=BoundaryCondition[], 
                           iterative=false,
                           reltol::Real = sqrt(eps()),
                           maxiter::Int = K.model.non * K.model.dim,
                           preconditioner = Identity(),
                           ordering=true)
    problem = K.model

    fixed = constrainedDoFs(problem, support)
    free = freeDoFs(problem, support)
    u = copy(f)
    fill!(u.a, 0.0)
    applyBoundaryConditions!(u, support)
    f_kin = K.A[:, fixed] * u.a[fixed]
    #u.a[free] = cholesky(Symmetric(K.A[free, free])) \ (f.a[free] - f_kin[free])
    if iterative
        u.a[free] = cg(K.A[free,free], f.a[free] - f_kin[free], Pl=preconditioner, reltol=reltol, maxiter=maxiter)
    elseif ordering == false
        u.a[free] = lu(K.A[free, free], q=nothing) \ (f.a[free] - f_kin[free])
    else
        u.a[free] = (K.A[free, free]) \ (f.a[free] - f_kin[free])
    end
    return u
end

"""
    gradMatrix(problem_u::Problem, problem_p::Problem)

Assembles the mixed gradient matrix G mapping a scalar field (pressure)
to a vector field (velocity), corresponding to the weak form

    (G p, v) = -∫_Ω p ∇·v dΩ

Notes
- Requires problem_u.pdim == problem_u.dim (vector field).
- Requires problem_p.pdim == 1 (scalar field).
- Requires problem_u.dim == problem_p.dim.
- The divergence matrix is obtained as D = -G'.
"""
function gradMatrix(problem_u::Problem, problem_p::Problem)

    @assert problem_u.dim == problem_p.dim
    @assert problem_u.pdim == problem_u.dim
    @assert problem_p.pdim == 1

    gmsh.model.setCurrent(problem_u.name)

    dim   = problem_u.dim
    pdim  = problem_u.pdim
    non_u = problem_u.non
    non_p = problem_p.non

    # Conservative IJV estimate
    lengthOfIJV = estimateLengthOfIJV(problem_u)

    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem_u.material)
        phName = problem_u.material[ipg].phName
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for (edim, etag) in dimTags
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for itype in 1:length(elemTypes)
                et = elemTypes[itype]

                _, _, order, numNodes::Int, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                intPoints, intWeights =
                    gmsh.model.mesh.getIntegrationPoints(
                        et, "Gauss" * string(2order + 1)
                    )

                numIntPoints = length(intWeights)

                # --- pressure basis (scalar) ---
                _, fun_p, _ =
                    gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h_p = reshape(fun_p, :, numIntPoints)

                # --- velocity basis (vector) ---
                _, dfun_u, _ =
                    gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h_u = reshape(dfun_u, :, numIntPoints)

                nnet = zeros(Int, length(elemTags[itype]), numNodes)

                invJac = zeros(3, 3numIntPoints)
                ∂h_u  = zeros(dim, numNodes * numIntPoints)

                Ke = zeros(pdim * numNodes, numNodes)

                for j in 1:length(elemTags[itype])
                    elem = elemTags[itype][j]

                    for a in 1:numNodes
                        nnet[j, a] =
                            elemNodeTags[itype][(j-1)*numNodes + a]
                    end

                    jac, jacDet, _ =
                        gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)

                    for k in 1:numIntPoints
                        invJac[1:3, 3k-2:3k] .=
                            inv(Jac[1:3, 3k-2:3k])'
                    end

                    fill!(∂h_u, 0.0)
                    for k in 1:numIntPoints, a in 1:numNodes
                        ∂h_u[1:dim, (k-1)*numNodes + a] .=
                            invJac[1:dim, 3k-2:3k-(3-dim)] *
                            ∇h_u[a*3-2:a*3-(3-dim), k]
                    end

                    fill!(Ke, 0.0)

                    for k in 1:numIntPoints
                        w = jacDet[k] * intWeights[k]
                        _kernel_grad_scalar_to_vector!(
                            Ke, w, k,
                            h_p, ∂h_u,
                            numNodes, numNodes,
                            pdim, dim
                        )
                    end

                    # scatter to global IJV
                    for a in 1:(pdim*numNodes)
                        na = div(a-1, pdim) + 1
                        Ia = (nnet[j, na]-1)*pdim + mod(a-1, pdim) + 1

                        for b in 1:numNodes
                            Ib = nnet[j, b]
                            I[pos] = Ia
                            J[pos] = Ib
                            V[pos] = Ke[a, b]
                            pos += 1
                        end
                    end
                end
            end
        end
    end

    resize!(I, pos-1)
    resize!(J, pos-1)
    resize!(V, pos-1)

    K = sparse(I, J, V, pdim*non_u, non_p)
    dropzeros!(K)

    return SystemMatrix(K, problem_p, problem_u)
end

@inline function _kernel_grad_scalar_to_vector!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h_p, ∂h_u,
    numNodes_p::Int,
    numNodes_u::Int,
    pdim_u::Int,
    dim::Int
)
    @inbounds for a in 1:numNodes_u
        @inbounds for b in 1:numNodes_p
            Nb = h_p[(k-1)*numNodes_p + b]
            @inbounds for α in 1:dim
                dNa = ∂h_u[α, (k-1)*numNodes_u + a]
                ia  = (a-1)*pdim_u + α
                ib  = b
                Ke[ia, ib] -= Nb * dNa * w
            end
        end
    end
end

# --- Navier–Stokes / Oseen advection (vector -> vector) -----------------------

"""
    navierStokesAdvectionMatrix(problem::Problem, u;
                                time_index::Int = 1)

Assembles the (linearized) Navier–Stokes convection operator in Oseen/Picard form:

    C(w,v) = ∫_Ω (u·∇)w · v dΩ

where `u` is a given velocity field (the coefficient), `w` is the unknown velocity,
and `v` is the test function.

Matrix size:
- `problem.pdim == problem.dim` is required (vector field).
- Returns a `SystemMatrix` of size (dim*N) × (dim*N).

Arguments
- `problem::Problem`:
  Must represent the velocity unknown (vector field).
- `u`:
  Given advection velocity in the SAME dof ordering as `VectorField` uses:
  `[u1x,u1y(,u1z), u2x,u2y(,u2z), ...]`.
  Accepted forms:
  - `VectorField` (nodal, same mesh)
  - `AbstractVector{<:Real}` (length = dim*problem.non)

Keyword arguments
- `time_index`:
  If `u` stores multiple time steps internally (e.g. as a matrix), you can select
  which column to use. For plain vectors it is ignored.

Notes
- This is NOT the scalar `advectionMatrix`. This operator uses the given velocity `u`
  inside the Gauss integration, assembled in ONE pass (no splitting to 3 scalar fields).
- The resulting operator is block-diagonal by velocity components (no component mixing),
  consistent with (u·∇) acting component-wise in Cartesian coordinates.
"""
function navierStokesAdvectionMatrix(
    problem::Problem,
    u;
    time_index::Int = 1
)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    @assert problem.pdim == problem.dim "navierStokesAdvectionMatrix: requires vector field (pdim == dim)."

    dim  = problem.dim
    pdim = problem.pdim
    N    = problem.non
    dof  = pdim * N

    # --- extract dof vector for u ------------------------------------------------
    uvec = _ns_get_uvec(u, dof, time_index)

    # prealloc (same strategy as others)
    lengthOfIJV = estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    # loop physical groups (same as poissonMatrix/advectionMatrix)
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        pos = _assemble_navier_stokes_advection!(I, J, V, pos, problem, phName, uvec)
    end

    resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)

    C = sparse(I, J, V, dof, dof)
    dropzeros!(C)
    return SystemMatrix(C, problem)
end

# --- internal: get u dofs as a vector ----------------------------------------

# Try to be robust without touching existing code:
# - If u is already a vector: use it.
# - If u is a VectorField: try DoFs(u), else u.a (vector), else u.a[:,time_index].
function _ns_get_uvec(u, dof::Int, time_index::Int)
    if u isa AbstractVector
        @assert length(u) == dof "navierStokesAdvectionMatrix: length(u) must be dim*non."
        return u
    end

    # VectorField path: prefer DoFs(u) if available
    try
        uv = DoFs(u)
        if uv isa AbstractVector
            @assert length(uv) == dof "navierStokesAdvectionMatrix: DoFs(u) has wrong length."
            return uv
        end
        # If DoFs(u) returns a matrix (time-dependent storage), take a column
        if uv isa AbstractMatrix
            @assert size(uv, 1) == dof "navierStokesAdvectionMatrix: DoFs(u) has wrong size."
            return view(uv, :, time_index)
        end
    catch
        # fallthrough
    end

    # fallback to field storage convention u.a
    if hasproperty(u, :a)
        A = getproperty(u, :a)
        if A isa AbstractVector
            @assert length(A) == dof "navierStokesAdvectionMatrix: u.a has wrong length."
            return A
        elseif A isa AbstractMatrix
            @assert size(A, 1) == dof "navierStokesAdvectionMatrix: u.a has wrong size."
            return view(A, :, time_index)
        end
    end

    error("navierStokesAdvectionMatrix: cannot extract velocity dofs from `u`. Provide a VectorField or a dof vector.")
end

# --- kernel -------------------------------------------------------------------

# Convection operator: ∫ Na * (u·∇Nb) * w dΩ, block-diagonal by component.
@inline function _kernel_ns_convection!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int,
    u_gp::AbstractVector{<:Real}
)
    @inbounds for a in 1:numNodes
        Na = h[(k-1)*numNodes + a]
        @inbounds for b in 1:numNodes
            adv = 0.0
            @inbounds for d in 1:dim
                adv += u_gp[d] * ∂h[d, (k-1)*numNodes + b]
            end
            _add_blockdiag!(Ke, pdim, a, b, Na * adv * w)
        end
    end
    return nothing
end

# --- assembly for one physical group -----------------------------------------

function _assemble_navier_stokes_advection!(
    I::Vector{Int}, J::Vector{Int}, V::Vector{Float64},
    pos0::Int,
    problem::Problem,
    phName::String,
    uvec::AbstractVector{<:Real}
)
    gmsh.model.setCurrent(problem.name)

    pdim = problem.pdim
    dim  = problem.dim
    pos  = pos0

    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

    for dimTag in dimTags
        edim, etag = dimTag
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

        for itype in 1:length(elemTypes)
            et = elemTypes[itype]
            _, _, order, numNodes::Int64, _, _ = gmsh.model.mesh.getElementProperties(et)

            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
            numIntPoints = length(intWeights)

            # basis and grad basis
            _, fun, _  = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h = reshape(fun, :, numIntPoints)

            _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
            ∇h = reshape(dfun, :, numIntPoints)   # (3*numNodes, numIntPoints)

            # connectivity
            nnet = zeros(Int, length(elemTags[itype]), numNodes)
            @inbounds for j in 1:length(elemTags[itype])
                for a in 1:numNodes
                    nnet[j, a] = elemNodeTags[itype][(j-1)*numNodes + a]
                end
            end

            invJac = zeros(3, 3numIntPoints)
            ∂h     = zeros(dim, numNodes * numIntPoints)
            Ke     = zeros(pdim * numNodes, pdim * numNodes)

            # local velocity nodal values: u_loc[d,a]
            u_loc  = zeros(Float64, dim, numNodes)
            u_gp   = zeros(Float64, dim)

            @inbounds for j in 1:length(elemTags[itype])
                elem = elemTags[itype][j]

                jac, jacDet, _ = gmsh.model.mesh.getJacobian(elem, intPoints)
                Jac = reshape(jac, 3, :)

                @inbounds for k in 1:numIntPoints
                    invJac[1:3, 3*k-2:3*k] .= inv(Jac[1:3, 3*k-2:3*k])'
                end

                # physical gradients of basis
                fill!(∂h, 0.0)
                @inbounds for k in 1:numIntPoints, a in 1:numNodes
                    invJk = invJac[1:dim, 3*k-2:3*k-(3-dim)]
                    gha  = ∇h[a*3-2 : a*3-(3-dim), k]
                    ∂h[1:dim, (k-1)*numNodes + a] .= invJk * gha
                end

                # gather u nodal values for this element
                @inbounds for a in 1:numNodes
                    node = nnet[j, a]
                    base = (node - 1) * pdim
                    for d in 1:dim
                        u_loc[d, a] = uvec[base + d]
                    end
                end

                fill!(Ke, 0.0)

                # integrate
                @inbounds for k in 1:numIntPoints
                    # interpolate u to Gauss point: u_gp[d] = Σ_a N_a * u_loc[d,a]
                    for d in 1:dim
                        u_gp[d] = dot(view(u_loc, d, :), view(h, :, k))
                    end

                    w = jacDet[k] * intWeights[k]
                    _kernel_ns_convection!(Ke, w, k, h, ∂h, numNodes, pdim, dim, u_gp)
                end

                # scatter to global I,J,V
                @inbounds for a in 1:(pdim*numNodes)
                    na = (div(a-1, pdim) + 1)
                    Ia_node = nnet[j, na]
                    Ia = (Ia_node - 1) * pdim + (mod(a-1, pdim) + 1)

                    @inbounds for b in 1:(pdim*numNodes)
                        nb = (div(b-1, pdim) + 1)
                        Jb_node = nnet[j, nb]
                        Jb = (Jb_node - 1) * pdim + (mod(b-1, pdim) + 1)

                        I[pos] = Ia
                        J[pos] = Jb
                        V[pos] = Ke[a, b]
                        pos += 1
                    end
                end
            end
        end
    end

    return pos
end

#=
@inline function _kernel_ns_convection!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int,
    u_gp::AbstractVector{<:Real}
)
    @inbounds for a in 1:numNodes
        Na = h[(k-1)*numNodes + a]
        @inbounds for b in 1:numNodes
            adv = 0.0
            @inbounds for d in 1:dim
                adv += u_gp[d] * ∂h[d, (k-1)*numNodes + b]
            end
            _add_blockdiag!(Ke, pdim, a, b, Na * adv * w)
        end
    end
    return nothing
end
=#

#=
function _assemble_ns_advection!(
    I::Vector{Int}, J::Vector{Int}, V::Vector{Float64},
    pos0::Int,
    problem::Problem,
    phName::String,
    uvec::AbstractVector{<:Real}
)
    gmsh.model.setCurrent(problem.name)

    pdim = problem.pdim
    dim  = problem.dim
    pos  = pos0

    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

    for (edim, etag) in dimTags
        elemTypes, elemTags, elemNodeTags =
            gmsh.model.mesh.getElements(edim, etag)

        for itype in 1:length(elemTypes)
            et = elemTypes[itype]
            _, _, order, numNodes::Int, _, _ =
                gmsh.model.mesh.getElementProperties(et)

            intPoints, intWeights =
                gmsh.model.mesh.getIntegrationPoints(
                    et, "Gauss" * string(2order + 1)
                )
            numIntPoints = length(intWeights)

            # basis
            _, fun, _ =
                gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h = reshape(fun, :, numIntPoints)

            _, dfun, _ =
                gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
            ∇h = reshape(dfun, :, numIntPoints)

            # connectivity
            nnet = zeros(Int, length(elemTags[itype]), numNodes)
            @inbounds for j in 1:length(elemTags[itype])
                for a in 1:numNodes
                    nnet[j, a] =
                        elemNodeTags[itype][(j-1)*numNodes + a]
                end
            end

            invJac = zeros(3, 3numIntPoints)
            ∂h     = zeros(dim, numNodes * numIntPoints)
            Ke     = zeros(pdim * numNodes, pdim * numNodes)

            u_loc = zeros(Float64, dim, numNodes)
            u_gp  = zeros(Float64, dim)

            @inbounds for j in 1:length(elemTags[itype])
                elem = elemTags[itype][j]

                jac, jacDet, _ =
                    gmsh.model.mesh.getJacobian(elem, intPoints)
                Jac = reshape(jac, 3, :)

                for k in 1:numIntPoints
                    invJac[1:3, 3k-2:3k] .=
                        inv(Jac[1:3, 3k-2:3k])'
                end

                # physical gradients
                fill!(∂h, 0.0)
                for k in 1:numIntPoints, a in 1:numNodes
                    ∂h[1:dim, (k-1)*numNodes + a] .=
                        invJac[1:dim, 3k-2:3k-(3-dim)] *
                        ∇h[a*3-2:a*3-(3-dim), k]
                end

                # gather nodal u
                for a in 1:numNodes
                    node = nnet[j, a]
                    base = (node-1)*pdim
                    for d in 1:dim
                        u_loc[d, a] = uvec[base + d]
                    end
                end

                fill!(Ke, 0.0)

                for k in 1:numIntPoints
                    for d in 1:dim
                        u_gp[d] = dot(view(u_loc, d, :), view(h, :, k))
                    end
                    w = jacDet[k] * intWeights[k]
                    _kernel_ns_convection!(
                        Ke, w, k, h, ∂h, numNodes, pdim, dim, u_gp
                    )
                end

                # scatter
                for a in 1:(pdim*numNodes)
                    na = div(a-1, pdim) + 1
                    Ia = (nnet[j, na]-1)*pdim + mod(a-1, pdim) + 1
                    for b in 1:(pdim*numNodes)
                        nb = div(b-1, pdim) + 1
                        Jb = (nnet[j, nb]-1)*pdim + mod(b-1, pdim) + 1
                        I[pos] = Ia
                        J[pos] = Jb
                        V[pos] = Ke[a, b]
                        pos += 1
                    end
                end
            end
        end
    end

    return pos
end
=#