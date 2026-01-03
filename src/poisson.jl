export stiffnessMatrixPoisson, convectionMatrixPoisson, massMatrixPoisson
export gradDivMatrix, symmetricGradientMatrix, curlCurlMatrix

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

# --- public API ---------------------------------------------------------------

"""
    stiffnessMatrixPoisson(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global stiffness matrix for a Poisson-type diffusion operator

```

K_ab = ∫_Ω (∇N_a · ∇N_b) α(x) dΩ

```

`coefficient` can be a constant (`Number`) or an elementwise `ScalarField`
(interpolated to Gauss points using the Lagrange basis, as in the original implementation).
"""
function stiffnessMatrixPoisson(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
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

"""
    convectionMatrixPoisson(problem::Problem; coefficient::Union{Number,ScalarField}=1.0, dir::Int=1)

Assembles the global convection (advection/drift) matrix

```

C_ab = ∫_Ω N_a (∂N_b/∂x_dir) β(x) dΩ

```

- `dir=1` corresponds to x, `dir=2` to y, `dir=3` to z (when available).
- `coefficient` can be a constant (`Number`) or an elementwise `ScalarField`,
  interpolated to Gauss points as in the stiffness implementation.
"""
function convectionMatrixPoisson(
    problem::Problem;
    coefficient::Union{Number,ScalarField}=1.0,
    dir::Int = 1)
    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    @assert 1 ≤ dir ≤ problem.dim "convectionMatrixPoisson: dir out of bounds"

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

"""
    massMatrixPoisson(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global (weighted) mass / reaction matrix

```

M_ab = ∫_Ω N_a N_b c(x) dΩ

```

`coefficient` can be a constant (`Number`) or an elementwise `ScalarField`,
interpolated to Gauss points using the Lagrange basis (same mechanism as stiffness).
"""
function massMatrixPoisson(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
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

"""
    gradDivMatrix(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

Assembles the global matrix corresponding to the grad-div bilinear form

```

G_ab = ∫_Ω (∇·u_h) (∇·v_h) α(x) dΩ

```

This is the weak form associated with the operator `∇(∇·u)` (up to sign conventions),
commonly appearing in linear elasticity and in grad-div stabilization.

Notes
- Requires `problem.pdim == problem.dim` (vector unknown with one component per spatial dimension).
- `coefficient` can be a constant (`Number`) or an elementwise `ScalarField`, interpolated to Gauss points
  using the Lagrange basis (same mechanism as in `stiffnessMatrixPoisson`).
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