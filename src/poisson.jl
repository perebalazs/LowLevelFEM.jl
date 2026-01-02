export stiffnessMatrixPoisson_old
export stiffnessMatrixPoisson, convectionMatrixPoisson, massMatrixPoisson, convectionMatrixReynoldsSkew
export couetteMatrixReynoldsSkew


"""
    stiffnessMatrixPoisson_old(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

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
"""
function stiffnessMatrixPoisson_old(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
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

function _assemble_reynolds_couette!(
    I::Vector{Int}, J::Vector{Int}, V::Vector{Float64},
    pos0::Int,
    problem::Problem,
    phName::String;
    velocity::Number,
    eta::Number = problem.material[1].η,
    dir::Int = 1,
    hfield::ScalarField = problem.geometry.h,        # nodal ScalarField
    dhdxfield::ScalarField = problem.geometry.dhdx   # nodal ScalarField (∂h/∂x_dir)
)
    gmsh.model.setCurrent(problem.name)

    pdim = problem.pdim
    dim  = problem.dim
    @assert 1 ≤ dir ≤ dim "Reynolds Couette: dir out of bounds"

    etaV = Float64(6.0 * eta * velocity / 2.0)   # matches your old code: 6*η*V/2

    pos = pos0
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

    for idm in 1:length(dimTags)
        edim, etag = dimTags[idm]
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

        for itype in 1:length(elemTypes)
            et = elemTypes[itype]
            _, _, order, numNodes::Int64, _, _ = gmsh.model.mesh.getElementProperties(et)

            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
            numIntPoints = length(intWeights)

            comp, fun, _  = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            hN = reshape(fun, :, numIntPoints)      # (numNodes, numIntPoints)

            comp, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
            ∇N = reshape(dfun, :, numIntPoints)     # (3*numNodes, numIntPoints)

            nnet   = zeros(Int, length(elemTags[itype]), numNodes)
            invJac = zeros(3, 3numIntPoints)
            ∂N     = zeros(dim, numNodes * numIntPoints)
            Ke     = zeros(pdim * numNodes, pdim * numNodes)

            # connectivity
            @inbounds for j in 1:length(elemTags[itype]), a in 1:numNodes
                nnet[j, a] = elemNodeTags[itype][(j-1)*numNodes + a]
            end

            @inbounds for j in 1:length(elemTags[itype])
                elem = elemTags[itype][j]

                # nodal h and dhdx for this element (from nodal fields)
                hn    = @view hfield.a[nnet[j, :]]
                dhdxn = @view dhdxfield.a[nnet[j, :]]

                jac, jacDet, _ = gmsh.model.mesh.getJacobian(elem, intPoints)
                Jac = reshape(jac, 3, :)

                # inv(J)' per GP
                @inbounds for k in 1:numIntPoints
                    invJac[1:3, 3*k-2:3*k] .= inv(Jac[1:3, 3*k-2:3*k])'
                end

                # physical grads of shape: ∂N
                fill!(∂N, 0.0)
                @inbounds for k in 1:numIntPoints, a in 1:numNodes
                    invJk = invJac[1:dim, 3*k-2:3*k-(3-dim)]
                    gha   = ∇N[a*3-2 : a*3-(3-dim), k]
                    ∂N[1:dim, (k-1)*numNodes + a] .= invJk * gha
                end

                fill!(Ke, 0.0)

                # integrate
                @inbounds for k in 1:numIntPoints
                    # interpolate h and dhdx to Gauss point using N
                    Nk = @view hN[:, k]
                    h_gp    = dot(hn, Nk)
                    dhdx_gp = dot(dhdxn, Nk)

                    w = jacDet[k] * intWeights[k]
                    _kernel_reynolds_couette!(Ke, w, k, hN, ∂N, h_gp, dhdx_gp, numNodes, pdim; dir = dir, etaV = etaV)
                end

                # scatter
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

# Convection (test-derivative):
# ∫ (∂Na/∂x_dir) * Nb * w dA
@inline function _kernel_convection_test!(
    Ke::Matrix{Float64},
    w::Float64, k::Int,
    h, ∂h,
    numNodes::Int, pdim::Int, dim::Int, elem::Integer;
    dir::Int = 1
)
    @inbounds for a in 1:numNodes
        dNa = ∂h[dir, (k-1)*numNodes + a]
        @inbounds for b in 1:numNodes
            Nb = h[(k-1)*numNodes + b]
            _add_blockdiag!(Ke, pdim, a, b, dNa * Nb * w)
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

# Reynolds–Couette kernel (skew/product-rule form):
# Ke_ab += (6ηV/2) * [ Na * (∂Nb/∂x_dir) * h  -  Na*Nb * (∂h/∂x_dir) ] * w
@inline function kernel_reynolds_couette!(
    Ke::Matrix{Float64},
    w::Float64,
    k::Int,
    hN,          # shape: (numNodes, numGP)
    dN,          # shape: (dim, numNodes*numGP)
    numNodes::Int,
    pdim::Int,
    dim::Int;
    h_gp::Float64,
    dhdx_gp::Float64,
    dir::Int = 1
)
    @inbounds for a in 1:numNodes
        Na = hN[(k-1)*numNodes + a]

        @inbounds for b in 1:numNodes
            Nb  = hN[(k-1)*numNodes + b]
            dNb = dN[dir, (k-1)*numNodes + b]

            # skew-adjoint Reynolds Couette term
            val = Na * dNb * h_gp - Na * Nb * dhdx_gp

            _add_blockdiag!(Ke, pdim, a, b, val * w)
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
    dir::Int = 1
)
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

function convectionMatrixReynoldsSkew(
    problem::Problem;
    coefficient::Union{Number,ScalarField}=1.0,
    dir::Int = 1
)
    gmsh.model.setCurrent(problem.name)
    @assert 1 ≤ dir ≤ problem.dim

    lengthOfIJV = estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)

    pos = 1
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName

        # trial-derivative part:  ∫ Na * ∂Nb
        pos1 = _assemble_poissonlike!(
            I, J, V, pos,
            problem, phName, coefficient,
            _kernel_convection!;
            dir = dir
        )

        # test-derivative part: ∫ ∂Na * Nb
        pos2 = _assemble_poissonlike!(
            I, J, V, pos,
            problem, phName, coefficient,
            _kernel_convection_test!;
            dir = dir
        )

        # antisymmetrize locally in V
        @inbounds for k in pos:pos2-1
            V[k] *= 0.5
        end
        @inbounds for k in pos1:pos2-1
            V[k] *= -0.5
        end

        pos = pos2
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
    couetteMatrixReynoldsSkew(problem; velocity, eta=problem.material[1].η, dir=1,
                              hfield=problem.geometry.h, dhdxfield=problem.geometry.dhdx)

Assembles the Reynolds–Couette (wedge) operator in skew/product-rule form:

K2_ab = (6ηV/2) ∫ [ N_a * h * ∂_{x_dir}N_b  -  N_a*N_b*∂_{x_dir}h ] dΩ

This reproduces the stable behavior of the original `systemMatrix` implementation.
"""
function couetteMatrixReynoldsSkew_old(
    problem::Problem;
    velocity::Number,
    dir::Int = 1
)
    gmsh.model.setCurrent(problem.name)

    h     = problem.geometry.h
    dhdx = problem.geometry.dhdx
    η     = problem.material[1].η

    h_e     = nodesToElements(h)
    dhdx_e  = nodesToElements(dhdx)

    _, pa_h     = _build_elemwise_coeff_dict(h_e)
    _, pa_dhdx  = _build_elemwise_coeff_dict(dhdx_e)

    coeff = 6.0 * η * velocity

    lengthOfIJV = estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName

        pos = _assemble_poissonlike!(
            I, J, V, pos,
            problem,
            phName,
            coeff,
            (Ke, w, k, hN, dN, numNodes, pdim, dim; dir=dir) -> begin
                elem = _assemble_poissonlike_current_element()  # lásd lent
                h_gp    = dot(pa_h[elem][:,1], view(hN, :, k))
                dhdx_gp = dot(pa_dhdx[elem][:,1], view(hN, :, k))

                kernel_reynolds_couette!(
                    Ke, w, k,
                    hN, dN,
                    numNodes, pdim, dim;
                    h_gp = h_gp,
                    dhdx_gp = dhdx_gp,
                    dir = dir
                )
            end;
            dir = dir
        )
    end

    resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)

    dof = problem.pdim * problem.non
    C = sparse(I, J, V, dof, dof)
    dropzeros!(C)

    return SystemMatrix(C, problem)
end

#=
### Használat Reynolds első 3 tagra (a te jelöléseiddel)

* stiffness: `coefficient = α*k*h^3`
* convection: `coefficient = 6ηz*V*h`, `dir=1` (x-irány)
* mass/reaction: `coefficient = 6ηz*∂x(h)`

```julia
Kdiff = stiffnessMatrixPoisson(problem; coefficient = α*k*h3)
Kconv = convectionMatrixPoisson(problem; coefficient = 6ηz*V*h, dir=1)
Kmass = massMatrixPoisson(problem; coefficient = 6ηz*∂x(h))
K = Kdiff + Kconv + Kmass
```
=#

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

function couetteMatrixReynoldsSkew(
    problem::Problem;
    velocity::Number,
    dir::Int = 1
)
    gmsh.model.setCurrent(problem.name)

    h     = problem.geometry.h
    dhdx = problem.geometry.dhdx
    η     = problem.material[1].η

    h_e    = nodesToElements(h)
    dhdx_e = nodesToElements(dhdx)

    pa_h    = Dict(zip(h_e.numElem, h_e.A))
    pa_dhdx = Dict(zip(dhdx_e.numElem, dhdx_e.A))

    pdim = problem.pdim
    dim  = problem.dim
    coef = 6.0 * η * velocity / 2.0

    lengthOfIJV = estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)

    pos = 1

    for ipg in eachindex(problem.material)
        phName = problem.material[ipg].phName
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for (edim, etag) in dimTags
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for itype in eachindex(elemTypes)
                et = elemTypes[itype]
                _, _, order, numNodes::Int, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                intPoints, intWeights =
                    gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numGP = length(intWeights)

                comp, fun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                N = reshape(fun, :, numGP)

                comp, dfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇N = reshape(dfun, :, numGP)

                for (eidx, elem) in enumerate(elemTags[itype])
                    nodes = elemNodeTags[itype][(eidx-1)*numNodes+1 : eidx*numNodes]
                    Ke = zeros(pdim*numNodes, pdim*numNodes)

                    jac, jacDet, _ =
                        gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)

                    for k in 1:numGP
                        invJ = inv(Jac[1:dim, 3k-2:3k-(3-dim)])
                        dN = zeros(dim, numNodes)
                        for a in 1:numNodes
                            dN[:,a] = invJ * ∇N[a*3-2:a*3-(3-dim), k]
                        end

                        Nk = N[:,k]
                        h_gp    = dot(pa_h[elem][:,1], Nk)
                        dhdx_gp = dot(pa_dhdx[elem][:,1], Nk)

                        w = coef * jacDet[k] * intWeights[k]

                        for a in 1:numNodes, b in 1:numNodes
                            val = Nk[a] * dN[dir,b] * h_gp -
                                  Nk[a] * Nk[b] * dhdx_gp
                            _add_blockdiag!(Ke, pdim, a, b, val * w)
                        end
                    end

                    for a in 1:pdim*numNodes
                        na = div(a-1,pdim)+1
                        Ia = (nodes[na]-1)*pdim + mod(a-1,pdim)+1
                        for b in 1:pdim*numNodes
                            nb = div(b-1,pdim)+1
                            Jb = (nodes[nb]-1)*pdim + mod(b-1,pdim)+1
                            I[pos]=Ia; J[pos]=Jb; V[pos]=Ke[a,b]
                            pos+=1
                        end
                    end
                end
            end
        end
    end

    resize!(I,pos-1); resize!(J,pos-1); resize!(V,pos-1)
    dof = pdim * problem.non
    C = sparse(I,J,V,dof,dof)
    dropzeros!(C)

    return SystemMatrix(C, problem)
end

