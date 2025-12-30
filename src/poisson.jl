export stiffnessMatrixPoisson

"""
    stiffnessMatrixPoisson(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)

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
function stiffnessMatrixPoisson(problem::Problem; coefficient::Union{Number,ScalarField}=1.0)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
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
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                @inbounds for k in 1:numIntPoints, l in 1:numNodes
                    val = h[(k-1)*numNodes + l]
                    for kk in 1:pdim
                        row = (k-1)*pdim + kk
                        col = (l-1)*pdim + kk
                        H[row, col] = val
                    end
                end
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

