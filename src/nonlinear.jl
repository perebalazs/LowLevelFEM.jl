export nodePositionVector, ∇, curl, rot, div, grad
export tangentMatrixConstitutive, tangentMatrixInitialStress
export equivalentNodalForce, nonFollowerLoadVector
export applyDeformationBoundaryConditions!, suppressDeformationAtBoundaries!, suppressDeformationAtBoundaries
export solveDeformation, showDeformationResults

"""
    nodePositionVector(problem)

Returns the position vectors of all mesh nodes as a `VectorField` (initial configuration).

Returns: `R`

Types:
- `problem`: Problem
- `R`: VectorField
"""
function nodePositionVector(problem)
    dim = problem.dim
    non = problem.non
    if dim != 3
        error("nodePositionVector: This function only works for 3D problems.")
    end
    r = zeros(non * dim)
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            #nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(edim, -1, true, false)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes()
            r[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            r[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            r[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
        end
    end
    return VectorField([], reshape(r, :,1), [0], [], 1, :v3D, problem)
end

"""
    ∇(rr::Union{VectorField, ScalarField, TensorField}; nabla=:grad)

Computes derivatives of `rr`.
- If `rr` is a `ScalarField` and `nabla == :grad`, returns the gradient (a `VectorField`).
- If `rr` is a `VectorField` and `nabla == :grad`, returns the gradient (a `TensorField`).
- If `rr` is a `VectorField` and `nabla == :curl`, returns the curl (a `VectorField`).
- If `rr` is a `VectorField` and `nabla == :div`, returns the divergence (a `ScalarField`).
- If `rr` is a `TensorField` and `nabla == :div`, returns the divergence (a `VectorField`).

Returns: `ScalarField`, `VectorField`, or `TensorField`

Types:
- `rr`: `ScalarField`, `VectorField`, or `TensorField`
- `nabla`: Symbol

# 3D Examples (assumes `problem` is set as in the ∇ doc setup)

```julia
# One-time 3D setup (assumes examples/Fields/cube.geo exists with physical group "body")
using LowLevelFEM
gmsh.initialize()
gmsh.open("examples/Fields/cube.geo")
mat = material("body", E=210e3, ν=0.3, ρ=7.85e-9)
problem = Problem([mat], type=:Solid)

# 1) Gradient of a 3D scalar field: ∇f → VectorField
f(X,Y,Z) = X^2 + Y*Z
S = scalarField(problem, [field("body", f=f)])
G = ∇(S)  # VectorField with 3 components

# 2) Curl of a 3D vector field: ∇ × v → VectorField
vx(X,Y,Z) = 0
vy(X,Y,Z) = X
vz(X,Y,Z) = 0
V = vectorField(problem, [field("body", fx=vx, fy=vy, fz=vz)])
C = ∇(V, nabla=:curl)  # approx (0, 0, 1) everywhere

# 3) Divergence of a 3D vector field: ∇ ⋅ v → ScalarField
v1(X,Y,Z) = X
v2(X,Y,Z) = Y
v3(X,Y,Z) = Z
V2 = vectorField(problem, [field("body", fx=v1, fy=v2, fz=v3)])
D = ∇(V2, nabla=:div)  # ≈ 3

# 4) Divergence of a 3D tensor field: ∇ · T → VectorField (if T is TensorField)
# For example, a diagonal tensor T with only Tzz = g(Z): div(T) = (0, 0, ∂g/∂Z)
g(Z) = 10 - Z
T = tensorField(problem, [field("body", fz=g)])
DV = ∇(T, nabla=:div)  # VectorField

# Symmetric displacement gradient via operators
# A = (u ∘ ∇ + ∇ ∘ u) / 2
gmsh.finalize()
```
"""
function ∇(rr::Union{VectorField, ScalarField, TensorField}; nabla=:grad)
    problem = rr.model
    gmsh.model.setCurrent(problem.name)
    if rr.a == [;;]
        r = elementsToNodes(rr)
    else
        r = rr
    end
    DoFResults=false
    
    nsteps = r.nsteps
    ε = []
    numElem = Int[]
    ncoord2 = zeros(3 * problem.non)
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    type = :e
    if DoFResults == true
        E1 = zeros(non * 9, nsteps)
        pcs = zeros(Int64, non * dim)
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        ν = problem.material[ipg].ν
        dim = 0
        if problem.dim == 3 && r isa VectorField && nabla == :grad
            dim = 3
            pdim = 3
            rowsOfB = 9
            b = 1
            sz = 9
        elseif problem.dim == 3 && r isa VectorField && nabla == :div
            dim = 3
            pdim = 3
            rowsOfB = 9
            b = 1
            sz = 1
        elseif problem.dim == 3 && r isa VectorField && nabla == :curl
            dim = 3
            pdim = 3
            rowsOfB = 9
            b = 1
            sz = 3
        elseif problem.dim == 3 && r isa ScalarField
            dim = 3
            pdim = 1
            rowsOfB = 3
            b = 1
            sz = 3
        elseif r isa TensorField && nabla == :div
            dim = 3
            pdim = 9
            rowsOfB = 3
            b = 1
            sz = 3
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
            dim = 2
            rowsOfB = 4
            b = 1
        else
            error("∇: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            #nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes()
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                #e0 = zeros(rowsOfB * numNodes)
                nodeCoord = zeros(numNodes * 3)
                for k in 1:dim, j = 1:numNodes
                    nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
                end
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
                ∇h = reshape(dfun, :, numNodes)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
                h = reshape(fun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numNodes)
                ∂h = zeros(3, numNodes * numNodes)
                B = zeros(rowsOfB * numNodes, pdim * numNodes) # spzeros????
                nn2 = zeros(Int, pdim * numNodes)
                rr = zeros(numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numNodes
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        rr[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numNodes, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
                    end
                    B .*= 0
                    if dim == 3 && r isa VectorField
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*9-(9-1), l*3-(3-1)] = B[k*9-(9-2), l*3-(3-2)] = B[k*9-(9-3), l*3-(3-3)] = ∂h[1, (k-1)*numNodes+l]
                            B[k*9-(9-4), l*3-(3-1)] = B[k*9-(9-5), l*3-(3-2)] = B[k*9-(9-6), l*3-(3-3)] = ∂h[2, (k-1)*numNodes+l]
                            B[k*9-(9-7), l*3-(3-1)] = B[k*9-(9-8), l*3-(3-2)] = B[k*9-(9-9), l*3-(3-3)] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 3 && r isa ScalarField
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*3-(3-1), l] = ∂h[1, (k-1)*numNodes+l]
                            B[k*3-(3-2), l] = ∂h[2, (k-1)*numNodes+l]
                            B[k*3-(3-3), l] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif r isa TensorField && nabla == :div
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*3-(3-1), l*9-(9-1)] = B[k*3-(3-2), l*9-(9-2)] = B[k*3-(3-3), l*9-(9-3)] = ∂h[1, (k-1)*numNodes+l]
                            B[k*3-(3-1), l*9-(9-4)] = B[k*3-(3-2), l*9-(9-5)] = B[k*3-(3-3), l*9-(9-6)] = ∂h[2, (k-1)*numNodes+l]
                            B[k*3-(3-1), l*9-(9-7)] = B[k*3-(3-2), l*9-(9-8)] = B[k*3-(3-3), l*9-(9-9)] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif dim == 2 && rowsOfB == 4
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*4-2, l*2-1] = rr[k] < 1e-10 ? 0 : h[l, k] / rr[k]
                        end
                    else
                        error("solveStrain: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    e = zeros(sz * numNodes, nsteps) # tensors have nine elements
                    for k in 1:numNodes
                        if rowsOfB == 9 && dim == 3 && r isa VectorField
                            B1 = B[k*9-8:k*9, 1:3*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    if nabla == :grad
                                        e[(k-1)*9+1:k*9, kk] = [e0[1], e0[2], e0[3],
                                            e0[4], e0[5], e0[6],
                                            e0[7], e0[8], e0[9]]
                                    elseif nabla == :div
                                        e[k, kk] = e0[1] + e0[5] + e0[9]
                                    elseif nabla == :curl
                                        e[(k-1)*3+1:k*3, kk] = [e0[6] - e0[8],
                                            e0[7] - e0[3],
                                            e0[2] - e0[4]]
                                    end
                                end
                                if DoFResults == true
                                    if nabla == :grad
                                        E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[2], e0[3], e0[4], e0[5], e0[6], e0[7], e0[8], e0[9]]
                                    elseif nabla == :div
                                        E1[nnet[j, k], kk] .+= e0[1] + e0[5] + e0[9]
                                    elseif nabla == :curl
                                        E1[3*nnet[j, k]-2:3*nnet[j,k], kk] .+= [e0[6] - e0[8], e0[7] - e0[3], e0[2] - e0[4]]
                                    end
                                end
                            end
                        elseif rowsOfB == 3 && dim == 3 && r isa ScalarField && nabla == :grad
                            B1 = B[k*3-2:k*3, 1:numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*3+1:k*3, kk] = [e0[1], e0[2], e0[3]]
                                end
                                if DoFResults == true
                                    E1[3*nnet[j, k]-2:3*nnet[j,k], kk] .+= [e0[1], e0[2], e0[3]]
                                end
                            end
                        elseif r isa TensorField && nabla == :div
                            B1 = B[k*3-2:k*3, 1:9*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*3+1:k*3, kk] = [e0[1], e0[2], e0[3]]
                                end
                                if DoFResults == true
                                    E1[3*nnet[j, k]-2:3*nnet[j,k], kk] .+= [e0[1], e0[2], e0[3]]
                                end
                            end
                        elseif rowsOfB == 4 && dim == 2 && problem.type == :AxiSymmetric
                            B1 = B[k*4-3:k*4, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[4]/2, 0,
                                        e0[4]/2, e0[3], 0,
                                        0, 0, e0[2]]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[4], 0, e0[4], e0[3], 0, 0, 0, e0[2]]
                                end
                            end
                        else
                            error("solveStrain: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
                        end
                    end
                    if DoFResults == true
                        pcs[nnet[j,1:numNodes]] .+= 1
                    end
                    if DoFResults == false
                        push!(ε, e)
                    end
                end
            end
        end
    end
    if DoFResults == true
        for k in 1:rowsOfB
            for l in 1:non
                E1[k + rowsOfB * l - rowsOfB, :] ./= pcs[l]
            end
        end
    end
    if DoFResults == true
        if r isa VectorField && nabla == :grad
            return TensorField([], E1, r.t, [], nsteps, :e, problem)
        elseif r isa VectorField && nabla == :div
            return ScalarField([], E1, r.t, [], nsteps, :scalar, problem)
        elseif r isa VectorField && nabla == :curl
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        elseif r isa ScalarField
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        elseif r isa TensorField && nabla == :div
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        end
    else
        if r isa VectorField && nabla == :grad
            type = :tensor
            return TensorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        elseif r isa VectorField && nabla == :div
            type = :scalar
            return ScalarField(ε, [;;], r.t, numElem, nsteps, type, problem)
        elseif r isa VectorField && nabla == :curl
            type = :v3D
            return VectorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        elseif r isa ScalarField
            type = :v3D
            return VectorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        elseif r isa TensorField && nabla == :div
            type = :v3D
            return VectorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        end
    end
end

"""
    curl(r::VectorField)

Solves the rotation of the vector field `r`.
An alternative way to solve `curl` is to use `∇` as a differencial operator.

Return: VectorField

Types:
- `r`: VectorField

# 3D Example (assumes `problem` is set as in the ∇ doc setup)

```julia
# Assumes a 3D mesh with physical group "body".
vx(X, Y, Z) = 0
vy(X, Y, Z) = X
vz(X, Y, Z) = 0
v = vectorField(problem, [field("body", fx=vx, fy=vy, fz=vz)])
D1 = curl(v)
D2 = ∇ × v
```
"""
function curl(r::VectorField)
    return ∇(r, nabla=:curl)
end

"""
    rot(r::VectorField)

Solves the rotation of the vector field `r`. In some countries "rot" denotes the English "curl".
(See the `curl` function.)

Return: VectorField

Types:
- `r`: VectorField
"""
function rot(r::VectorField)
    return ∇(r, nabla=:curl)
end

import Base.div

"""
    div(r::Union{VectorField,TensorField})

Solves the divergence of the vector field or tensor field `r`.
An alternative way to solve `div` is to use `∇` as a differencial operator.

Return: ScalarField or VectorField

Types:

- `r`: VectorField or TensorField

# 3D Examples (assumes `problem` is set as in the ∇ doc setup)

```julia
# Assumes a 3D mesh with physical group "body".

# 1) Divergence of a 3D vector field → ScalarField
v1(X,Y,Z) = X
v2(X,Y,Z) = Y
v3(X,Y,Z) = Z
v = vectorField(problem, [field("body", fx=v1, fy=v2, fz=v3)])
D1 = div(v)
D2 = ∇ ⋅ v

# 2) Divergence of a 3D tensor field → VectorField
fsz(X, Y, Z) = 10 - Z
S = tensorField(problem, [field("body", fz=fsz)])
b1 = -div(S)
b2 = -S ⋅ ∇
```
"""
function div(r::Union{VectorField,TensorField})
    return ∇(r, nabla=:div)
end

"""
    grad(r::Union{ScalarField,VectorField})

Solves the gradient of the scalar field or vector field `r`.
An alternative way to solve `grad` is to use `∇` as a differencial operator.

Return: VectorField or TensorField

Types:

- `r`: ScalarField or VectorField

# 3D Examples

```julia
# Assumes a 3D mesh with physical group "body".

# 1) Gradient of a 3D scalar field → VectorField
f(X,Y,Z) = X^2 + Y*Z
S = scalarField(problem, [field("body", f=f)])
G1 = grad(S)
G2 = ∇(S)

# 2) Gradient of a 3D vector field → TensorField
vx(X,Y,Z) = X
vy(X,Y,Z) = Y
vz(X,Y,Z) = Z
V = vectorField(problem, [field("body", fx=vx, fy=vy, fz=vz)])
T1 = grad(V)
T2 = V ∘ ∇
```
"""
function grad(r::Union{VectorField,ScalarField})
    return ∇(r)
end

"""
    tangentMatrixConstitutive(r::VectorField)

Solves the constitutive part of the tangent matrix (when solving large deformation problems).
(See [^6]) `r` is the position vector field in the current configuration.

Return: SystemMatrix

Types:
    - `r`: VectorField

[^6]: Javier Bonet, Richard D. Wood: *Nonlinear Continuum Mechanics for Finite Element Analysis*, 
    Cambridge University Press, 2008, <https://doi.org/10.1017/CBO9780511755446>
"""
function tangentMatrixConstitutive(r::VectorField)
    problem = r.model
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
    nn = []
    I = []
    J = []
    V = []
    V = convert(Vector{Float64}, V)
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)
    II1 = zeros(6, 6)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        mtype = problem.material[ipg].type
        μ = problem.material[ipg].μ
        λ = problem.material[ipg].λ
        κ = problem.material[ipg].κ
        dim = problem.dim
        pdim = problem.pdim
        if mtype == :StVenantKirchhoff
            rowsOfB = 6
            b = 1
        elseif mtype == :NeoHookeCompressible
            rowsOfB = 6
            b = 1
        else
            error("stiffnessMatrixLinear: dimension is $(problem.dim), material type is $(mtype).")
        end

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
                #comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                #h = reshape(fun, :, numIntPoints)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                #H = zeros(9 * numIntPoints, pdim * numNodes)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                ∂H = spzeros(9 * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                sizehint!(∂H, 9 * numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r.a[nn2]
                    B .*= 0
                    ∂H .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            end
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                ∂H[k*9-(9-1), l*pdim-(pdim-1)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-2), l*pdim-(pdim-2)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-3), l*pdim-(pdim-3)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-4), l*pdim-(pdim-1)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-5), l*pdim-(pdim-2)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-6), l*pdim-(pdim-3)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-7), l*pdim-(pdim-1)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-8), l*pdim-(pdim-2)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-9), l*pdim-(pdim-3)] = ∂h[3, (k-1)*numNodes+l]
                            end
                            ∂H1 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F = ∂H1 * r1
                            for l in 1:numNodes
                                B[k*rowsOfB-5, l*pdim-2] = F[1] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-2] = F[4] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-2] = F[7] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-2] = (F[4] * ∂h[1, (k-1)*numNodes+l] + F[1] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-2] = (F[7] * ∂h[2, (k-1)*numNodes+l] + F[4] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-2] = (F[1] * ∂h[3, (k-1)*numNodes+l] + F[7] * ∂h[1, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-5, l*pdim-1] = F[2] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-1] = F[5] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-1] = F[8] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-1] = (F[5] * ∂h[1, (k-1)*numNodes+l] + F[2] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-1] = (F[8] * ∂h[2, (k-1)*numNodes+l] + F[5] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-1] = (F[2] * ∂h[3, (k-1)*numNodes+l] + F[8] * ∂h[1, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-5, l*pdim-0] = F[3] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-0] = F[6] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-0] = F[9] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-0] = (F[6] * ∂h[1, (k-1)*numNodes+l] + F[3] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-0] = (F[9] * ∂h[2, (k-1)*numNodes+l] + F[6] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-0] = (F[3] * ∂h[3, (k-1)*numNodes+l] + F[9] * ∂h[1, (k-1)*numNodes+l]) / 2
                            end
                        end
                    else
                        error("stiffnessMatrixLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        if mtype == :StVenantKirchhoff
                            λ0 = λ
                            μ0 = μ
                            C1 = [λ0+2μ0 λ0 λ0 0 0 0;
                                λ0 λ0+2μ0 λ0 0 0 0;
                                λ0 λ0 λ0+2μ0 0 0 0;
                                0 0 0 μ0 0 0;
                                0 0 0 0 μ0 0;
                                0 0 0 0 0 μ0]
                        elseif mtype == :NeoHookeCompressible
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            J1 = det(F2)
                            C2 = F2' * F2
                            iC2 = inv(C2)
                            iC1 = reshape(iC2, 9, 1)
                            iCiC1 = iC1[[1, 5, 9, 4, 6, 3]] * iC1[[1, 5, 9, 4, 6, 3]]'
                            i1 = [1, 2, 3, 1, 2, 3]
                            j1 = [1, 2, 3, 2, 3, 1]
                            k1 = [1, 2, 3, 1, 2, 3]
                            l1 = [1, 2, 3, 2, 3, 1]
                            for i2 in 1:6, j2 in 1:6
                                II1[i2, j2] = (iC2[i1[i2], k1[j2]] * iC2[j1[i2], l1[j2]] + iC2[i1[i2], l1[j2]] * iC2[j1[i2], k1[j2]]) / 2
                            end
                            #J1 = J1 > 0 ? J1 : 1e-14
                            C1 = λ * iCiC1 + 2 * (μ - λ * log(J1)) * II1
                        else
                            error("loadVectorNonLinear: type '$(problem.type)' is not defined.")
                        end
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += B1' * C1 * B1 * b * jacDet[k] * intWeights[k]
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end

"""
    tangentMatrixInitialStress(r::VectorField)

Solves the initial stress part of the tangent matrix (when solving large deformation problems).
(See [^6]) `r` is the position vector field in the current configuration.

Return: SystemMatrix

Types:

- `r`: VectorField
"""
function tangentMatrixInitialStress(r::VectorField)
    problem = r.model
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
    nn = []
    I = []
    J = []
    V = []
    V = convert(Vector{Float64}, V)
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        mtype = problem.material[ipg].type
        μ = problem.material[ipg].μ
        λ = problem.material[ipg].λ
        κ = problem.material[ipg].κ
        dim = problem.dim
        pdim = problem.pdim
        if mtype == :StVenantKirchhoff
            rowsOfB = 6
            b = 1
        elseif mtype == :NeoHookeCompressible
            rowsOfB = 6
            b = 1
        else
            error("stiffnessMatrixNonLinear: dimension is $(problem.dim), material type is $(mtype).")
        end

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
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                H = spzeros(9 * numIntPoints, 9 * numNodes)
                sizehint!(H, 9 * numNodes * numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                ∂H = spzeros(9 * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, dim * numNodes)
                nn9 = zeros(Int, 9 * numNodes)
                sizehint!(∂H, 9 * numNodes * numIntPoints)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:9
                        H[k*9-(9-kk), l*9-(9-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    for k in 1:9
                        nn9[k:9:9*numNodes] = 9 * nnet[j, 1:numNodes] .- (9 - k)
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r.a[nn2]
                    B .*= 0
                    ∂H .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            end
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                ∂H[k*9-(9-1), l*pdim-(pdim-1)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-2), l*pdim-(pdim-2)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-3), l*pdim-(pdim-3)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-4), l*pdim-(pdim-1)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-5), l*pdim-(pdim-2)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-6), l*pdim-(pdim-3)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-7), l*pdim-(pdim-1)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-8), l*pdim-(pdim-2)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-9), l*pdim-(pdim-3)] = ∂h[3, (k-1)*numNodes+l]
                            end
                        end
                    else
                        error("stiffnessMatrixNonLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        if mtype == :StVenantKirchhoff
                            I3 = [1 0 0; 0 1 0; 0 0 1]
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            E = (F2' * F2 - I3) / 2
                            S2 = 2μ * E + λ * (E[1, 1] + E[2, 2] + E[3, 3]) * I3
                        elseif mtype == :NeoHookeCompressible
                            I3 = [1 0 0; 0 1 0; 0 0 1]
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            J1 = det(F2)
                            C = F2' * F2
                            iC = inv(C)
                            #J1 = J1 > 0 ? J1 : 1e-14
                            S2 = μ * (I3 - iC) + λ * log(J1) * iC
                        else
                            error("loadVectorNonLinear: type '$(problem.type)' is not defined.")
                        end
                        ∂H1 = ∂H[k*9-(9-1):3:k*9-(9-9), 1:pdim*numNodes]
                        ∂H2 = ∂H[k*9-(9-2):3:k*9-(9-9), 1:pdim*numNodes]
                        ∂H3 = ∂H[k*9-(9-3):3:k*9-(9-9), 1:pdim*numNodes]
                        K1 += (∂H1' * S2 * ∂H1 + ∂H2' * S2 * ∂H2 + ∂H3' * S2 * ∂H3) * b * jacDet[k] * intWeights[k]
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end

"""
    equivalentNodalForce(r::VectorField)

Solves the equivalent nodal force (when solving large deformation problems).
(See [^6]) `r` is the position vector field in the current configuration.

Return: VectorField

Types:

- `r`: VectorField
"""
function equivalentNodalForce(r::VectorField)
    problem = r.model
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    f = zeros(problem.non * problem.pdim)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        mtype = problem.material[ipg].type
        μ = problem.material[ipg].μ
        λ = problem.material[ipg].λ
        κ = problem.material[ipg].κ
        dim = problem.dim
        pdim = problem.pdim
        if mtype == :StVenantKirchhoff
            rowsOfB = 6
            b = 1
        elseif mtype == :NeoHookeCompressible
            rowsOfB = 6
            b = 1
        else
            error("loadVectorNonLinear: dimension is $(problem.dim), material type is $(problem.type).")
        end

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
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                H = spzeros(9 * numIntPoints, 9 * numNodes)
                sizehint!(H, 9 * numNodes * numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                ∂H = spzeros(9 * numIntPoints, pdim * numNodes)
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn9 = zeros(Int, 9 * numNodes)
                sizehint!(∂H, 9 * numNodes * numIntPoints)
                #for k in 1:numIntPoints, l in 1:numNodes
                #    for kk in 1:3:9, ll in 1:3
                #H[k*9-(9-kk)-1+ll, l*pdim-(pdim-ll)] = h[(k-1)*numNodes+l]
                #    end
                #end
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    #for k in 1:9
                    #    nn9[k:9:9*numNodes] = 9 * nnet[j, 1:numNodes] .- (9 - k)
                    #end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r.a[nn2]
                    B .*= 0
                    ∂H .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            end
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                ∂H[k*9-(9-1), l*pdim-(pdim-1)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-2), l*pdim-(pdim-2)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-3), l*pdim-(pdim-3)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-4), l*pdim-(pdim-1)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-5), l*pdim-(pdim-2)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-6), l*pdim-(pdim-3)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-7), l*pdim-(pdim-1)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-8), l*pdim-(pdim-2)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-9), l*pdim-(pdim-3)] = ∂h[3, (k-1)*numNodes+l]
                            end
                            ∂H1 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F = ∂H1 * r1
                            for l in 1:numNodes
                                B[k*rowsOfB-5, l*pdim-2] = F[1] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-2] = F[4] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-2] = F[7] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-2] = (F[4] * ∂h[1, (k-1)*numNodes+l] + F[1] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-2] = (F[7] * ∂h[2, (k-1)*numNodes+l] + F[4] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-2] = (F[1] * ∂h[3, (k-1)*numNodes+l] + F[7] * ∂h[1, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-5, l*pdim-1] = F[2] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-1] = F[5] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-1] = F[8] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-1] = (F[5] * ∂h[1, (k-1)*numNodes+l] + F[2] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-1] = (F[8] * ∂h[2, (k-1)*numNodes+l] + F[5] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-1] = (F[2] * ∂h[3, (k-1)*numNodes+l] + F[8] * ∂h[1, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-5, l*pdim-0] = F[3] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-0] = F[6] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-0] = F[9] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-0] = (F[6] * ∂h[1, (k-1)*numNodes+l] + F[3] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-0] = (F[9] * ∂h[2, (k-1)*numNodes+l] + F[6] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-0] = (F[3] * ∂h[3, (k-1)*numNodes+l] + F[9] * ∂h[1, (k-1)*numNodes+l]) / 2
                            end
                        end
                    else
                        error("loadVectorNonLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    f1 .*= 0
                    for k in 1:numIntPoints
                        if mtype == :StVenantKirchhoff
                            I3 = [1 0 0; 0 1 0; 0 0 1]
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            E = (F2' * F2 - I3) / 2
                            S2 = 2μ * E + λ * (E[1, 1] + E[2, 2] + E[3, 3]) * I3
                        elseif mtype == :NeoHookeCompressible
                            I3 = [1 0 0; 0 1 0; 0 0 1]
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            J1 = det(F2)
                            C = F2' * F2
                            iC = inv(C)
                            #J1 = J1 > 0 ? J1 : 1e-14
                            S2 = μ * (I3 - iC) + λ * log(J1) * iC
                        else
                            error("loadVectorNonLinear: type '$(problem.type)' is not defined.")
                        end
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        f1 += B1' * S2[[1, 5, 9, 4, 6, 3]] * b * jacDet[k] * intWeights[k]
                    end
                    f[nn2] += f1
                end
            end
        end
    end
    return VectorField([], reshape(f, :,1), [0], [], 1, :v3D, problem)
end

"""
    nonFollowerLoadVector(r::VectorField, load)

Solves the non-follower load vector (when solving large deformation problems).
`r` is the position vector field in the current configuration.

Return: VectorField

Types:

- `r`: VectorField
"""
function nonFollowerLoadVector(r::VectorField, loads)
    problem = r.model
    gmsh.model.setCurrent(problem.name)
    if !isa(loads, Vector)
        error("loadVector: loads are not arranged in a vector. Put them in [...]")
    end
    pdim = problem.pdim
    DIM = problem.dim
    b = problem.thickness
    non = problem.non
    dof = pdim * non
    fp = zeros(dof)
    ncoord2 = zeros(3 * problem.non)
    F0 = ∇(r)
    F0 = elementsToNodes(F0)
    for n in 1:length(loads)
        name, fx, fy, fz = loads[n]
        if pdim == 3
            f = [.0, .0, .0]
        elseif pdim == 2
            f = [.0, .0]
        elseif pdim == 1
            f = [.0]
        else
            error("loadVector: dimension of the problem is $(problem.dim).")
        end
        dimTags = gmsh.model.getEntitiesForPhysicalName(name)
        for i ∈ 1:length(dimTags)
            dimTag = dimTags[i]
            dim = dimTag[1]
            tag = dimTag[2]
            elementTypes, elementTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
            nodeTags::Vector{Int64}, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, tag, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for ii in 1:length(elementTypes)
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementTypes[ii])
                nnoe = reshape(elemNodeTags[ii], numNodes, :)'
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(2order+1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elementTags[ii]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                HH = zeros(9 * numIntPoints, 9 * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                        end
                        for l in 1:9
                            HH[j*9-(9-l), k*9-(9-l)] = h[(j-1)*numNodes+k]
                        end
                    end
                end
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn9 = zeros(Int, 9 * numNodes)
                ∂h = zeros(pdim, numNodes * numIntPoints)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:9
                        nn9[k:9:9*numNodes] = 9 * nnet[l, 1:numNodes] .- (9 - k)
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    f1 .*= 0
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r.a[nn2]
                    F9 = F0.a[nn9]
                    for j in 1:numIntPoints
                        H9 = HH[j*9-(9-1):j*9, 1:9*numNodes]
                        F1 = H9 * F9
                        F1 = reshape(F1, 3, 3)
                        iF1 =inv(F1)
                        x = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 2]
                        y = 0
                        z = 0
                        if isa(fx, Function) || isa(fy, Function) || isa(fz, Function)
                            y = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 1]
                            z = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 0]
                        end
                        if fz == 2im
                            if isa(fx, Function)
                                error("heatConvectionVector: h cannot be a function.")
                            end
                            f[1] = isa(fy, Function) ? fx * fy(x, y, z) : fx * fy
                        else
                            f[1] = isa(fx, Function) ? fx(x, y, z) : fx
                        end
                        if pdim > 1
                            f[2] = isa(fy, Function) ? fy(x, y, z) : fy
                        end
                        if pdim == 3
                            f[3] = isa(fz, Function) ? fz(x, y, z) : fz
                        end
                        H1 = H[j*pdim-(pdim-1):j*pdim, 1:pdim*numNodes] # H1[...] .= H[...] ????
                        ############### NANSON ######## 3D ###################################
                        if DIM == 3 && dim == 3
                            Ja = jacDet[j]
                        elseif DIM == 3 && dim == 2
                            xy = Jac[1, 3*j-2] * Jac[2, 3*j-1] - Jac[2, 3*j-2] * Jac[1, 3*j-1]
                            yz = Jac[2, 3*j-2] * Jac[3, 3*j-1] - Jac[3, 3*j-2] * Jac[2, 3*j-1]
                            zx = Jac[3, 3*j-2] * Jac[1, 3*j-1] - Jac[1, 3*j-2] * Jac[3, 3*j-1]
                            Ja = √(xy^2 + yz^2 + zx^2)
                        elseif DIM == 3 && dim == 1
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2 + (Jac[3, 3*j-2])^2)
                        elseif DIM == 3 && dim == 0
                            Ja = 1
                        ############ 2D #######################################################
                        elseif DIM == 2 && dim == 2 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = jacDet[j] * b
                        elseif DIM == 2 && dim == 2 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * jacDet[j] * x
                        elseif DIM == 2 && dim == 1 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * x
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
                        else
                            error("loadVector: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                        end
                        f1 += H1' * iF1 * f * Ja * intWeights[j]
                    end
                    fp[nn2] += f1
                end
            end
        end
    end
    return VectorField([], reshape(fp, :,1), [0], [], 1, :v3D, problem)
end

"""
    applyDeformationBoundaryConditions!(deformVec, supports)

Applies displacement boundary conditions `supports` on deformation vector `deformVec`.
Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Returns: nothing

Types:
- `deformVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyDeformationBoundaryConditions!(deformVec, supports; fact=1.0)
    problem = deformVec.model
    if !isa(supports, Vector)
        error("applyDeformationBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    dof, st = size(deformVec.a)
    pdim = problem.pdim

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsX
                jj += 1
                if isa(ux, Function)
                    deformVec.a[j] += uux[jj] * fact
                else
                    deformVec.a[j] += ux * fact
                end
            end
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsY
                jj += 1
                if isa(uy, Function)
                    deformVec.a[j] += uuy[jj] * fact
                else
                    deformVec.a[j] += uy * fact
                end
            end
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsZ
                jj += 1
                if isa(uz, Function)
                    deformVec.a[j] += uuz[jj] * fact
                else
                    deformVec.a[j] += uz * fact
                end
            end
        end
    end
end

"""
    suppressDeformationAtBoundaries!(stiffMat, loadVec, supports)

Suppresses the displacements given in `support` in `stiffMat` and `loadVec` 
so that it is only necessary to consider them once during iteration.
`stiffMat` is the stiffness matrix, `loadVec` is the load vector.
`supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Returns: nothing

Types:
- `stiffMat`: SystemMatrix 
- `loadVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function suppressDeformationAtBoundaries!(stiffMat, loadVec, supports)
#function applyDeformationBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports; fix=1)
    if stiffMat.model != loadVec.model
        error("suppressDeformationAtBoundaries!: system matrix and load vector does not belong to the same problem.")
    end
    problem = stiffMat.model
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    dof, dof = size(stiffMat.A)
    pdim = problem.pdim

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim - 1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
                f0 = stiffMat.A[:, nodeTagsX] * uux * 0
            else
                f0 = stiffMat.A[:, nodeTagsX] * ux * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim - 2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
                f0 = stiffMat.A[:, nodeTagsY] * uuy * 0
            else
                f0 = stiffMat.A[:, nodeTagsY] * uy * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
                f0 = stiffMat.A[:, nodeTagsZ] * uuz * 0
            else
                f0 = stiffMat.A[:, nodeTagsZ] * uz * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
    end

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsX
                jj += 1
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(ux, Function)
                    loadVec.a[j] = uux[jj] * 0
                else
                    loadVec.a[j] = ux * 0
                end
            end
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsY
                jj += 1
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(uy, Function)
                    loadVec.a[j] = uuy[jj] * 0
                else
                    loadVec.a[j] = uy * 0
                end
            end
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsZ
                jj += 1
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(uz, Function)
                    loadVec.a[j] = uuz[jj] * 0
                else
                    loadVec.a[j] = uz * 0
                end
            end
        end
    end

    dropzeros!(stiffMat.A)
    #dropzeros!(massMat)
    #dropzeros!(dampMat)
end

"""
    suppressDeformationAtBoundaries(stiffMat, loadVec, supports)

Suppresses the displacements given in `support` in `stiffMat` and `loadVec` 
so that it is only necessary to consider them once during iteration.
`stiffMat` is the stiffness matrix, `loadVec` is the load vector.
`supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Return: stiffMat1, loadVec1

Types:
- `stiffMat`: SystemMatrix 
- `loadVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
- `stiffMat1`: SystemMatrix 
- `loadVec1`: VectorField
"""
function suppressDeformationAtBoundaries(stiffMat, loadVec, supports)
    if stiffMat.model != loadVec.model
        error("suppressDeformationAtBoundaries!: system matrix and load vector does not belong to the same problem.")
    end
    problem = stiffMat.model
    K1 = copy(stiffMat)
    f1 = copy(loadVec)
    suppressDeformationAtBoundaries!(K1, f1, supports)
    return K1, f1
end

"""
    solveDeformation(problem::Problem, load, supp;
                    followerLoad=false,
                    loadSteps = 3,
                    rampedLoad = true,
                    rampedSupport = false,
                    maxIteration = 10,
                    saveSteps = false,
                    saveIterations = false,
                    plotConvergence = false,
                    relativeError = 1e-5,
                    initialDeformation=nodePositionVector(problem))

Solves the deformed shape of a non-linearly elastic body...
"""
function solveDeformation(problem, load, supp;
    followerLoad=false,
    loadSteps = 3,
    rampedLoad = true,
    rampedSupport = false,
    maxIteration = 10,
    saveSteps = false,
    saveIterations = false,
    plotConvergence = false,
    relativeError = 1e-5,
    initialDeformation=nodePositionVector(problem))

    r0 = initialDeformation
    err0 = abs(maximum(r0.a) - minimum(r0.a))
    f = loadVector(problem, load)
    r = []
    r1 = []
    push!(r, r0.a)
    r0 = copy(r0)
    e = []

    for j in range(1, loadSteps)
        fact = rampedLoad == true ? j / loadSteps : 1
        if rampedSupport == true
            applyDeformationBoundaryConditions!(r0, supp, fact=1 / loadSteps)
        elseif rampedSupport == false && j == 1
            applyDeformationBoundaryConditions!(r0, supp)
        end
        err = 1
        i = 0
        while err > relativeError && i < maxIteration
            i += 1
       
            Kl = tangentMatrixConstitutive(r0)
            Knl = tangentMatrixInitialStress(r0)
            if followerLoad == false
                f = nonFollowerLoadVector(r0, load)
            end
            fnl = equivalentNodalForce(r0)
            K1, f1 = suppressDeformationAtBoundaries(Kl + Knl, fact * f - fnl, supp)
            q = solveDisplacement(K1, f1)
            r0 += q
            if saveIterations == true
                push!(r, r0.a)
                r0 = copy(r0)
            end
            err = maximum(abs.(q.a)) / err0
            if plotConvergence == true
                append!(e, err)
            end
        end
        if saveSteps == true
            push!(r, r0.a)
            r0 = copy(r0)
        end
    end
    if saveIterations == true || saveSteps == true
        n = length(r)
        r1 = zeros(length(r0.a), n)
        for i in 1:n
            r1[:, i] = r[i]
        end
    else
        n = length(r)
        r1 = zeros(length(r0.a), 1)
        r1[:,1] = r[:,n]
    end
    r1 = VectorField([], r1, 1:size(r1, 2), [], size(r1, 2), :v3D, problem)
    if plotConvergence == true
        return r1, e
    else
        return r1
    end
end

"""
    showDeformationResults(r::VectorField, comp; name=String, visible=Boolean)

Shows deformation result, where `r` contains the position vectors of nodes 
in the *current configuration*.
"""
function showDeformationResults(r::VectorField, comp; name=comp, visible=false)
    problem = r.model
    if r.a == [;;]
        r = elementsToNodes(r)
    end
    r0 = nodePositionVector(problem)
    u = copy(r)
    for i in 1:size(r.a, 2)
        u.a[:, i] = r.a[:, i] - r0.a
    end
    return showDoFResults(u, comp, name=name, visible=visible, factor=1)
end

function showDeformationResults(r::VectorField; name="vector", visible=false)
    return showDeformationResults(r, :vector, name=name, visible=visible)
end