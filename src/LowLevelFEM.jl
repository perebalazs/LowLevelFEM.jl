module LowLevelFEM

using LinearAlgebra, SparseArrays
using Arpack
import gmsh_jll
include(gmsh_jll.gmsh_api)
import .gmsh
export gmsh

"""
    Problem(names; thickness=..., type=..., bandwidth=...)

A structure containing the most important data of the problem. 
- name of the model (in gmsh)
- type of the problem: 3D "Solid", "PlaneStrain" or "PlaneStress"
- bandwidth optimization using built-in `gmsh` function. Possibilities: "RCMK" (default), Hilbert" and "Metis"
- dimension of the problem, determined from `type`
- material constants: Physical group, Young's modulus, Poisson's ratio,
  mass density (in vector of tuples `names`)
- thickness of the plate
- number of nodes (non)

Types:
- `names`: Vector{Touple{String, Float64, Float64, Float64}}
- `type`: String
- `bandwidth`: String
- `dim`: Integer
- `thickness`: Float64
- `non`: Integer
"""
struct Problem
    name::String
    type::String
    dim::Int64
    material::Vector{Tuple{String, Float64, Float64, Float64}}
    thickness::Float64
    non::Int64
    function Problem(mat; thickness=1, type="Solid", bandwidth="RCMK")
        if type == "Solid"
            dim = 3
        elseif type == "PlaneStress"
            dim = 2
        elseif type == "PlaneStrain"
            dim = 2
        else
            error("Problem = $type ????")
        end
        name = gmsh.model.getCurrent()
        gmsh.option.setString("General.GraphicsFontEngine", "Cairo")
        gmsh.option.setString("View.Format", "%.6g")
        
        material = mat
        elemTags = []
        elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim,-1)
        for i in 1:length(elementTags)
            for j in 1:length(elementTags[i])
                push!(elemTags, elementTags[i][j])
	        end
        end

        if bandwidth != "RCMK" && bandwidth != "Hilbert" && bandwidth != "Metis"
            error("Problem: bandwidth can be 'Hilbert', 'Metis', or 'RCMK'. Now it is '$(bandwidth)'")
        end

        method = bandwidth == "none" ? "RCMK" : bandwidth
        oldTags, newTags = gmsh.model.mesh.computeRenumbering(method, elemTags)
        if bandwidth != "none"
            gmsh.model.mesh.renumberNodes(oldTags, newTags)
        end

        non = length(oldTags)
        return new(name, type, dim, material, thickness, non)
    end
end

"""
    StressField(sigma, numElem, nsteps)

A structure containing the data of a stress field. 
- sigma: vector of ElementNodeData type stress data (see gmsh.jl)
- numElem: vector of tags of elements
- nsteps: number of stress fields stored in sigma (for animations).

Types:
- `sigma`: Vector{Matrix{Float64}}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
"""
struct StressField
    sigma::Vector{Matrix{Float64}}
    numElem::Vector{Int}
    nsteps::Int
end

"""
    FEM.material(name; E=2.0e5, ν=0.3, ρ=7.85e-9)

Returns a tuple in which `name` is the name of a physical group, 
`E` is the modulus of elasticity, `ν` Poisson's ratio and `ρ` is
the mass density.

Return: mat

Types:
- `mat`: Tuple(String, Float64, Float64, Float64)
- `name`: String
- `E`: double
- `ν`: double
- `ρ`: double
"""
function material(name; E=2.0e5, ν=0.3, ρ=7.85e-9)
    mat = name, E, ν, ρ
    return mat
end

"""
    FEM.displacementConstraint(name; ux=..., uy=..., uz=...)

Gives the displacement constraints on `name` physical group. At least one `ux`, 
`uy` or `uz` value have to be given (depending on the dimension of the problem).

Return: none

Types:
- `name`: string
- `ux`: double
- `uy`: double
- `uz`: double
"""
function displacementConstraint(name; ux=1im, uy=1im, uz=1im)
    bc0 = name, ux, uy, uz
    return bc0
end

"""
    FEM.load(name; fx=..., fy=..., fz=...)

Gives the intensity of distributed load on `name` physical group. At least one `fx`, 
`fy` or `fz` value have to be given (depending on the dimension of the problem).

Return: none

Types:
- `name`: string
- `ux`: double
- `uy`: double
- `uz`: double
"""
function load(name; fx=0, fy=0, fz=0)
    ld0 = name, fx, fy, fz
    return ld0
end


"""
    FEM.generateMesh(problem, surf, elemSize; approxOrder=1, algorithm=6, quadrangle=0, internalNodes=0)

Obsolate, use gmsh script (.geo) instead.

Return: none

Types:
- ``: x
"""
function generateMesh(surf, elemSize; approxOrder=1, algorithm=6, quadrangle=0, internalNodes=0)
    all = gmsh.model.getEntities(0)
    gmsh.model.mesh.setSize(all, elemSize)
    gmsh.model.mesh.setAlgorithm(2, surf, algorithm)
    gmsh.model.mesh.generate(1)
    gmsh.model.mesh.generate(2)
    if quadrangle
        gmsh.model.mesh.recombine()
    end
    if internalNodes
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
    else
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
    end
    gmsh.model.mesh.setOrder(approxOrder)
end


"""
    FEM.stiffnessMatrix(problem)

Solves the stiffness matrix of the `problem`.

Return: `stiffMat`

Types:
- `problem`: Problem
- `stiffMat`: SparseMatrix
"""
function stiffnessMatrix(problem; elements=[])
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
        phName, E, ν, ρ = problem.material[ipg]
        dim = 0
        if problem.dim == 3 && problem.type == "Solid"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            dim = 3
            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            dim = 2
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            dim = 2
            rowsOfB = 3
            b = 1
        else
            error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
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
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * dim, numNodes * dim)
                Jidx = zeros(Int, numNodes * dim, numNodes * dim)
                for k in 1:numNodes*dim, l in 1:numNodes*dim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, dim * numNodes)
                K1 = zeros(dim * numNodes, dim * numNodes)
                nn2 = zeros(Int, dim * numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*3-0, l*2-0] = B[k*3-2, l*2-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*3-0, l*2-1] = B[k*3-1, l*2-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*6-5, l*3-2] = B[k*6-2, l*3-1] = B[k*6-0, l*3-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*6-4, l*3-1] = B[k*6-2, l*3-2] = B[k*6-1, l*3-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*6-3, l*3-0] = B[k*6-1, l*3-1] = B[k*6-0, l*3-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:dim*numNodes]
                        K1 += B1' * D * B1 * b * jacDet[k] * intWeights[k]
                    end
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.dim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return K
end

"""
    FEM.massMatrix(problem; lumped=...)

Solves the mass matrix of the `problem`. If `lumped` is true, solves lumped mass matrix.

Return: `massMat`

Types:
- `problem`: Problem
- `lumped`: Boolean
- `massMat`: SparseMatrix
"""
function massMatrix(problem; elements=[], lumped=true)
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
        phName, E, ν, ρ = problem.material[ipg]
        dim = 0
        if problem.dim == 3 && problem.type == "Solid"
            dim = 3
            rowsOfH = 3
            b = 1
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            dim = 2
            rowsOfH = 2
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            dim = 2
            rowsOfH = 2
            b = 1
        else
            error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
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
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                Iidx = zeros(Int, numNodes * dim, numNodes * dim)
                Jidx = zeros(Int, numNodes * dim, numNodes * dim)
                for k in 1:numNodes*dim, l in 1:numNodes*dim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                nn2 = zeros(Int, 2 * numNodes)
                H = zeros(rowsOfH * numIntPoints, dim * numNodes)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:dim
                        H[k*dim-(dim-kk), l*2-(dim-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                M1 = zeros(dim * numNodes, dim * numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    M1 .*= 0
                    for k in 1:numIntPoints
                        H1 = H[k*dim-(dim-1):k*dim, 1:dim*numNodes]
                        M1 += H1' * H1 * jacDet[k] * intWeights[k]
                    end
                    M1 *= ρ * b
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, M1[:])
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.dim * problem.non
    M = sparse(I, J, V, dof, dof)
    if lumped == true
        M = spdiagm(vec(sum(M, dims=2))) # lumped mass matrix
    end
    dropzeros!(M)
    return M
end

"""
    FEM.loadVector(problem, loads)

Solves a load vector of `problem`. `loads` is a tuple of name of physical group 
`name`, coordinates `fx`, `fy` and `fz` of the intensity of distributed force.
It can solve traction or body force depending on the problem.
In case of 2D problems and Line physical group means surface force.
In case of 2D problems and Surface physical group means body force.
In case of 3D problems and Line physical group means edge force.
In case of 3D problems and Surface physical group means surface force.
In case of 3D problems and Volume physical group means body force.

Return: `loadVec`

Types:
- `problem`: Problem
- `loads`: Vector{Tuple{String, Float64, Float64, Float64}}
- `loadVec`: Vector
"""
function loadVector(problem, loads)
    gmsh.model.setCurrent(problem.name)
    pdim = problem.dim
    b = problem.thickness
    non = problem.non
    dof = pdim * non
    fp = zeros(dof)
    for n in 1:length(loads)
        name, fx, fy, fz = loads[n]
        if problem.dim == 3
            f = [fx, fy, fz]
        elseif problem.dim == 2
            f = [fx, fy]
        else
            error("applyBoundaryConditions: dimension of the problem is $(problem.dim).")
        end
        dimTags = gmsh.model.getEntitiesForPhysicalName(name)
        for i ∈ 1:length(dimTags)
            dimTag = dimTags[i]
            dim = dimTag[1]
            tag = dimTag[2]
            elementTypes, elementTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
            for ii in 1:length(elementTypes)
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementTypes[ii])
                nnoe = reshape(elemNodeTags[ii], numNodes, :)'
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(order+1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                        end
                    end
                end
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    f1 .*= 0
                    for j in 1:numIntPoints
                        H1 = H[j*pdim-(pdim-1):j*pdim, 1:pdim*numNodes] # H1[...] .= H[...] ????
                        ############### NANSON ###########################################
                        if pdim == 3 && dim == 3
                            Ja = jacDet[j]
                        elseif pdim == 3 && dim == 2
                            xy = Jac[1, 3*j-2] * Jac[2, 3*j-1] - Jac[2, 3*j-2] * Jac[1, 3*j-1]
                            yz = Jac[2, 3*j-2] * Jac[3, 3*j-1] - Jac[3, 3*j-2] * Jac[2, 3*j-1]
                            zx = Jac[3, 3*j-2] * Jac[1, 3*j-1] - Jac[1, 3*j-2] * Jac[3, 3*j-1]
                            Ja = √(xy^2 + yz^2 + zx^2)
                        elseif pdim == 3 && dim == 1
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2 + (Jac[3, 3*j-2])^2)
                        elseif pdim == 2 && dim == 2
                            Ja = jacDet[j] * b
                        elseif pdim == 2 && dim == 1
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        else
                            error("applyBoundaryConditions: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                        end
                        f1 += H1' * f * Ja * intWeights[j]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                    end
                    fp[nn2] += f1
                end
            end
        end
    end
    return fp
end

"""
    FEM.applyBoundaryConditions!(problem, stiffMat, loadVec, supports)

Applies displacement boundary conditions `supports` on a stiffness matrix
`stiffMat` and load vector `loadVec`. Mesh details are in `problem`. `supports`
is a tuple of `name` of physical group and prescribed displacements `ux`, `uy`
and `uz`.

Return: none

Types:
- `problem`: Problem
- `stiffMat`: SparseMatrix 
- `loadVec`: Vector 
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions!(problem, stiffMat, loadVec, supports)
    dof, dof = size(stiffMat)
    massMat = spzeros(dof, dof)
    dampMat = spzeros(dof, dof)
    applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports)
    massMat = []
    dampMat = []
    return
end

"""
    FEM.applyBoundaryConditions(problem, stiffMat, loadVec, supports)

Applies displacement boundary conditions `supports` on a stiffness matrix
`stiffMat` and load vector `loadVec`. Mesh details are in `problem`. `supports`
is a tuple of `name` of physical group and prescribed displacements `ux`, `uy`
and `uz`. Creates a new stiffness matrix and load vector.

Return: `stiffMat`, `loadVec`

Types:
- `problem`: Problem
- `stiffMat`: SparseMatrix 
- `loadVec`: Vector 
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions(problem, stiffMat0, loadVec0, supports)
    dof, dof = size(stiffMat0)
    massMat = spzeros(dof, dof)
    dampMat = spzeros(dof, dof)
    stiffMat = copy(stiffMat0)
    loadVec = copy(loadVec0)
    applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports)
    massMat = []
    dampMat = []
    return stiffMat, loadVec
end

"""
    FEM.getTagForPhysicalName(name)

Returns `tags` of elements of physical group `name`.

Return: `tags`

Types:
- `name`: String
- `tags`: Vector{Integer}
"""
function getTagForPhysicalName(name)
    dimTags = gmsh.model.getPhysicalGroups(-1)
    i = 1
    while gmsh.model.getPhysicalName(dimTags[i][1], dimTags[i][2]) != name
        i += 1
        if i > length(dimTags)
            error("Physical name '$name' does not exist.")
        end
    end
    return dimTags[i][2]
end

"""
    FEM.applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports)

Applies displacement boundary conditions `supports` on a stiffness matrix
`stiffMat`, mass matrix `massMat`, damping matrix `dampMat` and load vector `loadVec`.
Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Return: none

Types:
- `problem`: Problem
- `stiffMat`: SparseMatrix 
- `massMat`: SparseMatrix 
- `dampMat`: SparseMatrix 
- `loadVec`: Vector{Float64}
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports)
    gmsh.model.setCurrent(problem.name)
    dof, dof = size(stiffMat)
    pdim = problem.dim

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim - 1)
            f0 = stiffMat[:, nodeTagsX] * ux
            f0 = sum(f0, dims=2)
            loadVec .-= f0
        end
        if uy != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim - 2)
            f0 = stiffMat[:, nodeTagsX] * uy
            f0 = sum(f0, dims=2)
            loadVec .-= f0
        end
        if pdim == 3 && uz != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= 3
            f0 = stiffMat[:, nodeTagsY] * uz
            f0 = sum(f0, dims=2)
            loadVec .-= f0
        end
    end

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
            for j ∈ nodeTagsX
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                loadVec[j] = ux
            end
        end
        if uy != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-2)
            for j ∈ nodeTagsX
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                loadVec[j] = uy
            end
        end
        if pdim == 3 && uz != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= 3
            for j ∈ nodeTagsY
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                loadVec[j] = uz
            end
        end
    end

    dropzeros!(stiffMat)
    dropzeros!(massMat)
    dropzeros!(dampMat)
end

"""
    FEM.solveDisplacement(K, q)

Solves the equation K*q=f for the displacement vector `q`. `K` is the stiffness Matrix,
`q` is the load vector.

Return: `q`

Types:
- `K`: SparseMatrix 
- `f`: Vector{Float64} 
- `q`: Vector{Float64}
"""
function solveDisplacement(K, f)
    return K \ f
end

"""
    FEM.solveStress(problem, q)

Solves the stress field `S` from displacement vector `q`. Stress field is given
per elements, so it usually contains jumps at the boundary of elements. Details
of mesh is available in `problem`.

Return: `S`

Types:
- `problem`: Problem
- `q`: Vector{Float64}
- `S`: StressField
"""
function solveStress(problem, q)
    gmsh.model.setCurrent(problem.name)

    nsteps = size(q, 2)
    σ = []
    numElem = Int[]

    for ipg in 1:length(problem.material)
        phName, E, ν, ρ = problem.material[ipg]
        dim = 0
        if problem.dim == 3 && problem.type == "Solid"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            dim = 3
            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            dim = 2
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            dim = 2
            rowsOfB = 3
            b = 1
        else
            error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
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
                s0 = zeros(rowsOfB * numNodes)
                nodeCoord = zeros(numNodes * 3)
                for k in 1:dim, j = 1:numNodes
                    nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
                end
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
                ∇h = reshape(dfun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numNodes)
                ∂h = zeros(3, numNodes * numNodes)
                B = zeros(rowsOfB * numNodes, dim * numNodes)
                nn2 = zeros(Int, dim * numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numNodes
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numNodes, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*3-0, l*2-0] = B[k*3-2, l*2-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*3-0, l*2-1] = B[k*3-1, l*2-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*6-5, l*3-2] = B[k*6-2, l*3-1] = B[k*6-0, l*3-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*6-4, l*3-1] = B[k*6-2, l*3-2] = B[k*6-1, l*3-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*6-3, l*3-0] = B[k*6-1, l*3-1] = B[k*6-0, l*3-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    s = zeros(9numNodes, nsteps) # tensors have nine elements
                    for k in 1:numNodes
                        if rowsOfB == 6 && dim == 3
                            B1 = B[k*6-5:k*6, 1:3*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q[nn2, kk]
                                s[(k-1)*9+1:k*9, kk] = [s0[1], s0[4], s0[6],
                                    s0[4], s0[2], s0[5],
                                    s0[6], s0[5], s0[3]]
                            end
                        elseif rowsOfB == 3 && dim == 2 && problem.type == "PlaneStress"
                            B1 = B[k*3-2:k*3, 1:2*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q[nn2, kk]
                                s[(k-1)*9+1:k*9, kk] = [s0[1], s0[3], 0,
                                    s0[3], s0[2], 0,
                                    0, 0, 0]
                            end
                        else
                            error("solveStress: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
                        end
                    end
                    push!(σ, s)
                end
            end
        end
    end
    sigma = StressField(σ, numElem, nsteps)
    return sigma
end

"""
    FEM.initialDisplacement!(problem, name, u0; ux=..., uy=..., uz=...)

Changes the displacement values `ux`, `uy` and `uz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Original values are in
displacement vector `u0`.

Return: none

Types:
- `problem`: Problem
- `name`: String 
- `u0`: Vector{Float64}
- `ux`: Float64 
- `uy`: Float64 
- `uz`: Float64 
"""
function initialDisplacement!(problem, name, u0; ux=1im, uy=1im, uz=1im)
    dim = problem.dim
    phg = getTagForPhysicalName(name)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
    if ux != 1im
        for i in 1:length(nodeTags)
            u0[nodeTags[i]*dim-(dim-1)] = ux
        end
    end
    if uy != 1im
        for i in 1:length(nodeTags)
            u0[nodeTags[i]*dim-(dim-2)] = uy
        end
    end
    if dim == 3 && uz != 1im
        for i in 1:length(nodeTags)
            u0[nodeTags[i]*dim] = uz
        end
    end
end

"""
    FEM.initialVelocity!(problem, name, v0; vx=..., vy=..., vz=...)

Changes the velocity values `vx`, `vy` and `vz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Original values are in
velocity vector `v0`.

Return: none

Types:
- `problem`: Problem
- `name`: String 
- `v0`: Vector{Float64}
- `vx`: Float64 
- `vy`: Float64 
- `vz`: Float64 
"""
function initialVelocity!(problem, name, v0; vx=1im, vy=1im, vz=1im)
    initialDisplacement!(problem, name, v0, ux=vx, uy=vy, uz=vz)
end

"""
    FEM.nodalForce!(problem, name, f0; fx=..., fy=..., fz=...)

Changes the force values `fx`, `fy` and `fz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Original values are in
load vector `f0`.

Return: none

Types:
- `problem`: Problem
- `name`: String 
- `f0`: Vector{Float64}
- `fx`: Float64 
- `fy`: Float64 
- `fz`: Float64 
"""
function nodalForce!(problem, name, f0; fx=1im, fy=1im, fz=1im)
    initialDisplacement!(problem, name, f0, ux=fx, uy=fy, uz=fz)
end

"""
    FEM.nodalAcceleration!(problem, name, a0; ax=..., ay=..., az=...)

Changes the acceleration values `ax`, `ay` and `az` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Original values are in
acceleration vector `a0`.

Return: none

Types:
- `problem`: Problem
- `name`: String 
- `a0`: Vector{Float64}
- `ax`: Float64
- `ay`: Float64
- `az`: Float64
"""
function nodalAcceleration!(problem, name, a0; ax=1im, ay=1im, az=1im)
    initialDisplacement!(problem, name, a0, ux=ax, uy=ay, uz=az)
end

"""
    FEM.smallestPeriodTime(K, M)

Solves the smallest period of time for a dynamic problem given by stiffness
matrix `K` and the mass matrix `M`.`

Return: `Δt`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `Δt`: Float64 
"""
function smallestPeriodTime(K, M)
    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LM)

    err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1]) / norm(K * ϕ[:,1])
    if err > 1e-3 # || true
        error("Túl nagy a hiba a legnagyobb sajátérték számításánál: $err")
    end
    Δt = 2π / √(real(abs(ω²[1])))
    return Δt
end

"""
    FEM.CDM(K, M, C, f, u0, v0, T, Δt)

Solves a transient dynamic problem using central difference method (explicit).
`K` is the stiffness Matrix, `M` is the mass matrix, `C` is the damping matrix,
`f` is the load vector, `u0` is the initial displacement, `v0` is the initial
velocity, `T` is the upper bound ot the time intervall (lower bound is zero)
and `Δt` is the time step size. Returns the displacement vectors and velocity
vectors in each time step arranged in the columns of the two matrices `u` and `v`
and a vector `t` of the time instants used.

Return: `u`, `v`, `t`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `C`: SparseMatrix
- `f`: Vector{Float64}
- `u0`: Vector{Float64}
- `v0`: Vector{Float64}
- `T`: Float64
- `Δt`: Float64 
- `u`: Matrix{Float64}
- `v`: Matrix{Float64}
- `t`: Vector{Float64}
"""
function CDM(K, M, C, f, u0, v0, T, Δt)
    invM = spdiagm(1 ./ diag(M))
    nsteps = ceil(Int64, T / Δt)
    dof, dof = size(K)

    u = zeros(dof, nsteps)
    v = zeros(dof, nsteps)
    #p = zeros(nsteps)
    t = zeros(nsteps)
    kene = zeros(nsteps)
    sene = zeros(nsteps)
    diss = zeros(nsteps)

    a0 = M \ (f - K * u0 - C * v0)
    u00 = u0 - v0 * Δt + a0 * Δt^2 / 2

    u[:, 1] = u0
    v[:, 1] = v0
    t[1] = 0
    kene[1] = dot(v0' * M, v0) / 2
    sene[1] = dot(u0' * K, u0) / 2

    for i in 2:nsteps
        u1 = 2.0 * u0 - u00 + Δt * Δt * invM * (f - K * u0) - Δt * invM * (C * (u0 - u00))
        u[:, i] = u1
        v1 = (u1 - u0) / Δt
        v[:, i] = v1
        t[i] = t[i-1] + Δt
        kene[i] = dot(v1' * M, v1) / 2
        sene[i] = dot(u1' * K, u1) / 2
        #diss[i] = dot(v1' * C, v1)
        u00 = u0
        u0 = u1
    end
    return u, v, t
end

"""
    FEM.showDoFResults(problem, q, comp; t=..., name=..., visible=...)

Loads nodal results into a View in gmsh. `q` is the field to show, `comp` is
the component of the field ("uvec", "ux", "uy", "uz", "vvec", "vx", "vy", "vz"),
`t` is a vector of time steps (same number of columns as `q`), `name` is a
title to display and `visible` is a true or false value to toggle on or off the 
initial visibility in gmsh. If `q` has more columns, then a sequence of results
will be shown (eg. as an animation). This function returns the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `q`: Vector{Matrix}
- `comp`: String
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showDoFResults(problem, q, comp; t=[0.0], name="u", visible=false)
    gmsh.model.setCurrent(problem.name)
    dim = problem.dim
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, -1, true)
    non = length(nodeTags)
    uvec = gmsh.view.add(name)
    if size(q, 2) != length(t)
        error("showDoFResults: number of time steps missmatch ($(size(q,2)) <==> $(length(t))).")
    end
    for j in 1:length(t)
        k = 1im
        if comp == "uvec" || comp == "vvec"
            nc = 3
            u = zeros(3 * non)
            for i in 1:length(nodeTags)
                u[3i-2] = q[dim*nodeTags[i]-(dim-1), j]
                u[3i-1] = q[dim*nodeTags[i]-(dim-2), j]
                u[3i-0] = dim == 3 ? q[dim*nodeTags[i]-(dim-3), j] : 0
            end
        else
            nc = 1
            if comp == "ux" || comp == "vx"
                k = 1
            elseif comp == "uy" || comp == "vy"
                k = 2
            elseif comp == "uz" || comp == "vz"
                k = 3
            else
                error("ShowDisplacementResults: component is $comp ????")
            end
            u = zeros(non)
            for i in 1:length(nodeTags)
                u[i] = dim == 2 && k == 3 ? 0 : q[dim*nodeTags[i]-(dim-k)]
            end
        end
        gmsh.view.addHomogeneousModelData(uvec, j-1, problem.name, "NodeData", nodeTags, u, t[j], nc)
    end

    gmsh.view.option.setNumber(uvec, "DisplacementFactor", 0)
    gmsh.view.option.setNumber(uvec, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(uvec, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uvec, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(uvec, "Visible", 0)
    end
    display("$comp..ok")
    return uvec
end

"""
    FEM.showStressResults(problem, S, comp; t=..., name=..., visible=..., smooth=...)

Loads stress results into a View in gmsh. `S` is a stress field to show, `comp` is
the component of the field ("s", "sx", "sy", "sz", "sxy", "syz", "szx"),
`t` is a vector of time steps (same length as the number of stress states),
`name` is a title to display, `visible` is a true or false value to toggle on or
off the initial visibility in gmsh and `smooth` is a true of false value to toggle
smoothing the stress field on or off. If length of `t` is more than one, then a 
sequence of results will be shown (eg. as an animation). This function returns
the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `S`: StressField
- `comp`: String
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStressResults(problem, S, comp; t=[0.0], name="σ", visible=false, smooth=true)
    gmsh.model.setCurrent(problem.name)
    dim = problem.dim
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    if S.nsteps != length(t)
        error("showStressResults: number of time steps missmatch ($(S.nsteps) <==> $(length(t))).")
    end
    SS = gmsh.view.add(name)
    σ = S.sigma
    numElem = S.numElem
    for jj in 1:length(t)

        k = 1im
        if comp == "s"
            σcomp = [σ[i][:,jj] for i in 1:length(S.numElem)]
            nc = 9
        else
            nc = 1
            if comp == "sx"
                k = 8
            elseif comp == "sy"
                k = 4
            elseif comp == "sz"
                k = 0
            elseif comp == "sxy" || comp == "syx"
                k = 7
            elseif comp == "syz" || comp == "szy"
                k = 3
            elseif comp == "szx" || comp == "sxz"
                k = 6
            else
                error("ShowStressResults: component is $comp ????")
            end
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                sx = zeros(div(size(σ[i], 1), 9))
                for j in 1:(div(size(σ[i], 1), 9))
                    sx[j] = σ[i][9j-k, jj]
                end
                push!(σcomp, sx)
            end
        end
        gmsh.view.addModelData(SS, jj-1, problem.name, "ElementNodeData", numElem, σcomp, t[jj], nc)
    end

    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end

    gmsh.view.option.setNumber(SS, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(SS, "TargetError", -1e-4)
    gmsh.view.option.setNumber(SS, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(SS, "Visible", 0)
    end
    display("$comp..ok")
    return SS
end

"""
    FEM.plotOnPath(problem, pathName, field, points; numOfSteps=..., name=..., visible=...)

Load a 2D plot on a path into a View in gmsh. `field` is the number of View in
gmsh from which the data of a field is imported. `pathName` is the name of a
physical group which contains a curve. The curve is devided into equal length
intervals with number of `points` points. The field is shown at this points.
`numOfSteps` is the sequence number of steps. `name` is the title of graph and
`visible` is a true or false value to toggle on or off the initial visibility 
in gmsh. This function returns the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `pathName`: String
- `field`: Integer
- `points`: Integer
- `numOfStep`: Integer
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function plotOnPath(problem, pathName, field, points; numOfStep=0, name="path", visible=false)
    gmsh.model.setCurrent(problem.name)
    dimTags = gmsh.model.getEntitiesForPhysicalName(pathName)
    i = 1
    while dimTags[i][1] != 1
        i += 1
        if i > length(dimTags)
            error("Physical name '$name' with dimension ONE does not exist.")
        end
    end
    path = dimTags[i][2]
    dataType, tags, data, time, numComponents = gmsh.view.getModelData(field, numOfStep)
    bounds = gmsh.model.getParametrizationBounds(1, path)
    bound1 = bounds[1][1]
    bound2 = bounds[2][1]
    step0 = (bound2 - bound1) / (points - 1)
    cv = zeros(4)
    CoordValue = []
    pt0 = gmsh.model.getValue(1, path, [bound1])
    for i in 1:points
        pt1 = gmsh.model.getValue(1, path, [bound1 + (i - 1) * step0])
        cv[1:3] = pt1 - pt0
	val, dis = gmsh.view.probe(field, pt1[1], pt1[2], pt1[3])
        if dis < 1e-5
            if numComponents == 1
                v = val[1]
            elseif numComponents == 3
                v = √(val[1]^2 + val[1]^2 + val[1]^2)
            elseif numComponents == 9
                v = √(0.5 * ((val[1] - val[5])^2 + (val[5] - val[9])^2 + (val[9] - val[1])^2 + 6 * (val[2]^2 + val[3]^2 + val[6]^2)))
            else
                error("Vagy skalás vagy vektor vagy tenzor...")
            end
        else
            v = 0
        end
        cv[4] = v
        append!(CoordValue, cv)
    end
    pathView = gmsh.view.add(name)
    gmsh.view.addListData(pathView, "SP", points, CoordValue)

    gmsh.view.option.setNumber(pathView, "Type", 2)
    gmsh.view.option.setNumber(pathView, "Axes", 1)

    if visible == false
        gmsh.view.option.setNumber(pathView, "Visible", 0)
    end
    return pathView
end

end #module
