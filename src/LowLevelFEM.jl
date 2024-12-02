module LowLevelFEM

using LinearAlgebra, SparseArrays
using Arpack
import gmsh_jll
include(gmsh_jll.gmsh_api)
import .gmsh
export gmsh

"""
    Material(phName, E, ν, ρ, k, c, α)

A structure containing the material constants. 
- E: elastic modulus,
- ν: Poisson's ratio,
- ρ: mass density,
- k: heat conductivity,
- c: specific heat,
- α: heat expansion coefficient
`phName` is the name of the physical group where the given material is used.

Types:
- `phName`: String
- `E`: Float64
- `ν`: Float64
- `ρ`: Float64
- `k`: Float64
- `c`: Float64
- `α`: Float64
"""
struct Material
    phName::String
    E::Float64
    ν::Float64
    ρ::Float64
    k::Float64
    c::Float64
    α::Float64
end


"""
    Problem(materials; thickness=..., type=..., bandwidth=...)

A structure containing the most important data of the problem. 
- Parts of the model with their material constants. More materials can be given. (see `material` function)
- type of the problem: 3D "Solid", "PlaneStrain", "PlaneStress", "AxiSymmetric",
  "PlaneHeatConduction", "HeatConduction", "AxiSymmetricHeatConduction".
  In the case of "AxiSymmetric", the axis of symmetry is the "y" axis, 
  while the geometry must be drawn in the positive "x" half-plane.
- bandwidth optimization using built-in `gmsh` function.
  Possibilities: "RCMK", "Hilbert", "Metis" or "none" (default)
- dimension of the problem, determined from `type`
- material constants: Young's modulus, Poisson's ratio,
  mass density, heat conduction corfficient, specific heat, heat 
  expansion coefficient (in a vector of material structure `materials`)
- thickness of the plate
- number of nodes (non)

Types:
- `materials`: Material
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
    pdim::Int64
    material::Vector{Material}
    thickness::Float64
    non::Int64
    function Problem(mat; thickness=1, type="Solid", bandwidth="none")
        if type == "Solid"
            dim = 3
            pdim = 3
        elseif type == "PlaneStress"
            dim = 2
            pdim = 2
        elseif type == "PlaneStrain"
            dim = 2
            pdim = 2
        elseif type == "AxiSymmetric"
            dim = 2
            pdim = 2
        elseif type == "PlaneHeatConduction"
            dim = 2
            pdim = 1
        elseif type == "HeatConduction"
            dim = 3
            pdim = 1
        elseif type == "AxiSymmetricHeatConduction"
            dim = 2
            pdim = 1
        else
            error("Problem type can be: 'Solid', PlaneStress', 'PlaneStrain', 'AxiSymmetric', 'PlaneHeatConduction', 'HeatConduction' or 'AxiSymmetricHeatConduction'. Now problem type = $type ????")
        end
        if !isa(mat, Vector)
            error("Problem: materials are not arranged in a vector. Put them in [...]")
        end
        name = gmsh.model.getCurrent()
        gmsh.option.setString("General.GraphicsFontEngine", "Cairo")
        gmsh.option.setString("View.Format", "%.6g")

        material = mat
        elemTags = []
        for ipg in 1:length(material)
            phName = material[ipg].phName
            dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
            for idm in 1:length(dimTags)
                dimTag = dimTags[idm]
                edim = dimTag[1]
                etag = dimTag[2]
                if edim != dim
                    error("Problem: dimension of the problem ($dim) is different than the dimension of finite elements ($edim).")
                end
                elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(edim, etag)
                for i in 1:length(elementTags)
                    for j in 1:length(elementTags[i])
                        push!(elemTags, elementTags[i][j])
                    end
                end
            end
        end

        if bandwidth != "RCMK" && bandwidth != "Hilbert" && bandwidth != "Metis" && bandwidth != "none"
            error("Problem: bandwidth can be 'Hilbert', 'Metis', 'RCMK' or 'none'. Now it is '$(bandwidth)'")
        end

        method = bandwidth == "none" ? "RCMK" : bandwidth
        oldTags, newTags = gmsh.model.mesh.computeRenumbering(method, elemTags)
        if bandwidth == "none"
            permOldTags = sortperm(oldTags)
            sortNewTags = 1:length(oldTags)
            newTags[permOldTags] = sortNewTags
        end
        gmsh.model.mesh.renumberNodes(oldTags, newTags)

        non = length(oldTags)
        return new(name, type, dim, pdim, material, thickness, non)
    end
end

"""
    VectorField(sigma, numElem, nsteps)

A structure containing the data of a heat flux field. 
- sigma: vector of ElementNodeData type heat flux data (see gmsh.jl)
- numElem: vector of tags of elements
- nsteps: number of stress fields stored in sigma (for animations).
- type: type of data (eg. stress "s", strain "s", heat flux "q")

Types:
- `sigma`: Vector{Matrix{Float64}}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: String
"""
struct VectorField
    sigma::Vector{Matrix{Float64}}
    numElem::Vector{Int}
    nsteps::Int
    type::String
end

"""
    TensorField(sigma, numElem, nsteps)

A structure containing the data of a stress or strain field. 
- sigma: vector of ElementNodeData type stress data (see gmsh.jl)
- numElem: vector of tags of elements
- nsteps: number of stress fields stored in sigma (for animations).
- type: type of data (eg. stress "s", strain "s", heat flux "q")

Types:
- `sigma`: Vector{Matrix{Float64}}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: String
"""
struct TensorField
    sigma::Vector{Matrix{Float64}}
    numElem::Vector{Int}
    nsteps::Int
    type::String
end

"""
    Modal(f, ϕ)

A structure containing the eigenfrequencies and eigen modes.
- f: eigenfrequencies
- ϕ: eigen modes

Types:
- `f`: Matrix{Float64}
- `ϕ`: Vector{Float64}
"""
struct Modal
    f::Vector{Float64}
    ϕ::Matrix{Float64}
end

"""
    FEM.material(name; E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5)

Returns a structure in which `name` is the name of a physical group, 
`E` is the modulus of elasticity, `ν` Poisson's ratio and `ρ` is
the mass density, `k` is the heat conductivity, `c` is the specific
heat, `α` is the coefficient of heat expansion.

Return: mat

Types:
- `mat`: Material
- `name`: String
- `E`: Float64
- `ν`: Float64
- `ρ`: Float64
- `k`: Float64
- `c`: Float64
- `α`: Float64
"""
function material(name; E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5)
    return Material(name, E, ν, ρ, k, c, α)
end

"""
    FEM.displacementConstraint(name; ux=..., uy=..., uz=...)

Gives the displacement constraints on `name` physical group. At least one `ux`, 
`uy` or `uz` value have to be given (depending on the dimension of the problem).
`ux`, `uy` or `uz` can be a constant value, or a function of `x`, `y` and `z`.
(E.g. `fn(x,y,z)=5*(5-x)); FEM.displacementConstraint("support1", ux=fn)`)

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `ux`: Float64 or Function
- `uy`: Float64 or Function
- `uz`: Float64 or Function
"""
function displacementConstraint(name; ux=1im, uy=1im, uz=1im)
    bc0 = name, ux, uy, uz
    return bc0
end

"""
    FEM.load(name; fx=..., fy=..., fz=...)

Gives the intensity of distributed load on `name` physical group. At least one `fx`, 
`fy` or `fz` value have to be given (depending on the dimension of the problem). `fx`, 
`fy` or `fz` can be a constant value, or a function of `x`, `y` and `z`.
(E.g. `fn(x,y,z)=5*(5-x)); FEM.load("load1", fx=fn)`)

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `fx`: Float64 or Function
- `fy`: Float64 or Function
- `fz`: Float64 or Function
"""
function load(name; fx=0, fy=0, fz=0)
    ld0 = name, fx, fy, fz
    return ld0
end

"""
    FEM.elasticSupport(name; kx=..., ky=..., kz=...)

Gives the distributed stiffness of the elastic support on `name` physical group.
`kx`, `ky` or `kz` can be a constant value, or a function of `x`, `y` and `z`.
(E.g. `fn(x,y,z)=5*(5-x)); FEM.elasticSupport("supp1", kx=fn)`)
Default values are 1.

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `kx`: Float64 or Function
- `ky`: Float64 or Function
- `kz`: Float64 or Function
"""
function elasticSupport(name; kx=0, ky=0, kz=0)
    es0 = name, kx, ky, kz
    return es0
end

"""
    FEM.temperatureConstraint(name; T=...)

Gives the temperature constraints on `name` physical group. 
`T` can be a constant value, or a function of `x`, `y` and `z`.
(E.g. `fn(x,y,z)=5*(5-x)); FEM.temperatureConstraint("surf1", T=fn)`)

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `T`: Float64 or Function
"""
function temperatureConstraint(name; T=1im)
    bc0 = name, T, 1im, 1im
    return bc0
end

"""
    FEM.heatFlux(name; qn=...)

Gives the heat flux normal to the surface of `name` physical group.
`qn` can be a constant value, or a function of `x`, `y` and `z`.
(E.g. `fn(x,y,z)=5*(5-x)); FEM.load("flux1", qn=fn)`)

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `qn`: Float64 or Function
"""
function heatFlux(name; qn=0)
    p1 =0
    p2 =0
    qn0 = -qn
    fl0 = name, qn0, p1, p2
    return fl0
end

"""
    FEM.heatSource(name; h=...)

Gives the body heat source in `name` physical group.
`h` can be a constant value, or a function of `x`, `y` and `z`.
(E.g. `fn(x,y,z)=5*(5-x)); FEM.load("source1", h=fn)`)

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `h`: Float64 or Function
"""
function heatSource(name; h=0)
    p1 =0
    p2 =0
    h0 = -h
    sr0 = name, h0, p1, p2
    return sr0
end

"""
    FEM.heatConvection(name; h=..., Tₐ=...)

Gives the heat convection of the surface given with `name` physical group.
`h` is the heat transfer coefficient of the surrounding media,
`Tₐ` is the ambient temperature.

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `h`: Float64
- `Tₐ`: Float64
"""
function heatConvection(name; h=10., Tₐ=20.)
    p = 2im
    hcv0 = name, h, Tₐ, p
    return hcv0
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
    if problem.type == "AxiSymmetric"
        return stiffnessMatrixAXI(problem, elements=elements)
    else
        return stiffnessMatrixSolid(problem, elements=elements)
    end
end

function stiffnessMatrixSolid(problem; elements=[])
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
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 3 && problem.type == "Solid"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
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
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
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
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += B1' * D * B1 * b * jacDet[k] * intWeights[k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
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
    return K
end

function stiffnessMatrixAXI(problem; elements=[])
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
    ncoord2 = zeros(3 * problem.non)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 2 && problem.type == "AxiSymmetric"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                ν 1-ν ν 0;
                ν ν 1-ν 0;
                0 0 0 (1-2ν)/2]

            rowsOfB = 4
        else
            error("stiffnessMatrixAxiSymmetric: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                r = zeros(numIntPoints)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 4
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-3, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = B[k*rowsOfB-0, l*pdim-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-2, l*pdim-1] = h[l, k] / r[k]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        #r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += 2π * B1' * D * B1 * r[k] * jacDet[k] * intWeights[k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
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
    return K
end

"""
    FEM.heatConductionMatrix(problem)

Solves the heat conduction matrix of the `problem`.

Return: `heatCondMat`

Types:
- `problem`: Problem
- `heatCondMat`: SparseMatrix
"""
function heatConductionMatrix(problem; elements=[])
    if problem.type == "AxiSymmetricHeatConduction"
        return heatCondMatrixAXI(problem, elements=elements)
    else
        return heatCondMatrixSolid(problem, elements=elements)
    end
end

function heatCondMatrixSolid(problem; elements=[])
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
        kk = problem.material[ipg].k
        dim = problem.dim
        pdim = problem.pdim
        b = problem.thickness
        if problem.type == "HeatConduction"
            rowsOfB = 3
        elseif problem.type == "PlaneHeatConduction"
            rowsOfB = 2
        else
            error("heatCondMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
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
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
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
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 2
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-2, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += B1' * B1 * kk * b * jacDet[k] * intWeights[k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
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
    return K
end

function heatCondMatrixAXI(problem; elements=[])
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
    ncoord2 = zeros(3 * problem.non)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        kk = problem.material[ipg].k
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 2 && problem.type == "AxiSymmetricHeatConduction"
            rowsOfB = 2
        else
            error("heatCondMatrixAXI: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                r = zeros(numIntPoints)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 2
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            #B[k*rowsOfB-2, l*pdim-1] = h[l, k] / r[k]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        #r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += 2π * B1' * B1 * kk * r[k] * jacDet[k] * intWeights[k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
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
    ncoord2 = zeros(3 * problem.non)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dim = problem.dim
        pdim = problem.pdim
        ρ = problem.material[ipg].ρ
        if problem.dim == 3 && problem.type == "Solid"
            rowsOfH = 3
            b = 1
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            rowsOfH = 2
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            rowsOfH = 2
            b = 1
        elseif problem.dim == 2 && problem.type == "AxiSymmetric"
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
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                nn2 = zeros(Int, pdim * numNodes)
                H = zeros(rowsOfH * numIntPoints, pdim * numNodes)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:pdim
                        H[k*pdim-(pdim-kk), l*pdim-(pdim-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                M1 = zeros(pdim * numNodes, pdim * numNodes)
                if problem.type != "AxiSymmetric"
                    for j in 1:length(elemTags[i])
                        elem = elemTags[i][j]
                        for k in 1:numNodes
                            nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                        end
                        jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                        M1 .*= 0
                        for k in 1:numIntPoints
                            H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                            M1 += H1' * H1 * jacDet[k] * intWeights[k]
                        end
                        M1 *= ρ * b
                        for k in 1:pdim
                            nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                        end
                        append!(I, nn2[Iidx[:]])
                        append!(J, nn2[Jidx[:]])
                        append!(V, M1[:])
                    end
                elseif problem.type == "AxiSymmetric"
                    for j in 1:length(elemTags[i])
                        elem = elemTags[i][j]
                        for k in 1:numNodes
                            nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                        end
                        jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                        M1 .*= 0
                        for k in 1:numIntPoints
                            r = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                            H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                            M1 += H1' * H1 * jacDet[k] * r * intWeights[k]
                        end
                        M1 *= 2π * ρ * b
                        for k in 1:pdim
                            nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                        end
                        append!(I, nn2[Iidx[:]])
                        append!(J, nn2[Jidx[:]])
                        append!(V, M1[:])
                    end
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.pdim * problem.non
    M = sparse(I, J, V, dof, dof)
    if lumped == true
        M = spdiagm(vec(sum(M, dims=2))) # lumped mass matrix
    end
    dropzeros!(M)
    return M
end

"""
    FEM.heatCapacityMatrix(problem; lumped=...)

Solves the heat capacity matrix of the `problem`. If `lumped` is true, solves lumped heat capacity matrix.

Return: `heatCapMat`

Types:
- `problem`: Problem
- `lumped`: Boolean
- `massMat`: SparseMatrix
"""
function heatCapacityMatrix(problem; elements=[], lumped=false)
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
    ncoord2 = zeros(3 * problem.non)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dim = problem.dim
        pdim = problem.pdim
        c = problem.material[ipg].c
        ρ = problem.material[ipg].ρ
        if problem.dim == 3 && problem.type == "Solid"
            rowsOfH = 3
            b = 1
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            rowsOfH = 2
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            rowsOfH = 2
            b = 1
        elseif problem.dim == 2 && problem.type == "AxiSymmetric"
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
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                nn2 = zeros(Int, pdim * numNodes)
                H = zeros(rowsOfH * numIntPoints, pdim * numNodes)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:pdim
                        H[k*pdim-(pdim-kk), l*pdim-(pdim-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                M1 = zeros(pdim * numNodes, pdim * numNodes)
                if problem.type != "AxiSymmetric"
                    for j in 1:length(elemTags[i])
                        elem = elemTags[i][j]
                        for k in 1:numNodes
                            nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                        end
                        jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                        M1 .*= 0
                        for k in 1:numIntPoints
                            H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                            M1 += H1' * H1 * jacDet[k] * intWeights[k]
                        end
                        M1 *= c * ρ * b
                        for k in 1:pdim
                            nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                        end
                        append!(I, nn2[Iidx[:]])
                        append!(J, nn2[Jidx[:]])
                        append!(V, M1[:])
                    end
                elseif problem.type == "AxiSymmetric"
                    for j in 1:length(elemTags[i])
                        elem = elemTags[i][j]
                        for k in 1:numNodes
                            nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                        end
                        jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                        M1 .*= 0
                        for k in 1:numIntPoints
                            r = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                            H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                            M1 += H1' * H1 * jacDet[k] * r * intWeights[k]
                        end
                        M1 *= 2π * c * ρ * b
                        for k in 1:pdim
                            nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                        end
                        append!(I, nn2[Iidx[:]])
                        append!(J, nn2[Jidx[:]])
                        append!(V, M1[:])
                    end
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.pdim * problem.non
    M = sparse(I, J, V, dof, dof)
    if lumped == true
        M = spdiagm(vec(sum(M, dims=2))) # lumped mass matrix
    end
    dropzeros!(M)
    return M
end

"""
    FEM.dampingMatrix(K, M, ωₘₐₓ; α=0.0, ξ=..., β=...)

Generates the damping matrix for proportional damping case. **C**=α**M**+β**K**
or **C**=α**M**+β₁**K**+β₂**KM⁻¹K**+β₃**KM⁻¹KM⁻¹K**+⋅⋅⋅. The latter corresponds 
to the damping characteristic characterized by a power series consisting of powers
of the natural frequencies with odd exponents. ξᵢ (`ξ` in the argument list) are the values ​​of the 
individual members of the series corresponding to the ωₘₐₓ value. βᵢ (`β` in the argument list) are the 
coefficients of the series. (see [^4]) Either `ξ` or `β` must be specified. `ξ` or `β` are scalars or 
vectors. `K` is the stiffness matrix, `M` is the mass matrix and `ωₘₐₓ` is the 
largest natural frequency.

Return: `dampingMatrix`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `ωₘₐₓ`: Float64
- `α`: Float64
- `ξ`: Float64 of Vector{Float64}
- `β`: Float64 of Vector{Float64}
- `dampingMatrix`: SparseMatrix
"""
function dampingMatrix(K, M, ωₘₐₓ; α=0.0, ξ=0.01, β=[2ξ[i]/(ωₘₐₓ)^(2i-1) for i in 1:length(ξ)])
    dof, dof = size(M)
    dof2, dof2 = size(K)
    if dof != nnz(M)
        error("dampingMatrix: M is not lumped!")
    end
    if dof != dof2
        error("dampingMatrix: sizes of M and K are not match: $dof <--> $dof2!")
    end
    invM = spdiagm(1 ./ diag(M))
    C = spzeros(dof, dof)
    MK = copy(K)
    iMK = invM * K
    C += α * M
    C += β[1] * MK
    for i in 2:length(β)
        MK *= iMK
        C += β[i] * MK
    end
    dropzeros!(C)
    return C
end

"""
    FEM.elasticSupportMatrix(problem, elSupp)

Solves the elastic support matrix of the `problem`. `elSupp` is a vector of elastic
supports defined in function `FEM.elasticSupport`. This matrix must be added
to the stiffness matrix.

Return: `elSuppMat`

Types:
- `problem`: Problem
- `elSupp`: Vector{Tuple{String, Float64, Float64, Float64}}
- `elSuppMat`: SparseMatrix
"""
function elasticSupportMatrix(problem, elSupports)
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

    if !isa(elSupports, Vector)
        error("elasticSupportMatrix: elastic supports are not arranged in a vector. Put them in [...]")
    end
    pdim = problem.pdim
    DIM = problem.dim
    b = problem.thickness
    non = problem.non
    dof = pdim * non
    ncoord2 = zeros(3 * problem.non)
    for n in 1:length(elSupports)
        name, kx, ky, kz = elSupports[n]
        if problem.pdim == 3
            f = [0, 0, 0]
        elseif problem.pdim == 2
            f = [0, 0]
        elseif problem.pdim == 1
            f = [0]
        else
            error("elasticSupportMatrix: dimension of the problem is $(problem.pdim).")
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
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(order+1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elementTags[ii]), numNodes)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                        end
                    end
                end
                C1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    C1 .*= 0
                    for j in 1:numIntPoints
                        x = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 2]
                        y = 0
                        z = 0
                        if isa(kx, Function) || isa(ky, Function) || isa(kz, Function)
                            y = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 1]
                            z = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 0]
                        end
                        f[1] = isa(kx, Function) ? kx(x, y, z) : kx
                        if problem.pdim > 1
                            f[2] = isa(ky, Function) ? ky(x, y, z) : ky
                        end
                        if problem.pdim == 3
                            f[3] = isa(kz, Function) ? kz(x, y, z) : kz
                        end
                        r = x
                        if pdim == 3
                            kk = [f[1] 0 0; 0 f[2] 0; 0 0 f[3]]
                        elseif pdim == 2
                            kk = [f[1] 0; 0 f[2]]
                        elseif pdim == 1
                            kk = [f[1]]
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
                        elseif DIM == 2 && dim == 2 && problem.type != "AxiSymmetric" && problem.type != "AxiSymmetricHeatConduction"
                            Ja = jacDet[j] * b
                        elseif DIM == 2 && dim == 2 && (problem.type == "AxiSymmetric" || problem.type == "AxiSymmetricHeatConduction")
                            Ja = 2π * jacDet[j] * r
                        elseif DIM == 2 && dim == 1 && problem.type != "AxiSymmetric" && problem.type != "AxiSymmetricHeatConduction"
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == "AxiSymmetric" || problem.type == "AxiSymmetricHeatConduction")
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * r
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
                        else
                            error("applyBoundaryConditions: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                        end
                        C1 += H1' * kk * H1 * Ja * intWeights[j]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, C1[:])
                end
            end
        end
    end

    C = sparse(I, J, V, dof, dof)
    dropzeros!(C)
    return C
end

"""
    FEM.heatConvectionMatrix(problem, heatConvection)

Solves the heat convection matrix of the `problem`. `heatConvection` 
is a vector of heat convection boundary condicions defined in function
`FEM.heatConduction`. This matrix must be substracted from the heat
conduction matrix.

Return: `heatConvMat`

Types:
- `problem`: Problem
- `heatConvection`: Vector{Tuple{String, Float64, Float64, Float64}}
- `heatConvMat`: SparseMatrix
"""
function heatConvectionMatrix(problem, heatConvection)
    if !isa(heatConvection, Vector)
        error("heatConvection: heat convections are not arranged in a vector. Put them in [...]")
    end
    return elasticSupportMatrix(problem, heatConvection)
end

"""
    FEM.loadVector(problem, loads)

Solves a load vector of `problem`. `loads` is a tuple of name of physical group 
`name`, coordinates `fx`, `fy` and `fz` of the intensity of distributed force.
It can solve traction or body force depending on the problem.
In case of 2D problems and Point physical group means concentrated force.
In case of 2D problems and Line physical group means surface force.
In case of 2D problems and Surface physical group means body force.
In case of 3D problems and Point physical group means concentrated force.
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
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(order+1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elementTags[ii]), numNodes)
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
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    f1 .*= 0
                    for j in 1:numIntPoints
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
                        r = x
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
                        elseif DIM == 2 && dim == 2 && problem.type != "AxiSymmetric" && problem.type != "AxiSymmetricHeatConduction"
                            Ja = jacDet[j] * b
                        elseif DIM == 2 && dim == 2 && (problem.type == "AxiSymmetric" || problem.type == "AxiSymmetricHeatConduction")
                            Ja = 2π * jacDet[j] * r
                        elseif DIM == 2 && dim == 1 && problem.type != "AxiSymmetric" && problem.type != "AxiSymmetricHeatConduction"
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == "AxiSymmetric" || problem.type == "AxiSymmetricHeatConduction")
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * r
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
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
    FEM.heatFluxVector(problem, heatFlux)

Solves a heat flux or heat source vector of `problem`. `heatFlux` is a tuple of name of physical group 
`name`, heat flux `qn` normal to the surface of the body. The outward direction is positive.
It can solve heat flux depending on the problem.
In case of 2D problems and Point physical group means concentrated heat flux.
In case of 2D problems and Line physical group means surface heat flux.
In case of 2D problems and Surface physical group means body heat source.
In case of 3D problems and Point physical group means concentrated heat flux.
In case of 3D problems and Line physical group means edge heat source.
In case of 3D problems and Surface physical group means surface heat flux.
In case of 3D problems and Volume physical group means body heat source.

Return: `heatFluxVec`

Types:
- `problem`: Problem
- `heatFlux`: Vector{Tuple{String, Float64, Float64, Float64}}
- `heatFluxVec`: Vector
"""
function heatFluxVector(problem, loads)
    if !isa(loads, Vector)
        error("heatFluxVector: heat fluxes are not arranged in a vector. Put them in [...]")
    end
    return loadVector(problem, loads)
end

"""
    FEM.heatSourceVector(problem, heatSource)

Solves a heat flux or heat source vector of `problem`. `heatSource` is a tuple of name of physical group 
`name`, heat flux `qn` normal to the surface of the body. The outward direction is positive.
It can solve heat flux depending on the problem.
In case of 2D problems and Point physical group means concentrated heat flux.
In case of 2D problems and Line physical group means surface heat flux.
In case of 2D problems and Surface physical group means body heat source.
In case of 3D problems and Point physical group means concentrated heat flux.
In case of 3D problems and Line physical group means edge heat source.
In case of 3D problems and Surface physical group means surface heat flux.
In case of 3D problems and Volume physical group means body heat source.
Same as the `heatFluxVector` function.

Return: `heatSourceVec`

Types:
- `problem`: Problem
- `heatSource`: Vector{Tuple{String, Float64, Float64, Float64}}
- `heatSourceVec`: Vector
"""
function heatSourceVector(problem, loads)
    if !isa(loads, Vector)
        error("heatSource: heat sources are not arranged in a vector. Put them in [...]")
    end
    return loadVector(problem, loads)
end

"""
    FEM.heatConvectionVector(problem, heatConvection)

Solves a heat convection vector of `problem`. `heatConvection` is a tuple of name of physical group 
`name`, coordinates `fx`, `fy` and `fz` of the intensity of distributed force.

Return: `heatConvVec`

Types:
- `problem`: Problem
- `heatConvection`: Vector{Tuple{String, Float64, Float64, Float64}}
- `heatConvVec`: Vector
"""
function heatConvectionVector(problem, heatConvection)
    if !isa(heatConvection, Vector)
        error("heatConvectionVector: heat convections are not arranged in a vector. Put them in [...]")
    end
    return loadVector(problem, heatConvection)
end

"""
    FEM.thermalLoadVector(problem, T; T₀=...)

Solves the thermal load vector from a temperature field `To` for problem `problem`.
`T₀` is the initial temperature field.

Return: `thermLoadVec`

Types:
- `problem`: Problem
- `T`: Vector{Float64}
- `T₀`: Vector{Float64}
- `thermLoadVec`: Vector{Float64}
"""
function thermalLoadVector(problem, T; T₀=1im)
    if problem.type == "AxiSymmetric"
        return thermalLoadVectorAXI(problem, T, T₀=T₀)
    else
        return thermalLoadVectorSolid(problem, T, T₀=T₀)
    end
end

function thermalLoadVectorSolid(problem, T; T₀=1im)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    dim = problem.dim
    pdim = problem.pdim
    dof = problem.non * pdim
    fT = zeros(dof)
    if T₀ == 1im
        T₀ = zeros(problem.non)
    end
    if size(T) != size(T₀) || size(T,1) != problem.non
        error("thermalLoadVectorSolid: size of T [$(size(T))] != size of T₀ [$(size(T₀))], non=$(problem.non)")
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        α = problem.material[ipg].α
        if problem.dim == 3 && problem.type == "Solid"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = α
            E0 = [1,1,1,0,0,0]
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness * α
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            rowsOfB = 3
            b = α
            E0 = [1,1,0]
        else
            error("thermalLoadVectorSolid: dimension is $(problem.dim), problem type is $(problem.type).")
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
                pdimT = 1
                H = zeros(pdimT * numIntPoints, pdimT * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdimT
                            H[j*pdimT-(pdimT-l), k*pdimT-(pdimT-l)] = h[k, j]
                        end
                    end
                end
                invJac = zeros(3, 3numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn1 = zeros(Int, pdimT * numNodes)
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
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:pdimT
                        nn1[k:pdimT:pdimT*numNodes] = pdimT * nnet[j, 1:numNodes] .- (pdimT - k)
                    end
                    f1 .*= 0
                    for k in 1:numIntPoints
                        H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        f1 += B1' * D * E0 * H1 * (T[nn1] - T₀[nn1]) * b * jacDet[k] * intWeights[k]
                    end
                    fT[nn2] += f1
                end
            end
        end
    end
    return fT
end

function thermalLoadVectorAXI(problem, T; T₀=1im)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    dim = problem.dim
    pdim = problem.pdim
    dof = problem.non * pdim
    fT = zeros(dof)
    if T₀ == 1im
        T₀ = zeros(problem.non)
    end
    if size(T) != size(T₀) || size(T,1) != problem.non
        error("thermalLoadVectorSolid: size of T [$(size(T))] != size of T₀ [$(size(T₀))], non=$(problem.non)")
    end
    ncoord2 = zeros(3 * problem.non)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        α = problem.material[ipg].α
        if problem.dim == 2 && problem.type == "AxiSymmetric"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                ν 1-ν ν 0;
                ν ν 1-ν 0;
                0 0 0 (1-2ν)/2]

            rowsOfB = 4
            b = α
            E0 = [1,1,1,0]
        else
            error("stiffnessMatrixAxiSymmetric: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                pdimT = 1
                H = zeros(pdimT * numIntPoints, pdimT * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdimT
                            H[j*pdimT-(pdimT-l), k*pdimT-(pdimT-l)] = h[k, j]
                        end
                    end
                end
                invJac = zeros(3, 3numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn1 = zeros(Int, pdimT * numNodes)
                r = zeros(numIntPoints)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 4
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-3, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = B[k*rowsOfB-0, l*pdim-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-2, l*pdim-1] = h[l, k] / r[k]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    f1 .*= 0
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:pdimT
                        nn1[k:pdimT:pdimT*numNodes] = pdimT * nnet[j, 1:numNodes] .- (pdimT - k)
                    end
                    for k in 1:numIntPoints
                        H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        f1 += 2π * B1' * D * E0  * H1 * (T[nn1] - T₀[nn1]) * b * r[k] * jacDet[k] * intWeights[k]
                    end
                    fT[nn2] += f1
                end
            end
        end
    end
    return fT
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
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
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
    if !isa(supports, Vector)
        error("applyBoundaryConditions: supports are not arranged in a vector. Put them in [...]")
    end
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
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    dof, dof = size(stiffMat)
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
                f0 = stiffMat[:, nodeTagsX] * uux
            else
                f0 = stiffMat[:, nodeTagsX] * ux
                f0 = sum(f0, dims=2)
            end
            loadVec .-= f0
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim - 2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
                f0 = stiffMat[:, nodeTagsY] * uuy
            else
                f0 = stiffMat[:, nodeTagsY] * uy
                f0 = sum(f0, dims=2)
            end
            loadVec .-= f0
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
                f0 = stiffMat[:, nodeTagsZ] * uuz
            else
                f0 = stiffMat[:, nodeTagsZ] * uz
                f0 = sum(f0, dims=2)
            end
            loadVec .-= f0
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
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                if isa(ux, Function)
                    loadVec[j] = uux[jj]
                else
                    loadVec[j] = ux
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
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                if isa(uy, Function)
                    loadVec[j] = uuy[jj]
                else
                    loadVec[j] = uy
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
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                if isa(uz, Function)
                    loadVec[j] = uuz[jj]
                else
                    loadVec[j] = uz
                end
            end
        end
    end

    dropzeros!(stiffMat)
    dropzeros!(massMat)
    dropzeros!(dampMat)
end

"""
    FEM.applyElasticSupport!(problem, stiffMat, elastSupp)

Applies elastic support boundary conditions `elastSupp` on a stiffness matrix
`stiffMat`. Mesh details are in `problem`. `elastSupp` is a tuple of `name`
of physical group and prescribed `kx`, `ky` and `kz` stiffnesses.

Return: none

Types:
- `problem`: Problem
- `stiffMat`: SparseMatrix 
- `elastSupp`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyElasticSupport!(problem, stiffMat, elastSupp)
    if !isa(elastSupp, Vector)
        error("applyElasticSupport!: elastic supports are not arranged in a vector. Put them in [...]")
    end
    C0 = elasticSupportMatrix(problem, elastSupp)
    stiffMat .+= C0
end

"""
    FEM.applyHeatConvection!(problem, heatCondMat, heatFluxVec, heatConv)

Applies heat convectiom boundary conditions `heatConv` on a heat conduction matrix
`heatCondMat` and heat flux vector `heatFluxVec`. Mesh details are in `problem`. `heatCond`
is a tuple of `name` of physical group and prescribed heat transfer coefficient `h`
and ambient temperature `Tₐ`.

Return: none

Types:
- `problem`: Problem
- `heatCondMat`: SparseMatrix 
- `heatFluxVec`: Vector 
- `heatConv`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyHeatConvection!(problem, heatCondMat, heatFluxVec, heatConv)
    if !isa(heatConv, Vector)
        error("applyHeatConvection!: heat convections are not arranged in a vector. Put them in [...]")
    end
    hf0 = heatConvectionVector(problem, heatConv)
    C0 = heatConvectionMatrix(problem, heatConv)
    heatCondMat .+= C0
    heatFluxVec .+= hf0
    return nothing
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
    FEM.solveTemperature(K, q)

Solves the equation K*T=q for the temperature vector `T`. `K` is the heat conduction matrix,
`q` is the heat flux vector.

Return: `T`

Types:
- `K`: SparseMatrix 
- `q`: Vector{Float64} 
- `T`: Vector{Float64}
"""
function solveTemperature(K, q)
    return K \ q
end

"""
    FEM.solveStrain(problem, q; DoFResults=false)

Solves the strain field `E` from displacement vector `q`. Strain field is given
per elements, so it usually contains jumps at the boundary of elements. Details
of mesh is available in `problem`. If `DoFResults` is true, `E` is a matrix with
nodal results. In this case `showDoFResults` can be used to show the results.

Return: `E`

Types:
- `problem`: Problem
- `q`: Vector{Float64}
- `E`: TensorField or Matrix{Float64}
"""
function solveStrain(problem, q; DoFResults=false)
    gmsh.model.setCurrent(problem.name)

    type = "e"
    nsteps = size(q, 2)
    ε = []
    numElem = Int[]
    ncoord2 = zeros(3 * problem.non)
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    if DoFResults == true
        E1 = zeros(non * 9, nsteps)
        pcs = zeros(Int64, non * dim)
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        ν = problem.material[ipg].ν
        dim = 0
        if problem.dim == 3 && problem.type == "Solid"
            dim = 3
            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            dim = 2
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            dim = 2
            rowsOfB = 3
            b = 1
        elseif problem.dim == 2 && problem.type == "AxiSymmetric"
            dim = 2
            rowsOfB = 4
            b = 1
        else
            error("solveStrain: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
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
                B = zeros(rowsOfB * numNodes, dim * numNodes)
                nn2 = zeros(Int, dim * numNodes)
                r = zeros(numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numNodes
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
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
                    elseif dim == 2 && rowsOfB == 4
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*4-2, l*2-1] = r[k] < 1e-10 ? 0 : h[l, k] / r[k]
                        end
                    else
                        error("solveStrain: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    e = zeros(9numNodes, nsteps) # tensors have nine elements
                    for k in 1:numNodes
                        if rowsOfB == 6 && dim == 3 && problem.type == "Solid"
                            B1 = B[k*6-5:k*6, 1:3*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * q[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[4], e0[6],
                                        e0[4], e0[2], e0[5],
                                        e0[6], e0[5], e0[3]]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[4], e0[6], e0[4], e0[2], e0[5], e0[6], e0[5], e0[3]]
                                end
                            end
                        elseif rowsOfB == 3 && dim == 2 && problem.type == "PlaneStress"
                            B1 = B[k*3-2:k*3, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * q[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[3], 0,
                                        e0[3], e0[2], 0,
                                        0, 0, ν/(ν-1)*(e0[1]+e0[2])]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[3], 0, e0[3], e0[2], 0, 0, 0, 0]
                                end
                            end
                        elseif rowsOfB == 3 && dim == 2 && problem.type == "PlaneStrain"
                            B1 = B[k*3-2:k*3, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * q[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[3], 0,
                                        e0[3], e0[2], 0,
                                        0, 0, 0]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[3], 0, e0[3], e0[2], 0, 0, 0, ν*(e0[1]+e0[2])]
                                end
                            end
                        elseif rowsOfB == 4 && dim == 2 && problem.type == "AxiSymmetric"
                            B1 = B[k*4-3:k*4, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * q[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[4], 0,
                                        e0[4], e0[3], 0,
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
        for k in 1:9
            for l in 1:non
                E1[k + 9 * l - 9, :] ./= pcs[l]
            end
        end
    end
    if DoFResults == true
        return E1
    else
        epsilon = TensorField(ε, numElem, nsteps, type)
        return epsilon
    end
end

"""
    FEM.solveStress(problem, q; DoFResults=false)

Solves the stress field `S` from displacement vector `q`. Stress field is given
per elements, so it usually contains jumps at the boundary of elements. Details
of mesh is available in `problem`. If `DoFResults` is true, `S` is a matrix with
nodal results. In this case `showDoFResults` can be used to show the results.

Return: `S`

Types:
- `problem`: Problem
- `q`: Vector{Float64}
- `S`: TensorField or Matrix{Float64}
"""
function solveStress(problem, q; T=1im, T₀=1im, DoFResults=false)
    gmsh.model.setCurrent(problem.name)

    type = "s"
    nsteps = size(q, 2)
    σ = []
    numElem = Int[]
    ncoord2 = zeros(3 * problem.non)
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    if T₀ == 1im
        T₀ = zeros(problem.non)
    end
    if DoFResults == true
        S1 = zeros(non * 9, nsteps)
        pcs = zeros(Int64, non * dim)
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        α = problem.material[ipg].α
        ακ = α * E / ν / (1 - 2ν)
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
            E0 = [1,1,1,0,0,0]
        elseif problem.dim == 2 && problem.type == "PlaneStress"
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            dim = 2
            rowsOfB = 3
            b = problem.thickness
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == "PlaneStrain"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            dim = 2
            rowsOfB = 3
            b = 1
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == "AxiSymmetric"
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                ν 1-ν ν 0;
                ν ν 1-ν 0;
                0 0 0 (1-2ν)/2]
            dim = 2
            rowsOfB = 4
            b = 1
        else
            error("solveStress: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                #s0 = zeros(rowsOfB * numNodes)
                nodeCoord = zeros(numNodes * 3)
                for k in 1:dim, j = 1:numNodes
                    nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
                end
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
                ∇h = reshape(dfun, :, numNodes)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
                h = reshape(fun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                pdimT = 1
                H = zeros(pdimT * numNodes, pdimT * numNodes)
                for j in 1:numNodes
                    for k in 1:numNodes
                        for l in 1:pdimT
                            H[j*pdimT-(pdimT-l), k*pdimT-(pdimT-l)] = h[k, j]
                        end
                    end
                end
                invJac = zeros(3, 3numNodes)
                ∂h = zeros(3, numNodes * numNodes)
                B = zeros(rowsOfB * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn1 = zeros(Int, pdimT * numNodes)
                r = zeros(numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numNodes
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numNodes, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
                    end
                    B .*= 0
                    if pdim == 2 && rowsOfB == 3
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif pdim == 3 && rowsOfB == 6
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif pdim == 2 && rowsOfB == 4
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*4-2, l*2-1] = r[k] < 1e-10 ? 0 : h[l, k] / r[k]
                        end
                    else
                        error("solveStress: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:pdimT
                        nn1[k:pdimT:pdimT*numNodes] = pdimT * nnet[j, 1:numNodes] .- (pdimT - k)
                    end
                    s = zeros(9numNodes, nsteps) # tensors have nine elements
                    for k in 1:numNodes
                        if rowsOfB == 6 && pdim == 3 && problem.type == "Solid"
                            H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                            B1 = B[k*rowsOfB-5:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q[nn2, kk]
                                if T != 1im
                                    s0 -= D * E0 * H1 * (T[nn1, kk] - T₀[nn1]) * α
                                end
                                if DoFResults == false
                                    s[(k-1)*9+1:k*9, kk] = [s0[1], s0[4], s0[6],
                                        s0[4], s0[2], s0[5],
                                        s0[6], s0[5], s0[3]]
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[4], s0[6], s0[4], s0[2], s0[5], s0[6], s0[5], s0[3]]
                                end
                            end
                        elseif rowsOfB == 3 && pdim == 2 && problem.type == "PlaneStress"
                            H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                            B1 = B[k*rowsOfB-2:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q[nn2, kk]
                                if T != 1im
                                    s0 -= D * E0 * H1 * (T[nn1, kk] - T₀[nn1]) * α
                                end
                                if DoFResults == false
                                    s[(k-1)*9+1:k*9, kk] = [s0[1], s0[3], 0,
                                        s0[3], s0[2], 0,
                                        0, 0, 0]
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[3], 0, s0[3], s0[2], 0, 0, 0, 0]
                                end
                            end
                        elseif rowsOfB == 3 && dim == 2 && problem.type == "PlaneStrain"
                            H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                            B1 = B[k*rowsOfB-2:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q[nn2, kk]
                                if T != 1im
                                    s0 -= D * E0 * H1 * (T[nn1, kk] - T₀[nn1]) * α
                                end
                                if DoFResults == false
                                    s[(k-1)*9+1:k*9, kk] = [s0[1], s0[3], 0,
                                        s0[3], s0[2], 0,
                                        0, 0, ν*(s0[1]+s0[2])]
                                    # PlaneStain: σz ≠ 0
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[3], 0, s0[3], s0[2], 0, 0, 0, ν*(s0[1]+s0[2])]
                                end
                            end
                        elseif rowsOfB == 4 && dim == 2 && problem.type == "AxiSymmetric"
                            B1 = B[k*4-3:k*4, 1:2*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q[nn2, kk]
                                if DoFResults == false
                                    s[(k-1)*9+1:k*9, kk] = [s0[1], s0[4], 0,
                                        s0[4], s0[3], 0,
                                        0, 0, s0[2]]
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[4], 0, s0[4], s0[3], 0, 0, 0, s0[2]]
                                end
                            end
                        else
                            error("solveStress: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
                        end
                    end
                    if DoFResults == true
                        pcs[nnet[j,1:numNodes]] .+= 1
                    end
                    if DoFResults == false
                        push!(σ, s)
                    end
                end
            end
        end
    end
    if DoFResults == true
        for k in 1:9
            for l in 1:non
                S1[k + 9 * l - 9, :] ./= pcs[l]
            end
        end
    end
    if DoFResults == true
        return S1
    else
        sigma = TensorField(σ, numElem, nsteps, type)
        return sigma
    end
end

"""
    FEM.solveHeatFlux(problem, T; DoFResults=false)

Solves the heat flux field `q` from temperature vector `T`. heat flux is given
per elements, so it usually contains jumps at the boundary of elements. Details
of mesh is available in `problem`. If `DoFResults` is true, `q` is a matrix with
nodal results. In this case `showDoFResults` can be used to show the results.

Return: `q`

Types:
- `problem`: Problem
- `T`: Vector{Float64}
- `q`: VectorField or Matrix{Float}
"""
function solveHeatFlux(problem, T; DoFResults=false)
    gmsh.model.setCurrent(problem.name)

    type = "q"
    nsteps = size(T, 2)
    σ = []
    numElem = Int[]
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    ncoord2 = zeros(3 * non)
    if DoFResults == true
        q1 = zeros(non * dim, nsteps)
        pcs = zeros(Int64, non * dim)
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        kT = -problem.material[ipg].k
        dim = 0
        if problem.dim == 3 && problem.type == "HeatConduction"
            dim = 3
            rowsOfB = 3
            b = 1
        elseif problem.dim == 2 && problem.type == "PlaneHeatConduction"
            dim = 2
            rowsOfB = 2
            b = 1
        elseif problem.dim == 2 && problem.type == "AxiSymmetricHeatConduction"
            dim = 2
            rowsOfB = 2
            b = 1
        else
            error("solveStress: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
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
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
                h = reshape(fun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                pdimT = 1
                H = zeros(pdimT * numNodes, pdimT * numNodes)
                for j in 1:numNodes
                    for k in 1:numNodes
                        for l in 1:pdimT
                            H[j*pdimT-(pdimT-l), k*pdimT-(pdimT-l)] = h[k, j]
                        end
                    end
                end
                invJac = zeros(3, 3numNodes)
                ∂h = zeros(3, numNodes * numNodes)
                B = zeros(rowsOfB * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                r = zeros(numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numNodes
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numNodes, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
                    end
                    B .*= 0
                    if pdim == 1 && rowsOfB == 2
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif pdim == 1 && rowsOfB == 3
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*rowsOfB-2, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("solveStress: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    s = zeros(3numNodes, nsteps) # vectors have three elements
                    for k in 1:numNodes
                        if rowsOfB == 3 && pdim == 1
                            B1 = B[k*rowsOfB-2:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = B1 * T[nn2, kk] * kT
                                s[(k-1)*3+1:k*3, kk] = [s0[1], s0[2], s0[3]]
                                if DoFResults == true
                                    #q1[dim*nnet[j, k]-2, kk] += s0[1]
                                    #q1[dim*nnet[j, k]-1, kk] += s0[2]
                                    #q1[dim*nnet[j, k]-0, kk] += s0[3]
                                    q1[dim*nnet[j, k]-2:dim*nnet[j,k], kk] .+= s0
                                end
                            end
                        elseif rowsOfB == 2 && pdim == 1
                            B1 = B[k*rowsOfB-1:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = B1 * T[nn2, kk] * kT
                                s[(k-1)*3+1:k*3, kk] = [s0[1], s0[2], 0]
                                if DoFResults == true
                                    #q1[dim*nnet[j, k]-1, kk] += s0[1]
                                    #q1[dim*nnet[j, k]-0, kk] += s0[2]
                                    q1[dim*nnet[j, k]-1:dim*nnet[j,k], kk] .+= s0
                                end
                            end
                        else
                            error("solveStress: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
                        end
                    end
                    if DoFResults == true
                        pcs[nnet[j,1:numNodes]] .+= 1
                    end
                    if DoFResults == false
                        push!(σ, s)
                    end
                end
            end
        end
    end
    if DoFResults == true
        for k in 1:dim
            for l in 1:non
                q1[k + dim * l - dim, :] ./= pcs[l]
            end
        end
    end
    if DoFResults == true
        return q1
    else
        sigma = VectorField(σ, numElem, nsteps, type)
        return sigma
    end
end

"""
    FEM.solveEigenModes(K, M; n=6, fₘᵢₙ=1.01)

Solves the eigen frequencies and mode shapes of a problem given by stiffness
matrix `K` and the mass matrix `M`. `n` is the number of eigenfrequencies to solve,
and solves the eigenfrequencies greater than `fₘᵢₙ`. Returns the struct of eigenfrequencies
and eigen modes.

Return: `modes`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `n`: Int64
- `fₘᵢₙ`: Float64
- `modes`: Modal 
"""
function solveEigenModes(K, M; n=6, fₘᵢₙ=1.01)
    ωₘᵢₙ² = (2π * fₘᵢₙ)^2
    ω², ϕ = Arpack.eigs(K, M, nev=n, which=:LR, sigma=ωₘᵢₙ², maxiter=10000)
    #if real(ω²[1]) > 0.999 && real(ω²[1]) < 1.001
    #    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    #end
    #err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1]) / norm(K * ϕ[:,1])
    #if err > 1e-3 # || true
    #    error("The error in the calculation of the smallest eigenvalue is too large: $err")
    #end
    f = sqrt.(abs.(real(ω²))) / 2π
    ϕ1 = real(ϕ)
    return Modal(f, ϕ1)
end

#=
"""
    FEM.resultant(problem, phName, field, component)

Solves a load vector of `problem`. `loads` is a tuple of name of physical group 
`name`, coordinates `fx`, `fy` and `fz` of the intensity of distributed force.
It can solve traction or body force depending on the problem.
In case of 2D problems and Point physical group means concentrated force.
In case of 2D problems and Line physical group means surface force.
In case of 2D problems and Surface physical group means body force.
In case of 3D problems and Point physical group means concentrated force.
In case of 3D problems and Line physical group means edge force.
In case of 3D problems and Surface physical group means surface force.
In case of 3D problems and Volume physical group means body force.

Return: `loadVec`

Types:
- `problem`: Problem
- `loads`: Vector{Tuple{String, Float64, Float64, Float64}}
- `loadVec`: Vector
"""
=#
function resultant(problem, field, phName)
    gmsh.model.setCurrent(problem.name)
    pdim = problem.pdim
    DIM = problem.dim
    b = problem.thickness
    non = problem.non
    dof = pdim * non
    fp = zeros(dof)
    ncoord2 = zeros(3 * problem.non)
    sum0 = 0
    dataType, tags, data, time, numComponents = gmsh.view.getModelData(field, 0)
    if numComponents != 1
        error("resultant: number of component of the field must be one.")
    end
    for n in 1:1#length(loads)
        #name, fx, fy, fz = loads[n]
        #if pdim == 3
        #    f = [.0, .0, .0]
        #elseif pdim == 2
        #    f = [.0, .0]
        #elseif pdim == 1
        #    f = [.0]
        #else
        #    error("resultant: DOF per nodes is $(pdim).")
        #end
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
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
                #nnoe = reshape(elemNodeTags[ii], numNodes, :)'
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(order+1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elementTags[ii]), numNodes)
                #H = zeros(pdim * numIntPoints, pdim * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            #H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                        end
                    end
                end
                #f1 = zeros(pdim * numNodes)
                #nn2 = zeros(Int, pdim * numNodes)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    s1 = 0
                    for j in 1:numIntPoints
                        x = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 2]
                        y = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 1]
                        z = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 0]
                        f, d = gmsh.view.probe(field, x, y, z)
                        r = x
                        #H1 = H[j*pdim-(pdim-1):j*pdim, 1:pdim*numNodes] # H1[...] .= H[...] ????
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
                        elseif DIM == 2 && dim == 2 && problem.type != "AxiSymmetric" && problem.type != "AxiSymmetricHeatConduction"
                            Ja = jacDet[j] * b
                        elseif DIM == 2 && dim == 2 && (problem.type == "AxiSymmetric" || problem.type == "AxiSymmetricHeatConduction")
                            Ja = 2π * jacDet[j] * r
                        elseif DIM == 2 && dim == 1 && problem.type != "AxiSymmetric" && problem.type != "AxiSymmetricHeatConduction"
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == "AxiSymmetric" || problem.type == "AxiSymmetricHeatConduction")
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * r
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
                        else
                            error("applyBoundaryConditions: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                        end
                        #f1 += H1' * f * Ja * intWeights[j]
                        s1 += f[1] * Ja * intWeights[j]
                    end
                    for k in 1:pdim
                        #nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                    end
                    sum0 += s1
                end
            end
        end
    end
    return sum0
end

"""
    FEM.resultant(field, phName; dim=2)

Solves the resultant of `field`` (load vector or heat flux vector) on `phName` physical group.
Dimension of the problem (DoF per node) van be given with `dim`. Return the resultant(s) in `tuple`.
Number of the members in `tuple` depends on `dim`.

Return: `res``

or

Return: `resx`, `resy`

or

Return: `resx`, `resy`, `resz`

Types:
- `field`: Vector{Float64}
- `phName`: String 
- `dim`: Int64
- `res`: Float64 
- `resx`: Float64 
- `resy`: Float64 
- `resz`: Float64 
"""
function resultant(field, phName; dim=2, axiSymmetric=false)
    if !isa(field, Vector) && !isa(field, Matrix)
        error("resultant: field must be a Vector or Matrix")
    end
    ph1 = getTagForPhysicalName(phName)
    nodes0, coords = gmsh.model.mesh.getNodesForPhysicalGroup(-1,ph1)
    nodes = Vector{Int64}(nodes0)
    s = [0.0, 0.0, 0.0]
    for i in 1:dim
        for j in 1:length(nodes)
            b = axiSymmetric == true ? 2π * coords[3j-2] : 1
            s[i] += field[dim * nodes[j] - (dim - i)] * b
        end
    end
    if dim == 1
        return s[1]
    elseif dim == 2
        return s[1], s[2]
    elseif dim == 3
        return s[1], s[2], s[3]
    end
end

#=
function elementToNode(problem, field)
    if field.type == "q"
        comp = 3
    elseif field.type == "s" || field.type == "e"
        comp = 9
    end
    dof = problem.non * comp
    vec = zeros(dof, field.nsteps)
    for i in 1:length(field.numElem)
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

    for i in 1:nsteps
        for j in 1:length(field.numElem)
            for k in 1:comp
end
=#

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
    FEM.largestPeriodTime(K, M)

Solves the largest period of time for a dynamic problem given by stiffness
matrix `K` and the mass matrix `M`.`

Return: `Δt`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `Δt`: Float64 
"""
function largestPeriodTime(K, M)
    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LR, sigma=0.01, maxiter=10000)
    if real(ω²[1]) > 0.999 && real(ω²[1]) < 1.001
        ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    end
    err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1]) / norm(K * ϕ[:,1])
    if err > 1e-3 # || true
        error("The error in the calculation of the smallest eigenvalue is too large: $err")
    end
    Δt = 2π / √(real(abs(ω²[1])))
    return Δt
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
        error("The error in the calculation of the largest eigenvalue is too large: $err")
    end
    Δt = 2π / √(real(abs(ω²[1])))
    return Δt
end

"""
    FEM.CDM(K, M, C, f, u0, v0, T, Δt)

Solves a transient dynamic problem using central difference method (CDM) (explicit).
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
    FEM.HHT(K, M, f, u0, v0, T, Δt; α=..., δ=..., γ=..., β=...)

Solves a transient dynamic problem using HHT-α method[^1] (implicit).
`K` is the stiffness Matrix, `M` is the mass matrix, `f` is the load vector, 
`u0` is the initial displacement, `v0` is the initial velocity, `T` is the 
upper bound ot the time intervall (lower bound is zero) and `Δt` is the time 
step size. Returns the displacement vectors and velocity vectors in each time 
step arranged in the columns of the two matrices `u` and `v` and a vector `t` 
of the time instants used. For the meaning of `α`, `β` and `γ` see [1]. If
`δ` is given, γ=0.5+δ and β=0.25⋅(0.5+γ)².

[^1]: Hilber, Hans M., Thomas JR Hughes, and Robert L. Taylor. "Improved 
    numerical dissipation for time integration algorithms in structural 
    dynamics." Earthquake Engineering & Structural Dynamics 5.3 (1977): 283-292.

Return: `u`, `v`, `t`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `f`: Vector{Float64}
- `u0`: Vector{Float64}
- `v0`: Vector{Float64}
- `T`: Float64
- `Δt`: Float64 
- `α`: Float64
- `β`: Float64
- `γ`: Float64
- `δ`: Float64
- `u`: Matrix{Float64}
- `v`: Matrix{Float64}
- `t`: Vector{Float64}
"""
function HHT(K, M, f, u0, v0, T, Δt; α=0.0, δ=0.0, γ=0.5 + δ, β=0.25 * (0.5 + γ)^2)
    nsteps = ceil(Int64, T / Δt)
    dof, dof = size(K)

    u = zeros(dof, nsteps)
    v = zeros(dof, nsteps)
    t = zeros(nsteps)
    kene = zeros(nsteps)
    sene = zeros(nsteps)
    diss = zeros(nsteps)

    dt = Δt
    dtdt = dt * dt

    c0 = 1.0 / (β * dtdt)
    c1 = γ / (β * dt)
    c2 = 1.0 / (β * dt)
    c3 = 0.5 / β - 1.0
    c4 = γ / β - 1.0
    c5 = dt / 2.0 * (γ / β - 2.0)
    c6 = dt * (1.0 - γ)
    c7 = γ * dt

    a0 = M \ (f - K * u0)

    u[:, 1] = u0
    v[:, 1] = v0
    t[1] = 0
    kene[1] = dot(v0' * M, v0) / 2
    sene[1] = dot(u0' * K, u0) / 2
    
    A = (α + 1) * K + M * c0
    AA = lu(A)

    for i in 2:nsteps
        b = f + M * (u0 * c0 + v0 * c2 + a0 * c3) + α * K * u0
        u1 = AA \ b
        u[:, i] = u1
        a1 = (u1 - u0) * c0 - v0 * c2 - a0 * c3
        v1 = v0 + a0 * c6 + a1 * c7
        v[:, i] = v1
        t[i] = t[i-1] + Δt
        kene[i] = dot(v1' * M, v1) / 2
        sene[i] = dot(u1' * K, u1) / 2
        #diss[i] = dot(v1' * C, v1)
        u0 = u1
        v0 = v1
        a0 = a1
    end
    return u, v, t
end

"""
    FEM.HHTaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, δ=0.0, γ=0.5 + δ, β=0.25 * (0.5 + γ)^2)

Gives some functions (graphs) for accuracy analysis of the HHT-α method[^1]. 
`ωₘᵢₙ` and `ωₘₐₓ` are the square root of smallest and largest eigenvalues of the
**Kϕ**=ω²**Mϕ** eigenvalue problem, `Δt` is the time step size. `type` is one of the
following values:
- "SR": spectral radius
- "ADR": algorithmic damping ratio
- "PE": period error
For details see [^2] and [^3]. 
`n` is the number of points in the graph. For the meaning of `α`, `β` and `γ`
see [1]. If `δ` is given, γ=0.5+δ and β=0.25⋅(0.5+γ)².
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

[^2]: Belytschko, Ted, and Thomas JR, Hughes: "Computational methods for 
    transient analysis", North-Holland, (1983).

[^3]: Serfőző, D., Pere, B.: A method to accurately define arbitrary algorithmic
    damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).
    <https://doi.org/10.1007/s00419-023-02454-9>

Return: `xy`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `ωₘᵢₙ`: Float64
- `ωₘₐₓ`: Float64
- `Δt`: Float64 
- `n`: Int64
- `α`: Float64
- `β`: Float64
- `γ`: Float64
- `δ`: Float64
- `xy`: Tuple{Vector{Float64},Vector{Float64}}
"""
function HHTaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, δ=0.0, γ=0.5 + δ, β=0.25 * (0.5 + γ)^2)
    x = zeros(n)
    y = similar(x)
    invT = range(ωₘᵢₙ/2π, length=n, stop=ωₘₐₓ/2π)
    for i ∈ 1:n
        ω = 2π * invT[i]
        A1 = [1 0 -Δt^2*β
            0 1 -Δt*γ
            (1+α)*ω^2 0 1]
        A2 = [1 Δt Δt^2*(0.5-β)
            0 1 Δt*(1-γ)
            α*ω^2 0 0]

        A = A1 \ A2

        eig = eigen(A)
        ρ, idx = findmax(abs, eig.values)
        λ = eig.values[idx]
        σ = real(λ)
        ε = imag(λ)
        if type == "SR"
            x[i] = log(invT[i] * Δt)
            y[i] = ρ
        elseif type == "ADR"
            x[i] = invT[i] * Δt
            Ω = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = -log(ρ) / 2Ω
            #y[i] = -log(ρ) / atan(ε, σ)
        elseif type == "PE"
            x[i] = invT[i] * Δt
            Ω = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = 1 - Ω/(2π*Δt*invT[i])
        else
            str1 = "HHTaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\nSR: spectral radius\n"
            str5 = "ADR: algorithmic damping ratio\n"
            str6 = "PE: period error\n"
            str7 = "\nFor details see Serfőző, D., Pere, B.: A method to accurately define arbitrary\n"
            str8 = "algorithmic damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).\n"
            str9 = "https://doi.org/10.1007/s00419-023-02454-9\n"
            error(str1*str2*str3*str5*str6*str7*str8*str9)
        end
    end
    return x, y
end

"""
    FEM.CDMaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=..., ξ=..., β=..., show_β=..., show_ξ=...)

Gives some functions (graphs) for accuracy analysis of the CDM method. 
`ωₘᵢₙ` and `ωₘₐₓ` are the square root of smallest and largest eigenvalues of the
**Kϕ**=ω²**Mϕ** eigenvalue problem, `Δt` is the time step size. `type` is one of the
following values:
- "SR": spectral radius
- "PDR": physical damping ratio
- "ADR": algorithmic damping ratio
- "PE": period error
For details see [^3]. 
`n` is the number of points in the graph. The damping matrix is assembled in the 
following ways: **C**=α**M**+β**K** or **C**=α**M**+β₁**K**+β₂**KM⁻¹K**+β₃**KM⁻¹KM⁻¹K**+⋅⋅⋅. 
The latter corresponds to the damping characteristic characterized by a power series consisting of powers
of the natural frequencies with odd exponents. ξᵢ (`ξ` in the argument list) are the values ​​of the 
individual members of the series corresponding to the ωₘₐₓ value. βᵢ (`β` in the argument list) are the 
coefficients of the series. (see [^4]) Either `ξ` or `β` must be specified. `ξ` or `β` are scalars or 
vectors. If `show_β` or `show_ξ` is `true`, the corresponding `β` or `ξ` values will be 
sent to the output.
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

[^4]: Serfőző, D., Pere, B.: An effective reduction method with Caughey damping for 
    spurious oscillations in dynamic problems, preprint, <https://doi.org/10.21203/rs.3.rs-3930320/v1>

Return: `xy`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `ωₘᵢₙ`: Float64
- `ωₘₐₓ`: Float64
- `Δt`: Float64 
- `n`: Int64
- `α`: Float64
- `β`: Float64 of Vector{Float64}
- `ξ`: Float64 of Vector{Float64}
- `show_β`: Boolean
- `show_ξ`: Boolean
- `xy`: Tuple{Vector{Float64},Vector{Float64}}
"""
function CDMaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, ξ=0.01, β=[2ξ[i]/(ωₘₐₓ)^(2i-1) for i in 1:length(ξ)], show_β=false, show_ξ=false)
    if show_β == true
        println("β = $β")
    end
    if show_ξ == true
        ξ = [β[i] / 2 * ωₘₐₓ^((2i-1)) for i in 1:length(β)]
        println("ξ = $ξ")
    end
    #Tₘᵢₙ /= √(1-ξₘₐₓ^2)
    x = zeros(n)
    y = zeros(n)
    ω = range(ωₘᵢₙ, length=n, stop=ωₘₐₓ)
    for i ∈ 1:n
        ξ = α/ω[i]
        for j in 1:length(β)
            ξ += β[j] / 2 * ω[i]^(2j-1)
        end
        Ω = Δt * ω[i]
        A = [2-2ξ*Ω-Ω^2 2ξ*Ω-1
            1 0]

        eig = eigen(A)
        ρ, idx = findmax(abs, eig.values)
        λ = eig.values[idx]
        σ = real(λ)
        ε = imag(λ)
        if type == "SR"
            x[i] = log((ω[i] / 2π) * Δt)
            y[i] = ρ
        elseif type == "ADR"
            x[i] = (ω[i] / 2π) * Δt
            Ω0 = √((log(ρ))^2 + (atan(ε,σ))^2 / 4)
            y[i] = -log(ρ) / 2Ω0
        elseif type == "PDR"
            x[i] = (ω[i] / 2π) * Δt
            for j in 1:length(β)
                y[i] += β[j] / 2 * (2π * x[i] / Δt) ^ (2j-1)
            end
        elseif type == "PE"
            x[i] = (ω[i] / 2π) * Δt
            Ω0 = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = 1 - Ω0/(Δt*ω[i])
        else
            str1 = "CDMaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\nSR: spectral radius\n"
            str4 = "PDR: physical damping ratio\n"
            str5 = "ADR: algorithmic damping ratio\n"
            str6 = "PE: period error\n"
            str7 = "\nFor details see Serfőző, D., Pere, B.: A method to accurately define arbitrary\n"
            str8 = "algorithmic damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).\n"
            str9 = "https://doi.org/10.1007/s00419-023-02454-9\n"
            error(str1*str2*str3*str4*str5*str6*str7*str8*str9)
        end
    end
    return x, y
end


"""
    FEM.showDoFResults(problem, q, comp; t=..., name=..., visible=...)

Loads nodal results into a View in gmsh. `q` is the field to show, `comp` is
the component of the field ("uvec", "ux", "uy", "uz", "vvec", "vx", "vy", "vz",
"qvec", "qx", "qy", "qz", "T", "p", "qn", "s", "sx", "sy", "sz", "sxy", "syx", "syz",
"szy", "szx", "sxz", "e", "ex", "ey", "ez", "exy", "eyx", "eyz", "ezy", "ezx", "exz", "seqv"),
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
function showDoFResults(problem, q, comp; t=[0.0], name=comp, visible=false, ff = 0)
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    pdim = problem.pdim
    pdim = div(size(q,1), problem.non) 
    nodeTags = []
    ##############################################################################
    if problem.type == "Reynolds"
        phName = problem.geometry.phName
        tag = getTagForPhysicalName(phName)
        nT, coords = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
        append!(nodeTags, nT)
    ##############################################################################
    else #########################################################################
        for ipg in 1:length(problem.material)
            phName = problem.material[ipg].phName
            tag = getTagForPhysicalName(phName)
            nT, coords = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
            append!(nodeTags, nT)
        end
    end #########################################################################

    #nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, -1, true)
    non = length(nodeTags)
    uvec = gmsh.view.add(name)
    if size(q, 2) != length(t)
        error("showDoFResults: number of time steps missmatch ($(size(q,2)) <==> $(length(t))).")
    end
    if comp[1] == 's' || comp[1] == 'e'
        pdim = 9
    end
    for j in 1:length(t)
        k = 1im
        if comp == "uvec" || comp == "vvec" || comp == "qvec"
            nc = 3
            u = zeros(3 * non)
            for i in 1:length(nodeTags)
                u[3i-2] = q[pdim*nodeTags[i]-(pdim-1), j]
                u[3i-1] = pdim > 1 ? q[pdim*nodeTags[i]-(pdim-2), j] : 0
                u[3i-0] = pdim == 3 ? q[pdim*nodeTags[i]-(pdim-3), j] : 0
            end
        elseif comp == "s" || comp == "e"
            nc = 9
            u = zeros(9 * non)
            for i in 1:length(nodeTags)
                u[9i-8] = q[pdim*nodeTags[i]-(pdim-1), j]
                u[9i-7] = q[pdim*nodeTags[i]-(pdim-2), j]
                u[9i-6] = q[pdim*nodeTags[i]-(pdim-3), j]
                u[9i-5] = q[pdim*nodeTags[i]-(pdim-4), j]
                u[9i-4] = q[pdim*nodeTags[i]-(pdim-5), j]
                u[9i-3] = q[pdim*nodeTags[i]-(pdim-6), j]
                u[9i-2] = q[pdim*nodeTags[i]-(pdim-7), j]
                u[9i-1] = q[pdim*nodeTags[i]-(pdim-8), j]
                u[9i-0] = q[pdim*nodeTags[i]-(pdim-9), j]
            end
        else
            nc = 1
            if comp == "ux" || comp == "vx" || comp == "p" || comp == "T" || comp == "qx"  || comp == "qn" || comp == "sx" || comp == "ex"
                k = 1
            elseif comp == "uy" || comp == "vy" || comp == "qy" || comp == "syx" || comp == "eyx"
                k = 2
            elseif comp == "uz" || comp == "vz" || comp == "qz" || comp == "szx" || comp == "ezx"
                k = 3
            elseif comp == "sxy" || comp == "exy"
                k = 4
            elseif comp == "sy" || comp == "ey"
                k = 5
            elseif comp == "szy" || comp == "ezy"
                k = 6
            elseif comp == "sxz" || comp == "exz"
                k = 7
            elseif comp == "syz" || comp == "eyz"
                k = 8
            elseif comp == "sz" || comp == "ez"
                k = 9
            elseif comp == "seqv"
                k = 10
            else
                error("ShowDisplacementResults: component is $comp ????")
            end
            u = zeros(non)
            for i in 1:length(nodeTags)
                null = k <= 9 ? 0 : √(0.5 * ((q[pdim*nodeTags[i]-(pdim-1), j]-q[pdim*nodeTags[i]-(pdim-5), j])^2+(q[pdim*nodeTags[i]-(pdim-5), j]-q[pdim*nodeTags[i]-(pdim-9), j])^2+(q[pdim*nodeTags[i]-(pdim-9), j]-q[pdim*nodeTags[i]-(pdim-1), j])^2 + 6*(q[pdim*nodeTags[i]-(pdim-2), j]^2+q[pdim*nodeTags[i]-(pdim-3), j]^2+q[pdim*nodeTags[i]-(pdim-6), j]^2)))
                u[i] = k > pdim ? null : q[pdim*nodeTags[i]-(pdim-k), j]
            end
        end
        gmsh.view.addHomogeneousModelData(uvec, j-1, problem.name, "NodeData", nodeTags, u, t[j], nc)
    end

    gmsh.view.option.setNumber(uvec, "DisplacementFactor", 0)
    gmsh.view.option.setNumber(uvec, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(uvec, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uvec, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(uvec, "Visible", 0)
    end
    if ff == 0 && length(t) > 1
        gmsh.view.option.setNumber(uvec, "ShowTime", 1)
    elseif ff == 1
        gmsh.view.option.setNumber(uvec, "ShowTime", 6)
    end
    return uvec
end

"""
    FEM.showModalResults(problem, Φ, name=..., visible=...)

Loads modal results into a View in gmsh. `Φ` is a struct of Modal. `name` is a
title to display and `visible` is a true or false value to toggle on or off the 
initial visibility in gmsh. Click on ▷| to change the results. This function 
returns the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `Φ`: Modal
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showModalResults(problem, Φ::Modal; name="modal", visible=false, ff=1)
    return showDoFResults(problem, Φ.ϕ, "p", t=Φ.f, name=name, visible=visible, ff=ff)
end

"""
    FEM.showStrainResults(problem, E, comp; t=..., name=..., visible=..., smooth=...)

Loads strain results into a View in gmsh. `E` is a strain field to show, `comp` is
the component of the field ("e", "ex", "ey", "ez", "exy", "eyz", "ezx"),
`t` is a vector of time steps (same length as the number of stress states),
`name` is a title to display, `visible` is a true or false value to toggle on or
off the initial visibility in gmsh and `smooth` is a true of false value to toggle
smoothing the stress field on or off. If length of `t` is more than one, then a 
sequence of results will be shown (eg. as an animation). This function returns
the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `E`: TensorField
- `comp`: String
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStrainResults(problem, E, comp; t=[0.0], name=comp, visible=false, smooth=true)
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    if E.nsteps != length(t)
        error("showStrainResults: number of time steps missmatch ($(S.nsteps) <==> $(length(t))).")
    end
    EE = gmsh.view.add(name)
    ε = E.sigma
    numElem = E.numElem
    for jj in 1:length(t)

        k = 1im
        if comp == "e"
            εcomp = [ε[i][:,jj] for i in 1:length(E.numElem)]
            nc = 9
        else
            nc = 1
            if comp == "ex"
                k = 8
            elseif comp == "ey"
                k = 4
            elseif comp == "ez"
                k = 0
            elseif comp == "exy" || comp == "eyx"
                k = 7
            elseif comp == "eyz" || comp == "ezy"
                k = 3
            elseif comp == "ezx" || comp == "exz"
                k = 6
            else
                error("ShowStressResults: component is $comp ????")
            end
            εcomp = []
            sizehint!(εcomp, length(numElem))
            for i in 1:length(E.numElem)
                ex = zeros(div(size(ε[i], 1), 9))
                for j in 1:(div(size(ε[i], 1), 9))
                    ex[j] = ε[i][9j-k, jj]
                end
                push!(εcomp, ex)
            end
        end
        gmsh.view.addModelData(EE, jj-1, problem.name, "ElementNodeData", numElem, εcomp, t[jj], nc)
    end

    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end

    gmsh.view.option.setNumber(EE, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(EE, "TargetError", -1e-4)
    gmsh.view.option.setNumber(EE, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(EE, "Visible", 0)
    end
    if length(t) > 1
        gmsh.view.option.setNumber(EE, "ShowTime", 1)
    end
    #display("$comp..ok")
    return EE
end

"""
    FEM.showElementResults(problem, F, comp; t=..., name=..., visible=..., smooth=...)

Same as `ShowStressResults` or `showStrainResults`, depending on the type of `F` data field.

Return: `tag`

Types:
- `problem`: Problem
- `F`: TensorField
- `comp`: String
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showElementResults(problem, F, comp; t=[0.0], name=comp, visible=false, smooth=true)
    if F.type == "e"
        return showStrainResults(problem, F, comp; t=t, name=comp, visible=false, smooth=smooth)
    elseif F.type == "s"
        return showStressResults(problem, F, comp; t=t, name=comp, visible=false, smooth=smooth)
    elseif F.type == "q"
        return showHeatFluxResults(problem, F, comp; t=t, name=comp, visible=false, smooth=smooth)
    else
        error("showElementResults: type is '$type'")
    end
end

"""
    FEM.showStressResults(problem, S, comp; t=..., name=..., visible=..., smooth=...)

Loads stress results into a View in gmsh. `S` is a stress field to show, `comp` is
the component of the field ("s", "sx", "sy", "sz", "sxy", "syz", "szx", "seqv"),
`t` is a vector of time steps (same length as the number of stress states),
`name` is a title to display, `visible` is a true or false value to toggle on or
off the initial visibility in gmsh and `smooth` is a true of false value to toggle
smoothing the stress field on or off. If length of `t` is more than one, then a 
sequence of results will be shown (eg. as an animation). This function returns
the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `S`: TensorField
- `comp`: String
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStressResults(problem, S, comp; t=[0.0], name=comp, visible=false, smooth=true)
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
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
        elseif comp == "seqv"
            nc = 1
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                seqv = zeros(div(size(σ[i], 1), 9))
                for j in 1:(div(size(σ[i], 1), 9))
                    sx = σ[i][9j-8, jj]
                    sy = σ[i][9j-4, jj]
                    sz = σ[i][9j-0, jj]
                    sxy = σ[i][9j-7, jj]
                    syz = σ[i][9j-3, jj]
                    szx = σ[i][9j-6, jj]
                    seqv[j] = √(((sx-sy)^2+(sy-sz)^2+(sz-sx)^2)/2. + 3*(sxy^2+syz^2+szx^2))
                end
                push!(σcomp, seqv)
            end
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

    gmsh.view.option.setNumber(SS, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(SS, "TargetError", -1e-4)
    gmsh.view.option.setNumber(SS, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(SS, "Visible", 0)
    end
    if length(t) > 1
        gmsh.view.option.setNumber(SS, "ShowTime", 1)
    end
    #display("$comp..ok")
    return SS
end

"""
    FEM.showHeatFluxResults(problem, Q, comp; t=..., name=..., visible=..., smooth=...)

Loads heat flux results into a View in gmsh. `Q` is a heat flux field to show, `comp` is
the component of the field ("qvec", "qx", "qy", "qz", "q"),
`t` is a vector of time steps (same length as the number of stress states),
`name` is a title to display, `visible` is a true or false value to toggle on or
off the initial visibility in gmsh and `smooth` is a true of false value to toggle
smoothing the stress field on or off. If length of `t` is more than one, then a 
sequence of results will be shown (eg. as an animation). This function returns
the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `S`: TensorField
- `comp`: String
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showHeatFluxResults(problem, S, comp; t=[0.0], name=comp, visible=false, smooth=true)
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    if S.nsteps != length(t)
        error("showStressResults: number of time steps missmatch ($(S.nsteps) <==> $(length(t))).")
    end
    SS = gmsh.view.add(name)
    σ = S.sigma
    numElem = S.numElem
    for jj in 1:length(t)

        k = 1im
        if comp == "qvec"
            σcomp = [σ[i][:,jj] for i in 1:length(S.numElem)]
            nc = 3
        elseif comp == "q"
            nc = 1
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                seqv = zeros(div(size(σ[i], 1), 3))
                for j in 1:(div(size(σ[i], 1), 3))
                    sx = σ[i][3j-2, jj]
                    sy = σ[i][3j-1, jj]
                    sz = σ[i][3j-0, jj]
                    seqv[j] = √(sx^2+sy^2+sz^2)
                end
                push!(σcomp, seqv)
            end
        else
            nc = 1
            if comp == "qx"
                k = 2
            elseif comp == "qy"
                k = 1
            elseif comp == "qz"
                k = 0
            else
                error("ShowHeatFluxResults: component is $comp ????")
            end
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                ss = zeros(div(size(σ[i], 1), 3))
                for j in 1:(div(size(σ[i], 1), 3))
                    ss[j] = σ[i][3j-k, jj]
                end
                push!(σcomp, ss)
            end
        end
        gmsh.view.addModelData(SS, jj-1, problem.name, "ElementNodeData", numElem, σcomp, t[jj], nc)
    end

    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end

    gmsh.view.option.setNumber(SS, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(SS, "TargetError", -1e-4)
    gmsh.view.option.setNumber(SS, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(SS, "Visible", 0)
    end
    if length(t) > 1
        gmsh.view.option.setNumber(SS, "ShowTime", 1)
    end
    #display("$comp..ok")
    return SS
end

"""
    FEM.plotOnPath(problem, pathName, field; points=100, step=..., plot=..., name=..., visible=...)

Load a 2D plot on a path into a View in gmsh. `field` is the number of View in
gmsh from which the data of a field is imported. `pathName` is the name of a
physical group which contains a curve. The curve is devided into equal length
intervals with number of `points` points. The field is shown at this points.
`step` is the sequence number of displayed step. If no step is given, shows all 
the aviable steps as an animation. If `plot` is true, additional return parameter, a tuple of
vectors is given back, in which `x` is a vector of values in horizontal axis, `y` is a vector
of values in vertical axis of a plot (see `Plots` package). `name` is the title of graph and
`visible` is a true or false value to toggle on or off the initial visibility 
in gmsh. This function returns the tag of View.

Return: `tag`

or

Return: `tag`, `xy`

Types:
- `problem`: Problem
- `pathName`: String
- `field`: Integer
- `points`: Integer
- `step`: Integer
- `plot`: Boolean
- `name`: String
- `visible`: Boolean
- `tag`: Integer
- `xy`: Tuples{Vector{Float64},Vector{Float64}}
"""
function plotOnPath(problem, pathName, field; points=100, step=1im, plot=false, name="field [$field] on $pathName", visible=false)
    gmsh.model.setCurrent(problem.name)
    dimTags = gmsh.model.getEntitiesForPhysicalName(pathName)
    if points < 2
        error("plotOnPath: number of points is less than two (points = $points)!")
    end
    CoordValue = []
    pt0 = [0, 0, 0]
    pt1 = pt0
    stepRange = 1:2
    for ii ∈ 1:length(dimTags)
        if dimTags[ii][1] != 1
            error("Physical name '$name' with dimension ONE does not exist.")
        end
        path = dimTags[ii][2]
        dataType, tags, data, time, numComponents = gmsh.view.getModelData(field, 0)
        bounds = gmsh.model.getParametrizationBounds(1, path)
        bound1 = bounds[1][1]
        bound2 = bounds[2][1]
        step0 = (bound2 - bound1) / (points - 1)
        nbstep = Int(gmsh.view.option.getNumber(field, "NbTimeStep"))
        if ii == 1
            pt0 = gmsh.model.getValue(1, path, [bound1])
        end
        if step == 1im
            stepRange = 1:nbstep
        else
            stepRange = step >= nbstep ? nbstep : step + 1
            if step >= nbstep
                @warn("plotOnPath: step is greater than max. number of steps (max. number is chosen)  $step <==> $(nbstep)!")
            end
        end
        cv = zeros(3 + length(stepRange))
        for i in 1:points
            pt1 = gmsh.model.getValue(1, path, [bound1 + (i - 1) * step0])
            cv[1:3] = pt1 - pt0
            for j in 1:length(stepRange)
                v = 0
                val, dis = gmsh.view.probe(field, pt1[1], pt1[2], pt1[3], stepRange[j] - 1)
                if dis < 1e-5
                    if numComponents == 1
                        v = val[1]
                    elseif numComponents == 3
                        v = √(val[1]^2 + val[2]^2 + val[3]^2)
                    elseif numComponents == 9
                        v = √(0.5 * ((val[1] - val[5])^2 + (val[5] - val[9])^2 + (val[9] - val[1])^2 + 6 * (val[2]^2 + val[3]^2 + val[6]^2)))
                    else
                        error("Vagy skalár vagy vektor vagy tenzor...")
                    end
                else
                    v = 0
                end
                cv[3+j] = v
            end
            append!(CoordValue, cv)
        end
    end
    pathView = gmsh.view.add(name)
    gmsh.view.addListData(pathView, "SP", points * length(dimTags), CoordValue)

    gmsh.view.option.setNumber(pathView, "Type", 2)
    gmsh.view.option.setNumber(pathView, "Axes", 1)

    if visible == false
        gmsh.view.option.setNumber(pathView, "Visible", 0)
    end

    if plot == true
        step = step == 1im ? 0 : step
        x = zeros(points * length(dimTags))
        y = zeros(points * length(dimTags))
        y[1] = CoordValue[4+step]
        x0 = 0
        y0 = 0
        z0 = 0
        for i in 2:points * length(dimTags)
            x1 = CoordValue[(length(stepRange)+3)*(i-1)+1]
            y1 = CoordValue[(length(stepRange)+3)*(i-1)+2]
            z1 = CoordValue[(length(stepRange)+3)*(i-1)+3]
            x[i] = x[i-1] + √((x1 - x0)^2 + (y1 - y0)^2 + (z1 - z0)^2)
            y[i] = CoordValue[(length(stepRange)+3)*(i-1)+4+step]
            x0 = x1
            y0 = y1
            z0 = z1
        end
        xy = x, y
        return pathView, xy
    else
        return pathView
    end
end

end #module
