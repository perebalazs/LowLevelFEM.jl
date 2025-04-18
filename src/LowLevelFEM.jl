module LowLevelFEM

using LinearAlgebra, SparseArrays
using Arpack
#using Base.Threads
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
- type of the problem: 3D `:Solid`, `:PlaneStrain`, `:PlaneStress`, `:AxiSymmetric`,
  `:PlaneHeatConduction`, `:HeatConduction`, `:AxiSymmetricHeatConduction`.
  In the case of `:AxiSymmetric`, the axis of symmetry is the "y" axis, 
  while the geometry must be drawn in the positive "x" half-plane.
- bandwidth optimization using built-in `gmsh` function.
  Possibilities: `:RCMK`, `:Hilbert`, `:Metis` or `:none` (default)
- dimension of the problem, determined from `type`
- material constants: Young's modulus, Poisson's ratio,
  mass density, heat conduction corfficient, specific heat, heat 
  expansion coefficient (in a vector of material structure `materials`)
- thickness of the plate
- number of nodes (non)

Types:
- `materials`: Material
- `type`: Symbol
- `bandwidth`: String
- `dim`: Integer
- `thickness`: Float64
- `non`: Integer
"""
struct Problem
    name::String
    type::Symbol
    dim::Int64
    pdim::Int64
    material::Vector{Material}
    thickness::Float64
    non::Int64
    function Problem(mat; thickness=1, type=:Solid, bandwidth=:none)
        if type == :Solid
            dim = 3
            pdim = 3
        elseif type == :PlaneStress
            dim = 2
            pdim = 2
        elseif type == :PlaneStrain
            dim = 2
            pdim = 2
        elseif type == :AxiSymmetric
            dim = 2
            pdim = 2
        elseif type == :PlaneHeatConduction
            dim = 2
            pdim = 1
        elseif type == :HeatConduction
            dim = 3
            pdim = 1
        elseif type == :AxiSymmetricHeatConduction
            dim = 2
            pdim = 1
        else
            error("Problem type can be: `:Solid`, `:PlaneStress`, `:PlaneStrain`, `:AxiSymmetric`, `:PlaneHeatConduction`, `:HeatConduction` or `:AxiSymmetricHeatConduction`. Now problem type = $type ????")
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
		    if length(elementTags[i]) == 0
		        error("Problem: No mesh in model '$name'.")
		    end
                    for j in 1:length(elementTags[i])
                        push!(elemTags, elementTags[i][j])
                    end
                end
            end
        end

        if bandwidth != :RCMK && bandwidth != :Hilbert && bandwidth != :Metis && bandwidth != :none
            error("Problem: bandwidth can be `:Hilbert`, `:Metis`, `:RCMK` or `:none`. Now it is `$(bandwidth)`")
        end

        method = bandwidth == :none ? :RCMK : bandwidth
        oldTags, newTags = gmsh.model.mesh.computeRenumbering(method, elemTags)
        if bandwidth == :none
            permOldTags = sortperm(oldTags)
            sortNewTags = 1:length(oldTags)
            newTags[permOldTags] = sortNewTags
        end
        gmsh.model.mesh.renumberNodes(oldTags, newTags)

        non = length(oldTags)
        return new(name, type, dim, pdim, material, thickness, non)
    end
end

using SparseArrays
struct Transformation
    T::SparseMatrixCSC{Float64}
    non::Int64
    dim::Int64
end

import Base.transpose
function transpose(A::Transformation)
    return Transformation(transpose(A.T), A.non, A.dim)
end

import Base.adjoint
function adjoint(A::Transformation)
    return Transformation(adjoint(A.T), A.non, A.dim)
end

import Base.*
function *(A::Transformation, B)
    n = size(B, 1)
    m = size(B, 2)
    non = A.non
    dim = A.dim
    if dim * non == n
        if issparse(B)
            return dropzeros(A.T * B)
        else
            return A.T * B
        end
    elseif 9non == n
        C = zeros(3non, 3)
        D = zeros(3non, 3)
        E = zeros(n, m)
        T = []
        I = []
        J = []
        V = Float64[]
        T1 = zeros(9)
        I0 = [1, 2, 3, 1, 2, 3, 1, 2, 3]
        J0 = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        if dim == 2
            for i in 1:non
                T1 = [A.T[2i-1, 2i-1], A.T[2i, 2i-1], 0, A.T[2i-1, 2i], A.T[2i, 2i], 0, 0, 0, 1]
                Idx = I0 .+ (3i - 3)
                Jdx = J0 .+ (3i - 3)
                append!(I, Idx)
                append!(J, Jdx)
                append!(V, T1)
            end
            fn(x, y) = y
            T = sparse(I, J, V, 3non, 3non, fn)
            dropzeros!(T)
        else
            T = A.T
        end
        for k in 1:m
            for i in 1:non
                for j = 1:3
                    C[3i-2:3i, j] = B[9i-9+3j-2:9i-9+3j, k]
                end
            end
            D = T * C
            for i in 1:non
                for j = 1:3
                    E[9i-9+3j-2:9i-9+3j, k] = D[3i-2:3i, j]
                end
            end
        end
        return E
    else
        error("*(A::Transformation, B): size missmatch dim * non = $dim * $non ≠ $n.")
    end
end

function *(B, A::Transformation)
    n = size(B, 1)
    m = size(B, 2)
    non = A.non
    dim = A.dim
    if dim * non == n
        if issparse(B)
            return dropzeros(A.T * B)
        else
            return A.T * B
        end
    elseif 9non == n
        C = zeros(3, 3non)
        D = zeros(3, 3non)
        E = zeros(n, m)
        T = []
        I = []
        J = []
        V = Float64[]
        T1 = zeros(9)
        I0 = [1, 2, 3, 1, 2, 3, 1, 2, 3]
        J0 = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        if dim == 2
            for i in 1:non
                T1 = [A.T[2i-1, 2i-1], A.T[2i, 2i-1], 0, A.T[2i-1, 2i], A.T[2i, 2i], 0, 0, 0, 1]
                Idx = I0 .+ (3i - 3)
                Jdx = J0 .+ (3i - 3)
                append!(I, Idx)
                append!(J, Jdx)
                append!(V, T1)
            end
            fn(x, y) = y
            T = sparse(I, J, V, 3non, 3non, fn)
            dropzeros!(T)
        else
            T = A.T
        end
        for k in 1:m
            for i in 1:non
                for j = 1:3
                    C[1:3, 3i-3+j] = B[9i-9+3j-2:9i-9+3j, k]
                end
            end
            D = C * T
            for i in 1:non
                for j = 1:3
                    E[9i-9+3j-2:9i-9+3j, k] = D[1:3, 3i-3+j]
                end
            end
        end
        return E
    else
        error("*(B, A::Transformation): size missmatch dim * non = $dim * $non ≠ $n.")
    end
end

function *(A::Transformation, B::Transformation)
    if A.non != B.non || A.dim != B.dim
        error("*(A::Transformation, B::Transformation): size missmatch non = $(A.non) ≠ $(B.non), dim = $(A.dim) ≠ $(B.dim).")
    end
    return Transformation(dropzeros(A.T * B.T), A.non, A.dim)
end

"""
    VectorField(sigma, numElem, nsteps)

A structure containing the data of a heat flux field. 
- sigma: vector of ElementNodeData type heat flux data (see gmsh.jl)
- numElem: vector of tags of elements
- nsteps: number of stress fields stored in sigma (for animations).
- type: type of data (eg. heat flux `:q`)

Types:
- `sigma`: Vector{Matrix{Float64}}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: Symbol
"""
struct VectorField
    sigma::Vector{Matrix{Float64}}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
end

"""
    TensorField(sigma, numElem, nsteps)

A structure containing the data of a stress or strain field. 
- sigma: vector of ElementNodeData type stress data (see gmsh.jl)
- numElem: vector of tags of elements
- nsteps: number of stress fields stored in sigma (for animations).
- type: type of data (eg. stress `:s` and strain `:e`)

Types:
- `sigma`: Vector{Matrix{Float64}}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: Symbol
"""
struct TensorField
    sigma::Vector{Matrix{Float64}}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
end
"""
    CoordinateSystem(vec1, vec2)

A structure containing the data of a coordinate system.
- `vec1`: direction of the new x axis.
- `vec2`: together with `vec1` determine the xy plane
If the problem is two dimensional, it is enough to give the first two
elements of `vec1`. Elements of `vec1` and `vec2` can be functions. In
3D case the functions have three arguments (x, y, and z coordinates),
otherwise (in 2D case) the number of arguments is two (x and y coordinates).

Types:
- `vec1`: Vector{Float64}
- `vec2`: Vector{Float64}

# Examples

```julia
# 2D case
nx(x, y) = x
ny(x, y) = y
cs = FEM.CoordinateSystem([nx, ny])
Q = FEM.rotateNodes(problem, "body", cs)
q2 = Q' * q1 # where `q1` is in Cartesian, `q2` is in Axisymmetric coordinate system and
             # `q1` is a nodal displacement vector.
S2 = Q' * S1 * Q # where `S1` is a stress field in Cartesian coordinate system while
                 # `S2` is in Axisymmetric coordinate system.

# 3D case
n1x(x, y, z) = x
n1y(x, y, z) = y
n2x(x, y, z) = -y
n2y = n1x
cs = FEM.CoordinateSystem([n1x, n1y, 0], [n2x, n2y, 0])
Q = FEM.rotateNodes(problem, "body", cs)
```
"""
struct CoordinateSystem
    vec1::Vector{Float64}
    vec2::Vector{Float64}
    vec1f::Vector{Function}
    vec2f::Vector{Function}
    i1::Vector{Int64}
    i2::Vector{Int64}
    function CoordinateSystem(vect1)
        f(x) = 0
        i1x = isa(vect1[1], Function) ? 1 : 0
        i1y = isa(vect1[2], Function) ? 1 : 0
        v1x = isa(vect1[1], Function) ? 1 : vect1[1]
        v1y = isa(vect1[2], Function) ? 1 : vect1[2]
        v1fx = isa(vect1[1], Function) ? vect1[1] : f
        v1fy = isa(vect1[2], Function) ? vect1[2] : f
        return new([v1x, v1y, 0], [0, 0, 0], [v1fx, v1fy, f], [f, f, f], [i1x, i1y, 0], [0, 0, 0])
    end
    function CoordinateSystem(vect1, vect2)
        f(x) = 0
        i1x = isa(vect1[1], Function) ? 1 : 0
        i1y = isa(vect1[2], Function) ? 1 : 0
        i1z = isa(vect1[3], Function) ? 1 : 0
        i2x = isa(vect2[1], Function) ? 1 : 0
        i2y = isa(vect2[2], Function) ? 1 : 0
        i2z = isa(vect2[3], Function) ? 1 : 0
        v1x = isa(vect1[1], Function) ? 1 : vect1[1]
        v1y = isa(vect1[2], Function) ? 1 : vect1[2]
        v1z = isa(vect1[3], Function) ? 1 : vect1[3]
        v2x = isa(vect2[1], Function) ? 1 : vect2[1]
        v2y = isa(vect2[2], Function) ? 1 : vect2[2]
        v2z = isa(vect2[3], Function) ? 1 : vect2[3]
        v1fx = isa(vect1[1], Function) ? vect1[1] : f
        v1fy = isa(vect1[2], Function) ? vect1[2] : f
        v1fz = isa(vect1[3], Function) ? vect1[3] : f
        v2fx = isa(vect2[1], Function) ? vect2[1] : f
        v2fy = isa(vect2[2], Function) ? vect2[2] : f
        v2fz = isa(vect2[3], Function) ? vect2[3] : f
        return new([v1x, v1y, v1z], [v2x, v2y, v2z], [v1fx, v1fy, v1fz], [v2fx, v2fy, v2fz], [i1x, i1y, i1z], [i2x, i2y, i2z])
    end
end

"""
    Eigen(f, ϕ)

A structure containing the eigenfrequencies and eigen modes.
- f: eigenfrequencies
- ϕ: eigen modes

Types:
- `f`: Matrix{Float64}
- `ϕ`: Vector{Float64}
"""
struct Eigen
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
`Tₐ` is the ambient temperature. The ambient temperature can be either
a constant or a function of x, y and z.

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `h`: Float64
- `Tₐ`: Float64 or Function
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
    if problem.type == :AxiSymmetric
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
        if problem.dim == 3 && problem.type == :Solid
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
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
                #appendlock = ReentrantLock()
                #Threads.@threads for j in 1:length(elemTags[i])
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
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
                    #Threads.lock(appendlock)
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                    #Threads.unlock(appendlock)
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
        if problem.dim == 2 && problem.type == :AxiSymmetric
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
    FEM.nonLinearStiffnessMatrix(problem, q)

Solves the nonlinear stiffness matrix of the `problem`. `q` is a
displacement field.

Return: `stiffMat`

Types:
- `problem`: Problem
- `q`: Vector{Float64}
- `stiffMat`: SparseMatrix
"""
function nonLinearStiffnessMatrix(problem, q; elements=[])
    if problem.type == :AxiSymmetric
        return nonLinearStiffnessMatrixAXI(problem, elements=elements)
    else
        return nonLinearStiffnessMatrixSolid(problem, q, elements=elements)
    end
end

function nonLinearStiffnessMatrixSolid(problem, q; elements=[])
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
    S1 = zeros(problem.dim, problem.dim)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 3 && problem.type == :Solid
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
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
                ∂h = zeros(dim, numNodes * numIntPoints)
                ∂H = zeros(dim * numIntPoints, numNodes)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                K0 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                #for k in 1:numIntPoints, l in 1:numNodes
                #    for kk in 1:pdim
                #        H[k*pdim-(pdim-kk), l*pdim-(pdim-kk)] = h[(k-1)*numNodes+l]
                #    end
                #end
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
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
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    B .*= 0
                    ∂H .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            ∂H[k*pdim*dim-3, l*pdim-1] = ∂H[k*pdim*dim-1, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            ∂H[k*pdim*dim-2, l*pdim-1] = ∂H[k*pdim*dim-0, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                            ∂H[k*dim-2, l] = ∂h[1, (k-1)*numNodes+l]
                            ∂H[k*dim-1, l] = ∂h[2, (k-1)*numNodes+l]
                            ∂H[k*dim-0, l] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("nonLinearStiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    q1 = q[nn2]
                    for k in 1:numIntPoints
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        ∂H1 = ∂H[k*dim-(dim-1):k*dim, 1:numNodes]
                        σ1 = D * B1 * q1
                        if problem.type == :Solid
                            S1[1,1] = σ1[1]
                            S1[2,2] = σ1[2]
                            S1[3,3] = σ1[3]
                            S1[1,2] = S1[2,1] = σ1[4]
                            S1[2,3] = S1[3,2] = σ1[5]
                            S1[3,1] = S1[1,3] = σ1[6]
                        else
                            error("nonLinearStiffnessMatrix: only 'Solid' is implemented")
                        end
                        K0 = ∂H1' * S1 * ∂H1 * b * jacDet[k] * intWeights[k]
                        for kk in 1:3
                            K1[kk:dim:dim*numNodes, kk:dim:dim*numNodes] += K0
                        end
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

function nonLinearStiffnessMatrixAXI(problem; elements=[])
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
        if problem.dim == 2 && problem.type == :AxiSymmetric
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
    if problem.type == :AxiSymmetricHeatConduction
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
        if problem.type == :HeatConduction
            rowsOfB = 3
        elseif problem.type == :PlaneHeatConduction
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
        if problem.dim == 2 && problem.type == :AxiSymmetricHeatConduction
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
        if problem.dim == 3 && problem.type == :Solid
            rowsOfH = 3
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            rowsOfH = 2
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            rowsOfH = 2
            b = 1
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
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
                if problem.type != :AxiSymmetric
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
                elseif problem.type == :AxiSymmetric
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
        if problem.dim == 3 && problem.type == :HeatConduction
            rowsOfH = 3
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneHeatConduction
            rowsOfH = 2
            b = 1
        elseif problem.dim == 2 && problem.type == :AxiSymmetricHeatConduction
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
                if problem.type != :AxiSymmetricHeatConduction
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
                elseif problem.type == :AxiSymmetricHeatConduction
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
    FEM.latentHeatMatrix(problem, u, v, T0)

Solves the latent heat matrix of the `problem`. With this matrix the generated heat due to deformations
(given with displacement field `u` and velocity field `v`) can be solved. `T0` is the current temperature
field which is given in absolute temperature scale (Kelvin).

Return: `latHeatMat`

Types:
- `problem`: Problem
- `u`: Vector{Float64}
- `v`: Vector{Float64}
- `T0`: Vector{Float64}
- `latHeatMat`: SparseMatrix
"""
function latentHeatMatrix(problem, u, v, T0)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
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
        α = problem.material[ipg].α
        κ = E / 3 / (1 - 2ν)
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 3 && problem.type == :HeatConduction
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneHeatConduction
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness
        else
            error("latentHeatMatrix: dimension is $(problem.dim), problem type is $(problem.type).")
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
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                M1 = zeros(pdim * numNodes, pdim * numNodes)
                ∇H = zeros(1 * numIntPoints, dim * numNodes)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, dim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn1 = zeros(Int, pdim * numNodes)
                nn2 = zeros(Int, dim * numNodes)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:pdim
                        H[k*pdim-(pdim-kk), l*pdim-(pdim-kk)] = h[(k-1)*numNodes+l]
                    end
                end
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
                            B[k*rowsOfB-0, l*dim-0] = B[k*rowsOfB-2, l*dim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*dim-1] = B[k*rowsOfB-1, l*dim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-5, l*dim-2] = B[k*rowsOfB-2, l*dim-1] = B[k*rowsOfB-0, l*dim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*dim-1] = B[k*rowsOfB-2, l*dim-2] = B[k*rowsOfB-1, l*dim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*dim-0] = B[k*rowsOfB-1, l*dim-1] = B[k*rowsOfB-0, l*dim-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    ∇H .*= 0
                    if dim == 2
                        for k in 1:numIntPoints, l in 1:numNodes
                            ∇H[k*1-0, l*dim-1] = ∂h[1, (k-1)*numNodes+l]
                            ∇H[k*1-0, l*dim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            ∇H[k*1-0, l*dim-2] = ∂h[1, (k-1)*numNodes+l]
                            ∇H[k*1-0, l*dim-1] = ∂h[2, (k-1)*numNodes+l]
                            ∇H[k*1-0, l*dim-0] = ∂h[3, (k-1)*numNodes+l]
                        end
                    end
                    nn1 = nnet[j, 1:numNodes]
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    M1 .*= 0
                    K1 .*= 0
                    for k in 1:numIntPoints
                        q1 = u[nn2]
                        dq1 = v[nn2]
                        T01 = T0[nn1]
                        H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                        ∇H1 = ∇H[k, 1:dim*numNodes]'
                        M1 += H1' * H1 * (∇H1 * dq1) * jacDet[k] * intWeights[k]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:dim*numNodes]
                        K1 += H1' * H1 * (q1' * B1' * D * B1 * dq1) / (H1 * T01)[1] * b * jacDet[k] * intWeights[k]
                    end
                    M1 *= κ * α * b
                    KM1 = K1 - M1
                    append!(I, nn1[Iidx[:]])
                    append!(J, nn1[Jidx[:]])
                    append!(V, KM1[:])
                end
            end
        end
    end
    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return K
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
supports defined in function `FEM.elasticSupport`. With the displacementent vector `q` in hand the
reaction force vector `fR` arising from the elastic support can be solved.
(`fR = heatConvMat * q`)

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
                        elseif DIM == 2 && dim == 2 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = jacDet[j] * b
                        elseif DIM == 2 && dim == 2 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * jacDet[j] * r
                        elseif DIM == 2 && dim == 1 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * r
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
                        else
                            error("elasticSupportMatrix: dimension of the problem is $(problem.dim), dimension of load is $dim.")
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
`FEM.heatConduction`. With the heat convection vector (see the `heatConvectionVector` function)
`heatConvVec`, temperature field vector `T` in hand the heat flux vector `qCV` arising from the
heat convection boundary condition can be solved. `qCV = heatConvMat * T - heatConvVec`

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
- In case of 2D problems and Point physical group means concentrated force.
- In case of 2D problems and Line physical group means surface force.
- In case of 2D problems and Surface physical group means body force.
- In case of 3D problems and Point physical group means concentrated force.
- In case of 3D problems and Line physical group means edge force.
- In case of 3D problems and Surface physical group means surface force.
- In case of 3D problems and Volume physical group means body force.

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
                        elseif DIM == 2 && dim == 2 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = jacDet[j] * b
                        elseif DIM == 2 && dim == 2 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * jacDet[j] * r
                        elseif DIM == 2 && dim == 1 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * r
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
                        else
                            error("loadVector: dimension of the problem is $(problem.dim), dimension of load is $dim.")
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
It can solve heat flux (or heat source) depending on the problem.
- In case of 2D problems and Point physical group means concentrated heat flux.
- In case of 2D problems and Line physical group means surface heat flux.
- In case of 2D problems and Surface physical group means body heat source.
- In case of 3D problems and Point physical group means concentrated heat flux.
- In case of 3D problems and Line physical group means edge heat source.
- In case of 3D problems and Surface physical group means surface heat flux.
- In case of 3D problems and Volume physical group means body heat source.

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
It can solve heat flux (or heat source) depending on the problem.
- In case of 2D problems and Point physical group means concentrated heat flux.
- In case of 2D problems and Line physical group means surface heat flux.
- In case of 2D problems and Surface physical group means body heat source.
- In case of 3D problems and Point physical group means concentrated heat flux.
- In case of 3D problems and Line physical group means edge heat source.
- In case of 3D problems and Surface physical group means surface heat flux.
- In case of 3D problems and Volume physical group means body heat source.
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

Solves a heat convection vector of `problem`. `heatConvection` 
is a vector of heat convection boundary condicions defined in function
`FEM.heatConduction`. With the heat convection matrix (see the `heatConvectionMatrix` function)
`heatConvMat`, temperature field vector `T` in hand the heat flux vector `qCV` arising from the
heat convection boundary condition can be solved. `qCV = heatConvMat * T - heatConvVec`

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

Solves the thermal load vector from a temperature field `T` for problem `problem`.
`T₀` is the initial temperature field.

Return: `thermLoadVec`

Types:
- `problem`: Problem
- `T`: Vector{Float64}
- `T₀`: Vector{Float64}
- `thermLoadVec`: Vector{Float64}
"""
function thermalLoadVector(problem, T; T₀=1im)
    if problem.type == :AxiSymmetric
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
        if problem.dim == 3 && problem.type == :Solid
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = α
            E0 = [1,1,1,0,0,0]
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness * α
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == :PlaneStrain
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
        if problem.dim == 2 && problem.type == :AxiSymmetric
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                ν 1-ν ν 0;
                ν ν 1-ν 0;
                0 0 0 (1-2ν)/2]

            rowsOfB = 4
            b = α
            E0 = [1,1,1,0]
        else
            error("thermalLoadVectorAxiSymmetric: dimension is $(problem.dim), problem type is $(problem.type).")
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
    FEM.applyBoundaryConditions!(problem, heatCondMat, heatCapMat, heatFluxVec, supports)

Applies boundary conditions `supports` on a heat conduction matrix
`heatCondMat`, heat capacity matrix `heatCapMat` and heat flux vector `heatFluxVec`. Mesh details are in `problem`. `supports`
is a tuple of `name` of physical group and prescribed temperature `T`.

Return: `stiffMat`, `loadVec`

Types:
- `problem`: Problem
- `stiffMat`: SparseMatrix 
- `loadVec`: Vector 
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions!(problem, heatCondMat, heatCapMat, heatFluxVec, supports; fix=1)
    if !isa(supports, Vector)
        error("applyBoundaryConditions: supports are not arranged in a vector. Put them in [...]")
    end
    dof, dof = size(heatCondMat)
    dampMat = spzeros(dof, dof)
    applyBoundaryConditions!(problem, heatCondMat, heatCapMat, dampMat, heatFluxVec, supports, fix=fix)
    dampMat = []
    return
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
function applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports; fix=1)
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
                stiffMat[j, j] = fix
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = fix
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
                stiffMat[j, j] = fix
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = fix
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
                stiffMat[j, j] = fix
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = fix
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
    FEM.applyBoundaryConditions!(problem, dispVec, supports)

Applies displacement boundary conditions `supports` on a displacement vector
`dispVec`. Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Return: none

Types:
- `problem`: Problem
- `dispVec`: Vector{Float64}
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions!(problem, dispVec, supports)
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
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
                dispVec[nodeTagsX,:] .= uux
            else
                dispVec[nodeTagsX,:] .= ux
            end
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim - 2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
                dispVec[nodeTagsY,:] .= uuy
            else
                dispVec[nodeTagsY,:] .= uy
            end
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
                dispVec[nodeTagsZ,:] .= uuz
            else
                dispVec[nodeTagsZ,:] .= uz
            end
        end
    end
end

"""
    FEM.field(name; f=..., fx=..., fy=..., fz=..., fxy=..., fyz=..., fzx=...)

Gives the value of scalar, vector or tensor field on `name` physical group. At least one `fx`, 
`fy` or `fz` etc. value have to be given (depending on the dimension of the problem). `fx`, 
`fy` or `fz` etc. can be a constant value, or a function of `x`, `y` and `z`.
(E.g. `fn(x,y,z)=5*(5-x)); FEM.field("surf1", fx=fn)`)

Return: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function,...x7}

Types:
- `name`: String
- `f`: Float64 or Function
- `fx`: Float64 or Function
- `fy`: Float64 or Function
- `fz`: Float64 or Function
- `fxy`: Float64 or Function
- `fyz`: Float64 or Function
- `fzx`: Float64 or Function

# Examples

```julia
f1(x, y, z) = y
f2 = FEM.field("face1", f=f1)
qq = FEM.scalarField(problem, [f2])
qqq = FEM.showDoFResults(problem, qq, :scalar)
```
"""
function field(name; f=:no, fx=:no, fy=:no, fz=:no, fxy=:no, fyx=:no, fyz=:no, fzy=:no, fzx=:no, fxz=:no)
    fxy = fyx != :no ? fyx : fxy
    fyz = fzy != :no ? fzy : fyz
    fzx = fxz != :no ? fxz : fzx
    ld0 = name, f, fx, fy, fz, fxy, fyz, fzx
    return ld0
end

"""
    FEM.scalarField(problem, dataField)

Defines a scalar field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`.

Return: Vector{Float64}

Types:
- `problem`: Problem
- `dataField`: Vector{Tuple{String, Float64,...}}

# Examples

```julia
f2 = FEM.field("face1", f=1)
qq = FEM.scalarField(problem, [f2])
qqq = FEM.showDoFResults(problem, qq, :scalar)
```
"""
function scalarField(problem, dataField)
    if !isa(dataField, Vector)
        error("applyBoundaryConditions!: dataField are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    pdim = 1
    non = problem.non
    field = zeros(non)

    for i in 1:length(dataField)
        name, f, fx, fy, fz, fxy, fyz, fzx = dataField[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(f, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if f != :no
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim - 1)
            if isa(f, Function)
                ff = f.(xx, yy, zz)
                field[nodeTagsX,:] .= ff
            else
                field[nodeTagsX,:] .= f
            end
        end
    end
    return field
end

"""
    FEM.vectorField(problem, dataField)

Defines a vector field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`.

Return: Vector{Float64}

Types:
- `problem`: Problem
- `dataField`: Vector{Tuple{String, Float64,...}}

# Examples

```julia
f1(x, y, z) = sin(x)
f2(x, y, z) = 5y
ff1 = FEM.field("face1", fx=f1, fy=f2, fz=0)
ff2 = FEM.field("face2", fx=f2, fy=f1, fz=1)
qq = FEM.vectorField(problem, [ff1, ff2])
qq0 = FEM.showDoFResults(problem, qq, :vector)
```
"""
function vectorField(problem, dataField)
    if !isa(dataField, Vector)
        error("applyBoundaryConditions!: dataField are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    pdim = problem.pdim
    non = problem.non
    field = zeros(non * pdim)

    for i in 1:length(dataField)
        name, f, fx, fy, fz, fxy, fyz, fzx = dataField[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(fx, Function) || isa(fy, Function) || isa(fz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if fx != :no
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim - 1)
            if isa(fx, Function)
                ffx = fx.(xx, yy, zz)
                field[nodeTagsX,:] .= ffx
            else
                field[nodeTagsX,:] .= fx
            end
        end
        if fy != :no
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim - 2)
            if isa(fy, Function)
                ffy = fy.(xx, yy, zz)
                field[nodeTagsY,:] .= ffy
            else
                field[nodeTagsY,:] .= fy
            end
        end
        if fz != :no && pdim == 3
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= pdim
            if isa(fz, Function)
                ffz = fz.(xx, yy, zz)
                field[nodeTagsZ,:] .= ffz
            else
                field[nodeTagsZ,:] .= fz
            end
        end
    end
    return field
end

"""
    FEM.constrainedDoFs(problem, supports)

Returns the serial numbers of constrained degrees of freedom. Support is a vector of boundary conditions given with the function `displacementConstraint`.

Return: `DoFs`

Types:
- `problem`: Problem
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
- `DoFs`: Vector{Int64}
"""
function constrainedDoFs(problem, supports)
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    #non = problem.non
    pdim = problem.pdim
    #dof = non * pdim
    cdofs = []

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags::Vector{Int64}, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        nodeTagsX = []
        nodeTagsY = []
        nodeTagsZ = []
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
        end
        cdofs = cdofs ∪ nodeTagsX ∪ nodeTagsY ∪ nodeTagsZ
    end

    return cdofs
end

"""
    FEM.freeDoFs(problem, supports)

Returns the serial numbers of unconstrained degrees of freedom. Support is a vector of boundary conditions given with the function `displacementConstraint`.

Return: `DoFs`

Types:
- `problem`: Problem
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
- `DoFs`: Vector{Int64}
"""
function freeDoFs(problem, supports)
    cdofs = constrainedDoFs(problem, supports)
    dof = problem.non * problem.pdim
    fdofs = setdiff(1:dof, cdofs)
    return fdofs
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
`heatCondMat` and heat flux vector `heatFluxVec`. Mesh details are in `problem`. `heatConv`
is a tuple of `name` of physical group and prescribed heat transfer coefficient `h`
and ambient temperature `Tₐ`. The ambient temperature can be either a constant or a 
function of x, y and z coordinates.

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
    FEM.solveDisplacement(problem, load, supp)

Solves the displacement vector `q` of `problem` with loads `load` and
supports `supp`.

Return: `q`

Types:
- `problem`: Problem 
- `load`: Vector{Tuple} 
- `supp`: Vector{Tuple}
- `q`: Vector{Float64}
"""
function solveDisplacement(problem, load, supp)
    K = stiffnessMatrix(problem)
    f = loadVector(problem, load)
    applyBoundaryConditions!(problem, K, f, supp)
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
    FEM.solveTemperature(problem, flux, temp)

Solves the temperature vector `T` of `problem` with given heat flux `flux` and
temperature `temp`.

Return: `T`

Types:
- `problem`: Problem 
- `flux`: Vector{Tuple} 
- `temp`: Vector{Tuple}
- `T`: Vector{Float64}
"""
function solveTemperature(problem, flux, temp)
    K = heatConductionMatrix(problem)
    q = heatFluxVector(problem, flux)
    #applyHeatConvection!(problem, K, q, temp)
    applyBoundaryConditions(problem, K, q, temp)
    return K \ q
end

"""
    FEM.solveTemperature(problem, flux, temp, heatconv)

Solves the temperature vector `T` of `problem` with given heat flux `flux`,
temperature `temp` and heat convection `heatconv`.

Return: `T`

Types:
- `problem`: Problem 
- `flux`: Vector{Tuple} 
- `temp`: Vector{Tuple}
- `heatconv`: Vector{Tuple}
- `T`: Vector{Float64}
"""
function solveTemperature(problem, flux, temp, heatconv)
    K = heatConductionMatrix(problem)
    q = heatFluxVector(problem, flux)
    applyHeatConvection!(problem, K, q, heatconv)
    applyBoundaryConditions(problem, K, q, temp)
    return K \ q
end

"""
    FEM.solveStrain(problem, q; DoFResults=false)

Solves the strain field `E` from displacement vector `q`. Strain field is given
per elements, so it usually contains jumps at the boundaries of elements. Details
of mesh is available in `problem`. If `DoFResults` is true, `E` is a matrix with
nodal results. In this case `showDoFResults` can be used to show the results 
(otherwise `showStrainResults` or `showElementResults`).

Return: `E`

Types:
- `problem`: Problem
- `q`: Vector{Float64}
- `E`: TensorField or Matrix{Float64}
"""
function solveStrain(problem, q; DoFResults=false)
    gmsh.model.setCurrent(problem.name)

    type = :e
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
        if problem.dim == 3 && problem.type == :Solid
            dim = 3
            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            dim = 2
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            dim = 2
            rowsOfB = 3
            b = 1
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
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
                        if rowsOfB == 6 && dim == 3 && problem.type == :Solid
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
                        elseif rowsOfB == 3 && dim == 2 && problem.type == :PlaneStress
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
                        elseif rowsOfB == 3 && dim == 2 && problem.type == :PlaneStrain
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
                        elseif rowsOfB == 4 && dim == 2 && problem.type == :AxiSymmetric
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
    FEM.solveStress(problem, q; T=..., T₀=..., DoFResults=false)

Solves the stress field `S` from displacement vector `q`. Stress field is given
per elements, so it usually contains jumps at the boundary of elements. Details
of mesh is available in `problem`. If `DoFResults` is true, `S` is a matrix with
nodal results. In this case `showDoFResults` can be used to show the results 
(otherwise `showStressResults` or `showElementResults`).
If the `T` temperature field (and `T₀` initial temperature field if it differs from zero) is given, the
function solves also the thermal stresses.

Return: `S`

Types:
- `problem`: Problem
- `q`: Vector{Float64}
- `T`: Vector{Float64}
- `T₀`: Vector{Float64}
- `S`: TensorField or Matrix{Float64}
"""
function solveStress(problem, q; T=1im, T₀=1im, DoFResults=false)
    gmsh.model.setCurrent(problem.name)

    type = :s
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
        if problem.dim == 3 && problem.type == :Solid
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
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            dim = 2
            rowsOfB = 3
            b = problem.thickness
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            dim = 2
            rowsOfB = 3
            b = 1
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
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
                        if rowsOfB == 6 && pdim == 3 && problem.type == :Solid
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
                        elseif rowsOfB == 3 && pdim == 2 && problem.type == :PlaneStress
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
                        elseif rowsOfB == 3 && dim == 2 && problem.type == :PlaneStrain
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
                                    # PlaneStrain: σz ≠ 0
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[3], 0, s0[3], s0[2], 0, 0, 0, ν*(s0[1]+s0[2])]
                                end
                            end
                        elseif rowsOfB == 4 && dim == 2 && problem.type == :AxiSymmetric
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
nodal results. In this case `showDoFResults` can be used to show the results
(otherwise `showHeatFluxResults` or `showElementResults`).

Return: `q`

Types:
- `problem`: Problem
- `T`: Vector{Float64}
- `q`: VectorField or Matrix{Float}
"""
function solveHeatFlux(problem, T; DoFResults=false)
    gmsh.model.setCurrent(problem.name)

    type = :q
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
        if problem.dim == 3 && problem.type == :HeatConduction
            dim = 3
            rowsOfB = 3
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneHeatConduction
            dim = 2
            rowsOfB = 2
            b = 1
        elseif problem.dim == 2 && problem.type == :AxiSymmetricHeatConduction
            dim = 2
            rowsOfB = 2
            b = 1
        else
            error("solveHeatFlux: dimension is $(problem.dim), problem type is $(problem.type).")
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
                        error("solveHeatFlux: rows of B is $rowsOfB, dimension of the problem is $dim.")
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
                            error("solveHeatFlux: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
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
    FEM.elementsToNodes(problem, T)

Solves the nodal results `F` from the elemental results `T`.
`T` can be tensor field or vector field.

Return: `F`

Types:
- `problem`: Problem
- `T`: TensorField or VectorField
- `F`: Matrix{Float64}
"""
function elementsToNodes(problem, S)
    gmsh.model.setCurrent(problem.name)

    type = S.type
    nsteps = S.nsteps
    numElem = S.numElem
    σ = S.sigma
    non = problem.non
    if type == :s || type == :e
        epn = 9
    elseif type == :q
        epn = 3
    else
        error("elementsToNodes: type is $type .")
    end
    s = zeros(non * epn, nsteps)
    pcs = zeros(Int64, non)

    for e in 1:length(numElem)
        elementType, nodeTags, dim, tag = gmsh.model.mesh.getElement(numElem[e])
        for i in 1:length(nodeTags)
            s[(nodeTags[i]-1) * epn + 1: nodeTags[i] * epn, :] .+= σ[e][(i-1)*epn+1:i*epn, :]
            pcs[nodeTags[i]] += 1
        end
    end
    for l in 1:non
        s[epn * (l - 1) + 1: epn * l, :] ./= pcs[l]
    end
    return s
end

"""
    FEM.fieldError(problem, F)

Solves the deviation of field results `F` (stresses, strains, heat flux components) at nodes, where the field has jumps.
The result can be displayed with the `showDoFResults` function.

Return: `e`

Types:
- `problem`: Problem
- `F`: TensorField or VectorField
- `e`: Matrix{Float64}
"""
function fieldError(problem, S)
    gmsh.model.setCurrent(problem.name)

    type = S.type
    nsteps = S.nsteps
    numElem = S.numElem
    σ = S.sigma
    non = problem.non
    if type == :s || type == :e
        epn = 9
    elseif type == :q
        epn = 3
    else
        error("fieldError: type is $type .")
    end
    avg = zeros(non * epn, nsteps)
    res = zeros(non * epn, nsteps)
    #m = zeros(nsteps)
    pcs = zeros(Int64, non)

    for e in 1:length(numElem)
        elementType, nodeTags, dim, tag = gmsh.model.mesh.getElement(numElem[e])
        for i in 1:length(nodeTags)
            avg[(nodeTags[i]-1) * epn + 1: nodeTags[i] * epn, :] .+= σ[e][(i-1)*epn+1:i*epn, :]
            pcs[nodeTags[i]] += 1
            #m = [max(avg[j], m[j]) for j in 1:nsteps]
        end
    end
    for l in 1:non
        avg[epn * (l - 1) + 1: epn * l, :] ./= pcs[l]
    end
    for e in 1:length(numElem)
        elementType, nodeTags, dim, tag = gmsh.model.mesh.getElement(numElem[e])
        for i in 1:length(nodeTags)
            res[(nodeTags[i]-1) * epn + 1: nodeTags[i] * epn, :] .+= (avg[(nodeTags[i]-1) * epn + 1: nodeTags[i] * epn, :] .- σ[e][(i-1)*epn+1:i*epn, :]).^2
        end
    end
    for l in 1:non
        res[epn * (l - 1) + 1: epn * l, :] /= pcs[l]
        res[epn * (l - 1) + 1: epn * l, :] .= sqrt.(res[epn * (l - 1) + 1: epn * l, :]) # ./ m
    end
    return res
end

"""
    FEM.solveEigenModes(K, M; n=6, fₘᵢₙ=1.01)

Solves the eigen frequencies and mode shapes of a problem given by stiffness
matrix `K` and the mass matrix `M`. `n` is the number of eigenfrequencies to solve,
and solves the eigenfrequencies greater than `fₘᵢₙ`. Returns the struct of eigenfrequencies
and eigen modes. Results can be presented by `showModalResults` function.

Return: `modes`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `n`: Int64
- `fₘᵢₙ`: Float64
- `modes`: Eigen 
"""
function solveEigenModes(K, M; n=6, fₘᵢₙ=0.01)
    ωₘᵢₙ² = (2π * fₘᵢₙ)^2
    ω², ϕ = Arpack.eigs(K, M, nev=n, which=:LR, sigma=ωₘᵢₙ², maxiter=10000)
    #if real(ω²[1]) > 0.999 && real(ω²[1]) < 1.001
    #    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    #end
    #err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1]) / norm(K * ϕ[:,1])
    #if err > 1e-3 # || true
    #    @warn("The error in the calculation of the smallest eigenvalue is too large: $err")
    #end
    f = sqrt.(abs.(real(ω²))) / 2π
    ϕ1 = real(ϕ)
    return Eigen(f, ϕ1)
end

"""
    FEM.solveBucklingModes(K, Knl; n=6)

Solves the critical force multipliers and buckling mode shapes of a problem given by stiffness
matrix `K` and the nonlinear stiffness matrix `Knl`. `n` is the number of buckling modes to solve.
Returns the struct of critical forces and buckling modes. Results can be presented by `showBucklingResults` function.

Return: `modes`

Types:
- `K`: SparseMatrix
- `Knl`: SparseMatrix
- `n`: Int64
- `modes`: Eigen 
"""
function solveBucklingModes(K, Knl; n=6)
    λ, ϕ = Arpack.eigs(K, -Knl, nev=n, which=:LR, sigma=0, maxiter=1000)
    #if real(ω²[1]) > 0.999 && real(ω²[1]) < 1.001
    #    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    #end
    #err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1]) / norm(K * ϕ[:,1])
    #if err > 1e-3 # || true
    #    error("The error in the calculation of the smallest eigenvalue is too large: $err")
    #end
    f = abs.(real(λ))
    ϕ1 = real(ϕ)
    return Eigen(f, ϕ1)
end

"""
    FEM.solveModalAnalysis(problem; constraints=[]; loads=[], n=6)

Solves the first `n` eigenfrequencies and the corresponding 
mode shapes for the `problem`, when `loads` and 
`constraints` are applied. `loads` and `contraints` are optional. 
Result can be presented by `showModalResults` function. 
`loads` and `constraints` can be defined by `load` and `displacementConstraint` functions,
respectively. If `loads` are given, it solves the eigenfrequencies of a prestressed structure.

Return: `modes`

Types:
- `problem`: Problem
- `loads`: Vector{tuples}
- `constraints`: Vector{tuples}
- `n`: Int64
- `modes`: Eigen
"""
function solveModalAnalysis(problem; constraints=[], loads=[], n=6, fₘᵢₙ=0.1)
    if !isa(loads, Vector)
        error("solveModalAnalysis: loads are not arranged in a vector. Put them in [...]")
    end
    if !isa(constraints, Vector)
        error("solveModalAnalysis: constraints are not arranged in a vector. Put them in [...]")
    end
    dof = problem.pdim * problem.non
    K = stiffnessMatrix(problem)
    M = massMatrix(problem)
    if length(constraints) == 0
        return solveEigenModes(K, M, n=n, fₘᵢₙ=fₘᵢₙ)
    elseif length(loads) == 0
        fdof = freeDoFs(problem, constraints)
        cdof = constrainedDoFs(problem, constraints)
        K = K[fdof, fdof]
        M = M[fdof, fdof]
        #f = zeros(dof)
        #applyBoundaryConditions!(problem, K, M, f, constraints, fix=fₘᵢₙ<1 ? fₘᵢₙ/10 : 0.1)
        mod = solveEigenModes(K, M, n=n, fₘᵢₙ=fₘᵢₙ)
        nn = length(mod.f)
        ϕ1 = zeros(dof, nn)
        #display("ϕ1[fdof,:]: $(size(ϕ1[fdof,:]))")
        #display("mod.ϕ: $(size(mod.ϕ))")
        ϕ1[fdof,:] .= mod.ϕ
        applyBoundaryConditions!(problem, ϕ1, constraints)
        return Eigen(mod.f, ϕ1)
    else
        fdof = freeDoFs(problem, constraints)
        cdof = constrainedDoFs(problem, constraints)
        ϕ1 = zeros(dof, n)
        f = loadVector(problem, loads)
        applyBoundaryConditions!(problem, K, f, constraints)
        q = solveDisplacement(K, f)

        err = 1
        count = 0
        while err > 1e-3 && count < 10
            count += 1
            q0 = copy(q)
            Knl = nonLinearStiffnessMatrix(problem, q)
            applyBoundaryConditions!(problem, Knl, f, constraints)
            q = solveDisplacement(K + Knl, f)
            err = sum(abs, q - q0) / (sum(abs, q0) == 0 ? 1 : sum(abs, q0))
        end
        if count == 10
            @warn("solveModalAnalysis: number of iterations is $count.")
        end
        Knl = nonLinearStiffnessMatrix(problem, q)

        #applyBoundaryConditions!(problem, K, M, Knl, f, constraints, fix=fₘᵢₙ<1 ? fₘᵢₙ/10 : 0.1)
        mod = solveEigenModes((K + Knl)[fdof,fdof], M[fdof,fdof], n=n, fₘᵢₙ=fₘᵢₙ)
        ϕ1[fdof,:] .= mod.ϕ
        applyBoundaryConditions!(problem, ϕ1, constraints)
        return Eigen(mod.f, ϕ1)
    end
end

"""
    FEM.solveBuckling(problem, loads, constraints; n=6)

Solves the multipliers for the first `n` critical forces and the corresponding 
buckling shapes for the instability of the `problem`, when `loads` and 
`constraints` are applied. Result can be presented by `showBucklingResults`
function. `loads` and `constraints` can be defined by `load` and `displacementConstraint` functions,
respectively.

Return: `buckling`

Types:
- `problem`: Problem
- `loads`: Vector{tuples}
- `constraints`: Vector{tuples}
- `n`: Int64
- `buckling`: Eigen 
"""
function solveBuckling(problem, loads, constraints; n=6)
    f = loadVector(problem, loads)
    K = stiffnessMatrix(problem)
    applyBoundaryConditions!(problem, K, f, constraints)
    q = solveDisplacement(K, f)

    err = 1
    count = 0
    while err > 1e-3 && count < 10
        count += 1
        q0 = copy(q)
        Knl = nonLinearStiffnessMatrix(problem, q)
        applyBoundaryConditions!(problem, Knl, f, constraints)
        q = solveDisplacement(K + Knl, f)
        err = sum(abs, q - q0) / (sum(abs, q0) == 0 ? 1 : sum(abs, q0))
    end
    if count == 10
        @warn("solveBuckling: number of iterations is $count.")
    end
    Knl = nonLinearStiffnessMatrix(problem, q)

    applyBoundaryConditions!(problem, Knl, f, constraints)
    return solveBucklingModes(K, Knl, n=n)
end

"""
    FEM.resultant(problem, field, phName; grad=false, component=:x)

Solves the resultant of `field` on `phName` physical group.
Return the resultant(s) in `tuple`.
Number of the members in `tuple` depends on the dimension of `problem`.
It can solve the resultant of a load vector (sum of the elements of the vector),
if `field` is a vector of floats. If `field` is a view (tag of a view in gmsh), then
the integral of the field is solved. `field` must have only one component.
If `grad` is `true`, then the gradient of the `field` will be evaluated and `component` of the gradient
(`:x`, `:y` or `:z`) will be used to solve the resultant.

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
function resultant(problem, field, phName; grad=false, component=:x, offsetX=0, offsetY=0, offsetZ=0)
    if !isa(field, Vector) && !isa(field, Matrix)
        return resultant2(problem, field, phName, grad, component, offsetX, offsetY, offsetZ)
    end
    if problem.type == :NavierStokes
        dim = problem.pdim + 1
    else
	dim = problem.pdim
    end
    axiSymmetric = false
    if problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction
        axiSymmetric = true
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

function resultant2(problem, field, phName, grad, component, offsetX, offsetY, offsetZ)
    gmsh.model.setCurrent(problem.name)
    pdim = problem.pdim
    DIM = problem.dim
    b = problem.thickness
    non = problem.non
    dof = pdim * non
    fp = zeros(dof)
    ncoord2 = zeros(3 * problem.non)
    sum0 = 0
    if component == :x
        comp0 = 1
    elseif component == :y
        comp0 = 2
    elseif component == :z
        comp0 = 3
    else
        error("resultant: invalid component '$component'")
    end
    dataType, tags, data, time, numComponents = gmsh.view.getModelData(field, 0)
    if numComponents != 1
        error("resultant: number of component of the field must be one.")
    end
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
                    f, d = gmsh.view.probe(field, x+offsetX, y+offsetY, z+offsetZ, -1, -1, grad, -1)
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
                    elseif DIM == 2 && dim == 2 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                        Ja = jacDet[j] * b
                    elseif DIM == 2 && dim == 2 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                        Ja = 2π * jacDet[j] * r
                    elseif DIM == 2 && dim == 1 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                        Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                    elseif DIM == 2 && dim == 1 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                        Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * r
                    elseif DIM == 2 && dim == 0
                        Ja = 1
                    ############ 1D #######################################################
                    elseif DIM == 1 && dim == 1
                        Ja = (Jac[1, 3*j-2]) * b
                    else
                        error("resultant: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                    end
                    #f1 += H1' * f * Ja * intWeights[j]
                    s1 += f[comp0] * Ja * intWeights[j]
                end
                for k in 1:pdim
                    #nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                end
                sum0 += s1
            end
        end
    end
    return sum0
end

"""
    FEM.initialDisplacement(problem, name; ux=..., uy=..., uz=...)

Sets the displacement values `ux`, `uy` and `uz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Returns the initial
displacement vector `u0`.

Return: u0

Types:
- `problem`: Problem
- `name`: String 
- `u0`: Vector{Float64}
- `ux`: Float64 
- `uy`: Float64 
- `uz`: Float64 
"""
function initialDisplacement(problem, name; ux=1im, uy=1im, uz=1im)
    dim = problem.dim
    u0 = zeros(problem.non * problem.dim)
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
    return u0
end

"""
    FEM.initialDisplacement!(problem, name, u0; ux=..., uy=..., uz=...)

Changes the displacement values to `ux`, `uy` and `uz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Original values are in
displacement vector `u0`.

Return: u0

Types:
- `problem`: Problem
- `name`: String 
- `ux`: Float64 
- `uy`: Float64 
- `uz`: Float64 
- `u0`: Vector{Float64}
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
    FEM.initialVelocity(problem, name; vx=..., vy=..., vz=...)

Sets the velocity values `vx`, `vy` and `vz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Returns the initial
velocity vector `v0`.

Return: v0

Types:
- `problem`: Problem
- `name`: String 
- `vx`: Float64 
- `vy`: Float64 
- `vz`: Float64 
- `v0`: Vector{Float64}
"""
function initialVelocity(problem, name; vx=1im, vy=1im, vz=1im)
    return initialDisplacement(problem, name, ux=vx, uy=vy, uz=vz)
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
    FEM.initialTemperature(problem, name; T=...)

Sets the temperature value `T` at nodes belonging to physical group `name`.
Returns the `T0` initial nodal temperature vector.

Return: T0

Types:
- `problem`: Problem
- `name`: String 
- `T`: Float64 
- `T0`: Vector{Float64}
"""
function initialTemperature(problem, name; T=1im)
    dim = problem.pdim
    T0 = zeros(problem.non * problem.pdim)
    phg = getTagForPhysicalName(name)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
    if T != 1im
        for i in 1:length(nodeTags)
            T0[nodeTags[i]*dim-(dim-1)] = T
        end
    end
    return T0
end

"""
    FEM.initialTemperature!(problem, name, T0; T=...)

Changes the tempetature value to `T` at nodes belonging to physical group `name`.
Original values are in temperature vector `T0`.

Return: none

Types:
- `problem`: Problem
- `name`: String 
- `T0`: Vector{Float64}
- `T`: Float64 
"""
function initialTemperature!(problem, name, T0; T=1im)
    dim = problem.pdim
    phg = getTagForPhysicalName(name)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
    if T != 1im
        for i in 1:length(nodeTags)
            T0[nodeTags[i]*dim-(dim-1)] = T
        end
    end
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
matrix `K` and the mass matrix `M`.

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
        @warn("The error in the calculation of the smallest eigenvalue is too large: $err")
    end
    Δt = 2π / √(abs(real(ω²[1])))
    return Δt
end

"""
    FEM.smallestPeriodTime(K, M)

Solves the smallest period of time for a dynamic problem given by stiffness
matrix `K` and the mass matrix `M`.

Return: `Δt`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `Δt`: Float64 
"""
function smallestPeriodTime(K, M)
    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LM, maxiter=10000)
    
    err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1]) / norm(K * ϕ[:,1])
    if err > 1e-3 # || true
        @warn("The error in the calculation of the largest eigenvalue is too large: $err")
    end
    Δt = 2π / √(abs(real(ω²[1])))
    return Δt
end

"""
    FEM.smallestEigenValue(K, M)

Solves the largest eigenvalue for a transient problem given by stiffness (heat conduction)
matrix `K` and the mass (heat capacity) matrix `M` (`C`).

Return: `λₘₐₓ`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `λₘₐₓ`: Float64 
"""
function smallestEigenValue(K, C)
    λ, ϕ = Arpack.eigs(K, C, nev=1, which=:LR, sigma=0.0001, maxiter=10000)
    if real(λ[1]) > 0.999 && real(λ[1]) < 1.001
        λ, ϕ = Arpack.eigs(K, C, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    end
    err = norm(K * ϕ[:,1] - λ[1] * C * ϕ[:,1]) / norm(K * ϕ[:,1])
    if err > 1e-3 # || true
        @warn("The error in the calculation of the largest eigenvalue is too large: $err")
    end
    λₘₐₓ = abs(real(λ[1]))
    return λₘₐₓ
end

"""
    FEM.largestEigenValue(K, M)

Solves the smallest eigenvalue for a transient problem given by stiffness (heat conduction)
matrix `K` and the mass (heat capacity) matrix `M` (`C`).

Return: `λₘᵢₙ`

Types:
- `K`: SparseMatrix
- `M`: SparseMatrix
- `λₘᵢₙ`: Float64 
"""
function largestEigenValue(K, C)
    λ, ϕ = Arpack.eigs(K, C, nev=1, which=:LM)

    err = norm(K * ϕ[:,1] - λ[1] * C * ϕ[:,1]) / norm(K * ϕ[:,1])
    if err > 1e-3 # || true
        @warn("The error in the calculation of the smallest eigenvalue is too large: $err")
    end
    λₘᵢₙ = abs(real(λ[1]))
    return λₘᵢₙ
end

"""
    FEM.CDM(K, M, C, f, u0, v0, T, Δt)

Solves a transient dynamic problem using central difference method (CDM) (explicit).
`K` is the stiffness Matrix, `M` is the mass matrix, `C` is the damping matrix,
`f` is the load vector, `u0` is the initial displacement, `v0` is the initial
velocity, `T` is the upper bound of the time intervall (lower bound is zero)
and `Δt` is the time step size. Returns the displacement vectors and velocity
vectors in each time step arranged in the columns of the two matrices `u` and `v`
and a vector `t` of the time instants used.

The critical (largest allowed) time step is `Δtₘₐₓ = Tₘᵢₙ / π * (√(1 + ξₘₐₓ^2) - ξₘₐₓ)`
where `Tₘᵢₙ` is the time period of the largest eigenfrequency and `ξₘₐₓ` is the largest
modal damping.

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

function CDM(K, M, f, u0, v0, T, Δt)
    dof, dof = size(K)
    C = zeros(dof, dof)
    return CDM(K, M, C, f, u0, v0, T, Δt)
end

"""
    FEM.HHT(K, M, f, u0, v0, T, Δt; α=..., δ=..., γ=..., β=...)

Solves a transient dynamic problem using HHT-α method[^1] (implicit).
`K` is the stiffness Matrix, `M` is the mass matrix, `f` is the load vector, 
`u0` is the initial displacement, `v0` is the initial velocity, `T` is the 
upper bound of the time intervall (lower bound is zero) and `Δt` is the time 
step size. Returns the displacement vectors and velocity vectors in each time 
step arranged in the columns of the two matrices `u` and `v` and a vector `t` 
of the time instants used. For the meaning of `α`, `β` and `γ` see [^1]. If
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
    FEM.FDM(K, C, q, T0, tₘₐₓ, Δt; ϑ=...)

Solves a transient heat conduction problem using Finite Difference Method (FDM).
Introducing a `ϑ` parameter, special cases can be used as the Forward Euler (explicit, ϑ=0),
Backward Euler (implicit, ϑ=1), Crank-Nicolson (ϑ=0.5) and intermediate cases (0<ϑ<1).
(This method is known as ϑ-method. See [^5].)
`K` is the heat conduction matrix, `C` is the heat capacity matrix,
`q` is the heat flux vector, `T0` is the initial temperature, `tₘₐₓ` is the upper 
bound of the time intervall (lower bound is zero)
and `Δt` is the time step size. Returns the nodal temperature vectors in each time 
step arranged in the columns of the matrix `T`
and a vector `t` of the time instants used.

The critical (largest allowed) time step is `Δtₘₐₓ = 2 / ((1-2ϑ)*λₘₐₓ)`
where `λₘₐₓ` is the largest eigenvalue of (**K**+λ**C**)**θ**=**0** 
eigenvalue problem and `ϑ` is the parameter of the ϑ-method. Default value of `ϑ`
is 1/2.

[^5]: Bathe, K. J.: Finite element procedures, Wiley, 1983, <https://doi.org/10.1002/nag.1610070412>

Return: `T`, `t`

Types:
- `K`: SparseMatrix
- `C`: SparseMatrix
- `q`: Vector{Float64}
- `T0`: Vector{Float64}
- `tₘₐₓ`: Float64
- `Δt`: Float64 
- `T`: Matrix{Float64}
- `t`: Vector{Float64}
"""
function FDM(K, C, q, T0, tₘₐₓ, Δt; ϑ=0.5)
    nnz = nonzeros(C)
    dof, dof = size(C)
    nsteps = ceil(Int64, tₘₐₓ / Δt)
    T = zeros(dof, nsteps)
    t = zeros(nsteps)
    T[:, 1] = T0
    t[1] = 0
    if ϑ == 0 && nnz == dof
        invC = spdiagm(1 ./ diag(C))
        for i in 2:nsteps
            T1 = T0 - (1 - ϑ) * Δt * invC * K * T0 + Δt * q
            T[:, i] = T1
            t[i] = t[i-1] + Δt
            T0 = T1
        end
    else
        A = C + ϑ * Δt * K
        b = C - (1 - ϑ) * Δt * K
        AA = lu(A)
        for i in 2:nsteps
            bb = b * T0 + Δt * q
            T1 = AA \ bb
            T[:, i] = T1
            t[i] = t[i-1] + Δt
            T0 = T1
        end
    end

    return T, t
end

"""
    FEM.CDMaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=..., ξ=..., β=..., show_β=..., show_ξ=...)

Gives some functions (graphs) for accuracy analysis of the CDM method. 
`ωₘᵢₙ` and `ωₘₐₓ` are the square root of smallest and largest eigenvalues of the
**Kϕ**=ω²**Mϕ** eigenvalue problem, `Δt` is the time step size. `type` is one of the
following values:
- `:SR`: spectral radius
- `:PDR`: physical damping ratio
- `:ADR`: algorithmic damping ratio
- `:PE`: period error
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
        if type == :SR
            x[i] = log((ω[i] / 2π) * Δt)
            y[i] = ρ
        elseif type == :ADR
            x[i] = (ω[i] / 2π) * Δt
            Ω0 = √((log(ρ))^2 + (atan(ε,σ))^2 / 4)
            y[i] = -log(ρ) / 2Ω0
        elseif type == :PDR
            x[i] = (ω[i] / 2π) * Δt
            for j in 1:length(β)
                y[i] += β[j] / 2 * (2π * x[i] / Δt) ^ (2j-1)
            end
        elseif type == :PE
            x[i] = (ω[i] / 2π) * Δt
            Ω0 = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = 1 - Ω0/(Δt*ω[i])
        else
            str1 = "CDMaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\n:SR - spectral radius\n"
            str4 = ":PDR - physical damping ratio\n"
            str5 = ":ADR - algorithmic damping ratio\n"
            str6 = ":PE - period error\n"
            str7 = "\nFor details see Serfőző, D., Pere, B.: A method to accurately define arbitrary\n"
            str8 = "algorithmic damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).\n"
            str9 = "https://doi.org/10.1007/s00419-023-02454-9\n"
            error(str1*str2*str3*str4*str5*str6*str7*str8*str9)
        end
    end
    return x, y
end

"""
    FEM.HHTaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, δ=0.0, γ=0.5 + δ, β=0.25 * (0.5 + γ)^2)

Gives some functions (graphs) for accuracy analysis of the HHT-α method[^1]. 
`ωₘᵢₙ` and `ωₘₐₓ` are the square root of smallest and largest eigenvalues of the
**Kϕ**=ω²**Mϕ** eigenvalue problem, `Δt` is the time step size. `type` is one of the
following values:
- `:SR`: spectral radius
- `:ADR`: algorithmic damping ratio
- `:PE`: period error
For details see [^2] and [^3]. 
`n` is the number of points in the graph. For the meaning of `α`, `β` and `γ`
see [^1]. If `δ` is given, γ=0.5+δ and β=0.25⋅(0.5+γ)².
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

[^2]: Belytschko, Ted, and Thomas JR, Hughes: "Computational methods for 
    transient analysis", North-Holland, (1983).

[^3]: Serfőző, D., Pere, B.: A method to accurately define arbitrary algorithmic
    damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).
    <https://doi.org/10.1007/s00419-023-02454-9>

Return: `xy`

Types:
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
        if type == :SR
            x[i] = log(invT[i] * Δt)
            y[i] = ρ
        elseif type == :ADR
            x[i] = invT[i] * Δt
            Ω = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = -log(ρ) / 2Ω
            #y[i] = -log(ρ) / atan(ε, σ)
        elseif type == :PE
            x[i] = invT[i] * Δt
            Ω = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = 1 - Ω/(2π*Δt*invT[i])
        else
            str1 = "HHTaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\n:SR - spectral radius\n"
            str5 = ":ADR - algorithmic damping ratio\n"
            str6 = ":PE - period error\n"
            str7 = "\nFor details see Serfőző, D., Pere, B.: A method to accurately define arbitrary\n"
            str8 = "algorithmic damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).\n"
            str9 = "https://doi.org/10.1007/s00419-023-02454-9\n"
            error(str1*str2*str3*str5*str6*str7*str8*str9)
        end
    end
    return x, y
end

"""
    FEM.FDMaccuracyAnalysis(λₘᵢₙ, λₘₐₓ, Δt; type, n=100, ϑ=...)

Gives a functions (graphs) for accuracy analysis of the ϑ-method[^5]. 
`λₘᵢₙ` and `λₘₐₓ` are the smallest and largest eigenvalues of the
**Kθ**=λ**Cθ** eigenvalue problem, `Δt` is the time step size. `type` is the
"SR" spectral radius.
`n` is the number of points in the graph. For the meaning of `ϑ` see [^5].
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

- "SR": spectral radius as a function of λ⋅Δt

Return: `xy`

Types:
- `λₘᵢₙ`: Float64
- `λₘₐₓ`: Float64
- `Δt`: Float64 
- `n`: Int64
- `ϑ`: Float64
- `xy`: Tuple{Vector{Float64},Vector{Float64}}
"""

function FDMaccuracyAnalysis(λₘᵢₙ, λₘₐₓ, Δt; type=:SR, n=100, ϑ=0.5)
    x = zeros(n)
    y = similar(x)
    Λ = range(λₘᵢₙ, length=n, stop=λₘₐₓ)
    for i ∈ 1:n
        A1 = 1 - (1 - ϑ) * Δt * Λ[i]
        A2 = 1 + ϑ * Δt * Λ[i]

        A = A1 \ A2

        ρ = abs(A)
        if type == :SR
            x[i] = Λ[i] * Δt
            y[i] = ρ
        else
            str1 = "FDMaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\nSR: spectral radius\n"
            error(str1*str2*str3)
        end
    end
    return x, y
end

"""
    FEM.rotateNodes(problem, phName, CoordSys)

Creates the `T` transformation matrix, which rotates the nodal coordinate system
of the nodes in `phName` physical group to the coordinate systen defined by `CoordSys`.
The mesh belongs to `problem`.

Return: `T`

Types:
- `problem`: Problem
- `phName`: String
- `CoordSys`: CoordinateSystem
- `T`: SparseMatrix
"""
function rotateNodes(problem, phName, CoordSys)
    dim = problem.dim
    non = problem.non
    dof = non * dim
    phg = getTagForPhysicalName(phName)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
    nonz = length(nodeTags) * dim * dim + (non - length(nodeTags)) * dim
    lengthOfIJV = length(nodeTags) * dim * dim + non * dim
    I = []
    J = []
    V = []
    V = convert(Vector{Float64}, V)
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)
    vec2 = [0.0, 0.0, 0.0]
    fn(x,y) = y
    IJ = 1:dof
    VV = ones(dof)
    append!(I, IJ)
    append!(J, IJ)
    append!(V, VV)
    
    if dim == 3
        I0 = [1, 2, 3, 1, 2, 3, 1, 2, 3]
        J0 = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        e1 = CoordSys.vec1
        e1 ./= √(dot(e1, e1))
        I3 = [1 0 0; 0 1 0; 0 0 1]
        e2 = (I3 - e1 * e1') * CoordSys.vec2
        e2 ./= √(dot(e2, e2))
        e3 = [e1[2] * e2[3] - e1[3] * e2[2], e1[3] * e2[1] - e1[1] * e2[3], e1[1] * e2[2] - e1[2] * e2[1]]
        T1 = [e1 e2 e3]

        if CoordSys.i1[1] == 1 || CoordSys.i1[2] == 1 || CoordSys.i1[3] == 1 || CoordSys.i2[1] == 1 || CoordSys.i2[2] == 1 || CoordSys.i2[3] == 1
            for i in 1:length(nodeTags)
                x = coord[i * 3 - 2]
                y = coord[i * 3 - 1]
                z = coord[i * 3]
                e1[1] = CoordSys.i1[1] == 1 ? CoordSys.vec1f[1](x, y, z) : CoordSys.vec1[1]
                e1[2] = CoordSys.i1[2] == 1 ? CoordSys.vec1f[2](x, y, z) : CoordSys.vec1[2]
                e1[3] = CoordSys.i1[3] == 1 ? CoordSys.vec1f[3](x, y, z) : CoordSys.vec1[3]
                vec2[1] = CoordSys.i2[1] == 1 ? CoordSys.vec2f[1](x, y, z) : CoordSys.vec2[1]
                vec2[2] = CoordSys.i2[2] == 1 ? CoordSys.vec2f[2](x, y, z) : CoordSys.vec2[2]
                vec2[3] = CoordSys.i2[3] == 1 ? CoordSys.vec2f[3](x, y, z) : CoordSys.vec2[3]
                e1 ./= √(dot(e1, e1))
                e2 = (I3 - e1 * e1') * vec2
                e2 ./= √(dot(e2, e2))
                e3 = [e1[2] * e2[3] - e1[3] * e2[2], e1[3] * e2[1] - e1[1] * e2[3], e1[1] * e2[2] - e1[2] * e2[1]]
                T1 = [e1 e2 e3]
                Iidx = (nodeTags[i] * 3 - 3) .+ I0
                Jidx = (nodeTags[i] * 3 - 3) .+ J0
                append!(I, Iidx)
                append!(J, Jidx)
                append!(V, T1[:])
            end
        else
            for i in 1:length(nodeTags)
                Iidx = (nodeTags[i] * 3 - 3) .+ I0
                Jidx = (nodeTags[i] * 3 - 3) .+ J0
                append!(I, Iidx)
                append!(J, Jidx)
                append!(V, T1[:])
            end
        end
    elseif dim == 2
        I0 = [1, 2, 1, 2]
        J0 = [1, 1, 2, 2]
        e1 = CoordSys.vec1[1:2]
        e1 ./= √(dot(e1, e1))
        e2 = [-e1[2], e1[1]]
        T1 = [e1 e2]

        if CoordSys.i1[1] == 1 || CoordSys.i1[2] == 1
            for i in 1:length(nodeTags)
                x = coord[i * 3 - 2]
                y = coord[i * 3 - 1]
                e1[1] = CoordSys.i1[1] == 1 ? CoordSys.vec1f[1](x, y) : CoordSys.vec1[1]
                e1[2] = CoordSys.i1[2] == 1 ? CoordSys.vec1f[2](x, y) : CoordSys.vec1[2]
                e1 ./= √(dot(e1, e1))
                e2 = [-e1[2], e1[1]]
                T1 = [e1 e2]
                Iidx = (nodeTags[i] * 2 - 2) .+ I0
                Jidx = (nodeTags[i] * 2 - 2) .+ J0
                append!(I, Iidx)
                append!(J, Jidx)
                append!(V, T1[:])
            end
        else
            for i in 1:length(nodeTags)
                Iidx = (nodeTags[i] * 2 - 2) .+ I0
                Jidx = (nodeTags[i] * 2 - 2) .+ J0
                append!(I, Iidx)
                append!(J, Jidx)
                append!(V, T1[:])
            end
        end
    else
        error("rotateNodes: dimension of the problem is 2 o 3, now it is $dim.")
    end
    T = sparse(I, J, V, dof, dof, fn)
    dropzeros!(T)
    return Transformation(T, non, dim)
end

"""
    FEM.showDoFResults(problem, q, comp; t=..., name=..., visible=...)

Loads nodal results into a View in gmsh. `q` is the field to show, `comp` is
the component of the field (:vector, :uvec, :ux, :uy, :uz, :vvec, :vx, :vy, :vz,
:qvec, :qx, :qy, :qz, :T, :p, :qn, :s, :sx, :sy, :sz, :sxy, :syx, :syz,
:szy, :szx, :sxz, :e, :ex, :ey, :ez, :exy, :eyx, :eyz, :ezy, :ezx, :exz, :seqv, :scalar),
`t` is a vector of time steps (same number of columns as `q`), `name` is a
title to display and `visible` is a true or false value to toggle on or off the 
initial visibility in gmsh. If `q` has more columns, then a sequence of results
will be shown (eg. as an animation). This function returns the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `q`: Vector{Matrix}
- `comp`: Symbol
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showDoFResults(problem, q, comp; t=[0.0], name=comp, visible=false, ff = 0)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    pdim = problem.pdim
    pdim = div(size(q,1), problem.non) 
    nodeTags = []
    ##############################################################################
    if problem.type == :Reynolds || problem.type == :NavierStokes
        phName = problem.material.phName
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
    if comp == :s || comp == :e
        pdim = 9
    end
    for j in 1:length(t)
        k = 1im
        if comp == :uvec || comp == :vvec || comp == :qvec || comp == :vector
            nc = 3
            u = zeros(3 * non)
            for i in 1:length(nodeTags)
                u[3i-2] = q[pdim*nodeTags[i]-(pdim-1), j]
                u[3i-1] = pdim > 1 ? q[pdim*nodeTags[i]-(pdim-2), j] : 0
                u[3i-0] = pdim == 3 ? q[pdim*nodeTags[i]-(pdim-3), j] : 0
            end
        elseif comp == :s || comp == :e
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
            if comp == :ux || comp == :vx || comp == :p || comp == :T || comp == :qx  || comp == :qn || comp == :sx || comp == :ex || comp == :scalar
                k = 1
            elseif comp == :uy || comp == :vy || comp == :qy || comp == :syx || comp == :eyx
                k = 2
            elseif comp == :uz || comp == :vz || comp == :qz || comp == :szx || comp == :ezx
                k = 3
            elseif comp == :sxy || comp == :exy
                k = 4
            elseif comp == :sy || comp == :ey
                k = 5
            elseif comp == :szy || comp == :ezy
                k = 6
            elseif comp == :sxz || comp == :exz
                k = 7
            elseif comp == :syz || comp == :eyz
                k = 8
            elseif comp == :sz || comp == :ez
                k = 9
            elseif comp == :seqv
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
    if ff == 1 || ff == 2
        gmsh.view.option.setNumber(uvec, "AdaptVisualizationGrid", 1)
    else
        gmsh.view.option.setNumber(uvec, "AdaptVisualizationGrid", 0)
    end
    gmsh.view.option.setNumber(uvec, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uvec, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(uvec, "Visible", 0)
    end
    if ff == 0 && length(t) > 1
        gmsh.view.option.setNumber(uvec, "ShowTime", 1)
    elseif ff == 1 || ff == 2
        gmsh.view.option.setNumber(uvec, "ShowTime", 6)
        gmsh.view.option.setNumber(uvec, "VectorType", 5)
    end
    return uvec
end

"""
    FEM.showModalResults(problem, Φ, name=..., visible=...)

Loads modal results into a View in gmsh. `Φ` is a struct of Eigen. `name` is a
title to display and `visible` is a true or false value to toggle on or off the 
initial visibility in gmsh. Click on ▷| to change the results. This function 
returns the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `Φ`: Eigen
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showModalResults(problem, Φ::Eigen; name=:modal, visible=false, ff=1)
    return showDoFResults(problem, Φ.ϕ, :uvec, t=Φ.f, name=name, visible=visible, ff=ff)
end

"""
    FEM.showBucklingResults(problem, Φ, name=..., visible=...)

Loads buckling results into a View in gmsh. `Φ` is a struct of Eigen. `name` is a
title to display and `visible` is a true or false value to toggle on or off the 
initial visibility in gmsh. Click on ▷| to change the results. This function 
returns the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `Φ`: Eigen
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showBucklingResults(problem, Φ::Eigen; name="buckling", visible=false, ff=2)
    return showDoFResults(problem, Φ.ϕ, :uvec, t=Φ.f, name=name, visible=visible, ff=ff)
end

"""
    FEM.showStrainResults(problem, E, comp; t=..., name=..., visible=..., smooth=...)

Loads strain results into a View in gmsh. `E` is a strain field to show, `comp` is
the component of the field (:e, :ex, :ey, :ez, :exy, :eyz, :ezx),
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
- `comp`: Symbol
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStrainResults(problem, E, comp; t=[0.0], name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
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
        if comp == :e
            εcomp = [ε[i][:,jj] for i in 1:length(E.numElem)]
            nc = 9
        else
            nc = 1
            if comp == :ex
                k = 8
            elseif comp == :ey
                k = 4
            elseif comp == :ez
                k = 0
            elseif comp == :exy || comp == :eyx
                k = 7
            elseif comp == :eyz || comp == :ezy
                k = 3
            elseif comp == :ezx || comp == :exz
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
- `comp`: Symbol
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showElementResults(problem, F, comp; t=[0.0], name=comp, visible=false, smooth=true)
    if F.type == :e
        return showStrainResults(problem, F, comp; t=t, name=comp, visible=false, smooth=smooth)
    elseif F.type == :s
        return showStressResults(problem, F, comp; t=t, name=comp, visible=false, smooth=smooth)
    elseif F.type == :q
        return showHeatFluxResults(problem, F, comp; t=t, name=comp, visible=false, smooth=smooth)
    else
        error("showElementResults: type is '$type'")
    end
end

"""
    FEM.showStressResults(problem, S, comp; t=..., name=..., visible=..., smooth=...)

Loads stress results into a View in gmsh. `S` is a stress field to show, `comp` is
the component of the field (:s, :sx, :sy, :sz, :sxy, :syz, :szx, :seqv),
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
- `comp`: Symbol
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStressResults(problem, S, comp; t=[0.0], name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
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
        if comp == :s
            σcomp = [σ[i][:,jj] for i in 1:length(S.numElem)]
            nc = 9
        elseif comp == :seqv
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
            if comp == :sx
                k = 8
            elseif comp == :sy
                k = 4
            elseif comp == :sz
                k = 0
            elseif comp == :sxy || comp == :syx
                k = 7
            elseif comp == :syz || comp == :szy
                k = 3
            elseif comp == :szx || comp == :sxz
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
the component of the field (:qvec, :qx, :qy, :qz, :q),
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
- `comp`: Symbol
- `t`: Vector{Float64}
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showHeatFluxResults(problem, S, comp; t=[0.0], name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
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
        if comp == :qvec
            σcomp = [σ[i][:,jj] for i in 1:length(S.numElem)]
            nc = 3
        elseif comp == :q
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
            if comp == :qx
                k = 2
            elseif comp == :qy
                k = 1
            elseif comp == :qz
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
    FEM.plotOnPath(problem, pathName, field; points=100, step=..., plot=..., name=..., visible=..., offsetX=..., offsetY=..., offsetZ=...)

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
function plotOnPath(problem, pathName, field; points=100, step=1im, plot=false, name="field [$field] on $pathName", visible=false, offsetX=0, offsetY=0, offsetZ=0)
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
                val, dis = gmsh.view.probe(field, pt1[1]+offsetX, pt1[2]+offsetY, pt1[3]+offsetZ, stepRange[j] - 1, -1, false, -1)
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
        step = 0
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

"""
    FEM.showOnSurface(field, phName; grad=false, component=:x, offsetX=0, offsetY=0, offsetZ=0, name=phName, visible=false)

Shows the values of a scalar field at a surface which has a physical name `phName`.
`field` is the tag of a view in GMSH. The values of the field are calculated at the
intersection with the surface. `grad` has a true or false value to toggle on or off
the gradient of the field. `component` is the component of the gradient of `field`
(:x, :y, :z) to be shown. `offsetX`, `offsetY`, `offsetZ` are the offsets in the
x, y and z directions where the values are picked from. `name` is a title to display
and `visible` is a true or false value to toggle on or off the initial visibility in gmsh.

Return: `tag`

Types:
- `field`: Integer
- `phName`: String
- `grad`: Boolean
- `component`: Symbol
- `offsetX`: Float64
- `offsetY`: Float64
- `offsetZ`: Float64
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showOnSurface(field, phName; grad=false, component=:x, offsetX=0, offsetY=0, offsetZ=0, name=phName, visible=false)
    SS = gmsh.view.add(name)
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
    nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(-1, -1, true, false)
    ncoord2 = similar(ncoord)
    ncoord2[nodeTags*3 .- 2] .= ncoord[1:3:length(ncoord)]
    ncoord2[nodeTags*3 .- 1] .= ncoord[2:3:length(ncoord)]
    ncoord2[nodeTags*3 .- 0] .= ncoord[3:3:length(ncoord)]
    if component == :x
        comp0 = 1
    elseif component == :y
        comp0 = 2
    elseif component == :z
        comp0 = 3
    else
        error("resultant: invalid component '$component'")
    end
    ret2 = []
    ret3 = []
    ret4 = []
    el2 = 0
    el3 = 0
    el4 = 0
    x = [0.0, 0.0, 0.0]
    for idm in 1:length(dimTags)
        dimTag = dimTags[idm]
        edim = dimTag[1]
        etag = dimTag[2]
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
        for i in 1:length(elemTypes)
            et = elemTypes[i]
            elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
            nn = zeros(numPrimaryNodes * 4)
            for j in 1:length(elemTags[i])
                elem = elemTags[i][j]
                for l in 1:numPrimaryNodes
                    for k in 1:3
                        x[k] = ncoord2[elemNodeTags[i][j*numNodes-(numNodes-l)]*3-(3-k)]
                        nn[k*numPrimaryNodes-numPrimaryNodes + l] = x[k]
                    end
                    k = 4
		    f, d = gmsh.view.probe(field, x[1]+offsetX, x[2]+offsetY, x[3]+offsetZ, -1, -1, grad, -1)
                    nn[4*numPrimaryNodes-numPrimaryNodes + l] = f[comp0]
                end
                if numPrimaryNodes == 3
                append!(ret3,nn)
                elseif numPrimaryNodes == 4
                    append!(ret4,nn)
                elseif numPrimaryNodes == 2
                    append!(ret2,nn)
                end
            end
            if numPrimaryNodes == 3
                el3 += length(elemTags[i])
            elseif numPrimaryNodes == 4
                el4 += length(elemTags[i])
            elseif numPrimaryNodes == 2
                el2 += length(elemTags[i])
            end
        end
    end
    if el3 > 0
        gmsh.view.addListData(SS, "ST", el3, ret3)
    end
    if el4 > 0
        gmsh.view.addListData(SS, "SQ", el4, ret4)
    end
    if el2 > 0
        gmsh.view.addListData(SS, "SL", el2, ret2)
    end
    if visible == false
        gmsh.view.option.setNumber(SS, "Visible", 0)
    end
    return SS
end


"""
    FEM.openPreProcessor(; openGL=...)

Launches the GMSH preprocessor window with openGL disabled by default.

Return: none

Types:
- `openGL`: Boolean
"""
function openPreProcessor(; openGL=false)
    if openGL == false
        ENV["LIBGL_ALWAYS_SOFTWARE"] = "true"
    else
        ENV["LIBGL_ALWAYS_SOFTWARE"] = "false"
    end
    gmsh.fltk.run()
end

"""
    FEM.openPostProcessor(; model=...)

Launches the GMSH postprocessor window with open postprocessor tree (of `model`).

Return: none

Types:
- `model`: Int64
"""
function openPostProcessor(; model=0)
    gmsh.fltk.openTreeItem(LazyString(model)*"Modules/Post-processing")
    gmsh.fltk.run()
end

"""
    FEM.setParameter(name, value)

Defines a parameter `name` and sets its value to `value`. 

Return: none

Types:
- `name`: String
- `value`: Float64
"""
function setParameter(name, value)
    gmsh.parser.setNumber(name, [value])
end

end #module
