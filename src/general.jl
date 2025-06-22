"""
    Material(phName, E, ν, ρ, k, c, α, λ, μ, κ)

A structure containing the material constants. 
- E: elastic modulus,
- ν: Poisson's ratio,
- ρ: mass density,
- k: heat conductivity,
- c: specific heat,
- α: heat expansion coefficient
- λ: Lamé parameter
- μ: Lamé parameter
- κ: Bulk modulus
`phName` is the name of the physical group where the given material is used.

Types:
- `phName`: String
- `E`: Float64
- `ν`: Float64
- `ρ`: Float64
- `k`: Float64
- `c`: Float64
- `α`: Float64
- `λ`: Float64
- `μ`: Float64
- `κ`: Float64
"""
struct Material
    phName::String
    E::Float64
    ν::Float64
    ρ::Float64
    k::Float64
    c::Float64
    α::Float64
    λ::Float64
    μ::Float64
    κ::Float64
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
        elseif type == :StVenantKirchhoff
            dim = 3
            pdim = 3
        elseif type == :NeoHookeCompressible
            dim = 3
            pdim = 3
        else
            error("Problem type can be: `:Solid`, `:PlaneStress`, `:PlaneStrain`, `:AxiSymmetric`, `:PlaneHeatConduction`, `:HeatConduction`, `:AxiSymmetricHeatConduction`,
            `StVenantKirchhoff` or `NeoHookeCompressible`. Now problem type = $type ????")
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
    ScalarField(A, a, t, numElem, nsteps, type)
    ScalarField(problem, dataField)

A structure containing the data of a heat flux field. 
- A: vector of ElementNodeData type scalar data (see gmsh.jl)
- numElem: vector of tags of elements
- nsteps: number of stress fields stored in `A` (for animations).
- type: type of data (eg. heat flux `:q`)

Types:
- `A`: Vector{Vector{Float64}}
- `a`: Matrix{Float64}
- `t`: Vector{Float64}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: Symbol
"""
struct ScalarField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
    function ScalarField(A0, a0, t0, numElem0, nsteps0, type0)
        return new(A0, a0, t0, numElem0, nsteps0, type0)
    end
    function ScalarField(problem, dataField)
        if !isa(dataField, Vector)
            error("ScalarField: dataField are not arranged in a vector. Put them in [...]")
        end
        gmsh.model.setCurrent(problem.name)

        type = :scalarInElements
        nsteps = 1
        A = []
        numElem = Int[]
        pdim = 1
        for i in 1:length(dataField)
            phName, f, fx, fy, fz, fxy, fyz, fzx = dataField[i]
            dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
            for idm in 1:length(dimTags)
                dimTag = dimTags[idm]
                edim = dimTag[1]
                etag = dimTag[2]
                elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
                for i in 1:length(elemTypes)
                    et = elemTypes[i]
                    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                    sc1 = zeros(numNodes, nsteps)
                    for j in 1:length(elemTags[i])
                        elem = elemTags[i][j]
                        push!(numElem, elem)
                        for k in 1:numNodes
                            nodeTag = elemNodeTags[i][(j-1)*numNodes+k]
                            if f != :no
                                if isa(f, Function)
                                    coord, parametricCoord, dim, tag = gmsh.model.mesh.getNode(nodeTag)
                                    x = coord[1]
                                    y = coord[2]
                                    z = coord[3]
                                    ff = f(x, y, z)
                                else
                                    ff = f
                                end
                                sc1[k] = f
                            end
                        end
                        push!(A, sc1)
                    end
                end
            end
        end
        a = [;;]
        t = []
        return new(A, a, t, numElem, nsteps, type)
    end
end

"""
    VectorField(A, a, t, numElem, nsteps, type)

A structure containing the data of a heat flux field. 
- A: vector of ElementNodeData type heat flux data (see gmsh.jl)
- numElem: vector of tags of elements
- nsteps: number of stress fields stored in `A` (for animations).
- type: type of data (eg. heat flux `:q`)

Types:
- `A`: Vector{Matrix{Float64}}
- `a`: Matrix{Float64}
- `t`: Vector{Float64}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: Symbol
"""
struct VectorField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
end

"""
    TensorField(A, a, t, numElem, nsteps, type)

A structure containing the data of a stress or strain field. 
- A: vector of ElementNodeData type stress data (see gmsh.jl)
- numElem: vector of tags of elements
- nsteps: number of stress fields stored in `A` (for animations).
- type: type of data (eg. stress `:s` and strain `:e`)

Types:
- `A`: Vector{Matrix{Float64}}
- `a`: Matrix{Float64}
- `t`: Vector{Float64}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: Symbol
"""
struct TensorField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
end

import Base.copy
function copy(A::ScalarField)
    a = copy(A.A)
    b = copy(A.a)
    c = copy(A.t)
    d = copy(A.numElem)
    e = copy(A.nsteps)
    f = A.type
    return ScalarField(a, b, c, d, e, f)
end

function copy(A::VectorField)
    a = copy(A.A)
    b = copy(A.a)
    c = copy(A.t)
    d = copy(A.numElem)
    e = copy(A.nsteps)
    f = A.type
    return VectorField(a, b, c, d, e, f)
end

function copy(A::TensorField)
    a = copy(A.A)
    b = copy(A.a)
    c = copy(A.t)
    d = copy(A.numElem)
    e = copy(A.nsteps)
    f = A.type
    return TensorField(a, b, c, d, e, f)
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
    FEM.material(name; E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5, λ=νE/(1+ν)/(1-2ν), μ=E/(1+ν)/2, κ=E/(1-2ν)/3)

Returns a structure in which `name` is the name of a physical group, 
`E` is the modulus of elasticity, `ν` Poisson's ratio and `ρ` is
the mass density, `k` is the heat conductivity, `c` is the specific
heat, `α` is the coefficient of heat expansion, `λ` and `μ` are the 
Lamé parameters, `κ` is the Bulk modulus.

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
- `λ`: Float64
- `μ`: Float64
- `κ`: Float64
"""
function material(name; E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5, μ=E/(1+ν)/2, λ=2μ*ν/(1-2ν), κ=2μ*(1+ν)/(1-2ν)/3)
    return Material(name, E, ν, ρ, k, c, α, λ, μ, κ)
end

import Base.*
function *(A::TensorField, B::TensorField)
    if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
        if length(A.A) != length(B.A)
            error("*(A::TensoeField, B::TensorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
        end
        if A.numElem != B.numElem
            error("*(A::TensoeField, B::TensorField): tensor fields are not compatible.")
        end
        nsteps = A.nsteps
        nsteps2 = B.nsteps
        if nsteps != nsteps2
            error("*(A::TensoeField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
        end
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            m = length(B.A[i]) ÷ 9
            if n != m
                error("*(A::TensoeField, B::TensorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
            end
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = reshape(reshape(A.A[i][9j-8:9j, k], 3, 3) * reshape(B.A[i][9j-8:9j, k], 3, 3), 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
    else
        error("*(A::TensorField, B::TensorField): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function *(A::ScalarField, B::TensorField)
    if length(A.A) != 0 && length(B.A) != 0
        sz = 0
        nsteps = B.nsteps
        sec = intersect(B.numElem, A.numElem)
        indS = []
        indT = []
        sizehint!(indS, length(sec))
        sizehint!(indT, length(sec))
        for i in sec
            append!(indS, findall(j -> j == i, A.numElem))
            append!(indT, findall(j -> j == i, B.numElem))
        end
        C = []
        num = []
        sizehint!(C, length(sec))
        sizehint!(num, length(sec))
        D = []
        for i in eachindex(sec)
            n = length(B.A[i]) ÷ 9
            if n != sz
                D = zeros(9n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][9j-8:9j, k]
                end
            end
            append!(num, sec[i])
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, B.t, num, B.nsteps, :e)
    else
        error("*(ScalarField, TensorField): data at nodes is not yet implemented.")
    end
end

function *(B::TensorField, A::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
        sz = 0
        nsteps = B.nsteps
        sec = intersect(B.numElem, A.numElem)
        indS = []
        indT = []
        sizehint!(indS, length(sec))
        sizehint!(indT, length(sec))
        for i in sec
            append!(indS, findall(j -> j == i, A.numElem))
            append!(indT, findall(j -> j == i, B.numElem))
        end
        C = []
        num = []
        sizehint!(C, length(sec))
        sizehint!(num, length(sec))
        D = []
        for i in eachindex(sec)
            n = length(B.A[i]) ÷ 9
            if n != sz
                D = zeros(9n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][9j-8:9j, k]
                end
            end
            append!(num, sec[i])
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, B.t, num, B.nsteps, :e)
    else
        error("*(TensorField, ScalarField): data at nodes is not yet implemented.")
    end
end

import Base./
function /(B::TensorField, A::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
        sz = 0
        nsteps = B.nsteps
        sec = intersect(B.numElem, A.numElem)
        indS = []
        indT = []
        sizehint!(indS, length(sec))
        sizehint!(indT, length(sec))
        for i in sec
            append!(indS, findall(j -> j == i, A.numElem))
            append!(indT, findall(j -> j == i, B.numElem))
        end
        C = []
        num = []
        sizehint!(C, length(sec))
        sizehint!(num, length(sec))
        D = []
        for i in eachindex(sec)
            n = length(B.A[i]) ÷ 9
            if n != sz
                D = zeros(9n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = B.A[indT[i]][9j-8:9j, k] / A.A[indS[i]][j, k]
                end
            end
            append!(num, sec[i])
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, B.t, num, B.nsteps, :e)
    else
        error("/(TensorField, ScalarField): data at nodes is not yet implemented.")
    end
end

function *(A::ScalarField, B::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
        sz = 0
        if A.nsteps != B.nsteps
            error("*(ScalarField, ScalarField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
        end
        nsteps = B.nsteps
        sec = intersect(B.numElem, A.numElem)
        indS = []
        indT = []
        sizehint!(indS, length(sec))
        sizehint!(indT, length(sec))
        for i in sec
            append!(indS, findall(j -> j == i, A.numElem))
            append!(indT, findall(j -> j == i, B.numElem))
        end
        C = []
        num = []
        sizehint!(C, length(sec))
        sizehint!(num, length(sec))
        D = []
        for i in eachindex(sec)
            n = length(B.A[i])
            if n != sz
                D = zeros(n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    D[j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][j, k]
                end
            end
            append!(num, sec[i])
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, num, A.nsteps, :e)
    else
        error("*(ScalarField, ScalarField): data at nodes is not yet implemented.")
    end
end

function /(A::ScalarField, B::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
        sz = 0
        if A.nsteps != B.nsteps
            error("/(ScalarField, ScalarField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
        end
        nsteps = B.nsteps
        sec = intersect(B.numElem, A.numElem)
        indS = []
        indT = []
        sizehint!(indS, length(sec))
        sizehint!(indT, length(sec))
        for i in sec
            append!(indS, findall(j -> j == i, A.numElem))
            append!(indT, findall(j -> j == i, B.numElem))
        end
        C = []
        num = []
        sizehint!(C, length(sec))
        sizehint!(num, length(sec))
        D = []
        for i in eachindex(sec)
            n = length(B.A[i])
            if n != sz
                D = zeros(n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    D[j, k] = A.A[indS[i]][j, k] / B.A[indT[i]][j, k]
                end
            end
            append!(num, sec[i])
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, num, A.nsteps, :e)
    else
        error("/(ScalarField, ScalarField): data at nodes is not yet implemented.")
    end
end

import Base.transpose
function transpose(A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                D = zeros(9n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        D[9j-8:9j, k] = reshape(transpose(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                    end
                end
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("transpose(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
        end
    else
        error("transpose(TensorField): data at nodes is not yet implemented.")
    end
end

import Base.adjoint
function adjoint(A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                D = zeros(9n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        D[9j-8:9j, k] = reshape(adjoint(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                    end
                end
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("adjoint(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
        end
    else
        error("adjoint(TensorField): data at nodes is not yet implemented.")
    end
end

function unitTensor(A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                D = zeros(9n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        D[9j-8:9j, k] = reshape([1 0 0; 0 1 0; 0 0 1], 9, 1)
                    end
                end
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("unit(A::TensorField): TensorField type ($(A.type) is not yet implemented.")
        end
    else
        error("transpose(TensorField): data at nodes is not yet implemented.")
    end
end

function trace(A::TensorField)
    if length(A.A) != 0
        sz = 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            D = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                if sz != n
                    D = zeros(n, nsteps)
                    sz = n
                end
                for j in 1:n
                    for k in 1:nsteps
                        trace = A.A[i][9j-8, k] + A.A[i][9j-4, k] + A.A[i][9j, k]
                        D[j, k] = trace
                    end
                end
                push!(C, D)
            end
            a = [;;]
            return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("trace(A::TensorField): TensorField type ($(A.type) is not yet implemented.")
        end
    else
        error("trace(TensorField): data at nodes is not yet implemented.")
    end
end

import LinearAlgebra.det
function det(A::TensorField)
    if length(A.A) != 0
        sz = 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            D = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                if sz != n
                    D = zeros(n, nsteps)
                    sz = n
                end
                for j in 1:n
                    for k in 1:nsteps
                        d = LinearAlgebra.det(reshape(A.A[i][9j-8:9j, k], 3, 3))
                        D[j, k] = d
                    end
                end
                push!(C, D)
            end
            a = [;;]
            return ScalarField(C, a, A.t, A.numElem, A.nsteps, :sc)
        else
            error("det(A::TensorField): TensorField type ($(A.type) is not yet implemented.")
        end
    else
        error("det(TensorField): data at nodes is not yet implemented.")
    end
end

import Base.+
function +(A::ScalarField, B::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
        if A.type == B.type
            if length(A.A) != length(B.A)
                error("+(A::ScalarField, B::ScalarField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("+(A::ScalarField, B::ScalarField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            for i in eachindex(sec)
                #n = length(A.A[i]) ÷ 9
                #m = length(B.A[i]) ÷ 9
                #if n != m
                #    error("+(A::VectorField, B::VectorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
                #end
                D = A.A[i] + B.A[i]
                append!(num, sec[i])
                push!(C, D)
            end
            for i in eachindex(dif1)
                D = A.A[i]
                append!(num, dif1[i])
                push!(C, D)
            end
            for i in eachindex(dif2)
                D = B.A[i]
                append!(num, dif2[i])
                push!(C, D)
            end
            a = [;;]
            return ScalarField(C, a, A.t, num, A.nsteps, A.type)
        else
            error("+(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if A.type == B.type
            return ScalarField([], A.a + B.a, A.t, [], A.nsteps, A.type)
        else
            error("+(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("+(ScalarField, ScalarField): internal error")
    end
end

import Base.-
function -(A::ScalarField, B::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
        if A.type == B.type
            if length(A.A) != length(B.A)
                error("-(A::ScalarField, B::ScalarField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("-(A::ScalarField, B::ScalarField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            for i in eachindex(sec)
                #n = length(A.A[i]) ÷ 9
                #m = length(B.A[i]) ÷ 9
                #if n != m
                #    error("+(A::VectorField, B::VectorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
                #end
                D = A.A[i] - B.A[i]
                append!(num, sec[i])
                push!(C, D)
            end
            for i in eachindex(dif1)
                D = A.A[i]
                append!(num, dif1[i])
                push!(C, D)
            end
            for i in eachindex(dif2)
                D = -B.A[i]
                append!(num, dif2[i])
                push!(C, D)
            end
            a = [;;]
            return ScalarField(C, a, A.t, num, A.nsteps, A.type)
        else
            error("+(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if A.type == B.type
            return ScalarField([], A.a - B.a, A.t, [], A.nsteps, A.type)
        else
            error("+(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("+(ScalarField, ScalarField): internal error")
    end
end

import Base.+
function +(A::VectorField, B::VectorField)
    if length(A.A) != 0 && length(B.A) != 0
        if (A.type == :u3D && B.type == :u3D) || (A.type == :u2D && B.type == :u2D) || (A.type == :f3D && B.type == :f3D) || (A.type == :f2D && B.type == :f2D)
            if length(A.A) != length(B.A)
                error("+(A::VectorField, B::VectorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("+(A::VectorField, B::VectorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            for i in eachindex(sec)
                #n = length(A.A[i]) ÷ 9
                #m = length(B.A[i]) ÷ 9
                #if n != m
                #    error("+(A::VectorField, B::VectorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
                #end
                D = A.A[i] + B.A[i]
                append!(num, sec[i])
                push!(C, D)
            end
            for i in eachindex(dif1)
                D = A.A[i]
                append!(num, dif1[i])
                push!(C, D)
            end
            for i in eachindex(dif2)
                D = B.A[i]
                append!(num, dif2[i])
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, num, A.nsteps, A.type)
        else
            error("+(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if (A.type == :u3D && B.type == :u3D) || (A.type == :u2D && B.type == :u2D) || (A.type == :f3D && B.type == :f3D) || (A.type == :f2D && B.type == :f2D)
            return VectorField([], A.a + B.a, A.t, [], A.nsteps, A.type)
        else
            error("+(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("+(VectorField, VectorField): internal error")
    end
end

import Base.-
function -(A::VectorField, B::VectorField)
    if length(A.A) != 0 && length(B.A) != 0
        if (A.type == :u3D && B.type == :u3D) || (A.type == :u2D && B.type == :u2D) || (A.type == :f3D && B.type == :f3D) || (A.type == :f2D && B.type == :f2D)
            if length(A.A) != length(B.A)
                error("-(A::VectorField, B::VectorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("-(A::VectorField, B::VectorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            for i in eachindex(sec)
                D = A.A[i] - B.A[i]
                append!(num, sec[i])
                push!(C, D)
            end
            for i in eachindex(dif1)
                D = A.A[i]
                append!(num, dif1[i])
                push!(C, D)
            end
            for i in eachindex(dif2)
                D = -B.A[i]
                append!(num, dif2[i])
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, num, A.nsteps, :e)
        else
            error("-(A::VectorField, B::VectorField): VectorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if (A.type == :u3D && B.type == :u3D) || (A.type == :u2D && B.type == :u2D) || (A.type == :f3D && B.type == :f3D) || (A.type == :f2D && B.type == :f2D)
            return VectorField([], A.a - B.a, A.t, [], A.nsteps, A.type)
        else
            error("-(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("-(VectorField, VectorField): internal error")
    end
end

function *(A::VectorField, b)
    if length(A.A) != 0
        if A.type == :u3D || A.type == :u2D || A.type == :f3D || A.type == :f2D
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                D = A.A[i] * b
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("*(A::VectorField, b): VectorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0
        return VectorField([], A.a .* b, A.t, [], A.nsteps, A.type)
    else
        error("*(VectorField, b): internal error")
    end
end

function *(b, A::VectorField)
    if length(A.A) != 0
        if A.type == :u3D || A.type == :u2D || A.type == :f3D || A.type == :f2D
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                D = A.A[i] * b
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("*(A::VectorField, b): VectorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0
        return VectorField([], A.a .* b, A.t, [], A.nsteps, A.type)
    else
        error("*(b, VectorField): internal error")
    end
end

function /(A::VectorField, b)
    if length(A.A) != 0
        if A.type == :u3D || A.type == :u2D || A.type == :f3D || A.type == :f2D
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                D = A.A[i] / b
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("/(A::VectorField, b): VectorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0
        return VectorField([], A.a / b, A.t, [], A.nsteps, A.type)
    else
        error("/(VectorField, b): internal error")
    end
end

function +(A::TensorField, B::TensorField)
    if length(A.A) != 0 && length(B.A) != 0
        if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
            if length(A.A) != length(B.A)
                error("+(A::TensoeField, B::TensorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("+(A::TensoeField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            for i in eachindex(sec)
                #n = length(A.A[i]) ÷ 9
                #m = length(B.A[i]) ÷ 9
                #if n != m
                #    error("+(A::TensorField, B::TensorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
                #end
                D = A.A[i] + B.A[i]
                append!(num, sec[i])
                push!(C, D)
            end
            for i in eachindex(dif1)
                D = A.A[i]
                append!(num, dif1[i])
                push!(C, D)
            end
            for i in eachindex(dif2)
                D = B.A[i]
                append!(num, dif2[i])
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, num, A.nsteps, :e)
        else
            error("+(A::TensorField, B::TensorField): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    else
        error("+(TensorField, TensorField): data at nodes is not yet implemented.")
    end
end

import Base.-
function -(A::TensorField, B::TensorField)
    if length(A.A) != 0 && length(B.A) != 0
        if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
            if length(A.A) != length(B.A)
                error("-(A::TensoeField, B::TensorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("-(A::TensoeField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            for i in eachindex(sec)
                D = A.A[i] - B.A[i]
                append!(num, sec[i])
                push!(C, D)
            end
            for i in eachindex(dif1)
                D = A.A[i]
                append!(num, dif1[i])
                push!(C, D)
            end
            for i in eachindex(dif2)
                D = -B.A[i]
                append!(num, dif2[i])
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, num, A.nsteps, :e)
        else
            error("-(A::TensorField, B::TensorField): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    else
        error("-(TensorField, TensorField): data at nodes is not yet implemented.")
    end
end

function *(A::TensorField, b)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                D = A.A[i] * b
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("*(A::TensorField, b): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    else
        error("*(TensorField, Any): data at nodes is not yet implemented.")
    end
end

function *(b, A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                D = A.A[i] * b
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("*(A::TensorField, b): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    else
        error("*(Any, TensorField): data at nodes is not yet implemented.")
    end
end

function *(A::ScalarField, b)
    if length(A.A) != 0
        C = []
        for i in 1:length(A.A)
            D = A.A[i] * b
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e)
    else
        error("*(ScalarField, Any): data at nodes is not yet implemented.")
    end
end

function *(b, A::ScalarField)
    if length(A.A) != 0
        C = []
        for i in 1:length(A.A)
            D = A.A[i] * b
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e)
    else
        error("*(Any, ScalarField): data at nodes is not yet implemented.")
    end
end

function /(A::TensorField, b)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            C = []
            for i in 1:length(A.A)
                D = A.A[i] / b
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("/(A::TensorField, b): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    else
        error("/(TensorField, Any): data at nodes is not yet implemented.")
    end
end

function /(A::ScalarField, b)
    if length(A.A) != 0
        C = []
        for i in 1:length(A.A)
            D = A.A[i] / b
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e)
    else
        error("/(ScalarField, Any): data at nodes is not yet implemented.")
    end
end

function /(b, A::ScalarField)
    if length(A.A) != 0
        C = []
        for i in 1:length(A.A)
            D = b ./ A.A[i]
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e)
    else
        error("/(Any, ScalarField): data at nodes is not yet implemented.")
    end
end

import Base.inv
function inv(A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                D = zeros(9n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        D[9j-8:9j, k] = reshape(inv(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                    end
                end
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("inv(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
        end
    else
        error("inv(TensorField): data at nodes is not yet implemented.")
    end
end

import Base.sqrt
function sqrt(A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                D = zeros(9n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        D[9j-8:9j, k] = reshape(sqrt(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                    end
                end
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("sqrt(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
        end
    else
        error("sqrt(TensorField): data at nodes is not yet implemented.")
    end
end

import Base.log
function log(A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                D = zeros(9n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        D[9j-8:9j, k] = reshape(log(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                    end
                end
                push!(C, D)
            end
            a = [;;]
            return TensorField(C, a, A.t, A.numElem, A.nsteps, :e)
        else
            error("log(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
        end
    else
        error("log(TensorField): data at nodes is not yet implemented.")
    end
end

import Base.log
function log(A::ScalarField)
    if length(A.A) != 0
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i])
            D = zeros(n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[j, k] = log(A.A[i][j, k])
                end
            end
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e)
    else
        error("log(ScalarField): data at nodes is not yet implemented.")
    end
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
        error("scalarField: dataField are not arranged in a vector. Put them in [...]")
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
    return ScalarField([], reshape(field, :, 1), [0], [], 1, :scalarInNodes)
end

"""
    FEM.vectorField(problem, dataField; type=...)

Defines a vector field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`. `type` can be an arbitrary `Symbol`,
eg. `:u` or `:f`.

Return: VectorField

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
function vectorField(problem, dataField; type=:u)
    if !isa(dataField, Vector)
        error("applyBoundaryConditions!: dataField are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    pdim = problem.dim
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
    if pdim == 3
        type = Symbol(String(type) * "3D")
    elseif pdim == 2
        type = Symbol(String(type) * "2D")
    else
        error("vectorField: dimension is $pdim")
    end
    return VectorField([], reshape(field, :,1), [0], [], 1, type)
end

"""
    FEM.tensorField(problem, dataField; type=...)

Defines a vector field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`. `type` can be an arbitrary `Symbol`,
eg. `:u` or `:f`.

Return: TensorField

Types:
- `problem`: Problem
- `dataField`: Vector{Tuple{String, Float64,...}}

# Examples

```julia
f1(x, y, z) = sin(x)
f2(x, y, z) = 5y
ff1 = FEM.field("face1", fx=f1, fy=f2, fz=0, fxy=1, fyz=1, fzx=f2)
ff2 = FEM.field("face2", fx=f2, fy=f1, fz=1, fxy=1, fyz=f1, fzx=1)
qq = FEM.tensorField(problem, [ff1, ff2])
qq0 = FEM.showDoFResults(problem, qq, :tensor)
```
"""
function tensorField(problem, dataField; type=:e)
    if !isa(dataField, Vector)
        error("applyBoundaryConditions!: dataField are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    pdim = 9
    non = problem.non
    field = zeros(non * pdim)

    for i in 1:length(dataField)
        name, f, fx, fy, fz, fxy, fyz, fzx = dataField[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(fx, Function) || isa(fy, Function) || isa(fz, Function) || isa(fxy, Function) || isa(fyz, Function) || isa(fzx, Function)
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
            nodeTagsY .-= (pdim - 5)
            if isa(fy, Function)
                ffy = fy.(xx, yy, zz)
                field[nodeTagsY,:] .= ffy
            else
                field[nodeTagsY,:] .= fy
            end
        end
        if fz != :no
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= pdim
            if isa(fz, Function)
                ffz = fz.(xx, yy, zz)
                field[nodeTagsZ,:] .= ffz
            else
                field[nodeTagsZ,:] .= fz
            end
        end
        if fxy != :no
            nodeTagsXY = copy(nodeTags)
            nodeTagsXY *= pdim
            nodeTagsXY .-= (pdim - 4)
            if isa(fxy, Function)
                ffxy = fxy.(xx, yy, zz)
                field[nodeTagsXY,:] .= ffxy
            else
                field[nodeTagsXY,:] .= fxy
            end
        end
        if fxy != :no
            nodeTagsYX = copy(nodeTags)
            nodeTagsYX *= pdim
            nodeTagsYX .-= (pdim - 2)
            if isa(fxy, Function)
                ffxy = fxy.(xx, yy, zz)
                field[nodeTagsYX,:] .= ffxy
            else
                field[nodeTagsYX,:] .= fxy
            end
        end
        if fyz != :no
            nodeTagsYZ = copy(nodeTags)
            nodeTagsYZ *= pdim
            nodeTagsYZ .-= (pdim - 8)
            if isa(fyz, Function)
                ffyz = fyz.(xx, yy, zz)
                field[nodeTagsYZ,:] .= ffyz
            else
                field[nodeTagsYZ,:] .= fyz
            end
        end
        if fyz != :no
            nodeTagsZY = copy(nodeTags)
            nodeTagsZY *= pdim
            nodeTagsZY .-= (pdim - 6)
            if isa(fyz, Function)
                ffyz = fyz.(xx, yy, zz)
                field[nodeTagsZY,:] .= ffyz
            else
                field[nodeTagsZY,:] .= fyz
            end
        end
        if fzx != :no
            nodeTagsZX = copy(nodeTags)
            nodeTagsZX *= pdim
            nodeTagsZX .-= (pdim - 3)
            if isa(fyz, Function)
                ffyz = fyz.(xx, yy, zz)
                field[nodeTagsZX,:] .= ffyz
            else
                field[nodeTagsZX,:] .= fyz
            end
        end
        if fzx != :no
            nodeTagsXZ = copy(nodeTags)
            nodeTagsXZ *= pdim
            nodeTagsXZ .-= (pdim - 7)
            if isa(fyz, Function)
                ffyz = fyz.(xx, yy, zz)
                field[nodeTagsXZ,:] .= ffyz
            else
                field[nodeTagsXZ,:] .= fyz
            end
        end
    end
    return TensorField([], reshape(field, :,1), [0], [], 1, type)
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
    FEM.elementsToNodes(problem, T)

Solves the nodal results `F` from the elemental results `T`.
`T` can be tensor field or vector field.

Return: `F`

Types:
- `problem`: Problem
- `T`: TensorField or VectorField
- `F`: TensorField or VectorField
"""
function elementsToNodes(problem, S)
    gmsh.model.setCurrent(problem.name)

    type = S.type
    nsteps = S.nsteps
    numElem = S.numElem
    σ = S.A
    non = problem.non
    if type == :s || type == :e || type == :F
        epn = 9
    elseif type == :q3D
        epn = 3
    elseif type == :q2D
        epn = 2
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
    if type == :q3D || type == :q2D
        return VectorField([], s, S.t, [], S.nsteps, type)
    elseif type == :e || type == :s
        return TensorField([], s, S.t, [], S.nsteps, type)
    else
        error("elementsToNodes: internal error, type=$type.")
    end
end

"""
    FEM.fieldError(problem, F)

Solves the deviation of field results `F` (stresses, strains, heat flux components) at nodes, where the field has jumps.
The result can be displayed with the `showDoFResults` function.

Return: `e`

Types:
- `problem`: Problem
- `F`: TensorField or VectorField
- `e`: ScalarField
"""
function fieldError(problem, S)
    gmsh.model.setCurrent(problem.name)

    type = S.type
    nsteps = S.nsteps
    numElem = S.numElem
    σ = S.A
    non = problem.non
    if type == :s || type == :e
        epn = 9
    elseif type == :q3D
        epn = 3
    elseif type == :q2D
        epn = 2
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
    if type == :q3D || type == :q2D
        return VectorField([], res, S.t, [], S.nsteps, type)
    elseif type == :e || type == :s
        return TensorField([], res, S.t, [], S.nsteps, type)
    else
        error("elementsToNodes: internal error, type=$type.")
    end
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

Return: `res`

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
    FEM.showDoFResults(problem, q, comp; name=..., visible=...)

Loads nodal results into a View in gmsh. `q` is the field to show, `comp` is
the component of the field (:vector, :uvec, :ux, :uy, :uz, :vvec, :vx, :vy, :vz,
:qvec, :qx, :qy, :qz, :T, :p, :qn, :s, :sx, :sy, :sz, :sxy, :syx, :syz,
:szy, :szx, :sxz, :e, :ex, :ey, :ez, :exy, :eyx, :eyz, :ezy, :ezx, :exz, :seqv, :scalar),
`name` is a title to display and `visible` is a true or false value to toggle on or off the 
initial visibility in gmsh. If `q` has more columns, then a sequence of results
will be shown (eg. as an animation). This function returns the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `q`: ScalarField, VectorField or TensorField
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
    t = q.t
    if q.a == [;;]
        error("showDoFResults: No data")
    end
    dim = problem.dim
    pdim = problem.pdim
    pdim = div(size(q.a,1), problem.non) 
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
    if size(q.a, 2) != length(t)
        error("showDoFResults: number of time steps missmatch ($(size(q.a,2)) <==> $(length(t))).")
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
                u[3i-2] = q.a[pdim*nodeTags[i]-(pdim-1), j]
                u[3i-1] = pdim > 1 ? q.a[pdim*nodeTags[i]-(pdim-2), j] : 0
                u[3i-0] = pdim == 3 ? q.a[pdim*nodeTags[i]-(pdim-3), j] : 0
            end
        elseif comp == :s || comp == :e || comp == :tensor
            nc = 9
            u = zeros(9 * non)
            for i in 1:length(nodeTags)
                u[9i-8] = q.a[pdim*nodeTags[i]-(pdim-1), j]
                u[9i-7] = q.a[pdim*nodeTags[i]-(pdim-2), j]
                u[9i-6] = q.a[pdim*nodeTags[i]-(pdim-3), j]
                u[9i-5] = q.a[pdim*nodeTags[i]-(pdim-4), j]
                u[9i-4] = q.a[pdim*nodeTags[i]-(pdim-5), j]
                u[9i-3] = q.a[pdim*nodeTags[i]-(pdim-6), j]
                u[9i-2] = q.a[pdim*nodeTags[i]-(pdim-7), j]
                u[9i-1] = q.a[pdim*nodeTags[i]-(pdim-8), j]
                u[9i-0] = q.a[pdim*nodeTags[i]-(pdim-9), j]
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
                null = k <= 9 ? 0 : √(0.5 * ((q.a[pdim*nodeTags[i]-(pdim-1), j]-q[pdim*nodeTags[i]-(pdim-5), j])^2+(q[pdim*nodeTags[i]-(pdim-5), j]-q[pdim*nodeTags[i]-(pdim-9), j])^2+(q[pdim*nodeTags[i]-(pdim-9), j]-q[pdim*nodeTags[i]-(pdim-1), j])^2 + 6*(q[pdim*nodeTags[i]-(pdim-2), j]^2+q[pdim*nodeTags[i]-(pdim-3), j]^2+q[pdim*nodeTags[i]-(pdim-6), j]^2)))
                u[i] = k > pdim ? null : q.a[pdim*nodeTags[i]-(pdim-k), j]
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
    return showDoFResults(problem, VectorField([], Φ.ϕ, Φ.f, [], length(Φ.f), :u3D), :uvec, name=name, visible=visible, ff=ff)
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
    return showDoFResults(problem, VectorField([], Φ.ϕ, Φ.f, [], length(Φ.f), :u3D), :uvec, name=name, visible=visible, ff=ff)
end

"""
    FEM.showStrainResults(problem, E, comp; name=..., visible=..., smooth=...)

Loads strain results into a View in gmsh. `E` is a strain field to show, `comp` is
the component of the field (:e, :ex, :ey, :ez, :exy, :eyz, :ezx),
`name` is a title to display, `visible` is a true or false value to toggle on or
off the initial visibility in gmsh and `smooth` is a true of false value to toggle
smoothing the stress field on or off. If `E` contains more than one time steps, then a 
sequence of results will be shown (eg. as an animation). This function returns
the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `E`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStrainResults(problem, E, comp; name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    t = E.t
    if E.nsteps != length(t)
        error("showStrainResults: number of time steps missmatch ($(E.nsteps) <==> $(length(t))).")
    end
    EE = gmsh.view.add(name)
    if E.A == []
        error("showStrainResults: No data")
    end
    ε = E.A
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
    FEM.showElementResults(problem, F, comp; name=..., visible=..., smooth=...)

Same as `ShowStressResults` or `showStrainResults`, depending on the type of `F` data field.

Return: `tag`

Types:
- `problem`: Problem
- `F`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showElementResults(problem, F, comp; name=comp, visible=false, smooth=true)
    if F.type == :e
        return showStrainResults(problem, F, comp; name=comp, visible=false, smooth=smooth)
    elseif F.type == :s
        return showStressResults(problem, F, comp; name=comp, visible=false, smooth=smooth)
    elseif F.type == :q2D || F.type == :q3D
        return showHeatFluxResults(problem, F, comp; name=comp, visible=false, smooth=smooth)
    else
        error("showElementResults: type is '$(F.type)'")
    end
end

"""
    FEM.showStressResults(problem, S, comp; name=..., visible=..., smooth=...)

Loads stress results into a View in gmsh. `S` is a stress field to show, `comp` is
the component of the field (:s, :sx, :sy, :sz, :sxy, :syz, :szx, :seqv),
`name` is a title to display, `visible` is a true or false value to toggle on or
off the initial visibility in gmsh and `smooth` is a true of false value to toggle
smoothing the stress field on or off. If `S` contains more than one time steps, then a 
sequence of results will be shown (eg. as an animation). This function returns
the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `S`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStressResults(problem, S, comp; name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    t = S.t
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    if S.nsteps != length(t)
        error("showStressResults: number of time steps missmatch ($(S.nsteps) <==> $(length(t))).")
    end
    SS = gmsh.view.add(name)
    if S.A == []
        error("showStressResults: No data")
    end
    σ = S.A
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
    FEM.showHeatFluxResults(problem, Q, comp; name=..., visible=..., smooth=...)

Loads heat flux results into a View in gmsh. `Q` is a heat flux field to show, `comp` is
the component of the field (:qvec, :qx, :qy, :qz, :q),
`name` is a title to display, `visible` is a true or false value to toggle on or
off the initial visibility in gmsh and `smooth` is a true of false value to toggle
smoothing the stress field on or off. If `Q` contains more than one time steps, then a 
sequence of results will be shown (eg. as an animation). This function returns
the tag of View.

Return: `tag`

Types:
- `problem`: Problem
- `S`: VectorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showHeatFluxResults(problem, S, comp; name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    t = S.t
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    if S.nsteps != length(t)
        error("showStressResults: number of time steps missmatch ($(S.nsteps) <==> $(length(t))).")
    end
    SS = gmsh.view.add(name)
    if S.A == []
        error("showHeatFluxResults: No data")
    end
    σ = S.A
    numElem = S.numElem
    for jj in 1:length(t)

        k = 1im
        if comp == :qvec
            σcomp = []
            for i in eachindex(S.numElem)
                es = length(σ[i][:,jj]) ÷ dim
                aa = zeros(3es)
                aa[1:3:3es] .= σ[i][1:dim:dim*es, jj]
                aa[2:3:3es] .= σ[i][2:dim:dim*es, jj]
                aa[2:3:3es] .= dim == 3 ? σ[i][3:dim:dim*es, jj] : zeros(es,1)
                push!(σcomp, aa)
            end
            #σcomp = [σ[i][:,jj] for i in 1:length(S.numElem)]
            nc = 3
        elseif comp == :q
            nc = 1
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                seqv = zeros(div(size(σ[i], 1), dim))
                for j in 1:(div(size(σ[i], 1), dim))
                    sx = σ[i][dim*j-(dim-1), jj]
                    sy = σ[i][dim*j-(dim-2), jj]
                    sz = dim == 3 ? σ[i][dim*j-(dim-3), jj] : 0
                    seqv[j] = √(sx^2+sy^2+sz^2)
                end
                push!(σcomp, seqv)
            end
        else
            nc = 1
            if comp == :qx
                k = 1
            elseif comp == :qy
                k = 2
            elseif comp == :qz && dim == 3
                k = 3
            else
                error("ShowHeatFluxResults: component is $comp ????")
            end
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                ss = zeros(div(size(σ[i], 1), dim))
                for j in 1:(div(size(σ[i], 1), dim))
                    ss[j] = σ[i][dim*j-(dim-k), jj]
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

