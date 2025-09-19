export Problem, material, getEigenVectors, getEigenValues
export displacementConstraint, load, elasticSupport
export temperatureConstraint, heatFlux, heatSource, heatConvection
export field, scalarField, vectorField, tensorField
export constrainedDoFs, freeDoFs
export elementsToNodes, nodesToElements, projectTo2D, expandTo3D, isNodal, isElementwise
export fieldError, resultant
export rotateNodes
export showDoFResults, showModalResults, showBucklingResults
export showStrainResults, showStressResults, showElementResults, showHeatFluxResults
export plotOnPath, showOnSurface
export openPreProcessor, openPostProcessor, setParameter, setParameters
export probe
export saveField, loadField, isSaved

"""
    Material(phName, type, E, ν, ρ, k, c, α, λ, μ, κ)

Structure containing the material type and constants.
- type: constitutive law (`:Hooke`, `:StVenantKirchhoff`, `:NeoHookeCompressible`)
- E: elastic modulus
- ν: Poisson's ratio
- ρ: mass density
- k: thermal conductivity
- c: specific heat
- α: thermal expansion coefficient
- λ: Lamé parameter
- μ: Lamé parameter
- κ: bulk modulus
`phName` is the name of the physical group where the material is used.

Types:
- `phName`: String
- `type`: Symbol
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
    type::Symbol
    E::Float64
    ν::Float64
    ρ::Float64
    k::Float64
    c::Float64
    α::Float64
    λ::Float64
    μ::Float64
    κ::Float64
    A::Float64
    Material() = new()
    Material(name) = new(name, :none, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    Material(name, type, E, ν, ρ, k, c, α, λ, μ, κ, A) = new(name, type, E, ν, ρ, k, c, α, λ, μ, κ, A)
end


"""
    Problem(materials; thickness=..., type=..., bandwidth=...)

Structure containing key data for a problem.
- Parts of the model with their material constants. Multiple materials can be provided (see `material`).
- Problem type: `:Solid`, `:PlaneStrain`, `:PlaneStress`, `:AxiSymmetric`, `:HeatConduction`, `:PlaneHeatConduction`, 
`:AxiSymmetricHeatConduction`, `:Truss`.
  For `:AxiSymmetric`, the symmetry axis is `y`, and the geometry must be drawn in the positive `x` half-plane.
- Bandwidth optimization using Gmsh's built-in reordering. Options: `:RCMK`, `:Hilbert`, `:Metis`, or `:none` (default).
- Dimension of the problem, determined by `type`.
- Material constants (vector of `Material`; see the `Material` struct).
- Plate thickness (for 2D plate problems).
- Number of nodes (`non`).
- Geometry dimension.
- Problem dimension (e.g., a 3D heat conduction problem is a 1D problem).
- In case of 2D truss displacements have to be fixed in the third direction.

Types:
- `materials`: Material
- `type`: Symbol
- `bandwidth`: Symbol
- `dim`: Integer
- `thickness`: Float64
- `non`: Integer
- `dim`: Integer
- `pdim`: Integer
"""
struct Problem
    name::String
    type::Symbol
    dim::Int64
    pdim::Int64
    material::Vector{Material}
    thickness::Float64
    non::Int64
    Problem() = new()
    Problem(name, type, dim, pdim, material, thickness, non) = new(name, type, dim, pdim, material, thickness, non)
    function Problem(mat; thickness=1, type=:Solid, bandwidth=:none)
        
        if Sys.CPU_THREADS != Threads.nthreads()
            @warn "Number of threads($(Threads.nthreads())) ≠ logical threads in CPU($(Sys.CPU_THREADS))."
        end
        
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
        elseif type == :Truss
            dim = 3
            pdim = 3
        else
            error("Problem type can be: `:Solid`, `:PlaneStress`, `:PlaneStrain`, `:AxiSymmetric`, `:PlaneHeatConduction`, `:HeatConduction`, `:AxiSymmetricHeatConduction`.
            Now problem type = $type ????")
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
                if edim != dim && (type != :Truss || edim != 1)
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

"""
    Transformation(T::SparseMatrixCSC{Float64}, non::Int64, dim::Int64)

Structure containing the transformation matrix `T` at each node in the FE mesh, the number of
nodes `non`, and the problem dimension `dim`.

Types:
- `T`: SparseMatrixCSC{Float64}
- `non`: Int64
- `dim`: Int64
"""
struct Transformation
    T::SparseMatrixCSC{Float64}
    non::Int64
    dim::Int64
end

"""
    SystemMatrix(A::SparseMatrixCSC{Float64}, model::Problem)

Structure containing the stiffness/mass/heat conduction/heat capacity/latent heat/... matrix and the associated `Problem`.

Types:
- `A`: SparseMatrixCSC{Float64}
- `model`: Problem
"""
struct SystemMatrix
    A::SparseMatrixCSC
    model::Problem
end

"""
    Base.show(io::IO, M::SystemMatrix)

Internal function to display `SystemMatrix` as a `SparseMatrixCSC{Float64}`.
"""
function Base.show(io::IO, M::SystemMatrix)
    display(M.A)
end

import Base.copy

"""
    Base.copy(A::SystemMatrix)

Internal function to copy the whole content of a `SystemMatrix`.
"""
function copy(A::SystemMatrix)
    return SystemMatrix(copy(A.A), A.model)
end

abstract type AbstractField end

"""
    ScalarField(A, a, t, numElem, nsteps, type, model)
    ScalarField(problem, dataField)

Structure containing all data of a scalar field (e.g., temperature).
- A: vector of element-wise scalar data
- a: matrix of nodal values of the scalar field
- numElem: vector of element tags
- nsteps: number of time steps stored in `A` (for animations)
- type: type of data (e.g., `:T` for temperature)
- model: associated `Problem`

Types:
- `A`: Vector{Vector{Float64}}
- `a`: Matrix{Float64}
- `t`: Vector{Float64}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: Symbol
- `model`: Problem

# Example

```julia
s(x,y,z) = 2x + 3y
fs = field("body", f=s)
S = ScalarField(problem, [fs])
```
Here `S` is defined element-wise.
"""
struct ScalarField <: AbstractField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
    model::Problem
    function ScalarField(A0, a0, t0, numElem0, nsteps0, type0, model)
        return new(A0, a0, t0, numElem0, nsteps0, type0, model)
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
                                sc1[k, 1] = f
                            end
                        end
                        push!(A, sc1)
                    end
                end
            end
        end
        a = [;;]
        t = []
        return new(A, a, t, numElem, nsteps, type, problem)
    end
end

"""
    VectorField(A, a, t, numElem, nsteps, type, model)

Structure containing the data of a vector field (e.g., displacement or heat flux).
- A: vector of element-wise vector data
- a: matrix of nodal values of the vector field
- numElem: vector of element tags
- nsteps: number of time steps stored in `A` (for animations)
- type: type of data (e.g., `:u`, `:q`)
- model: associated `Problem`

Types:
- `A`: Vector{Matrix{Float64}}
- `a`: Matrix{Float64}
- `t`: Vector{Float64}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: Symbol
- `model`: Problem
"""
struct VectorField <: AbstractField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
    model::Problem
end

"""
    TensorField(A, a, t, numElem, nsteps, type, model)

Structure containing the data of a tensor field (e.g., stress or strain).
- A: vector of element-wise tensor data
- a: matrix of nodal values of the tensor field
- numElem: vector of element tags
- nsteps: number of time steps stored in `A` (for animations)
- type: type of data (e.g., `:s`, `:e`)
- model: associated `Problem`

Types:
- `A`: Vector{Matrix{Float64}}
- `a`: Matrix{Float64}
- `t`: Vector{Float64}
- `numElem`: Vector{Integer}
- `nsteps`: Integer
- `type`: Symbol
- `model`: Problem
"""
struct TensorField <: AbstractField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
    model::Problem
end

function copy(A::AbstractField)
    T = typeof(A)
    
    a = copy(A.A)
    b = copy(A.a)
    c = copy(A.t)
    d = copy(A.numElem)
    e = copy(A.nsteps)
    f = A.type
    g = A.model
    return T(a, b, c, d, e, f, g)
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
cs = CoordinateSystem([nx, ny])
Q = rotateNodes(problem, "body", cs)
q2 = Q' * q1 # where `q1` is in Cartesian, `q2` is in Axisymmetric coordinate system and
             # `q1` is a nodal displacement vector.
S2 = Q' * S1 * Q # where `S1` is a stress field in Cartesian coordinate system while
                 # `S2` is in Axisymmetric coordinate system.

# 3D case
n1x(x, y, z) = x
n1y(x, y, z) = y
n2x(x, y, z) = -y
n2y = n1x
cs = CoordinateSystem([n1x, n1y, 0], [n2x, n2y, 0])
Q = rotateNodes(problem, "body", cs)
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
    Eigen(f, ϕ, model)

A structure containing the eigenfrequencies and eigen modes.
- f: eigenfrequencies
- ϕ: eigen modes
- model: same as `Problem`

Types:
- `f`: Matrix{Float64}
- `ϕ`: Vector{Float64}
- `model`: Problem
"""
struct Eigen
    f::Vector{Float64}
    ϕ::Matrix{Float64}
    model::Problem
end

"""
    material(name; type=:Hooke, E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5, λ=νE/(1+ν)/(1-2ν), μ=E/(1+ν)/2, κ=E/(1-2ν)/3)

Returns a structure in which `name` is the name of a physical group,
`type` is the name of the constitutive law (e.g., `:Hooke`),
`E` is the modulus of elasticity, `ν` is Poisson's ratio, and `ρ` is
the mass density. `k` is the thermal conductivity, `c` is the specific
heat, `α` is the coefficient of thermal expansion, `λ` and `μ` are the 
Lamé parameters, and `κ` is the bulk modulus.

Returns: `mat`

# Examples

```julia
mat = material("body", E=210e3, ν=0.3, ρ=7.85e-9)
```

Types:
- `mat`: Material
- `name`: String
- `type`: Symbol
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
function material(name; type=:Hooke, E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5, μ=E/(1+ν)/2, λ=2μ*ν/(1-2ν), κ=2μ*(1+ν)/(1-2ν)/3, A=1)
    if type != :Hooke &&
        type != :StVenantKirchhoff &&
        type != :NeoHookeCompressible
        error("material: type can be :Hooke, :StVenantKirchhoff or :NeoHookeCompressible. Now type is $type.")
    end
    return Material(name, type, E, ν, ρ, k, c, α, λ, μ, κ, A)
end

"""
    getEigenVectors(A::TensorField)

A function to extract the columns of a tensor field to separate vector fields.

Return: N1, N2, N3

Types:
- `A`: TensorField
- `N1`: VectorField
- `N2`: VectorField
- `N3`: VectorField

# Examples:

```julia
using LinearAlgebra
λ, Q = eigen(S)
N1, N2, N2 = getEigenVectors(Q)
```
"""
function getEigenVectors(A::TensorField)
    if length(A.A) != 0
        if A isa TensorField
            nsteps = A.nsteps
            c1 = []
            c2 = []
            c3 = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 9
                h1 = zeros(3n, nsteps)
                h2 = zeros(3n, nsteps)
                h3 = zeros(3n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        E = reshape(A.A[i][9j-8:9j, k], 3, 3)
                        h1[3j-2:3j, k] = reshape(E[:,1], 3, 1)
                        h2[3j-2:3j, k] = reshape(E[:,2], 3, 1)
                        h3[3j-2:3j, k] = reshape(E[:,3], 3, 1)
                    end
                end
                push!(c1, h1)
                push!(c2, h2)
                push!(c3, h3)
            end
            a = [;;]
            return VectorField(c1, a, A.t, A.numElem, A.nsteps, :v3D, A.model), 
            VectorField(c2, a, A.t, A.numElem, A.nsteps, :v3D, A.model), 
            VectorField(c3, a, A.t, A.numElem, A.nsteps, :v3D, A.model)
        else
            error("getEigenVectors(A::TensorField): A is not a TensorField.")
        end
    else
        error("getEigenVectors(TensorField): data at nodes is not yet implemented.")
    end
end

"""
    getEigenValues(A::VectorField)

A function to extract the elements of a vector field to separate scalar fields.

Return: λ1, λ2, λ3

Types:
- `A`: VectorField
- `λ1`: ScalarField
- `λ2`: ScalarField
- `λ3`: ScalarField

# Examples:

```julia
using LinearAlgebra
λ, Q = eigen(S)
λ1, λ2, λ2 = getEigenValues(λ)
```
"""
function getEigenValues(A::VectorField)
    if length(A.A) != 0
        if A isa VectorField
            nsteps = A.nsteps
            c1 = []
            c2 = []
            c3 = []
            for i in 1:length(A.A)
                n = length(A.A[i]) ÷ 3
                h1 = zeros(n, nsteps)
                h2 = zeros(n, nsteps)
                h3 = zeros(n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        E = reshape(A.A[i][3j-2:3j, k], 3, 1)
                        h1[j, k] = E[1,1]
                        h2[j, k] = E[2,1]
                        h3[j, k] = E[3,1]
                    end
                end
                push!(c1, h1)
                push!(c2, h2)
                push!(c3, h3)
            end
            a = [;;]
            return ScalarField(c1, a, A.t, A.numElem, A.nsteps, :scalar, A.model), 
            ScalarField(c2, a, A.t, A.numElem, A.nsteps, :scalar, A.model), 
            ScalarField(c3, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
        else
            error("getEigenValues(A::VectorField): A is not a VectorField.")
        end
    else
        error("getEigenValues(VectorField): data at nodes is not yet implemented.")
    end
end

"""
    displacementConstraint(name; ux=..., uy=..., uz=...)

Specifies displacement constraints on the physical group `name`. At least one of `ux`, `uy`,
or `uz` must be provided (depending on the problem dimension). `ux`, `uy`, or `uz` can be a
constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); displacementConstraint("support1", ux=fn)`)

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

# Examples

```julia
hc = heatConvection("outer", h=10.0, Tₐ=20.0)
```

# Examples

```julia
src = heatSource("body", h=1.0e6)
```

# Examples

```julia
q = heatFlux("out", qn=500.0)
```

# Examples

```julia
bcT = temperatureConstraint("hot_face", T=100.0)
```

# Examples

```julia
ld = load("load", fy=-1.0)
```

# Examples

```julia
bc = displacementConstraint("supp", ux=0, uy=0)
```

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
    load(name; fx=..., fy=..., fz=...)

Specifies a distributed load on the physical group `name`. At least one of `fx`, `fy`, or `fz`
must be provided (depending on the problem dimension). `fx`, `fy`, or `fz` can be a constant
or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); load("load1", fx=fn)`)

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

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
    elasticSupport(name; kx=..., ky=..., kz=...)

Specifies distributed stiffness for an elastic support on the physical group `name`.
`kx`, `ky`, or `kz` can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); elasticSupport("supp1", kx=fn)`)
Default values are 1.

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

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
    temperatureConstraint(name; T=...)

Specifies temperature constraints on the physical group `name`.
`T` can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); temperatureConstraint("surf1", T=fn)`)

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `T`: Float64 or Function
"""
function temperatureConstraint(name; T=1im)
    bc0 = name, T, 1im, 1im
    return bc0
end

"""
    heatFlux(name; qn=...)

Specifies the heat flux normal to the surface of the physical group `name`.
`qn` can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); load("flux1", qn=fn)`)

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

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
    heatSource(name; h=...)

Specifies the volumetric heat source in the physical group `name`.
`h` can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); load("source1", h=fn)`)

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

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
    heatConvection(name; h=..., Tₐ=...)

Specifies convective boundary conditions on the surface in the physical group `name`.
`h` is the heat transfer coefficient of the surrounding medium; `Tₐ` is the ambient temperature.
`Tₐ` can be either a constant or a function of `x`, `y`, and `z`.

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

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
    generateMesh(problem, surf, elemSize; approxOrder=1, algorithm=6, quadrangle=0, internalNodes=0)

Obsolete; use a Gmsh script (.geo) instead.

Returns: nothing
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
    field(name; f=..., fx=..., fy=..., fz=..., fxy=..., fyz=..., fzx=...)

Specifies the value of a scalar, vector, or tensor field on the physical group `name`.
At least one of `fx`, `fy`, or `fz` etc. must be provided (depending on the problem dimension).
Each component can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); field("surf1", fx=fn)`)

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function, ... x7}

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
f2 = field("face1", f=f1)
qq = scalarField(problem, [f2])
qqq = showDoFResults(qq, :scalar)
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
    scalarField(problem, dataField)

Defines a scalar field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`.

Return: Vector{Float64}

Types:
- `problem`: Problem
- `dataField`: Vector{Tuple{String, Float64,...}}

# Examples

```julia
f2 = field("face1", f=1)
qq = scalarField(problem, [f2])
qqq = showDoFResults(qq, :scalar)
```

Here ScalarField is defined in nodes.
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
    return ScalarField([], reshape(field, :, 1), [0], [], 1, :scalar, problem)
end

"""
    vectorField(problem, dataField; type=...)

Defines a vector field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`. `type` can be an arbitrary `Symbol`,
e.g., `:u` or `:f`.

Return: VectorField

Types:
- `problem`: Problem
- `dataField`: Vector{Tuple{String, Float64,...}}

# Examples

```julia
f1(x, y, z) = sin(x)
f2(x, y, z) = 5y
ff1 = field("face1", fx=f1, fy=f2, fz=0)
ff2 = field("face2", fx=f2, fy=f1, fz=1)
qq = vectorField(problem, [ff1, ff2])
qq0 = showDoFResults(qq, :vector)
```

Here VectorField is defined in nodes.
"""
function vectorField(problem, dataField)
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
        type = :v3D #Symbol(String(type) * "3D")
    elseif pdim == 2
        type = :v2D #Symbol(String(type) * "2D")
    else
        error("vectorField: dimension is $pdim")
    end
    return VectorField([], reshape(field, :,1), [0], [], 1, type, problem)
end

"""
    tensorField(problem, dataField; type=...)

Defines a vector field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`. `type` can be an arbitrary `Symbol`,
e.g., `:u` or `:f`.

Return: TensorField

Types:
- `problem`: Problem
- `dataField`: Vector{Tuple{String, Float64,...}}

# Examples

```julia
f1(x, y, z) = sin(x)
f2(x, y, z) = 5y
ff1 = field("face1", fx=f1, fy=f2, fz=0, fxy=1, fyz=1, fzx=f2)
ff2 = field("face2", fx=f2, fy=f1, fz=1, fxy=1, fyz=f1, fzx=1)
qq = tensorField(problem, [ff1, ff2])
qq0 = showDoFResults(qq, :tensor)
```

Here TensorField is defined in nodes.
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
    return TensorField([], reshape(field, :,1), [0], [], 1, type, problem)
end

"""
    constrainedDoFs(problem, supports)

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
    freeDoFs(problem, supports)

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
    elementsToNodes(T)

Solves the nodal results `F` from the elemental results `T`.
`T` can be ScalarField, VectorField or TensorField.

Return: `F`

Types:
- `T`: ScalarField, VectorField or TensorField
- `F`: ScalarField, VectorField or TensorField
"""
function elementsToNodes(S)
    problem = S.model
    gmsh.model.setCurrent(problem.name)
    
    if S.a != [;;]
        return S
    end
    T = typeof(S)
    type = S.type
    nsteps = S.nsteps
    numElem = S.numElem
    σ = S.A
    non = problem.non
    if S isa TensorField
        epn = 9
    elseif (S.type == :v3D || S.model.dim == 3) && S isa VectorField
        epn = 3
    elseif (S.type ==:v2D || S.model.dim == 2) && S isa VectorField
        epn = 2
    elseif S isa ScalarField #type == :scalar
        epn = 1
    else
        error("elementsToNodes: type is $type .")
    end
    s = zeros(non * epn, nsteps)
    pcs = zeros(Int64, non)
    #display("epn = $epn")
    #display("size of s = $(size(s))")
    
    for e in 1:length(numElem)
        #display("e=$e")
        elementType, nodeTags, dim, tag = gmsh.model.mesh.getElement(numElem[e])
        #display("nodeTags = $nodeTags")
        for i in 1:length(nodeTags)
            #display("size of σ[$(e)] = $(size(σ[e]))")
            s[(nodeTags[i]-1) * epn + 1: nodeTags[i] * epn, :] .+= 
            σ[e][(i-1)*epn+1:i*epn, :]
            pcs[nodeTags[i]] += 1
        end
    end
    for l in 1:non
        s[epn * (l - 1) + 1: epn * l, :] ./= pcs[l]
    end
    return T([], s, S.t, [], S.nsteps, type, problem)
end

"""
    nodesToElements(T)

Solves the element results `F` from the nodal results `T`.
`T` can be ScalarField, VectorField or TensorField.

Return: `F`

Types:
- `T`: ScalarField, VectorField or TensorField
- `F`: ScalarField, VectorField or TensorField
"""
function nodesToElements(r::Union{ScalarField,VectorField,TensorField})
    problem = r.model
    gmsh.model.setCurrent(problem.name)
    
    if r.A != []
        return r
    end
    T = typeof(r)
    if r isa ScalarField
        size = 1
    elseif r isa VectorField
        size = r.model.dim
    elseif r isa TensorField
        size = 9
    end
    type = r.type
    
    nsteps = r.nsteps
    ε = []
    numElem = Int[]
    ncoord2 = zeros(3 * problem.non)
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        #ν = problem.material[ipg].ν
        dim = 0
        if problem.dim == 3 && problem.type == :Solid
            dim = 3
            rowsOfH = 3
        elseif problem.dim == 2 && problem.type == :PlaneStress
            dim = 2
            rowsOfH = 2
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            dim = 2
            rowsOfH = 2
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
            dim = 2
            rowsOfH = 2
        else
            error("deformationGradient: dimension is $(problem.dim), problem type is $(problem.type).")
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
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
                h = reshape(fun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                nn2 = zeros(Int, size * numNodes)
                H = zeros(size * numNodes, size * numNodes)
                for k in 1:numNodes, l in 1:numNodes
                    for kk in 1:size
                        H[k*size-(size-kk), l*size-(size-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    push!(numElem, elem)
                    for k in 1:size
                        nn2[k:size:size*numNodes] = size * nnet[j, 1:numNodes] .- (size - k)
                    end
                    e = zeros(size*numNodes, nsteps)
                    for k in 1:numNodes
                        H1 = H[k*size-(size-1):k*size, 1:size*numNodes]
                        for kk in 1:nsteps
                            e0 = H1 * r.a[nn2, kk]
                            e[(k-1)*size+1:k*size, kk] = [e0[ll] for ll in 1:size]
                        end
                    end
                    push!(ε, e)
                end
            end
        end
    end
    return T(ε, [;;], r.t, numElem, nsteps, type, problem)
end

"""
    isNodal(field)

Check if a given field is defined at nodes (nodal quantity).

Nodal quantities are associated with mesh nodes, for example
displacements, nodal forces, or nodal temperatures.

# Examples
```julia
isNodal(displacement_field)   # returns true
isNodal(strain_field)         # returns false
```
"""
function isNodal(a::Union{ScalarField,VectorField,TensorField})
    if a.a != [;;] && a.A == []
        return true
    elseif a.a == [;;] && a.A != []
        return false
    else
        error("isNodal: internal error.")
    end
end

"""
isElementwise(field)

Check if a given field is defined per element (elementwise quantity).

Elementwise quantities are associated with finite elements as a whole,
for example stresses, strains, or energy densities evaluated inside elements.

Examples
```julia
isElementwise(displacement_field)   # returns false
isElementwise(strain_field)         # returns true
```
"""
function isElementwise(a::Union{ScalarField,VectorField,TensorField})
    if a.a == [;;] && a.A != []
        return true
    elseif a.a != [;;] && a.A == []
        return false
    else
        error("isElementwise: internal error.")
    end
end

"""
    expandTo3D(v2D::VectorField)

Expand a 2D vector field into 3D by adding a zero z-component.

return: VectorField

# Examples
```julia
V3D = expandTo3D(V2D)
"""
function expandTo3D(a::VectorField)
    problem = a.model
    if a.type != :v2D
        error("planeToSpace: argument must be a 2D VectorField.")
    end
    if a.a == [;;]
        b = []
        for i in 1:length(a.numElem)
            B = zeros(length(a.A[i]) ÷ 2 * 3, a.nsteps)
            B[1:3:length(B), :] .= a.A[i][1:2:length(a.A[i]), :]
            B[2:3:length(B), :] .= a.A[i][2:2:length(a.A[i]), :]
            push!(b, B)
        end
        return VectorField(b, [;;], a.t, a.numElem, a.nsteps, :v3D, a.model)
    else
        non = problem.non
        nsteps = a.nsteps
        b = zeros(3non, nsteps)
        b[1:3:3non, :] = a.a[1:2:2non, :]
        b[2:3:3non, :] = a.a[2:2:2non, :]
        return VectorField([], b, a.t, [], nsteps, :v3D, a.model)
    end
end

"""
    projectTo2D(v3D::VectorField)

Project a 3D vector field into 2D by adding a zero z-component.

return: VectorField

# Examples
```julia
V2D = expandTo3D(V3D)
"""
function projectTo2D(a::VectorField)
    problem = a.model
    if a.type != :v3D
        error("planeToSpace: argument must be a 3D VectorField.")
    end
    if a.a == [;;]
        b = []
        for i in 1:length(a.numElem)
            B = zeros(length(a.A[i]) ÷ 3 * 2, a.nsteps)
            B[1:2:length(B), :] .= a.A[i][1:3:length(a.A[i]), :]
            B[2:2:length(B), :] .= a.A[i][2:3:length(a.A[i]), :]
            push!(b, B)
        end
        return VectorField(b, [;;], a.t, a.numElem, a.nsteps, :v2D, a.model)
    else
        non = problem.non
        nsteps = a.nsteps
        b = zeros(2non, nsteps)
        b[1:2:2non, :] = a.a[1:3:3non, :]
        b[2:2:2non, :] = a.a[2:3:3non, :]
        return VectorField([], b, a.t, [], nsteps, :v2D, a.model)
    end
end

"""
    fieldError(F)

Computes the deviation of field results `F` (stresses, strains, heat flux components) at nodes
where the field has jumps. The result can be displayed with the `showDoFResults` function.

Returns: `e`

Types:
- `F`: VectorField or TensorField
- `e`: ScalarField
"""
function fieldError(S)
    problem = S.model
    gmsh.model.setCurrent(problem.name)
    
    type = S.type
    nsteps = S.nsteps
    numElem = S.numElem
    σ = S.A
    non = problem.non
    if type == :s || type == :e
        epn = 9
    elseif type == :v3D
        epn = 3
    elseif type == :v2D
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
    if type == :v3D || type == :v2D
        return VectorField([], res, S.t, [], S.nsteps, type, problem)
    elseif type == :e || type == :s
        return TensorField([], res, S.t, [], S.nsteps, type, problem)
    else
        error("fieldError: internal error, type=$type.")
    end
end

"""
    resultant(problem, field, phName; grad=false, component=:x)

Computes the resultant of `field` on the physical group `phName`.
Returns the resultant(s) in a `tuple`. The number of elements in the tuple depends on the
dimension of `problem`.
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
    if !isa(field, VectorField) && !isa(field, Matrix)
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
            s[i] += field.a[dim * nodes[j] - (dim - i)] * b
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
                    ff = isnan(f[comp0]) ? 0 : f[comp0]
                    s1 += ff * Ja * intWeights[j]
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
    rotateNodes(problem, phName, CoordSys)

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
    showDoFResults(q, comp; name=..., visible=...)

Loads nodal results into a View in Gmsh. `q` is the field to show, `comp` is
the component of the field (:vector, :uvec, :ux, :uy, :uz, :vvec, :vx, :vy, :vz,
:qvec, :qx, :qy, :qz, :T, :p, :qn, :s, :sx, :sy, :sz, :sxy, :syx, :syz,
:szy, :szx, :sxz, :e, :ex, :ey, :ez, :exy, :eyx, :eyz, :ezy, :ezx, :exz, :seqv, :scalar, :tensor),
`name` is a title to display and `visible` is a Boolean to toggle the initial visibility in Gmsh on or off.
If `q` has more columns, then a sequence of results will be shown (e.g., as an animation).
This function returns the tag of the View.

Returns: `tag`

Types:
- `q`: ScalarField, VectorField or TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showDoFResults(q, comp; name=comp, visible=false, ff = 0, factor=0)
    problem = q.model
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    t = q.t
    #display(q)
    if q.a == [;;]
        q = elementsToNodes(q)
    end
    if q.a == [;;]
        error("showDoFResults: No data")
    end
    edim = problem.dim
    if problem.type == :Truss
        edim = 1
    end
    pdim = problem.pdim
    pdim = div(size(q.a,1), problem.non)
    nodeTags = []
    ##############################################################################
    #if problem.type == :Reynolds || problem.type == :NavierStokes
    #    phName = problem.material.phName
    #    tag = getTagForPhysicalName(phName)
    #    nT, coords = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
    #    append!(nodeTags, nT)
    ##############################################################################
    #else #########################################################################
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        tag = getTagForPhysicalName(phName)
        nT, coords = gmsh.model.mesh.getNodesForPhysicalGroup(edim, tag)
        append!(nodeTags, nT)
    end
    #end #########################################################################
    
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
    elseif visible == true
        gmsh.view.option.setNumber(uvec, "Visible", 1)
    end
    gmsh.view.option.setNumber(uvec, "ShowTime", 0)
    if ff == 0 && length(t) > 1
        gmsh.view.option.setNumber(uvec, "ShowTime", 1)
    elseif ff == 1 || ff == 2
        gmsh.view.option.setNumber(uvec, "ShowTime", 6)
        gmsh.view.option.setNumber(uvec, "VectorType", 5)
    end
    gmsh.view.option.setNumber(uvec, "DisplacementFactor", factor)
    return uvec
end

"""
    showModalResults(Φ, name=..., visible=...)

Loads modal results into a View in Gmsh. `Φ` is an `Eigen` struct. `name` is a title to display and
`visible` is a Boolean to toggle the initial visibility in Gmsh on or off. Click on ▷| to change
the results. This function returns the tag of the View.

Returns: `tag`

Types:
- `Φ`: Eigen
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showModalResults(Φ::Eigen; name=:modal, visible=false, ff=1)
    return showDoFResults(VectorField([], Φ.ϕ, Φ.f, [], length(Φ.f), :v3D, Φ.model), :uvec, name=name, visible=visible, ff=ff)
end

"""
    showBucklingResults(Φ, name=..., visible=...)

Loads buckling results into a View in Gmsh. `Φ` is an `Eigen` struct. `name` is a title to display and
`visible` is a Boolean to toggle the initial visibility in Gmsh on or off. Click on ▷| to change
the results. This function returns the tag of the View.

Returns: `tag`

Types:
- `Φ`: Eigen
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showBucklingResults(Φ::Eigen; name="buckling", visible=false, ff=2)
    return showDoFResults(VectorField([], Φ.ϕ, Φ.f, [], length(Φ.f), :v3D, Φ.model), :uvec, name=name, visible=visible, ff=ff)
end

"""
    showStrainResults(E, comp; name=..., visible=..., smooth=...)

Loads strain results into a View in Gmsh. `E` is a strain field to show, `comp` is
the component of the field (:e, :ex, :ey, :ez, :exy, :eyz, :ezx),
`name` is a title to display, `visible` is a Boolean to toggle the initial visibility in Gmsh on or
off and `smooth` is a Boolean to toggle
smoothing the stress field on or off. If `E` contains more than one time steps, then a 
sequence of results will be shown (e.g., as an animation). This function returns
the tag of the View.

Returns: `tag`

Types:
- `E`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStrainResults(E, comp; name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    problem = E.model
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
    gmsh.view.option.setNumber(EE, "ShowTime", 0)
    if length(t) > 1
        gmsh.view.option.setNumber(EE, "ShowTime", 1)
    end
    #display("$comp..ok")
    return EE
end

"""
    showElementResults(F, comp; name=..., visible=..., smooth=...)

Same as `ShowStressResults` or `showStrainResults`, depending on the type of `F` data field.

Return: `tag`

Types:
- `F`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showElementResults(F, comp; name=comp, visible=false, smooth=false, factor=0)
    if F.type == :e
        return showStrainResults(F, comp, name=name, visible=visible, smooth=smooth)
    elseif F.type == :s
        return showStressResults(F, comp, name=name, visible=visible, smooth=smooth)
    elseif F isa VectorField && F.A != []
        return showHeatFluxResults(F, comp, name=name, visible=visible, smooth=smooth, factor=factor)
    elseif F isa ScalarField && F.A != []
        return showScalarResults(F, name=name, visible=visible, smooth=smooth, factor=factor)
    else
        error("showElementResults: type is '$(F.type)'")
    end
end

"""
    showStressResults(S, comp; name=..., visible=..., smooth=...)

Loads stress results into a View in Gmsh. `S` is a stress field to show, `comp` is
the component of the field (:s, :sx, :sy, :sz, :sxy, :syz, :szx, :seqv),
`name` is a title to display, `visible` is a Boolean to toggle the initial visibility in Gmsh on or
off, and `smooth` is a Boolean to toggle smoothing the stress field on or off. If `S` contains more
than one time step, a sequence of results will be shown (e.g., as an animation). This function returns
the tag of the View.

Returns: `tag`

Types:
- `S`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStressResults(S, comp; name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    problem = S.model
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
    gmsh.view.option.setNumber(SS, "ShowTime", 0)
    if length(t) > 1
        gmsh.view.option.setNumber(SS, "ShowTime", 1)
    end
    #display("$comp..ok")
    return SS
end

"""
    showHeatFluxResults(Q, comp; name=..., visible=..., smooth=...)

Loads heat flux results into a View in Gmsh. `Q` is a heat flux field to show, `comp` is
the component of the field (:qvec, :qx, :qy, :qz, :q),
`name` is a title to display, `visible` is a Boolean to toggle the initial visibility in Gmsh on or
off, and `smooth` is a Boolean to toggle smoothing on or off. If `Q` contains more than one time step,
a sequence of results will be shown (e.g., as an animation). This function returns the tag of the View.

Returns: `tag`

Types:
- `S`: VectorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showHeatFluxResults(S, comp; name=comp, visible=false, smooth=true, factor=0)
    problem = S.model
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
        if comp == :qvec || comp == :vector
            σcomp = []
            for i in eachindex(S.numElem)
                es = length(σ[i][:,jj]) ÷ dim
                aa = zeros(3es)
                aa[1:3:3es] .= σ[i][1:dim:dim*es, jj]
                aa[2:3:3es] .= σ[i][2:dim:dim*es, jj]
                aa[3:3:3es] .= dim == 3 ? σ[i][3:dim:dim*es, jj] : zeros(es,1)
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
    gmsh.view.option.setNumber(SS, "ShowTime", 0)
    if length(t) > 1
        gmsh.view.option.setNumber(SS, "ShowTime", 1)
    end
    gmsh.view.option.setNumber(SS, "DisplacementFactor", factor)
    #display("$comp..ok")
    return SS
end

function showScalarResults(S; name="scalar", visible=false, smooth=false, factor=0)
    problem = S.model
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    dim = 1
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
        
        nc = 1
        σcomp = []
        sizehint!(σcomp, length(numElem))
        for i in 1:length(S.numElem)
            sc = zeros(div(size(σ[i], 1), dim))
            for j in 1:(div(size(σ[i], 1), dim))
                sc[j] = σ[i][j, jj]
            end
            push!(σcomp, sc)
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
    gmsh.view.option.setNumber(SS, "ShowTime", 0)
    if length(t) > 1
        gmsh.view.option.setNumber(SS, "ShowTime", 1)
    end
    gmsh.view.option.setNumber(SS, "DisplacementFactor", factor)
    #display("$comp..ok")
    return SS
end

"""
    plotOnPath(problem, pathName, field; points=100, step=..., plot=..., name=..., visible=..., offsetX=..., offsetY=..., offsetZ=...)

Loads a 2D plot along a path into a View in Gmsh. `field` is the View id in
Gmsh from which the field data is imported. `pathName` is the name of a
physical group that contains a curve. The curve is divided into equal-length
segments with `points` sampling points. The field is shown at these points.
`step` is the sequence number of the displayed step. If no step is given, it shows all
available steps as an animation. If `plot` is true, an additional return parameter (a tuple of
vectors) is returned, where `x` is the horizontal axis and `y` is the vertical axis of the plot
(see the `Plots` package). `name` is the title of the graph, and
`visible` is a Boolean to toggle the initial visibility in Gmsh on or off.
This function returns the tag of the View.

Returns: `tag`

or

Returns: `tag`, `xy`

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
function plotOnPath(pathName, field; points=100, step=1im, plot=false, name="field [$field] on $pathName", visible=false, offsetX=0, offsetY=0, offsetZ=0)
    #gmsh.model.setCurrent(problem.name)
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
    gmsh.view.option.setNumber(pathView, "AdaptVisualizationGrid", 0)
    
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
    showOnSurface(field, phName; grad=false, component=:x, offsetX=0, offsetY=0, offsetZ=0, name=phName, visible=false)

Shows the values of a scalar field on a surface with physical name `phName`.
`field` is the tag of a View in Gmsh. The values of the field are calculated at the
intersection with the surface. `grad` is a Boolean to toggle the gradient of the field on or off.
`component` is the component of the gradient of `field` (:x, :y, :z) to be shown. `offsetX`, `offsetY`, `offsetZ`
are the offsets in the x, y, and z directions where the values are sampled. `name` is a title to display,
and `visible` is a Boolean to toggle the initial visibility in Gmsh on or off.

Returns: `tag`

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
function showOnSurface(phName, field; grad=false, component=:x, offsetX=0, offsetY=0, offsetZ=0, name=phName, visible=false)
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
    openPreProcessor(; openGL=...)

Launches the Gmsh preprocessor window with OpenGL disabled by default.

Returns: nothing

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
    openPostProcessor(; model=...)

Launches the Gmsh postprocessor window with the postprocessor tree opened (of `model`).

Returns: nothing

Types:
- `model`: Int64
"""
function openPostProcessor(; model=0)
    gmsh.fltk.openTreeItem(LazyString(model)*"Modules/Post-processing")
    gmsh.fltk.run()
end

"""
    setParameter(name, value)

Defines a parameter `name` and sets its value to `value`. 

Returns: nothing

Types:
- `name`: String
- `value`: Float64
"""
function setParameter(name, value)
    gmsh.parser.setNumber(name, [value])
end

"""
    setParameters(name, value)

Defines a parameter `name` and sets its value to `value`, which is a `Vector{Float64}`. 

Returns: nothing

Types:
- `name`: String
- `value`: Vector{Float64}
"""
function setParameters(name, value)
    gmsh.parser.setNumber(name, value)
end

"""
    probe(A::Union{ScalarField,VectorField,TensorField}, x::Number, y::Number, z::Number; step=Int)

Get the value of the field `A` at point coordinates `x`, `y`, `z` at time step `step`.

Returns: Float64 or Vector{Float64} or Matrix{Float64}

Types:

- `A`: ScalarField or VectorField or TensorField
- `x`: Number
- `y`: Number
- `z`: Number
- `step`: Int
"""
function probe(A::TensorField, x, y, z; step=1)
    elementTag, elementType, nodeTags, u, v, w = gmsh.model.mesh.getElementByCoordinates(x, y, z, -1, false)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementType)
    comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementType, [u, v, w], "Lagrange")
    SS = [0.0, 0, 0, 0, 0, 0, 0, 0, 0]
    if A.a == [;;]
        ind = findfirst(i -> i == elementTag, A.numElem)
        for i in range(1, 9)
            SS[i] = round(fun' * A.A[ind][i:9:9numNodes, step], digits=10)
        end
    elseif A.A == []
        for i in range(1, 9)
            SS[i] = round(fun' * A.a[9nodeTags.-(9-i), step], digits=10)
        end
    end
    return reshape(SS, 3, 3)
end

function probe(A::VectorField, x, y, z; step=1)
    elementTag, elementType, nodeTags, u, v, w = gmsh.model.mesh.getElementByCoordinates(x, y, z, -1, false)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementType)
    comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementType, [u, v, w], "Lagrange")
    dim = 0
    SS = []
    if A.model.dim == 2
        dim = 2
        SS = [0.0, 0]
    elseif A.model.dim == 3
        dim = 3
        SS = [0.0, 0, 0]
    else
        error("probe: dimension cannot be determined.")
    end
    if A.a == [;;]
        ind = findfirst(i -> i == elementTag, A.numElem)
        for i in range(1, dim)
            SS[i] = round(fun' * A.A[ind][i:dim:dim*numNodes, step], digits=10)
        end
    elseif A.A == []
        for i in range(1, dim)
            SS[i] = round(fun' * A.a[dim*nodeTags.-(dim-i), step], digits=10)
        end
    end
    return SS
end

function probe(A::ScalarField, x, y, z; step=1)
    elementTag, elementType, nodeTags, u, v, w = gmsh.model.mesh.getElementByCoordinates(x, y, z, -1, false)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementType)
    comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementType, [u, v, w], "Lagrange")
    dim = 0
    SS = 0
    if A.a == [;;]
        ind = findfirst(i -> i == elementTag, A.numElem)
        SS = round(fun' * A.A[ind][1:numNodes, step], digits=10)
    elseif A.A == []
        SS = round(fun' * A.a[nodeTags, step], digits=10)
    end
    return SS
end

"""
    probe(A::Union{ScalarField,VectorField,TensorField}, s::String; step=Int)

Get the value of the field `A` at a point given by its physical name in Gmsh at time step `step`.

Returns: Float64 or Vector{Float64} or Matrix{Float64}

Types:

- `A`: ScalarField or VectorField or TensorField
- `x`: Number
- `y`: Number
- `z`: Number
- `step`: Int
"""
function probe(A::Union{ScalarField,VectorField,TensorField}, name::String; step=1)
    phtag = getTagForPhysicalName(name)
    pttag = gmsh.model.getEntitiesForPhysicalGroup(0, phtag)[1]
    coord = gmsh.model.getValue(0, pttag, [])
    return probe(A, coord[1], coord[2], coord[3], step=step)
end

"""
    saveField(fileName::String, variable::Union{ScalarField,VectorField,TensorField,Number})

Saves `variable` of type ScalarField, VectorField, or TensorField to a file named `fileName`.
The name of the file will be complemented with the string "-LLF-Data.jld2"

Returns: nothing

Types:

- `fileName`: String
- `variable`: ScalarField, VectorField or TensorField
"""
function saveField(fileName::String, variable::Union{ScalarField,VectorField,TensorField,Number})
    @save fileName * "-LLF-Data.jld2" variable
end

"""
    loadField(fileName::String)

Loads a ScalarField, VectorField, or TensorField from a file named `fileName` (without "-LLF-Data.jld2"). 

Returns: variable

Types:

- `fileName`: String
- `variable`: ScalarField, VectorField or TensorField
"""
function loadField(fileName::String)
    @load fileName * "-LLF-Data.jld2" variable
    return variable
end

"""
    isSaved(fileName::String)

Checks whether a variable has been saved or not.

Returns: Boolean

Types:

- `fileName`: String
"""
function isSaved(fileName::String)
    return isfile(fileName * "-LLF-Data.jld2")
end
