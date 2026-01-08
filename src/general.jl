export Problem, Material, getEigenVectors, getEigenValues
export displacementConstraint, load, elasticSupport, BoundaryCondition
export temperatureConstraint, heatFlux, heatSource, heatConvection
export field, scalarField, vectorField, tensorField, ScalarField, VectorField, TensorField
export constrainedDoFs, freeDoFs, DoFs
export elementsToNodes, nodesToElements, projectTo2D, expandTo3D, isNodal, isElementwise
export fieldError, resultant, integrate, ∫, normalVector, tangentVector
export rotateNodes, CoordinateSystem
export showDoFResults, showModalResults, showBucklingResults
export showStrainResults, showStressResults, showElementResults, showHeatFluxResults
export plotOnPath, showOnSurface
export openPreProcessor, openPostProcessor, setParameter, setParameters
export probe
export saveField, loadField, isSaved
export ∂x, ∂y, ∂z, ∂t

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
    λ::Float64
    μ::Float64
    κ::Float64
    ρ::Float64
    k::Float64
    c::Float64
    α::Float64
    η::Float64
    p₀::Float64
    A::Float64

    function Material(
        phName::String;
        type = :Hooke,
        E = nothing, ν = nothing, λ = nothing, μ = nothing, κ = nothing,
        ρ = 7.85e-9,
        k = 45.0,
        c = 4.2e8,
        α = 1.2e-5,
        η = 1e-7,
        p₀ = 0.1,
        A = 1.0
    )
        type ∈ (:Hooke, :StVenantKirchhoff, :NeoHookeCompressible, :JFO) ||
            error("Invalid material type: $type")

        ec = elastic_constants(E=E, ν=ν, λ=λ, μ=μ, κ=κ)

        return new(
            phName,
            type,
            ec.E, ec.ν, ec.λ, ec.μ, ec.κ,
            ρ, k, c, α,
            η, p₀, A
        )
    end
end

function elastic_constants(; E=nothing, ν=nothing, λ=nothing, μ=nothing, κ=nothing)
    given = Dict(:E=>E, :ν=>ν, :λ=>λ, :μ=>μ, :κ=>κ)
    specified = filter(kv -> kv[2] !== nothing, given)
    n = length(specified)
    DEFAULT_POISSON = 0.3
    DEFAULT_YOUNG = 2.0e5

    if n == 0
        E = DEFAULT_YOUNG
        ν = DEFAULT_POISSON
        μ = E / (2*(1+ν))
        λ = 2μ*ν/(1-2ν)
        κ = E / (3*(1-2ν))

    elseif n == 1 && haskey(specified, :E)
        E = specified[:E]
        ν = DEFAULT_POISSON
        μ = E / (2*(1+ν))
        λ = 2μ*ν/(1-2ν)
        κ = E / (3*(1-2ν))

    elseif n == 2 && haskey(specified, :E) && haskey(specified, :ν)
        E, ν = specified[:E], specified[:ν]
        μ = E / (2*(1+ν))
        λ = 2μ*ν/(1-2ν)
        κ = E / (3*(1-2ν))

    elseif n == 2 && haskey(specified, :μ) && haskey(specified, :ν)
        μ, ν = specified[:μ], specified[:ν]
        E = 2μ*(1+ν)
        λ = 2μ*ν/(1-2ν)
        κ = E / (3*(1-2ν))

    elseif n == 2 && haskey(specified, :λ) && haskey(specified, :μ)
        λ, μ = specified[:λ], specified[:μ]
        ν = λ / (2*(λ+μ))
        E = μ*(3λ+2μ)/(λ+μ)
        κ = λ + 2μ/3

    elseif n == 2 && haskey(specified, :κ) && haskey(specified, :μ)
        κ, μ = specified[:κ], specified[:μ]
        λ = κ - 2μ/3
        ν = λ / (2*(λ+μ))
        E = 2μ*(1+ν)

    elseif n < 2
        error("Insufficient elastic parameters. Specify at least E+ν, μ+ν, λ+μ, κ+μ, or E alone.")

    else
        error("Elastic constants are overdetermined: $(keys(specified))")
    end

    return (; E, ν, λ, μ, κ)
end

#=
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
    η::Float64
    p₀::Float64
    A::Float64
    Material() = new()
    Material(name) = new(name, :none, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    Material(name, type, E, ν, ρ, k, c, α, λ, μ, κ, η, p₀, A) = new(name, type, E, ν, ρ, k, c, α, λ, μ, κ, η, p₀, A)
end
=#

mutable struct Geometry
    nameGap::String
    nameVolume::String
    dim::Int64
    tagTop::Int64
    h::Any
    dhdx::Any
    hh::Any
    p::Any
    Geometry() = new("", "", 0, 0, nothing, nothing, nothing, nothing)
    function Geometry(nameGap, nameVolume)
        dimTags = gmsh.model.getEntitiesForPhysicalName(nameGap)
        dim = dimTags[1][1]
        for i in 1:length(dimTags)
            if dim ≠ dimTags[i][1]
                error("Geometry: different surface dimensions $dim <==> $(dimTags[i][1])")
            end
        end
        #tagTop = showGapThickness(nameGap)
        #tagBottom = 0 #showGapThickness(nameGeo)
        return new(nameGap, nameVolume, dim, 0, nothing, nothing, nothing, nothing)
    end
end

"""
    Problem(materials; thickness=..., type=..., bandwidth=..., dim=..., fdim=...)

Structure containing key data for a problem.
- Parts of the model with their material constants. Multiple materials can be provided (see `material`).
- Problem type: `:Solid`, `:PlaneStrain`, `:PlaneStress`, `:AxiSymmetric`, `:HeatConduction`, `:PlaneHeatConduction`, 
`:AxiSymmetricHeatConduction`, `:Truss`, `:General`.
  For `:AxiSymmetric`, the symmetry axis is `y`, and the geometry must be drawn in the positive `x` half-plane.
- Bandwidth optimization using Gmsh's built-in reordering. Options: `:RCMK`, `:Hilbert`, `:Metis`, or `:none` (default).
- Dimension of the problem, determined by `type`.
- Material constants (vector of `Material`; see the `Material` struct).
- Plate thickness (for 2D plate problems).
- Number of nodes (`non`).
- Geometry dimension.
- Problem dimension (e.g., a 3D heat conduction problem is a 1D problem).
- In case of 2D truss displacements have to be fixed in the third direction.
- `dim` is number of dimensions in space (in case of :General)
- `fdim` is the number of unknown fields (eg. scalar->1, 2D vector->2, 3D vector->3 or more) (in case of :General)

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
    geometry::Geometry
    Problem() = new()
    Problem(name, type, dim, pdim, material, thickness, non, geometry) = new(name, type, dim, pdim, material, thickness, non, geometry)
    function Problem(mat; thickness=1.0, type=:Solid, bandwidth=:none, nameTopSurface="", nameVolume="", dim=3, fdim=dim)
        if type == :dummy
            return new("dummy", :dummy, 0, 0, mat, 0, 0)
        end
        pdim = 3
        dim0 = dim
        pdim0 = fdim

        #if Sys.CPU_THREADS != Threads.nthreads()
        #    @warn "Number of threads($(Threads.nthreads())) ≠ logical threads in CPU($(Sys.CPU_THREADS))."
        #end
        geometry = Geometry()

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
        elseif type == :Poisson3D || type == :Poisson
            dim = 3
            pdim = 1
        elseif type == :Poisson2D
            dim = 2
            pdim = 1
        elseif type == :Poisson1D
            dim = 1
            pdim = 1
        elseif type == :Reynolds
            geometry = Geometry(nameTopSurface, nameVolume)
            dim = geometry.dim
            pdim = 1
        elseif type == :General
            dim = dim0
            pdim = pdim0
        else
            error("Problem type can be: `:Solid`, `:PlaneStress`, `:PlaneStrain`, `:AxiSymmetric`, `:PlaneHeatConduction`, `:HeatConduction`, `:AxiSymmetricHeatConduction`, `:Poisson1D`, `:Poisson2D`, `:Poisson3D`, `:General`.
            Now problem type = $type ????")
        end
        if !isa(mat, Vector)
            error("Problem: materials are not arranged in a vector. Put them in [...]")
        end
        name = gmsh.model.getCurrent()
        gmsh.option.setString("General.GraphicsFontEngine", "Cairo")
        gmsh.option.setString("View.Format", "%.6g")
        
        material = mat
        
        if bandwidth ≠ :none
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
            
            #method = bandwidth == :none ? :RCMK : bandwidth
            oldTags, newTags = gmsh.model.mesh.computeRenumbering(bandwidth, elemTags)
            #permOldTags = sortperm(oldTags)
            #sortNewTags = 1:length(oldTags)
            #newTags[permOldTags] = sortNewTags
            gmsh.model.mesh.renumberNodes(oldTags, newTags)
        end
        
        nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()
        non = length(nodeTags)
        return new(name, type, dim, pdim, material, thickness, non, geometry)
    end
end

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
    show(io, M.A)
end

function Base.show(io::IO, ::MIME"text/plain", M::SystemMatrix)
    show(io, M.A)
    #println(io, "SystemMatrix ($(size(M.A,1)) × $(size(M.A,2))), nnz = $(nnz(M.A))")
end

import Base:display
function display(M::SystemMatrix)
    display(M.A)
end

abstract type AbstractField end

"""
    ScalarField(A, a, t, numElem, nsteps, type, model)
    ScalarField(problem::Problem, dataField::Vector; steps=1, tmin=0.0, tmax=tmin+(steps-1))
    ScalarField(problem::Problem, phName::String, data::Union{Number,Function};
                steps=1, tmin=0.0, tmax=tmin+(steps-1))
    ScalarField(s::ScalarField; steps=1, tmin=0.0, tmax=tmin+(steps-1), step=1)

Container for a time-dependent scalar field defined on a finite element mesh
(e.g. temperature, pressure, potential).

A `ScalarField` can store the field either
- **element-wise** (values at the element nodes, stored in `A`), or
- **nodally** (values at global mesh nodes, stored in `a`).

Time dependence is handled by storing multiple time steps for each spatial
degree of freedom.

---

## Stored data

- `A` : Vector of matrices holding **element-wise** values  
  (`A[e][k,i]` = value at local node `k` of element `e` at time step `i`)
- `a` : Matrix of **nodal** values  
  (`a[n,i]` = value at node `n` at time step `i`)
- `t` : Vector of time instants corresponding to the stored time steps
- `numElem` : Vector of element tags associated with `A`
- `nsteps` : Number of stored time steps
- `type` : Symbol identifying the physical meaning of the field
- `model` : Associated `Problem`

At a given time, either `A` or `a` is typically populated, depending on whether
the field is element-wise or nodal.

---

## Constructors

### Low-level constructor
```julia
ScalarField(A, a, t, numElem, nsteps, type, model)
```

Directly constructs a `ScalarField` from preallocated data arrays.

---

### From spatial field definitions

```julia
ScalarField(problem, dataField; steps, tmin, tmax)
```

Constructs an **element-wise** scalar field from a vector of field definitions
(e.g. as returned by `field(...)`).
The scalar value may be:

* a constant,
* a spatial function `f(x,y,z)`,
* or a space–time function `f(x,y,z,t)`.

Values are evaluated at the element nodes for each requested time step.

---

### From a physical group

```julia
ScalarField(problem, phName, data; steps, tmin, tmax)
```

Convenience constructor for defining a scalar field on a single physical group.
Equivalent to calling `field(phName, f=data)` internally.

---

### From an existing ScalarField

```julia
ScalarField(s; steps, tmin, tmax, step)
```

Creates a new `ScalarField` by **replicating a selected time step** of an
existing field.
This is useful for initializing a time-dependent field from a static solution.

* `step` selects the time index of `s` to be replicated.
* The new field contains `steps` identical time slices.

---

## Notes

* Time steps do not need to be uniformly spaced.
* No assumptions are made about governing equations; this is a pure data container.
* Spatial interpolation and projections (e.g. nodal ↔ element-wise) are handled
  by separate utility functions.

---

## Example

```julia
s(x,y,z) = 2x + 3y
fs = field("body", f=s)

S1 = ScalarField(problem, [fs])
S2 = ScalarField(problem, "body", s)
```

Here `S1` and `S2` define equivalent element-wise scalar fields.
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
    function ScalarField(problem::Problem, dataField::Vector; steps=1, tmin=0.0, tmax=tmin+(steps-1))
        gmsh.model.setCurrent(problem.name)
        nodeTags, coords, _ = gmsh.model.mesh.getNodes()
        type = :scalar
        nsteps = steps
        t = range(start=tmin, stop=tmax, length=steps)
        A = Vector{Matrix{Float64}}()
        numElem = Int[]

        @inbounds for data in dataField
            #phName, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = data
            phName = data.phName
            f = data.f
            dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

            for (edim, etag) in dimTags
                elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
                @inbounds for (i, et) in enumerate(elemTypes)
                    _, _, _, numNodes, _, _ = gmsh.model.mesh.getElementProperties(et)
                    for (j, elem) in enumerate(elemTags[i])
                        sc1 = zeros(numNodes, nsteps)
                        push!(numElem, elem)
                        if f !== nothing
                            for it in 1:steps
                                for k in 1:numNodes
                                    nodeTag = elemNodeTags[i][(j-1)*numNodes + k]
                                    if f isa Function
                                        x = coords[3*nodeTag - 2]
                                        y = coords[3*nodeTag - 1]
                                        z = coords[3*nodeTag]
                                        if applicable(f,1,2,3,4)
                                            sc1[k, it] = f(x, y, z, t[it])
                                        elseif applicable(f,1,2,3)
                                            sc1[k, it] = f(x, y, z)
                                        else
                                            error("ScalarField: function must have 3 or 4 arguments.")
                                        end
                                    else
                                        sc1[k, it] = f
                                    end
                                end
                            end
                        end
                        push!(A, sc1)
                    end
                end
            end
        end

        return ScalarField(A, [;;], t, numElem, nsteps, type, problem)
    end
    function ScalarField(problem::Problem, phName::String, data::Union{Number,Function}; steps=1, tmin=0.0, tmax=tmin+(steps-1))
        f = field(phName, f=data)
        return ScalarField(problem, [f], steps=steps, tmin=tmin, tmax=tmax)
    end
    function ScalarField(s::ScalarField; steps=1, tmin=0.0, tmax=tmin+(steps-1), step=1)
        if step > s.nsteps || step < 1
            error("ScalarField: `step` is out of range - 0 < $step < $(s.nsteps)")
        end
        A = Vector{Matrix{Float64}}(undef, length(s.numElem))
        s1 = nodesToElements(s)
        for j = 1:length(s1.numElem)
            A[j] = zeros(size(s1.A[j], 1), steps)
            for i in 1:steps
                A[j][:, i] = s1.A[j][:, step]
            end
        end
        t = range(start=tmin, stop=tmax, length=steps)
        ScalarField(A, [;;], t, s1.numElem, steps, s1.type, s1.model)
    end
end

"""
    VectorField(A, a, t, numElem, nsteps, type, model)
    VectorField(problem::Problem, dataField::Vector)
    VectorField(problem::Problem, phName::String, data::Vector)
    VectorField(problem::Problem, phName::String, func::Function)
    VectorField(comps::Vector{ScalarField})

Container for a (possibly time-dependent) vector field defined on a finite
element mesh (e.g. displacement, velocity, heat flux).

A `VectorField` stores a 3D vector field either
- **element-wise**, with values given at element nodes (`A`), or
- **nodally**, with values given at global mesh nodes (`a`).

Time dependence is handled by storing multiple time steps for each spatial
degree of freedom, analogously to `ScalarField`.

---

## Stored data

- `A` : Vector of matrices holding **element-wise** vector values  
  (`A[e][3k-2:3k, i]` = vector value at local node `k` of element `e`
  at time step `i`)
- `a` : Matrix of **nodal** vector values  
  (layout follows the same interleaved component ordering)
- `t` : Vector of time instants corresponding to the stored time steps
- `numElem` : Vector of element tags associated with `A`
- `nsteps` : Number of stored time steps
- `type` : Symbol identifying the physical meaning of the vector field
  (e.g. `:v3D`)
- `model` : Associated `Problem`

Vector components are stored in **interleaved form**
(`x₁,y₁,z₁,x₂,y₂,z₂,…`) for each element or node.

---

## Constructors

### Low-level constructor
```julia
VectorField(A, a, t, numElem, nsteps, type, model)
```

Directly constructs a `VectorField` from preallocated data arrays.

---

### From element-wise field definitions

```julia
VectorField(problem, dataField::Vector)
```

Constructs an **element-wise** vector field from a vector of field definitions
(e.g. as returned by `field(...)`).
Each component may be specified as:

* a constant,
* a spatial function `f(x,y,z)`.

Values are evaluated at the element nodes.

---

### From a physical group and component data

```julia
VectorField(problem, phName, data::Vector)
```

Convenience constructor for defining a vector field on a single physical group,
where `data = [fx, fy, fz]` specifies the three components (constants or
functions).

---

### From a vector-valued function

```julia
VectorField(problem, phName, func::Function)
```

Constructs an element-wise vector field by evaluating a function
`func(x,y,z) -> (vx,vy,vz)` at the element nodes.

---

### From scalar components

```julia
VectorField(comps::Vector{ScalarField})
```

Assembles a 3D vector field from exactly three compatible `ScalarField`s.

Requirements:

* exactly three components,
* same `Problem`,
* identical element numbering,
* identical time discretization.

Each scalar field provides one vector component.

---

## Notes

* Time steps do not need to be uniformly spaced.
* No governing equations are implied; this is a pure data container.
* Spatial operations (gradient, divergence, projections, etc.) are handled
  by separate utility functions.
* Element-wise storage is the primary representation; nodal storage may be
  empty depending on construction.

---

## Example

```julia
vx(x,y,z) = x + y
vy(x,y,z) = z

fv = field("body", fx=vx, fy=vy, fz=3)
V1 = VectorField(problem, [fv])

V2 = VectorField(problem, "body", [vx, vy, 3])
```

Here `V1` and `V2` define equivalent element-wise vector fields.
"""
struct VectorField <: AbstractField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
    model::Problem
    function VectorField(A0, a0, t0, numElem0, nsteps0, type0, model)
        return new(A0, a0, t0, numElem0, nsteps0, type0, model)
    end
    function VectorField(comps::Vector{ScalarField})
        s1 = nodesToElements(comps[1])
        s2 = nodesToElements(comps[2])
        s3 = nodesToElements(comps[3])
        comps1 = [s1, s2, s3]
        n = length(comps1)
        @assert n == 3 "VectorField must have exactly 3 scalar components."

        # check same problem
        prob = comps1[1].model
        for s in comps1
            @assert s.model === prob "All scalar fields must belong to the same problem."
        end

        # check element ordering
        elems = comps1[1].numElem
        nsteps = comps1[1].nsteps
        for s in comps1
            @assert s.numElem == elems "ScalarField element numbering mismatch."
            @assert s.nsteps == nsteps "ScalarField time steps mismatch."
        end

        #@disp comps1[1]
        #@disp s1

        # assemble elementwise vector field
        ε0 = Vector{Matrix{Float64}}(undef, length(elems))
        #ε = zeros(Matrix{Float64}, length(elems))
        for i in 1:length(elems)
            # comps1[j].a[i] is a (numNodes × nsteps) matrix
            #blocks = [comps1[j].A[i] for j in 1:3]
            block = zeros(size(comps1[1].A[i],1) * 3, size(comps1[1].A[i],2))
            for l = 1:3
                block[l:3:end, :] = comps1[l].A[i]
            end
            #ε0[i] = vcat(blocks...)   # 3*numNodes × nsteps
            ε0[i] = block
        end
        #@disp ε0

        a = [;;]
        return new(copy(ε0), a, comps1[1].t, elems, comps1[1].nsteps, :v3D, prob)
    end
    function VectorField(problem::Problem, dataField)
        if !isa(dataField, Vector)
            error("VectorField: dataField are not arranged in a vector. Put them in [...]")
        end
        gmsh.model.setCurrent(problem.name)
        
        type = :v3D
        nsteps = 1
        A = []
        numElem = Int[]
        ff = 0
        pdim = 1
        for i in 1:length(dataField)
            #phName, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = dataField[i]
            phName = dataField[i].phName
            fx = dataField[i].fx
            fy = dataField[i].fy
            fz = dataField[i].fz
            dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
            for idm in 1:length(dimTags)
                dimTag = dimTags[idm]
                edim = dimTag[1]
                etag = dimTag[2]
                elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
                for i in 1:length(elemTypes)
                    et = elemTypes[i]
                    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                    for j in 1:length(elemTags[i])
                        vec1 = zeros(3numNodes, nsteps)
                        elem = elemTags[i][j]
                        push!(numElem, elem)
                        for k in 1:numNodes
                            nodeTag = elemNodeTags[i][(j-1)*numNodes+k]
                            coord, parametricCoord, dim, tag = gmsh.model.mesh.getNode(nodeTag)
                            x = coord[1]
                            y = coord[2]
                            z = coord[3]
                            if fx !== nothing
                                if fx isa Function
                                    ff = fx(x, y, z)
                                else
                                    ff = fx
                                end
                                vec1[3k-2, 1] = ff
                            end
                            if fy !== nothing
                                if fy isa Function
                                    ff = fy(x, y, z)
                                else
                                    ff = fy
                                end
                                vec1[3k-1, 1] = ff
                            end
                            if fz !== nothing
                                if fz isa Function
                                    ff = fz(x, y, z)
                                else
                                    ff = fz
                                end
                                vec1[3k, 1] = ff
                            end
                        end
                        push!(A, vec1)
                    end
                end
            end
        end
        a = [;;]
        t = [0.0]
        return new(A, a, t, numElem, nsteps, type, problem)
    end
    function VectorField(problem::Problem, phName::String, data::Vector)
        if size(data) == (3,)
            f = field(phName, fx=data[1], fy=data[2], fz=data[3])
            return VectorField(problem, [f])
        else
            error("VectorField: size of data is $(size(data)).")
        end
    end
    function VectorField(problem::Problem, phName::String, func::Function)
        gmsh.model.setCurrent(problem.name)
        
        type = :v3D
        nsteps = 1
        A = []
        numElem = Int[]
        ff = 0
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                for j in 1:length(elemTags[i])
                    vec1 = zeros(3numNodes, nsteps)
                    elem = elemTags[i][j]
                    push!(numElem, elem)
                    for k in 1:numNodes
                        nodeTag = elemNodeTags[i][(j-1)*numNodes+k]
                        coord, parametricCoord, dim, tag = gmsh.model.mesh.getNode(nodeTag)
                        x = coord[1]
                        y = coord[2]
                        z = coord[3]
                        vec1[3k-2:3k] = func(x, y, z)
                    end
                    push!(A, vec1)
                end
            end
        end
        a = [;;]
        t = [0.0]
        return new(A, a, t, numElem, nsteps, type, problem)
    end
end

"""
    TensorField(A, a, t, numElem, nsteps, type, model)
    TensorField(problem::Problem, dataField::Vector)
    TensorField(problem::Problem, phName::String, data::Matrix)
    TensorField(comps::Matrix{ScalarField})

Container for a (possibly time-dependent) second-order tensor field defined on a
finite element mesh (e.g. stress, strain, conductivity tensor).

A `TensorField` stores a full 3×3 tensor at each spatial location, either
- **element-wise**, with values given at element nodes (`A`), or
- **nodally**, with values given at global mesh nodes (`a`).

Time dependence is handled by storing multiple time steps for each tensor
component.

---

## Stored data

- `A` : Vector of matrices holding **element-wise** tensor values  
  (`A[e][9k-8:9k, i]` = tensor components at local node `k` of element `e`
  at time step `i`)
- `a` : Matrix of **nodal** tensor values  
  (same component ordering as element-wise storage)
- `t` : Vector of time instants corresponding to the stored time steps
- `numElem` : Vector of element tags associated with `A`
- `nsteps` : Number of stored time steps
- `type` : Symbol identifying the physical meaning of the tensor field
  (e.g. `:e`, `:stress`, `:tensor`)
- `model` : Associated `Problem`

Tensor components are stored in **row-major, interleaved form**
for each node:
```

(xx, yx, zx,
xy, yy, zy,
xz, yz, zz)

```
repeated for all nodes and time steps.

---

## Constructors

### Low-level constructor
```julia
TensorField(A, a, t, numElem, nsteps, type, model)
```

Directly constructs a `TensorField` from preallocated data arrays.

---

### From element-wise field definitions

```julia
TensorField(problem, dataField::Vector)
```

Constructs an **element-wise** tensor field from a vector of field definitions
(e.g. as returned by `field(...)`).

Each tensor component may be specified as:

* a constant,
* a spatial function `f(x,y,z)`.

Values are evaluated at the element nodes.

---

### From a constant or functional tensor

```julia
TensorField(problem, phName, data::Matrix)
```

Convenience constructor for defining a tensor field on a single physical group,
where `data` is a `3×3` matrix whose entries are constants or functions
`f(x,y,z)`.

---

### From scalar components

```julia
TensorField(comps::Matrix{ScalarField})
```

Assembles a 3×3 tensor field from a `3×3` matrix of compatible `ScalarField`s.

Requirements:

* exactly a 3×3 matrix of scalar fields,
* identical `Problem`,
* identical element numbering,
* identical time discretization.

Each scalar field provides one tensor component.

---

## Notes

* Time steps do not need to be uniformly spaced.
* No symmetry is assumed; all 9 tensor components are stored explicitly.
* This is a pure data container; no constitutive or kinematic assumptions
  are implied.
* Spatial operations (e.g. divergence, invariants, projections) are handled
  by separate utility functions.
* Element-wise storage is the primary representation; nodal storage may be
  empty depending on construction.

---

## Example

```julia
tx(x,y,z)  = x + y
txy(x,y,z) = z

ft = field("body", fx=tx, fxy=txy, fz=3)
T1 = TensorField(problem, [ft])

T2 = TensorField(problem, "body", [1 0 0;
                                  0 2 0;
                                  0 0 3])
```

Here `T1` and `T2` define equivalent element-wise tensor fields.
"""
struct TensorField <: AbstractField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
    model::Problem
    function TensorField(A0, a0, t0, numElem0, nsteps0, type0, model)
        return new(A0, a0, t0, numElem0, nsteps0, type0, model)
    end
    function TensorField(problem, dataField)
        if !isa(dataField, Vector)
            error("TensorField: dataField are not arranged in a vector. Put them in [...]")
        end
        gmsh.model.setCurrent(problem.name)
        
        type = :e
        nsteps = 1
        A = []
        numElem = Int[]
        ff = 0
        pdim = 1
        for i in 1:length(dataField)
            #phName, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = dataField[i]
            phName = dataField[i].phName
            fx = dataField[i].fx
            fy = dataField[i].fy
            fz = dataField[i].fz
            fxy = dataField[i].fxy
            fyz = dataField[i].fyz
            fzx = dataField[i].fzx
            fyx = dataField[i].fyx
            fzy = dataField[i].fzy
            fxz = dataField[i].fxz
            dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
            for idm in 1:length(dimTags)
                dimTag = dimTags[idm]
                edim = dimTag[1]
                etag = dimTag[2]
                elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
                for i in 1:length(elemTypes)
                    et = elemTypes[i]
                    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                    for j in 1:length(elemTags[i])
                        ten1 = zeros(9numNodes, nsteps)
                        elem = elemTags[i][j]
                        push!(numElem, elem)
                        for k in 1:numNodes
                            nodeTag = elemNodeTags[i][(j-1)*numNodes+k]
                            coord, parametricCoord, dim, tag = gmsh.model.mesh.getNode(nodeTag)
                            x = coord[1]
                            y = coord[2]
                            z = coord[3]
                            if fx !== nothing
                                if fx isa Function
                                    ff = fx(x, y, z)
                                else
                                    ff = fx
                                end
                                ten1[9k-8, 1] = ff
                            end
                            if fy !== nothing
                                if fy isa Function
                                    ff = fy(x, y, z)
                                else
                                    ff = fy
                                end
                                ten1[9k-4, 1] = ff
                            end
                            if fz !== nothing
                                if fz isa Function
                                    ff = fz(x, y, z)
                                else
                                    ff = fz
                                end
                                ten1[9k, 1] = ff
                            end
                            if fxy !== nothing
                                if fxy isa Function
                                    ff = fxy(x, y, z)
                                else
                                    ff = fxy
                                end
                                ten1[9k-5, 1] = ff
                            end
                            if fyz !== nothing
                                if fyz isa Function
                                    ff = fyz(x, y, z)
                                else
                                    ff = fyz
                                end
                                ten1[9k-1, 1] = ff
                            end
                            if fzx !== nothing
                                if fzx isa Function
                                    ff = fzx(x, y, z)
                                else
                                    ff = fzx
                                end
                                ten1[9k-6, 1] = ff
                            end
                            if fyx !== nothing
                                if fyx isa Function
                                    ff = fyx(x, y, z)
                                else
                                    ff = fyx
                                end
                                ten1[9k-7, 1] = ff
                            end
                            if fzy !== nothing
                                if fzy isa Function
                                    ff = fzy(x, y, z)
                                else
                                    ff = fzy
                                end
                                ten1[9k-3, 1] = ff
                            end
                            if fxz !== nothing
                                if fxz isa Function
                                    ff = fxz(x, y, z)
                                else
                                    ff = fxz
                                end
                                ten1[9k-2, 1] = ff
                            end
                        end
                        push!(A, ten1)
                    end
                end
            end
        end
        a = [;;]
        t = [0.0]
        return new(A, a, t, numElem, nsteps, type, problem)
    end
    function TensorField(problem::Problem, phName::String, data::Matrix)
        if size(data) == (3,3)
            g = [x isa Function ? x : ((_,_,_)->x) for x in data]
            f = field(phName, fx=g[1], fy=g[5], fz=g[9], fxy=g[4], fyz=g[8], fzx=g[3], fyx=g[2], fzy=g[6], fxz=g[7])
            #f = field(phName, fx=data[1], fy=data[5], fz=data[9], fxy=data[2], fyz=data[6], fzx=data[3])
            return TensorField(problem, [f])
        else
            error("TensorField: size of data is $(size(data)).")
        end
    end
    function TensorField(comps::Matrix{ScalarField})
        @assert size(comps) == (3,3) "TensorField requires a 3×3 matrix of ScalarFields."
   
        # --- 1) Konvertálás elementwise alakra (nodesToElements) ---
        comps2 = Matrix{ScalarField}(undef, 3, 3)
        @inbounds for i in 1:3, j in 1:3
            comps2[i,j] = nodesToElements(comps[i,j])
        end
   
        # --- 2) Ellenőrzés: minden komponens ugyanahhoz a problémához tartozik ---
        prob   = comps2[1,1].model
        numElem = comps2[1,1].numElem
        nsteps = comps2[1,1].nsteps
        t      = comps2[1,1].t
   
        @inbounds for i in 1:3, j in 1:3
            s = comps2[i,j]
            @assert s.model === prob        "TensorField: ScalarFields must share same Problem."
            @assert s.numElem == numElem    "TensorField: element list mismatch."
            @assert s.nsteps  == nsteps     "TensorField: time step mismatch."
        end
   
        # --- 3) Elemwise tensor összeállítása ---
        A = Vector{Matrix{Float64}}(undef, length(numElem))
   
        for e in 1:length(numElem)
            #blocks = Matrix{Float64}[]
            block = zeros(size(comps2[1].A[e],1) * 9, size(comps2[1].A[e],2))
            # 9 komponens 3×3 mátrixból sorfolytonosan
            @inbounds for i in 1:3, j in 1:3
                #push!(blocks, comps2[i,j].A[e])
                block[(j-1)*3+i:9:end, :] = comps2[i,j].A[e]
            end
            #A[e] = vcat(blocks...)   # méret: (9*numNodes × nsteps)
            A[e] = block
        end
   
        # --- 4) TensorField visszaadása ---
        return TensorField(A, [;;], t, numElem, nsteps, :tensor, prob)
    end
end

struct BoundaryCondition
    phName::String

    # Dirichlet-type (primary variables)
    ux::Union{Nothing, Number, Function, ScalarField}
    uy::Union{Nothing, Number, Function, ScalarField}
    uz::Union{Nothing, Number, Function, ScalarField}
    p ::Union{Nothing, Number, Function, ScalarField}   # pressure / scalar field
    T ::Union{Nothing, Number, Function, ScalarField}   # temperature

    # Neumann-type (fluxes / forces)
    f::Union{Nothing, Number, Function, ScalarField}
    fx::Union{Nothing, Number, Function, ScalarField}
    fy::Union{Nothing, Number, Function, ScalarField}
    fz::Union{Nothing, Number, Function, ScalarField}
    fxy::Union{Nothing, Number, Function, ScalarField}
    fyz::Union{Nothing, Number, Function, ScalarField}
    fzx::Union{Nothing, Number, Function, ScalarField}
    fyx::Union{Nothing, Number, Function, ScalarField}
    fzy::Union{Nothing, Number, Function, ScalarField}
    fxz::Union{Nothing, Number, Function, ScalarField}
    #qx ::Union{Nothing, Number, Function}   # heat flux / generic scalar flux
    #qy ::Union{Nothing, Number, Function}   # heat flux / generic scalar flux
    #qz ::Union{Nothing, Number, Function}   # heat flux / generic scalar flux
    kx ::Union{Nothing, Number, Function, ScalarField}
    ky ::Union{Nothing, Number, Function, ScalarField}
    kz ::Union{Nothing, Number, Function, ScalarField}
    q ::Union{Nothing, Number, Function, ScalarField}   # heat flux / generic scalar flux
    qn ::Union{Nothing, Number, Function, ScalarField}   # heat flux / generic scalar flux
    h ::Union{Nothing, Number, Function, ScalarField}

    # Stress / traction components
    sx::Union{Nothing, Number, Function, ScalarField}
    sy::Union{Nothing, Number, Function, ScalarField}
    sz::Union{Nothing, Number, Function, ScalarField}
    sxy::Union{Nothing, Number, Function, ScalarField}
    syz::Union{Nothing, Number, Function, ScalarField}
    szx::Union{Nothing, Number, Function, ScalarField}

    function BoundaryCondition(phName::String; kwargs...)
        fields = fieldnames(BoundaryCondition)
        values = map(fields) do f
            f == :phName ? phName : get(kwargs, f, nothing)
        end
        return new(values...)
    end
end

function BC_specified_fields(bc::BoundaryCondition)
    return filter(f -> getfield(bc, f) !== nothing,
                  fieldnames(BoundaryCondition))
end

using SparseArrays

"""
    DoFs(f::AbstractField)

Return the vector of degrees of freedom (DOFs) associated with a field.

This corresponds to the nodal representation used in linear algebra
operations (assembly, solving, post-processing).
"""
@inline DoFs(f::AbstractField) = f.a

import Base:copy
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

import Base:show
"""
    Base.show(io::IO, ::MIME"text/plain", M::Union{ScalarField,VectorField,TensorField})

Plain-text display in REPL and notebooks.
"""
function Base.show(io::IO, ::MIME"text/plain", M::Union{ScalarField,VectorField,TensorField})
    # fejléc külön sorban (szebb notebookban)
    println(io, nameof(typeof(M)))
    
    io2 = IOContext(io, :compact => false)
    
    if isNodal(M)
        Base.show(io2, M.a)
    else
        Base.show(io2, M.A)
    end
end

import Base:display
function display(M::Union{ScalarField,VectorField,TensorField})
    if isNodal(M)
        display("$(nameof(typeof(M))) $(M.a)")
    else
        display("$(nameof(typeof(M))) $(M.A)")
    end
end

"""
    ∂(r::ScalarField, dir::Int) -> ScalarField

Compute the spatial derivative of a scalar field in the global `x`, `y`, or `z`
direction. The input field may be nodal or elementwise; nodal fields are
automatically converted to elementwise form using `nodesToElements(r)`.

# Arguments
- `r::ScalarField`: scalar field to differentiate.
- `dir::Int`: spatial direction:
    * `1` → ∂/∂x  
    * `2` → ∂/∂y  
    * `3` → ∂/∂z

# Returns
A new `ScalarField` in **elementwise representation**, containing  
the gradient component ∂r/∂(dir) at all nodes of each element.

# Notes
- Works for 1D, 2D, and 3D meshes.
- Uses the element Jacobian and shape-function derivatives in global coordinates.
- The global problem mesh is accessed via `r.model`.
"""
function ∂(r::ScalarField, dir::Int)
    #@info "∂: nodal"
    problem = r.model
    gmsh.model.setCurrent(problem.name)

    @assert dir ∈ 1:3 "d(r,dir): dir must be 1 (x), 2 (y), or 3 (z)."

    # ha nodális tér, konvertáljuk elemisre
    #r = rr.a == [;;] ? elementsToNodes(rr) : rr
    nsteps = r.nsteps

    # globális node koordináták gyors elérése
    nodeTags, ncoord, _ = gmsh.model.mesh.getNodes()
    non = problem.non

    ncoord2 = zeros(3 * non)
    @inbounds begin
        ncoord2[nodeTags .* 3 .- 2] = ncoord[1:3:end]
        ncoord2[nodeTags .* 3 .- 1] = ncoord[2:3:end]
        ncoord2[nodeTags .* 3 .- 0] = ncoord[3:3:end]
    end

    ε = Vector{Matrix{Float64}}()    # elemi eredmények
    numElem = Int[]                  # elem azonosítók

    # anyagcsoportok végigjárása
    for ipg in 0:length(problem.material)
        phName = ""
        if ipg == 0
            if r.model.geometry.nameVolume ≠ ""
                phName = r.model.geometry.nameVolume
            else
                continue
            end
        else
            phName = problem.material[ipg].phName
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        @inbounds for (edim, etag) in dimTags
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            @inbounds for it in 1:length(elemTypes)
                et = Int(elemTypes[it])

                # elem adatok cache-ből
                dim_et, numNodes, nodeCoord = _get_props_cached(et)
                ∇h_all, h_all = _get_basis_cached(et, nodeCoord, numNodes)

                # munkatömbök
                invJac = Matrix{Float64}(undef, 3, 3*numNodes)
                ∂h    = Matrix{Float64}(undef, 1, numNodes*numNodes)  # ∂/∂dir
                nn2   = Vector{Int}(undef, numNodes)
                nnet  = Matrix{Int}(undef, length(elemTags[it]), numNodes)

                dim_et = 3
                # elem ciklus
                for j in 1:length(elemTags[it])
                    elem = elemTags[it][j]

                    # node tags → gyors index
                    @inbounds for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[it][(j-1)*numNodes + k]
                    end

                    # Jacobian
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)

                    # invJac minden lokális csomópontra
                    for k in 1:numNodes
                        Jk = @view Jac[1:dim_et, (k-1)*3+1 : k*3]
                        invJk = pinv(Matrix(Jk'))       # stabil

                        fill!(@view(invJac[1:3, 3*k-2:3*k]), 0.0)
                        @views invJac[1:dim_et, 3*k-2:3*k-3+dim_et] .= invJk
                    end

                    # --- ∂h = ∂/∂dir (shape-derivált) ---
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numNodes, l in 1:numNodes
                        g = @view ∇h_all[(l-1)*3+1 : (l-1)*3+dim_et, k]
                        Jslice = @view invJac[:, 3*k-2 : 3*k-2+dim_et-1]

                        dphi = (Jslice * g)[dir]   # dir=1,2,3 → x,y,z irány
                        ∂h[1, (k-1)*numNodes + l] = dphi
                    end

                    # globális scalar DoF indexek
                    @inbounds for k in 1:numNodes
                        nn2[k] = nnet[j, k]        # scalar → csak 1 szabadságfok
                    end

                    # elem eredmény (numNodes × nsteps)
                    e = Matrix{Float64}(undef, numNodes, nsteps)

                    @inbounds for k in 1:numNodes
                        idx = (k-1)*numNodes
                        for kk in 1:nsteps
                            s = 0.0
                            for l in 1:numNodes
                                s += ∂h[1, idx + l] * r.a[nn2[l], kk]
                            end
                            e[k, kk] = s
                        end
                    end

                    push!(ε, e)
                    push!(numElem, elem)
                end
            end
        end
    end

    return ScalarField(ε, [;;], r.t, numElem, nsteps, :scalar, problem)
end

"""
    ∂e(r::ScalarField, dir::Int) -> ScalarField

Compute the spatial derivative of an **elementwise** scalar field purely at the
element level, without nodal conversion. This avoids mixing contributions from
neighboring elements and is suitable for discontinuous fields or post-processing
of elementwise quantities.

# Arguments
- `r::ScalarField`: must already be elementwise (`isElementwise(r) == true`).
- `dir::Int`: spatial direction (1 → x, 2 → y, 3 → z).

# Returns
An elementwise `ScalarField` containing the derivative ∂r/∂(dir) at the local
nodes of each element.

# Notes
- Uses the element-local Jacobian per integration point.
- Safer for interfaces or discontinuities than nodal differentiation.
- Only elementwise → elementwise transformations are performed.
"""
function ∂e(rr::ScalarField, dir::Int)
    #@info "∂e: elementwise scalar derivative"
    problem = rr.model
    gmsh.model.setCurrent(problem.name)

    dim = problem.dim              # 1, 2 vagy 3
    @assert 1 ≤ dir ≤ dim "∂e(rr,dir): dir must be between 1 and problem.dim."

    # elementwise scalar field kell
    @assert rr.A != [] "∂e: rr must be elementwise ScalarField (A ≠ [])."

    nsteps = rr.nsteps

    # elemTag → index in rr.A / rr.numElem
    elem_index = Dict{Int,Int}()
    @inbounds for (i, e) in enumerate(rr.numElem)
        elem_index[e] = i
    end

    ε       = Vector{Matrix{Float64}}()  # elementwise results
    numElem = Int[]                      # elemTags az eredményhez

    # --- iterate physical groups (ugyanúgy, mint a többi függvényben) ---
    for ipg in 0:length(problem.material)
        phName =
            ipg == 0 ?
                (rr.model.geometry.nameVolume ≠ "" ? rr.model.geometry.nameVolume : "") :
                problem.material[ipg].phName

        if phName == ""
            continue
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        @inbounds for (edim, etag) in dimTags
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            @inbounds for it in 1:length(elemTypes)
                et = Int(elemTypes[it])

                # cached FE data
                dim_et, numNodes, nodeCoord = _get_props_cached(et)
                ∇h_all, h_all = _get_basis_cached(et, nodeCoord, numNodes)

                # csak olyan elemtípussal foglalkozzunk, ahol az rr-ben van elem
                # (elem_index úgyis szűr, itt nem muszáj, de nem árt)

                # invJac: dim × (dim_et * numNodes) – minden csomóponthoz egy dim×dim_et blokk
                invJac = Matrix{Float64}(undef, dim, dim_et * numNodes)
                ∂h     = Matrix{Float64}(undef, 1, numNodes * numNodes)   # csak a kiválasztott dir komponens

                # elem ciklus
                @inbounds for j in 1:length(elemTags[it])
                    elem = elemTags[it][j]

                    # csak azok az elemek, amelyekhez rr-ben tényleg van adat
                    eidx = get(elem_index, elem, 0)
                    if eidx == 0
                        continue
                    end

                    # Jacobianok az elem csomópontjain
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)   # gmsh mindig 3 sorral jön, de csak 1..dim kell

                    # invJac blokkok minden lokális csomópontra
                    @inbounds for k in 1:numNodes
                        # Jk: dim × dim_et
                        Jk = @view Jac[1:dim, (k-1)*3+1 : (k-1)*3+dim_et]
                        invJk = pinv(Matrix(Jk'))   # dim_et × dim  →  dim × dim_et

                        @views invJac[:, (k-1)*dim_et+1 : k*dim_et] .= invJk
                    end

                    # ∂h = ∂φ_l/∂x_dir globális koordinátákban
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numNodes, l in 1:numNodes
                        # lokális grad (dim_et komponens)
                        g = @view ∇h_all[(l-1)*3+1 : (l-1)*3+dim_et, k]

                        # dim × dim_et blokk a k-adik csomópontnál
                        Jslice = @view invJac[:, (k-1)*dim_et+1 : k*dim_et]

                        # globális grad (dim komponens) → ebből választjuk a dir-ediket
                        dphi_vec = Jslice * g
                        dphi = dphi_vec[dir]

                        ∂h[1, (k-1)*numNodes + l] = dphi
                    end

                    # csomóponti értékek ezen az elemen (numNodes × nsteps)
                    vals = rr.A[eidx]
                    @assert size(vals, 1) == numNodes "∂e: numNodes mismatch with rr.A."

                    # elem eredmény: numNodes × nsteps
                    e = Matrix{Float64}(undef, numNodes, nsteps)

                    # ∑_l ∂φ_l * value_l
                    @inbounds for k in 1:numNodes
                        idx = (k-1)*numNodes
                        for kk in 1:nsteps
                            s = 0.0
                            for l in 1:numNodes
                                s += ∂h[1, idx + l] * vals[l, kk]
                            end
                            e[k, kk] = s
                        end
                    end

                    push!(ε, e)
                    push!(numElem, elem)
                end
            end
        end
    end

    return ScalarField(ε, [;;], rr.t, numElem, nsteps, :scalar, problem)
end

"""
    ∂x(r::ScalarField)
    ∂y(r::ScalarField)
    ∂z(r::ScalarField)

Compute directional derivatives of a scalar field `r` with respect to the
global coordinates x, y, or z.

These functions are convenience wrappers that automatically choose between  
the nodal differentiation (`∂`) and the elementwise differentiation (`∂e`)
depending on the representation of the input field.

# Dispatch logic
- If `isNodal(r)` is `true`:  
    - `∂x(r)` → `∂(r, 1)`  
    - `∂y(r)` → `∂(r, 2)`  
    - `∂z(r)` → `∂(r, 3)`  
- Otherwise (elementwise field):  
    - `∂x(r)` → `∂e(r, 1)`  
    - `∂y(r)` → `∂e(r, 2)`  
    - `∂z(r)` → `∂e(r, 3)`

# Coordinate conventions
- In **3D**: `(x, y, z)` are the global Cartesian coordinates.
- In **2D (plane stress/strain)**:  
  - `∂x` and `∂y` operate in the plane.  
  - `∂z` is mathematically allowed and returns zero for planar meshes.

# Returns
A `ScalarField` storing the derivative in **elementwise** representation.

# Examples

## 1) Derivative of a nodal field (automatic nodal→element conversion)

```julia
s = ScalarField(prob, "vol", (x,y,z) -> x^2 + y)

dxs = ∂x(s)      # computes 2x
dys = ∂y(s)      # computes 1
dzs = ∂z(s)      # computes 0
```

## 2) Derivative of an elementwise field

```julia
se = elementsToNodes(s) |> r -> ∂e(r, 1)  # manually
dxs = ∂x(se)                              # same result, auto-detected
```

## 3) Works on 2D meshes as well

```julia
s = ScalarField(prob2D, "surf", (x,y,z) -> x - y)

dxs = ∂x(s)      # = 1 everywhere
dys = ∂y(s)      # = -1 everywhere
dzs = ∂z(s)      # = zero field (2D mesh)
```

## 4) Using derivatives to build a gradient vector field

```julia
s = ScalarField(prob, "vol", (x,y,z) -> sin(x)*y)

gx = ∂x(s)       # y*cos(x)
gy = ∂y(s)       # sin(x)
gz = ∂z(s)       # 0

grad_s = VectorField([gx, gy, gz])
```

# Notes

* The output is always **elementwise**, enabling consistent post-processing.
* Direction indices follow FEM convention:
  1 → x, 2 → y, 3 → z.
"""
∂x(r::ScalarField) = isNodal(r) ? ∂(r, 1) : ∂e(r, 1)
∂y(r::ScalarField) = isNodal(r) ? ∂(r, 2) : ∂e(r, 2)
∂z(r::ScalarField) = isNodal(r) ? ∂(r, 3) : ∂e(r, 3)

"""
    ∂t(s::Union{ScalarField,VectorField,TensorField})

Compute the time derivative of a time-dependent field using second-order
finite differences in time.

The field values are assumed to be stored at discrete time instants `s.t`
and are differentiated independently at each spatial degree of freedom.
The result has the same type and layout (nodal or elementwise) as the input
field.

### Time discretization
- For interior time steps, a second-order central difference is used:
```

∂s/∂t(tᵢ) ≈ (s(tᵢ₊₁) − s(tᵢ₋₁)) / (tᵢ₊₁ − tᵢ₋₁)

```
- At the first and last time steps, the derivative is obtained by
second-order extrapolation to preserve overall accuracy.
- If only two time steps are present, a first-order difference is used
and assigned to both time levels.
- If `s.nsteps == 1`, a zero field is returned.

The time grid `s.t` is not required to be uniform.

### Notes
- The operation is purely algebraic in time and does not involve any
spatial operators or time-integration schemes.
- The derivative is computed in-place on newly allocated storage and
does not modify the original field.
- For elementwise fields, the time derivative is computed independently
for each element.

### Returns
A new field of the same concrete type as `s` (`ScalarField`, `VectorField`,
or `TensorField`) containing the time derivative.
"""

function ∂t(s::Union{ScalarField,VectorField,TensorField})
    if s.nsteps == 1
        return 0s
    end
    T = typeof(s)
    if isNodal(s)
        a = similar(s.a)
        time_derivative!(a, s.t, s.a)
        return T([], a, s.t, [], s.nsteps, s.type, s.model)
    elseif isElementwise(s)
        A = [similar(s.A[e]) for e in eachindex(s.A)]
        for e in eachindex(s.A)
            time_derivative!(A[e], s.t, s.A[e])
        end
        return T(A, [;;], s.t, s.numElem, s.nsteps, s.type, s.model)
    else
        error("∂t: internal error")
    end
end

function time_derivative!(a, t, src)
    nsteps = length(t)
    if nsteps == 2
        @inbounds begin
            a[:,1] .= (src[:,2] - src[:,1]) / (t[2] - t[1])
            a[:,2] .= a[:,1]
        end
    else
        @inbounds begin
            for i in 2:nsteps-1
                a[:,i] .= (src[:,i+1] - src[:,i-1]) / (t[i+1] - t[i-1])
            end
            a[:,1]   .= 2*(src[:,2]-src[:,1])/(t[2]-t[1]) - a[:,2]
            a[:,end] .= 2*(src[:,end]-src[:,end-1])/(t[end]-t[end-1]) - a[:,end-1]
        end
    end
end

"""
    v[k] -> ScalarField

Extract the k-th component (1,2,3) of a `VectorField` as a `ScalarField`.

# Arguments
- `k::Int`: component index (1 → x, 2 → y, 3 → z).

# Returns
A `ScalarField` in the same representation as the parent field (elementwise or nodal).

# Errors
Raises an error for any non-scalar indexing:
- `v[1:2]`    → ERROR  
- `v[[1,3]]`  → ERROR  
- `v[:]`      → ERROR

Only single-component extraction `v[k]` is allowed.
"""
function Base.getindex(v::VectorField, k::Int)
    @assert 1 ≤ k ≤ 3 "VectorField has exactly 3 components. Use v[1], v[2], v[3]."

    # ha elementwise mező
    if v.A != []
        numElem = length(v.A)
        numNodes = size(v.A[1], 1) ÷ 3   # 3 komponens miatt
        nsteps = v.nsteps

        Aout = Vector{Matrix{Float64}}(undef, numElem)

        @inbounds for i in 1:numElem
            Ai = v.A[i]
            rows = k:3:3*numNodes
            Aout[i] = Ai[rows, :]
        end

        return ScalarField(Aout, [;;], v.t, v.numElem, v.nsteps, :scalar, v.model)
    end

    # ha nodal mező
    if v.a != [;;]
        numNodes = size(v.a, 1) ÷ 3
        rows = k:3:3*numNodes
        aout = v.a[rows, :]

        return ScalarField([], aout, v.t, v.numElem, v.nsteps, :scalar, v.model)
    end

    error("VectorField: no data found in A or a.")
end

"""
    T[i,j] -> ScalarField

Extract a single tensor component (i,j) from a 3×3 `TensorField`.

# Arguments
- `i::Int`, `j::Int`: tensor indices in 1:3.

# Returns
A `ScalarField` corresponding to component Tᵢⱼ.

# Notes
Block extraction is done elementwise or nodally depending on the representation
of `T`.

# Errors
Raises an error if i or j is outside 1:3.
"""
function Base.getindex(T::TensorField, i::Int, j::Int)
    @assert 1 ≤ i ≤ 3 && 1 ≤ j ≤ 3 "TensorField indexing must be T[i,j], with i,j ∈ 1:3"

    block = (j-1)*3 + i   # which tensor component
    # --- elementwise ---
    if T.A != []
        numElem  = length(T.A)
        numNodes = size(T.A[1], 1) ÷ 9     # 9 tensor components
        Aout = Vector{Matrix{Float64}}(undef, numElem)

        @inbounds for e in 1:numElem
            Ae = T.A[e]
            rows = block:9:9*numNodes
            Aout[e] = Ae[rows, :]
        end

        return ScalarField(Aout, [;;], T.t, T.numElem, T.nsteps, :scalar, T.model)
    end

    # --- nodal ---
    if T.a != [;;]
        numNodes = size(T.a, 1) ÷ 9
        rows = block:9:9*numNodes
        aout = T.a[rows, :]

        return ScalarField([], aout, T.t, T.numElem, T.nsteps, :scalar, T.model)
    end

    error("TensorField: no data in A or a.")
end

"""
    T[i,:] -> VectorField

Extract the i-th row of a 3×3 `TensorField` as a 3-component `VectorField`.

Equivalent to: [T[i,1], T[i,2], T[i,3]]

# Arguments
- `i::Int`: row index (1 ≤ i ≤ 3).

# Returns
A `VectorField` of size 3.

# Notes
Colon indexing is required; `T[i]` is not allowed.
"""
function Base.getindex(T::TensorField, i::Int, ::Colon)
    @assert 1 ≤ i ≤ 3 "TensorField row index must be 1..3"

    comps = ScalarField[]
    for j in 1:3
        push!(comps, T[i,j])   # already returns ScalarField
    end
    return VectorField(comps)
end

"""
    T[:,j] -> VectorField

Extract the j-th column of a 3×3 `TensorField` as a 3-component `VectorField`.

Equivalent to: [T[1,j], T[2,j], T[3,j]]

# Arguments
- `j::Int`: column index (1 ≤ j ≤ 3).

# Returns
A `VectorField` of size 3.

# Notes
Colon indexing is required; `T[j]` is not allowed.
"""
function Base.getindex(T::TensorField, ::Colon, j::Int)
    @assert 1 ≤ j ≤ 3 "TensorField column index must be 1..3"

    comps = ScalarField[]
    for i in 1:3
        push!(comps, T[i,j])
    end
    return VectorField(comps)
end

"""
Invalid index pattern for `TensorField`.

Allowed:
- `T[i,j]`    → ScalarField
- `T[i,:]`    → VectorField (row)
- `T[:,j]`    → VectorField (column)

Not allowed:
- `T[i]`
- `T[:]`
- `T[1:2,1]`
- `T[:, :]`
- `T[[1,3],2]`
"""
function Base.getindex(T::TensorField, I...)
    error("""
    Invalid TensorField indexing.

    Allowed:
      T[i,j]    → ScalarField
      T[i,:]    → VectorField (row)
      T[:,j]    → VectorField (column)

    Everything else is invalid.
    """)
end

"""
    v[k] = s

Assign a scalar component into a `VectorField`.

This operation overwrites the k-th component (k = 1, 2, 3) of the vector field.

## Requirements
- `v` must be a `VectorField`
- `s` must be a `ScalarField`
- both fields must belong to the **same problem**
- both must have the **same number of time steps**
- storage layout must match:
    - If `v` is **elementwise** (`v.A ≠ []`), then `s` must also be elementwise
      and must have the same number of elements (`numElem`).
    - If `v` is **nodal** (`v.a ≠ [;;]`), then `s` must also be nodal and must have
      the same nodal size.

Only the rows belonging to component `k` are modified.  
Other components remain unchanged.

## Example
```julia
v = VectorField([s1, s2, s3])

v[1] = s4       # overwrite first component
v[3] = 2s1      # scale and assign into third component
```

## Errors

An informative error is thrown if:

* k ∉ 1:3
* the fields belong to different problems
* storage layouts do not match (nodal vs. elementwise)
* nodal or element counts differ
"""
function Base.setindex!(v::VectorField, s::ScalarField, k::Int)
    @assert 1 ≤ k ≤ 3 "VectorField has exactly 3 components."

    # check problem consistency
    @assert s.model === v.model "ScalarField and VectorField belong to different problems."
    @assert s.nsteps == v.nsteps "Time step mismatch in VectorField assignment."

    # elementwise case
    if v.A != []
        @assert s.A != [] "Cannot assign nodal field into elementwise VectorField."
        @assert length(s.A) == length(v.A) "numElem mismatch in VectorField assignment."

        numNodes = size(v.A[1], 1) ÷ 3
        @assert size(s.A[1], 1) == numNodes "Node count mismatch in VectorField assignment."

        @inbounds for e in 1:length(v.A)
            # overwrite only the block belonging to component k
            rows = k:3:3*numNodes
            v.A[e][rows, :] .= s.A[e]
        end
        return v
    end

    # nodal case
    if v.a != [;;]
        @assert s.a != [;;] "Cannot assign elementwise field into nodal VectorField."

        numNodes = size(v.a, 1) ÷ 3
        @assert size(s.a, 1) == numNodes "Node count mismatch in VectorField assignment."

        rows = k:3:3*numNodes
        v.a[rows, :] .= s.a
        return v
    end

    error("VectorField: no data to assign into (empty A and a).")
end

"""
    T[i, j] = s
    T[i, :] = v
    T[:, j] = v

Assign a scalar or vector component into a `TensorField`.

## Supported forms

### 1. Assign scalar component
`T[i, j] = s`

- Inserts the scalar field `s` into the tensor component `(i, j)`  
  (with i,j ∈ 1:3).
- Only the rows of the `(i,j)` block are overwritten.

### 2. Assign row vector
`T[i, :] = v`

- `v` must be a `VectorField` with 3 components.
- Replaces the entire i-th tensor row with the components of `v`.

### 3. Assign column vector
`T[:, j] = v`

- `v` must be a `VectorField`.
- Replaces the entire j-th tensor column with the components of `v`.

## Requirements
- Assigned field(s) must belong to the **same problem** as `T`.
- The number of time steps must match.
- Storage mode must match:
    - Elementwise tensor → only elementwise scalar/vector may be assigned.
    - Nodal tensor → only nodal scalar/vector may be assigned.
- Element counts or nodal sizes must be identical.

## Examples
```julia
# scalar insertion
T[1,1] = s1
T[2,3] = s2

# row insertion
T[1,:] = v1    # v1 is VectorField([s11, s12, s13])

# column insertion
T[:,3] = v2
```

## Errors

Errors are thrown when:

* indices are out of bounds
* shapes do not match (nodal vs elementwise mismatch)
* component sizes differ
* attempting unsupported forms like `T[1] = ...` or `T[:, :] = ...`
"""
function Base.setindex!(T::TensorField, s::ScalarField, i::Int, j::Int)
    @assert 1 ≤ i ≤ 3 && 1 ≤ j ≤ 3 "TensorField indexing must be T[i,j] with i,j ∈ 1:3."

    @assert s.model === T.model "TensorField and ScalarField belong to different problems."
    @assert s.nsteps == T.nsteps "Time step mismatch in TensorField assignment."

    block = (j-1)*3 + i  # component index

    # elementwise
    if T.A != []
        @assert s.A != [] "Cannot assign nodal field into elementwise TensorField."
        @assert length(s.A) == length(T.A) "numElem mismatch in TensorField assignment."

        numNodes = size(T.A[1], 1) ÷ 9
        @assert size(s.A[1], 1) == numNodes "Node count mismatch in TensorField assignment."

        rows = block : 9 : 9*numNodes

        @inbounds for e in 1:length(T.A)
            T.A[e][rows, :] .= s.A[e]
        end
        return T
    end

    # nodal
    if T.a != [;;]
        @assert s.a != [;;] "Cannot assign elementwise field into nodal TensorField."

        numNodes = size(T.a, 1) ÷ 9
        @assert size(s.a, 1) == numNodes

        rows = block : 9 : 9*numNodes
        T.a[rows, :] .= s.a
        return T
    end

    error("TensorField: empty storage (A and a both empty).")
end

function Base.setindex!(T::TensorField, v::VectorField, ::Colon, j::Int)
    @assert 1 ≤ j ≤ 3 "TensorField column index must be 1..3."
    @assert v.model === T.model "Model mismatch in TensorField assignment."
    @assert v.nsteps == T.nsteps "Time step mismatch."

    # expand column assignment
    for i in 1:3
        T[i, j] = v[i]   # reuses the ScalarField assignment above
    end

    return T
end

function Base.setindex!(T::TensorField, v::VectorField, i::Int, ::Colon)
    @assert 1 ≤ i ≤ 3 "TensorField row index must be 1..3."
    @assert v.model === T.model "Model mismatch in TensorField assignment."
    @assert v.nsteps == T.nsteps "Time step mismatch."

    # expand row assignment
    for j in 1:3
        T[i, j] = v[j]   # reuses ScalarField setter
    end

    return T
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
    getTagForPhysicalName(name)
                            
Returns `tag` of elements of physical group `name`.
                            
Returns: `tag`
                            
Types:
- `name`: String
- `tag`: Integer
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
    getDimForPhysicalName(name)
                            
Returns `dim` of elements of physical group `name`.
                            
Returns: `dim`
                            
Types:
- `name`: String
- `dim`: Integer
"""
function getDimForPhysicalName(name)
    dimTags = gmsh.model.getPhysicalGroups(-1)
    i = 1
    while gmsh.model.getPhysicalName(dimTags[i][1], dimTags[i][2]) != name
        i += 1
        if i > length(dimTags)
            error("Physical name '$name' does not exist.")
        end
    end
    return dimTags[i][1]
end

#=
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
function material(name; type=:Hooke, E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5, μ=E/(1+ν)/2, λ=2μ*ν/(1-2ν), κ=2μ*(1+ν)/(1-2ν)/3, η=1e-7, p₀=0.1, A=1)
    if type != :Hooke &&
        type != :StVenantKirchhoff &&
        type != :NeoHookeCompressible &&
        type ≠ :JFO
        error("material: type can be :Hooke, :StVenantKirchhoff or :NeoHookeCompressible. Now type is $type.")
    end
    return Material(name, type, E, ν, ρ, k, c, α, λ, μ, κ, η, p₀, A)
end
=#

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
function displacementConstraint(name; ux=nothing, uy=nothing, uz=nothing)
    ux === nothing && uy === nothing && uz === nothing &&
        error("displacementConstraint: at least one of ux, uy, uz must be specified.")
    return BoundaryCondition(name, ux=ux, uy=uy, uz=uz)
end
#=
function displacementConstraint(name; ux=1im, uy=1im, uz=1im)
    bc0 = name, ux, uy, uz
    return bc0
end
=#

"""
    load(name; fx=..., fy=..., fz=...)

Specifies a distributed load on the physical group `name`. At least one of `fx`, `fy`, or `fz`
must be provided (depending on the problem dimension). `fx`, `fy`, or `fz` can be a constant
or a function of `x`, `y`, and `z` or `ScalarField`.
(e.g., `fn(x,y,z) = 5*(5-x); load("load1", fx=fn)`)

Returns: Tuple{String, Float64 or Function, Float64 or Function, Float64 or Function}

Types:
- `name`: String
- `fx`: Float64 or Function
- `fy`: Float64 or Function
- `fz`: Float64 or Function
"""
function load(name; fx=nothing, fy=nothing, fz=nothing, p=nothing)
    fx === nothing && fy === nothing && fz === nothing && p === nothing &&
        error("load: at least one of fx, fy, fz, p must be specified.")
    if p !== nothing
        if p isa ScalarField
            pp = elementsToNodes(p)
            return BoundaryCondition(name, p=pp)
        end
        return BoundaryCondition(name, p=p)
    end
    if fx !== nothing
        ffx = fx isa ScalarField ? elementsToNodes(fx) : fx
    else
        ffx = fx
    end
    if fy !== nothing
        ffy = fy isa ScalarField ? elementsToNodes(fy) : fy
    else
        ffy = fy
    end
    if fz !== nothing
        ffz = fz isa ScalarField ? elementsToNodes(fz) : fz
    else
        ffz = fz
    end
    return BoundaryCondition(name, fx=ffx, fy=ffy, fz=ffz)
end
#=
function load(name; fx::Union{Number, Function, ScalarField}=0, fy::Union{Number, Function, ScalarField}=0, fz::Union{Number, Function, ScalarField}=0, p::Union{Number, Function, ScalarField}=0)
    if fx isa ScalarField
        if isElementwise(fx)
            fx = elementsToNodes(fx)
        end
    end
    if fy isa ScalarField
        if isElementwise(fy)
            fy = elementsToNodes(fy)
        end
    end
    if fz isa ScalarField
        if isElementwise(fz)
            fz = elementsToNodes(fz)
        end
    end
    if p ≠ 0
        if p isa ScalarField
            if isElementwise(p)
                p = elementsToNodes(p)
            end
        end
        fx = p
        fy = 1im
    end
    ld0 = name, fx, fy, fz
    return ld0
end
=#

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
function elasticSupport(name; kx=nothing, ky=nothing, kz=nothing)
    kx === nothing && ky === nothing && kz === nothing &&
        error("elasticSupport: at least one of kx, ky, kz, p must be specified.")
    return BoundaryCondition(name, kx=kx, ky=ky, kz=kz)
end
#=
function elasticSupport(name; kx=0, ky=0, kz=0)
    es0 = name, kx, ky, kz
    return es0
end
=#

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
function temperatureConstraint(name; T=nothing)
    T === nothing &&
        error("temperatureConstraint: T must be specified.")
    return BoundaryCondition(name, T=T)
end
#=
function temperatureConstraint(name; T=1im)
    bc0 = name, T, 1im, 1im
    return bc0
end
=#

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
function heatFlux(name; qn=nothing)
    qn === nothing &&
        error("heatFlux: qn must be specified.")
    return BoundaryCondition(name, qn=qn)
    #=
    p1 =0
    p2 =0
    qn0 = -qn
    fl0 = name, qn0, p1, p2
    return fl0
    =#
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
function heatSource(name; h=nothing)
    h === nothing &&
        error("heatSource: h must be specified.")
    return BoundaryCondition(name, h=h)
    #=
    p1 =0
    p2 =0
    h0 = -h
    sr0 = name, h0, p1, p2
    return sr0
    =#
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
function heatConvection(name; h=nothing, T∞=nothing)
    h === nothing && T∞ === nothing
        error("heatSource: h and T∞ must be specified.")
    return BoundaryCondition(name, h=h, T=T∞)
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
function field(name; f=nothing, fx=nothing, fy=nothing, fz=nothing, fxy=nothing, fyx=nothing, fyz=nothing, fzy=nothing, fzx=nothing, fxz=nothing)
    fyx = fyx !== nothing ? fyx : fxy
    fzy = fzy !== nothing ? fzy : fyz
    fxz = fxz !== nothing ? fxz : fzx
    return BoundaryCondition(name, f=f, fx=fx, fy=fy, fz=fz, fxy=fxy, fyz=fyz, fzx=fzx, fyx=fyx, fzy=fzy, fxz=fxz)
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
qqq = showDoFResults(qq)

qq2 = scalarField(problem, "body", 2.0)
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
        #name, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = dataField[i]
        name = dataField[i].phName
        f = dataField[i].f
        if f === nothing
            error("scalarField: f is not defined")
        end
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(f, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if f !== nothing
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
    
function scalarField(problem::Problem, phName::String, data::Union{Number,Function})
    f = field(phName, f=data)
    return scalarField(problem, [f])
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
qq0 = showDoFResults(qq)

qq2 = vectorField(problem, "body", [1, 2, ff2])
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
        #name, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = dataField[i]
        name = dataField[i].phName
        fx = dataField[i].fx
        fy = dataField[i].fy
        fz = dataField[i].fz
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(fx, Function) || isa(fy, Function) || isa(fz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if fx !== nothing
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
        if fy !== nothing
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
        if fz !== nothing && pdim == 3
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

function vectorField(problem::Problem, phName::String, data::Vector)
    if size(data) == (3,)
        f = field(phName, fx=data[1], fy=data[2], fz=data[3])
        return vectorField(problem, [f])
    else
        error("vectorField: size of data is $(size(data)).")
    end
end

vectorField(sm::Vector{ScalarField}) = elementsToNodes(VectorField([sm[1], sm[2], sm[3]]))
tensorField(sm) = elementsToNodes(TensorField(sm))

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
qq0 = showDoFResults(qq)

qq2 = tensorField(problem, "body", [1 0 0; 0 2 0; 0 0 f1])
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
        #name, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = dataField[i]
        name = dataField[i].phName
        fx = dataField[i].fx
        fy = dataField[i].fy
        fz = dataField[i].fz
        fxy = dataField[i].fxy
        fyz = dataField[i].fyz
        fzx = dataField[i].fzx
        fyx = dataField[i].fyx
        fzy = dataField[i].fzy
        fxz = dataField[i].fxz
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(fx, Function) || isa(fy, Function) || isa(fz, Function) || isa(fxy, Function) || isa(fyz, Function) || isa(fzx, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if fx !== nothing
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
        if fy !== nothing
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
        if fz !== nothing
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= pdim
            if isa(fz, Function)
                ffz = fz.(xx, yy, zz)
                field[nodeTagsZ,:] .= ffz
            else
                field[nodeTagsZ,:] .= fz
            end
        end
        if fxy !== nothing
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
        if fyx !== nothing
            nodeTagsYX = copy(nodeTags)
            nodeTagsYX *= pdim
            nodeTagsYX .-= (pdim - 2)
            if isa(fyx, Function)
                ffyx = fyx.(xx, yy, zz)
                field[nodeTagsYX,:] .= ffyx
            else
                field[nodeTagsYX,:] .= fyx
            end
        end
        if fyz !== nothing
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
        if fzy !== nothing
            nodeTagsZY = copy(nodeTags)
            nodeTagsZY *= pdim
            nodeTagsZY .-= (pdim - 6)
            if isa(fzy, Function)
                ffzy = fzy.(xx, yy, zz)
                field[nodeTagsZY,:] .= ffzy
            else
                field[nodeTagsZY,:] .= fzy
            end
        end
        if fzx !== nothing
            nodeTagsZX = copy(nodeTags)
            nodeTagsZX *= pdim
            nodeTagsZX .-= (pdim - 3)
            if isa(fzx, Function)
                ffzx = fzx.(xx, yy, zz)
                field[nodeTagsZX,:] .= ffzx
            else
                field[nodeTagsZX,:] .= fzx
            end
        end
        if fxz !== nothing
            nodeTagsXZ = copy(nodeTags)
            nodeTagsXZ *= pdim
            nodeTagsXZ .-= (pdim - 7)
            if isa(fxz, Function)
                ffxz = fxz.(xx, yy, zz)
                field[nodeTagsXZ,:] .= ffxz
            else
                field[nodeTagsXZ,:] .= fxz
            end
        end
    end
    return TensorField([], reshape(field, :,1), [0], [], 1, type, problem)
end

function tensorField(problem::Problem, phName::String, data::Matrix)
    if size(data) == (3,3)
        f = field(phName, fx=data[1], fy=data[5], fz=data[9], fxy=data[4], fyz=data[8], fzx=data[3], fyx=data[2], fzy=data[6], fxz=data[7])
        return tensorField(problem, [f])
    else
        error("tensorField: size of data is $(size(data)).")
    end
end

"""
    normalVector(problem::Problem, phName::String) -> VectorField

Compute outward unit normal vectors for all nodes belonging to a surface-type
physical group in a 3D Gmsh model.

For curve-type physical groups, the function returns normal curvature vectors:
the direction corresponds to the unit normal within the curve’s osculating
plane, and the vector magnitude equals the local curvature of the curve.

# Arguments
- `problem::Problem`: A `Problem` structure containing the current Gmsh model 
  (name, dimension, number of nodes, etc.).
- `phName::String`: The name of a physical surface group in Gmsh for which the normal
  vectors are computed.

# Description
The function sets the current model, queries the elements and nodes that belong to the 
given physical surface group, and evaluates the gradients of the Lagrange basis functions 
to compute local tangent vectors of the surface. Normal vectors are obtained as cross 
products of these tangents and normalized to unit length.

Each node belonging to the physical surface group is assigned its corresponding 
3D unit normal vector.

# Returns
- `VectorField`: A `VectorField` structure that stores the nodal normal vectors 
  on the given physical surface.

# Errors
- Throws an error if the provided physical group is not of surface type (`dim != 2`).

# Example
```julia
using LowLevelFEM

# Load a 3D geometry and mesh it in Gmsh
gmsh.initialize()
gmsh.open("box_with_surface.msh")

# Define a 3D model on the volume physical group "body"
mat = material("body")
prob = Problem([mat])

# Compute nodal normals on a physical surface named "leftWall"
nv = normalVector(problem, "leftWall")

# Show the normal vectors on the model
showDoFResults(nv)
openPostProcessor()
```
"""
function normalVector(problem::Problem, phName::String)
    gmsh.model.setCurrent(problem.name)

    non = problem.non
    ncoord2 = zeros(3 * non)
    numElem = Int[]
    normalvectors = Vector{Matrix{Float64}}()

    # Fizikai csoport lekérdezése (pl. felület)
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

    for dimTag in dimTags
        edim, etag = dimTag
        if edim == 1
            normalVector1D(problem, phName)
        elseif edim ≠ 2
            error("normalVector: physical group is niether a curve (1D) nor a surface (2D).")
        end

        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
        nodeTags, ncoord, _ = gmsh.model.mesh.getNodes(edim, etag, true, false)

        # Globális csomópontkoordináták betöltése
        ncoord2[nodeTags*3 .- 2] = ncoord[1:3:end]
        ncoord2[nodeTags*3 .- 1] = ncoord[2:3:end]
        ncoord2[nodeTags*3 .- 0] = ncoord[3:3:end]

        # --- Elemciklus ---
        for i in 1:length(elemTypes)
            et = elemTypes[i]
            elementName, eldim, order, numNodes, localNodeCoord, _ =
                gmsh.model.mesh.getElementProperties(et)

            # Lokális koordináták átalakítása 3D formába
            localNodeCoord2 = zeros(length(localNodeCoord) ÷ 2 * 3)
            localNodeCoord2[1:3:end] .= localNodeCoord[1:2:end]
            localNodeCoord2[2:3:end] .= localNodeCoord[2:2:end]

            # Bázisfüggvények deriváltjai a referenciadoménben
            _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, localNodeCoord2, "GradLagrange")
            ∇h = reshape(dfun, :, numNodes)  # (3*numEval × numNodes)

            nnet = zeros(Int, length(elemTags[i]), numNodes)

            # --- Elemeken végig ---
            for j in 1:length(elemTags[i])
                elem = elemTags[i][j]
                @inbounds for k in 1:numNodes
                    nnet[j, k] = elemNodeTags[i][(j-1)*numNodes + k]
                end

                xx = ncoord2[nnet[j, :]*3 .- 2]
                yy = ncoord2[nnet[j, :]*3 .- 1]
                zz = ncoord2[nnet[j, :]*3 .- 0]

                normvec_elem = zeros(3 * numNodes)

                # --- Csomópontonként ---
                for k in 1:numNodes
                    # Kiválasztjuk a lokális ∂N/∂ξ és ∂N/∂η értékeket a k-adik pontban
                    dN_dξ = ∇h[1:3:end, k]
                    dN_dη = ∇h[2:3:end, k]

                    # Jacobi-mátrix komponensek (a k-adik pontban)
                    dx_dξ = dot(dN_dξ, xx)
                    dy_dξ = dot(dN_dξ, yy)
                    dz_dξ = dot(dN_dξ, zz)
                    dx_dη = dot(dN_dη, xx)
                    dy_dη = dot(dN_dη, yy)
                    dz_dη = dot(dN_dη, zz)

                    # Kereszt-szorzat (∂x/∂ξ × ∂x/∂η)
                    nx = dy_dξ * dz_dη - dz_dξ * dy_dη
                    ny = dz_dξ * dx_dη - dx_dξ * dz_dη
                    nz = dx_dξ * dy_dη - dy_dξ * dx_dη

                    vnorm = √(nx^2 + ny^2 + nz^2)
                    normvec_elem[3*k-2] = nx / vnorm
                    normvec_elem[3*k-1] = ny / vnorm
                    normvec_elem[3*k-0] = nz / vnorm
                end

                push!(numElem, elem)
                push!(normalvectors, reshape(normvec_elem, :, 1))
            end
        end
    end

    return VectorField(normalvectors, [;;], [0.0], numElem, 1, :v3D, problem)
end

"""
    tangentVector(problem::Problem, phName::String) -> VectorField

Compute unit tangent vectors for all nodes of a curve-type physical group
in a 3D Gmsh model.

# Arguments
- `problem::Problem`: A `Problem` structure containing the current Gmsh model
  (name, dimension, number of nodes, etc.).
- `phName::String`: The name of a physical curve group in Gmsh for which the
  tangent vectors are computed.

# Description
The function sets the current model, queries the elements and nodes that belong
to the given 1D physical group, and evaluates the gradients of the Lagrange
basis functions to determine the local mapping derivative ∂x/∂ξ. Since this
derivative represents the geometric tangent direction of the curve at each
evaluation point, it is normalized to produce a unit tangent vector.

Each node belonging to the physical curve group is assigned the corresponding
3D unit tangent vector aligned with the parametric direction of the curve.

# Returns
- `VectorField`: A `VectorField` structure that stores the nodal tangent vectors
  on the given physical curve.

# Errors
- Throws an error if the provided physical group is not of curve type (`dim != 1`).

# Example
```julia
using LowLevelFEM

# Load a 3D geometry containing a curved edge
gmsh.initialize()
gmsh.open("curve_geometry.msh")

# Create a problem structure
mat = material("body")
prob = Problem([mat])

# Compute tangent vectors along the curve named "leadingEdge"
tv = tangentVector(prob, "leadingEdge")

# Visualize the tangent field
showDoFResults(tv)
openPostProcessor()
```
"""
function tangentVector(problem::Problem, phName::String)
    gmsh.model.setCurrent(problem.name)

    non = problem.non
    ncoord2 = zeros(3 * non)
    numElem = Int[]
    tangentvectors = Vector{Matrix{Float64}}()

    # Fizikai csoport (1D görbe)
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

    for (edim, etag) in dimTags
        if edim != 1
            error("tangentVector: physical group is not 1D (curve).")
        end

        elemTypes, elemTags, elemNodeTags =
            gmsh.model.mesh.getElements(edim, etag)

        # Csomópontok és 3D koordináták
        nodeTags, ncoord, _ =
            gmsh.model.mesh.getNodes(edim, etag, true, false)

        ncoord2[nodeTags*3 .- 2] = ncoord[1:3:end]
        ncoord2[nodeTags*3 .- 1] = ncoord[2:3:end]
        ncoord2[nodeTags*3 .- 0] = ncoord[3:3:end]

        # ----- Elemciklus -----
        for i in 1:length(elemTypes)
            et = elemTypes[i]

            elementName, eldim, order, numNodes, localNodeCoord, _ =
                gmsh.model.mesh.getElementProperties(et)

            # 1D → nincs 3D kiterjesztés, közvetlenül használható
            localNodeCoord2 = zeros(length(localNodeCoord) * 3)
            localNodeCoord2[1:3:end] .= localNodeCoord[1:end]
            #localNodeCoord2[2:3:end] .= localNodeCoord[2:end]

            # ∂N/∂ξ a lokális pontokban
            _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, localNodeCoord2, "GradLagrange")
            dN_dξ = reshape(dfun, :, numNodes)  # (1*numEval × numNodes)

            nnet = zeros(Int, length(elemTags[i]), numNodes)

            for j in 1:length(elemTags[i])
                elem = elemTags[i][j]

                @inbounds for k in 1:numNodes
                    nnet[j,k] = elemNodeTags[i][(j-1)*numNodes + k]
                end

                xx = ncoord2[nnet[j,:]*3 .- 2]
                yy = ncoord2[nnet[j,:]*3 .- 1]
                zz = ncoord2[nnet[j,:]*3 .- 0]

                tang_elem = zeros(3 * numNodes)

                # ----- Csomópontonként -----
                for k in 1:numNodes
                    dNk = dN_dξ[1:3:end,k]    # ∂Nk/∂ξ (numNodes hosszú sor)

                    # ∂x/∂ξ, ∂y/∂ξ, ∂z/∂ξ
                    dx_dξ = dot(dNk, xx)
                    dy_dξ = dot(dNk, yy)
                    dz_dξ = dot(dNk, zz)

                    # normalizált érintő
                    vnorm = √(dx_dξ^2 + dy_dξ^2 + dz_dξ^2)
                    if vnorm > 0
                        dx = dx_dξ / vnorm
                        dy = dy_dξ / vnorm
                        dz = dz_dξ / vnorm
                    else
                        dx = dy = dz = 0.0
                    end

                    tang_elem[3*k-2] = dx
                    tang_elem[3*k-1] = dy
                    tang_elem[3*k]   = dz
                end

                push!(numElem, elem)
                push!(tangentvectors, reshape(tang_elem, :, 1))
            end
        end
    end

    return VectorField(tangentvectors, [;;], [0.0], numElem, 1, :v3D, problem)
end

function normalVector1D(problem::Problem, phName::String)
    gmsh.model.setCurrent(problem.name)

    non = problem.non
    ncoord2 = zeros(3 * non)
    numElem = Int[]
    normalvectors = Vector{Matrix{Float64}}()

    # Fizikai csoport (1D görbe)
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

    for (edim, etag) in dimTags
        if edim != 1
            error("tangentVector: physical group is not 1D (curve).")
        end

        elemTypes, elemTags, elemNodeTags =
            gmsh.model.mesh.getElements(edim, etag)

        # Csomópontok és 3D koordináták
        nodeTags, ncoord, _ =
            gmsh.model.mesh.getNodes(edim, etag, true, false)

        ncoord2[nodeTags*3 .- 2] = ncoord[1:3:end]
        ncoord2[nodeTags*3 .- 1] = ncoord[2:3:end]
        ncoord2[nodeTags*3 .- 0] = ncoord[3:3:end]

        # ----- Elemciklus -----
        for i in 1:length(elemTypes)
            et = elemTypes[i]

            elementName, eldim, order, numNodes, localNodeCoord, _ =
                gmsh.model.mesh.getElementProperties(et)

            # 1D → nincs 3D kiterjesztés, közvetlenül használható
            localNodeCoord2 = zeros(length(localNodeCoord) * 3)
            localNodeCoord2[1:3:end] .= localNodeCoord[1:end]
            #localNodeCoord2[2:3:end] .= localNodeCoord[2:end]

            # ∂N/∂ξ a lokális pontokban
            _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, localNodeCoord2, "GradLagrange")
            dN_dξ = reshape(dfun, :, numNodes)  # (1*numEval × numNodes)

            nnet = zeros(Int, length(elemTags[i]), numNodes)

            for j in 1:length(elemTags[i])
                elem = elemTags[i][j]

                @inbounds for k in 1:numNodes
                    nnet[j,k] = elemNodeTags[i][(j-1)*numNodes + k]
                end

                xx = ncoord2[nnet[j,:]*3 .- 2]
                yy = ncoord2[nnet[j,:]*3 .- 1]
                zz = ncoord2[nnet[j,:]*3 .- 0]

                tang_elem = zeros(3 * numNodes)
                norm_elem = zeros(3 * numNodes)

                # ----- Csomópontonként -----
                for k in 1:numNodes
                    dNk = dN_dξ[1:3:end,k]    # ∂Nk/∂ξ (numNodes hosszú sor)

                    # ∂x/∂ξ, ∂y/∂ξ, ∂z/∂ξ
                    dx_dξ = dot(dNk, xx)
                    dy_dξ = dot(dNk, yy)
                    dz_dξ = dot(dNk, zz)

                    # normalizált érintő
                    vnorm = √(dx_dξ^2 + dy_dξ^2 + dz_dξ^2)
                    if vnorm > 0
                        dx = dx_dξ / vnorm
                        dy = dy_dξ / vnorm
                        dz = dz_dξ / vnorm
                    else
                        dx = dy = dz = 0.0
                    end

                    tang_elem[3*k-2] = dx
                    tang_elem[3*k-1] = dy
                    tang_elem[3*k]   = dz
                end

                for k in 1:numNodes
                    dNk = dN_dξ[1:3:end,k]    # ∂Nk/∂ξ (numNodes hosszú sor)

                    # ∂x/∂ξ, ∂y/∂ξ, ∂z/∂ξ
                    dtx_dξ = dot(dNk, tang_elem[1:3:end])
                    dty_dξ = dot(dNk, tang_elem[2:3:end])
                    dtz_dξ = dot(dNk, tang_elem[3:3:end])
                    dx_dξ = dot(dNk, xx)
                    dy_dξ = dot(dNk, yy)
                    dz_dξ = dot(dNk, zz)

                    # normalizált érintő
                    vnorm = √(dx_dξ^2 + dy_dξ^2 + dz_dξ^2)
                    if vnorm > 0
                        dx = dtx_dξ / vnorm
                        dy = dty_dξ / vnorm
                        dz = dtz_dξ / vnorm
                    else
                        dx = dy = dz = 0.0
                    end

                    norm_elem[3*k-2] = dx
                    norm_elem[3*k-1] = dy
                    norm_elem[3*k]   = dz
                end

                push!(numElem, elem)
                push!(normalvectors, reshape(norm_elem, :, 1))
            end
        end
    end

    return VectorField(normalvectors, [;;], [0.0], numElem, 1, :v3D, problem)
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
        #name, ux, uy, uz = supports[i]
        name = supports[i].phName
        ux = supports[i].ux
        uy = supports[i].uy
        uz = supports[i].uz
        T = supports[i].T
        sx = supports[i].sx
        sy = supports[i].sy
        sz = supports[i].sz
        sxy = supports[i].sxy
        syz = supports[i].syz
        szx = supports[i].szx
        
        T !== nothing &&
            (ux !== nothing || uy !== nothing || uz !== nothing) &&
            (sx !== nothing || sy !== nothing || sz !== nothing || sxy !== nothing || syz !== nothing || szx !== nothing) &&
            error("applyBoundaryConditions: only T or ux/uy/uz or sx/sy/sz/sxy/syz/szx can be given as BoundaryCondition.")
        phg = getTagForPhysicalName(name)
        nodeTags::Vector{Int64}, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        nodeTagsX = []
        nodeTagsY = []
        nodeTagsZ = []
        nodeTagsXY = []
        nodeTagsYZ = []
        nodeTagsZX = []
        nodeTagsYX = []
        nodeTagsZY = []
        nodeTagsXZ = []
        if T !== nothing && pdim == 1
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
        end
        if ux !== nothing
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
        end
        if uy !== nothing
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
        end
        if pdim == 3 && uz !== nothing
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= pdim
            nodeTagsZ .-= (pdim-0)
        end
        if pdim == 9 && sx !== nothing
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim - 1)
        end
        if pdim == 9 && sy !== nothing
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim - 5)
        end
        if pdim == 9 && sz !== nothing
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= pdim
            nodeTagsZ .-= (pdim - 9)
        end
        if pdim == 9 && sxy !== nothing
            nodeTagsXY = copy(nodeTags)
            nodeTagsXY *= pdim
            nodeTagsXY .-= (pdim - 4)
        end
        if pdim == 9 && sxy !== nothing
            nodeTagsYX = copy(nodeTags)
            nodeTagsYX *= pdim
            nodeTagsYX .-= (pdim - 2)
        end
        if pdim == 9 && syz !== nothing
            nodeTagsYZ = copy(nodeTags)
            nodeTagsYZ *= pdim
            nodeTagsYZ .-= (pdim - 8)
        end
        if pdim == 9 && syz !== nothing
            nodeTagsZY = copy(nodeTags)
            nodeTagsZY *= pdim
            nodeTagsZY .-= (pdim - 6)
        end
        if pdim == 9 && szx !== nothing
            nodeTagsZX = copy(nodeTags)
            nodeTagsZX *= pdim
            nodeTagsZX .-= (pdim - 3)
        end
        if pdim == 9 && szx !== nothing
            nodeTagsXZ = copy(nodeTags)
            nodeTagsXZ *= pdim
            nodeTagsXZ .-= (pdim - 7)
        end
        cdofs = unique(cdofs ∪ nodeTagsX ∪ nodeTagsY ∪ nodeTagsZ ∪ nodeTagsXY ∪ nodeTagsYZ ∪ nodeTagsZX ∪ nodeTagsYX ∪ nodeTagsZY ∪ nodeTagsXZ)
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
    
    for e in 1:length(numElem)
        elementType, nodeTags, dim, tag = gmsh.model.mesh.getElement(numElem[e])
        for i in 1:length(nodeTags)
            s[(nodeTags[i]-1) * epn + 1: nodeTags[i] * epn, :] .+= σ[e][(i-1)*epn+1:i*epn, :]
            pcs[nodeTags[i]] += 1
        end
    end
    for l in 1:non
        s[epn * (l - 1) + 1: epn * l, :] ./= pcs[l] == 0 ? 1 : pcs[l]
    end
    return T([], s, S.t, [], S.nsteps, type, problem)
end

"""
    nodesToElements(T, onPhysicalGroup="")

Solves the element results `F` from the nodal results `T`.
`T` can be ScalarField, VectorField or TensorField.
If `onPhysicalGroup` is an existing physical group in gmsh, field `T`
will be solve only on elements belonging to that physical group. Dimension
of physical group can be different than the dimension of the problem.
(eg. mapping from 3D volume to a 2D surface)

Return: `F`

Types:
- `T`: ScalarField, VectorField or TensorField
- `F`: ScalarField, VectorField or TensorField
- `onPhysicalGroup`: String
"""
function nodesToElements(r::Union{ScalarField,VectorField,TensorField}; onPhysicalGroup="")
    problem = r.model
    gmsh.model.setCurrent(problem.name)
    
    if isElementwise(r)
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
    nom = length(problem.material)
    phName = ""
    if onPhysicalGroup ≠ ""
        nom = 1
        phName = onPhysicalGroup
    end
    
    for ipg in 1:nom
        if onPhysicalGroup == ""
            phName = problem.material[ipg].phName
        end
        #ν = problem.material[ipg].ν
        dim = problem.dim
        #=
        if problem.dim == 3 && problem.type == :Solid
            dim = 3
            rowsOfH = 3
        elseif problem.dim == 3 && problem.type == :Truss
            dim = 3
            #rowsOfH = 3
        elseif problem.dim == 2 # && problem.type == :PlaneStress
            dim = 2
            rowsOfH = 2
        #elseif problem.dim == 2 && problem.type == :PlaneStrain
        #    dim = 2
        #    rowsOfH = 2
        #elseif problem.dim == 2 && problem.type == :AxiSymmetric
        #    dim = 2
        #    rowsOfH = 2
        else
            error("nodesToElements: dimension is $(problem.dim), problem type is $(problem.type).")
        end
        =#
        
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
    return !(a.a == [;;])
    if a.a != [;;] && a.A == []
        return true
    elseif a.a == [;;] && a.A != []
        return false
    else
        error("isNodal: Scalar field is empty.")
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
    return !(a.A == [])
    if a.a == [;;] && a.A != []
        return true
    elseif a.a != [;;] && a.A == []
        return false
    else
        error("isElementwise: Scalar field is empty.")
    end
end

"""
    expandTo3D(v2D::VectorField)

Expand a 2D vector field into 3D by adding a zero z-component.

return: VectorField

# Examples
```julia
V3D = expandTo3D(V2D)
```
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

Project a 3D vector field onto the xy-plane by discarding the z-component.

return: VectorField

# Examples
```julia
V2D = expandTo3D(V3D)
```
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
    resultant(field, phName)

Computes the resultant of vector field `field` on the physical group `phName`.
Returns the resultant(s) in a `tuple`. The number of elements in the tuple depends on the
dimension of problem (dimension of `field`).
It can solve for example the resultant of a load vector (sum of the elements of the vector).

Return: `resx`

or

Return: `resx`, `resy`

or

Return: `resx`, `resy`, `resz`

Types:
- `field`: VectorField
- `phName`: String 
- `resx`: Float64 
- `resy`: Float64 
- `resz`: Float64 
"""
function resultant(field::VectorField, phName::String)
    axiSymmetric = false
    problem = field.model
    if isElementwise(field)
        error("resultant: only nodal results are permitted.")
    end
    if problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction
        axiSymmetric = true
    end
    ph1 = getTagForPhysicalName(phName)
    nodes0, coords = gmsh.model.mesh.getNodesForPhysicalGroup(-1,ph1)
    nodes = Vector{Int64}(nodes0)
    dim = problem.dim
    s = [0.0, 0.0, 0.0]
    for i in 1:dim
        for j in 1:length(nodes)
            b = axiSymmetric == true ? 2π * coords[3j-2] : 1
            s[i] += field.a[dim * nodes[j] - (dim - i)] * b
        end
    end
    if dim == 1
        return [s[1]]
    elseif dim == 2
        return [s[1], s[2]]
    elseif dim == 3
        return [s[1], s[2], s[3]]
    end
end

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
    integrate(problem::Problem, phName::String, f::Union{Function,ScalarField}; step::Int64=1, order::Int64=...)
    ∫(problem::Problem, phName::String, f::Union{Function,ScalarField}; step::Int64=1, order::Int64=...)

Integrates the function or scalar field `f` over the physical group `phName` defined in the geometry of `problem`.
If `f` is a `ScalarField`, the time step `step` will be integrated. The optional parameter `order` controls the 
numerical integration rule. If `order > 0`, it is used as a hint for selecting the Gauss quadrature order
(i.e. the number of integration points) employed during integration. The exactness of the integration depends 
on the element geometry and the regularity of the integrand.

Returns: integral

Types:
- `problem`: Problem
- `phName`: String
- `f`: Function (of x, y and z) or ScalarField
- `integral`: Number

Examples:

```julia
f(x, y, z) = x^2 + y^2
Iz = integrate(prob, "body", f)
```
"""
function integrate(problem::Problem, phName::String, f::Union{Function,ScalarField}; step::Int64=1, order::Int64=0)
    gmsh.model.setCurrent(problem.name)
    f2 = 0
    if f isa ScalarField
        f2 = elementsToNodes(f)
    end
    DIM = problem.dim
    b = problem.thickness
    ncoord2 = zeros(3 * problem.non)
    sum0 = 0
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
    for i ∈ 1:length(dimTags)
        dimTag = dimTags[i]
        dim = dimTag[1]
        tag = dimTag[2]
        elementTypes, elementTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
        nodeTags::Vector{Int64}, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, tag, true, false)
        #ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
        #ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
        #ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
        @inbounds for (local_idx, nd) in enumerate(nodeTags)
            ci = 3*(nd-1)
            cn = 3*(local_idx-1)
            ncoord2[ci+1] = ncoord[cn+1]
            ncoord2[ci+2] = ncoord[cn+2]
            ncoord2[ci+3] = ncoord[cn+3]
        end

        for ii in 1:length(elementTypes)
            elementName, dim, order1, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementTypes[ii])
            qorder = order <= 0 ? order1 : order
            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(qorder))
            numIntPoints = length(intWeights)
            comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
            h = reshape(fun, :, numIntPoints)
            nnet = zeros(Int, length(elementTags[ii]), numNodes)
            @inbounds for l in 1:length(elementTags[ii])
                elem = elementTags[ii][l]
                #for k in 1:numNodes
                #    nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                #end
                nodes = elemNodeTags[ii]
                offset = (l - 1) * numNodes
                nodeids = nodes[offset + 1 : offset + numNodes]
                jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                Jac = reshape(jac, 3, :)
                s1 = 0
                first = (l - 1) * numNodes + 1
                @inbounds for j in 1:numIntPoints
                    x = 0.0
                    y = 0.0
                    z = 0.0
                    @inbounds for nk in 1:numNodes
                        nd    = nodes[first + nk - 1]      # globális csomópont index
                        baseg = 3 * (nd - 1)
                        hk    = h[nk, j]
                        x += hk * ncoord2[baseg + 1]
                        y += hk * ncoord2[baseg + 2]
                        z += hk * ncoord2[baseg + 3]
                    end
                    ff = 0
                    if f isa Function
                        ff = f(x, y, z)
                    elseif f isa ScalarField
                        #ff = h[:, j]' * f2.a[nnet[l, :]]
                        ff = h[:, j]' * f2.a[nodeids, step]
                    else
                        error("integrate: 3rd argument must be a Function or a ScalarField")
                    end
                    r = x
                    ############### NANSON ######## 3D ###################################
                    if (DIM == 3 || DIM == 2) && dim == 3
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
                        error("integrate: dimension of the problem is $(problem.dim), dimension of ? is $dim.")
                    end
                    #@disp ff
                    s1 += ff * Ja * intWeights[j]
                end
                sum0 += s1
            end
        end
    end
    return sum0
end

∫(problem::Problem, phName::String, f::Union{Function,ScalarField}; step::Int64=1) = integrate(problem, phName, f, step=step)

"""
    time_integral!(s, t, d; s0 = 0.0)

Reconstruct a time-dependent signal `s` from its time derivative `d`
using the inverse of the second-order central difference scheme.

The reconstructed signal is unique up to an additive time-independent
constant, specified by `s0`.
"""

function time_integral!(s, t, d; s0=0.0)
    nsteps = length(t)

    # reference value
    s[:, 1] .= s0

    if nsteps == 1
        return s
    end

    if nsteps == 2
        # linear reconstruction
        dt = t[2] - t[1]
        s[:, 2] .= s[:, 1] .+ dt .* d[:, 1]
        return s
    end

    @inbounds begin
        # bootstrap second step using forward Euler
        dt = t[2] - t[1]
        s[:, 2] .= s[:, 1] .+ dt .* d[:, 1]

        # central inverse recursion
        for i in 2:nsteps-1
            s[:, i+1] .= s[:, i-1] .+ (t[i+1] - t[i-1]) .* d[:, i]
        end
    end

    return s
end

"""
    integrate(s::Union{ScalarField,VectorField,TensorField})
    ∫(s::Union{ScalarField,VectorField,TensorField})

Compute the time integral of a time-dependent field using a discrete
inverse of the central finite-difference time derivative.

The result reproduces the original field (up to a time-independent
constant) when applied to a field obtained by `∂t`.
"""
function integrate(s::Union{ScalarField,VectorField,TensorField})
    s0 = 0.0
    T = typeof(s)

    if isNodal(s)
        a = similar(s.a)
        time_integral!(a, s.t, s.a; s0=s0)
        return T([], a, s.t, [], s.nsteps, s.type, s.model)

    elseif isElementwise(s)
        A = [similar(s.A[e]) for e in eachindex(s.A)]
        for e in eachindex(s.A)
            time_integral!(A[e], s.t, s.A[e]; s0=s0)
        end
        return T(A, [;;], s.t, s.numElem, s.nsteps, s.type, s.model)

    else
        error("∫t: internal error")
    end
end

∫(s::Union{ScalarField,VectorField,TensorField}) = integrate(s)

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
    showDoFResults(q, comp; name=..., visible=..., factor=0)

Loads nodal results into a View in Gmsh. `q` is the field to show, `comp` is
the component of the field (:vector, :uvec, :ux, :uy, :uz, :vvec, :vx, :vy, :vz,
:qvec, :qx, :qy, :qz, :T, :p, :qn, :s, :sx, :sy, :sz, :sxy, :syx, :syz,
:szy, :szx, :sxz, :e, :ex, :ey, :ez, :exy, :eyx, :eyz, :ezy, :ezx, :exz, :seqv, :scalar, :tensor),
`name` is a title to display and `visible` is a Boolean to toggle the initial visibility in Gmsh on or off.
If `q` has more columns, then a sequence of results will be shown (e.g., as an animation).
`factor` multiplies the DoF result to increase for better visibility.
This function returns the tag of the View.

Returns: `tag`

Types:
- `q`: ScalarField, VectorField or TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showDoFResults(q::Union{ScalarField,VectorField,TensorField}, comp::Symbol; name=comp, visible=false, ff = 0, factor=0)
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
    #display("showDoFResults: dim=$pdim, type=$(q.type).")
    nodeTags = []
    u = []
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
    if problem.type == :Reynolds
        if problem.geometry.nameVolume ≠ ""
            phName = problem.geometry.nameVolume
            tag = getTagForPhysicalName(phName)
            nT, coords = gmsh.model.mesh.getNodesForPhysicalGroup(edim + 1, tag)
            append!(nodeTags, nT)
        end
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
    
    gmsh.view.option.setNumber(uvec, "GlyphLocation", 2)
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

function showDoFResults(q::Union{ScalarField,VectorField,TensorField}; name=q.type, visible=false, ff = 0, factor=0)
    if q isa ScalarField
        showDoFResults(q, :scalar, name=name, visible=visible, factor=factor)
    elseif q isa VectorField
        showDoFResults(q, :vector, name=name, visible=visible, factor=factor)
    elseif q isa TensorField
        showDoFResults(q, :tensor, name=name, visible=visible, factor=factor)
    else
        error("showDoFResults: argument must be a ScalarField, VectorField or TensorField.")
    end
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
        if comp == :e || comp == :tensor
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
            elseif comp == :exy
                k = 5
            elseif comp == :eyz
                k = 1
            elseif comp == :ezx
                k = 6
            elseif comp == :eyx
                k = 7
            elseif comp == :ezy
                k = 3
            elseif comp == :exz
                k = 2
            else
                error("ShowStrainResults: component is $comp ????")
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
function showElementResults(F::Union{ScalarField,VectorField,TensorField}, comp; name=comp, visible=false, smooth=false, factor=0)
    if F.type == :e || F.type == :tensor && isElementwise(F)
        return showStrainResults(F, comp, name=name, visible=visible, smooth=smooth)
    elseif F.type == :s && isElementwise(F)
        return showStressResults(F, comp, name=name, visible=visible, smooth=smooth)
    elseif F isa VectorField && isElementwise(F)
        return showHeatFluxResults(F, comp, name=name, visible=visible, smooth=smooth, factor=factor)
    elseif F isa ScalarField && isElementwise(F)
        return showScalarResults(F, name=name, visible=visible, smooth=smooth, factor=factor)
    elseif isNodal(F)
        return showDoFResults(F, comp, name=name, visible=visible, smooth=smooth, factor=factor)
    else
        error("showElementResults: type is '$(F.type)'")
    end
end

function showElementResults(q::Union{ScalarField,VectorField,TensorField}; name="__field__", visible=false, smooth=false, ff = 0, factor=0)
    if q isa ScalarField
        showElementResults(q, :scalar, name=name == "__field__" ? :scalar : name, visible=visible, smooth=smooth)
    elseif q isa VectorField
        showElementResults(q, :vector, name=name == "__field__" ? :vector : name, visible=visible, smooth=smooth, factor=factor)
    elseif q isa TensorField
        showElementResults(q, :tensor, name=name == "__field__" ? :tensor : name, visible=visible, smooth=smooth)
    else
        error("showElementResults: argument must be a ScalarField, VectorField or TensorField.")
    end
end

function showStressResults(q::TensorField; name="StressField", visible=false, smooth=false, ff = 0, factor=0)
    if q isa TensorField
        showStressResults(q, :s, name=name, visible=visible, smooth=smooth)
    else
        error("showStressResults: argument must be a TensorField.")
    end
end

function showStrainResults(q::TensorField; name="StrainField", visible=false, ff = 0, factor=0)
    if q isa TensorField
        showStrainResults(q, :e, name=name, visible=visible, smooth=smooth)
    else
        error("showStrainResults: argument must be a TensorField.")
    end
end

function showHeatFluxResults(q::VectorField; name="VectorField", visible=false, ff = 0, factor=0)
    if q isa VectorField
        showHeatFluxResults(q, :vector, name=name, visible=visible, smooth=smooth)
    else
        error("showHeatFluxResults: argument must be a VectorField.")
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
function showStressResults(S::TensorField, comp; name=comp, visible=false, smooth=false)
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
        if comp == :s || comp == :tensor
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
function showHeatFluxResults(S::VectorField, comp; name=comp, visible=false, smooth=true, factor=0)
    problem = S.model
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    if S.type == :v3D
        dim = 3
    elseif S.type == :v2D
        dim = 2
    else
        error("showHeatFluxResults: type is $(S.type)")
    end
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
    
    gmsh.view.option.setNumber(SS, "GlyphLocation", 2)
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

function showScalarResults(S; name="ScalarField", visible=false, smooth=false, factor=0)
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
    plotOnPath(pathName, field; points=100, step=..., plot=..., name=..., visible=..., offsetX=..., offsetY=..., offsetZ=...)

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
function showOnSurface(field::Number, phName::String; grad=false, component=:x, offsetX=0, offsetY=0, offsetZ=0, name=phName, visible=false)
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

function showOnSurface(field::ScalarField, phName::String; name=phName)
    s = elementsToNodes(field)
    s = nodesToElements(s, onPhysicalGroup=phName)
    view = showElementResults(s, name=name)
    return view
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
