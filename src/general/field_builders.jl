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
struct Eigen{T}
    f::Vector{T}
    ϕ::Matrix{T}
    model::Union{Problem,Nothing}
    problems::Union{Vector{Problem},Nothing}
    offsets::Union{Vector{Int},Nothing}
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
    material(name; type=:Hooke, E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5, μ=E/(1+ν)/2, λ=2μ*ν/(1-2ν), κ=2μ*(1+ν)/(1-2ν)/3, η=1e-7, p₀=0.1, A=1)

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

#=
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
=#

"""
    displacementConstraint(name; ux=nothing, uy=nothing, uz=nothing)

Specifies displacement constraints on the physical group `name`. At least one of `ux`, `uy`,
or `uz` must be provided (depending on the problem dimension). `ux`, `uy`, or `uz` can be a
constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); displacementConstraint("support1", ux=fn)`)

Returns: BoundaryCondition

# Examples

```julia
hc = heatConvection("outer", h=10.0, T∞=20.0)
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
    load(name; fx=nothing, fy=nothing, fz=nothing, p=nothing) -> BoundaryCondition

Create engineering distributed-load data on physical group `name`. Supply at least
one Cartesian component or pressure `p`. Values may be numbers, spatial
functions `(x,y,z)->value`, or scalar fields (converted to nodal storage).
If `p` is supplied it takes precedence over Cartesian component keywords.

This declares load data consumed by the engineering assembly interface; it does
not itself assemble a DSL RHS. The engineering API name is intended to remain
when its backend migrates to DSL-backed wrappers. For direct weak forms use
`f = ∫(U ⋅ traction; Γ="right")`.
"""
function load(name; fx=nothing, fy=nothing, fz=nothing, p=nothing)
    fx === nothing && fy === nothing && fz === nothing && p === nothing &&
        error("load: at least one of fx, fy, fz, p must be specified.")
    if p !== nothing
        if p isa ScalarField
            pp = elementsToNodes(p)
            return LoadCondition(name, p=pp)
        end
        return LoadCondition(name, p=p)
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
    return LoadCondition(name, fx=ffx, fy=ffy, fz=ffz)
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
    elasticSupport(name; kx=nothing, ky=nothing, kz=nothing)

Specifies distributed stiffness for an elastic support on the physical group `name`.
`kx`, `ky`, or `kz` can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); elasticSupport("supp1", kx=fn)`)
Default values are 1.

Returns: BoundaryCondition

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
    temperatureConstraint(name; T=nothing)

Specifies temperature constraints on the physical group `name`.
`T` can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); temperatureConstraint("surf1", T=fn)`)

Returns: BoundaryCondition

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
    heatFlux(name; qn=nothing)

Specifies the heat flux normal to the surface of the physical group `name`.
`qn` can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); load("flux1", qn=fn)`)

Returns: BoundaryCondition

Types:
- `name`: String
- `qn`: Float64 or Function
"""
function heatFlux(name; qn=nothing)
    qn === nothing &&
        error("heatFlux: qn must be specified.")
    return LoadCondition(name, qn=qn)
    #=
    p1 =0
    p2 =0
    qn0 = -qn
    fl0 = name, qn0, p1, p2
    return fl0
    =#
end

"""
    heatSource(name; h=nothing)

Specifies the volumetric heat source in the physical group `name`.
`h` can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); load("source1", h=fn)`)

Returns: BoundaryCondition

Types:
- `name`: String
- `h`: Float64 or Function
"""
function heatSource(name; h=nothing)
    h === nothing &&
        error("heatSource: h must be specified.")
    return LoadCondition(name, h=h)
    #=
    p1 =0
    p2 =0
    h0 = -h
    sr0 = name, h0, p1, p2
    return sr0
    =#
end

"""
    heatConvection(name; h=nothing, T∞=nothing)

Specifies convective boundary conditions on the surface in the physical group `name`.
`h` is the heat transfer coefficient of the surrounding medium; `Tₐ` is the ambient temperature.
`Tₐ` can be either a constant or a function of `x`, `y`, and `z`.

Returns: BoundaryCondition

Types:
- `name`: String
- `h`: Float64
- `Tₐ`: Float64 or Function
"""
function heatConvection(name; h=nothing, T∞=nothing)
    h === nothing || T∞ === nothing &&
        error("heatSource: h and T∞ must be specified.")
    return BoundaryCondition(name, h=h, T∞=T∞)
end

"""
    generateMesh(surf, elemSize; approxOrder=1, algorithm=6, quadrangle=0, internalNodes=0)

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
    field(name; f=nothing, fx=nothing, fy=nothing, fz=nothing, 
                fxy=nothing, fyx=nothing, fyz=nothing, 
                fzy=nothing, fzx=nothing, fxz=nothing)

Specifies the value of a scalar, vector, or tensor field on the physical group `name`.
At least one of `fx`, `fy`, or `fz` etc. must be provided (depending on the problem dimension).
Each component can be a constant or a function of `x`, `y`, and `z`.
(e.g., `fn(x,y,z) = 5*(5-x); field("surf1", fx=fn)`)

Returns: BoundaryCondition

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
- `dataField`: Vector{BoundaryCondition}

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
        #f = dataField[i].f
        vals = dataField[i].values
        f = get(vals, :f, nothing)
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
    vectorField(problem::Problem, dataField::Vector{BoundaryCondition})

Defines a vector field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`. `type` can be an arbitrary `Symbol`,
e.g., `:u` or `:f`.

Return: VectorField

Types:
- `problem`: Problem
- `dataField`: Vector{BoundaryCondition}

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
        error("vectorField: dataField are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    pdim = problem.dim
    non = problem.non
    field = zeros(non * pdim)
    
    for i in 1:length(dataField)
        #name, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = dataField[i]
        name = dataField[i].phName
        vals = dataField[i].values
        fx = get(vals, :fx, nothing)
        fy = get(vals, :fy, nothing)
        fz = get(vals, :fz, nothing)
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
    tensorField(problem::Problem, dataField::Vector{BoundaryCondition}; type=:e)

Defines a vector field from `dataField`, which is a tuple of `name` of physical group and
prescribed values or functions. Mesh details are in `problem`. `type` can be an arbitrary `Symbol`,
e.g., `:u` or `:f`.

Return: TensorField

Types:
- `problem`: Problem
- `dataField`: Vector{BoundaryCondition}

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
        error("tensorField!: dataField are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    pdim = 9
    non = problem.non
    field = zeros(non * pdim)
    
    for i in 1:length(dataField)
        #name, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = dataField[i]
        name = dataField[i].phName
        vals = dataField[i].values
        fx = get(vals, :fx, nothing)
        fy = get(vals, :fy, nothing)
        fz = get(vals, :fz, nothing)
        fxy = get(vals, :fxy, nothing)
        fyz = get(vals, :fyz, nothing)
        fzx = get(vals, :fzx, nothing)
        fyx = get(vals, :fyx, nothing)
        fzy = get(vals, :fzy, nothing)
        fxz = get(vals, :fxz, nothing)
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

