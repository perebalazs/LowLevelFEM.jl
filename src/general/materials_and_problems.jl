"""
    Material(phName; type=:Hooke, E=nothing, ν=nothing, λ=nothing,
             μ=nothing, κ=nothing, ρ=7.85e-9, k=45.0, c=4.2e8,
             α=1.2e-5, η=1e-7, p₀=0.1, A=1.0)

Material parameters assigned to a physical group.

This structure stores mechanical, thermal and fluid-related material
properties used during assembly. The material is associated with a
Gmsh physical group via `phName`.

# Arguments

- `phName::String`  
  Name of the physical group where the material is applied.

# Keyword Arguments

## Elastic properties
With no elastic constants, use `E=2e5` and `ν=0.3`. Supplying only `E`
keeps `ν=0.3`; supplying only `κ` keeps `E=2e5`. Supported pairs are
`(E,ν)`, `(μ,ν)`, `(λ,μ)`, and `(κ,μ)`. Other combinations are rejected;
missing elastic constants are computed from the selected combination.

`type` selects `:Hooke`, `:StVenantKirchhoff`, `:NeoHookeCompressible`, or
`:JFO`. It is material metadata; the selected assembly or constitutive routine
must support that material law. Parameters use the user's consistent unit
system; the defaults are not converted automatically.

- `E::Float64`  
  Young's modulus
- `ν::Float64`  
  Poisson's ratio
- `λ::Float64`  
  First Lamé parameter
- `μ::Float64`  
  Second Lamé parameter (shear modulus)
- `κ::Float64`  
  Bulk modulus

## Density
- `ρ::Float64`  
  Mass density

## Thermal properties
- `k::Float64`  
  Thermal conductivity
- `c::Float64`  
  Specific heat capacity
- `α::Float64`  
  Thermal expansion coefficient

## Fluid / tribology-related properties
- `η::Float64`  
  Dynamic viscosity
- `p₀::Float64`  
  Cavitation pressure

## Structural property
- `A::Float64`  
  Cross-sectional area (used e.g. in truss elements)

# Notes

- Elastic constants are internally completed to a consistent set
  using `elastic_constants`.
- Default values are provided for all optional parameters.
- The material definition is independent of the field type
  (scalar, vector, tensor) and can be used in multifield problems.

# Examples

```julia
mat = Material("solid", E=210e3, ν=0.3)

mat = Material("fluid",
    η=1e-6,
    p₀=0.1
)

mat = Material("truss",
    E=210e3,
    A=1.0
)
```
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

    elseif n == 1 && haskey(specified, :κ)
        κ = specified[:κ]
        E = DEFAULT_YOUNG
        ν = (1 - E/3κ) / 2
        μ = E / (2*(1+ν))
        λ = 2μ*ν/(1-2ν)
        #κ = E / (3*(1-2ν))

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

"""
    material(phName::String; kwargs...) -> Material

Compatibility spelling of [`Material`](@ref), forwarding all material keywords
and defaults unchanged. Prefer `Material` in new models.
"""
function material(  phName::String;
        type = :Hooke,
        E = nothing, ν = nothing, λ = nothing, μ = nothing, κ = nothing,
        ρ = 7.85e-9,
        k = 45.0,
        c = 4.2e8,
        α = 1.2e-5,
        η = 1e-7,
        p₀ = 0.1,
        A = 1.0)
    return Material(phName, type=type, E=E, ν=ν, λ=λ, μ=μ, κ=κ, ρ=ρ, k=k, c=c, α=α, η=η, p₀=p₀, A=A)
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

"""
    renumberNodes!(material::Material; method::Symbol=:RCMK)
    renumberNodes!(materials::Vector{Material}; method::Symbol=:RCMK)

Renumber the Gmsh mesh nodes associated with the given material domain(s).

The renumbering is computed from the finite elements belonging to the
physical groups referenced by the supplied `Material` objects. The supported
methods are `:RCMK`, `:Hilbert`, and `:Metis`.

If several material domains belong to the same finite-element problem, it is
recommended to pass all of them together in a vector so that the renumbering
is computed from the complete element connectivity of the problem.

# Warning

Node renumbering modifies the node tags of the Gmsh mesh and must therefore
be performed **before creating any `Problem` or `Field` objects** that use
the mesh.

Do not call `renumberNodes!` after a `Problem` or `Field` has already been
created, because existing node/DoF mappings may then become inconsistent
with the modified Gmsh mesh.

# Examples

```julia
mat = Material("body")

renumberNodes!(mat; method=:RCMK)

u = Field([mat], type=:VectorField, dim=3, field=:u)
```

For a problem containing several material domains:

```julia
mat1 = Material("body1")
mat2 = Material("body2")

renumberNodes!([mat1, mat2]; method=:RCMK)

u = Field([mat1, mat2], type=:VectorField, dim=3, field=:u)
```

"""
function renumberNodes!(
    materials::Vector{Material};
    method::Symbol=:RCMK
    )
    method in (:RCMK, :Hilbert, :Metis) ||
    error(
        "renumberNodes!: method can be `:RCMK`, `:Hilbert` or `:Metis`. " *
        "Now method = `$method`."
    )

    isempty(materials) &&
        error("renumberNodes!: material vector must not be empty.")

    elemTags = Int[]

    for material in materials
        phName = material.phName
   
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
   
        isempty(dimTags) &&
            error(
                "renumberNodes!: physical group '$phName' was not found."
            )
   
        for (edim, etag) in dimTags
            _, elementTags, _ =
                gmsh.model.mesh.getElements(edim, etag)
   
            for tags in elementTags
                append!(elemTags, tags)
            end
        end
    end

    isempty(elemTags) &&
        error(
            "renumberNodes!: no finite elements were found in the supplied material domains."
        )

    oldTags, newTags =
        gmsh.model.mesh.computeRenumbering(method, elemTags)

    gmsh.model.mesh.renumberNodes(oldTags, newTags)

    return nothing
end

function renumberNodes!(
    material::Material;
    method::Symbol=:RCMK
    )
    return renumberNodes!([material]; method=method)
end

function ensureContiguousNodeTags!()
    nodeTags, _, _ = gmsh.model.mesh.getNodes()

    isempty(nodeTags) &&
        error("ensureContiguousNodeTags!: the current Gmsh model contains no nodes.")

    n = length(nodeTags)

    if minimum(nodeTags) != 1 || maximum(nodeTags) != n
        oldTags = sort(nodeTags)
        newTags = collect(1:n)

        gmsh.model.mesh.renumberNodes(oldTags, newTags)

        @info "Non-contiguous Gmsh node tags detected. Nodes were automatically renumbered to 1:$n."
    end

    return nothing
end

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
    ApproximationData()

Mutable cache attached to a [`Problem`](@ref) for solver-space approximation
transformations. `T` maps solver-space coefficients to the ordinary global
LowLevelFEM/Gmsh storage space, while `R` maps representable global fields
back to the solver space. `metadata` stores basis-specific auxiliary data
(e.g. GLL topology for spectral elements).

For the ordinary full-order Lagrange basis all three fields remain `nothing`,
so no identity matrix is stored.
"""
mutable struct ApproximationData
    T::Union{Nothing,SparseMatrixCSC{Float64,Int}}
    R::Union{Nothing,SparseMatrixCSC{Float64,Int}}
    metadata::Any
end

ApproximationData() = ApproximationData(nothing, nothing, nothing)


"""
    Problem(materials; kwargs...)

Defines a finite element problem on the current Gmsh model.

The `Problem` structure collects all information required for assembly,
including geometry, discretization, unknown field structure and associated
materials. It serves as the central object for operator-based formulations.

# Arguments

- `materials::Vector{Material}`  
  List of materials assigned to physical groups in the mesh.

# Keyword Arguments

## Geometry and discretization

- `dim::Int`  
  Spatial dimension of the problem (1, 2 or 3).
- `thickness::Float64`  
  Thickness parameter for 2D formulations.

## Field definition

- `fieldName::Symbol`  
  Name of the unknown field (e.g. `:u`, `:T`, `:p`).
- `rhsName::Symbol`  
  Name of the right-hand side field.
- `reducedOrder::Bool`

  If `true`, the active approximation order is reduced from `p` to `p - 1`,
  while assembly and field storage remain on the original order-`p` mesh.
  This option is independent of the solver basis. Default: `false`.
- `basis::Symbol`

  Solver-space basis. Currently `:lagrange` and `:spectral` are supported.
  The global LowLevelFEM/Gmsh storage representation remains Lagrange in both
  cases. Default: `:lagrange`.

## Problem configuration

- `type::Symbol`  
  Predefined problem setup. Supported values:

  - `:Solid`  
  - `:PlaneStress`  
  - `:PlaneStrain`  
  - `:AxiSymmetric`  
  - `:HeatConduction`  
  - `:PlaneHeatConduction`  
  - `:AxiSymmetricHeatConduction`  
  - `:Truss`  
  - `:ScalarField`  
  - `:VectorField`  
  - `:TensorField`  
  - `:Reynolds`

  These presets define default values for:
  - spatial dimension (`dim`)
  - number of unknown components (`pdim`)
  - field naming (`field`, `rhs_field`)

- `bandwidth::Symbol`  
  Optional node renumbering strategy for bandwidth reduction.  
  Options: `:RCMK`, `:Hilbert`, `:Metis`, `:none` (default).

## Geometry names (used in special cases)

- `nameTopSurface::String`  
  Name of the top surface (used e.g. in Reynolds-type problems)
- `nameVolume::String`  
  Name of the volume region

# Internal Fields

- `dim::Int`  
  Spatial dimension of the mesh
- `pdim::Int`  
  Number of unknown components per node  
  (e.g. 1 = scalar, 2/3 = vector, 9 = tensor)
- `non::Int`  
  Number of mesh nodes
- `geometry::Geometry`  
  Geometry information extracted from Gmsh
- `material::Vector{Material}`  
  Material definitions assigned to physical groups
- `approximation::ApproximationData`
  Cached solver-space transformation and basis-specific auxiliary data.

# Notes

- The problem is defined on the currently active Gmsh model.
- The total number of degrees of freedom is `non * pdim`.
- The structure is compatible with scalar, vector and tensor fields,
  and supports multifield formulations.

# Examples

```julia
# 3D elasticity
mat = Material("solid", E=210e3, ν=0.3)
P = Problem([mat], type=:Solid)

# Scalar field in 2D
P = Problem([mat], type=:ScalarField, dim=2)

# Vector field in 3D
P = Problem([mat], type=:VectorField, dim=3)

# Tensor field
P = Problem([mat], type=:TensorField, dim=3)

# Reynolds-type problem
P = Problem([mat],
    type=:Reynolds,
    nameTopSurface="top",
    nameVolume="fluid"
)
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
    field::Symbol
    rhs_field::Symbol
    reducedOrder::Bool
    basis::Symbol
    approximation::ApproximationData

    Problem() = new()

    function Problem(
        name, type, dim, pdim, material, thickness, non, geometry,
        field::Symbol, rhs_field::Symbol,
        reducedOrder::Bool=false, basis::Symbol=:lagrange
        )
        P = new(
            name, type, dim, pdim, material, thickness, non, geometry,
            field, rhs_field, reducedOrder, basis, ApproximationData()
        )
        _initialize_approximation_data!(P)
        return P
    end

    function Problem(
        name, type, dim, pdim, material, thickness, non, geometry,
        field::Symbol, rhs_field::Symbol,
        reducedOrder::Bool, basis::Symbol, approx::ApproximationData
        )
        P = new(
            name, type, dim, pdim, material, thickness, non, geometry,
            field, rhs_field, reducedOrder, basis, approx
        )
        return P
    end

    function Problem(
        name, type, dim, pdim, material, thickness, non, geometry,
        reducedOrder::Bool=false, basis::Symbol=:lagrange
        )
        P = new(
            name, type, dim, pdim, material, thickness, non, geometry,
            :unknown, :rhs, reducedOrder, basis, ApproximationData()
        )
        _initialize_approximation_data!(P)
        return P
    end

    function Problem(mat; thickness=1.0, type=:Solid, bandwidth=:none,
        nameTopSurface=nothing, nameVolume=nothing, dim::Int=3,
        fieldName::Symbol=:unknown, rhsName::Symbol=:rhs,
        field::Symbol=fieldName, rhs_field::Symbol=rhsName,
        reducedOrder::Bool=false, basis::Symbol=:lagrange)
        if type == :dummy
            return new(
                "dummy", :dummy, 0, 0, mat, 0, 0, Geometry(),
                field, rhs_field, reducedOrder, basis, ApproximationData()
            )
        end
        basis in (:lagrange, :spectral) ||
            error(
                "Problem: basis must be :lagrange or :spectral."
            )

        basis === :spectral && reducedOrder &&
            error(
                "Problem: basis=:spectral combined with reducedOrder=true " *
                "is part of the intended API but is not yet implemented."
            )
        
        pdim = 3
        dim0 = dim

        #if Sys.CPU_THREADS != Threads.nthreads()
        #    @warn "Number of threads($(Threads.nthreads())) ≠ logical threads in CPU($(Sys.CPU_THREADS))."
        #end
        geometry = Geometry()

        if type == :Solid
            dim = 3
            pdim = 3
            field = :u
            rhs_field = :f
        elseif type == :PlaneStress
            dim = 2
            pdim = 2
            field = :u
            rhs_field = :f
        elseif type == :PlaneStrain
            dim = 2
            pdim = 2
            field = :u
            rhs_field = :f
        elseif type == :AxiSymmetric
            dim = 2
            pdim = 2
            field = :u
            rhs_field = :f
        elseif type == :PlaneHeatConduction
            dim = 2
            pdim = 1
            field = :T
            rhs_field = :q
        elseif type == :HeatConduction
            dim = 3
            pdim = 1
            field = :T
            rhs_field = :q
        elseif type == :AxiSymmetricHeatConduction
            dim = 2
            pdim = 1
            field = :T
            rhs_field = :q
        elseif type == :Truss
            dim = 3
            pdim = 3
            field = :u
            rhs_field = :f
        elseif type == :ScalarField
            dim0 < 1 && error("Problem: dimension of a :ScalarField problem must be one, two or three.")
            dim0 > 3 && error("Problem: dimension of a :ScalarField problem must be equal or less than three.")
            dim = dim0
            pdim = 1
            if nameTopSurface !== nothing && nameVolume !== nothing
                geometry = Geometry(nameTopSurface, nameVolume)
                dim = geometry.dim
            end
        elseif type == :VectorField
            dim0 < 1 && error("Problem: dimension of a :VectorField problem must be one, two or three.")
            dim0 == 1 && error("Problem: dimension of a :VectorField problem must be greater than one.")
            dim0 > 3 && error("Problem: dimension of a :VectorField problem must be two or three.")
            dim = dim0
            pdim = dim0
        elseif type == :TensorField
            dim0 < 1 && error("Problem: dimension of a :TensorField problem must be one, two or three.")
            dim0 == 1 && error("Problem: dimension of a :TensorField problem must be greater than one.")
            dim0 > 3 && error("Problem: dimension of a :TensorField problem must be two or three.")
            dim = dim0
            pdim = 9
        elseif type == :Reynolds
            geometry = Geometry(nameTopSurface, nameVolume)
            dim = geometry.dim
            pdim = 1
            field = :p
            rhs_field = :q
        else
            error("Problem type can be: `:Solid`, `:PlaneStress`, `:PlaneStrain`, `:AxiSymmetric`, `:PlaneHeatConduction`, `:HeatConduction`, `:AxiSymmetricHeatConduction`, `:Truss`, `:ScalarField`, `VectorField`.
            Now problem type = $type ????")
        end
        if !isa(mat, Vector)
            error("Problem: materials are not arranged in a vector. Put them in [...]")
        end
        name = gmsh.model.getCurrent()
        ensureContiguousNodeTags!()
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

            if bandwidth != :none
                oldTags, newTags = gmsh.model.mesh.computeRenumbering(bandwidth, elemTags)
                gmsh.model.mesh.renumberNodes(oldTags, newTags)
            end
        end

        nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()
        non = length(nodeTags)
        if nameTopSurface !== nothing && nameVolume !== nothing
            initialize(geometry, mat, non, field, rhs_field)
        end
        P = new(
            name, type, dim, pdim, material, thickness, non, geometry,
            field, rhs_field, reducedOrder, basis, ApproximationData()
        )

        _initialize_approximation_data!(P)

        return P
    end
end

"""
    _has_approximation(P)

Return `true` when `P` uses an intrinsic solver-space approximation
transformation (reduced polynomial order or a non-Lagrange basis).
"""
@inline _has_approximation(P::Problem) =
    P.reducedOrder || P.basis !== :lagrange

"""
    _initialize_approximation_data!(P)

Build and cache the intrinsic approximation transformation of `P` when the
corresponding implementation is already loaded. In a fully loaded
LowLevelFEM module this is called eagerly by the `Problem` constructor.

The guarded lookup also keeps source include order robust: if a specialized
builder is defined later during module loading, `_ensure_approximation_data!`
will initialize the cache before first use.
"""
function _initialize_approximation_data!(P::Problem)
    P.type === :dummy && return P.approximation

    if P.basis === :spectral
        P.reducedOrder &&
            error(
                "Problem: basis=:spectral combined with reducedOrder=true " *
                "is not yet implemented."
            )

        if isdefined(@__MODULE__, :_initialize_spectral_approximation!)
            getfield(@__MODULE__, :_initialize_spectral_approximation!)(P)
        end

    elseif P.basis === :lagrange
        if P.reducedOrder &&
           isdefined(@__MODULE__, :_initialize_reduced_approximation!)
            getfield(@__MODULE__, :_initialize_reduced_approximation!)(P)
        end

    else
        error("Problem: unsupported basis $(P.basis).")
    end

    return P.approximation
end

"""
    _ensure_approximation_data!(P)

Ensure that the intrinsic approximation transformation required by `P` has
been built and cached. Full-order Lagrange fields intentionally keep an empty
cache and need no transformation matrix.
"""
function _ensure_approximation_data!(P::Problem)
    _has_approximation(P) || return P.approximation

    if P.approximation.T === nothing || P.approximation.R === nothing
        _initialize_approximation_data!(P)
    end

    P.approximation.T === nothing &&
        error(
            "Approximation transformation for basis=$(P.basis), " *
            "reducedOrder=$(P.reducedOrder) is not available."
        )

    P.approximation.R === nothing &&
        error(
            "Approximation restriction for basis=$(P.basis), " *
            "reducedOrder=$(P.reducedOrder) is not available."
        )

    return P.approximation
end

"""
    _approximation_matrices(P)

Return the cached intrinsic prolongation and restriction matrices of `P`.
"""
function _approximation_matrices(P::Problem)
    data = _ensure_approximation_data!(P)
    return data.T, data.R
end



"""
    Field

Alias for [`Problem`](@ref).

`Field` can be used instead of `Problem` when defining unknown fields in
multifield formulations. Both names refer to exactly the same type.
"""
const Field = Problem

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
    ndofs(problem::Problem)

Return total number of dofs for a single-field problem.
"""
ndofs(P::Problem) = P.non * P.pdim

