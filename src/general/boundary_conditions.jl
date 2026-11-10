"""
    CoordinateSystem(e1::VectorField; handedness=:right)
    CoordinateSystem(e1::VectorField, e2::VectorField; handedness=:right)

Define a local orthonormal coordinate system from direction `VectorField`s.

The direction fields may contain constant or spatially varying components and
their finite-element support determines where the local coordinate system is
active. The same `CoordinateSystem` object can be used by [`solveField`](@ref),
[`rotateResult`](@ref), and [`rotateNodes`](@ref).

In 2D, only `e1` is supplied. It defines the first local basis direction. The
second direction is generated automatically. For `handedness=:right`,

```math
\\mathbf e_2 = (-e_{1y}, e_{1x}),
```

while `handedness=:left` reverses this direction.

In 3D, both `e1` and `e2` are required. `e1` is normalized and `e2` defines the
local 1-2 plane. During construction of the transformation, `e2` is
orthogonalized against `e1` and normalized. The third basis direction is

```math
\\mathbf e_3 = \\mathbf e_1 \\times \\mathbf e_2
```

for `handedness=:right`, and the cross-product order is reversed for
`handedness=:left`.

# Arguments

- `e1::VectorField`: first local basis direction.
- `e2::VectorField`: second reference direction; required in 3D and not allowed
  in 2D.
- `handedness::Symbol`: `:right` (default) or `:left`.

# Notes

The direction fields must contain exactly one time step and must belong to the
same `Problem`. Multiple coordinate systems with disjoint support can be
combined in `solveField` and `rotateResult`.

# Examples

```julia
# 2D polar basis on a region called "body"
er = VectorField(U, "body", [(x, y, z) -> x, (x, y, z) -> y, 0.0])
cs = CoordinateSystem(er)

# 3D cylindrical basis
ex = VectorField(U, "body", [(x, y, z) -> x, (x, y, z) -> y, 0.0])
ez = VectorField(U, "body", [0.0, 0.0, 1.0])
cs = CoordinateSystem(ex, ez)

Q = rotateNodes(U, cs)
```
"""
struct CoordinateSystem
    e1::VectorField
    e2::Union{VectorField,Nothing}
    sig::Float64
end

@inline function _coordinate_handedness(handedness::Symbol)
    handedness === :right && return 1.0
    handedness === :left  && return -1.0
    error("CoordinateSystem: handedness must be :right or :left.")
end

function CoordinateSystem(
    e1::VectorField;
    handedness::Symbol=:right
)
    dim = e1.model.dim

    dim == 2 ||
        error(
            "CoordinateSystem: the one-direction form is only valid in 2D; " *
            "a second direction field is required in 3D."
        )

    return CoordinateSystem(
        e1,
        nothing,
        _coordinate_handedness(handedness)
    )
end

function CoordinateSystem(
    e1::VectorField,
    e2::VectorField;
    handedness::Symbol=:right
)
    e1.model === e2.model ||
        error("CoordinateSystem: e1 and e2 must belong to the same Problem.")

    e1.model.dim == 3 ||
        error("CoordinateSystem: in 2D only one direction field must be given.")

    return CoordinateSystem(
        e1,
        e2,
        _coordinate_handedness(handedness)
    )
end

"""
    BoundaryCondition(phName::String; field=nothing, problem=field, kwargs...)

General container type for prescribing boundary conditions in LowLevelFEM.

- `phName` is a Gmsh physical group name.
- `problem` optionally binds the condition to a specific field (required for multi-field).
- `field` alias for `problem`
- `values` stores component-wise boundary condition data, such as `ux`, `uy`, `T`, `exx`, etc.,
  as `Number`, `Function`, or `ScalarField`.

`BoundaryCondition` unifies all Dirichlet- and Neumann-type boundary data
(mechanical, thermal, elastic, convective, etc.) into a single, extensible
structure. A boundary condition always refers to a *physical group* of the mesh
(e.g. `"left"`, `"right"`, `"wall"`) and may prescribe several components of its associated field.

The same `BoundaryCondition` type is used consistently across:
- static and transient solid mechanics
- modal and buckling analysis
- heat conduction and thermo-mechanical problems

Boundary conditions are typically passed as `Vector{BoundaryCondition}` to
assembly, load, and solver routines.

# Prescribed fields

## Mechanical (solid mechanics)

### Dirichlet-type (essential)
- `ux, uy, uz` : prescribed displacement components
- `p`          : prescribed pressure (if applicable)

### Neumann-type (natural)
- `fx, fy, fz` : prescribed surface forces / tractions

### Elastic (Robin-type)
- `kx, ky, kz` : elastic support (spring stiffness per displacement component)

## Thermal

### Dirichlet-type
- `T`          : prescribed temperature

### Neumann-type
- `q`          : prescribed heat flux
- `qn`         : prescribed normal heat flux

### Convection (Robin-type)
- `h`          : heat transfer coefficient
- `T`          : ambient temperature (used together with `h`)

## Poisson-type problems

### Dirichlet-type
- `p`          : prescribed field

### Neumann-type
- `q`          : prescribed force

For heat convection, the field `T` is interpreted as *ambient temperature*.
This is intentional and handled consistently by the heat assembly routines.

# Value types

Each field may be specified as:
- a scalar number (`Number`)
- a `ScalarField`
- a function `(x,y,z) -> value`

depending on the context and the assembly routine using it.

# Usage pattern

It is recommended to:
- use **separate `BoundaryCondition` objects** for conceptually different
  constraints (e.g. displacement support vs. load),
- and collect them in a vector:
```julia
bc = [
    displacementConstraint("left", ux=0.0, uy=0.0),
    load("right", fx=1.0),
    temperatureConstraint("wall", T=300.0)
]
```

# Notes

* Only the fields relevant to a given analysis are accessed; unused fields
  may safely remain `nothing`.
* The container stores keyword/value pairs; it does not validate every
  physical combination. The consuming solver or assembly routine determines
  which components are meaningful.
* For generic DSL fields, component keys follow `P.field` and `P.rhs_field`.
  Use `field=P` in a mixed problem to identify the constrained unknown.
* DSL natural loads are assembled as linear forms with `∫`; storing a load
  in a `BoundaryCondition` alone does not assemble it.
* With local solver coordinates, prescribed components are local. Use
  [`Global`](@ref) to mark a Cartesian right-hand side.

See also:

* `displacementConstraint`
* `load`
* `elasticSupport`
* `temperatureConstraint`
* `heatFlux`
* `heatConvection`
"""
struct BoundaryCondition
    phName::String
    problem::Union{Problem,Nothing}
    values::Dict{Symbol,Union{Number,Function,ScalarField}}

    function BoundaryCondition(phName::String; field=nothing, problem=field, kwargs...)
        vals = Dict{Symbol,Union{Number,Function,ScalarField}}()
        for (k, v) in kwargs
            if v === nothing
                continue
            end
            vals[k] = v
        end
        return new(phName, problem, vals)
    end
end
#=
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
    kx ::Union{Nothing, Number, Function}
    ky ::Union{Nothing, Number, Function}
    kz ::Union{Nothing, Number, Function}
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

    CS::Union{Nothing, CoordinateSystem}
    Q::Union{Nothing,Transformation}

    function BoundaryCondition(phName::String; kwargs...)
        fields = fieldnames(BoundaryCondition)
        values = map(fields) do f
            f == :phName ? phName : get(kwargs, f, nothing)
        end
        return new(values...)
    end
end
=#

"""
    BoundaryConditionFields(bc::BoundaryCondition)

Return the defined field keys of a boundary condition.

# Arguments
- `bc::BoundaryCondition`: Boundary condition object with named values.

# Returns
- An iterator over `Symbol` keys present in `bc.values`.

# Example
```julia
keys = BoundaryConditionFields(bc)
```
"""
function BoundaryConditionFields(bc::BoundaryCondition)
    return keys(bc.values)
end
#=
function BoundaryConditionFields(bc::BoundaryCondition)
    return filter(f -> getfield(bc, f) !== nothing,
                  fieldnames(BoundaryCondition))
end
=#

"""
    LoadCondition(phName; field=nothing, problem=field, kwargs...)

Generic load specification for assembling right-hand side vectors.

- `phName` is a Gmsh physical group name.
- `problem` optionally binds the load to a specific field (required for multi-field).
- `field` alias for `problem`
- `values` stores component-wise load data, such as `fx`, `fy`, `q`, `σxx`, etc.,
  as `Number`, `Function`, or `ScalarField`.

The interpretation (surface traction vs body force vs flux, etc.) is decided
by the assembly routine, not by this struct.
"""
struct LoadCondition
    phName::String
    problem::Union{Problem,Nothing}
    values::Dict{Symbol,Union{Number,Function,ScalarField,Nothing}}

    function LoadCondition(name, prob, vals)
        return new(name, prob, vals)        
    end
    function LoadCondition(phName::String; field=nothing, problem=field, kwargs...)
        vals = Dict{Symbol,Union{Number,Function,ScalarField,Nothing}}()
        for (k, v) in kwargs
            v = v isa ScalarField ? elementsToNodes(v) : v
            vals[k] = v
        end
        return new(phName, problem, vals)
    end
end

function loads_for_problem(loads::Vector{LoadCondition}, P::Problem)
    return filter(lc -> lc.problem === nothing || lc.problem === P, loads)
end

"""
    _extract_rhs_components(problem::Problem, vals::Dict)

Extract RHS components using `problem.rhs_field` prefix.
Returns a Dict{String,Any} mapping component string -> value.
"""
function _extract_rhs_components(problem::Problem, vals::Dict)
    prefix = String(problem.rhs_field)
    comps = Dict{String,Union{Number,Function,ScalarField,Nothing}}()

    for (sym, val) in vals
        s = String(sym)
        if startswith(s, prefix)
            comp = String(chopprefix(s,prefix))
            comps[comp] = val
        end
    end
    return comps
end

"""
    _check_load_keys(problem::Problem, vals::Dict)

Validate load keys. Throws error if invalid combination is detected.
"""
function _check_load_keys(problem::Problem, vals::Dict)

    rhs_comps = _extract_rhs_components(problem, vals)

    p = get(vals, :p, nothing)
    T = get(vals, :T, nothing)
    h = get(vals, :h, nothing)
    qn = get(vals, :qn, nothing)
    hs = get(vals, :hs, nothing)
    q = get(vals, :q, nothing)

    # --- Solid-like problems ---
    if problem.type in (:Solid, :PlaneStress, :PlaneStrain, :AxiSymmetric)

        if p !== nothing && !isempty(rhs_comps)
            error("loadVector: pressure and RHS components cannot be given together.")
        end

        # --- Heat problems ---
    elseif problem.type in (:HeatConduction,
        :PlaneHeatConduction,
        :AxiSymmetricHeatConduction)

        count = 0
        for x in (hs, qn, h, q)
            count += x !== nothing ? 1 : 0
        end

        if count > 1
            error("loadVector: multiple heat-type loads defined in same LoadCondition.")
        end
    end

    # unknown keys
    allowed = Set([:p, :T, :h, :qn, :hs, :q, :T∞])
    for (sym, _) in vals
        s = String(sym)
        if startswith(s, String(problem.rhs_field))
            continue
        elseif sym in allowed
            continue
        else
            error("loadVector: unknown load key '$sym'.")
        end
    end
end

"""
    MultiPointConstraint(; master::String, slave::String,
                         field=nothing, problem=field, kwargs...)
    MPC(; master::String, slave::String,
        field=nothing, problem=field, kwargs...)

Define a multi-point constraint between master and slave physical groups.

`MPC` is a shorthand alias for `MultiPointConstraint`.

The constraint identifies corresponding degrees of freedom on the `master`
and `slave` physical groups. The constrained components can be selected
individually using Boolean keyword arguments.

If no component keyword is specified, all components of the associated field
are constrained.

# Arguments

- `master::String`:
  Name of the Gmsh physical group containing the master node(s).

- `slave::String`:
  Name of the Gmsh physical group containing the slave node(s).

# Keyword Arguments

- `field`:
  Field (`Problem`) to which the constraint belongs. This is required when
  MPCs are used in multifield systems.

- `problem`:
  Backward-compatible alias for `field`.

- `kwargs...`:
  Boolean component selectors. Component names follow the field name and its
  components, for example `ux`, `uy`, `uz`, `φ`, `φx`, `φy`, `φz`.

  A value of `true` constrains the corresponding component, while `false`
  leaves it unconstrained.

# Behavior

Two types of coupling are supported.

## Single-master coupling

If the master physical group contains a single node, all selected slave
degrees of freedom are tied to that master node. This can be used for
remote-point-type constraints.

## Periodic coupling

If the master and slave physical groups contain multiple nodes, the
correspondence is obtained from the periodic node mapping stored in the
Gmsh model.

Periodic pairing must therefore be defined in Gmsh before constructing the
constraint.

# Examples

Tie all displacement components of a boundary to a remote point:

```julia
mpc = MPC(
    master="remote",
    slave="right",
    field=U
)
```

Constrain only selected components:

```julia
mpc = MPC(
    master="remote",
    slave="right",
    field=U,
    ux=true,
    uy=true,
    uz=false
)
```

Periodic displacement coupling:

```julia
periodic = MPC(
    master="right",
    slave="left",
    field=U
)

u = solveField(
    K,
    f;
    support=[bc],
    mpc=[periodic]
)
```

# Notes

MPC relations may be chained. During solution, chained constraints are  
resolved so that each constrained degree of freedom points directly to its  
final master degree of freedom.  
"""
struct MultiPointConstraint
    master::String
    slave::String
    problem::Union{Problem,Nothing}
    components::Dict{Symbol,Bool}

    function MultiPointConstraint(;master::String, slave::String, field=nothing, problem=field, kwargs...)
        comps = Dict{Symbol,Bool}()

        for (k, v) in kwargs
            v === nothing && continue
            v isa Bool || error("MPC: component '$k' must be Bool.")
            comps[k] = v
        end

        return new(master, slave, problem, comps)
    end
end

"""
    MPC

Alias for [`MultiPointConstraint`](@ref).

`MPC` can be used as a shorter name when defining multi-point constraints.
Both names refer to exactly the same type.
"""
const MPC = MultiPointConstraint

function mpcRepresentativeMap(
    problem::Problem,
    mpcs::Vector{MPC}
    )

    ndof = problem.non * problem.pdim

    # ----------------------------------------------------------
    # Union-find structure
    #
    # parent[i] gives the current representative tree.
    # When two equivalence classes are merged, the class
    # containing the master DOF remains the representative.
    # ----------------------------------------------------------

    parent = collect(1:ndof)

    function find_root!(dof::Int)

        root = dof

        # Find root.
        while parent[root] != root
            root = parent[root]
        end

        # Path compression.
        current = dof

        while parent[current] != current
            next = parent[current]
            parent[current] = root
            current = next
        end

        return root
    end

    function union_to_master!(
        masterDof::Int,
        slaveDof::Int
        )

        masterRoot = find_root!(masterDof)
        slaveRoot  = find_root!(slaveDof)

        masterRoot == slaveRoot &&
            return nothing

        # Preserve master-side representative.
        parent[slaveRoot] = masterRoot

        return nothing
    end

    gmsh.model.setCurrent(problem.name)

    comp_map = _bc_component_map(problem)
    prefix = String(problem.field)

    for mpc in mpcs

        # ------------------------------------------------------
        # Multifield filter
        # ------------------------------------------------------

        if mpc.problem !== nothing &&
           mpc.problem !== problem

            continue
        end

        masterTag =
            getTagForPhysicalName(mpc.master)

        slaveTag =
            getTagForPhysicalName(mpc.slave)

        masterNodes, _ =
            gmsh.model.mesh.getNodesForPhysicalGroup(
                -1,
                masterTag
            )

        isempty(masterNodes) &&
            error(
                "MPC: master physical group " *
                "'$(mpc.master)' contains no nodes."
            )

        # ------------------------------------------------------
        # Active components
        # ------------------------------------------------------

        active_components =
            trues(problem.pdim)

        for (sym, active) in mpc.components

            s = String(sym)

            startswith(s, prefix) ||
                error(
                    "MPC: unknown component '$sym'."
                )

            comp =
                String(chopprefix(s, prefix))

            haskey(comp_map, comp) ||
                error(
                    "MPC: invalid component '$comp' " *
                    "for field $(problem.field)."
                )

            active_components[
                comp_map[comp]
            ] = active
        end

        components =
            findall(active_components)

        # ======================================================
        # CASE 1:
        # One master node -> all slave nodes tied to one master
        # ======================================================

        if length(masterNodes) == 1

            slaveNodes, _ =
                gmsh.model.mesh.getNodesForPhysicalGroup(
                    -1,
                    slaveTag
                )

            masterNode =
                Int(masterNodes[1])

            for node in slaveNodes

                slaveNode = Int(node)

                slaveNode == masterNode &&
                    continue

                for comp in components

                    masterDof =
                        problem.pdim *
                        (masterNode - 1) +
                        comp

                    slaveDof =
                        problem.pdim *
                        (slaveNode - 1) +
                        comp

                    union_to_master!(
                        masterDof,
                        slaveDof
                    )
                end
            end

            continue
        end

        # ======================================================
        # CASE 2:
        # Multiple master nodes -> Gmsh periodic pairing
        # ======================================================

        bdim = problem.dim - 1

        masterEntities =
            Set(
                gmsh.model.getEntitiesForPhysicalGroup(
                    bdim,
                    masterTag
                )
            )

        slaveEntities =
            gmsh.model.getEntitiesForPhysicalGroup(
                bdim,
                slaveTag
            )

        isempty(slaveEntities) &&
            error(
                "MPC: periodic slave physical group " *
                "'$(mpc.slave)' contains no boundary entities."
            )

        for slaveEntity in slaveEntities

            tagMaster,
            slaveNodes,
            periodicMasterNodes,
            affineTransform =
                gmsh.model.mesh.getPeriodicNodes(
                    bdim,
                    slaveEntity,
                    true
                )

            tagMaster < 0 &&
                error(
                    "MPC: physical group '$(mpc.slave)' " *
                    "is not defined as periodic in Gmsh."
                )

            tagMaster in masterEntities ||
                error(
                    "MPC: Gmsh periodic master of " *
                    "'$(mpc.slave)' does not belong to " *
                    "'$(mpc.master)'."
                )

            length(slaveNodes) ==
            length(periodicMasterNodes) ||
                error(
                    "MPC: inconsistent periodic node mapping."
                )

            for k in eachindex(slaveNodes)

                slaveNode =
                    Int(slaveNodes[k])

                masterNode =
                    Int(periodicMasterNodes[k])

                for comp in components

                    slaveDof =
                        problem.pdim *
                        (slaveNode - 1) +
                        comp

                    masterDof =
                        problem.pdim *
                        (masterNode - 1) +
                        comp

                    union_to_master!(
                        masterDof,
                        slaveDof
                    )
                end
            end
        end
    end

    # ----------------------------------------------------------
    # Convert union-find structure to the canonical representative
    # map expected by mpcTransformation().
    # ----------------------------------------------------------

    rep = Vector{Int}(undef, ndof)

    for dof in 1:ndof
        rep[dof] = find_root!(dof)
    end

    return rep
end

#=
function mpcRepresentativeMap(
    problem::Problem,
    mpcs::Vector{MPC}
    )

    ndof = problem.non * problem.pdim
    rep = collect(1:ndof)

    gmsh.model.setCurrent(problem.name)

    comp_map = _bc_component_map(problem)
    prefix = String(problem.field)

    for mpc in mpcs

        # ------------------------------------------------------
        # Multifield filter
        # ------------------------------------------------------

        if mpc.problem !== nothing &&
           mpc.problem !== problem

            continue
        end

        masterTag =
            getTagForPhysicalName(mpc.master)

        slaveTag =
            getTagForPhysicalName(mpc.slave)

        masterNodes, _ =
            gmsh.model.mesh.getNodesForPhysicalGroup(
                -1,
                masterTag
            )

        # ------------------------------------------------------
        # Active components
        # ------------------------------------------------------

        active_components =
            trues(problem.pdim)

        for (sym, active) in mpc.components

            s = String(sym)

            startswith(s, prefix) ||
                error(
                    "MPC: unknown component '$sym'."
                )

            comp =
                String(chopprefix(s, prefix))

            haskey(comp_map, comp) ||
                error(
                    "MPC: invalid component '$comp' " *
                    "for field $(problem.field)."
                )

            active_components[
                comp_map[comp]
            ] = active
        end

        components =
            findall(active_components)

        # ======================================================
        # CASE 1:
        # One master node -> all slave nodes tied to one master
        # ======================================================

        if length(masterNodes) == 1

            slaveNodes, _ =
                gmsh.model.mesh.getNodesForPhysicalGroup(
                    -1,
                    slaveTag
                )

            masterNode =
                Int(masterNodes[1])

            for node in slaveNodes

                slaveNode = Int(node)

                slaveNode == masterNode &&
                    continue

                for comp in components

                    masterDof =
                        problem.pdim *
                        (masterNode - 1) +
                        comp

                    slaveDof =
                        problem.pdim *
                        (slaveNode - 1) +
                        comp

                    rep[slaveDof] =
                        masterDof
                end
            end

            continue
        end

        # ======================================================
        # CASE 2:
        # Multiple master nodes -> Gmsh periodic pairing
        # ======================================================

        bdim = problem.dim - 1

        masterEntities =
            Set(
                gmsh.model.getEntitiesForPhysicalGroup(
                    bdim,
                    masterTag
                )
            )

        slaveEntities =
            gmsh.model.getEntitiesForPhysicalGroup(
                bdim,
                slaveTag
            )

        isempty(slaveEntities) &&
            error(
                "MPC: periodic slave physical group " *
                "'$(mpc.slave)' contains no boundary entities."
            )

        for slaveEntity in slaveEntities

            tagMaster,
            slaveNodes,
            periodicMasterNodes,
            affineTransform =
                gmsh.model.mesh.getPeriodicNodes(
                    bdim,
                    slaveEntity,
                    true
                )

            tagMaster < 0 &&
                error(
                    "MPC: physical group '$(mpc.slave)' " *
                    "is not defined as periodic in Gmsh."
                )

            tagMaster in masterEntities ||
                error(
                    "MPC: Gmsh periodic master of " *
                    "'$(mpc.slave)' does not belong to " *
                    "'$(mpc.master)'."
                )

            length(slaveNodes) ==
            length(periodicMasterNodes) ||
                error(
                    "MPC: inconsistent periodic node mapping."
                )

            for k in eachindex(slaveNodes)

                slaveNode =
                    Int(slaveNodes[k])

                masterNode =
                    Int(periodicMasterNodes[k])

                for comp in components

                    slaveDof =
                        problem.pdim *
                        (slaveNode - 1) +
                        comp

                    masterDof =
                        problem.pdim *
                        (masterNode - 1) +
                        comp

                    # Detect contradictory MPC definitions.
                    if rep[slaveDof] != slaveDof &&
                       rep[slaveDof] != masterDof

                        error(
                            "MPC: DOF $slaveDof has " *
                            "conflicting master DOFs."
                        )
                    end

                    rep[slaveDof] =
                        masterDof
                end
            end
        end
    end

    _canonicalize_mpc_rep!(rep)

    return rep
end
=#

function mpcTransformation(rep::Vector{Int})

    p = unique(rep)

    representative_to_reduced = zeros(Int, length(rep))

    for (j, dof) in enumerate(p)
        representative_to_reduced[dof] = j
    end

    full_to_reduced =
        [representative_to_reduced[rep[i]] for i in eachindex(rep)]

    I = collect(eachindex(rep))
    J = full_to_reduced
    V = ones(Float64, length(rep))

    T = sparse(I, J, V, length(rep), length(p))

    return T, p, full_to_reduced
end

function mpcReducedBCData(full_to_reduced, fixed, uD)

    nred = maximum(full_to_reduced)

    fixed_r = Int[]
    uD_r = zeros(Float64, nred)
    assigned = falses(nred)

    for i in fixed

        j = full_to_reduced[i]
        value = uD[i]

        if assigned[j]

            isapprox(uD_r[j], value) ||
                error(
                    "MPC: incompatible prescribed values on tied DOFs."
                )

        else

            uD_r[j] = value
            assigned[j] = true
            push!(fixed_r, j)

        end
    end

    sort!(unique!(fixed_r))

    free_r = setdiff(1:nred, fixed_r)

    return free_r, fixed_r, uD_r
end

"""
    mpcReducedBCData(full_to_reduced, fixed, uD::AbstractMatrix)

Map time-dependent prescribed full-space values to the MPC-reduced space.
"""
function mpcReducedBCData(
    full_to_reduced,
    fixed,
    uD::AbstractMatrix
    )

    nsteps = size(uD, 2)

    free_r, fixed_r, uD1 =
        mpcReducedBCData(
            full_to_reduced,
            fixed,
            @view(uD[:, 1])
        )

    uD_r =
        zeros(
            Float64,
            length(uD1),
            nsteps
        )

    uD_r[:, 1] .= uD1

    for i in 2:nsteps

        free_i, fixed_i, uDi =
            mpcReducedBCData(
                full_to_reduced,
                fixed,
                @view(uD[:, i])
            )

        free_i == free_r ||
            error(
                "MPC: reduced free DOFs change between time steps."
            )

        fixed_i == fixed_r ||
            error(
                "MPC: reduced constrained DOFs change between time steps."
            )

        uD_r[:, i] .= uDi
    end

    return free_r, fixed_r, uD_r
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
