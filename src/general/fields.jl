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

`kind=:field` (default) marks ordinary field data. For `:eigen_real` and
`:eigen_imag`, columns are eigenmodes and `t` contains the corresponding real
or imaginary eigenvalue components rather than physical time instants.


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
    kind::Symbol
    function ScalarField(A0, a0, t0, numElem0, nsteps0, type0, model, kind::Symbol=:field)
        kind in (:field, :eigen_real, :eigen_imag) ||
            throw(ArgumentError("ScalarField: kind must be :field, :eigen_real or :eigen_imag."))
        return new(A0, a0, t0, numElem0, nsteps0, type0, model, kind)
    end
    function ScalarField(problem::Problem, dataField::Vector; steps=1, tmin=0.0, tmax=tmin+(steps-1))
        gmsh.model.setCurrent(problem.name)
        nodeTags, coords, _ = gmsh.model.mesh.getNodes()
        ###################################################
        nodeIndex = Dict{Int,Int}()
        for (i, tag) in enumerate(nodeTags)
            nodeIndex[tag] = i
        end
        ###################################################
        type = :scalar
        nsteps = steps
        t = range(start=tmin, stop=tmax, length=steps)
        A = Vector{Matrix{Float64}}()
        numElem = Int[]

        @inbounds for data in dataField
            #phName, f, fx, fy, fz, fxy, fyz, fzx, fyx, fzy, fxz = data
            phName = data.phName
            vals = data.values
            f = get(vals, :f, nothing)
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
                                    ii = nodeIndex[nodeTag]
                                    if f isa Function
                                        x = coords[3*ii - 2]
                                        y = coords[3*ii - 1]
                                        z = coords[3*ii]
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
        t = s.kind === :field ?
            collect(range(start=tmin, stop=tmax, length=steps)) :
            fill(s.t[step], steps)
        ScalarField(A, [;;], t, s1.numElem, steps, s1.type, s1.model, s.kind)
    end
end

"""
    VectorField(A, a, t, numElem, nsteps, type, model)
    VectorField(problem::Problem, dataField::Vector)
    VectorField(problem::Problem, phName::String, data::Vector)
    VectorField(problem::Problem, phName::String, func::Function)
    VectorField(comps::Vector)

Container for a (possibly time-dependent) vector field defined on a finite
element mesh (e.g. displacement, velocity, heat flux).

A `VectorField` stores a 3D vector field either
- **element-wise**, with values given at element nodes (`A`), or
- **nodally**, with values given at global mesh nodes (`a`).

Time dependence is handled by storing multiple time steps for each spatial
degree of freedom, analogously to `ScalarField`.

---

## Stored data

`kind=:field` (default) marks ordinary field data. For `:eigen_real` and
`:eigen_imag`, columns are eigenmodes and `t` contains the corresponding real
or imaginary eigenvalue components rather than physical time instants.


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
(`x₁,y₁,z₁,x₂,y₂,z₂,…`) for element-wise `:v3D` data. A nodal
solution or assembled load instead has `model.pdim` components per node;
for a 2D vector unknown this means `(x₁,y₁,x₂,y₂,…)`.

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

Assembles a 3D vector field from exactly three components. At least one
component must be a `ScalarField`; the others may be compatible scalar fields
or numbers. Numbers are expanded over the reference field's element support
and columns. Nodal scalar components are converted to element-wise storage.

Requirements:

* exactly three components, at least one a scalar field,
* all field components share the same `Problem`,
* identical element numbering,
* identical time discretization.

Each scalar field provides one vector component.

---

### From an existing VectorField

```julia
VectorField(s; steps, tmin, tmax, step)
```

Creates a new `VectorField` by **replicating a selected time step** of an
existing field.
This is useful for initializing a time-dependent field from a static solution.

* `step` selects the time index of `s` to be replicated.
* The new field contains `steps` identical time slices.

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
    kind::Symbol
    function VectorField(A0, a0, t0, numElem0, nsteps0, type0, model, kind::Symbol=:field)
        kind in (:field, :eigen_real, :eigen_imag) ||
            throw(ArgumentError("VectorField: kind must be :field, :eigen_real or :eigen_imag."))
        return new(A0, a0, t0, numElem0, nsteps0, type0, model, kind)
    end
    function VectorField(comps::AbstractVector)

        length(comps) == 3 ||
            error("VectorField must have exactly 3 components.")

        all(c -> c isa ScalarField || c isa Number, comps) ||
            error(
                "VectorField components must be ScalarField or Number."
            )

        # At least one ScalarField is required to determine the model,
        # element support and time discretization.
        iref = findfirst(c -> c isa ScalarField, comps)

        iref === nothing &&
            error(
                "VectorField: at least one component must be a ScalarField " *
                "when using VectorField([c1, c2, c3])."
            )

        ref = nodesToElements(comps[iref])

        prob = ref.model
        elems = ref.numElem
        nsteps = ref.nsteps
        t = ref.t

        comps1 = Vector{ScalarField}(undef, 3)

        for i in 1:3

            c = comps[i]

            if c isa ScalarField

                s = nodesToElements(c)

                s.model === prob ||
                    error(
                        "VectorField: all ScalarField components must belong " *
                        "to the same Problem."
                    )

                s.numElem == elems ||
                    error(
                        "VectorField: ScalarField element numbering mismatch."
                    )

                s.nsteps == nsteps ||
                    error(
                        "VectorField: ScalarField time-step mismatch."
                    )

                comps1[i] = s

            else
                # Constant component with the same element support and
                # time discretization as the reference ScalarField.
                A = Vector{Matrix{Float64}}(undef, length(elems))

                for e in eachindex(elems)
                    A[e] = fill(
                        Float64(c),
                        size(ref.A[e], 1),
                        nsteps
                    )
                end

                comps1[i] = ScalarField(
                    A,
                    [;;],
                    copy(t),
                    copy(elems),
                    nsteps,
                    :scalar,
                    prob
                )
            end
        end

        # Assemble interleaved vector components.
        A = Vector{Matrix{Float64}}(undef, length(elems))

        for e in eachindex(elems)

            numNodes = size(comps1[1].A[e], 1)

            block = zeros(
                3 * numNodes,
                nsteps
            )

            for c in 1:3
                block[c:3:end, :] .= comps1[c].A[e]
            end

            A[e] = block
        end

        return VectorField(
            A,
            [;;],
            copy(t),
            copy(elems),
            nsteps,
            :v3D,
            prob
        )
    end
#= cab be deleted if the previous constructor is working
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
        return new(copy(ε0), a, comps1[1].t, elems, comps1[1].nsteps, :v3D, prob, :field)
    end
=#
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
            vals = dataField[i].values
            fx = get(vals, :fx, nothing)
            fy = get(vals, :fy, nothing)
            fz = get(vals, :fz, nothing)
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
        return new(A, a, t, numElem, nsteps, type, problem, :field)
    end
    function VectorField(problem::Problem, phName::String, data::Vector)
        @assert length(data) == 3 "VectorField must have exactly 3 components."

        comps = Vector{ScalarField}(undef, 3)

        for i in 1:3
            if data[i] isa ScalarField
                @assert data[i].model === problem "VectorField: ScalarFields must belong to the same problem."
                comps[i] = data[i]

            elseif data[i] isa Number || data[i] isa Function
                comps[i] = ScalarField(problem, phName, data[i])

            else
                error("Unsupported type in VectorField component: $(typeof(data[i]))")
            end
        end

        return VectorField(comps)
    end
    #function VectorField(problem::Problem, phName::String, data::Vector)
    #    if size(data) == (3,)
    #        f = field(phName, fx=data[1], fy=data[2], fz=data[3])
    #        return VectorField(problem, [f])
    #    else
    #        error("VectorField: size of data is $(size(data)).")
    #    end
    #end
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
        return new(A, a, t, numElem, nsteps, type, problem, :field)
    end
    function VectorField(s::VectorField; steps=1, tmin=0.0, tmax=tmin+(steps-1), step=1)
        if step > s.nsteps || step < 1
            error("VectorField: `step` is out of range - 0 < $step < $(s.nsteps)")
        end
        A = Vector{Matrix{Float64}}(undef, length(s.numElem))
        s1 = nodesToElements(s)
        for j = 1:length(s1.numElem)
            A[j] = zeros(size(s1.A[j], 1), steps)
            for i in 1:steps
                A[j][:, i] = s1.A[j][:, step]
            end
        end
        t = s.kind === :field ?
            collect(range(start=tmin, stop=tmax, length=steps)) :
            fill(s.t[step], steps)
        VectorField(A, [;;], t, s1.numElem, steps, s1.type, s1.model, s.kind)
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

`kind=:field` (default) marks ordinary field data. For `:eigen_real` and
`:eigen_imag`, columns are eigenmodes and `t` contains the corresponding real
or imaginary eigenvalue components rather than physical time instants.


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

Tensor components are stored in **column-major, interleaved form**
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

Both examples construct element-wise tensor fields; their values differ.
Matrix-based construction also accepts `ScalarField` entries belonging to the
same `Problem`, and mixed scalar components can be used with
`TensorField(comps::Matrix)`. Nodal components are converted to element-wise
storage; component fields must have compatible support and column counts.
"""
struct TensorField <: AbstractField
    A::Vector{Matrix{Float64}}
    a::Matrix{Float64}
    t::Vector{Float64}
    numElem::Vector{Int}
    nsteps::Int
    type::Symbol
    model::Problem
    kind::Symbol
    function TensorField(A0, a0, t0, numElem0, nsteps0, type0, model, kind::Symbol=:field)
        kind in (:field, :eigen_real, :eigen_imag) ||
            throw(ArgumentError("TensorField: kind must be :field, :eigen_real or :eigen_imag."))
        return new(A0, a0, t0, numElem0, nsteps0, type0, model, kind)
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
        return new(A, a, t, numElem, nsteps, type, problem, :field)
    end
    function TensorField(problem::Problem, phName::String, data::Matrix)
        @assert size(data) == (3,3) "TensorField must be a 3×3 matrix."

        comps = Matrix{ScalarField}(undef, 3, 3)

        for i in 1:3, j in 1:3
            val = data[i,j]

            if val isa ScalarField
                @assert val.model === problem "TensorField: ScalarFields must belong to the same problem."
                comps[i,j] = val

            elseif val isa Number || val isa Function
                comps[i,j] = ScalarField(problem, phName, val)

            else
                error("Unsupported type in TensorField component: $(typeof(val))")
            end
        end

        return TensorField(comps)
    end
    #function TensorField(problem::Problem, phName::String, data::Matrix)
    #    if size(data) == (3,3)
    #        g = [x isa Function ? x : ((_,_,_)->x) for x in data]
    #        f = field(phName, fx=g[1], fy=g[5], fz=g[9], fxy=g[4], fyz=g[8], fzx=g[3], fyx=g[2], fzy=g[6], fxz=g[7])
    #        #f = field(phName, fx=data[1], fy=data[5], fz=data[9], fxy=data[2], fyz=data[6], fzx=data[3])
    #        return TensorField(problem, [f])
    #    else
    #        error("TensorField: size of data is $(size(data)).")
    #    end
    #end
    function TensorField(comps::AbstractMatrix)

        size(comps) == (3, 3) ||
            error("TensorField requires a 3×3 matrix.")

        all(c -> c isa ScalarField || c isa Number, comps) ||
            error("TensorField components must be ScalarField or Number.")

        # At least one ScalarField is required to determine the mesh support.
        iref = findfirst(c -> c isa ScalarField, comps)

        iref === nothing &&
            error(
                "TensorField: at least one component must be a ScalarField " *
                "when using TensorField(matrix)."
            )

        ref = nodesToElements(comps[iref])

        prob    = ref.model
        numElem = ref.numElem
        nsteps  = ref.nsteps
        t       = ref.t

        comps2 = Matrix{ScalarField}(undef, 3, 3)

        @inbounds for i in 1:3, j in 1:3

            c = comps[i, j]

            if c isa ScalarField

                s = nodesToElements(c)

                s.model === prob ||
                    error("TensorField: ScalarFields must share the same Problem.")

                s.numElem == numElem ||
                    error("TensorField: element list mismatch.")

                s.nsteps == nsteps ||
                    error("TensorField: time step mismatch.")

                comps2[i, j] = s

            else
                A = Vector{Matrix{Float64}}(undef, length(numElem))

                for e in eachindex(numElem)
                    A[e] = fill(
                        Float64(c),
                        size(ref.A[e], 1),
                        nsteps
                    )
                end

                comps2[i, j] = ScalarField(
                    A,
                    [;;],
                    copy(t),
                    copy(numElem),
                    nsteps,
                    :scalar,
                    prob
                )
            end
        end

        A = Vector{Matrix{Float64}}(undef, length(numElem))

        for e in eachindex(numElem)

            numNodes = size(comps2[1, 1].A[e], 1)

            block = zeros(
                9 * numNodes,
                nsteps
            )

            @inbounds for i in 1:3, j in 1:3
                block[(j-1)*3+i:9:end, :] .= comps2[i, j].A[e]
            end

            A[e] = block
        end

        return TensorField(
            A,
            [;;],
            copy(t),
            copy(numElem),
            nsteps,
            :tensor,
            prob
        )
    end
    #=
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
    =#
end

"""
    SystemVector

Global right-hand-side container. `.a::AbstractMatrix{Float64}` has one row
per DOF and one column per right-hand side or stored step.

A single-field vector has `.model` set and `.problems === nothing`. A
multifield vector has `.model === nothing`, an ordered `.problems` vector,
and zero-based starting `.offsets` into the concatenated rows. Its field order
must match the corresponding [`SystemMatrix`](@ref).

Construct through `SystemVector(field)` or `SystemVector([f1, f2, ...])`.
"""
struct SystemVector
    a::AbstractMatrix{Float64}

    # ---- Single-field metadata ----
    model::Union{Problem,Nothing}

    # ---- Block metadata ----
    problems::Union{Vector{Problem},Nothing}
    offsets::Union{Vector{Int},Nothing}
end

"""
    SystemVector(field::Union{ScalarField,VectorField,TensorField})

Wrap the nodal RHS matrix `field.a` without dropping any columns. The matrix
is shared with the input; the constructor does not convert element-wise data.
Use [`elementsToNodes`](@ref) first when nodal field values are intended.

For an assembled load from [`∫`](@ref), all columns are retained as independent
right-hand sides. See also [`SystemVector`](@ref).
"""
function SystemVector(field::Union{ScalarField,VectorField,TensorField})

    # assume single load step
    a = field.a

    return SystemVector(a, field.model, nothing, nothing)
end

"""
    SystemVector(fields::Vector)

Concatenate nodal RHS matrices in the supplied field order for a multifield
system. All columns are preserved; fields must have the same column count,
unique `Problem`s, and the correct nodal DOF counts. An empty list is rejected.

The order must match the block order of [`SystemMatrix`](@ref). For example,
`F = SystemVector([fu, fp])` combines already assembled velocity and pressure
loads, and `u, p = solveField(K, F; support=bc)` reconstructs their solutions.
"""
function SystemVector(fields::Vector)

    isempty(fields) && error("SystemVector: empty field list.")

    # ----------------------------------------------------------
    # 1) Extract problems in GIVEN order
    # ----------------------------------------------------------
    problems = Problem[]

    for f in fields
        P = f.model
        push!(problems, P)
    end

    # check uniqueness (no duplicate fields)
    length(unique(problems)) == length(problems) ||
        error("SystemVector: duplicate Problems in field list.")

    # ----------------------------------------------------------
    # 2) Compute offsets (trial order)
    # ----------------------------------------------------------
    offsets = Vector{Int}(undef, length(problems))
    offsets[1] = 0

    for i in 2:length(problems)
        offsets[i] = offsets[i-1] + ndofs(problems[i-1])
    end

    total_dofs = offsets[end] + ndofs(problems[end])
    nsteps = size((fields[1]).a, 2)
    for f in fields
        size(f.a, 2) == nsteps ||
            error("SystemVector: all fields must have the same number of steps.")
    end
    a_global = zeros(total_dofs, nsteps)

    # ----------------------------------------------------------
    # 3) Fill global vector
    # ----------------------------------------------------------
    for (idx, f) in enumerate(fields)

        P = f.model
        offset = offsets[idx]
        nloc = ndofs(P)

        # single load step assumption
        a_local = f.a

        size(a_local, 1) == nloc ||
            error("SystemVector: local RHS size mismatch for problem $(P.name).")

        a_global[offset+1:offset+nloc, :] .= a_local
    end

    return SystemVector(a_global, nothing, problems, offsets)
end

