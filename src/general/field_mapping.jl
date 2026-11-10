function parse_component(sym::Symbol, problem::Problem)
    s = String(sym)
    pf = String(problem.field)

    if !startswith(s, pf)
        return nothing
    end

    comp = s[length(pf)+1:end]   # lehet "", "x", "xy", stb.

    return comp
end

function component_index(problem::Problem, suffix::String)

    if problem.pdim == 1
        return 1

    elseif problem.pdim == 2 || problem.pdim == 3
        if suffix == "x"
            return 1
        elseif suffix == "y"
            return 2
        elseif suffix == "z"
            return 3
        else
            error("Unknown component '$suffix'")
        end

    elseif problem.pdim == 9
        tensor_map = Dict(
            "xx"=>1, "yx"=>2, "zx"=>3,
            "xy"=>4, "yy"=>5, "zy"=>6,
            "xz"=>7, "yz"=>8, "zz"=>9,
            "x"=>1, "y"=>5, "z"=>9
        )

        if haskey(tensor_map, suffix)
            return tensor_map[suffix]
        else
            error("Unknown tensor component '$suffix'")
        end

    else
        error("Unsupported pdim $(problem.pdim)")
    end
end

@inline function nodeToDof(problem::Problem, node::Int, comp::Int)
    return (node - 1) * problem.pdim + comp
end

function nodesOnPhysicalGroup(problem::Problem, phName::String)

    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

    isempty(dimTags) && error("Physical group '$phName' not found.")

    nodes = Int[]

    for dimTag in dimTags
        edim = dimTag[1]
        etag = dimTag[2]

        # csak valódi perem entitás
        if edim < problem.dim || true
            elementTypes, elementTags, nodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for nt in nodeTags
                append!(nodes, nt)
            end
        else
            error("Cannot apply Dirichlet BC on entity of dimension $edim in $(problem.dim)D problem.")
        end
    end

    return unique(nodes)
end

"""
    constrainedDoFs(problem::Problem, supports::Vector{BoundaryCondition})

Returns the serial numbers of constrained degrees of freedom. Support is a vector of boundary conditions given with the function `displacementConstraint`.

Return: `DoFs`

Types:
- `problem`: Problem
- `supports`: Vector{BoundaryCondition}
- `DoFs`: Vector{Int64}
"""
function constrainedDoFs(problem::Problem,
    supports::Vector{BoundaryCondition})

    dofs = Int[]

    for bc in supports

        # multi-field check
        if bc.problem !== nothing && bc.problem !== problem
            continue
        end

        node_ids = nodesOnPhysicalGroup(problem, bc.phName)

        for (sym, val) in bc.values
            comp = String(sym)
            fld = String(problem.field)

            if !startswith(comp, fld)
                continue
            end

            suffix = String(chopprefix(comp, fld))

            cidx = component_index(problem, suffix)

            for node in node_ids
                push!(dofs, (node - 1) * problem.pdim + cidx)
            end
        end
    end

    return unique(dofs)
end
#=
function constrainedDoFs(problem, supports)
    if !isa(supports, Vector)
        error("constrainedDoFs: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    #non = problem.non
    pdim = problem.pdim
    #dof = non * pdim
    cdofs = Int[]
    
    for i in 1:length(supports)
        name = supports[i].phName
        ux = supports[i].ux
        uy = supports[i].uy
        uz = supports[i].uz
        T = supports[i].T
        p = supports[i].p
        sx = supports[i].sx
        sy = supports[i].sy
        sz = supports[i].sz
        sxy = supports[i].sxy
        syz = supports[i].syz
        szx = supports[i].szx
        
        if (T !== nothing && (ux !== nothing || uy !== nothing || uz !== nothing)) ||
            (T !== nothing && (sx !== nothing || sy !== nothing || sz !== nothing ||
            sxy !== nothing || syz !== nothing || szx !== nothing)) ||
            ((ux !== nothing || uy !== nothing || uz !== nothing) &&
            (sx !== nothing || sy !== nothing || sz !== nothing ||
            sxy !== nothing || syz !== nothing || szx !== nothing))
            
            error("BoundaryCondition: T, displacement and stress BCs are mutually exclusive.")
        end

        phg = getTagForPhysicalName(name)
        nodeTags::Vector{Int64}, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        nodeTagsX = Int[]
        nodeTagsY = Int[]
        nodeTagsZ = Int[]
        nodeTagsXY = Int[]
        nodeTagsYZ = Int[]
        nodeTagsZX = Int[]
        nodeTagsYX = Int[]
        nodeTagsZY = Int[]
        nodeTagsXZ = Int[]
        if T !== nothing && pdim == 1
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
        end
        if p !== nothing && pdim == 1
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
        end
        if ux !== nothing && (pdim == 1 || pdim == 2 || pdim == 3)
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
        end
        if uy !== nothing && (pdim == 2 || pdim == 3)
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
        end
        if pdim == 3 && uz !== nothing
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= pdim
            nodeTagsZ .-= (pdim-3)
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
=#

"""
    allDoFs(problem)

Returns the serial numbers of degrees of freedom of the whole model.

Return: `DoFs`

Types:
- `problem`: Problem
- `DoFs`: Vector{Int64}
"""
function allDoFs(problem)
    nodes0 = Int[]
    nodes = Int[]
    pdim = problem.pdim
    for mat in problem.material
        dimTags = gmsh.model.getEntitiesForPhysicalName(mat.phName)
        for (edim, etag) in dimTags
            nodeTags, _, _ = gmsh.model.mesh.getNodes(edim, etag, true)
            nodes0 = unique(nodes0 ∪ nodeTags)
        end
    end
    for p in 1:pdim
        nodes = unique(nodes ∪ (pdim * nodes0 .- (pdim - p)))
    end
    return nodes
end

"""
    freeDoFs(problem, supports)

Returns the serial numbers of unconstrained degrees of freedom. Support is 
a vector of boundary conditions given with the function `displacementConstraint`.

Return: `DoFs`

Types:
- `problem`: Problem
- `supports`: Vector{BoundaryCondition}
- `DoFs`: Vector{Int64}
"""
function freeDoFs(problem, supports)
    cdofs = constrainedDoFs(problem, supports)
    adofs = allDoFs(problem)
    fdofs = setdiff(adofs, cdofs)
    return fdofs
end

#=
"""
    elementsToNodes(T)

Solves the nodal results `F` from the elemental results `T`.
`T` can be ScalarField, VectorField or TensorField.

Return: `F`

Types:
- `T`: ScalarField, VectorField or TensorField
- `F`: ScalarField, VectorField or TensorField
"""
=#

"""
    elementsToNodes(S::Union{ScalarField,VectorField,TensorField})

Convert an elementwise field to a nodal field by averaging the values
of all elements connected to each node.

Only elements contained in `S.numElem` contribute to the nodal average.
For discontinuous elementwise fields, this operation smooths jumps across
element boundaries by arithmetic averaging.
"""
function elementsToNodes(
    S::Union{ScalarField,VectorField,TensorField}
    )
    isNodal(S) && return S

    problem = S.model
    gmsh.model.setCurrent(problem.name)

    T = typeof(S)
    nsteps = S.nsteps
    numElem = S.numElem
    σ = S.A
    non = problem.non

    epn =
        S isa ScalarField ? 1 :
        S isa VectorField ? (S.type == :v2D ? 2 : 3) :
        S isa TensorField ? 9 :
        error("elementsToNodes: unsupported field type $(typeof(S)).")

    s = zeros(Float64, non * epn, nsteps)
    pcs = zeros(Int, non)

    isempty(numElem) &&
        return T(
            Matrix{Float64}[],
            s,
            S.t,
            Int[],
            nsteps,
            S.type,
            problem,
            S.kind
        )

    # Map each Gmsh element tag to the corresponding position in S.A.
    elementIndex = Dict{Int,Int}()
    sizehint!(elementIndex, length(numElem))

    for (fieldIndex, elementTag) in pairs(numElem)
        haskey(elementIndex, elementTag) &&
            error(
                "elementsToNodes: duplicate element tag $elementTag " *
                "in the field."
            )

        elementIndex[elementTag] = fieldIndex
    end

    # Retrieve the complete mesh connectivity in a single Gmsh call.
    elementTypes, elementTags, elementNodeTags =
        gmsh.model.mesh.getElements(-1, -1)

    found = 0

    for block in eachindex(elementTypes)
        tags = elementTags[block]
        connectivity = elementNodeTags[block]
        numberOfElements = length(tags)

        numberOfElements == 0 && continue

        numberOfNodes =
            length(connectivity) ÷ numberOfElements

        @inbounds for localElement in 1:numberOfElements
            elementTag = Int(tags[localElement])
            fieldIndex = get(elementIndex, elementTag, 0)

            fieldIndex == 0 && continue

            values = σ[fieldIndex]

            size(values, 1) == epn * numberOfNodes ||
                error(
                    "elementsToNodes: incompatible number of local " *
                    "values for element $elementTag."
                )

            nodeOffset = (localElement - 1) * numberOfNodes

            for localNode in 1:numberOfNodes
                nodeTag =
                    Int(connectivity[nodeOffset + localNode])

                pcs[nodeTag] += 1

                localOffset = (localNode - 1) * epn
                nodalOffset = (nodeTag - 1) * epn

                for step in 1:nsteps
                    for component in 1:epn
                        s[nodalOffset + component, step] +=
                            values[localOffset + component, step]
                    end
                end
            end

            found += 1
        end

        found == length(numElem) && break
    end

    found == length(numElem) ||
        error(
            "elementsToNodes: only $found of $(length(numElem)) " *
            "field elements were found in the current Gmsh model."
        )

    @inbounds for node in 1:non
        count = pcs[node]
        count == 0 && continue

        factor = 1.0 / count
        nodalOffset = (node - 1) * epn

        for step in 1:nsteps
            for component in 1:epn
                s[nodalOffset + component, step] *= factor
            end
        end
    end

    return T(
        Matrix{Float64}[],
        s,
        S.t,
        Int[],
        nsteps,
        S.type,
        problem,
        S.kind
    )
end

"""
    nodesToElements(T::Union{ScalarField,VectorField,TensorField}; 
                    onPhysicalGroup="")

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
function nodesToElements(
    r::Union{ScalarField,VectorField,TensorField};
    onPhysicalGroup=""
    )
    if isElementwise(r)
        return r
    end

    problem = r.model
    gmsh.model.setCurrent(problem.name)

    T = typeof(r)

    ncomp =
        r isa ScalarField ? 1 :
        r isa VectorField ? (r.type == :v2D ? 2 : 3) :
        9

    nsteps = r.nsteps

    ε = Vector{Matrix{Float64}}()
    numElem = Int[]

    # r.a follows the ordering returned by getNodes(), while element
    # connectivity contains Gmsh node tags. Build the mapping only when
    # node tags cannot be used directly as indices.
    nodeTags, _, _ = gmsh.model.mesh.getNodes()

    contiguousNodeTags =
        length(nodeTags) == problem.non &&
        all(i -> nodeTags[i] == i, eachindex(nodeTags))

    nodeIndex = contiguousNodeTags ?
        nothing :
        Dict{Int,Int}(Int(tag) => i for (i, tag) in enumerate(nodeTags))

    physicalGroups = onPhysicalGroup == "" ?
        (material.phName for material in problem.material) :
        (onPhysicalGroup,)

    for phName in physicalGroups
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for (edim, etag) in dimTags
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for i in eachindex(elemTypes)
                nelem = length(elemTags[i])
                nelem == 0 && continue

                connectivity = elemNodeTags[i]
                numNodes = length(connectivity) ÷ nelem

                sizehint!(ε, length(ε) + nelem)
                sizehint!(numElem, length(numElem) + nelem)

                for j in 1:nelem
                    push!(numElem, Int(elemTags[i][j]))

                    e = Matrix{Float64}(
                        undef,
                        ncomp * numNodes,
                        nsteps
                    )

                    offset = (j - 1) * numNodes

                    @inbounds for it in 1:nsteps
                        for l in 1:numNodes
                            nodeTag = Int(connectivity[offset + l])

                            idx = contiguousNodeTags ?
                                nodeTag :
                                nodeIndex[nodeTag]

                            localOffset = (l - 1) * ncomp
                            globalOffset = (idx - 1) * ncomp

                            for component in 1:ncomp
                                e[localOffset + component, it] =
                                    r.a[globalOffset + component, it]
                            end
                        end
                    end

                    push!(ε, e)
                end
            end
        end
    end

    return T(
        ε,
        [;;],
        r.t,
        numElem,
        nsteps,
        r.type,
        problem,
        r.kind
    )
end

"""
    elementsToElements(
        S::Union{ScalarField,VectorField,TensorField};
        onPhysicalGroup::String,
        fromPhysicalGroup::String=""
    )

Transfer or restrict an elementwise field to the elements of another physical
group while preserving elementwise discontinuities.

Two operations are supported:

1. same-dimensional restriction,
2. transfer to codimension-one boundary or interface elements.

No nodal averaging is performed.

# Arguments

- `S`:
  Elementwise `ScalarField`, `VectorField` or `TensorField`.

# Keyword Arguments

- `onPhysicalGroup::String`:
  Target Gmsh physical group.

- `fromPhysicalGroup::String=""`:
  Optional source physical group used to disambiguate codimension-one
  transfers on internal interfaces.

# Same-dimensional restriction

If the target physical group has the same dimension as the source elements,
only source elements belonging to the target group are retained.

Supported cases are

```text
1D -> 1D
2D -> 2D
3D -> 3D
```

The element values are copied directly and their elementwise representation  
is preserved.

This is useful, for example, when an elementwise result is defined on an  
entire beam structure but only one beam segment should be extracted.

# Codimension-one transfer

If the target physical group has dimension one less than the source elements,  
each target element is matched to a side of a source element through common  
primary Gmsh node tags.

Supported cases are

```text
1D -> 0D
2D -> 1D
3D -> 2D
```

Values are copied from the matched source element without nodal averaging.  
Target nodes that are not present in the matched source element are assigned  
zero.

If a target element is adjacent to more than one eligible source element,  
the transfer is ambiguous and an error is thrown. Use `fromPhysicalGroup`  
to select the required source side of an internal interface.

# Returns

A field of the same type as `S` containing values on the target elements.

# Examples

Restrict an elementwise beam result to one beam segment:

```julia
M_ABC = elementsToElements(
    M;
    onPhysicalGroup="ABC"
)
```

Transfer a volume field to a boundary surface:

```julia
surface_field = elementsToElements(
    volume_field;
    onPhysicalGroup="surface"
)
```

Select one side of an internal interface:

```julia
interface_field = elementsToElements(
    volume_field;
    onPhysicalGroup="interface",
    fromPhysicalGroup="body_1"
)
```

# Notes

The source field must be elementwise.

For codimension-one transfer, the source and target meshes must be conforming.

To obtain nodally averaged values instead of preserving elementwise  
discontinuities, use

```julia
nodesToElements(
    elementsToNodes(S);
    onPhysicalGroup="..."
)
```

"""
function elementsToElements(
    S::Union{ScalarField,VectorField,TensorField};
    onPhysicalGroup::String,
    fromPhysicalGroup::String=""
    )
    isElementwise(S) ||
        error("elementsToElements: the source field must be elementwise.")

    isempty(onPhysicalGroup) &&
        error("elementsToElements: onPhysicalGroup must not be empty.")

    length(S.A) == length(S.numElem) ||
        error(
            "elementsToElements: the number of local value arrays does not " *
            "match the number of source element tags."
        )

    isempty(S.numElem) &&
        error("elementsToElements: the source field contains no elements.")

    problem = S.model
    gmsh.model.setCurrent(problem.name)

    ncomp =
        S isa ScalarField ? 1 :
        S isa VectorField && S.type == :v2D ? 2 :
        S isa VectorField && S.type == :v3D ? 3 :
        S isa TensorField ? 9 :
        error(
            "elementsToElements: unsupported field type $(typeof(S)) " *
            "with type $(S.type)."
        )

    # Element tags are the stable link between the field data and the Gmsh
    # connectivity. No assumption is made about the order of S.numElem.
    sourceIndex = Dict{Int,Int}()
    sizehint!(sourceIndex, length(S.numElem))

    for (fieldIndex, elementTag0) in pairs(S.numElem)
        elementTag = Int(elementTag0)
        haskey(sourceIndex, elementTag) &&
            error(
                "elementsToElements: duplicate source element tag " *
                "$elementTag."
            )
        sourceIndex[elementTag] = fieldIndex
    end

    nsource = length(S.numElem)
    # Keep indices into the bulk Gmsh connectivity instead of allocating a
    # separate node-tag vector for every source element.
    sourceBlock = zeros(Int, nsource)
    sourceLocalElement = zeros(Int, nsource)
    sourceNumberOfNodes = zeros(Int, nsource)
    sourceDimension = fill(-1, nsource)
    sourceElementTypes = Set{Int}()
    foundSourceElements = 0

    elementTypes, elementTags, elementNodeTags =
        gmsh.model.mesh.getElements(-1, -1)

    for block in eachindex(elementTypes)
        elementType = Int(elementTypes[block])
        tags = elementTags[block]
        connectivity = elementNodeTags[block]
        numberOfElements = length(tags)
        numberOfElements == 0 && continue

        _, elementDimension, _, numberOfNodes, _, _ =
            gmsh.model.mesh.getElementProperties(elementType)

        @inbounds for localElement in 1:numberOfElements
            elementTag = Int(tags[localElement])
            fieldIndex = get(sourceIndex, elementTag, 0)
            fieldIndex == 0 && continue

            firstNode = (localElement - 1) * numberOfNodes + 1
            lastNode = localElement * numberOfNodes
            sourceBlock[fieldIndex] = block
            sourceLocalElement[fieldIndex] = localElement
            sourceNumberOfNodes[fieldIndex] = numberOfNodes
            sourceDimension[fieldIndex] = elementDimension
            push!(sourceElementTypes, elementType)
            foundSourceElements += 1

            values = S.A[fieldIndex]
            size(values, 1) == ncomp * numberOfNodes ||
                error(
                    "elementsToElements: incompatible local data size for " *
                    "source element $elementTag."
                )
            size(values, 2) == S.nsteps ||
                error(
                    "elementsToElements: incompatible number of steps for " *
                    "source element $elementTag."
                )
        end
    end

    foundSourceElements == nsource ||
        error(
            "elementsToElements: only $foundSourceElements of $nsource " *
            "source elements were found in the current Gmsh model."
        )

    sourceDim = sourceDimension[1]
    all(==(sourceDim), sourceDimension) ||
        error(
            "elementsToElements: all source elements must have the same " *
            "dimension."
        )

    sourceDim in (1, 2, 3) ||
        error(
            "elementsToElements: source element dimension $sourceDim is " *
            "not supported."
        )

    # fromPhysicalGroup restricts the source candidates before sides are
    # indexed. This resolves the two-sided trace on internal interfaces.
    allowedSource = trues(nsource)

    if !isempty(fromPhysicalGroup)
        fill!(allowedSource, false)
        fromEntities =
            gmsh.model.getEntitiesForPhysicalName(fromPhysicalGroup)

        isempty(fromEntities) &&
            error(
                "elementsToElements: source physical group " *
                "\"$fromPhysicalGroup\" contains no entities."
            )

        fromDimensions = unique(first.(fromEntities))
        fromDimensions == [sourceDim] ||
            error(
                "elementsToElements: source physical group " *
                "\"$fromPhysicalGroup\" must have dimension $sourceDim."
            )

        for (entityDimension, entityTag) in fromEntities
            _, tagsByType, _ =
                gmsh.model.mesh.getElements(entityDimension, entityTag)

            for tags in tagsByType
                for elementTag0 in tags
                    fieldIndex =
                        get(sourceIndex, Int(elementTag0), 0)
                    fieldIndex == 0 && continue
                    allowedSource[fieldIndex] = true
                end
            end
        end

        any(allowedSource) ||
            error(
                "elementsToElements: source physical group " *
                "\"$fromPhysicalGroup\" contains none of the elements " *
                "stored in the source field."
            )
    end

    @inline function sideKey(nodes)
        numberOfNodes = length(nodes)
        numberOfNodes in 1:4 ||
            error(
                "elementsToElements: only point, line, triangular and " *
                "quadrangular sides are supported."
            )

        a = Int(nodes[1])
        numberOfNodes == 1 && return (1, a, 0, 0, 0)

        b = Int(nodes[2])
        a > b && ((a, b) = (b, a))
        numberOfNodes == 2 && return (2, a, b, 0, 0)

        c = Int(nodes[3])
        b > c && ((b, c) = (c, b))
        a > b && ((a, b) = (b, a))
        numberOfNodes == 3 && return (3, a, b, c, 0)

        d = Int(nodes[4])
        c > d && ((c, d) = (d, c))
        b > c && ((b, c) = (c, b))
        a > b && ((a, b) = (b, a))
        return (4, a, b, c, d)
    end

    targetDim = getDimForPhysicalName(onPhysicalGroup)
    targetTag = getTagForPhysicalName(onPhysicalGroup)

    targetEntities =
        gmsh.model.getEntitiesForPhysicalGroup(targetDim, targetTag)

    isempty(targetEntities) &&
        error(
            "elementsToElements: target physical group " *
            "\"$onPhysicalGroup\" contains no entities."
        )

    #targetDimensions = unique(first.(targetEntities))

    # ------------------------------------------------------------------
    # Same-dimensional restriction: 1D→1D, 2D→2D, 3D→3D
    # ------------------------------------------------------------------
    if targetDim == sourceDim

        Aout = Matrix{Float64}[]
        targetElementTags = Int[]

        for entityTag in targetEntities
            _, tagsByType, _ =
                gmsh.model.mesh.getElements(targetDim, entityTag)

            for tags in tagsByType
                for elementTag0 in tags
                    elementTag = Int(elementTag0)

                    fieldIndex = get(sourceIndex, elementTag, 0)
                    fieldIndex == 0 && continue
                    allowedSource[fieldIndex] || continue

                    push!(Aout, copy(S.A[fieldIndex]))
                    push!(targetElementTags, elementTag)
                end
            end
        end

        isempty(targetElementTags) &&
            error(
                "elementsToElements: target physical group " *
                "\"$onPhysicalGroup\" contains none of the elements " *
                "stored in the source field."
            )

        return typeof(S)(
            Aout,
            [;;],
            S.t,
            targetElementTags,
            S.nsteps,
            S.type,
            problem
        )
    end
    targetDim in (sourceDim, sourceDim - 1) ||
        error(
            "elementsToElements: target physical group " *
            "\"$onPhysicalGroup\" must have dimension $sourceDim or " *
            "$(sourceDim - 1)."
        )

    # Only sides that occur in the target group need to be retained. This
    # avoids a large dictionary containing every internal source side.
    targetSideKeys = Set{NTuple{5,Int}}()
    for entityTag in targetEntities
        targetDimension = targetDim
        targetTypes, tagsByType, nodeTagsByType =
            gmsh.model.mesh.getElements(targetDimension, entityTag)

        for block in eachindex(targetTypes)
            targetType = Int(targetTypes[block])
            tags = tagsByType[block]
            connectivity = nodeTagsByType[block]
            numberOfElements = length(tags)
            numberOfElements == 0 && continue

            _, _, _, numberOfNodes, _, numberOfPrimaryNodes =
                gmsh.model.mesh.getElementProperties(targetType)
            sizehint!(
                targetSideKeys,
                length(targetSideKeys) + numberOfElements
            )

            @inbounds for localElement in 1:numberOfElements
                firstNode = (localElement - 1) * numberOfNodes + 1
                primaryNodes = view(
                    connectivity,
                    firstNode:firstNode + numberOfPrimaryNodes - 1
                )
                push!(targetSideKeys, sideKey(primaryNodes))
            end
        end
    end

    # A positive value identifies the unique source element. Zero marks an
    # ambiguous side shared by multiple eligible source elements.
    sideSources = Dict{NTuple{5,Int},Int}()

    function addSide!(nodes, fieldIndex)
        key = sideKey(nodes)
        key in targetSideKeys || return nothing
        previous = get(sideSources, key, -1)

        if previous == -1
            sideSources[key] = fieldIndex
        elseif previous != fieldIndex
            sideSources[key] = 0
        end

        return nothing
    end

    # Gmsh returns side nodes element by element, in the same element order as
    # getElementsByType. Primary nodes are sufficient for identifying a side.
    for elementType in sourceElementTypes
        tagsByType, _ =
            gmsh.model.mesh.getElementsByType(elementType)
        numberOfElements = length(tagsByType)
        numberOfElements == 0 && continue

        if sourceDim == 1
            _, _, _, _, _, numberOfPrimaryNodes =
                gmsh.model.mesh.getElementProperties(elementType)
            numberOfPrimaryNodes == 2 ||
                error(
                    "elementsToElements: unsupported one-dimensional " *
                    "source element type $elementType."
                )

            _, connectivity =
                gmsh.model.mesh.getElementsByType(elementType)
            _, _, _, numberOfNodes, _, _ =
                gmsh.model.mesh.getElementProperties(elementType)

            @inbounds for localElement in 1:numberOfElements
                fieldIndex =
                    get(sourceIndex, Int(tagsByType[localElement]), 0)
                fieldIndex == 0 && continue
                allowedSource[fieldIndex] || continue

                firstNode = (localElement - 1) * numberOfNodes + 1
                addSide!(view(connectivity, firstNode:firstNode), fieldIndex)
                addSide!(
                    view(connectivity, firstNode + 1:firstNode + 1),
                    fieldIndex
                )
            end
        elseif sourceDim == 2
            edgeNodes = gmsh.model.mesh.getElementEdgeNodes(
                elementType,
                -1,
                true
            )

            length(edgeNodes) % (2 * numberOfElements) == 0 ||
                error(
                    "elementsToElements: invalid edge-node data returned " *
                    "by Gmsh for element type $elementType."
                )

            numberOfEdges =
                length(edgeNodes) ÷ (2 * numberOfElements)

            @inbounds for localElement in 1:numberOfElements
                fieldIndex =
                    get(sourceIndex, Int(tagsByType[localElement]), 0)
                fieldIndex == 0 && continue
                allowedSource[fieldIndex] || continue

                for localEdge in 1:numberOfEdges
                    firstNode =
                        ((localElement - 1) * numberOfEdges +
                         localEdge - 1) * 2 + 1
                    addSide!(
                        view(edgeNodes, firstNode:firstNode + 1),
                        fieldIndex
                    )
                end
            end
        else
            for faceType in (3, 4)
                faceNodes = gmsh.model.mesh.getElementFaceNodes(
                    elementType,
                    faceType,
                    -1,
                    true
                )
                isempty(faceNodes) && continue

                length(faceNodes) % (faceType * numberOfElements) == 0 ||
                    error(
                        "elementsToElements: invalid face-node data " *
                        "returned by Gmsh for element type $elementType."
                    )

                numberOfFaces =
                    length(faceNodes) ÷
                    (faceType * numberOfElements)

                @inbounds for localElement in 1:numberOfElements
                    fieldIndex =
                        get(sourceIndex, Int(tagsByType[localElement]), 0)
                    fieldIndex == 0 && continue
                    allowedSource[fieldIndex] || continue

                    for localFace in 1:numberOfFaces
                        firstNode =
                            ((localElement - 1) * numberOfFaces +
                             localFace - 1) * faceType + 1
                        addSide!(
                            view(
                                faceNodes,
                                firstNode:firstNode + faceType - 1
                            ),
                            fieldIndex
                        )
                    end
                end
            end
        end
    end

    Aout = Matrix{Float64}[]
    targetElementTags = Int[]

    for entityTag in targetEntities
        targetDimension = targetDim
        targetTypes, tagsByType, nodeTagsByType =
            gmsh.model.mesh.getElements(targetDimension, entityTag)

        for block in eachindex(targetTypes)
            targetType = Int(targetTypes[block])
            tags = tagsByType[block]
            connectivity = nodeTagsByType[block]
            numberOfElements = length(tags)
            numberOfElements == 0 && continue

            _, elementDimension, _, numberOfNodes, _, numberOfPrimaryNodes =
                gmsh.model.mesh.getElementProperties(targetType)

            elementDimension == sourceDim - 1 ||
                error(
                    "elementsToElements: unexpected target element " *
                    "dimension $elementDimension."
                )

            sizehint!(Aout, length(Aout) + numberOfElements)
            sizehint!(
                targetElementTags,
                length(targetElementTags) + numberOfElements
            )

            @inbounds for localElement in 1:numberOfElements
                elementTag = Int(tags[localElement])
                firstNode = (localElement - 1) * numberOfNodes + 1
                lastNode = localElement * numberOfNodes
                targetNodes = view(connectivity, firstNode:lastNode)
                key = sideKey(
                    view(targetNodes, 1:numberOfPrimaryNodes)
                )
                fieldIndex = get(sideSources, key, -1)

                if fieldIndex == 0
                    qualifier = isempty(fromPhysicalGroup) ?
                        " Specify fromPhysicalGroup to select the source " *
                        "side of the interface." :
                        " The selected fromPhysicalGroup still contains " *
                        "more than one matching source element."
                    error(
                        "elementsToElements: target element $elementTag " *
                        "matches multiple source elements." *
                        qualifier
                    )
                end

                targetValues = zeros(
                    Float64,
                    ncomp * numberOfNodes,
                    S.nsteps
                )

                if fieldIndex > 0
                    sourceValues = S.A[fieldIndex]
                    sourceConnectivityBlock = sourceBlock[fieldIndex]
                    localSourceElement =
                        sourceLocalElement[fieldIndex]
                    numberOfSourceNodes =
                        sourceNumberOfNodes[fieldIndex]
                    sourceFirstNode =
                        (localSourceElement - 1) * numberOfSourceNodes + 1
                    sourceNodes = view(
                        elementNodeTags[sourceConnectivityBlock],
                        sourceFirstNode:
                        sourceFirstNode + numberOfSourceNodes - 1
                    )

                    for targetLocalNode in 1:numberOfNodes
                        nodeTag = Int(targetNodes[targetLocalNode])
                        sourceLocalNode = 0
                        for localSourceNode in 1:numberOfSourceNodes
                            if Int(sourceNodes[localSourceNode]) == nodeTag
                                sourceLocalNode = localSourceNode
                                break
                            end
                        end
                        sourceLocalNode == 0 && continue

                        targetOffset =
                            (targetLocalNode - 1) * ncomp
                        sourceOffset =
                            (sourceLocalNode - 1) * ncomp
                        for step in 1:S.nsteps
                            for component in 1:ncomp
                                targetValues[
                                    targetOffset + component,
                                    step
                                ] = sourceValues[
                                    sourceOffset + component,
                                    step
                                ]
                            end
                        end
                    end
                end

                push!(Aout, targetValues)
                push!(targetElementTags, elementTag)
            end
        end
    end

    isempty(targetElementTags) &&
        error(
            "elementsToElements: target physical group " *
            "\"$onPhysicalGroup\" contains no mesh elements."
        )

    return typeof(S)(
        Aout,
        [;;],
        S.t,
        targetElementTags,
        S.nsteps,
        S.type,
        problem,
        :field
    )
end

"""
    isNodal(field::Union{ScalarField,VectorField,TensorField})

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
    return !isempty(a.a)
    if a.a != [;;] && a.A == []
        return true
    elseif a.a == [;;] && a.A != []
        return false
    else
        error("isNodal: Scalar field is empty.")
    end
end

"""
    isElementwise(field::Union{ScalarField,VectorField,TensorField})

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
    return !isempty(a.A)
    if a.a == [;;] && a.A != []
        return true
    elseif a.a != [;;] && a.A == []
        return false
    else
        error("isElementwise: Scalar field is empty.")
    end
end

"""
    expandTo3D(V2D::VectorField)

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
    projectTo2D(V3D::VectorField)

Project a 3D vector field onto the xy-plane by discarding the z-component.

return: VectorField

# Examples
```julia
V2D = projectTo2D(V3D)
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

