"""
    openPreProcessor(; openGL=true)

Launches the Gmsh preprocessor window with OpenGL disabled by default.

Returns: nothing

Types:
- `openGL`: Boolean
"""
function openPreProcessor(; openGL=true)
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
    getParameter(name)

Gets the value of the parameter `name`. 

Returns: value

Types:
- `name`: String
- `value`: Float64
"""
function getParameter(name)
    return gmsh.parser.getNumber(name)[1]
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
    probe(A::Union{ScalarField,VectorField,TensorField}, 
          x::Number, y::Number, z::Number; step=Int)

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
        ind = _probe_element_index(A, elementTag)
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
        ind = _probe_element_index(A, elementTag)
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
        ind = _probe_element_index(A, elementTag)
        SS = round(fun' * A.A[ind][1:numNodes, step], digits=10)
    elseif A.A == []
        SS = round(fun' * A.a[nodeTags, step], digits=10)
    end
    return SS
end

function _probe_element_index(A, elementTag)
    ind = findfirst(==(elementTag), A.numElem)

    ind === nothing && error(
        "probe: no elementwise value is available at this point. " *
        "The point may lie on a lower-dimensional boundary element. " *
        "Convert the field to a nodal field with elementsToNodes(...) and try again."
    )

    return ind
end


"""
    probe(A::Union{ScalarField,VectorField,TensorField}, name::String; step=1)

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

function _problem_without_approximation(P::Problem)
    return Problem(
        P.name,
        P.type,
        P.dim,
        P.pdim,
        P.material,
        P.thickness,
        P.non,
        P.geometry,
        P.field,
        P.rhs_field,
        P.reducedOrder,
        P.basis,
        ApproximationData()
    )
end

"""
    saveField(fileName::String, variable::Union{ScalarField,VectorField,TensorField,Number})

Saves `variable` of type ScalarField, VectorField, or TensorField to a file named `fileName`.
The name of the file will be complemented with the string "-LLF-Data.jld2".
The field `kind` is preserved when saving and loading.

Returns: nothing

Types:

- `fileName`: String
- `variable`: ScalarField, VectorField or TensorField
"""
function saveField(
    fileName::String,
    var::Union{ScalarField,VectorField,TensorField}
    )
    P = _problem_without_approximation(var.model)
    T = typeof(var)

    variable = T(
        var.A,
        var.a,
        var.t,
        var.numElem,
        var.nsteps,
        var.type,
        P,
        var.kind
    )

    @save fileName * "-LLF-Data.jld2" variable

    return nothing
end

function saveField(fileName::String, variable::Number)
    @save fileName * "-LLF-Data.jld2" variable
    return nothing
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

