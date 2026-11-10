"""
    fieldError(F::Union{ScalarField,VectorField,TensorField})

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
    resultant(field::VectorField, phName::String) -> Vector{Float64}
    resultant(problem, view_or_field, phName; grad=false, component=:x,
              offsetX=0, offsetY=0, offsetZ=0)

Sum nodal vector entries on the named physical group. The two-argument form
requires nodal input and returns a length-`problem.dim` vector. Only the first
stored column is used; this is suitable for a resultant of assembled nodal
loads, not quadrature of an arbitrary sampled traction. Axisymmetric presets
multiply each selected entry by the local `2πr` factor.

The legacy three-argument form returns a scalar/tuple for vector fields.
For a scalar Gmsh view tag, it instead integrates probed values over the group:
`grad=true` requests the view gradient, `component=:x/:y/:z` selects the
probed component, and `offsetX/Y/Z` shift the probing coordinates. These
keywords apply to the view-integration branch, not to nodal summation.
The legacy `Matrix` branch expects field metadata and is not an ordinary
array-based interface; prefer a `VectorField`.

Use [`integrate`](@ref) for explicit spatial quadrature of a function or
`ScalarField`, selecting `step` there when several columns are present.
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
    is_elementwise = f isa ScalarField && f.a == [;;]
    # In the usual case, the element order in the field is identical to the
    # Gmsh element order. Build a tag-to-field map only as a fallback for
    # reordered or partially overlapping elementwise fields.
    elemToField = Dict{Int,Int}()
    DIM = problem.dim
    b = problem.thickness
    #ncoord2 = zeros(3 * problem.non)
    sum0 = 0.0
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
    for i ∈ 1:length(dimTags)
        dimTag = dimTags[i]
        dim = dimTag[1]
        tag = dimTag[2]
        elementTypes, elementTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
        #nodeTags::Vector{Int64}, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, tag, true, false)
        #ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
        #ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
        #ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
        #@inbounds for (local_idx, nd) in enumerate(nodeTags)
        #    ci = 3*(nd-1)
        #    cn = 3*(local_idx-1)
        #    ncoord2[ci+1] = ncoord[cn+1]
        #    ncoord2[ci+2] = ncoord[cn+2]
        #    ncoord2[ci+3] = ncoord[cn+3]
        #end

        for ii in 1:length(elementTypes)
            et = elementTypes[ii]
            elementName, dim, order1, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementTypes[ii])
            qorder = order <= 0 ? order1 : order
            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(qorder))
            numIntPoints = length(intWeights)
            comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h = reshape(fun, :, numIntPoints)
            #nnet = zeros(Int, length(elementTags[ii]), numNodes)
            # Batched geometric data for the current geometric entity only.
            # The returned data follow the same element order as `getElements`
            # above, so no global element-to-local-index map is needed.
            jacAll, detAll, coordAll = gmsh.model.mesh.getJacobians(et, intPoints, tag)

            numElements = length(elementTags[ii])
            @assert length(detAll) == numElements * numIntPoints
            @assert length(jacAll) == 9 * numElements * numIntPoints
            @assert length(coordAll) == 3 * numElements * numIntPoints

            fieldStart = 0
            if is_elementwise && numElements > 0
                firstElem = elementTags[ii][1]
                firstIdx = findfirst(==(firstElem), f.numElem)

                if firstIdx !== nothing && firstIdx + numElements - 1 <= length(f.numElem)
                    fieldStart = firstIdx
                    @inbounds for k in 1:numElements
                        if f.numElem[fieldStart + k - 1] != elementTags[ii][k]
                            fieldStart = 0
                            break
                        end
                    end
                end

                if fieldStart == 0 && isempty(elemToField)
                    sizehint!(elemToField, length(f.numElem))
                    @inbounds for k in eachindex(f.numElem)
                        elemToField[f.numElem[k]] = k
                    end
                end
            end

            nodes = elemNodeTags[ii]
            @inbounds for l in 1:length(elementTags[ii])
                elem = elementTags[ii][l]
                #for k in 1:numNodes
                #    nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                #end
                offset = (l - 1) * numNodes
                fieldIdx = 0
                if is_elementwise
                    fieldIdx = fieldStart > 0 ? fieldStart + l - 1 : get(elemToField, elem, 0)
                    fieldIdx == 0 && continue
                end

                fieldValues = is_elementwise ? f.A[fieldIdx] : nothing

                gidx = l - 1

                jacOffset = gidx * 9 * numIntPoints
                detOffset = gidx * numIntPoints
                coordOffset = gidx * 3 * numIntPoints

                Jac = reshape(@view(jacAll[jacOffset + 1 : jacOffset + 9 * numIntPoints]), 3, :)

                jacDet = @view detAll[detOffset + 1 : detOffset + numIntPoints]

                coord = @view coordAll[coordOffset + 1 : coordOffset + 3 * numIntPoints]

                #jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                #Jac = reshape(jac, 3, :)
                s1 = 0.0
                #first = (l - 1) * numNodes + 1
                @inbounds for j in 1:numIntPoints
                    x = coord[3j - 2]
                    y = coord[3j - 1]
                    z = coord[3j]
                    #x = 0.0
                    #y = 0.0
                    #z = 0.0
                    #@inbounds for nk in 1:numNodes
                    #    nd    = nodes[first + nk - 1]      # globális csomópont index
                    #    baseg = 3 * (nd - 1)
                    #    hk    = h[nk, j]
                    #    x += hk * ncoord2[baseg + 1]
                    #    y += hk * ncoord2[baseg + 2]
                    #    z += hk * ncoord2[baseg + 3]
                    #end
                    ff = 0.0
                    if f isa Function
                        ff = f(x, y, z)
                    elseif is_elementwise
                        @simd for k in 1:numNodes
                            ff += h[k, j] * fieldValues[k, step]
                        end
                    elseif f isa ScalarField
                        @simd for k in 1:numNodes
                            ff += h[k, j] * f.a[nodes[offset + k], step]
                        end
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

∫(problem::Problem, phName::String, f::Union{Function,ScalarField}; step::Int64=1, order::Int64=0) = integrate(problem, phName, f, step=step, order=order)

"""
    time_integral!(s, t, d; s0=0.0)

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
    rotateNodes(problem::Problem, coordSys::CoordinateSystem)
    rotateNodes(problem::Problem, phName::String, coordSys::CoordinateSystem)

Build the nodal transformation matrix associated with `coordSys`.

The returned [`Transformation`](@ref) maps local vector components to global
Cartesian components at the active nodes. The coordinate system is evaluated
with the same backend used by [`solveField`](@ref) and [`rotateResult`](@ref).

In the two-argument form, the active nodes are determined directly by the
support of the direction fields stored in `coordSys`.

The three-argument form is retained for compatibility with existing
`rotateNodes` workflows. In this form `phName` is checked against the support of
the coordinate system; the physical group and the coordinate-system support
must contain the same mesh nodes.

# Arguments

- `problem::Problem`: finite-element field whose nodal components are transformed.
- `phName::String`: optional Gmsh physical-group name used to verify the active
  node set.
- `coordSys::CoordinateSystem`: local coordinate-system definition.

# Returns

A [`Transformation`](@ref) object `Q` such that

```math
\\mathbf v_\\mathrm{global} = Q\\,\\mathbf v_\\mathrm{local}.
```

# Example

```julia
er = VectorField(U, "body", [(x, y, z) -> x, (x, y, z) -> y, 0.0])
cs = CoordinateSystem(er)
Q = rotateNodes(U, cs)
```
"""
function rotateNodes(
    problem::Problem,
    coordSys::CoordinateSystem
)
    Q, _ = _coordinate_transformation(problem, coordSys)
    return Q
end

function rotateNodes(
    problem::Problem,
    phName::String,
    coordSys::CoordinateSystem
)
    Q, coord_nodes = _coordinate_transformation(problem, coordSys)
    group_nodes = _nodesOnPhysicalGroup(problem, phName)

    coord_nodes == group_nodes ||
        error(
            "rotateNodes: coordinate-system support and physical group " *
            "'$phName' do not contain the same nodes. " *
            "Define the CoordinateSystem direction field(s) on '$phName', " *
            "or use rotateNodes(problem, coordSys)."
        )

    return Q
end

