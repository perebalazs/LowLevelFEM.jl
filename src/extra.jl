export pressureConstraint, flowRate
export solvePressure, pressureInVolume, solveShearStress


function showGapThickness(phName; name=phName, visible=false)
    SS = gmsh.view.add(name)
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
    nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(-1, -1, true, false)
    ncoord2 = similar(ncoord)
    ncoord2[nodeTags*3 .- 2] .= ncoord[1:3:length(ncoord)]
    ncoord2[nodeTags*3 .- 1] .= ncoord[2:3:length(ncoord)]
    ncoord2[nodeTags*3 .- 0] .= ncoord[3:3:length(ncoord)]
    ret2 = []
    ret3 = []
    ret4 = []
    c = 0
    for idm in 1:length(dimTags)
        dimTag = dimTags[idm]
        edim = dimTag[1]
        etag = dimTag[2]
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
        @inbounds for i in 1:length(elemTags)
            c += length(elemTags[i])
        end
    end
    sizehint!(ret2, c)
    sizehint!(ret3, c)
    sizehint!(ret4, c)

    el2 = 0
    el3 = 0
    el4 = 0
    for idm in 1:length(dimTags)
        dimTag = dimTags[idm]
        edim = dimTag[1]
        etag = dimTag[2]
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
        for i in 1:length(elemTypes)
            et = elemTypes[i]
            elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
            nn = zeros(numPrimaryNodes * 4)
            @inbounds for j in 1:length(elemTags[i])
                #elem = elemTags[i][j]
                @inbounds for l in 1:numPrimaryNodes
                    for k in 1:2
                        nn[k*numPrimaryNodes-numPrimaryNodes + l] = ncoord2[elemNodeTags[i][j*numNodes-(numNodes-l)]*3-(3-k)]
                    end
                    k = 3
                    nn[k*numPrimaryNodes-numPrimaryNodes + l] = 0
                    k = 4
                    k2 = 3
                    nn[4*numPrimaryNodes-numPrimaryNodes + l] = ncoord2[elemNodeTags[i][j*numNodes-(numNodes-l)]*3-(3-k2)]
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

function pressureConstraint(name::String; p=1im)
    bc0 = name, p, 1im, 1im
    return bc0
end

function flowRate(name::String; q=0)
    ld0 = name, q
    return ld0
end

function initialize(problem::Problem)
    h0 = ScalarField(problem, problem.geometry.nameGap, (x, y, z) -> z)
    h0 = elementsToNodes(h0)
    h = projectScalarField(h0, from=problem.geometry.nameGap, to=problem.material[1].phName, gap=true)
    h0 = nodesToElements(h, onPhysicalGroup=problem.material[1].phName)
    dhdx = elementsToNodes(∂x(h0))
    dhdx = elementsToNodes(dhdx)
    problem.geometry.h = h
    problem.geometry.dhdx = dhdx
end

"""
    projectScalarField(p::ScalarField; from="", to="", gap=false, binSize=0.7)

Projects a scalar field defined on a 2D surface onto the nodes of another mesh
(typically a 3D volume or a parallel surface), using a robust AABB-based
spatial binning algorithm.

This function is primarily intended for use in the Reynolds equation solver,
where a pressure field defined on a lubricating surface must be transferred
to another surface or to the surrounding 3D mesh.

---

## Functionality

The projection is performed by:
1. Collecting all 2D finite elements belonging to the physical group `from`.
2. Assigning each source element to all spatial bins intersected by its
   axis-aligned bounding box (AABB) in the *(x,y)* plane.
3. Assigning target nodes (from the physical group `to`) to spatial bins
   based on their *(x,y)* coordinates.
4. For each target node, testing candidate source elements from the same bin:
   - computing local coordinates using Gmsh,
   - performing an inside-element test (triangular or quadrilateral),
   - interpolating the scalar field using Gmsh-provided Lagrange basis functions.

The algorithm is deterministic, robust with respect to bin size, and works
for arbitrary element order supported by Gmsh.

---

## Arguments

- `p::ScalarField`  
  Scalar field defined on the source surface (`from`).
  In the Reynolds context, this typically represents the pressure field.

---

## Keyword Arguments

- `from::String`  
  Name of the physical group defining the **source 2D surface**
  on which the scalar field `p` is defined.

- `to::String`  
  Name of the physical group defining the **target mesh**
  (surface or volume) whose nodes receive the projected values.

- `gap::Bool = false`  
  If `true`, the source surface nodes are temporarily projected onto the
  *z = 0* plane during the projection.
  
  This option is useful in Reynolds-type lubrication problems, where:
  - the height of the lubricat is defined on a geometrically curved upper surface,
  - but the projection should be performed purely in *(x,y)* coordinates.
  
  After the projection, the original node coordinates are restored.

- `binSize::Real = 0.7`  
  Controls the spatial bin size relative to the smallest source element size.
  
  Smaller values increase the number of bins (finer spatial partitioning),
  potentially improving performance for large meshes.
  
  Due to the AABB-based binning strategy, **correctness does not depend on
  the bin size**; this parameter only affects performance.

---

## Returns

- `ScalarField`  
  A new scalar field defined on the target mesh nodes, containing the
  interpolated values of `p`.

---

## Notes

- This function relies exclusively on the Gmsh API for element localization
  and basis function evaluation; no custom shape functions are implemented.
- Both triangular and quadrilateral elements are supported, for arbitrary
  polynomial order.
- The AABB-based binning guarantees that no valid node–element associations
  are missed, even for small bin sizes or highly distorted elements.
- The function is optimized for single-threaded execution; parallel scaling
  may be limited by Gmsh API overhead.

---

## Typical Use Case (Reynolds Equation)

- Solve the Reynolds equation on a lubricating surface (`from`),
  obtaining a pressure field `p`.
- Project the pressure field onto a surrounding 3D mesh (`to`)
  for use in structural, thermal, or multiphysics coupling.

---
"""
function projectScalarField(pp::ScalarField; from="", to="", gap=false, binSize=0.7)
    problem = pp.model
    if isempty(from) || isempty(to)
        error("projectScalarField: physical groups must be given.")
    end
    if isElementwise(pp)
        p = elementsToNodes(p)
    else
        p = pp
    end
    dimF = 0
    if gap == true
        tag = getTagForPhysicalName(from)
        dimF = getDimForPhysicalName(from)
        nodeTagsUpper, coordUpper = gmsh.model.mesh.getNodesForPhysicalGroup(dimF, tag)
        #h0 = elementsToNodes(ScalarField(problem, from, (x, y, z) -> z))
    end

    # --- 1) 2D elemek összegyűjtése ---
    dimTags2D = gmsh.model.getEntitiesForPhysicalName(from)

    elems2D = Int[]
    conn2D = Dict{Int,Vector{Int}}()
    etype2D = Dict{Int,Int}()
    prin2D = Dict{Int,Int}()

    for (edim, etag) in dimTags2D
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
        for i in eachindex(elemTypes)
            et = Int(elemTypes[i])
            _, dim, _, numNodes, _, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
            dim == 2 || error("projectScalarField: 3D mesh as a source is not yet implemented.")

            tags = elemTags[i]
            nodes = elemNodeTags[i]

            @inbounds for j in eachindex(tags)
                elem = Int(tags[j])
                push!(elems2D, elem)

                i1 = (j - 1) * numNodes + 1
                conn = Vector{Int}(undef, numNodes)
                @inbounds for k in 1:numNodes
                    conn[k] = Int(nodes[i1+k-1])
                end

                conn2D[elem] = conn
                etype2D[elem] = et
                prin2D[elem] = numPrimaryNodes
            end
        end
    end

    # ---- Elem indexelés: elemTag -> idx, és tömbök a gyors belső ciklushoz
    ne = length(elems2D)
    elem2idx = Dict{Int,Int}()
    sizehint!(elem2idx, ne)
    et_of = Vector{Int}(undef, ne)
    pn_of = Vector{Int}(undef, ne)
    conn_of = Vector{Vector{Int}}(undef, ne)

    @inbounds for idx in 1:ne
        e = elems2D[idx]
        elem2idx[e] = idx
        et_of[idx] = etype2D[e]
        pn_of[idx] = prin2D[e]
        conn_of[idx] = conn2D[e]
    end

    # --- 2) Bounding box + minimális elemsize ---
    bbxmin = Inf
    bbxmax = -Inf
    bbymin = Inf
    bbymax = -Inf
    esizemin = Inf

    nodeCoordCache = Dict{Int,NTuple{3,Float64}}()

    dimTags = gmsh.model.getEntitiesForPhysicalName(from)
    for (edim, etag) in dimTags
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(edim, etag)
        bbxmin = min(bbxmin, xmin)
        bbxmax = max(bbxmax, xmax)
        bbymin = min(bbymin, ymin)
        bbymax = max(bbymax, ymax)

        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
        for i in eachindex(elemTypes)
            _, _, _, numNodes, _, _ = gmsh.model.mesh.getElementProperties(elemTypes[i])

            @inbounds for j in eachindex(elemTags[i])
                xemin = Inf
                xemax = -Inf
                yemin = Inf
                yemax = -Inf

                nodes = elemNodeTags[i][(j-1)*numNodes+1:j*numNodes]
                @inbounds for n in nodes
                    coord = get!(nodeCoordCache, Int(n)) do
                        c, _, _, _ = gmsh.model.mesh.getNode(n)
                        (c[1], c[2], c[3])
                    end
                    xemin = min(xemin, coord[1])
                    xemax = max(xemax, coord[1])
                    yemin = min(yemin, coord[2])
                    yemax = max(yemax, coord[2])
                end
                esizemin = min(esizemin, max(xemax - xemin, yemax - yemin))
            end
        end
    end

    bbx = bbxmax - bbxmin
    bby = bbymax - bbymin

    nbbx = Int(bbx ÷ (binSize * esizemin))
    nbby = Int(bby ÷ (binSize * esizemin))

    inv_dx = nbbx / bbx
    inv_dy = nbby / bby

    ebins = [Int[] for _ in 1:(nbbx*nbby)]
    nbins = [Int[] for _ in 1:(nbbx*nbby)]

    # --- 3) Elembinning  (>>> AABB BINNING) ---
    # Elem AABB (x,y) alapján minden lefedett binbe betesszük az elemet.
    # Megjegyzés: clamp + floor ugyanazt a bin-index definíciót használja, mint a node binning.
    for (edim, etag) in dimTags
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
        for i in eachindex(elemTypes)
            _, _, _, numNodes, _, _ = gmsh.model.mesh.getElementProperties(elemTypes[i])

            for j in eachindex(elemTags[i])
                elem_tag = Int(elemTags[i][j])
                idx = elem2idx[elem_tag]

                nodes = elemNodeTags[i][(j-1)*numNodes+1:j*numNodes]

                # AABB a csomópont-koordinátákból (x,y)
                xemin = Inf
                xemax = -Inf
                yemin = Inf
                yemax = -Inf
                @inbounds for n in nodes
                    c = nodeCoordCache[Int(n)]
                    xemin = min(xemin, c[1])
                    xemax = max(xemax, c[1])
                    yemin = min(yemin, c[2])
                    yemax = max(yemax, c[2])
                end

                # lefedett bin tartomány
                nxmin = clamp(Int(floor((xemin - bbxmin) * inv_dx)), 0, nbbx - 1)
                nxmax = clamp(Int(floor((xemax - bbxmin) * inv_dx)), 0, nbbx - 1)
                nymin = clamp(Int(floor((yemin - bbymin) * inv_dy)), 0, nbby - 1)
                nymax = clamp(Int(floor((yemax - bbymin) * inv_dy)), 0, nbby - 1)

                # minden érintett binbe: sorrend determinisztikus (nx outer, ny inner)
                @inbounds for nx in nxmin:nxmax
                    base = nx * nbby
                    @inbounds for ny in nymin:nymax
                        push!(ebins[base+ny+1], idx)
                    end
                end
            end
        end
    end

    # --- 4) Node binning ---
    tag = getTagForPhysicalName(to)
    dimT = getDimForPhysicalName(to)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(dimT, tag)
    coord = reshape(coord, 3, :)

    @inbounds for i in eachindex(nodeTags)
        nx = clamp(Int(floor((coord[1, i] - bbxmin) * inv_dx)), 0, nbbx - 1)
        ny = clamp(Int(floor((coord[2, i] - bbymin) * inv_dy)), 0, nbby - 1)
        push!(nbins[nx*nbby+ny+1], Int(nodeTags[i]))
    end

    # --- 5) Interpoláció ---
    a = zeros(problem.non)
    pvals = vec(p.a)
    uvw = Vector{Float64}(undef, 3)

    if gap == true
        for (ni, nn) in enumerate(nodeTagsUpper)
            gmsh.model.mesh.setNode(nn, [coordUpper[3ni-2], coordUpper[3ni-1], 0.0], [0, 0, 0])
        end
    end

    for b in eachindex(ebins)
        elems = ebins[b]   # idx-ek
        nodes = nbins[b]
        isempty(elems) && continue
        isempty(nodes) && continue

        for node in nodes
            coordN = get!(nodeCoordCache, node) do
                c, _, _, _ = gmsh.model.mesh.getNode(node)
                (c[1], c[2], c[3])
            end

            for idx in elems
                et = et_of[idx]
                pn = pn_of[idx]
                conn = conn_of[idx]

                u, v, w = gmsh.model.mesh.getLocalCoordinatesInElement(
                    elems2D[idx], coordN[1], coordN[2], coordN[3]
                )

                if (pn == 4 && -1.0001 ≤ u ≤ 1.0001 && -1.0001 ≤ v ≤ 1.0001) ||
                   (pn == 3 && -0.01 ≤ u && -0.01 ≤ v && u + v ≤ 1.01)

                    uvw[1] = u
                    uvw[2] = v
                    uvw[3] = w
                    _, bf, _ = gmsh.model.mesh.getBasisFunctions(et, uvw, "Lagrange")

                    val = 0.0
                    @inbounds for i in eachindex(conn)
                        val += bf[i] * pvals[conn[i]]
                    end

                    a[node] = val
                    break
                end
            end
        end
    end

    if gap == true
        for (ni, nn) in enumerate(nodeTagsUpper)
            gmsh.model.mesh.setNode(nn, [coordUpper[3ni-2], coordUpper[3ni-1], coordUpper[3ni]], [0, 0, 0])
        end
    end

    ret = ScalarField([], reshape(a, :, 1), [0.0], [], 1, :scalar, problem)
    if isNodal(pp)
        return ret
    else
        return nodesToElements(ret)
    end
end

function pressureInVolume(p::ScalarField)
    if p.model.geometry.nameVolume == ""
        error("initializePressure: no volume for lubricant has been defined.")
    end

    pp = projectScalarField(p, from=p.model.material[1].phName, to=p.model.geometry.nameVolume)

    if isnothing(p.model.geometry.hh)
        hh = projectScalarField(p.model.geometry.h, from=p.model.material[1].phName, to=p.model.geometry.nameVolume)
        hh = nodesToElements(hh, onPhysicalGroup=p.model.geometry.nameVolume)
        p.model.geometry.hh = hh
    end

    return nodesToElements(pp, onPhysicalGroup=p.model.geometry.nameVolume)
end

#function systemMatrix(problem, αInNodes::ScalarField, velocity::Number, height::ScalarField)
function systemMatrix_old(problem, velocity::Number)
    gmsh.model.setCurrent(problem.name)
    if problem.type ≠ :Reynolds && problem.material[1].type ≠ :JFO
        error("systemMatrix: only Reynolds equation with JFO cavitation model is implemented.")
    end
    if length(problem.material) ≠ 1
        error("systemMatrix: only one material is permitted.")
    end
    if problem.geometry.h == nothing
        #error("systemMatrix: Problem must be initialized with \"initialize(problem)\"")
        initialize(problem)
    end

    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
    nn = []
    I = []
    J = []
    V = []
    V2 = []
    V = convert(Vector{Float64}, V)
    V2 = convert(Vector{Float64}, V2)
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)
    sizehint!(V2, lengthOfIJV)
    ncoord2 = zeros(3 * problem.non)
    phName = problem.material[1].phName
    nameGap = problem.geometry.nameGap
    η = problem.material[1].η
    κ = problem.material[1].κ

    for ipg in 1:1#length(problem.material)
        #η = problem.material.η
        dim = problem.dim
        pdim = problem.pdim
        if dim == 2
            rowsOfB = 2
        elseif dim == 1
            rowsOfB = 1
        else
            error("systemMatrix: rowsOfB = $rowsOfB")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            tag = getTagForPhysicalName(phName)
            nodeTags, ncoord = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                ∂h = zeros(dim, numNodes * numIntPoints)
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                K2 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:pdim
                        H[k*pdim-(pdim-kk), l*pdim-(pdim-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                x = zeros(numIntPoints)
                y = dim == 2 ? zeros(numIntPoints) : 0
                #α = zeros(numIntPoints)
                hh = zeros(numIntPoints)
                dhhdx = zeros(numIntPoints)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        x[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                        if dim == 2
                            y[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 1]
                        end
                        #α[k] = h[:, k]' * αInNodes.a[nnet[j, :]]
                        hh[k] = h[:, k]' * problem.geometry.h.a[nnet[j, :]]
                        dhhdx[k] = h[:, k]' * problem.geometry.dhdx.a[nnet[j, :]]
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 2
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 1 && rowsOfB == 1
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    K2 .*= 0
                    #hh = 0
                    #dhdx = 0
                    for k in 1:numIntPoints
                        #g = α[k] < 1.0 ? 0.0 : 1.0
                        #if dim == 1
                        #    hh = gmsh.view.probe(problem.geometry.tagTop, x[k], 0, 0, -1, -1, false, -1)[1][1]
                        #    hhB = gmsh.view.probe(problem.geometry.tagBottom, x[k], 0, 0, -1, -1, false, -1)[1][1]
                        #    hh -= hhB
                        #    dhdx = gmsh.view.probe(problem.geometry.tagTop, x[k], 0, 0,-1,-1,true, -1)[1][1]
                        #elseif dim == 2
                        #    hh = gmsh.view.probe(problem.geometry.tagTop, x[k], y[k], 0, -1, -1, false, -1)[1][1]
                        #    hhB = gmsh.view.probe(problem.geometry.tagBottom, x[k], y[k], 0, -1, -1, false, -1)[1][1]
                        #    hh -= hhB
                        #    dhdx = gmsh.view.probe(problem.geometry.tagTop, x[k], y[k], 0,-1,-1,true, -1)[1][1]
                        #else
                        #    error("systemMatrixReynoldsCavitation: dim = $dim.")
                        #end
                        H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        dHdx1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB-(rowsOfB-1), 1:pdim*numNodes]
                        K1 += (B1' * B1) * hh[k]^3 * 1 * κ * jacDet[k] * intWeights[k]
                        K2 -= (H1' * H1) * 6.0 * dhhdx[k] * velocity * η / 2 * jacDet[k] * intWeights[k]
                        K2 += (H1' * dHdx1) * 6.0 * hh[k] * velocity * η / 2 * jacDet[k] * intWeights[k]
                        #K2 += jacDet[k] * intWeights[k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                    append!(V2, K2[:])
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    KK = sparse(I, J, V2, dof, dof)
    dropzeros!(K)
    dropzeros!(KK)
    return SystemMatrix(K, problem), SystemMatrix(KK, problem)
    GC.gc()
end

function systemMatrix(problem, velocity::Number)
    gmsh.model.setCurrent(problem.name)

    if problem.type ≠ :Reynolds && problem.material[1].type ≠ :JFO
        error("systemMatrix: only Reynolds equation with JFO cavitation model is implemented.")
    end
    if length(problem.material) ≠ 1
        error("systemMatrix: only one material is permitted.")
    end
    if problem.geometry.h === nothing
        initialize(problem)
    end

    # --- becslés az I,J,V,V2 méretére ---
    elemTypes_all, elemTags_all, elemNodeTags_all = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum((div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.dim)^2 * length(elemTags_all[i]) for i in 1:length(elemTags_all))

    # --- index és érték vektorok, típusosan ---
    I  = Int[]
    J  = Int[]
    V  = Float64[]
    V2 = Float64[]
    sizehint!(I,  lengthOfIJV)
    sizehint!(J,  lengthOfIJV)
    sizehint!(V,  lengthOfIJV)
    sizehint!(V2, lengthOfIJV)

    # csak akkor kell, ha tényleg használod valahol:
    nn = Vector{Matrix{Int}}()

    ncoord2 = zeros(3 * problem.non)
    phName  = problem.material[1].phName
    η       = problem.material[1].η
    κ       = problem.material[1].κ

    dim  = problem.dim
    pdim = problem.pdim
    rowsOfB = dim  # 1D → 1, 2D → 2; az eredeti logika is ez volt

    # (most csak egy anyagot engedünk továbbra is)
    for ipg in 1:1
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        tag     = getTagForPhysicalName(phName)

        # node koordináták egyszer / physical group-onként
        nodeTags, ncoord = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
        @inbounds begin
            ncoord2[nodeTags .* 3 .- 2] .= ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags .* 3 .- 1] .= ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags .* 3 .- 0] .= ncoord[3:3:length(ncoord)]
        end

        for idm in 1:length(dimTags)
            edim, etag = dimTags[idm]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

            for i in 1:length(elemTypes)
                et = elemTypes[i]

                elementName, dim_el, order, numNodes::Int, localNodeCoord, numPrimaryNodes =
                    gmsh.model.mesh.getElementProperties(et)

                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)

                # GradLagrange
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)

                # Lagrange
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)

                # --- preallocáció elem-típus szinten ---
                numElem_i = length(elemTags[i])
                nnet  = Matrix{Int}(undef, numElem_i, numNodes)
                invJac = Matrix{Float64}(undef, 3, 3 * numIntPoints)

                nloc  = numNodes * pdim
                Iidx  = Vector{Int}(undef, nloc * nloc)
                Jidx  = Vector{Int}(undef, nloc * nloc)

                # Iidx, Jidx mátrix-séma (egyszer per elem-típus)
                idx = 1
                @inbounds for k in 1:nloc, l in 1:nloc
                    Iidx[idx] = l
                    Jidx[idx] = k
                    idx += 1
                end

                ∂h  = Matrix{Float64}(undef, dim, numNodes * numIntPoints)
                H   = Matrix{Float64}(undef, pdim * numIntPoints, pdim * numNodes)
                B   = Matrix{Float64}(undef, rowsOfB * numIntPoints, pdim * numNodes)
                K1  = Matrix{Float64}(undef, nloc, nloc)
                K2  = Matrix{Float64}(undef, nloc, nloc)
                nn2 = Vector{Int}(undef, nloc)

                x     = Vector{Float64}(undef, numIntPoints)
                y     = Vector{Float64}(undef, numIntPoints)  # dim==1 esetén nem használjuk
                hh    = Vector{Float64}(undef, numIntPoints)
                dhhdx = Vector{Float64}(undef, numIntPoints)

                # H mátrix: blokk-diagonális ismétlés (eredeti logika tisztábban)
                @inbounds for k in 1:numIntPoints, l in 1:numNodes
                    val = h[(k-1)*numNodes + l]
                    for kk in 1:pdim
                        row = (k-1)*pdim + kk
                        col = (l-1)*pdim + kk
                        H[row, col] = val
                    end
                end

                # --- elemenkénti ciklus ugyanazon elem-típusra ---
                for j in 1:numElem_i
                    elem = elemTags[i][j]

                    @inbounds for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes + k]
                    end

                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)

                    # inverz Jacobian és x,y,hh,dhhdx
                    @inbounds for k in 1:numIntPoints
                        # 3×3 blokk
                        jstart = 3k - 2
                        @views Jblk = Jac[1:3, jstart:jstart+2]
                        @views invJblk = inv(Jblk)'   # ha kell, később lehet kézi 3×3 inverzre cserélni
                        invJac[1:3, jstart:jstart+2] .= invJblk

                        # x, y, hh, dhhdx interpoláció
                        # shape: h[k] ∈ R^{numNodes}
                        # nnet[j, :] node indexek
                        xk = 0.0
                        yk = 0.0
                        hhk = 0.0
                        dhdxk = 0.0
                        for l in 1:numNodes
                            w = h[l, k]
                            nd = nnet[j, l]
                            idx3 = 3nd
                            xk    += w * ncoord2[idx3 - 2]
                            if dim == 2
                                yk += w * ncoord2[idx3 - 1]
                            end
                            hhk   += w * problem.geometry.h.a[nd]
                            dhdxk += w * problem.geometry.dhdx.a[nd]
                        end
                        x[k]     = xk
                        if dim == 2
                            y[k] = yk
                        end
                        hh[k]    = hhk
                        dhhdx[k] = dhdxk
                    end

                    # ∂h kinullázása és kiszámítása
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numIntPoints, l in 1:numNodes
                        # grad shape in fizikai tér: invJac * ∇h
                        # ∇h indexeléssel óvatosan: [l*3-2 : l*3-(3-dim)] × k
                        col_grad_start = 3l - 2
                        col_grad_end   = 3l - (3-dim)
                        jstart = 3k - 2
                        jend   = 3k - (3-dim)
                        @views g_ref = ∇h[col_grad_start:col_grad_end, k]
                        @views invJd = invJac[1:dim, jstart:jend]
                        ∂h[:, (k-1)*numNodes + l] .= invJd * g_ref
                    end

                    # B kinullázása és kitöltése
                    fill!(B, 0.0)
                    if dim == 2 && rowsOfB == 2
                        @inbounds for k in 1:numIntPoints, l in 1:numNodes
                            row1 = (k-1)*rowsOfB + 1
                            row2 = row1 + 1
                            col  = l*pdim
                            B[row1, col] = ∂h[1, (k-1)*numNodes + l]
                            B[row2, col] = ∂h[2, (k-1)*numNodes + l]
                        end
                    elseif dim == 1 && rowsOfB == 1
                        @inbounds for k in 1:numIntPoints, l in 1:numNodes
                            row = k
                            col = l*pdim
                            B[row, col] = ∂h[1, (k-1)*numNodes + l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end

                    # K1, K2 kinullázása
                    fill!(K1, 0.0)
                    fill!(K2, 0.0)

                    # integráció pontonként
                    @inbounds for k in 1:numIntPoints
                        jac_k = jacDet[k]
                        w_k   = intWeights[k]

                        # lokális blokkok
                        @views H1    = H[(k-1)*pdim+1 : k*pdim, :]
                        @views B1    = B[(k-1)*rowsOfB+1 : k*rowsOfB, :]
                        @views dHdx1 = B[(k-1)*rowsOfB+1 : (k-1)*rowsOfB+1, :]  # eredetiben ez volt

                        # skálák
                        fac1 = hh[k]^3 * κ * jac_k * w_k
                        fac2 = 6.0 * velocity * η / 2 * jac_k * w_k

                        # K1 += (B1' * B1) * fac1
                        K1 .+= (B1' * B1) .* fac1

                        # K2 -= (H1' * H1) * fac2 * dhhdx[k]
                        # K2 += (H1' * dHdx1) * fac2 * hh[k]
                        K2 .-= (H1' * H1)    .* (fac2 * dhhdx[k])
                        K2 .+= (H1' * dHdx1) .* (fac2 * hh[k])
                    end

                    # globális szabadsági fok indexek (nn2)
                    @inbounds for k in 1:pdim
                        # k-adik komponens: k, k+pdim, ...
                        start = k
                        for (idx_node, nd) in enumerate(nnet[j, :])
                            nn2[start + (idx_node-1)*pdim] = pdim * nd - (pdim - k)
                        end
                    end

                    # I, J, V, V2 bővítése
                    append!(I, nn2[Iidx])
                    append!(J, nn2[Jidx])
                    append!(V,  K1[:])
                    append!(V2, K2[:])
                end

                push!(nn, nnet)
            end
        end
    end

    dof = problem.pdim * problem.non
    K   = sparse(I, J, V,  dof, dof)
    KK  = sparse(I, J, V2, dof, dof)
    dropzeros!(K)
    dropzeros!(KK)

    return SystemMatrix(K, problem), SystemMatrix(KK, problem)
end

function flowRateVector(problem, loads)
    gmsh.model.setCurrent(problem.name)
    dim0 = problem.dim
    pdim = problem.pdim
    non = problem.non
    dof = pdim * non
    #η = problem.material[1].η
    fp = zeros(dof)
    ncoord2 = zeros(3 * problem.non)
    for n in 1:length(loads)
        name, ff = loads[n]
        f = ff
        dimTags = gmsh.model.getEntitiesForPhysicalName(name)
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
                nnoe = reshape(elemNodeTags[ii], numNodes, :)'
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(order+1))
                numIntPoints = length(intWeights)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elementTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                        end
                    end
                end
                ∂h = zeros(dim0, numNodes * numIntPoints)
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[i][(l-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim0, (k-1)*numNodes+l] .= invJac[1:dim0, k*3-2:k*3-(3-dim0)] * ∇h[l*3-2:l*3-(3-dim0), k]
                    end
                    f1 .*= 0
                    for j in 1:numIntPoints
                        x = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 2]
                        if isa(ff, Function)
                            y = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 1]
                            f = ff(x, y)
                        end
                        H1 = H[j*pdim-(pdim-1):j*pdim, 1:pdim*numNodes] # H1[...] .= H[...] ????
                        ############### NANSON ###########################################
                        if pdim == 1 && dim == 2
                            Ja = jacDet[j]
                        elseif pdim == 1 && dim == 1
                            Ja = jacDet[j]
                            #Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2)
                            #Ja = 1
                        elseif pdim == 1 && dim == 0
                            Ja = 1
                        else
                            error("massFlowVectorReynolds: dimension of the problem is $(problem.pdim), dimension of load is $dim.")
                        end
                        f1 += H1' * f * Ja * intWeights[j]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                    end
                    fp[nn2] += f1
                end
            end
        end
    end
    type = :none
    if problem.dim == 1
        type = :scalar
    elseif problem.dim == 2
        type = :scalar
    else
        error("problem.dim = $(problem.dim)")
    end
    return ScalarField([], reshape(fp, :,1), [0.0], [], 1, type, problem)
end

#=
function solvePressure(problem, load, BC, V; cav=false, periodicSlave="", periodicMaster="")
    if problem.type == :dummy
        return nothing
    end
    if problem.type != :Reynolds
        error("solvePressure: bad problem type '$type'.")
    end
    p0 = problem.material[1].p₀
    κ = problem.material[1].κ
    if problem.geometry.h == nothing
        initialize(problem)
    end
    fluid = constrainedDoFs(problem, [pressureConstraint(problem.material[1].phName, p=0)])
    bc = constrainedDoFs(problem, BC)
    free = setdiff(fluid, bc)
    #one = zeros(problem.non * problem.pdim) .+ 1
    #type = :none
    #if problem.dim == 1
    #    type = :scalar
    #elseif problem.dim == 2
    #    type = :scalar
    #else
    #    error("problem.dim = $(problem.dim)")
    #end
    #α = ScalarField([], reshape(one, :,1), [0.0], [], 1, type, problem)
    α = scalarField(problem, problem.material[1].phName, 1)
    α0 = copy(α)
    f0 = flowRateVector(problem, load)
    f = copy(f0)
    err = 1
    c = 0
    if problem.dim == 1
        if length(periodicSlave) != 0
            pbc1 = getTagForPhysicalName(periodicSlave)
            nbc1::Int64 = gmsh.model.mesh.getNodesForPhysicalGroup(0, pbc1)[1][1]
            pbc2 = getTagForPhysicalName(periodicMaster)
            nbc2::Int64 = gmsh.model.mesh.getNodesForPhysicalGroup(0, pbc2)[1][1]
        
            G = zeros(1, problem.non)
            #GT = zeros(problem.non)

            G[nbc1] = 1
            G[nbc2] = -1
            GT = G'

            while err > 1e-10 && c < 101
                c += 1
                α0 = copy(α)
                fill!(α.a, 0.0)
                K, KK = systemMatrix(problem, -V)
                f = copy(f0)
                f -= KK * α0
                applyBoundaryConditions!(α, BC)
                α.a[bc] .= exp.((α.a[bc] .- p0) ./ κ)
                p1 = K.A[:, bc] * α.a[bc]
                #Reynolds.applyBoundaryConditions!(problem, K, f, BC)
                #a0 = similar(α0.a)
                a0 = K.A[free, free] \ (f.a[free] - p1[free])
                KG = K \ GT
                λ = (G * KG) \ G * a0
                #a1 = ScalarField([], -KG * λ, [0.0], [], 1, a0.type, a0.model)
                a1 = -KG * λ
                α.a[free] = a0 + a1
                
                err = sum(abs, α.a - α0.a) / sum(abs, α0.a)
                if c == 100
                    @warn("solvePressure: number of iterations reached step $c.")
                end
            end
        elseif length(periodicSlave) == 0
            while err > 1e-10 && c < 101
                c += 1
                α0 = copy(α)
                K, KK = systemMatrix(problem, -V)
                f = copy(f0)
                f -= KK * α
                #f .= f0
                applyBoundaryConditions!(problem, K, f, BC)
                α = K \ f
                err = sum(abs, α.a - α0.a) / sum(abs, α0.a)
                if c == 100
                    @warn("solvePressure: number of iterations reached step $c.")
                end
            end
        end
    elseif problem.dim == 2
        if length(periodicSlave) != 0
            pbc1 = getTagForPhysicalName(periodicSlave) # right0: slave
            cv1 = (gmsh.model.getEntitiesForPhysicalGroup(1, pbc1))[1]
            tagMaster, slave::Vector{Int64}, master::Vector{Int64}, affineTransform = gmsh.model.mesh.getPeriodicNodes(1, cv1, true)
            @info "master = $tagMaster, slave = $cv1"

            G = zeros(length(master), problem.non)
            GT = zeros(problem.non, length(master))

            for i in 1:length(master)
                G[i, master[i]] = 1
                G[i, slave[i]] = -1

                GT[master[i], i] = 1
                GT[slave[i], i] = -1
            end
            
            while err > 1e-10 && c < 101
                c += 1
                α0 = copy(α)
                K, KK = systemMatrix(problem, -V)
                f = copy(f0)
                f -= KK * α0
                applyBoundaryConditions!(α, BC)
                α.a[bc] .= exp.((α.a[bc] .- p0) ./ κ)
                p1 = K.A[:, bc] * α.a[bc]
                a0 = K.A[free, free] \ (f.a[free] - p1[free]) # free × 1
                KG = K.A[free, free] \ GT[free, :] # free × master
                λ = (G[:, free] * KG) \ G[:, free] * a0 # master × 1
                a1 = -KG * λ # free × 1
                α.a[free] = a0 + a1
                
                err = sum(abs, α.a - α0.a) / sum(abs, α0.a)
                if c == 100
                    @warn("solvePressure: number of iterations reached step $c.")
                end
            end

        elseif length(periodicSlave) == 0

            while err > 1e-10 && c < 101
                c += 1
                α0 = copy(α)
                K, KK = systemMatrix(problem, -V)
                f = copy(f0)
                f -= KK * α
                applyBoundaryConditions!(α, BC)
                α.a[bc] .= exp.((α.a[bc] .- p0) ./ κ)
                p1 = K.A[:, bc] * α.a[bc]
                a0 = K.A[free, free] \ (f.a[free,1] - p1[free])
                α.a[free] = a0
                
                err = sum(abs, α.a - α0.a) / sum(abs, α0.a)
                if c == 100
                    @warn("solvePressure: number of iterations reached step $c.")
                end
            end
        end
    end
    @info "solvePressure: number of iterations is $c."
    pret = zeros(problem.non, 1)
    αcav = zeros(problem.non, 1)
    g = [i < 1.0 ? 0.0 : 1.0 for i in α.a[fluid]]
    αcav[fluid] = 1 .- [i < 1.0 ? i : 1.0 for i in α.a[fluid]]
    α = [i < 1 ? 1 : i for i in α.a[fluid]]
    if cav == false
        pret[fluid] =  p0 .+ κ * g .* log.(α)
        return ScalarField([], reshape(pret, :,1), [0], [], 1, :other, problem)
    else
        #return α, αcav # log.(α), αcav
        pret[fluid] =  p0 .+ κ * g .* log.(α)
        return ScalarField([], reshape(pret, :,1), [0], [], 1, :other, problem), ScalarField([], reshape(αcav, :,1), [0], [], 1, :p, problem)
    end
end
=#

function solvePressure(problem, load, BC, V; cav=false, periodicSlave="", periodicMaster="")
    if problem.type == :dummy
        return nothing
    end
    if problem.type != :Reynolds
        error("solvePressure: bad problem type '$type'.")
    end
    p0 = problem.material[1].p₀
    κ = problem.material[1].κ
    if problem.geometry.h == nothing
        initialize(problem)
    end
    fluid = constrainedDoFs(problem, [pressureConstraint(problem.material[1].phName, p=0)])
    bc = constrainedDoFs(problem, BC)
    free = setdiff(fluid, bc)
    #one = zeros(problem.non * problem.pdim) .+ 1
    #type = :none
    #if problem.dim == 1
    #    type = :scalar
    #elseif problem.dim == 2
    #    type = :scalar
    #else
    #    error("problem.dim = $(problem.dim)")
    #end
    #α = ScalarField([], reshape(one, :,1), [0.0], [], 1, type, problem)
    α = scalarField(problem, problem.material[1].phName, 1)
    α0 = copy(α)
    f0 = LowLevelFEM.flowRateVector(problem, load)
    f = copy(f0)
    err = 1
    c = 0
    if problem.dim == 1
        if length(periodicSlave) != 0
            pbc1 = LowLevelFEM.getTagForPhysicalName(periodicSlave)
            nbc1::Int64 = gmsh.model.mesh.getNodesForPhysicalGroup(0, pbc1)[1][1]
            pbc2 = LowLevelFEM.getTagForPhysicalName(periodicMaster)
            nbc2::Int64 = gmsh.model.mesh.getNodesForPhysicalGroup(0, pbc2)[1][1]

            G = zeros(1, problem.non)
            #GT = zeros(problem.non)

            G[nbc1] = 1
            G[nbc2] = -1
            GT = G'

            while err > 1e-10 && c < 101
                c += 1
                α0 = copy(α)
                fill!(α.a, 0.0)
                K, KK = LowLevelFEM.systemMatrix(problem, -V)
                f = copy(f0)
                f -= KK * α0
                applyBoundaryConditions!(α, BC)
                α.a[bc] .= exp.((α.a[bc] .- p0) ./ κ)
                p1 = K.A[:, bc] * α.a[bc]
                #Reynolds.applyBoundaryConditions!(problem, K, f, BC)
                #a0 = similar(α0.a)
                a0 = K.A[free, free] \ (f.a[free] - p1[free])
                KG = K \ GT
                λ = (G * KG) \ G * a0
                #a1 = ScalarField([], -KG * λ, [0.0], [], 1, a0.type, a0.model)
                a1 = -KG * λ
                α.a[free] = a0 + a1

                err = sum(abs, α.a - α0.a) / sum(abs, α0.a)
                if c == 100
                    @warn("solvePressure: number of iterations reached step $c.")
                end
            end
        elseif length(periodicSlave) == 0
            while err > 1e-10 && c < 101
                c += 1
                α0 = copy(α)
                K, KK = LowLevelFEM.systemMatrix(problem, -V)
                f = copy(f0)
                f -= KK * α
                #f .= f0
                applyBoundaryConditions!(problem, K, f, BC)
                α = K \ f
                err = sum(abs, α.a - α0.a) / sum(abs, α0.a)
                if c == 100
                    @warn("solvePressure: number of iterations reached step $c.")
                end
            end
        end
    elseif problem.dim == 2
        if length(periodicSlave) != 0
            pbc1 = LowLevelFEM.getTagForPhysicalName(periodicSlave) # right0: slave
            cv1 = (gmsh.model.getEntitiesForPhysicalGroup(1, pbc1))[1]
            tagMaster, slave::Vector{Int64}, master::Vector{Int64}, affineTransform = gmsh.model.mesh.getPeriodicNodes(1, cv1, true)
            @info "master = $tagMaster, slave = $cv1"

            G = zeros(length(master), problem.non)
            #GT = zeros(problem.non, length(master))

            for i in 1:length(master)
                G[i, master[i]] = 1
                G[i, slave[i]] = -1

                #GT[master[i], i] = 1
                #GT[slave[i], i] = -1
            end

            indG = findall(x -> x in free, master)
            #indG2 = findall(x -> x in free, slave)

            while err > 1e-10 && c < 101
                c += 1
                α0 = copy(α)
                K, KK = LowLevelFEM.systemMatrix(problem, -V)
                f = copy(f0)
                f -= KK * α0
                applyBoundaryConditions!(α, BC)
                α.a[bc] .= exp.((α.a[bc] .- p0) ./ κ)
                p1 = K.A[:, bc] * α.a[bc]
                a0 = K.A[free, free] \ (f.a[free] - p1[free]) # free × 1
                KG = K.A[free, free] \ G'[free, indG] # free × master
                #@disp fluid
                #@disp bc
                #@disp free
                #@disp master
                #@disp indG
                #@disp indG2
                #display(G[indG, master])
                #display(K.A[free, free])
                #display(KG)
                #@disp GT
                λ = (G[indG, free] * KG) \ G[indG, free] * a0 # master × 1
                a1 = -KG * λ # free × 1
                α.a[free] = a0 + a1

                err = sum(abs, α.a - α0.a) / sum(abs, α0.a)
                if c == 100
                    @warn("solvePressure: number of iterations reached step $c.")
                end
            end

        elseif length(periodicSlave) == 0

            while err > 1e-10 && c < 101
                c += 1
                α0 = copy(α)
                K, KK = LowLevelFEM.systemMatrix(problem, -V)
                f = copy(f0)
                f -= KK * α
                applyBoundaryConditions!(α, BC)
                α.a[bc] .= exp.((α.a[bc] .- p0) ./ κ)
                p1 = K.A[:, bc] * α.a[bc]
                a0 = K.A[free, free] \ (f.a[free, 1] - p1[free])
                α.a[free] = a0

                err = sum(abs, α.a - α0.a) / sum(abs, α0.a)
                if c == 100
                    @warn("solvePressure: number of iterations reached step $c.")
                end
            end
        end
    end
    @info "solvePressure: number of iterations is $c."
    pret = zeros(problem.non, 1)
    αcav = zeros(problem.non, 1)
    g = [i < 1.0 ? 0.0 : 1.0 for i in α.a[fluid]]
    αcav[fluid] = 1 .- [i < 1.0 ? i : 1.0 for i in α.a[fluid]]
    α = [i < 1 ? 1 : i for i in α.a[fluid]]
    if cav == false
        pret[fluid] = p0 .+ κ * g .* log.(α)
        return ScalarField([], reshape(pret, :, 1), [0], [], 1, :other, problem)
    else
        #return α, αcav # log.(α), αcav
        pret[fluid] = p0 .+ κ * g .* log.(α)
        return ScalarField([], reshape(pret, :, 1), [0], [], 1, :other, problem), ScalarField([], reshape(αcav, :, 1), [0], [], 1, :p, problem)
    end
end

function solveShearStress(p, V; component=:xz)
    grad_p = grad_xy(p)
    grad_p = expandTo3D(grad_p)

    ex = VectorField(p.model, p.model.material[1].phName, [1, 0, 0])
    ey = VectorField(p.model, p.model.material[1].phName, [0, 1, 0])

    px = grad_p ⋅ ex
    py = grad_p ⋅ ey

    h = p.model.geometry.h;

    η = p.model.material[1].η

    if component == :xz
        τxz = -px * h / 2 - η / h * V
        return τxz
    elseif component == :yz
        τyz = -py * h / 2
        return τyz
    else
        error("solveShearStress: wrong component name ($(component)), use :xz or :yz.")
    end
end
