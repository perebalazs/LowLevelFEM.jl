export projectScalarField, fieldsToVolume


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

#function pressureConstraint(name::String; p=1im)
#    bc0 = name, p, 1im, 1im
#    return bc0
#end

#function flowRate(name::String; q=0)
#    ld0 = name, q
#    return ld0
#end

function initialize(problem::Problem)
    h0 = ScalarField(problem, problem.geometry.nameGap, (x, y, z) -> z)
    h0 = elementsToNodes(h0)
    h = projectScalarField(h0, from=problem.geometry.nameGap, to=problem.material[1].phName, gap=true)
    h0 = nodesToElements(h, onPhysicalGroup=problem.material[1].phName)
    dhdx = elementsToNodes(∂x(h0))
    #dhdx = elementsToNodes(dhdx)
    problem.geometry.h = h
    problem.geometry.dhdx = dhdx
end

"""
    projectScalarField(p::Union{ScalarField,Vector{ScalarField}}; from="", to="", gap=false, binSize=0.7)

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

- `ScalarField` or `Vector{ScalarField}`
  A new scalar field(s) defined on the target mesh nodes, containing the
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
function projectScalarField(pp::Union{ScalarField,Vector{ScalarField}}; from="", to="", gap=false, binSize=0.7)
    problem = pp isa Vector ? pp[1].model : pp.model
    if isempty(from) || isempty(to)
        error("projectScalarField: physical groups must be given.")
    end
    rows = pp isa Vector ? length(pp) : 1
    p = pp isa Vector ? [elementsToNodes(pp[i]) for i in eachindex(pp)] : [elementsToNodes(pp)]
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

    nbbx = max(1, Int(bbx ÷ (binSize * esizemin)))
    nbby = max(1, Int(bby ÷ (binSize * esizemin)))

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
    a = zeros(problem.non, rows)
    pvals = Matrix{Float64}(undef, problem.non, rows)
    for j in 1:rows
        pvals[:, j] .= vec(p[j].a)
    end

    #pvals = [vec(p[i].a) for i in eachindex(p)]
    uvw = Vector{Float64}(undef, 3)

    if gap == true
        for (ni, nn) in enumerate(nodeTagsUpper)
            gmsh.model.mesh.setNode(nn, [coordUpper[3ni-2], coordUpper[3ni-1], 0.0], [0, 0, 0])
        end
    end

    val = Vector{Float64}(undef, rows)
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

                    fill!(val, 0.0)
                    @inbounds for ii in eachindex(conn)
                        nid = conn[ii]
                        bfi = bf[ii]
                        for j in 1:rows
                            val[j] += bfi * pvals[nid, j]
                        end
                    end

                    #fill!(val, 0.0)
                    #@inbounds for i in eachindex(conn), j in 1:rows
                    #    val[j] += bf[i] * pvals[conn[i], j]
                    #end

                    @inbounds for i in 1:rows
                        a[node, i] = val[i]
                    end
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

    if pp isa Vector
        ret = [ScalarField([], reshape(a[:,i], :, 1), [0.0], [], 1, :scalar, problem) for i in 1:rows]
        ret = [isNodal(pp[i]) ? ret[i] : nodesToElements(ret[i], onPhysicalGroup=to) for i in 1:rows]
    else
        ret = ScalarField([], reshape(a[:, 1], :, 1), [0.0], [], 1, :scalar, problem)
        ret = isNodal(pp) ? ret : nodesToElements(ret, onPhysicalGroup=to)
    end
    return ret
end

function fieldsToVolume(p0::Union{ScalarField,Vector{ScalarField}})
    if p0 isa Vector
        p = p0
    else
        p = [p0]
    end
    if p[1].model.geometry.nameVolume == ""
        error("initializePressure: no volume for lubricant has been defined.")
    end

    if isnothing(p[1].model.geometry.hh)
        ret = projectScalarField(vcat(p[1].model.geometry.h, p), from=p[1].model.material[1].phName, to=p[1].model.geometry.nameVolume)
        hh = nodesToElements(ret[1], onPhysicalGroup=p[1].model.geometry.nameVolume)
        p[1].model.geometry.hh = hh
        pp = ret[2:end]
        pp = [nodesToElements(pp[i], onPhysicalGroup=p[1].model.geometry.nameVolume) for i in eachindex(pp)]
        return pp
    else
        pp = projectScalarField(p, from=p[1].model.material[1].phName, to=p[1].model.geometry.nameVolume)
        pp = [nodesToElements(pp[i], onPhysicalGroup=p[1].model.geometry.nameVolume) for i in eachindex(pp)]
        return pp
    end
end
