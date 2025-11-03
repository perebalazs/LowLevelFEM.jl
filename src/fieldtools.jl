module FieldTools

using LinearAlgebra: dot
using ..LowLevelFEM

export probe_field_bulk

"""
    build_bbox(rr::ScalarField)

Előkészíti az elemek *bounding boxait* gyors térbeli szűréshez.
Minden elemhez `(xmin, xmax, ymin, ymax, zmin, zmax)` értékeket tárol,
amelyeket a `probe_field_bulk` használ majd a keresés gyorsítására.
"""
function build_bbox(rr::ScalarField)
    gmsh.model.setCurrent(rr.model.name)

    elem_bboxes = Vector{NTuple{6,Float64}}(undef, length(rr.numElem))
    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    X, Y, Z = coords[1:3:end], coords[2:3:end], coords[3:3:end]
    _, _, dim, _ = gmsh.model.mesh.getElement(rr.numElem[1])
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)

    counter = 1
    for (i, et) in enumerate(elemTypes)
        _, dim, order, numNodes, _, _ = gmsh.model.mesh.getElementProperties(et)
        for j in 1:length(elemTags[i])
            nodes = elemNodeTags[i][(j-1)*numNodes+1:j*numNodes]
            xs, ys, zs = X[nodes], Y[nodes], Z[nodes]
            elem_bboxes[counter] = (
                minimum(xs), maximum(xs),
                minimum(ys), maximum(ys),
                minimum(zs), maximum(zs)
            )
            counter += 1
        end
    end

    return elem_bboxes
end



"""
    probe_field_bulk(rr::ScalarField, points::Matrix{Float64})

Gyors, tisztán Julia-alapú tömeges lekérdezés (`(x,y,z)` pontok → mezőértékek).
A `points` mátrix mérete `N×3` legyen (egy sor egy pontot jelent).
Előfeltétel: a `rr.bbox` mező a `build_bbox(rr)` eredménye legyen.
"""
function probe_field_bulk(rr::ScalarField, points::Matrix{Float64})
    gmsh.model.setCurrent(rr.model.name)
    npts = size(points, 1)
    vals = fill(NaN, npts)

    rr_bbox = build_bbox(rr)

    # 1. Lekérjük a hálóhoz tartozó adatokat egyszer
    _, _, dim, _ = gmsh.model.mesh.getElement(rr.numElem[1])
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    xnode, ynode, znode = coords[1:3:end], coords[2:3:end], coords[3:3:end]

    # 2. Végigmegyünk az összes kiértékelendő ponton
    @inbounds for ipt in 1:npts
        x, y, z = points[ipt, 1], points[ipt, 2], points[ipt, 3]

        # 2a. Gyors bbox szűrés
        for (i, bbox) in enumerate(rr_bbox)
            (xmin, xmax, ymin, ymax, zmin, zmax) = bbox
            if !(x ≥ xmin && x ≤ xmax && y ≥ ymin && y ≤ ymax && z ≥ zmin && z ≤ zmax)
                continue
            end

            # 2b. Megvizsgáljuk az elemen belüli interpolációt
            _, _, dim, _ = gmsh.model.mesh.getElement(rr.numElem[i])
            et = elemTypes[1]  # feltételezve azonos típus a mezőben
            _, dim, order, numNodes, _, _ = gmsh.model.mesh.getElementProperties(et)
            nodes = elemNodeTags[1][(i-1)*numNodes+1:i*numNodes]
            xs, ys, zs = xnode[nodes], ynode[nodes], znode[nodes]

            # Természetes koordináták (ξ, η, ζ)
            ξηζ = natural_coordinates_generic(x, y, z, xs, ys, zs, et)
            if ξηζ === nothing
                continue
            end
            ξ, η, ζ = ξηζ

            # Lagrange-alakfüggvények és interpoláció
            N = lagrange_shape_generic(ξ, η, ζ, et)
            if length(N) != numNodes
                continue
            end

            vals[ipt] = dot(N, rr.A[i][:, 1])
            break  # kilépés az első megfelelő elemmel
        end
    end

    return vals
end

"""
    probe_field(rr::ScalarField, x, y, z)

Gyors, tisztán Julia-alapú lekérdezés (Gmsh nélkül), amely visszaadja
a `ScalarField` értékét az `(x, y, z)` pontban.

Támogatott elemtípusok: `tet4`, `tet10`, `hex8`, `hex20`, `tri3`, `tri6`.
"""
function probe_field(rr::ScalarField, x::Float64, y::Float64, z::Float64)
    gmsh.model.setCurrent(rr.model.name)

    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    X = coords[1:3:end]; Y = coords[2:3:end]; Z = coords[3:3:end]

    if isempty(rr.numElem)
        error("probe_field: field has no associated elements.")
    end
    _, _, dim, _ = gmsh.model.mesh.getElement(rr.numElem[1])
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)

    @inbounds for (i, et) in enumerate(elemTypes)
        elementName, dim, order, numNodes::Int, localCoord, numPrimary =
            gmsh.model.mesh.getElementProperties(et)

        for (j, elem) in enumerate(elemTags[i])
            if !(elem in rr.numElem)
                continue
            end

            node_ids = elemNodeTags[i][(j-1)*numNodes+1 : j*numNodes]
            xnodes = X[node_ids]; ynodes = Y[node_ids]; znodes = Z[node_ids]

            # gyors bbox szűrés
            if x < minimum(xnodes)-1e-8 || x > maximum(xnodes)+1e-8 ||
               y < minimum(ynodes)-1e-8 || y > maximum(ynodes)+1e-8 ||
               z < minimum(znodes)-1e-8 || z > maximum(znodes)+1e-8
                continue
            end

            ξηζ = natural_coordinates_generic(x, y, z, xnodes, ynodes, znodes, et)
            if ξηζ === nothing
                continue
            end
            ξ, η, ζ = ξηζ

            N = lagrange_shape_generic(ξ, η, ζ, et)
            if length(N) != numNodes
                continue
            end

            ind = findfirst(isequal(elem), rr.numElem)
            vals = rr.A[ind][:, 1]
            return dot(N, vals)
        end
    end
    return NaN
end


# ---------------------------------------------------------------------
# Shape függvények
# ---------------------------------------------------------------------

function lagrange_shape_generic(ξ, η, ζ, et)
    if et == 4       # tet4
        return [1-ξ-η-ζ, ξ, η, ζ]

    elseif et == 5   # hex8
        return 0.125 .* [
            (1-ξ)*(1-η)*(1-ζ),
            (1+ξ)*(1-η)*(1-ζ),
            (1+ξ)*(1+η)*(1-ζ),
            (1-ξ)*(1+η)*(1-ζ),
            (1-ξ)*(1-η)*(1+ζ),
            (1+ξ)*(1-η)*(1+ζ),
            (1+ξ)*(1+η)*(1+ζ),
            (1-ξ)*(1+η)*(1+ζ)
        ]

    elseif et == 11  # tet10
        L1 = 1-ξ-η-ζ; L2=ξ; L3=η; L4=ζ
        return [
            L1*(2L1-1),
            L2*(2L2-1),
            L3*(2L3-1),
            L4*(2L4-1),
            4*L1*L2,
            4*L2*L3,
            4*L3*L1,
            4*L1*L4,
            4*L2*L4,
            4*L3*L4
        ]

    elseif et == 17  # hex20
        N = Vector{Float64}(undef, 20)
        ξi = [-1,1,1,-1,-1,1,1,-1]
        ηi = [-1,-1,1,1,-1,-1,1,1]
        ζi = [-1,-1,-1,-1,1,1,1,1]
        for n in 1:8
            N[n] = 0.125*(1+ξi[n]*ξ)*(1+ηi[n]*η)*(1+ζi[n]*ζ)*(ξi[n]*ξ+ηi[n]*η+ζi[n]*ζ-2)
        end
        N[9]  = 0.25*(1-ξ^2)*(1-η)*(1-ζ)
        N[10] = 0.25*(1-η^2)*(1+ξ)*(1-ζ)
        N[11] = 0.25*(1-ξ^2)*(1+η)*(1-ζ)
        N[12] = 0.25*(1-η^2)*(1-ξ)*(1-ζ)
        N[13] = 0.25*(1-ξ^2)*(1-η)*(1+ζ)
        N[14] = 0.25*(1-η^2)*(1+ξ)*(1+ζ)
        N[15] = 0.25*(1-ξ^2)*(1+η)*(1+ζ)
        N[16] = 0.25*(1-η^2)*(1-ξ)*(1+ζ)
        N[17] = 0.25*(1-ζ^2)*(1-ξ)*(1-η)
        N[18] = 0.25*(1-ζ^2)*(1+ξ)*(1-η)
        N[19] = 0.25*(1-ζ^2)*(1+ξ)*(1+η)
        N[20] = 0.25*(1-ζ^2)*(1-ξ)*(1+η)
        return N

    elseif et == 2   # tri3
        return [1-ξ-η, ξ, η]

    elseif et == 9   # tri6
        L1=1-ξ-η; L2=ξ; L3=η
        return [L1*(2L1-1), L2*(2L2-1), L3*(2L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1]

    else
        error("lagrange_shape_generic: elemtype $et not supported yet")
    end
end


# ---------------------------------------------------------------------
# Természetes koordináták (ξ, η, ζ)
# ---------------------------------------------------------------------

function natural_coordinates_generic(x, y, z, xnodes, ynodes, znodes, et)
    if et in (4, 11)  # tetra
        A = [
            xnodes[2]-xnodes[1] xnodes[3]-xnodes[1] xnodes[4]-xnodes[1];
            ynodes[2]-ynodes[1] ynodes[3]-ynodes[1] ynodes[4]-ynodes[1];
            znodes[2]-znodes[1] znodes[3]-znodes[1] znodes[4]-znodes[1];
        ]
        b = [x - xnodes[1], y - ynodes[1], z - znodes[1]]
        ξηζ = A \ b
        ξ, η, ζ = ξηζ
        if ξ >= -1e-8 && η >= -1e-8 && ζ >= -1e-8 && ξ + η + ζ <= 1.0 + 1e-8
            return (ξ, η, ζ)
        else
            return nothing
        end

    elseif et in (5, 17)  # hex8, hex20
        ξ = η = ζ = 0.0
        for it in 1:5
            N = lagrange_shape_generic(ξ, η, ζ, et)
            xe = sum(N .* xnodes)
            ye = sum(N .* ynodes)
            ze = sum(N .* znodes)
            dx = x - xe; dy = y - ye; dz = z - ze
            if abs(dx)+abs(dy)+abs(dz) < 1e-10
                break
            end
            ξ += dx / (maximum(xnodes)-minimum(xnodes))
            η += dy / (maximum(ynodes)-minimum(ynodes))
            ζ += dz / (maximum(znodes)-minimum(znodes))
        end
        if abs(ξ)<=1.01 && abs(η)<=1.01 && abs(ζ)<=1.01
            return (ξ, η, ζ)
        else
            return nothing
        end

    elseif et in (2, 9)  # tri3, tri6
        A = [
            xnodes[2]-xnodes[1] xnodes[3]-xnodes[1];
            ynodes[2]-ynodes[1] ynodes[3]-ynodes[1];
        ]
        b = [x - xnodes[1], y - ynodes[1]]
        ξη = A \ b
        ξ, η = ξη
        if ξ >= -1e-8 && η >= -1e-8 && ξ + η <= 1.0 + 1e-8
            return (ξ, η, 0.0)
        else
            return nothing
        end

    else
        return nothing
    end
end

end # module

