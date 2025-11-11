module FieldTools

using LinearAlgebra: dot
using ..LowLevelFEM

export probe_field_bulk, build_surface_grid
export probe_field_bulk_fallback

#using LinearAlgebra

const SUPPORTED_ELEMENT_TYPES = Set([1, 8, 2, 9, 3, 10])

using Polyester
using Base: @propagate_inbounds

"""
    build_surface_grid(rr; nx::Int=0, ny::Int=0)

2D-s (alsó/felső) felülethez rács-alapú térindex. A rács celláiba az elemek
xy-bbox alapján kerülnek. Visszaad egy `grid` named tuple-t.
"""
function build_surface_grid(rr; nx::Int=0, ny::Int=0)
    gmsh.model.setCurrent(rr.model.name)

    # -- egyszeri beolvasás Gmsh-ból
    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    Xf = Float32.(coords[1:3:end]); Yf = Float32.(coords[2:3:end]); Zf = Float32.(coords[3:3:end])
    node_index = Dict{Int32,Int32}(Int32(tag)=>Int32(i) for (i,tag) in enumerate(nodeTags))

    # csak a 2D-s mező elemei (tri/quad)
    elem_tags = Int.(collect(rr.numElem))
    isempty(elem_tags) && error("surface grid: üres mező")
    _, _, dim, _ = gmsh.model.mesh.getElement(elem_tags[1])
    dim == 2 || error("surface grid: a rr mező 2D felületen legyen (dim=2)")

    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)

    # --- elemlisták (csak azok, amelyek a rr.numElem-ben vannak)
    elem_nodes = Vector{Vector{Int32}}()
    elem_types = Int[]
    elem_coords = Vector{NTuple{3,Vector{Float32}}}()
    elem_bboxes = NTuple{4,Float32}[]

    # bbox globális tartomány XY-ben
    gxmin = typemax(Float32); gxmax = typemin(Float32)
    gymin = typemax(Float32); gymax = typemin(Float32)

    rrset = Set(rr.numElem)
    for (t, et) in enumerate(elemTypes)
        # csak tri3/tri6/quad4/quad9
        et in (2,9,3,10) || continue
        # elem tulajdonság: csomópont/elem
        _, _, _, nper, _, _ = gmsh.model.mesh.getElementProperties(et)
        tags = elemTags[t]; nlist = elemNodeTags[t]
        @inbounds for j in 1:length(tags)
            tag = tags[j]
            tag ∈ rrset || continue
            offs = (j-1)*nper + 1
            nodes = Int32.(nlist[offs:offs+nper-1])
            xs = Vector{Float32}(undef, nper)
            ys = similar(xs); zs = similar(xs)
            for (k,nid) in enumerate(nodes)
                idx = node_index[nid]
                xs[k] = Xf[idx]; ys[k] = Yf[idx]; zs[k] = Zf[idx]
            end
            xmin = minimum(xs); xmax = maximum(xs)
            ymin = minimum(ys); ymax = maximum(ys)
            push!(elem_nodes, nodes)
            push!(elem_types, et)
            push!(elem_coords, (xs,ys,zs))
            push!(elem_bboxes, (xmin, xmax, ymin, ymax))
            gxmin = min(gxmin, xmin); gxmax = max(gxmax, xmax)
            gymin = min(gymin, ymin); gymax = max(gymax, ymax)
        end
    end

    ne = length(elem_nodes)
    ne == length(rr.A) || error("surface grid: eltér az elemek és rr.A száma")

    # --- rács felbontás (heuriszta: ~10-20 elem/bin átlag)
    if nx == 0 || ny == 0
        sidex = gxmax - gxmin
        sidey = gymax - gymin
        # cél elemsűrűség ~ 16 elem/bin
        bins = max(1, round(Int, sqrt(ne/16)))
        ar = sidex / max(sidey, eps(Float32))
        nx = max(8, round(Int, bins*sqrt(Float64(ar))))
        ny = max(8, round(Int, bins/sqrt(Float64(ar))))
    end
    dx = (gxmax - gxmin + eps(Float32)) / nx
    dy = (gymax - gymin + eps(Float32)) / ny

    # --- bin lista: minden cellához elemindexek
    bins = [Int32[] for _ in 1:(nx*ny)]

    @inline cellid(ix,iy) = (iy-1)*nx + ix

    # elemek bepakolása a lefedett bin(ek)be
    for i in 1:ne
        (xmin,xmax,ymin,ymax) = elem_bboxes[i]
        ix0 = max(1, floor(Int, (xmin-gxmin)/dx)+1)
        ix1 = min(nx,  ceil(Int, (xmax-gxmin)/dx)+1)
        iy0 = max(1, floor(Int, (ymin-gymin)/dy)+1)
        iy1 = min(ny,  ceil(Int, (ymax-gymin)/dy)+1)
        for iy in iy0:iy1, ix in ix0:ix1
            push!(bins[cellid(ix,iy)], Int32(i))
        end
    end

    return (
        elem_nodes=elem_nodes,
        elem_types=elem_types,
        elem_coords=elem_coords,
        elem_bboxes=elem_bboxes,
        # rács
        bins=bins, nx=nx, ny=ny,
        gxmin=gxmin, gymin=gymin, dx=dx, dy=dy
    )
end


"""
    probe_field_bulk_binned(rr, points, grid) -> Vector{Float64}

Rács-indexelt, többszálú lekérdezés. `points` alakja: N×3 (Float64).
Csak xy-ban keresünk (vékony szilárd), majd a rr (2D) alakfüggvényeivel interpolálunk.
"""
function probe_field_bulk_binned(rr, points::AbstractMatrix{<:Real}, grid)
    gmsh.model.setCurrent(rr.model.name)
    npts = size(points,1)
    vals = fill(NaN, npts)

    elem_nodes   = grid.elem_nodes
    elem_types   = grid.elem_types
    elem_coords  = grid.elem_coords
    bins         = grid.bins
    nx=grid.nx; ny=grid.ny
    gxmin=grid.gxmin; gymin=grid.gymin; dx=grid.dx; dy=grid.dy

    @inline function bin_of_point(x::Float64, y::Float64)
        ix = clamp(Int(floor((x - gxmin) / dx)) + 1, 1, nx)
        iy = clamp(Int(floor((y - gymin) / dy)) + 1, 1, ny)
        return ix, iy
    end

    #@inline function bin_of_point(x::Float64, y::Float64)
    #    ix = clamp(fld(Int, floor((x - gxmin)/dx), 1) + 1, 1, nx)
    #    iy = clamp(fld(Int, floor((y - gymin)/dy), 1) + 1, 1, ny)
    #    return ix, iy
    #end
    @inline function neigh_bins(ix::Int, iy::Int)
        # aktuális + szomszédos 8 cella (robosztus a bbox határán)
        I = (max(1,ix-1)):(min(nx,ix+1))
        J = (max(1,iy-1)):(min(ny,iy+1))
        return I, J
    end
    @inline cellid(ix,iy) = (iy-1)*nx + ix

    # --- TÖBBSZÁLÚ kiértékelés
    @batch for ipt in 1:npts
        x = Float64(points[ipt,1]); y = Float64(points[ipt,2]); z = Float64(points[ipt,3])
        ix, iy = bin_of_point(x,y)
        I, J = neigh_bins(ix,iy)

        found = false
        @inbounds for j in J, i in I
            for ei in bins[cellid(i,j)]
                et = elem_types[ei]
                xs, ys, zs = elem_coords[ei]
                ξηζ = natural_coordinates_generic(x, y, z, xs, ys, zs, et)
                ξηζ === nothing && continue
                ξ, η, ζ = ξηζ
                N = lagrange_shape_generic(ξ, η, ζ, et)
                length(N) == length(elem_nodes[ei]) || continue
                # mezőérték az adott elemhez rr.A[ei]
                vals[ipt] = dot(N, rr.A[ei][:,1])    # Float64 eredmény
                found = true
                break
            end
            found && break
        end
        # ha nem talált elemet: vals[ipt] marad NaN
    end

    return vals
end

"""
    build_bbox(rr::ScalarField)

Előkészíti az elemek xy-síkbeli bounding boxait, és visszaad egy cache-t,
amelyet a `probe_field_bulk` többször is fel tud használni.
A visszatérési érték egy named tuple a következő kulcsokkal:
    * `bbox` – kibővített (xmin, xmax, ymin, ymax) tartományok
    * `elem_nodes` – a mezőhöz tartozó elemek csomópont-listái
    * `elem_types` – a hozzájuk tartozó Gmsh elemtípusok
    * `elem_coords` – előre kigyűjtött (x, y, z) csomópont-koordináták
    * `elem_dims` – az elem topológiai dimenziója (1 vagy 2)
"""
function build_bbox(rr::ScalarField)
    gmsh.model.setCurrent(rr.model.name)

    elem_tags = Int.(collect(rr.numElem))
    if isempty(elem_tags)
        return (
            bbox=NTuple{4,Float64}[],
            elem_nodes=Vector{Vector{Int}}(),
            elem_types=Int[],
            elem_coords=Vector{NTuple{3,Vector{Float64}}}(),
            elem_dims=Int[],
        )
    end

    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    xcoords = coords[1:3:end]
    ycoords = coords[2:3:end]
    zcoords = coords[3:3:end]
    node_index = Dict{Int,Int}()
    for (pos, tag_raw) in enumerate(nodeTags)
        node_index[Int(tag_raw)] = pos
    end

    elem_nodes = Vector{Vector{Int}}(undef, length(elem_tags))
    elem_types = Vector{Int}(undef, length(elem_tags))
    elem_coords = Vector{NTuple{3,Vector{Float64}}}(undef, length(elem_tags))
    elem_dims = Vector{Int}(undef, length(elem_tags))
    elem_bboxes = Vector{NTuple{4,Float64}}(undef, length(elem_tags))

    for (i, tag) in enumerate(elem_tags)
        et, node_ids, dim, _ = gmsh.model.mesh.getElement(tag)
        et in SUPPORTED_ELEMENT_TYPES || error("build_bbox: element type $et not supported (only line/tri/quad)")
        dim in (1, 2) || error("build_bbox: expected 1D/2D mesh, got dim=$dim for element $tag")

        node_vec = Int.(node_ids)
        elem_nodes[i] = node_vec
        elem_types[i] = et
        elem_dims[i] = dim

        xs = Vector{Float64}(undef, length(node_vec))
        ys = similar(xs)
        zs = similar(xs)
        for (k, nid) in enumerate(node_vec)
            idx = get(node_index, nid, nothing)
            idx === nothing && error("build_bbox: node $nid missing from gmsh node list")
            xs[k] = xcoords[idx]
            ys[k] = ycoords[idx]
            zs[k] = zcoords[idx]
        end
        elem_coords[i] = (xs, ys, zs)

        xmin = minimum(xs)
        xmax = maximum(xs)
        ymin = minimum(ys)
        ymax = maximum(ys)
        span = max(xmax - xmin, ymax - ymin)
        pad = max(1e-9, 1e-6 * span)
        elem_bboxes[i] = (xmin - pad, xmax + pad, ymin - pad, ymax + pad)
    end

    return (
        bbox=elem_bboxes,
        elem_nodes=elem_nodes,
        elem_types=elem_types,
        elem_coords=elem_coords,
        elem_dims=elem_dims,
    )
end

"""
    probe_field_bulk_fallback(rr::ScalarField, points::AbstractMatrix{Float64})

Gyors, tisztán Julia-alapú tömeges lekérdezés (`(x,y,z)` pontok → mezőértékek)
1D (line2/line3) és 2D (tri3/tri6/quad4/quad9) elemekre.
A kereséshez csak az `(x, y)` koordinátákat használjuk, a mezőt ezt
követően interpoláljuk a lokális természetes koordináták alapján.
"""
function probe_field_bulk_fallback(rr::ScalarField, points::AbstractMatrix{Float64})
    gmsh.model.setCurrent(rr.model.name)
    npts = size(points, 1)
    vals = fill(NaN, npts)

    cache = build_bbox(rr)
    elem_bboxes = cache.bbox
    elem_types = cache.elem_types
    elem_coords = cache.elem_coords
    elem_nodes = cache.elem_nodes

    nelem = length(elem_bboxes)
    nelem == length(rr.A) || error("probe_field_bulk: element/DOF mismatch")

    for (i, nodes) in enumerate(elem_nodes)
        size(rr.A[i], 1) == length(nodes) || error("probe_field_bulk: coefficient size mismatch at element index $i")
    end

    @inbounds for ipt in 1:npts
        x = points[ipt, 1]
        y = points[ipt, 2]
        z = points[ipt, 3]

        for i in 1:nelem
            xmin, xmax, ymin, ymax = elem_bboxes[i]
            if x < xmin || x > xmax || y < ymin || y > ymax
                continue
            end

            xs, ys, zs = elem_coords[i]
            ξηζ = natural_coordinates_generic(x, y, z, xs, ys, zs, elem_types[i])
            ξηζ === nothing && continue
            ξ, η, ζ = ξηζ

            N = lagrange_shape_generic(ξ, η, ζ, elem_types[i])
            length(N) == length(elem_nodes[i]) || continue

            vals[ipt] = dot(N, rr.A[i][:, 1])
            break
        end
    end

    return vals
end



"""
    probe_field(rr::ScalarField, x, y, z)

Gyors, tisztán Julia-alapú lekérdezés (Gmsh nélkül), amely visszaadja
a `ScalarField` értékét az `(x, y, z)` pontban.

Támogatott elemtípusok: `line2`, `line3`, `tri3`, `tri6`, `quad4`, `quad9`, `tet4`, `tet10`, `hex8`, `hex20`.
"""
function probe_field(rr::ScalarField, x::Float64, y::Float64, z::Float64)
    gmsh.model.setCurrent(rr.model.name)

    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    X = coords[1:3:end]
    Y = coords[2:3:end]
    Z = coords[3:3:end]

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

            node_ids = elemNodeTags[i][(j-1)*numNodes+1:j*numNodes]
            xnodes = X[node_ids]
            ynodes = Y[node_ids]
            znodes = Z[node_ids]

            # gyors bbox szűrés
            if x < minimum(xnodes) - 1e-8 || x > maximum(xnodes) + 1e-8 ||
               y < minimum(ynodes) - 1e-8 || y > maximum(ynodes) + 1e-8 ||
               z < minimum(znodes) - 1e-8 || z > maximum(znodes) + 1e-8
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

"""
    probe_field_bulk(rr::ScalarField, points::Matrix{<:Real}; grid=nothing)

Ha a `grid` argumentum nincs megadva, a hagyományos (lassabb) módszert használja,
amely minden elemet végigvizsgál.  
Ha `grid` meg van adva (`build_surface_grid` kimenete), akkor a gyors,
rács-indexelt (Polyester) verziót hívja.
"""
function probe_field_bulk(rr::ScalarField, points::AbstractMatrix{<:Real}; grid=nothing)
    if grid === nothing
    @info "probe_field_bulk: running legacy slow mode"

    # --- építsük fel a bbox adatokat (mint eddig a build_bbox-ban) ---
    bboxdata = build_bbox(rr)
    elem_nodes = bboxdata.elem_nodes
    elem_types = bboxdata.elem_types
    elem_coords = bboxdata.elem_coords

    npts = size(points, 1)
    vals = fill(NaN, npts)

    for ipt in 1:npts
        x, y, z = points[ipt, 1], points[ipt, 2], points[ipt, 3]
        found = false
        for (ei, et) in enumerate(elem_types)
            xs, ys, zs = elem_coords[ei]
            ξηζ = natural_coordinates_generic(x, y, z, xs, ys, zs, et)
            ξηζ === nothing && continue
            ξ, η, ζ = ξηζ
            N = lagrange_shape_generic(ξ, η, ζ, et)
            if length(N) == length(elem_nodes[ei])
                vals[ipt] = dot(N, rr.A[ei][:,1])
                found = true
                break
            end
        end
        found || (vals[ipt] = NaN)
    end

    return vals
    #=
    if grid === nothing
        # --- fallback a régi, lassabb módszerre ---
        gmsh.model.setCurrent(rr.model.name)
        npts = size(points, 1)
        vals = fill(NaN, npts)
        elem_nodes, elem_types, elem_coords = rr.elem_nodes, rr.elem_types, rr.elem_coords

        @info "probe_field_bulk: running legacy slow mode on $(length(elem_nodes)) elements"

        for ipt in 1:npts
            x, y, z = points[ipt, 1], points[ipt, 2], points[ipt, 3]
            found = false
            for (ei, et) in enumerate(elem_types)
                xs, ys, zs = elem_coords[ei]
                ξηζ = natural_coordinates_generic(x, y, z, xs, ys, zs, et)
                ξηζ === nothing && continue
                ξ, η, ζ = ξηζ
                N = lagrange_shape_generic(ξ, η, ζ, et)
                if length(N) == length(elem_nodes[ei])
                    vals[ipt] = dot(N, rr.A[ei][:,1])
                    found = true
                    break
                end
            end
            found || (vals[ipt] = NaN)
        end

        return vals
        =#

    else
        # --- gyors, binninges verzió ---
        @info "probe_field_bulk: using binned grid with $(length(grid.elem_nodes)) elements and $(length(grid.bins)) bins"
        return probe_field_bulk_binned(rr, points, grid)
    end
end

# ---------------------------------------------------------------------
# Shape függvények
# ---------------------------------------------------------------------

function lagrange_shape_generic(ξ, η, ζ, et)
    if et == 1       # line2
        return 0.5 .* [(1 - ξ), (1 + ξ)]

    elseif et == 8   # line3
        return [
            0.5 * ξ * (ξ - 1),
            1 - ξ^2,
            0.5 * ξ * (ξ + 1)
        ]

    elseif et == 4   # tet4
        return [1 - ξ - η - ζ, ξ, η, ζ]

    elseif et == 5   # hex8
        return 0.125 .* [
            (1 - ξ) * (1 - η) * (1 - ζ),
            (1 + ξ) * (1 - η) * (1 - ζ),
            (1 + ξ) * (1 + η) * (1 - ζ),
            (1 - ξ) * (1 + η) * (1 - ζ),
            (1 - ξ) * (1 - η) * (1 + ζ),
            (1 + ξ) * (1 - η) * (1 + ζ),
            (1 + ξ) * (1 + η) * (1 + ζ),
            (1 - ξ) * (1 + η) * (1 + ζ)
        ]

    elseif et == 11  # tet10
        L1 = 1 - ξ - η - ζ
        L2 = ξ
        L3 = η
        L4 = ζ
        return [
            L1 * (2L1 - 1),
            L2 * (2L2 - 1),
            L3 * (2L3 - 1),
            L4 * (2L4 - 1),
            4 * L1 * L2,
            4 * L2 * L3,
            4 * L3 * L1,
            4 * L1 * L4,
            4 * L2 * L4,
            4 * L3 * L4
        ]

    elseif et == 17  # hex20
        N = Vector{Float64}(undef, 20)
        ξi = [-1, 1, 1, -1, -1, 1, 1, -1]
        ηi = [-1, -1, 1, 1, -1, -1, 1, 1]
        ζi = [-1, -1, -1, -1, 1, 1, 1, 1]
        for n in 1:8
            N[n] = 0.125 * (1 + ξi[n] * ξ) * (1 + ηi[n] * η) * (1 + ζi[n] * ζ) * (ξi[n] * ξ + ηi[n] * η + ζi[n] * ζ - 2)
        end
        N[9] = 0.25 * (1 - ξ^2) * (1 - η) * (1 - ζ)
        N[10] = 0.25 * (1 - η^2) * (1 + ξ) * (1 - ζ)
        N[11] = 0.25 * (1 - ξ^2) * (1 + η) * (1 - ζ)
        N[12] = 0.25 * (1 - η^2) * (1 - ξ) * (1 - ζ)
        N[13] = 0.25 * (1 - ξ^2) * (1 - η) * (1 + ζ)
        N[14] = 0.25 * (1 - η^2) * (1 + ξ) * (1 + ζ)
        N[15] = 0.25 * (1 - ξ^2) * (1 + η) * (1 + ζ)
        N[16] = 0.25 * (1 - η^2) * (1 - ξ) * (1 + ζ)
        N[17] = 0.25 * (1 - ζ^2) * (1 - ξ) * (1 - η)
        N[18] = 0.25 * (1 - ζ^2) * (1 + ξ) * (1 - η)
        N[19] = 0.25 * (1 - ζ^2) * (1 + ξ) * (1 + η)
        N[20] = 0.25 * (1 - ζ^2) * (1 - ξ) * (1 + η)
        return N

    elseif et == 2   # tri3
        return [1 - ξ - η, ξ, η]

    elseif et == 9   # tri6
        L1 = 1 - ξ - η
        L2 = ξ
        L3 = η
        return [L1 * (2L1 - 1), L2 * (2L2 - 1), L3 * (2L3 - 1), 4 * L1 * L2, 4 * L2 * L3, 4 * L3 * L1]

    elseif et == 3   # quad4
        return 0.25 .* [
            (1 - ξ) * (1 - η),
            (1 + ξ) * (1 - η),
            (1 + ξ) * (1 + η),
            (1 - ξ) * (1 + η)
        ]

    elseif et == 10  # quad9
        return [
            0.25 * ξ * (ξ - 1) * η * (η - 1),
            0.25 * ξ * (ξ + 1) * η * (η - 1),
            0.25 * ξ * (ξ + 1) * η * (η + 1),
            0.25 * ξ * (ξ - 1) * η * (η + 1),
            0.5 * (1 - ξ^2) * η * (η - 1),
            0.5 * ξ * (ξ + 1) * (1 - η^2),
            0.5 * (1 - ξ^2) * η * (η + 1),
            0.5 * ξ * (ξ - 1) * (1 - η^2),
            (1 - ξ^2) * (1 - η^2)
        ]

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
            znodes[2]-znodes[1] znodes[3]-znodes[1] znodes[4]-znodes[1]
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
            dx = x - xe
            dy = y - ye
            dz = z - ze
            if abs(dx) + abs(dy) + abs(dz) < 1e-10
                break
            end
            ξ += dx / (maximum(xnodes) - minimum(xnodes) + eps())
            η += dy / (maximum(ynodes) - minimum(ynodes) + eps())
            ζ += dz / (maximum(znodes) - minimum(znodes) + eps())
        end
        if abs(ξ) <= 1.01 && abs(η) <= 1.01 && abs(ζ) <= 1.01
            return (ξ, η, ζ)
        else
            return nothing
        end

    elseif et in (2, 9)  # tri3, tri6
        A = [
            xnodes[2]-xnodes[1] xnodes[3]-xnodes[1];
            ynodes[2]-ynodes[1] ynodes[3]-ynodes[1]
        ]
        b = [x - xnodes[1], y - ynodes[1]]
        ξη = A \ b
        ξ, η = ξη
        if ξ >= -1e-8 && η >= -1e-8 && ξ + η <= 1.0 + 1e-8
            return (ξ, η, 0.0)
        else
            return nothing
        end

    elseif et in (1, 8)  # line2, line3
        x1 = xnodes[1]
        y1 = ynodes[1]
        x2 = xnodes[end]
        y2 = ynodes[end]
        dx = x2 - x1
        dy = y2 - y1
        len2 = dx * dx + dy * dy
        len2 < 1e-20 && return nothing
        t = ((x - x1) * dx + (y - y1) * dy) / len2
        if t < -1e-8 || t > 1 + 1e-8
            return nothing
        end
        ξ = 2t - 1
        return (ξ, 0.0, 0.0)

    elseif et in (3, 10)  # quad4, quad9
        return quad_reference_coordinates(x, y, xnodes, ynodes, et)

    else
        return nothing
    end
end


function quad_reference_coordinates(x, y, xnodes, ynodes, et)
    ξ = 0.0
    η = 0.0
    tol = 1e-10
    fd = 1e-6

    for _ in 1:15
        N = lagrange_shape_generic(ξ, η, 0.0, et)
        xe = dot(N, xnodes)
        ye = dot(N, ynodes)
        rx = x - xe
        ry = y - ye
        if abs(rx) + abs(ry) < tol
            break
        end

        dN_dξ = (lagrange_shape_generic(ξ + fd, η, 0.0, et) .- lagrange_shape_generic(ξ - fd, η, 0.0, et)) ./ (2fd)
        dN_dη = (lagrange_shape_generic(ξ, η + fd, 0.0, et) .- lagrange_shape_generic(ξ, η - fd, 0.0, et)) ./ (2fd)

        dxdξ = dot(dN_dξ, xnodes)
        dydξ = dot(dN_dξ, ynodes)
        dxdη = dot(dN_dη, xnodes)
        dydη = dot(dN_dη, ynodes)

        detJ = dxdξ * dydη - dxdη * dydξ
        abs(detJ) < 1e-14 && break

        Δξ = (rx * dydη - ry * dxdη) / detJ
        Δη = (-rx * dydξ + ry * dxdξ) / detJ

        ξ += Δξ
        η += Δη

        if abs(Δξ) + abs(Δη) < tol
            break
        end
    end

    if abs(ξ) <= 1.0 + 1e-6 && abs(η) <= 1.0 + 1e-6
        return (ξ, η, 0.0)
    else
        return nothing
    end
end










#####################################################################################################
#=
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
=#

end # module

