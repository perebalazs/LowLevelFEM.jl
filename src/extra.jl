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
            for j in 1:length(elemTags[i])
                elem = elemTags[i][j]
                for l in 1:numPrimaryNodes
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
    tag = getTagForPhysicalName(problem.material[1].phName)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(problem.dim, tag)
    coord = reshape(coord, 3, :)
    #val = probe_field_bulk_fallback(h0, coord')
    grid = build_surface_grid(h0, nx=500, ny=500)
    val = probe_field_bulk_walking(h0, coord', grid)
    h1 = scalarField(problem, problem.material[1].phName, (x,y,z)->z)
    h1.a[nodeTags] .-= val
    h1.a .*= -1.0
    h = h1 # nodesToElements(h1, onPhysicalGroup="bottom")

    #fh(x, y, z) = gmsh.view.probe(problem.geometry.tagTop, x, y, 0)[1][1] - z
    #h = scalarField(problem, problem.material[1].phName, fh)
    grad_h = ∇(h)
    ex = VectorField(problem, problem.material[1].phName, [1,0,0])
    dhdx = grad_h * ex
    problem.geometry.h = h
    problem.geometry.dhdx = elementsToNodes(dhdx)
    #gmsh.view.remove(problem.geometry.tagTop)
end

function pressureInVolume(p::ScalarField)
    if p.model.geometry.nameVolume == ""
        error("initializePressure: no volume for lubricant has been defined.")
    end

    tag = LowLevelFEM.getTagForPhysicalName(p.model.geometry.nameVolume)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(p.model.dim + 1, tag)
    coord = reshape(coord, 3, :)
    coord_mat = Matrix(coord')
    p1 = nodesToElements(p)
    grid = build_surface_grid(p1, nx=500, ny=500)
    #val = probe_field_bulk(p1, coord_mat, grid=grid)
    #val = probe_field_bulk_fallback(p1, coord_mat)
    val = probe_field_bulk_walking(p1, coord_mat, grid)
    pp = scalarField(p.model, p.model.geometry.nameVolume, 0)
    pp.a[nodeTags] .= val

    if isnothing(p.model.geometry.hh)
        #tag = LowLevelFEM.getTagForPhysicalName(p.model.geometry.nameVolume)
        #nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(p.model.dim + 1, tag)
        #coord = reshape(coord, 3, :)
        h_top = ScalarField(p.model, p.model.geometry.nameGap, (x, y, z) -> z)
        grid = build_surface_grid(h_top, nx=500, ny=500)
        #val = probe_field_bulk(h_top, coord_mat, grid=grid)
        #val = probe_field_bulk_fallback(h_top, coord_mat)
        val = probe_field_bulk_walking(h_top, coord_mat, grid)
        hh = scalarField(p.model, p.model.geometry.nameVolume, 0)
        hh.a[nodeTags] .= val
        hh = nodesToElements(hh, onPhysicalGroup=p.model.geometry.nameVolume)
        p.model.geometry.hh = nodesToElements(hh)
    end

    return pp

    #=
    return nodesToElements(pp, onPhysicalGroup=p.model.geometry.nameVolume)

    if p.model.geometry.nameVolume == ""
        error("initializePressure: no volume for lubricant has been defined.")
    end
    if p.model.geometry.hh == nothing
        view_h = showDoFResults(p.model.geometry.h)
        fh(x, y, z) = gmsh.view.probe(view_h, x, y, 0, -1, -1, false, -1)[1][1]
        p.model.geometry.hh = ScalarField(p.model, [field(p.model.geometry.nameVolume, f=fh)])
        gmsh.view.remove(view_h)
    end
    view_p = showDoFResults(p)
    fp(x, y, z) = gmsh.view.probe(view_p, x, y, 0, -1, -1, false, -1)[1][1]
    pp = ScalarField(p.model, [field(p.model.geometry.nameVolume, f=fp)])
    gmsh.view.remove(view_p)
    return pp
    =#
end

#=
function initializePressure(p::ScalarField)
    if p.model.geometry.nameVolume == ""
        error("initializePressure: no volume for lubricant has been defined.")
    end
    if p.model.geometry.hh == nothing
        view_h = showDoFResults(p.model.geometry.h)
        fh(x, y, z) = gmsh.view.probe(view_h, x, y, 0, -1, -1, false, -1)[1][1]
        p.model.geometry.hh = ScalarField(p.model, [field(p.model.geometry.nameVolume, f=fh)])
        gmsh.view.remove(view_h)
    end
    view_p = showDoFResults(p)
    fp(x, y, z) = gmsh.view.probe(view_p, x, y, 0, -1, -1, false, -1)[1][1]
    p.model.geometry.p = ScalarField(p.model, [field(p.model.geometry.nameVolume, f=fp)])
    gmsh.view.remove(view_p)
end
=#

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
