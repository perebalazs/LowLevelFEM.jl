export pressureConstraint, flowRate
export solvePressure


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
    fh(x, y, z) = gmsh.view.probe(problem.geometry.tagTop, x, y, 0)[1][1] - z
    h = scalarField(problem, problem.material[1].phName, fh)
    grad_h = ∇(h)
    ex = VectorField(problem, problem.material[1].phName, [1,0,0])
    dhdx = grad_h * ex
    problem.geometry.h = h
    problem.geometry.dhdx = elementsToNodes(dhdx)
    gmsh.view.remove(problem.geometry.tagTop)
end

#function systemMatrix(problem, αInNodes::ScalarField, velocity::Number, height::ScalarField)
function systemMatrix(problem, velocity::Number)
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
