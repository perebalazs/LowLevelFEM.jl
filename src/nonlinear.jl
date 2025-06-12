function nodePositionVector(problem)
    dim = problem.dim
    non = problem.non
    if dim != 3
        error("nodePositionVector: This function only works for 3D problems.")
    end
    r = zeros(non * dim)
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            r[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            r[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            r[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
        end
    end
    return r
end

function deformationGradient(problem, r)
    gmsh.model.setCurrent(problem.name)

    type = :F
    nsteps = size(r, 2)
    F = []
    numElem = Int[]
    ncoord2 = zeros(3 * problem.non)
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim != 3 || problem.type != :Solid
            error("deformationGradient: dimension is $(problem.dim), problem type is $(problem.type). (They must be 3 and :Solid.)")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                F1 = zeros(9numNodes, nsteps)
                nodeCoord = zeros(numNodes * 3)
                for k in 1:dim, j = 1:numNodes
                    nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
                end
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
                ∇h = reshape(dfun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numNodes)
                ∂h = zeros(dim, numNodes * numNodes)
                ∂H = spzeros(9 * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                sizehint!(∂H, 9 * numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numNodes
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numNodes, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r[nn2]
                    ∂H .*= 0
                    for k in 1:numNodes
                        for l in 1:numNodes
                            ∂H[k*9-(9-1), l*pdim-(pdim-1)] = ∂h[1, (k-1)*numNodes+l]
                            ∂H[k*9-(9-2), l*pdim-(pdim-1)] = ∂h[2, (k-1)*numNodes+l]
                            ∂H[k*9-(9-3), l*pdim-(pdim-1)] = ∂h[3, (k-1)*numNodes+l]
                            ∂H[k*9-(9-4), l*pdim-(pdim-2)] = ∂h[1, (k-1)*numNodes+l]
                            ∂H[k*9-(9-5), l*pdim-(pdim-2)] = ∂h[2, (k-1)*numNodes+l]
                            ∂H[k*9-(9-6), l*pdim-(pdim-2)] = ∂h[3, (k-1)*numNodes+l]
                            ∂H[k*9-(9-7), l*pdim-(pdim-3)] = ∂h[1, (k-1)*numNodes+l]
                            ∂H[k*9-(9-8), l*pdim-(pdim-3)] = ∂h[2, (k-1)*numNodes+l]
                            ∂H[k*9-(9-9), l*pdim-(pdim-3)] = ∂h[3, (k-1)*numNodes+l]
                        end
                    end
                    push!(numElem, elem)
                    for k in 1:numNodes
                        ∂H1 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                        for kk in 1:nsteps
                            F0 = ∂H1 * r[nn2, kk]
                            F1[(k-1)*9+1:k*9, kk] = [F0[1], F0[4], F0[7],
                                F0[2], F0[5], F0[8],
                                F0[3], F0[6], F0[9]]
                        end
                    end
                    push!(F, F1)
                end
            end
        end
    end
    sigma = TensorField(F, numElem, nsteps, type)
    return sigma
end

function stiffnessMatrixLinear(problem, r)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
    nn = []
    I = []
    J = []
    V = []
    V = convert(Vector{Float64}, V)
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 3 && problem.type == :Solid
            #=D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]
                =#
            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            rowsOfB = 3
            b = 1
        else
            error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                #comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                #h = reshape(fun, :, numIntPoints)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                #H = zeros(9 * numIntPoints, pdim * numNodes)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                ∂H = spzeros(9 * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                sizehint!(∂H, 9 * numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r[nn2]
                    B .*= 0
                    ∂H .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            end
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                ∂H[k*9-(9-1), l*pdim-(pdim-1)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-2), l*pdim-(pdim-2)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-3), l*pdim-(pdim-3)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-4), l*pdim-(pdim-1)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-5), l*pdim-(pdim-2)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-6), l*pdim-(pdim-3)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-7), l*pdim-(pdim-1)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-8), l*pdim-(pdim-2)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-9), l*pdim-(pdim-3)] = ∂h[3, (k-1)*numNodes+l]
                            end
                            ∂H1 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F = ∂H1 * r1
                            for l in 1:numNodes
                                B[k*rowsOfB-5, l*pdim-2] = F[1] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-2] = F[4] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-2] = F[7] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-2] = (F[4] * ∂h[1, (k-1)*numNodes+l] + F[1] * ∂h[2, (k-1)*numNodes+l]) / 2.0
                                B[k*rowsOfB-1, l*pdim-2] = (F[7] * ∂h[2, (k-1)*numNodes+l] + F[4] * ∂h[3, (k-1)*numNodes+l]) / 2.0
                                B[k*rowsOfB-0, l*pdim-2] = (F[1] * ∂h[3, (k-1)*numNodes+l] + F[7] * ∂h[1, (k-1)*numNodes+l]) / 2.0
                                B[k*rowsOfB-5, l*pdim-1] = F[2] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-1] = F[5] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-1] = F[8] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-1] = (F[5] * ∂h[1, (k-1)*numNodes+l] + F[2] * ∂h[2, (k-1)*numNodes+l]) / 2.0
                                B[k*rowsOfB-1, l*pdim-1] = (F[8] * ∂h[2, (k-1)*numNodes+l] + F[5] * ∂h[3, (k-1)*numNodes+l]) / 2.0
                                B[k*rowsOfB-0, l*pdim-1] = (F[2] * ∂h[3, (k-1)*numNodes+l] + F[8] * ∂h[1, (k-1)*numNodes+l]) / 2.0
                                B[k*rowsOfB-5, l*pdim-0] = F[3] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-0] = F[6] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-0] = F[9] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-0] = (F[6] * ∂h[1, (k-1)*numNodes+l] + F[3] * ∂h[2, (k-1)*numNodes+l]) / 2.0
                                B[k*rowsOfB-1, l*pdim-0] = (F[9] * ∂h[2, (k-1)*numNodes+l] + F[6] * ∂h[3, (k-1)*numNodes+l]) / 2.0
                                B[k*rowsOfB-0, l*pdim-0] = (F[1] * ∂h[3, (k-1)*numNodes+l] + F[9] * ∂h[1, (k-1)*numNodes+l]) / 2.0
                            end
                        end
                    else
                        error("stiffnessMatrixLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        Ex = 10
                        νxy = 0.3
                        λ = Ex * νxy / ((1 + νxy) * (1 - 2νxy))
                        μ = Ex / (2 * (1 + νxy))
                        #∂H1 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                        #F1 = ∂H1 * r1
                        #F2 = reshape(F1, 3, 3)
                        #F2 = [1.0 0 0; 0 1.0 0; 0 0 1.0]
                        #J1 = det(F2)
                        #C2 = F2' * F2
                        #iC2 = inv(C2)
                        #iC1 = reshape(iC2, 9, 1)
                        #iCiC1 = iC1[[1, 5, 9, 4, 6, 3]] * iC1[[1, 5, 9, 4, 6, 3]]'
                        #i1 = [1, 2, 3, 1, 2, 3]
                        #j1 = [1, 2, 3, 2, 3, 1]
                        #k1 = [1, 2, 3, 1, 2, 3]
                        #l1 = [1, 2, 3, 2, 3, 1]
                        #II1 = zeros(6, 6)
                        #for i2 in 1:6, j2 in 1:6
                        #    II1[i2, j2] = (iC2[i1[i2], k1[j2]] * iC2[j1[i2], l1[j2]] + iC2[i1[i2], l1[j2]] * iC2[j1[i2], k1[j2]]) / 2
                        #end
                        #C1 = λ * iCiC1 + 2 * (μ - λ * log(J1)) * II1
                        λ0 = λ
                        μ0 = μ
                        C1 = [λ0+2μ0 λ0 λ0 0 0 0;
                            λ0 λ0+2μ0 λ0 0 0 0;
                            λ0 λ0 λ0+2μ0 0 0 0;
                            0 0 0 μ0 0 0;
                            0 0 0 0 μ0 0;
                            0 0 0 0 0 μ0]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += B1' * C1 * B1 * b * jacDet[k] * intWeights[k]
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return K
end

function stiffnessMatrixNonLinear(problem, r)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
    nn = []
    I = []
    J = []
    V = []
    V = convert(Vector{Float64}, V)
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 3 && problem.type == :Solid
            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            rowsOfB = 3
            b = 1
        else
            error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                H = spzeros(9 * numIntPoints, 9 * numNodes)
                sizehint!(H, 9 * numNodes * numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                ∂H = spzeros(9 * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, dim * numNodes)
                nn9 = zeros(Int, 9 * numNodes)
                sizehint!(∂H, 9 * numNodes * numIntPoints)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:9
                        H[k*9-(9-kk), l*9-(9-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    for k in 1:9
                        nn9[k:9:9*numNodes] = 9 * nnet[j, 1:numNodes] .- (9 - k)
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r[nn2]
                    #S1 = S[nn9]
                    B .*= 0
                    ∂H .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            end
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                ∂H[k*9-(9-1), l*pdim-(pdim-1)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-2), l*pdim-(pdim-2)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-3), l*pdim-(pdim-3)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-4), l*pdim-(pdim-1)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-5), l*pdim-(pdim-2)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-6), l*pdim-(pdim-3)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-7), l*pdim-(pdim-1)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-8), l*pdim-(pdim-2)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-9), l*pdim-(pdim-3)] = ∂h[3, (k-1)*numNodes+l]
                            end
                        end
                    else
                        error("stiffnessMatrixLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        Ex = 10
                        νxy = 0.3
                        λ = Ex * νxy / ((1 + νxy) * (1 - 2νxy))
                        μ = Ex / (2 * (1 + νxy))
                        I3 = [1 0 0; 0 1 0; 0 0 1]
                        #H1 = H[k*9-(9-1):k*9, 1:9*numNodes]
                        ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                        F1 = ∂H0 * r1
                        F2 = reshape(F1, 3, 3)
                        #J1 = det(F2)
                        E = (F2' * F2 - I3) / 2
                        #iC = inv(C)
                        S2 = 2μ * E + λ * (E[1, 1] + E[2, 2] + E[3, 3]) * I3
                        #S0 = reshape(H1 * S1, 3, 3)
                        ∂H1 = ∂H[k*9-(9-1):3:k*9-(9-9), 1:pdim*numNodes]
                        ∂H2 = ∂H[k*9-(9-2):3:k*9-(9-9), 1:pdim*numNodes]
                        ∂H3 = ∂H[k*9-(9-3):3:k*9-(9-9), 1:pdim*numNodes]
                        K1 += (∂H1' * S2 * ∂H1 + ∂H2' * S2 * ∂H2 + ∂H3' * S2 * ∂H3) * b * jacDet[k] * intWeights[k]
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return K
end

function loadVectorNonLinear(problem, r)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    f = zeros(problem.non * problem.pdim)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 3 && problem.type == :Solid
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            rowsOfB = 3
            b = 1
        else
            error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                H = spzeros(9 * numIntPoints, 9 * numNodes)
                sizehint!(H, 9 * numNodes * numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                ∂H = spzeros(9 * numIntPoints, pdim * numNodes)
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn9 = zeros(Int, 9 * numNodes)
                sizehint!(∂H, 9 * numNodes * numIntPoints)
                #for k in 1:numIntPoints, l in 1:numNodes
                #    for kk in 1:3:9, ll in 1:3
                #H[k*9-(9-kk)-1+ll, l*pdim-(pdim-ll)] = h[(k-1)*numNodes+l]
                #    end
                #end
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:9
                        nn9[k:9:9*numNodes] = 9 * nnet[j, 1:numNodes] .- (9 - k)
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r[nn2]
                    #S1 = S[nn9]
                    B .*= 0
                    ∂H .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            end
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints
                            for l in 1:numNodes
                                ∂H[k*9-(9-1), l*pdim-(pdim-1)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-2), l*pdim-(pdim-2)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-3), l*pdim-(pdim-3)] = ∂h[1, (k-1)*numNodes+l]
                                ∂H[k*9-(9-4), l*pdim-(pdim-1)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-5), l*pdim-(pdim-2)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-6), l*pdim-(pdim-3)] = ∂h[2, (k-1)*numNodes+l]
                                ∂H[k*9-(9-7), l*pdim-(pdim-1)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-8), l*pdim-(pdim-2)] = ∂h[3, (k-1)*numNodes+l]
                                ∂H[k*9-(9-9), l*pdim-(pdim-3)] = ∂h[3, (k-1)*numNodes+l]
                            end
                            ∂H1 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F = ∂H1 * r1
                            for l in 1:numNodes
                                B[k*rowsOfB-5, l*pdim-2] = F[1] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-2] = F[4] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-2] = F[7] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-2] = (F[4] * ∂h[1, (k-1)*numNodes+l] + F[1] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-2] = (F[7] * ∂h[2, (k-1)*numNodes+l] + F[4] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-2] = (F[1] * ∂h[3, (k-1)*numNodes+l] + F[7] * ∂h[1, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-5, l*pdim-1] = F[2] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-1] = F[5] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-1] = F[8] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-1] = (F[5] * ∂h[1, (k-1)*numNodes+l] + F[2] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-1] = (F[8] * ∂h[2, (k-1)*numNodes+l] + F[5] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-1] = (F[2] * ∂h[3, (k-1)*numNodes+l] + F[8] * ∂h[1, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-5, l*pdim-0] = F[3] * ∂h[1, (k-1)*numNodes+l]
                                B[k*rowsOfB-4, l*pdim-0] = F[6] * ∂h[2, (k-1)*numNodes+l]
                                B[k*rowsOfB-3, l*pdim-0] = F[9] * ∂h[3, (k-1)*numNodes+l]
                                B[k*rowsOfB-2, l*pdim-0] = (F[6] * ∂h[1, (k-1)*numNodes+l] + F[3] * ∂h[2, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-1, l*pdim-0] = (F[9] * ∂h[2, (k-1)*numNodes+l] + F[6] * ∂h[3, (k-1)*numNodes+l]) / 2
                                B[k*rowsOfB-0, l*pdim-0] = (F[1] * ∂h[3, (k-1)*numNodes+l] + F[9] * ∂h[1, (k-1)*numNodes+l]) / 2
                            end
                        end
                    else
                        error("stiffnessMatrixLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    f1 .*= 0
                    for k in 1:numIntPoints
                        Ex = 10
                        νxy = 0.3
                        λ = Ex * νxy / ((1 + νxy) * (1 - 2νxy))
                        μ = Ex / (2 * (1 + νxy))
                        I3 = [1 0 0; 0 1 0; 0 0 1]
                        #H1 = H[k*9-(9-1):k*9, 1:9*numNodes]
                        ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                        F1 = ∂H0 * r1
                        F2 = reshape(F1, 3, 3)
                        #J1 = det(F2)
                        E = (F2' * F2 - I3) / 2
                        #iC = inv(C)
                        S2 = 2μ * E + λ * (E[1, 1] + E[2, 2] + E[3, 3]) * I3
                        #S2 = [1 0 0; 0 0 0; 0 0 0]
                        #S2 = -I3
                        #S3 = reshape(S2, 9, 1)
                        #H1 = H[k*9-(9-1):k*9, 1:9*numNodes]
                        #S0 = H1 * S1
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        f1 += B1' * S2[[1, 5, 9, 4, 6, 3]] * b * jacDet[k] * intWeights[k]
                    end
                    f[nn2] += f1
                end
            end
        end
    end
    return f
end

function followerLoadVector(problem, r, loads)
    gmsh.model.setCurrent(problem.name)
    if !isa(loads, Vector)
        error("loadVector: loads are not arranged in a vector. Put them in [...]")
    end
    pdim = problem.pdim
    DIM = problem.dim
    b = problem.thickness
    non = problem.non
    dof = pdim * non
    fp = zeros(dof)
    ncoord2 = zeros(3 * problem.non)
    F0 = deformationGradient(problem, r)
    F0 = elementsToNodes(problem, F0)
    for n in 1:length(loads)
        name, fx, fy, fz = loads[n]
        if pdim == 3
            f = [.0, .0, .0]
        elseif pdim == 2
            f = [.0, .0]
        elseif pdim == 1
            f = [.0]
        else
            error("loadVector: dimension of the problem is $(problem.dim).")
        end
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
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elementTags[ii]), numNodes)
                invJac = zeros(3, 3numIntPoints)
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                HH = zeros(9 * numIntPoints, 9 * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                        end
                        for l in 1:9
                            HH[j*9-(9-l), k*9-(9-l)] = h[(j-1)*numNodes+k]
                        end
                    end
                end
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn9 = zeros(Int, 9 * numNodes)
                ∂h = zeros(pdim, numNodes * numIntPoints)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:9
                        nn9[k:9:9*numNodes] = 9 * nnet[l, 1:numNodes] .- (9 - k)
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    f1 .*= 0
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r[nn2]
                    F9 = F0[nn9]
                    for j in 1:numIntPoints
                        H9 = HH[j*9-(9-1):j*9, 1:9*numNodes]
                        F1 = H9 * F9
                        F1 = reshape(F1, 3, 3)
                        x = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 2]
                        y = 0
                        z = 0
                        if isa(fx, Function) || isa(fy, Function) || isa(fz, Function)
                            y = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 1]
                            z = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 0]
                        end
                        if fz == 2im
                            if isa(fx, Function)
                                error("heatConvectionVector: h cannot be a function.")
                            end
                            f[1] = isa(fy, Function) ? fx * fy(x, y, z) : fx * fy
                        else
                            f[1] = isa(fx, Function) ? fx(x, y, z) : fx
                        end
                        if pdim > 1
                            f[2] = isa(fy, Function) ? fy(x, y, z) : fy
                        end
                        if pdim == 3
                            f[3] = isa(fz, Function) ? fz(x, y, z) : fz
                        end
                        H1 = H[j*pdim-(pdim-1):j*pdim, 1:pdim*numNodes] # H1[...] .= H[...] ????
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
                            Ja = 2π * jacDet[j] * x
                        elseif DIM == 2 && dim == 1 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * x
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
                        else
                            error("loadVector: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                        end
                        f1 += H1' * F1 * f * Ja * intWeights[j]
                    end
                    fp[nn2] += f1
                end
            end
        end
    end
    return fp
end


