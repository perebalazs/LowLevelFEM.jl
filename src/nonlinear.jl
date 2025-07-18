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
    return VectorField([], reshape(r, :,1), [0], [], 1, :u3D, problem)
end

function deformationGradient(problem, r; DoFResults=false)
    gmsh.model.setCurrent(problem.name)

    type = :F
    
    nsteps = r.nsteps
    ε = []
    numElem = Int[]
    ncoord2 = zeros(3 * problem.non)
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    type = :e
    if DoFResults == true
        E1 = zeros(non * 9, nsteps)
        pcs = zeros(Int64, non * dim)
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        ν = problem.material[ipg].ν
        dim = 0
        if problem.dim == 3 && problem.type == :Solid
            dim = 3
            rowsOfB = 9
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            dim = 2
            rowsOfB = 3
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            dim = 2
            rowsOfB = 3
            b = 1
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
            dim = 2
            rowsOfB = 4
            b = 1
        else
            error("solveStrain: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                #e0 = zeros(rowsOfB * numNodes)
                nodeCoord = zeros(numNodes * 3)
                for k in 1:dim, j = 1:numNodes
                    nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
                end
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
                ∇h = reshape(dfun, :, numNodes)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
                h = reshape(fun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                invJac = zeros(3, 3numNodes)
                ∂h = zeros(3, numNodes * numNodes)
                B = zeros(rowsOfB * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                rr = zeros(numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numNodes
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        rr[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numNodes, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*3-0, l*2-0] = B[k*3-2, l*2-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*3-0, l*2-1] = B[k*3-1, l*2-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 9
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*9-(9-1), l*3-(3-1)] = B[k*9-(9-2), l*3-(3-2)] = B[k*9-(9-3), l*3-(3-3)] = ∂h[1, (k-1)*numNodes+l]
                            B[k*9-(9-4), l*3-(3-1)] = B[k*9-(9-5), l*3-(3-2)] = B[k*9-(9-6), l*3-(3-3)] = ∂h[2, (k-1)*numNodes+l]
                            B[k*9-(9-7), l*3-(3-1)] = B[k*9-(9-8), l*3-(3-2)] = B[k*9-(9-9), l*3-(3-3)] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif dim == 2 && rowsOfB == 4
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*4-2, l*2-1] = rr[k] < 1e-10 ? 0 : h[l, k] / rr[k]
                        end
                    else
                        error("solveStrain: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    e = zeros(9numNodes, nsteps) # tensors have nine elements
                    for k in 1:numNodes
                        if rowsOfB == 9 && dim == 3 && problem.type == :Solid
                            B1 = B[k*9-8:k*9, 1:3*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[2], e0[3],
                                        e0[4], e0[5], e0[6],
                                        e0[7], e0[8], e0[9]]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[2], e0[3], e0[4], e0[5], e0[6], e0[7], e0[8], e0[9]]
                                end
                            end
                        elseif rowsOfB == 3 && dim == 2 && problem.type == :PlaneStress
                            B1 = B[k*3-2:k*3, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[3]/2, 0,
                                        e0[3]/2, e0[2], 0,
                                        0, 0, ν/(ν-1)*(e0[1]+e0[2])]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[3], 0, e0[3], e0[2], 0, 0, 0, ν/(ν-1)*(e0[1]+e0[2])]
                                end
                            end
                        elseif rowsOfB == 3 && dim == 2 && problem.type == :PlaneStrain
                            B1 = B[k*3-2:k*3, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[3]/2, 0,
                                        e0[3]/2, e0[2], 0,
                                        0, 0, 0]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[3], 0, e0[3], e0[2], 0, 0, 0, 0]
                                end
                            end
                        elseif rowsOfB == 4 && dim == 2 && problem.type == :AxiSymmetric
                            B1 = B[k*4-3:k*4, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[4]/2, 0,
                                        e0[4]/2, e0[3], 0,
                                        0, 0, e0[2]]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[4], 0, e0[4], e0[3], 0, 0, 0, e0[2]]
                                end
                            end
                        else
                            error("solveStrain: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
                        end
                    end
                    if DoFResults == true
                        pcs[nnet[j,1:numNodes]] .+= 1
                    end
                    if DoFResults == false
                        push!(ε, e)
                    end
                end
            end
        end
    end
    if DoFResults == true
        for k in 1:9
            for l in 1:non
                E1[k + 9 * l - 9, :] ./= pcs[l]
            end
        end
    end
    if DoFResults == true
        epsilon = TensorField([], E1, r.t, [], nsteps, type, problem)
        return epsilon
    else
        epsilon = TensorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        return epsilon
    end
end

function deformationGradient(r; DoFResults=false)
    return deformationGradient(r.model, r, DoFResults=DoFResults)
end

function gradientOf(r; DoFResults=false)
    return deformationGradient(r.model, r, DoFResults=DoFResults)
end

function ∇(r::VectorField; DoFResults=false)
    return deformationGradient(r.model, r, DoFResults=DoFResults)
end



#=
function deformationGradient(problem, r)
    gmsh.model.setCurrent(problem.name)

    type = :F
    nsteps = size(r.a, 2)
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
                    for k in 1:numNodes
                        for l in 1:numNodes
                            ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                        end
                    end
                    #r1 = r[nn2]
                    ∂H .*= 0
                    for k in 1:numNodes
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
                    push!(numElem, elem)
                    for k in 1:numNodes
                        ∂H1 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                        for kk in 1:nsteps
                            F0 = ∂H1 * r.a[nn2, kk]
                            F1[(k-1)*9+1:k*9, kk] = [F0[1], F0[2], F0[3],
                                F0[4], F0[5], F0[6],
                                F0[7], F0[8], F0[9]]
                        end
                    end
                    push!(F, F1)
                end
            end
        end
    end
    a = [;;]
    sigma = TensorField(F, a, r.t, numElem, nsteps, type)
    return sigma
end
=#

function tangentMatrixConstitutive(problem, r)
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
    II1 = zeros(6, 6)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        mtype = problem.material[ipg].type
        μ = problem.material[ipg].μ
        λ = problem.material[ipg].λ
        κ = problem.material[ipg].κ
        dim = problem.dim
        pdim = problem.pdim
        if mtype == :StVenantKirchhoff
            rowsOfB = 6
            b = 1
        elseif mtype == :NeoHookeCompressible
            rowsOfB = 6
            b = 1
        else
            error("stiffnessMatrixLinear: dimension is $(problem.dim), material type is $(mtype).")
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
                    r1 = r.a[nn2]
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
                                B[k*rowsOfB-0, l*pdim-0] = (F[3] * ∂h[3, (k-1)*numNodes+l] + F[9] * ∂h[1, (k-1)*numNodes+l]) / 2
                            end
                        end
                    else
                        error("stiffnessMatrixLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        if mtype == :StVenantKirchhoff
                            λ0 = λ
                            μ0 = μ
                            C1 = [λ0+2μ0 λ0 λ0 0 0 0;
                                λ0 λ0+2μ0 λ0 0 0 0;
                                λ0 λ0 λ0+2μ0 0 0 0;
                                0 0 0 μ0 0 0;
                                0 0 0 0 μ0 0;
                                0 0 0 0 0 μ0]
                        elseif mtype == :NeoHookeCompressible
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            J1 = det(F2)
                            C2 = F2' * F2
                            iC2 = inv(C2)
                            iC1 = reshape(iC2, 9, 1)
                            iCiC1 = iC1[[1, 5, 9, 4, 6, 3]] * iC1[[1, 5, 9, 4, 6, 3]]'
                            i1 = [1, 2, 3, 1, 2, 3]
                            j1 = [1, 2, 3, 2, 3, 1]
                            k1 = [1, 2, 3, 1, 2, 3]
                            l1 = [1, 2, 3, 2, 3, 1]
                            for i2 in 1:6, j2 in 1:6
                                II1[i2, j2] = (iC2[i1[i2], k1[j2]] * iC2[j1[i2], l1[j2]] + iC2[i1[i2], l1[j2]] * iC2[j1[i2], k1[j2]]) / 2
                            end
                            #J1 = J1 > 0 ? J1 : 1e-14
                            C1 = λ * iCiC1 + 2 * (μ - λ * log(J1)) * II1
                        else
                            error("loadVectorNonLinear: type '$(problem.type)' is not defined.")
                        end
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

function tangentMatrixConstitutive(r)
    return tangentMatrixConstitutive(r.model, r)
end

function tangentMatrixInitialStress(problem, r)
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
        mtype = problem.material[ipg].type
        μ = problem.material[ipg].μ
        λ = problem.material[ipg].λ
        κ = problem.material[ipg].κ
        dim = problem.dim
        pdim = problem.pdim
        if mtype == :StVenantKirchhoff
            rowsOfB = 6
            b = 1
        elseif mtype == :NeoHookeCompressible
            rowsOfB = 6
            b = 1
        else
            error("stiffnessMatrixNonLinear: dimension is $(problem.dim), material type is $(mtype).")
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
                    r1 = r.a[nn2]
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
                        error("stiffnessMatrixNonLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        if mtype == :StVenantKirchhoff
                            I3 = [1 0 0; 0 1 0; 0 0 1]
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            E = (F2' * F2 - I3) / 2
                            S2 = 2μ * E + λ * (E[1, 1] + E[2, 2] + E[3, 3]) * I3
                        elseif mtype == :NeoHookeCompressible
                            I3 = [1 0 0; 0 1 0; 0 0 1]
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            J1 = det(F2)
                            C = F2' * F2
                            iC = inv(C)
                            #J1 = J1 > 0 ? J1 : 1e-14
                            S2 = μ * (I3 - iC) + λ * log(J1) * iC
                        else
                            error("loadVectorNonLinear: type '$(problem.type)' is not defined.")
                        end
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

function tangentMatrixInitialStress(r)
    return tangentMatrixInitialStress(r.model, r)
end

function equivalentNodalForce(problem, r)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    f = zeros(problem.non * problem.pdim)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        mtype = problem.material[ipg].type
        μ = problem.material[ipg].μ
        λ = problem.material[ipg].λ
        κ = problem.material[ipg].κ
        dim = problem.dim
        pdim = problem.pdim
        if mtype == :StVenantKirchhoff
            rowsOfB = 6
            b = 1
        elseif mtype == :NeoHookeCompressible
            rowsOfB = 6
            b = 1
        else
            error("loadVectorNonLinear: dimension is $(problem.dim), material type is $(problem.type).")
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
                    #for k in 1:9
                    #    nn9[k:9:9*numNodes] = 9 * nnet[j, 1:numNodes] .- (9 - k)
                    #end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    r1 = r.a[nn2]
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
                                B[k*rowsOfB-0, l*pdim-0] = (F[3] * ∂h[3, (k-1)*numNodes+l] + F[9] * ∂h[1, (k-1)*numNodes+l]) / 2
                            end
                        end
                    else
                        error("loadVectorNonLinear: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    f1 .*= 0
                    for k in 1:numIntPoints
                        if mtype == :StVenantKirchhoff
                            I3 = [1 0 0; 0 1 0; 0 0 1]
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            E = (F2' * F2 - I3) / 2
                            S2 = 2μ * E + λ * (E[1, 1] + E[2, 2] + E[3, 3]) * I3
                        elseif mtype == :NeoHookeCompressible
                            I3 = [1 0 0; 0 1 0; 0 0 1]
                            ∂H0 = ∂H[k*9-(9-1):k*9, 1:pdim*numNodes]
                            F1 = ∂H0 * r1
                            F2 = reshape(F1, 3, 3)
                            J1 = det(F2)
                            C = F2' * F2
                            iC = inv(C)
                            #J1 = J1 > 0 ? J1 : 1e-14
                            S2 = μ * (I3 - iC) + λ * log(J1) * iC
                        else
                            error("loadVectorNonLinear: type '$(problem.type)' is not defined.")
                        end
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        f1 += B1' * S2[[1, 5, 9, 4, 6, 3]] * b * jacDet[k] * intWeights[k]
                    end
                    f[nn2] += f1
                end
            end
        end
    end
    return VectorField([], reshape(f, :,1), [0], [], 1, :f3D, problem)
end

function equivalentNodalForce(r)
    return equivalentNodalForce(r.model, r)
end

function nonFollowerLoadVector(problem, r, loads)
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
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(2order+1))
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
                    r1 = r.a[nn2]
                    F9 = F0.a[nn9]
                    for j in 1:numIntPoints
                        H9 = HH[j*9-(9-1):j*9, 1:9*numNodes]
                        F1 = H9 * F9
                        F1 = reshape(F1, 3, 3)
                        iF1 =inv(F1)
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
                        f1 += H1' * iF1 * f * Ja * intWeights[j]
                    end
                    fp[nn2] += f1
                end
            end
        end
    end
    return VectorField([], reshape(fp, :,1), [0], [], 1, :f3D, problem)
end

function nonFollowerLoadVector(r, loads)
    return nonFollowerLoadVector(r.model, r, loads)
end

"""
    FEM.applyDeformationBoundaryConditions!(problem, deformVec, supports)

Applies displacement boundary conditions `supports` on deformation vector `deformVec`.
Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Return: none

Types:
- `problem`: Problem
- `deformVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyDeformationBoundaryConditions!(problem, deformVec, supports; fact=1.0)
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    dof, st = size(deformVec.a)
    pdim = problem.pdim

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsX
                jj += 1
                if isa(ux, Function)
                    deformVec.a[j] += uux[jj] * fact
                else
                    deformVec.a[j] += ux * fact
                end
            end
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsY
                jj += 1
                if isa(uy, Function)
                    deformVec.a[j] += uuy[jj] * fact
                else
                    deformVec.a[j] += uy * fact
                end
            end
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsZ
                jj += 1
                if isa(uz, Function)
                    deformVec.a[j] += uuz[jj] * fact
                else
                    deformVec.a[j] += uz * fact
                end
            end
        end
    end
end

function applyDeformationBoundaryConditions!(deformVec, supports; fact=1.0)
    return applyDeformationBoundaryConditions!(deformVec.model, deformVec, supports, fact=fact)
end

"""
    FEM.suppressDeformationAtBoundaries!(problem, stiffMat, massMat, dampMat, loadVec, supports)

Applies displacement boundary conditions `supports` on a stiffness matrix
`stiffMat`, mass matrix `massMat`, damping matrix `dampMat` and load vector `loadVec`.
Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Return: none

Types:
- `problem`: Problem
- `stiffMat`: SparseMatrix 
- `massMat`: SparseMatrix 
- `dampMat`: SparseMatrix 
- `loadVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function suppressDeformationAtBoundaries!(problem, stiffMat, loadVec, supports)
#function applyDeformationBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports; fix=1)
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    dof, dof = size(stiffMat)
    pdim = problem.pdim

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim - 1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
                f0 = stiffMat[:, nodeTagsX] * uux * 0
            else
                f0 = stiffMat[:, nodeTagsX] * ux * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim - 2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
                f0 = stiffMat[:, nodeTagsY] * uuy * 0
            else
                f0 = stiffMat[:, nodeTagsY] * uy * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
                f0 = stiffMat[:, nodeTagsZ] * uuz * 0
            else
                f0 = stiffMat[:, nodeTagsZ] * uz * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
    end

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsX
                jj += 1
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(ux, Function)
                    loadVec.a[j] = uux[jj] * 0
                else
                    loadVec.a[j] = ux * 0
                end
            end
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsY
                jj += 1
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(uy, Function)
                    loadVec.a[j] = uuy[jj] * 0
                else
                    loadVec.a[j] = uy * 0
                end
            end
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
            end
            jj = 0
            for j ∈ nodeTagsZ
                jj += 1
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(uz, Function)
                    loadVec.a[j] = uuz[jj] * 0
                else
                    loadVec.a[j] = uz * 0
                end
            end
        end
    end

    dropzeros!(stiffMat)
    #dropzeros!(massMat)
    #dropzeros!(dampMat)
end

function suppressDeformationAtBoundaries!(stiffMat, loadVec, supports)
    return suppressDeformationAtBoundaries!(loadVec.model, stiffMat, loadVec, supports)
end

"""
    FEM.suppressDeformationAtBoundaries(problem, stiffMat, loadVec, supports)

Applies displacement boundary conditions `supports` on a stiffness matrix
`stiffMat`, mass matrix `massMat`, damping matrix `dampMat` and load vector `loadVec`.
Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Return: none

Types:
- `problem`: Problem
- `stiffMat`: SparseMatrix 
- `loadVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function suppressDeformationAtBoundaries(problem, stiffMat, loadVec, supports)
    K1 = copy(stiffMat)
    f1 = copy(loadVec)
    suppressDeformationAtBoundaries!(problem, K1, f1, supports)
    return K1, f1
end

function suppressDeformationAtBoundaries(stiffMat, loadVec, supports)
    return suppressDeformationAtBoundaries(loadVec.model, stiffMat, loadVec, supports)
end

function solveDeformation(problem, load, supp;
    followerLoad=false,
    loadSteps = 3,
    rampedLoad = true,
    rampedSupport = false,
    maxIteration = 10,
    saveSteps = false,
    saveIterations = false,
    plotConvergence = false,
    relativeError = 1e-5,
    initialDeformation=nodePositionVector(problem))

    r0 = initialDeformation
    err0 = abs(maximum(r0.a) - minimum(r0.a))
    f = loadVector(problem, load)
    r = []
    r1 = []
    push!(r, r0.a)
    r0 = copy(r0)
    e = []

    for j in range(1, loadSteps)
        fact = rampedLoad == true ? j / loadSteps : 1
        if rampedSupport == true
            applyDeformationBoundaryConditions!(problem, r0, supp, fact=1 / loadSteps)
        elseif rampedSupport == false && j == 1
            applyDeformationBoundaryConditions!(problem, r0, supp)
        end
        err = 1
        i = 0
        while err > relativeError && i < maxIteration
            i += 1
       
            Kl = tangentMatrixConstitutive(problem, r0)
            Knl = tangentMatrixInitialStress(problem, r0)
            if followerLoad == false
                f = nonFollowerLoadVector(problem, r0, load)
            end
            fnl = equivalentNodalForce(problem, r0)
            K1, f1 = suppressDeformationAtBoundaries(problem, Kl + Knl, fact * f - fnl, supp)
            q = solveDisplacement(K1, f1)
            r0 += q
            if saveIterations == true
                push!(r, r0.a)
                r0 = copy(r0)
            end
            err = maximum(abs.(q.a)) / err0
            if plotConvergence == true
                append!(e, err)
            end
        end
        if saveSteps == true
            push!(r, r0.a)
            r0 = copy(r0)
        end
    end
    if saveIterations == true || saveSteps == true
        n = length(r)
        r1 = zeros(length(r0.a), n)
        for i in 1:n
            r1[:, i] = r[i]
        end
    else
        n = length(r)
        r1 = zeros(length(r0.a), 1)
        r1[:,1] = r[:,n]
    end
    r1 = VectorField([], r1, 1:size(r1, 2), [], size(r1, 2), :u3D, problem)
    if plotConvergence == true
        return r1, e
    else
        return r1
    end
end

function showDeformationResults(problem, r, comp; name=comp, visible=false)
    r0 = nodePositionVector(problem)
    u = copy(r)
    for i in 1:size(r.a, 2)
        u.a[:, i] = r.a[:, i] - r0.a
    end
    return showDoFResults(problem, u, comp, name=name, visible=visible)
end

function showDeformationResults(r, comp; name=comp, visible=false)
    return showDeformationResults(r.model, r, comp, name=name, visible=visible)
end
