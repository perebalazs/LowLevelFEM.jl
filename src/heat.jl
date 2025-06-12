"""
    FEM.heatConductionMatrix(problem)

Solves the heat conduction matrix of the `problem`.

Return: `heatCondMat`

Types:
- `problem`: Problem
- `heatCondMat`: SparseMatrix
"""
function heatConductionMatrix(problem; elements=[])
    if problem.type == :AxiSymmetricHeatConduction
        return heatCondMatrixAXI(problem, elements=elements)
    else
        return heatCondMatrixSolid(problem, elements=elements)
    end
end

function heatCondMatrixSolid(problem; elements=[])
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
        kk = problem.material[ipg].k
        dim = problem.dim
        pdim = problem.pdim
        b = problem.thickness
        if problem.type == :HeatConduction
            rowsOfB = 3
        elseif problem.type == :PlaneHeatConduction
            rowsOfB = 2
        else
            error("heatCondMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
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
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
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
                    elseif dim == 3 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-2, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += B1' * B1 * kk * b * jacDet[k] * intWeights[k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
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

function heatCondMatrixAXI(problem; elements=[])
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
    ncoord2 = zeros(3 * problem.non)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        kk = problem.material[ipg].k
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 2 && problem.type == :AxiSymmetricHeatConduction
            rowsOfB = 2
        else
            error("heatCondMatrixAXI: dimension is $(problem.dim), problem type is $(problem.type).")
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
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                r = zeros(numIntPoints)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
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
                            #B[k*rowsOfB-2, l*pdim-1] = h[l, k] / r[k]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        #r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += 2π * B1' * B1 * kk * r[k] * jacDet[k] * intWeights[k]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
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

"""
    FEM.heatCapacityMatrix(problem; lumped=...)

Solves the heat capacity matrix of the `problem`. If `lumped` is true, solves lumped heat capacity matrix.

Return: `heatCapMat`

Types:
- `problem`: Problem
- `lumped`: Boolean
- `massMat`: SparseMatrix
"""
function heatCapacityMatrix(problem; elements=[], lumped=false)
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
    ncoord2 = zeros(3 * problem.non)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dim = problem.dim
        pdim = problem.pdim
        c = problem.material[ipg].c
        ρ = problem.material[ipg].ρ
        if problem.dim == 3 && problem.type == :HeatConduction
            rowsOfH = 3
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneHeatConduction
            rowsOfH = 2
            b = 1
        elseif problem.dim == 2 && problem.type == :AxiSymmetricHeatConduction
            rowsOfH = 2
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
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                nn2 = zeros(Int, pdim * numNodes)
                H = zeros(rowsOfH * numIntPoints, pdim * numNodes)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:pdim
                        H[k*pdim-(pdim-kk), l*pdim-(pdim-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                M1 = zeros(pdim * numNodes, pdim * numNodes)
                if problem.type != :AxiSymmetricHeatConduction
                    for j in 1:length(elemTags[i])
                        elem = elemTags[i][j]
                        for k in 1:numNodes
                            nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                        end
                        jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                        M1 .*= 0
                        for k in 1:numIntPoints
                            H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                            M1 += H1' * H1 * jacDet[k] * intWeights[k]
                        end
                        M1 *= c * ρ * b
                        for k in 1:pdim
                            nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                        end
                        append!(I, nn2[Iidx[:]])
                        append!(J, nn2[Jidx[:]])
                        append!(V, M1[:])
                    end
                elseif problem.type == :AxiSymmetricHeatConduction
                    for j in 1:length(elemTags[i])
                        elem = elemTags[i][j]
                        for k in 1:numNodes
                            nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                        end
                        jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                        M1 .*= 0
                        for k in 1:numIntPoints
                            r = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                            H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                            M1 += H1' * H1 * jacDet[k] * r * intWeights[k]
                        end
                        M1 *= 2π * c * ρ * b
                        for k in 1:pdim
                            nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                        end
                        append!(I, nn2[Iidx[:]])
                        append!(J, nn2[Jidx[:]])
                        append!(V, M1[:])
                    end
                end
                push!(nn, nnet)
            end
        end
    end
    dof = problem.pdim * problem.non
    M = sparse(I, J, V, dof, dof)
    if lumped == true
        M = spdiagm(vec(sum(M, dims=2))) # lumped mass matrix
    end
    dropzeros!(M)
    return M
end

"""
    FEM.latentHeatMatrix(problem, u, v, T0)

Solves the latent heat matrix of the `problem`. With this matrix the generated heat due to deformations
(given with displacement field `u` and velocity field `v`) can be solved. `T0` is the current temperature
field which is given in absolute temperature scale (Kelvin).

Return: `latHeatMat`

Types:
- `problem`: Problem
- `u`: Vector{Float64}
- `v`: Vector{Float64}
- `T0`: Vector{Float64}
- `latHeatMat`: SparseMatrix
"""
function latentHeatMatrix(problem, u, v, T0)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * problem.dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
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
        α = problem.material[ipg].α
        κ = E / 3 / (1 - 2ν)
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 3 && problem.type == :HeatConduction
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneHeatConduction
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness
        else
            error("latentHeatMatrix: dimension is $(problem.dim), problem type is $(problem.type).")
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
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                M1 = zeros(pdim * numNodes, pdim * numNodes)
                ∇H = zeros(1 * numIntPoints, dim * numNodes)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, dim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn1 = zeros(Int, pdim * numNodes)
                nn2 = zeros(Int, dim * numNodes)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:pdim
                        H[k*pdim-(pdim-kk), l*pdim-(pdim-kk)] = h[(k-1)*numNodes+l]
                    end
                end
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-0, l*dim-0] = B[k*rowsOfB-2, l*dim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*dim-1] = B[k*rowsOfB-1, l*dim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-5, l*dim-2] = B[k*rowsOfB-2, l*dim-1] = B[k*rowsOfB-0, l*dim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*dim-1] = B[k*rowsOfB-2, l*dim-2] = B[k*rowsOfB-1, l*dim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*dim-0] = B[k*rowsOfB-1, l*dim-1] = B[k*rowsOfB-0, l*dim-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    ∇H .*= 0
                    if dim == 2
                        for k in 1:numIntPoints, l in 1:numNodes
                            ∇H[k*1-0, l*dim-1] = ∂h[1, (k-1)*numNodes+l]
                            ∇H[k*1-0, l*dim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            ∇H[k*1-0, l*dim-2] = ∂h[1, (k-1)*numNodes+l]
                            ∇H[k*1-0, l*dim-1] = ∂h[2, (k-1)*numNodes+l]
                            ∇H[k*1-0, l*dim-0] = ∂h[3, (k-1)*numNodes+l]
                        end
                    end
                    nn1 = nnet[j, 1:numNodes]
                    for k in 1:dim
                        nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
                    end
                    M1 .*= 0
                    K1 .*= 0
                    for k in 1:numIntPoints
                        q1 = u[nn2]
                        dq1 = v[nn2]
                        T01 = T0[nn1]
                        H1 = H[k*pdim-(pdim-1):k*pdim, 1:pdim*numNodes]
                        ∇H1 = ∇H[k, 1:dim*numNodes]'
                        M1 += H1' * H1 * (∇H1 * dq1) * jacDet[k] * intWeights[k]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:dim*numNodes]
                        K1 += H1' * H1 * (q1' * B1' * D * B1 * dq1) / (H1 * T01)[1] * b * jacDet[k] * intWeights[k]
                    end
                    M1 *= κ * α * b
                    KM1 = K1 - M1
                    append!(I, nn1[Iidx[:]])
                    append!(J, nn1[Jidx[:]])
                    append!(V, KM1[:])
                end
            end
        end
    end
    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return K
end

"""
    FEM.heatConvectionMatrix(problem, heatConvection)

Solves the heat convection matrix of the `problem`. `heatConvection` 
is a vector of heat convection boundary condicions defined in function
`FEM.heatConduction`. With the heat convection vector (see the `heatConvectionVector` function)
`heatConvVec`, temperature field vector `T` in hand the heat flux vector `qCV` arising from the
heat convection boundary condition can be solved. `qCV = heatConvMat * T - heatConvVec`

Return: `heatConvMat`

Types:
- `problem`: Problem
- `heatConvection`: Vector{Tuple{String, Float64, Float64, Float64}}
- `heatConvMat`: SparseMatrix
"""
function heatConvectionMatrix(problem, heatConvection)
    if !isa(heatConvection, Vector)
        error("heatConvection: heat convections are not arranged in a vector. Put them in [...]")
    end
    return elasticSupportMatrix(problem, heatConvection)
end

"""
    FEM.heatFluxVector(problem, heatFlux)

Solves a heat flux or heat source vector of `problem`. `heatFlux` is a tuple of name of physical group 
`name`, heat flux `qn` normal to the surface of the body. The outward direction is positive.
It can solve heat flux (or heat source) depending on the problem.
- In case of 2D problems and Point physical group means concentrated heat flux.
- In case of 2D problems and Line physical group means surface heat flux.
- In case of 2D problems and Surface physical group means body heat source.
- In case of 3D problems and Point physical group means concentrated heat flux.
- In case of 3D problems and Line physical group means edge heat source.
- In case of 3D problems and Surface physical group means surface heat flux.
- In case of 3D problems and Volume physical group means body heat source.

Return: `heatFluxVec`

Types:
- `problem`: Problem
- `heatFlux`: Vector{Tuple{String, Float64, Float64, Float64}}
- `heatFluxVec`: Vector
"""
function heatFluxVector(problem, loads)
    if !isa(loads, Vector)
        error("heatFluxVector: heat fluxes are not arranged in a vector. Put them in [...]")
    end
    return loadVector(problem, loads)
end

"""
    FEM.heatSourceVector(problem, heatSource)

Solves a heat flux or heat source vector of `problem`. `heatSource` is a tuple of name of physical group 
`name`, heat flux `qn` normal to the surface of the body. The outward direction is positive.
It can solve heat flux (or heat source) depending on the problem.
- In case of 2D problems and Point physical group means concentrated heat flux.
- In case of 2D problems and Line physical group means surface heat flux.
- In case of 2D problems and Surface physical group means body heat source.
- In case of 3D problems and Point physical group means concentrated heat flux.
- In case of 3D problems and Line physical group means edge heat source.
- In case of 3D problems and Surface physical group means surface heat flux.
- In case of 3D problems and Volume physical group means body heat source.
Same as the `heatFluxVector` function.

Return: `heatSourceVec`

Types:
- `problem`: Problem
- `heatSource`: Vector{Tuple{String, Float64, Float64, Float64}}
- `heatSourceVec`: Vector
"""
function heatSourceVector(problem, loads)
    if !isa(loads, Vector)
        error("heatSource: heat sources are not arranged in a vector. Put them in [...]")
    end
    return loadVector(problem, loads)
end

"""
    FEM.heatConvectionVector(problem, heatConvection)

Solves a heat convection vector of `problem`. `heatConvection` 
is a vector of heat convection boundary condicions defined in function
`FEM.heatConduction`. With the heat convection matrix (see the `heatConvectionMatrix` function)
`heatConvMat`, temperature field vector `T` in hand the heat flux vector `qCV` arising from the
heat convection boundary condition can be solved. `qCV = heatConvMat * T - heatConvVec`

Return: `heatConvVec`

Types:
- `problem`: Problem
- `heatConvection`: Vector{Tuple{String, Float64, Float64, Float64}}
- `heatConvVec`: Vector
"""
function heatConvectionVector(problem, heatConvection)
    if !isa(heatConvection, Vector)
        error("heatConvectionVector: heat convections are not arranged in a vector. Put them in [...]")
    end
    return loadVector(problem, heatConvection)
end

"""
    FEM.thermalLoadVector(problem, T; T₀=...)

Solves the thermal load vector from a temperature field `T` for problem `problem`.
`T₀` is the initial temperature field.

Return: `thermLoadVec`

Types:
- `problem`: Problem
- `T`: Vector{Float64}
- `T₀`: Vector{Float64}
- `thermLoadVec`: Vector{Float64}
"""
function thermalLoadVector(problem, T; T₀=1im)
    if problem.type == :AxiSymmetric
        return thermalLoadVectorAXI(problem, T, T₀=T₀)
    else
        return thermalLoadVectorSolid(problem, T, T₀=T₀)
    end
end

function thermalLoadVectorSolid(problem, T; T₀=1im)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    dim = problem.dim
    pdim = problem.pdim
    dof = problem.non * pdim
    fT = zeros(dof)
    if T₀ == 1im
        T₀ = zeros(problem.non)
    end
    if size(T) != size(T₀) || size(T,1) != problem.non
        error("thermalLoadVectorSolid: size of T [$(size(T))] != size of T₀ [$(size(T₀))], non=$(problem.non)")
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        α = problem.material[ipg].α
        if problem.dim == 3 && problem.type == :Solid
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            rowsOfB = 6
            b = α
            E0 = [1,1,1,0,0,0]
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            rowsOfB = 3
            b = problem.thickness * α
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            rowsOfB = 3
            b = α
            E0 = [1,1,0]
        else
            error("thermalLoadVectorSolid: dimension is $(problem.dim), problem type is $(problem.type).")
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
                pdimT = 1
                H = zeros(pdimT * numIntPoints, pdimT * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdimT
                            H[j*pdimT-(pdimT-l), k*pdimT-(pdimT-l)] = h[k, j]
                        end
                    end
                end
                invJac = zeros(3, 3numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn1 = zeros(Int, pdimT * numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:pdimT
                        nn1[k:pdimT:pdimT*numNodes] = pdimT * nnet[j, 1:numNodes] .- (pdimT - k)
                    end
                    f1 .*= 0
                    for k in 1:numIntPoints
                        H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        f1 += B1' * D * E0 * H1 * (T[nn1] - T₀[nn1]) * b * jacDet[k] * intWeights[k]
                    end
                    fT[nn2] += f1
                end
            end
        end
    end
    return fT
end

function thermalLoadVectorAXI(problem, T; T₀=1im)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    dim = problem.dim
    pdim = problem.pdim
    dof = problem.non * pdim
    fT = zeros(dof)
    if T₀ == 1im
        T₀ = zeros(problem.non)
    end
    if size(T) != size(T₀) || size(T,1) != problem.non
        error("thermalLoadVectorSolid: size of T [$(size(T))] != size of T₀ [$(size(T₀))], non=$(problem.non)")
    end
    ncoord2 = zeros(3 * problem.non)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        α = problem.material[ipg].α
        if problem.dim == 2 && problem.type == :AxiSymmetric
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                ν 1-ν ν 0;
                ν ν 1-ν 0;
                0 0 0 (1-2ν)/2]

            rowsOfB = 4
            b = α
            E0 = [1,1,1,0]
        else
            error("thermalLoadVectorAxiSymmetric: dimension is $(problem.dim), problem type is $(problem.type).")
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
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, :, numIntPoints)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                pdimT = 1
                H = zeros(pdimT * numIntPoints, pdimT * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdimT
                            H[j*pdimT-(pdimT-l), k*pdimT-(pdimT-l)] = h[k, j]
                        end
                    end
                end
                invJac = zeros(3, 3numIntPoints)
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                nn1 = zeros(Int, pdimT * numNodes)
                r = zeros(numIntPoints)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 4
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-3, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = B[k*rowsOfB-0, l*pdim-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-2, l*pdim-1] = h[l, k] / r[k]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    f1 .*= 0
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:pdimT
                        nn1[k:pdimT:pdimT*numNodes] = pdimT * nnet[j, 1:numNodes] .- (pdimT - k)
                    end
                    for k in 1:numIntPoints
                        H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        f1 += 2π * B1' * D * E0  * H1 * (T[nn1] - T₀[nn1]) * b * r[k] * jacDet[k] * intWeights[k]
                    end
                    fT[nn2] += f1
                end
            end
        end
    end
    return fT
end

"""
    FEM.applyHeatConvection!(problem, heatCondMat, heatFluxVec, heatConv)

Applies heat convectiom boundary conditions `heatConv` on a heat conduction matrix
`heatCondMat` and heat flux vector `heatFluxVec`. Mesh details are in `problem`. `heatConv`
is a tuple of `name` of physical group and prescribed heat transfer coefficient `h`
and ambient temperature `Tₐ`. The ambient temperature can be either a constant or a 
function of x, y and z coordinates.

Return: none

Types:
- `problem`: Problem
- `heatCondMat`: SparseMatrix 
- `heatFluxVec`: Vector 
- `heatConv`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyHeatConvection!(problem, heatCondMat, heatFluxVec, heatConv)
    if !isa(heatConv, Vector)
        error("applyHeatConvection!: heat convections are not arranged in a vector. Put them in [...]")
    end
    hf0 = heatConvectionVector(problem, heatConv)
    C0 = heatConvectionMatrix(problem, heatConv)
    heatCondMat .+= C0
    heatFluxVec .+= hf0
    return nothing
end

"""
    FEM.solveTemperature(K, q)

Solves the equation K*T=q for the temperature vector `T`. `K` is the heat conduction matrix,
`q` is the heat flux vector.

Return: `T`

Types:
- `K`: SparseMatrix 
- `q`: Vector{Float64} 
- `T`: Vector{Float64}
"""
function solveTemperature(K, q)
    return K \ q
end

"""
    FEM.solveTemperature(problem, flux, temp)

Solves the temperature vector `T` of `problem` with given heat flux `flux` and
temperature `temp`.

Return: `T`

Types:
- `problem`: Problem 
- `flux`: Vector{Tuple} 
- `temp`: Vector{Tuple}
- `T`: Vector{Float64}
"""
function solveTemperature(problem, flux, temp)
    K = heatConductionMatrix(problem)
    q = heatFluxVector(problem, flux)
    #applyHeatConvection!(problem, K, q, temp)
    applyBoundaryConditions(problem, K, q, temp)
    return K \ q
end

"""
    FEM.solveTemperature(problem, flux, temp, heatconv)

Solves the temperature vector `T` of `problem` with given heat flux `flux`,
temperature `temp` and heat convection `heatconv`.

Return: `T`

Types:
- `problem`: Problem 
- `flux`: Vector{Tuple} 
- `temp`: Vector{Tuple}
- `heatconv`: Vector{Tuple}
- `T`: Vector{Float64}
"""
function solveTemperature(problem, flux, temp, heatconv)
    K = heatConductionMatrix(problem)
    q = heatFluxVector(problem, flux)
    applyHeatConvection!(problem, K, q, heatconv)
    applyBoundaryConditions(problem, K, q, temp)
    return K \ q
end

"""
    FEM.solveHeatFlux(problem, T; DoFResults=false)

Solves the heat flux field `q` from temperature vector `T`. heat flux is given
per elements, so it usually contains jumps at the boundary of elements. Details
of mesh is available in `problem`. If `DoFResults` is true, `q` is a matrix with
nodal results. In this case `showDoFResults` can be used to show the results
(otherwise `showHeatFluxResults` or `showElementResults`).

Return: `q`

Types:
- `problem`: Problem
- `T`: Vector{Float64}
- `q`: VectorField or Matrix{Float}
"""
function solveHeatFlux(problem, T; DoFResults=false)
    gmsh.model.setCurrent(problem.name)

    type = :q
    nsteps = size(T, 2)
    σ = []
    numElem = Int[]
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    ncoord2 = zeros(3 * non)
    if DoFResults == true
        q1 = zeros(non * dim, nsteps)
        pcs = zeros(Int64, non * dim)
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        kT = -problem.material[ipg].k
        dim = 0
        if problem.dim == 3 && problem.type == :HeatConduction
            dim = 3
            rowsOfB = 3
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneHeatConduction
            dim = 2
            rowsOfB = 2
            b = 1
        elseif problem.dim == 2 && problem.type == :AxiSymmetricHeatConduction
            dim = 2
            rowsOfB = 2
            b = 1
        else
            error("solveHeatFlux: dimension is $(problem.dim), problem type is $(problem.type).")
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
                s0 = zeros(rowsOfB * numNodes)
                nodeCoord = zeros(numNodes * 3)
                for k in 1:dim, j = 1:numNodes
                    nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
                end
                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
                ∇h = reshape(dfun, :, numNodes)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
                h = reshape(fun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                pdimT = 1
                H = zeros(pdimT * numNodes, pdimT * numNodes)
                for j in 1:numNodes
                    for k in 1:numNodes
                        for l in 1:pdimT
                            H[j*pdimT-(pdimT-l), k*pdimT-(pdimT-l)] = h[k, j]
                        end
                    end
                end
                invJac = zeros(3, 3numNodes)
                ∂h = zeros(3, numNodes * numNodes)
                B = zeros(rowsOfB * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                r = zeros(numNodes)
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)
                    for k in 1:numNodes
                        invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                        r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    end
                    ∂h .*= 0
                    for k in 1:numNodes, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
                    end
                    B .*= 0
                    if pdim == 1 && rowsOfB == 2
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif pdim == 1 && rowsOfB == 3
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*rowsOfB-2, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("solveHeatFlux: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    s = zeros(3numNodes, nsteps) # vectors have three elements
                    for k in 1:numNodes
                        if rowsOfB == 3 && pdim == 1
                            B1 = B[k*rowsOfB-2:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = B1 * T[nn2, kk] * kT
                                s[(k-1)*3+1:k*3, kk] = [s0[1], s0[2], s0[3]]
                                if DoFResults == true
                                    #q1[dim*nnet[j, k]-2, kk] += s0[1]
                                    #q1[dim*nnet[j, k]-1, kk] += s0[2]
                                    #q1[dim*nnet[j, k]-0, kk] += s0[3]
                                    q1[dim*nnet[j, k]-2:dim*nnet[j,k], kk] .+= s0
                                end
                            end
                        elseif rowsOfB == 2 && pdim == 1
                            B1 = B[k*rowsOfB-1:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = B1 * T[nn2, kk] * kT
                                s[(k-1)*3+1:k*3, kk] = [s0[1], s0[2], 0]
                                if DoFResults == true
                                    #q1[dim*nnet[j, k]-1, kk] += s0[1]
                                    #q1[dim*nnet[j, k]-0, kk] += s0[2]
                                    q1[dim*nnet[j, k]-1:dim*nnet[j,k], kk] .+= s0
                                end
                            end
                        else
                            error("solveHeatFlux: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
                        end
                    end
                    if DoFResults == true
                        pcs[nnet[j,1:numNodes]] .+= 1
                    end
                    if DoFResults == false
                        push!(σ, s)
                    end
                end
            end
        end
    end
    if DoFResults == true
        for k in 1:dim
            for l in 1:non
                q1[k + dim * l - dim, :] ./= pcs[l]
            end
        end
    end
    if DoFResults == true
        return q1
    else
        sigma = VectorField(σ, numElem, nsteps, type)
        return sigma
    end
end

"""
    FEM.initialTemperature(problem, name; T=...)

Sets the temperature value `T` at nodes belonging to physical group `name`.
Returns the `T0` initial nodal temperature vector.

Return: T0

Types:
- `problem`: Problem
- `name`: String 
- `T`: Float64 
- `T0`: Vector{Float64}
"""
function initialTemperature(problem, name; T=1im)
    dim = problem.pdim
    T0 = zeros(problem.non * problem.pdim)
    phg = getTagForPhysicalName(name)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
    if T != 1im
        for i in 1:length(nodeTags)
            T0[nodeTags[i]*dim-(dim-1)] = T
        end
    end
    return T0
end

"""
    FEM.initialTemperature!(problem, name, T0; T=...)

Changes the tempetature value to `T` at nodes belonging to physical group `name`.
Original values are in temperature vector `T0`.

Return: none

Types:
- `problem`: Problem
- `name`: String 
- `T0`: Vector{Float64}
- `T`: Float64 
"""
function initialTemperature!(problem, name, T0; T=1im)
    dim = problem.pdim
    phg = getTagForPhysicalName(name)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
    if T != 1im
        for i in 1:length(nodeTags)
            T0[nodeTags[i]*dim-(dim-1)] = T
        end
    end
end

"""
    FEM.FDM(K, C, q, T0, tₘₐₓ, Δt; ϑ=...)

Solves a transient heat conduction problem using Finite Difference Method (FDM).
Introducing a `ϑ` parameter, special cases can be used as the Forward Euler (explicit, ϑ=0),
Backward Euler (implicit, ϑ=1), Crank-Nicolson (ϑ=0.5) and intermediate cases (0<ϑ<1).
(This method is known as ϑ-method. See [^5].)
`K` is the heat conduction matrix, `C` is the heat capacity matrix,
`q` is the heat flux vector, `T0` is the initial temperature, `tₘₐₓ` is the upper 
bound of the time intervall (lower bound is zero)
and `Δt` is the time step size. Returns the nodal temperature vectors in each time 
step arranged in the columns of the matrix `T`
and a vector `t` of the time instants used.

The critical (largest allowed) time step is `Δtₘₐₓ = 2 / ((1-2ϑ)*λₘₐₓ)`
where `λₘₐₓ` is the largest eigenvalue of (**K**+λ**C**)**θ**=**0** 
eigenvalue problem and `ϑ` is the parameter of the ϑ-method. Default value of `ϑ`
is 1/2.

[^5]: Bathe, K. J.: Finite element procedures, Wiley, 1983, <https://doi.org/10.1002/nag.1610070412>

Return: `T`, `t`

Types:
- `K`: SparseMatrix
- `C`: SparseMatrix
- `q`: Vector{Float64}
- `T0`: Vector{Float64}
- `tₘₐₓ`: Float64
- `Δt`: Float64 
- `T`: Matrix{Float64}
- `t`: Vector{Float64}
"""
function FDM(K, C, q, T0, tₘₐₓ, Δt; ϑ=0.5)
    nnz = nonzeros(C)
    dof, dof = size(C)
    nsteps = ceil(Int64, tₘₐₓ / Δt)
    T = zeros(dof, nsteps)
    t = zeros(nsteps)
    T[:, 1] = T0
    t[1] = 0
    if ϑ == 0 && nnz == dof
        invC = spdiagm(1 ./ diag(C))
        for i in 2:nsteps
            T1 = T0 - (1 - ϑ) * Δt * invC * K * T0 + Δt * invC * q
            T[:, i] = T1
            t[i] = t[i-1] + Δt
            T0 = T1
        end
    else
        A = C + ϑ * Δt * K
        b = C - (1 - ϑ) * Δt * K
        AA = lu(A)
        for i in 2:nsteps
            bb = b * T0 + Δt * q
            T1 = AA \ bb
            T[:, i] = T1
            t[i] = t[i-1] + Δt
            T0 = T1
        end
    end

    return T, t
end

"""
    FEM.FDMaccuracyAnalysis(λₘᵢₙ, λₘₐₓ, Δt; type, n=100, ϑ=...)

Gives a functions (graphs) for accuracy analysis of the ϑ-method[^5]. 
`λₘᵢₙ` and `λₘₐₓ` are the smallest and largest eigenvalues of the
**Kθ**=λ**Cθ** eigenvalue problem, `Δt` is the time step size. `type` is the
"SR" spectral radius.
`n` is the number of points in the graph. For the meaning of `ϑ` see [^5].
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

- "SR": spectral radius as a function of λ⋅Δt

Return: `xy`

Types:
- `λₘᵢₙ`: Float64
- `λₘₐₓ`: Float64
- `Δt`: Float64 
- `n`: Int64
- `ϑ`: Float64
- `xy`: Tuple{Vector{Float64},Vector{Float64}}
"""

function FDMaccuracyAnalysis(λₘᵢₙ, λₘₐₓ, Δt; type=:SR, n=100, ϑ=0.5)
    x = zeros(n)
    y = similar(x)
    Λ = range(λₘᵢₙ, length=n, stop=λₘₐₓ)
    for i ∈ 1:n
        A1 = 1 - (1 - ϑ) * Δt * Λ[i]
        A2 = 1 + ϑ * Δt * Λ[i]

        A = A1 \ A2

        ρ = abs(A)
        if type == :SR
            x[i] = Λ[i] * Δt
            y[i] = ρ
        else
            str1 = "FDMaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\nSR: spectral radius\n"
            error(str1*str2*str3)
        end
    end
    return x, y
end

