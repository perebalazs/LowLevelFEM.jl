function systemMatrixPoisson(problem::Problem, multiplier::Union{Number,ScalarField})
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

    if multiplier isa Number
        p = multiplier
    else
        p = nodesToElements(multiplier)
        pa = Dict(zip(p.numElem, p.A))
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        #kk = problem.material[ipg].k
        dim = problem.dim
        pdim = problem.pdim
        b = problem.thickness
        if problem.type == :Poisson
            rowsOfB = 3
        elseif problem.type == :Poisson2D
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
                ∂h = zeros(dim, numNodes * numIntPoints)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                @inbounds for k in 1:numIntPoints, l in 1:numNodes
                    val = h[(k-1)*numNodes + l]
                    for kk in 1:pdim
                        row = (k-1)*pdim + kk
                        col = (l-1)*pdim + kk
                        H[row, col] = val
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
                    pa0 = p isa ScalarField ? pa[elem][:,1] : p
                    for k in 1:numIntPoints
                        val = pa0 isa Vector ? dot(pa0, h[:, k]) : pa0
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += B1' * B1 * val * b * jacDet[k] * intWeights[k]
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
    return SystemMatrix(K, problem)
end

function systemMatrixPoisson2(problem::Problem, multiplier::Union{Number,ScalarField})
    gmsh.model.setCurrent(problem.name)

    # --- estimate nnz length for I/J/V (pdim=1) ---
    elemTypes0, elemTags0, elemNodeTags0 = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = 0
    @inbounds for i in 1:length(elemTags0)
        ne = length(elemTags0[i])
        ne == 0 && continue
        # numNodes for this element type group (assumes uniform within group)
        numNodes = div(length(elemNodeTags0[i]), ne)
        lengthOfIJV += (numNodes * numNodes) * ne
    end

    # --- preallocate triplets (fast path: fill by index) ---
    I = Vector{Int}(undef, lengthOfIJV)
    J = Vector{Int}(undef, lengthOfIJV)
    V = Vector{Float64}(undef, lengthOfIJV)
    pos = 1

    # --- multiplier handling (kept inside this single function) ---
    pa = nothing
    constp = 1.0
    use_pa = false
    if multiplier isa Number
        constp = Float64(multiplier)
    else
        p = nodesToElements(multiplier)  # elementwise field
        pa = Dict{Int, Matrix{Float64}}(zip(p.numElem, p.A))
        use_pa = true
    end

    @inline mult_at(elem, hcol) = begin
        if use_pa
            v = get(pa, elem, nothing)
            v === nothing ? 1.0 : dot(v, hcol)
        else
            constp
        end
    end

    dim = problem.dim
    b   = problem.thickness
    @assert problem.pdim == 1 "This optimized version assumes pdim=1."

    # --- main loop over physical groups ---
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for idm in 1:length(dimTags)
            edim, etag = dimTags[idm]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

            for itype in 1:length(elemTypes)
                et = elemTypes[itype]
                tags = elemTags[itype]
                ntags = elemNodeTags[itype]
                ne = length(tags)
                ne == 0 && continue

                elementName, dim_e, order1, numNodes::Int64, localNodeCoord, numPrimaryNodes =
                    gmsh.model.mesh.getElementProperties(et)

                # Integration points: keep your original choice; fix "2order" -> 2*order1
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order1 + 1))
                numIntPoints = length(intWeights)

                # Basis functions
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
                h = reshape(fun, numNodes, numIntPoints)  # (numNodes, numIntPoints)

                comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                ∇h = reshape(dfun, 3*numNodes, numIntPoints)  # (3*numNodes, numIntPoints)

                # Buffers (allocated once per element type)
                Ke   = zeros(Float64, numNodes, numNodes)
                dNdx = zeros(Float64, dim, numNodes)
                nn2  = Vector{Int}(undef, numNodes)

                # Loop over elements in this type group
                @inbounds for je in 1:ne
                    elem = tags[je]

                    # connectivity (pdim=1 => dof index is node index; keep original semantics)
                    base = (je-1)*numNodes
                    for a in 1:numNodes
                        nn2[a] = ntags[base + a]
                    end

                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)  # (3, 3*numIntPoints)

                    fill!(Ke, 0.0)

                    @inbounds for gp in 1:numIntPoints
                        # inv(J)' for this Gauss point (3x3)
                        jcol = 3*(gp-1)
                        @views Jgp = Jac[:, jcol+1:jcol+3]
                        invJt = inv(Jgp)'  # small (3x3) -> ok, but still general

                        # compute dNdx[:, a] for all nodes a at this gp
                        for a in 1:numNodes
                            ia = 3*(a-1)
                            # reference gradients from gmsh: ∂Na/∂(ξ,η,ζ) in slots 1..3
                            for d in 1:dim
                                s = 0.0
                                for e in 1:dim
                                    s += invJt[d, e] * ∇h[ia + e, gp]
                                end
                                dNdx[d, a] = s
                            end
                        end

                        # multiplier at this (elem,gp)
                        @views α = mult_at(elem, h[:, gp])

                        w = α * b * jacDet[gp] * intWeights[gp]

                        # Ke[a,b] += (∇Na ⋅ ∇Nb) * w
                        for a in 1:numNodes
                            for b2 in 1:numNodes
                                s = 0.0
                                for d in 1:dim
                                    s += dNdx[d, a] * dNdx[d, b2]
                                end
                                Ke[a, b2] += s * w
                            end
                        end
                    end

                    # write triplets (no append!, no Any, no Iidx/Jidx)
                    for a in 1:numNodes
                        Ia = nn2[a]
                        for b2 in 1:numNodes
                            I[pos] = Ia
                            J[pos] = nn2[b2]
                            V[pos] = Ke[a, b2]
                            pos += 1
                        end
                    end
                end
            end
        end
    end

    # In case lengthOfIJV was an over-estimate (shouldn't be, but safe):
    if pos <= length(I)
        resize!(I, pos-1); resize!(J, pos-1); resize!(V, pos-1)
    end

    dof = problem.non  # pdim=1
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end
