export stiffnessMatrix, nonLinearStiffnessMatrix, massMatrix, dampingMatrix
export elasticSupportMatrix
export loadVector
export applyBoundaryConditions, applyBoundaryConditions!, applyElasticSupport!
export solveDisplacement, solveStrain, solveStress, solveAxialForce
export solveEigenModes, solveBucklingModes
export solveModalAnalysis, solveBuckling
export initialDisplacement, initialDisplacement!, initialVelocity, initialVelocity!
export nodalForce!, nodalAcceleration!
export largestPeriodTime, smallestPeriodTime, largestEigenValue, smallestEigenValue
export CMD, HHT, CDMaccuracyAnalysis, HHTaccuracyAnalysis

# linear.jl tetején
mutable struct SolidWS
    dH::Matrix{Float64}
    B::Matrix{Float64}
    nnet::Vector{Int}
    K1::Matrix{Float64}
    nn2::Vector{Int}
end

function make_ws(dim::Int, numNodes::Int, numIntPoints::Int, rowsOfB::Int, dof_elem::Int)
    SolidWS(
        zeros(dim, numNodes*numIntPoints),
        zeros(rowsOfB*numIntPoints, dof_elem),
        zeros(Int, numNodes),
        zeros(dof_elem, dof_elem),
        zeros(Int, dof_elem)
    )
end

mutable struct AxisWS
    dH::Matrix{Float64}   # dim × (numNodes*numIntPoints)
    B::Matrix{Float64}    # (rowsOfB*numIntPoints) × dof_elem
    nnet::Vector{Int}     # numNodes
    K1::Matrix{Float64}   # dof_elem × dof_elem
    nn2::Vector{Int}      # dof_elem
    r::Vector{Float64}    # numIntPoints
end

function make_ws_axi(dim::Int, numNodes::Int, numIntPoints::Int, rowsOfB::Int, dof_elem::Int)
    AxisWS(
        zeros(dim, numNodes*numIntPoints),
        zeros(rowsOfB*numIntPoints, dof_elem),
        zeros(Int, numNodes),
        zeros(dof_elem, dof_elem),
        zeros(Int, dof_elem),
        zeros(numIntPoints)
    )
end

"""
    stiffnessMatrix(problem)

Solves the stiffness matrix of the `problem`.

Returns: `stiffMat`

Types:
- `problem`: Problem
- `stiffMat`: SystemMatrix

# Examples

```julia
K = stiffnessMatrix(problem)
```
"""
function stiffnessMatrix(problem; elements=[])
    if problem.type == :AxiSymmetric && Threads.nthreads() == 1
        return stiffnessMatrixAXI(problem, elements=elements)
    elseif problem.type == :AxiSymmetric && Threads.nthreads() > 1
        return stiffnessMatrixAXIParallel(problem, elements=elements)
        #return stiffnessMatrixAXI(problem, elements=elements)
    elseif problem.type == :Truss
        return stiffnessMatrixTruss(problem, elements=elements)
    elseif (problem.type != :AxiSymmetric || problem.type != :Truss) && Threads.nthreads() == 1
        return stiffnessMatrixSolid(problem, elements=elements)
    elseif (problem.type != :AxiSymmetric || problem.type != :Truss) && Threads.nthreads() > 1
        return stiffnessMatrixSolidParallel(problem, elements=elements)
        #return stiffnessMatrixSolid(problem, elements=elements)
    end
end

@inline function inv3x3_fast(J::SMatrix{3,3,Float64})
    # kofaktorok + det, majd 1/det * adj(J)
    J11, J12, J13 = J[1,1], J[1,2], J[1,3]
    J21, J22, J23 = J[2,1], J[2,2], J[2,3]
    J31, J32, J33 = J[3,1], J[3,2], J[3,3]

    c11 =  (J22*J33 - J23*J32)
    c12 = -(J21*J33 - J23*J31)
    c13 =  (J21*J32 - J22*J31)

    c21 = -(J12*J33 - J13*J32)
    c22 =  (J11*J33 - J13*J31)
    c23 = -(J11*J32 - J12*J31)

    c31 =  (J12*J23 - J13*J22)
    c32 = -(J11*J23 - J13*J21)
    c33 =  (J11*J22 - J12*J21)

    detJ = J11*c11 + J12*c12 + J13*c13
    invdet = 1.0 / detJ

    # inv(J) = adj(J)/det, majd transzponálás, ha kell
    return @SMatrix [c11*invdet  c21*invdet  c31*invdet;
                     c12*invdet  c22*invdet  c32*invdet;
                     c13*invdet  c23*invdet  c33*invdet]
end

@inline function inv2x2_fast(A11::Float64, A12::Float64, A21::Float64, A22::Float64)
    detA = A11*A22 - A12*A21
    invdet = 1.0 / detA
    # [A11 A12; A21 A22]^{-1} = (1/det)[ A22 -A12; -A21 A11]
    return @SMatrix [ A22*invdet  -A12*invdet;
                     -A21*invdet   A11*invdet]
end

function stiffnessMatrixSolid(problem; elements=[])
    gmsh.model.setCurrent(problem.name)

    # --- előszámítás ---
    elemTypes_all, elemTags_all, elemNodeTags_all =
        gmsh.model.mesh.getElements(problem.dim, -1)

    lengthOfIJV = sum([
        (div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.dim)^2 *
        length(elemTags_all[i]) for i in 1:length(elemTags_all)
    ])

    nn = Vector{Any}()
    I  = Int[]
    J  = Int[]
    V  = Float64[]
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)

    # --- anyagcsoportonként ---
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim, pdim = problem.dim, problem.pdim

        # anyagmátrix D
        if problem.dim == 3 && problem.type == :Solid
            D = E / ((1 + ν) * (1 - 2ν)) * @SMatrix [
                1-ν  ν    ν    0            0            0;
                ν    1-ν  ν    0            0            0;
                ν    ν    1-ν  0            0            0;
                0    0    0    (1-2ν)/2     0            0;
                0    0    0    0            (1-2ν)/2     0;
                0    0    0    0            0            (1-2ν)/2
            ]
            rowsOfB, b = 6, 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * @SMatrix [
                1   ν   0;
                ν   1   0;
                0   0  (1-ν)/2
            ]
            rowsOfB, b = 3, problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            D = E / ((1 + ν) * (1 - 2ν)) * @SMatrix [
                1-ν   ν    0;
                ν    1-ν   0;
                0     0   (1-2ν)/2
            ]
            rowsOfB, b = 3, 1
        else
            error("stiffnessMatrixSolid: dimension=$(problem.dim), type=$(problem.type)")
        end

        # --- fizikai csoport entitásai ---
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # --- érintett elemtípusok ---
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x -> push!(types_in_group, x), t)
        end

        # --- típusonként ---
        for et in types_in_group
            elementName, edim, order, numNodes::Int, localNodeCoord, numPrimaryNodes =
                gmsh.model.mesh.getElementProperties(et)

            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order+1))
            numIntPoints = length(intWeights)

            comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
            ∇h_ref = reshape(dfun, :, numIntPoints)

            elemTags_global, elemNodeTags_global = gmsh.model.mesh.getElementsByType(et)

            jac_all, det_all, coords_all = gmsh.model.mesh.getJacobians(et, intPoints)
            nip = length(det_all) ÷ length(elemTags_global)

            # tag → globális index map
            idxmap = Dict{Int,Int}()
            sizehint!(idxmap, length(elemTags_global))
            @inbounds for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            # --- csoportonként ---
            for (edim, etag) in dimTags
                if edim != problem.dim; continue; end
                elemTags_local, elemNodeTags_local = gmsh.model.mesh.getElementsByType(et, etag)
                if isempty(elemTags_local); continue; end

                nloc = length(elemTags_local)
                nnet = zeros(Int, nloc, numNodes)

                ∂h = zeros(dim, numNodes * numIntPoints)
                B  = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim*numNodes, pdim*numNodes)
                nn2 = zeros(Int, pdim*numNodes)

                Iidx = [l for k in 1:pdim*numNodes, l in 1:pdim*numNodes]
                Jidx = [k for k in 1:pdim*numNodes, l in 1:pdim*numNodes]

                @inbounds for j in 1:nloc
                    # node-tag-ek
                    for k in 1:numNodes
                        nnet[j,k] = elemNodeTags_local[(j-1)*numNodes+k]
                    end

                    gidx = idxmap[elemTags_local[j]]
                    jac_slice = view(jac_all, gidx*9*nip + 1 : (gidx+1)*9*nip)
                    det_slice = view(det_all, gidx*nip + 1   : (gidx+1)*nip)

                    fill!(∂h, 0.0)

                    @inbounds for k in 1:numIntPoints
                        Jk = SMatrix{3,3}(jac_slice[(k-1)*9+1 : k*9])
                        invJT = inv(Jk)'

                        for l in 1:numNodes
                            ∂h[1:dim, (k-1)*numNodes+l] =
                                invJT[1:dim,1:dim] * ∇h_ref[l*3-2 : l*3-(3-dim), k]
                        end
                    end

                    fill!(B, 0.0)
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[1,(k-1)*numNodes+l]
                            B[k*rowsOfB-2, l*pdim-1] = ∂h[1,(k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-1] = ∂h[2,(k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[2,(k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-5, l*pdim-2] = ∂h[1,(k-1)*numNodes+l]
                            B[k*rowsOfB-2, l*pdim-1] = ∂h[1,(k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-0] = ∂h[1,(k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*pdim-1] = ∂h[2,(k-1)*numNodes+l]
                            B[k*rowsOfB-2, l*pdim-2] = ∂h[2,(k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = ∂h[2,(k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*pdim-0] = ∂h[3,(k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-1] = ∂h[3,(k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-2] = ∂h[3,(k-1)*numNodes+l]
                        end
                    end

                    fill!(K1, 0.0)
                    @inbounds for k in 1:numIntPoints
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 .+= B1' * D * B1 * b * det_slice[k] * intWeights[k]
                    end

                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] =
                            pdim * nnet[j,1:numNodes] .- (pdim - k)
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

function stiffnessMatrixSolidParallel(problem; elements=[])
    # --- BLAS szálak mentése + átállítása ---
    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    gmsh.model.setCurrent(problem.name)

    # Előszámítás I,J,V méretére
    elemTypes_all, elemTags_all, elemNodeTags_all =
        gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([
        (div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.dim)^2 *
        length(elemTags_all[i]) for i in 1:length(elemTags_all)
    ])

    I  = Int[]
    J  = Int[]
    V  = Float64[]
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)

    # Anyagcsoportonként
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim, pdim = problem.dim, problem.pdim

        # anyagmátrix
        if problem.dim == 3 && problem.type == :Solid
            D = E / ((1 + ν) * (1 - 2ν)) * @SMatrix [
                1-ν  ν    ν    0            0            0;
                ν    1-ν  ν    0            0            0;
                ν    ν    1-ν  0            0            0;
                0    0    0    (1-2ν)/2     0            0;
                0    0    0    0            (1-2ν)/2     0;
                0    0    0    0            0            (1-2ν)/2
            ]
            rowsOfB, b = 6, 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * @SMatrix [
                1   ν   0;
                ν   1   0;
                0   0  (1-ν)/2
            ]
            rowsOfB, b = 3, problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            D = E / ((1 + ν) * (1 - 2ν)) * @SMatrix [
                1-ν   ν    0;
                ν    1-ν   0;
                0     0   (1-2ν)/2
            ]
            rowsOfB, b = 3, 1
        else
            error("stiffnessMatrixSolidParallel: dimension=$(problem.dim), type=$(problem.type)")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # érintett elemtípusok
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x -> push!(types_in_group, x), t)
        end

        for et in types_in_group
            elementName, edim, order, numNodes::Int, _, _ =
                gmsh.model.mesh.getElementProperties(et)

            intPoints, intWeights =
                gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order+1))
            numIntPoints = length(intWeights)

            comp, dfun, ori =
                gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
            ∇h_ref = reshape(dfun, :, numIntPoints)

            elemTags_global, elemNodeTags_global =
                gmsh.model.mesh.getElementsByType(et)
            jac_all, det_all, coords_all =
                gmsh.model.mesh.getJacobians(et, intPoints)
            nip = length(det_all) ÷ length(elemTags_global)

            # gyors tag→index map
            idxmap = Dict{Int,Int}()
            sizehint!(idxmap, length(elemTags_global))
            @inbounds for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            # csoportonként
            for (edim, etag) in dimTags
                if edim != problem.dim; continue; end
                elemTags_local, elemNodeTags_local =
                    gmsh.model.mesh.getElementsByType(et, etag)
                if isempty(elemTags_local); continue; end

                nloc     = length(elemTags_local)
                dof_elem = pdim * numNodes

                # szálankénti puffer
                nt  = Threads.nthreads()
                Is  = [Int[] for _ in 1:nt]
                Js  = [Int[] for _ in 1:nt]
                Vs  = [Float64[] for _ in 1:nt]
                for t in 1:nt
                    sizehint!(Is[t], dof_elem^2 * cld(nloc, nt))
                    sizehint!(Js[t], dof_elem^2 * cld(nloc, nt))
                    sizehint!(Vs[t], dof_elem^2 * cld(nloc, nt))
                end

                Iidx = [l for k in 1:dof_elem, l in 1:dof_elem]
                Jidx = [k for k in 1:dof_elem, l in 1:dof_elem]

                # workspace-ek létrehozása
                W = [make_ws(dim, numNodes, numIntPoints, rowsOfB, dof_elem) for _ in 1:nt]

                @batch per=thread for j in 1:nloc
                    tid = Threads.threadid()
                    w   = W[tid]
                    Iᵗ, Jᵗ, Vᵗ = Is[tid], Js[tid], Vs[tid]

                    fill!(w.K1, 0.0)
                    fill!(w.nn2, 0)

                    # node-tag-ek
                    for k in 1:numNodes
                        w.nnet[k] = elemNodeTags_local[(j-1)*numNodes+k]
                    end

                    # globális index + jac szelet
                    gidx      = idxmap[elemTags_local[j]]
                    jac_slice = view(jac_all, gidx*9*nip + 1 : (gidx+1)*9*nip)
                    det_slice = view(det_all, gidx*nip + 1   : (gidx+1)*nip)

                    fill!(w.dH, 0.0)
                    for k in 1:numIntPoints
                        Jk    = SMatrix{3,3,Float64}(jac_slice[(k-1)*9+1 : k*9])
                        invJT = inv(Jk)'

                        for l in 1:numNodes
                            w.dH[1:dim, (k-1)*numNodes+l] =
                                invJT[1:dim,1:dim] * ∇h_ref[l*3-2 : l*3-(3-dim), k]
                        end
                    end

                    fill!(w.B, 0.0)
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            w.B[k*rowsOfB-0, l*pdim-0] = w.dH[1,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-2, l*pdim-1] = w.dH[1,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-0, l*pdim-1] = w.dH[2,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-1, l*pdim-0] = w.dH[2,(k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            w.B[k*rowsOfB-5, l*pdim-2] = w.dH[1,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-2, l*pdim-1] = w.dH[1,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-0, l*pdim-0] = w.dH[1,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-4, l*pdim-1] = w.dH[2,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-2, l*pdim-2] = w.dH[2,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-1, l*pdim-0] = w.dH[2,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-3, l*pdim-0] = w.dH[3,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-1, l*pdim-1] = w.dH[3,(k-1)*numNodes+l]
                            w.B[k*rowsOfB-0, l*pdim-2] = w.dH[3,(k-1)*numNodes+l]
                        end
                    end

                    for k in 1:numIntPoints
                        B1 = view(w.B, k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:dof_elem)
                        w.K1 .+= B1' * D * B1 * b * det_slice[k] * intWeights[k]
                    end

                    for c in 1:pdim
                        w.nn2[c:pdim:dof_elem] = pdim .* w.nnet .- (pdim - c)
                    end

                    append!(Iᵗ, w.nn2[Iidx[:]])
                    append!(Jᵗ, w.nn2[Jidx[:]])
                    append!(Vᵗ, w.K1[:])
                end

                append!(I, vcat(Is...))
                append!(J, vcat(Js...))
                append!(V, vcat(Vs...))
            end
        end
    end

    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)

    # --- BLAS szálak visszaállítása ---
    BLAS.set_num_threads(old_blas_threads)

    return SystemMatrix(K, problem)
end

function stiffnessMatrixSolidOld(problem; elements=[])
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
                #display("numIntPoints = $numIntPoints")
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
                #appendlock = ReentrantLock()
                #Threads.@threads for j in 1:length(elemTags[i])
                @inbounds for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    @inbounds for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    jac, jacDet, coord = @inline gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    @inbounds for k in 1:numIntPoints
                        invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    end
                    ∂h .*= 0
                    @inbounds for k in 1:numIntPoints, l in 1:numNodes
                        ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    end
                    B .*= 0
                    if dim == 2 && rowsOfB == 3
                        @inbounds for k in 1:numIntPoints, l in 1:numNodes
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
                    K1 .*= 0
                    @inbounds for k in 1:numIntPoints
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += B1' * D * B1 * b * jacDet[k] * intWeights[k]
                    end
                    @inbounds for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    #Threads.lock(appendlock)
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                    #Threads.unlock(appendlock)
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

function stiffnessMatrixAXI(problem; elements=[])
    gmsh.model.setCurrent(problem.name)

    # I,J,V becslés
    elemTypes_all, elemTags_all, elemNodeTags_all = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum(((div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.dim)^2) *
                      length(elemTags_all[i]) for i in 1:length(elemTags_all))

    I = Int[]; J = Int[]; V = Float64[]
    sizehint!(I, lengthOfIJV); sizehint!(J, lengthOfIJV); sizehint!(V, lengthOfIJV)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim, pdim = problem.dim, problem.pdim

        @assert dim == 2 && problem.type == :AxiSymmetric "stiffnessMatrixAXI: dim=2 és :AxiSymmetric kell."

        # Anyagmátrix (4×4) – [εrr, εzz, εθθ, γrz] sorrenddel
        D = E / ((1 + ν) * (1 - 2ν)) * @SMatrix [
            1-ν  ν    ν    0;
            ν    1-ν  ν    0;
            ν    ν    1-ν  0;
            0    0     0   (1-2ν)/2
        ]
        rowsOfB = 4
        b = 1.0     # axi: “vastagság” tényező normálisan 1

        # a teljes 2D háló csomópont-koordinátái (R, Z a 3D X,Y-ból: R=X, Z=Y)
        nodeTags, ncoord, _ = gmsh.model.mesh.getNodes(dim, -1, true, false)
        ncoord2 = zeros(3 * problem.non)
        ncoord2[nodeTags .* 3 .- 2] .= ncoord[1:3:length(ncoord)]   # R (X)
        ncoord2[nodeTags .* 3 .- 1] .= ncoord[2:3:length(ncoord)]   # Z (Y)
        ncoord2[nodeTags .* 3 .- 0] .= ncoord[3:3:length(ncoord)]   # (unused)

        # A fizikai csoport entitásai
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # érintett elemtípusok
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x->push!(types_in_group, x), t)
        end

        # típusonként
        for et in types_in_group
            elementName, edim, order, numNodes::Int, _, _ =
                gmsh.model.mesh.getElementProperties(et)

            # kvadratúra + alakfüggvények
            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order + 1))
            numIntPoints = length(intWeights)

            _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
            ∇h_ref = reshape(dfun, :, numIntPoints)

            _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h_ref = reshape(fun, :, numIntPoints)

            # típushoz tartozó összes elem jacobiánjai (előre)
            elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)
            jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)
            nip = length(det_all) ÷ length(elemTags_global)

            # gyors elemTag→glob.index map
            idxmap = Dict{Int,Int}(); sizehint!(idxmap, length(elemTags_global))
            @inbounds for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            # a csoport elemei
            for (edim, etag) in dimTags
                edim == dim || continue
                elemTags_local, elemNodeTags_local = gmsh.model.mesh.getElementsByType(et, etag)
                isempty(elemTags_local) && continue

                nloc     = length(elemTags_local)
                dof_elem = pdim * numNodes

                # munkaterek (újrafelhasználható tömbök)
                dH = zeros(dim, numNodes * numIntPoints)
                B  = zeros(rowsOfB * numIntPoints, dof_elem)
                K1 = zeros(dof_elem, dof_elem)
                nn2 = zeros(Int, dof_elem)
                r  = zeros(numIntPoints)
                nnet = zeros(Int, nloc, numNodes)

                # indexlisták (I,J) – kész mintaelemhez
                Iidx = repeat(1:dof_elem, inner=dof_elem)
                Jidx = repeat(1:dof_elem, outer=dof_elem)

                @inbounds for j in 1:nloc
                    # lokális csomópont-címkék
                    for k in 1:numNodes
                        nnet[j,k] = elemNodeTags_local[(j-1)*numNodes+k]
                    end

                    # jacobián-szeletek
                    gidx      = idxmap[elemTags_local[j]]
                    jac_slice = @view jac_all[gidx*9*nip + 1 : (gidx+1)*9*nip]
                    det_slice = @view det_all[gidx*nip + 1   : (gidx+1)*nip]

                    fill!(dH, 0.0)
                    fill!(r, 0.0)

                    @inbounds for k in 1:numIntPoints
                        # 3×3 inv(J)^T
                        Jk    = SMatrix{3,3,Float64}(jac_slice[(k-1)*9+1 : k*9])
                        invJT = inv(Jk)'

                        # ∂h = invJT[1:2,1:2]*∇h_ref(•)
                        @inbounds for l in 1:numNodes
                            dH[1:dim, (k-1)*numNodes+l] =
                                invJT[1:dim,1:dim] * ∇h_ref[l*3-2 : l*3-(3-dim), k]
                        end

                        # sugár r(k) = Σ_l h_l(k) * R_l
                        # (R_l = ncoord2[ node*3-2 ])
                        r[k] = dot(@view(h_ref[:,k]), ncoord2[nnet[j,:] .* 3 .- 2])
                    end

                    # B mátrix felépítése
                    fill!(B, 0.0)
                    @inbounds for k in 1:numIntPoints, l in 1:numNodes
                        # ε_rr
                        B[k*rowsOfB-3, l*pdim-1] = dH[1,(k-1)*numNodes+l]
                        # ε_zz
                        B[k*rowsOfB-2, l*pdim-0] = dH[2,(k-1)*numNodes+l]
                        # ε_θθ = u_r / r
                        B[k*rowsOfB-1, l*pdim-1] = h_ref[l,k] / r[k]
                        # γ_rz
                        B[k*rowsOfB-0, l*pdim-1] = dH[2,(k-1)*numNodes+l]
                        B[k*rowsOfB-0, l*pdim-0] = dH[1,(k-1)*numNodes+l]
                    end

                    # elem-mátrix
                    fill!(K1, 0.0)
                    @inbounds for k in 1:numIntPoints
                        @views B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:dof_elem]
                        K1 .+= (B1' * D * B1) * (2π * r[k] * det_slice[k] * intWeights[k] * b)
                    end

                    # globális indexek (nn2)
                    @inbounds for c in 1:pdim
                        nn2[c:pdim:dof_elem] = pdim .* nnet[j,1:numNodes] .- (pdim - c)
                    end

                    # illesztés
                    append!(I, nn2[Iidx]); append!(J, nn2[Jidx]); append!(V, K1[:])
                end
            end
        end
    end

    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end

function stiffnessMatrixAXIParallel(problem; elements=[])
    # BLAS szálak kikapcs. szálonkénti szorzások miatt
    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    gmsh.model.setCurrent(problem.name)

    # I,J,V becslés
    elemTypes_all, elemTags_all, elemNodeTags_all = gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum(((div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.dim)^2) *
                      length(elemTags_all[i]) for i in 1:length(elemTags_all))
    I = Int[]; J = Int[]; V = Float64[]
    sizehint!(I, lengthOfIJV); sizehint!(J, lengthOfIJV); sizehint!(V, lengthOfIJV)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim, pdim = problem.dim, problem.pdim

        @assert dim == 2 && problem.type == :AxiSymmetric "stiffnessMatrixAXIParallel: dim=2 és :AxiSymmetric kell."

        D = E / ((1 + ν) * (1 - 2ν)) * @SMatrix [
            1-ν  ν    ν    0;
            ν    1-ν  ν    0;
            ν    ν    1-ν  0;
            0     0    0   (1-2ν)/2
        ]
        rowsOfB = 4
        b = 1.0

        # teljes 2D háló csomópont-koordináták (R = X)
        nodeTags, ncoord, _ = gmsh.model.mesh.getNodes(dim, -1, true, false)
        ncoord2 = zeros(3 * problem.non)
        ncoord2[nodeTags .* 3 .- 2] .= ncoord[1:3:length(ncoord)]
        ncoord2[nodeTags .* 3 .- 1] .= ncoord[2:3:length(ncoord)]
        ncoord2[nodeTags .* 3 .- 0] .= ncoord[3:3:length(ncoord)]

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # érintett elemtípusok
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x->push!(types_in_group, x), t)
        end

        for et in types_in_group
            elementName, edim, order, numNodes::Int, _, _ =
                gmsh.model.mesh.getElementProperties(et)

            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order + 1))
            numIntPoints = length(intWeights)

            _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
            ∇h_ref = reshape(dfun, :, numIntPoints)

            _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h_ref = reshape(fun, :, numIntPoints)

            elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)
            jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)
            nip = length(det_all) ÷ length(elemTags_global)

            idxmap = Dict{Int,Int}(); sizehint!(idxmap, length(elemTags_global))
            @inbounds for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            for (edim, etag) in dimTags
                edim == dim || continue
                elemTags_local, elemNodeTags_local = gmsh.model.mesh.getElementsByType(et, etag)
                isempty(elemTags_local) && continue

                nloc     = length(elemTags_local)
                dof_elem = pdim * numNodes

                # szálankénti gyűjtők
                nt = Threads.nthreads()
                Is = [Int[] for _ in 1:nt]
                Js = [Int[] for _ in 1:nt]
                Vs = [Float64[] for _ in 1:nt]
                for t in 1:nt
                    sizehint!(Is[t], dof_elem^2 * cld(nloc, nt))
                    sizehint!(Js[t], dof_elem^2 * cld(nloc, nt))
                    sizehint!(Vs[t], dof_elem^2 * cld(nloc, nt))
                end

                Iidx = repeat(1:dof_elem, inner=dof_elem)
                Jidx = repeat(1:dof_elem, outer=dof_elem)

                # workspace-ek
                W = [make_ws_axi(dim, numNodes, numIntPoints, rowsOfB, dof_elem) for _ in 1:nt]

                @batch per=thread for j in 1:nloc
                    tid = Threads.threadid()
                    w = W[tid]
                    Iᵗ, Jᵗ, Vᵗ = Is[tid], Js[tid], Vs[tid]

                    fill!(w.K1, 0.0); fill!(w.nn2, 0)

                    # csomópont-címkék az elemhez
                    @inbounds for k in 1:numNodes
                        w.nnet[k] = elemNodeTags_local[(j-1)*numNodes+k]
                    end

                    # jac szeletek
                    gidx      = idxmap[elemTags_local[j]]
                    jac_slice = @view jac_all[gidx*9*nip + 1 : (gidx+1)*9*nip]
                    det_slice = @view det_all[gidx*nip + 1   : (gidx+1)*nip]

                    fill!(w.dH, 0.0); fill!(w.r, 0.0)
                    @inbounds for k in 1:numIntPoints
                        Jk    = SMatrix{3,3,Float64}(jac_slice[(k-1)*9+1 : k*9])
                        invJT = inv(Jk)'
                        @inbounds for l in 1:numNodes
                            w.dH[1:dim, (k-1)*numNodes+l] =
                                invJT[1:dim,1:dim] * ∇h_ref[l*3-2 : l*3-(3-dim), k]
                        end
                        w.r[k] = dot(@view(h_ref[:,k]), ncoord2[w.nnet .* 3 .- 2])
                    end

                    fill!(w.B, 0.0)
                    @inbounds for k in 1:numIntPoints, l in 1:numNodes
                        # ε_rr
                        w.B[k*rowsOfB-3, l*pdim-1] = w.dH[1,(k-1)*numNodes+l]
                        # ε_zz
                        w.B[k*rowsOfB-2, l*pdim-0] = w.dH[2,(k-1)*numNodes+l]
                        # ε_θθ
                        w.B[k*rowsOfB-1, l*pdim-1] = h_ref[l,k] / w.r[k]
                        # γ_rz
                        w.B[k*rowsOfB-0, l*pdim-1] = w.dH[2,(k-1)*numNodes+l]
                        w.B[k*rowsOfB-0, l*pdim-0] = w.dH[1,(k-1)*numNodes+l]
                    end

                    @inbounds for k in 1:numIntPoints
                        @views B1 = w.B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:dof_elem]
                        w.K1 .+= (B1' * D * B1) * (2π * w.r[k] * det_slice[k] * intWeights[k] * b)
                    end

                    @inbounds for c in 1:pdim
                        w.nn2[c:pdim:dof_elem] = pdim .* w.nnet .- (pdim - c)
                    end

                    append!(Iᵗ, w.nn2[Iidx]); append!(Jᵗ, w.nn2[Jidx]); append!(Vᵗ, w.K1[:])
                end

                append!(I, vcat(Is...)); append!(J, vcat(Js...)); append!(V, vcat(Vs...))
            end
        end
    end

    dof = problem.pdim * problem.non
    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)

    BLAS.set_num_threads(old_blas_threads)
    return SystemMatrix(K, problem)
end

function stiffnessMatrixAXIOld(problem; elements=[])
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
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 2 && problem.type == :AxiSymmetric
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                ν 1-ν ν 0;
                ν ν 1-ν 0;
                0 0 0 (1-2ν)/2]

            rowsOfB = 4
        else
            error("stiffnessMatrixAxiSymmetric: dimension is $(problem.dim), problem type is $(problem.type).")
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
                    if dim == 2 && rowsOfB == 4
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-3, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = B[k*rowsOfB-0, l*pdim-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-2, l*pdim-1] = h[l, k] / r[k]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        #r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += 2π * B1' * D * B1 * r[k] * jacDet[k] * intWeights[k]
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

function stiffnessMatrixTruss(problem; elements=[])
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(1, -1)
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
        E = problem.material[ipg].E
        A = problem.material[ipg].A
        dim = problem.dim
        pdim = problem.pdim
        
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            if edim ≠ 1
                error("stiffnessMatrixTruss: not 1D elements.")
            end
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(1, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                if numNodes ≠ 2
                    error("stiffnessMatrixTruss: truss element must have exactly two nodes.")
                end
                #intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                #numIntPoints = length(intWeights)
                #comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                #∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                #invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                #∂h = zeros(dim, numNodes * numIntPoints)
                #B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                #appendlock = ReentrantLock()
                #Threads.@threads for j in 1:length(elemTags[i])
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    #jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    #Jac = reshape(jac, 3, :)
                    #for k in 1:numIntPoints
                    #    invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    #end
                    #∂h .*= 0
                    #for k in 1:numIntPoints, l in 1:numNodes
                    #    ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    #end
                    #B .*= 0
                    #if dim == 2 && rowsOfB == 3
                    #    for k in 1:numIntPoints, l in 1:numNodes
                    #        B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                    #    end
                    #elseif dim == 3 && rowsOfB == 6
                    #    for k in 1:numIntPoints, l in 1:numNodes
                    #        B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                    #    end
                    #else
                    #    error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    #end
                    #K1 .*= 0
                    #for k in 1:numIntPoints
                    #    B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                    #    K1 += B1' * D * B1 * b * jacDet[k] * intWeights[k]
                    #end
                    K0 = [1 -1; -1 1]
                    x1 = ncoord2[nnet[j, 1] * 3 .- 2]
                    y1 = ncoord2[nnet[j, 1] * 3 .- 1]
                    z1 = ncoord2[nnet[j, 1] * 3 .- 0]
                    x2 = ncoord2[nnet[j, 2] * 3 .- 2]
                    y2 = ncoord2[nnet[j, 2] * 3 .- 1]
                    z2 = ncoord2[nnet[j, 2] * 3 .- 0]
                    L = √((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
                    c1 = (x2 - x1) / L
                    c2 = (y2 - y1) / L
                    c3 = (z2 - z1) / L
                    T = [c1 c2 c3 0 0 0; 0 0 0 c1 c2 c3]
                    k1 = A * E / L
                    K1 = k1 * T' * K0 * T
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    #Threads.lock(appendlock)
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, K1[:])
                    #Threads.unlock(appendlock)
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

"""
    nonLinearStiffnessMatrix(problem, q)

Solves the nonlinear stiffness matrix of the `problem`. `q` is a displacement field.

Returns: `stiffMat`

Types:
- `problem`: Problem
- `q`: VectorField
- `stiffMat`: SystemMatrix
"""
function nonLinearStiffnessMatrix(q; elements=[])
    problem = q.model
    if problem.type == :AxiSymmetric
        return nonLinearStiffnessMatrixAXI(problem, elements=elements)
    else
        return nonLinearStiffnessMatrixSolid(q, elements=elements)
    end
end

function nonLinearStiffnessMatrixSolid(q; elements=[])
    problem = q.model
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
    S1 = zeros(problem.dim, problem.dim)

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
                ∂h = zeros(dim, numNodes * numIntPoints)
                ∂H = zeros(dim * numIntPoints, numNodes)
                B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                K1 = zeros(pdim * numNodes, pdim * numNodes)
                K0 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                #for k in 1:numIntPoints, l in 1:numNodes
                #    for kk in 1:pdim
                #        H[k*pdim-(pdim-kk), l*pdim-(pdim-kk)] = h[(k-1)*numNodes+l]
                #    end
                #end
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
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
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    B .*= 0
                    ∂H .*= 0
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            ∂H[k*pdim*dim-3, l*pdim-1] = ∂H[k*pdim*dim-1, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            ∂H[k*pdim*dim-2, l*pdim-1] = ∂H[k*pdim*dim-0, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                            ∂H[k*dim-2, l] = ∂h[1, (k-1)*numNodes+l]
                            ∂H[k*dim-1, l] = ∂h[2, (k-1)*numNodes+l]
                            ∂H[k*dim-0, l] = ∂h[3, (k-1)*numNodes+l]
                        end
                    else
                        error("nonLinearStiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    q1 = q.a[nn2]
                    for k in 1:numIntPoints
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        ∂H1 = ∂H[k*dim-(dim-1):k*dim, 1:numNodes]
                        σ1 = D * B1 * q1
                        if problem.type == :Solid
                            S1[1,1] = σ1[1]
                            S1[2,2] = σ1[2]
                            S1[3,3] = σ1[3]
                            S1[1,2] = S1[2,1] = σ1[4]
                            S1[2,3] = S1[3,2] = σ1[5]
                            S1[3,1] = S1[1,3] = σ1[6]
                        else
                            error("nonLinearStiffnessMatrix: only 'Solid' is implemented")
                        end
                        K0 = ∂H1' * S1 * ∂H1 * b * jacDet[k] * intWeights[k]
                        for kk in 1:3
                            K1[kk:dim:dim*numNodes, kk:dim:dim*numNodes] += K0
                        end
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

function nonLinearStiffnessMatrixAXI(problem; elements=[])
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
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.pdim
        if problem.dim == 2 && problem.type == :AxiSymmetric
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                ν 1-ν ν 0;
                ν ν 1-ν 0;
                0 0 0 (1-2ν)/2]

            rowsOfB = 4
        else
            error("stiffnessMatrixAxiSymmetric: dimension is $(problem.dim), problem type is $(problem.type).")
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
                    if dim == 2 && rowsOfB == 4
                        for k in 1:numIntPoints, l in 1:numNodes
                            B[k*rowsOfB-3, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-1, l*pdim-0] = B[k*rowsOfB-0, l*pdim-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-2, l*pdim-1] = h[l, k] / r[k]
                        end
                    else
                        error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    K1 .*= 0
                    for k in 1:numIntPoints
                        #r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                        B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                        K1 += 2π * B1' * D * B1 * r[k] * jacDet[k] * intWeights[k]
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

"""
    massMatrix(problem; lumped=...)

Solves the mass matrix of the `problem`. If `lumped` is true, computes the lumped mass matrix.

Returns: `massMat`

Types:
- `problem`: Problem
- `lumped`: Boolean
- `massMat`: SystemMatrix

# Examples

```julia
M = massMatrix(problem; lumped=true)
```
"""
function massMatrix(problem; lumped=true)
    if problem.type == :Truss
        return massMatrixTruss(problem, lumped=lumped)
    elseif Threads.nthreads() == 1
        return massMatrixSolid(problem, lumped=lumped)
    elseif Threads.nthreads() > 1
        return massMatrixSolidParallel(problem, lumped=lumped)
        #return massMatrixSolid(problem, lumped=lumped)
    end
end

function massMatrixSolid(problem; elements = [], lumped::Bool = true)
    gmsh.model.setCurrent(problem.name)

    # --- I,J,V becslés a jobb allokációkhoz ---
    elemTypes_all, elemTags_all, elemNodeTags_all =
        gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum(((div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.pdim)^2) *
                      length(elemTags_all[i]) for i in 1:length(elemTags_all))

    I = Int[]; J = Int[]; V = Float64[]
    sizehint!(I, lengthOfIJV); sizehint!(J, lengthOfIJV); sizehint!(V, lengthOfIJV)

    # --- anyagcsoportonként ---
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        ρ   = problem.material[ipg].ρ
        dim = problem.dim
        pdim = problem.pdim

        rowsOfH, b =
            if     dim == 3 && problem.type == :Solid       ; (3, 1.0)
            elseif dim == 2 && problem.type == :PlaneStress ; (2, problem.thickness)
            elseif dim == 2 && problem.type == :PlaneStrain ; (2, 1.0)
            elseif dim == 2 && problem.type == :AxiSymmetric; (2, 1.0)
            else
                error("massMatrixSolid: unsupported dim=$(problem.dim), type=$(problem.type)")
            end

        # AXI-hoz szükséges a R koordináta (R=X)
        ncoord2 = zeros(3 * problem.non)
        if problem.type == :AxiSymmetric
            nodeTags, ncoord, _ = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags .* 3 .- 2] .= ncoord[1:3:length(ncoord)] # R
            ncoord2[nodeTags .* 3 .- 1] .= ncoord[2:3:length(ncoord)] # Z
            ncoord2[nodeTags .* 3 .- 0] .= ncoord[3:3:length(ncoord)]
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # érintett elemtípusok a csoportban
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x->push!(types_in_group, x), t)
        end

        # --- típusonként ---
        for et in types_in_group
            elementName, edim, order, numNodes::Int, _, _ =
                gmsh.model.mesh.getElementProperties(et)

            # kvadratúra + Lagrange alakfüggvény
            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order + 1))
            numIntPoints = length(intWeights)

            _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h_ref = reshape(fun, :, numIntPoints)

            # H mátrix (csak a h értékeket ismételjük a pdim komponensre)
            dof_elem = pdim * numNodes
            H = zeros(rowsOfH * numIntPoints, dof_elem)
            @inbounds for k in 1:numIntPoints, l in 1:numNodes
                val = h_ref[l, k]
                for c in 1:pdim
                    # 2D-ben rowsOfH==pdim (PlaneStress/Strain/Axi); 3D-ben rowsOfH=3, pdim=3 → stimmel
                    H[k*pdim-(pdim-c), l*pdim-(pdim-c)] = val
                end
            end

            # Jacobianok előre a típus összes elemére
            elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)
            jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)
            nip = length(det_all) ÷ length(elemTags_global)

            # elemTag -> globális index
            idxmap = Dict{Int,Int}(); sizehint!(idxmap, length(elemTags_global))
            @inbounds for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            # --- a csoport elemei ---
            for (edim, etag) in dimTags
                edim == dim || continue
                elemTags_local, elemNodeTags_local =
                    gmsh.model.mesh.getElementsByType(et, etag)
                isempty(elemTags_local) && continue

                nloc = length(elemTags_local)
                nnet = zeros(Int, nloc, numNodes)

                # fix indexminták (I,J)
                Iidx = repeat(1:dof_elem, inner=dof_elem)
                Jidx = repeat(1:dof_elem, outer=dof_elem)

                # munka pufferek
                M1  = zeros(dof_elem, dof_elem)
                nn2 = zeros(Int, dof_elem)
                r   = (problem.type == :AxiSymmetric) ? zeros(numIntPoints) : Float64[]

                @inbounds for j in 1:nloc
                    # node címkék
                    for k in 1:numNodes
                        nnet[j,k] = elemNodeTags_local[(j-1)*numNodes+k]
                    end

                    # Jacobian szeletek
                    gidx      = idxmap[elemTags_local[j]]
                    det_slice = @view det_all[gidx*nip + 1 : (gidx+1)*nip]

                    fill!(M1, 0.0)

                    if problem.type != :AxiSymmetric
                        @inbounds for k in 1:numIntPoints
                            @views H1 = H[k*pdim-(pdim-1):k*pdim, 1:dof_elem]
                            # ρ*b * ∫ HᵀH detJ w
                            M1 .+= (H1' * H1) * (det_slice[k] * intWeights[k])
                        end
                        M1 .*= (ρ * b)
                    else
                        # axi: r(k) kiszámítása
                        @inbounds for k in 1:numIntPoints
                            r[k] = dot(@view(h_ref[:,k]), ncoord2[nnet[j,:] .* 3 .- 2]) # R=X
                        end
                        @inbounds for k in 1:numIntPoints
                            @views H1 = H[k*pdim-(pdim-1):k*pdim, 1:dof_elem]
                            # 2π ρ b ∫ HᵀH r detJ w
                            M1 .+= (H1' * H1) * (r[k] * det_slice[k] * intWeights[k])
                        end
                        M1 .*= (2π * ρ * b)
                    end

                    # globális indexelés
                    @inbounds for c in 1:pdim
                        nn2[c:pdim:dof_elem] = pdim .* nnet[j,1:numNodes] .- (pdim - c)
                    end

                    append!(I, nn2[Iidx]); append!(J, nn2[Jidx]); append!(V, M1[:])
                end
            end
        end
    end

    dof = problem.pdim * problem.non
    M = sparse(I, J, V, dof, dof)
    if lumped
        M = spdiagm(vec(sum(M, dims = 2)))  # diagonál (row-sum)
    end
    dropzeros!(M)
    return SystemMatrix(M, problem)
end

function massMatrixSolid0(problem; elements=[], lumped=true)
    gmsh.model.setCurrent(problem.name)

    # felső becslés I,J,V hosszára
    elemTypes_all, elemTags_all, elemNodeTags_all =
        gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.dim)^2 *
                       length(elemTags_all[i]) for i in 1:length(elemTags_all)])

    I  = Int[]
    J  = Int[]
    V  = Float64[]
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)

    # globális koordináták (axi esethez kell R=x)
    nodeTags_all, ncoord_all, _ = gmsh.model.mesh.getNodes(problem.dim, -1, true, false)
    ncoord2 = zeros(3 * problem.non)
    ncoord2[nodeTags_all .* 3 .- 2] = ncoord_all[1:3:end]
    ncoord2[nodeTags_all .* 3 .- 1] = ncoord_all[2:3:end]
    ncoord2[nodeTags_all .* 3 .- 0] = ncoord_all[3:3:end]

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        ρ      = problem.material[ipg].ρ
        dim    = problem.dim
        pdim   = problem.pdim

        rowsOfH, b = if problem.dim == 3 && problem.type == :Solid
            (3, 1.0)
        elseif problem.dim == 2 && problem.type == :PlaneStress
            (2, problem.thickness)
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            (2, 1.0)
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
            (2, 1.0)
        else
            error("massMatrixSolid: dimension=$(problem.dim), type=$(problem.type)")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # érintett elemtípusok
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x->push!(types_in_group, x), t)
        end

        for et in types_in_group
            elementName, edim, order, numNodes::Int, _, _ =
                gmsh.model.mesh.getElementProperties(et)

            # kvadratúra + alakfüggvények
            intPoints, intWeights =
                gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order + 1))
            numIntPoints = length(intWeights)

            _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h_ref = reshape(fun, :, numIntPoints)   # (numNodes × numIntPoints)

            # globális (típusonként) elemlisták + jacobianok bulkban
            elemTags_global, elemNodeTags_global = gmsh.model.mesh.getElementsByType(et)
            jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)
            nip = length(det_all) ÷ length(elemTags_global)

            # elemTag → 0-based index a jac_all/det_all szeleteléshez
            idxmap = Dict{Int,Int}()
            for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            # végig a fizikai csoport elemein
            for (edim2, etag) in dimTags
                if edim2 != problem.dim; continue; end
                elemTags_local, elemNodeTags_local =
                    gmsh.model.mesh.getElementsByType(et, etag)
                if isempty(elemTags_local); continue; end

                nloc     = length(elemTags_local)
                dof_elem = pdim * numNodes

                # scatter indexek
                Iidx = [l for k in 1:dof_elem, l in 1:dof_elem]
                Jidx = [k for k in 1:dof_elem, l in 1:dof_elem]

                # H mátrix előkészítése (shape-függvény × dof)
                H = zeros(pdim*numIntPoints, dof_elem)
                for k in 1:numIntPoints, l in 1:numNodes
                    for kk in 1:pdim
                        H[(k-1)*pdim + kk, l*pdim - (pdim-kk)] = h_ref[l,k]
                    end
                end

                @inbounds for j in 1:nloc
                    # elem csomópont tagjei
                    nnet = zeros(Int, numNodes)
                    @inbounds for k in 1:numNodes
                        nnet[k] = elemNodeTags_local[(j-1)*numNodes + k]
                    end

                    # jacobik szeletei
                    gidx      = idxmap[elemTags_local[j]]
                    jac_slice = @view jac_all[gidx*9*nip + 1 : (gidx+1)*9*nip]
                    det_slice = @view det_all[gidx*nip + 1   : (gidx+1)*nip]

                    # helyi M1
                    M1 = MMatrix{dof_elem,dof_elem,Float64}(undef)
                    fill!(M1, 0.0)

                    if problem.type != :AxiSymmetric
                        @inbounds for k in 1:numIntPoints
                            H1 = SMatrix{pdim,dof_elem,Float64}(H[(k-1)*pdim+1:k*pdim, 1:dof_elem])
                            G  = H1' * H1 * det_slice[k] * intWeights[k]
                            M1 .+= G
                        end
                        M1 .*= ρ * b
                    else
                        @inbounds for k in 1:numIntPoints
                            # r(IP) = Σ_l h_l(IP) * R_l
                            acc = 0.0
                            @inbounds for l in 1:numNodes
                                Rl = ncoord2[nnet[l]*3 - 2]
                                acc += h_ref[l,k] * Rl
                            end
                            r = acc
                            H1 = SMatrix{pdim,dof_elem,Float64}(H[(k-1)*pdim+1:k*pdim, 1:dof_elem])
                            G  = H1' * H1 * det_slice[k] * intWeights[k] * r
                            M1 .+= G
                        end
                        M1 .*= 2π * ρ * b
                    end

                    # globális dof-indexek
                    nn2 = MVector{dof_elem,Int}(undef)
                    @inbounds for c in 1:pdim
                        nn2[c:pdim:dof_elem] = pdim .* nnet .- (pdim - c)
                    end

                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, M1[:])
                end
            end
        end
    end

    dof = problem.pdim * problem.non
    M = sparse(I, J, V, dof, dof)
    if lumped
        M = spdiagm(vec(sum(M, dims=2)))
    end
    dropzeros!(M)
    return SystemMatrix(M, problem)
end

function massMatrixSolidParallel(problem; elements = [], lumped::Bool = true)
    # BLAS szálak mentése + 1 szálra állítás
    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    gmsh.model.setCurrent(problem.name)

    # --- I,J,V becslés ---
    elemTypes_all, elemTags_all, elemNodeTags_all =
        gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum(((div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.pdim)^2) *
                      length(elemTags_all[i]) for i in 1:length(elemTags_all))

    I = Int[]; J = Int[]; V = Float64[]
    sizehint!(I, lengthOfIJV); sizehint!(J, lengthOfIJV); sizehint!(V, lengthOfIJV)

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        ρ   = problem.material[ipg].ρ
        dim = problem.dim
        pdim = problem.pdim

        rowsOfH, b =
            if     dim == 3 && problem.type == :Solid       ; (3, 1.0)
            elseif dim == 2 && problem.type == :PlaneStress ; (2, problem.thickness)
            elseif dim == 2 && problem.type == :PlaneStrain ; (2, 1.0)
            elseif dim == 2 && problem.type == :AxiSymmetric; (2, 1.0)
            else
                error("massMatrixSolidParallel: unsupported dim=$(problem.dim), type=$(problem.type)")
            end

        # AXI-hoz kell a R koordináta
        ncoord2 = zeros(3 * problem.non)
        if problem.type == :AxiSymmetric
            nodeTags, ncoord, _ = gmsh.model.mesh.getNodes(dim, -1, true, false)
            ncoord2[nodeTags .* 3 .- 2] .= ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags .* 3 .- 1] .= ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags .* 3 .- 0] .= ncoord[3:3:length(ncoord)]
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # érintett elemtípusok
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x->push!(types_in_group, x), t)
        end

        for et in types_in_group
            elementName, edim, order, numNodes::Int, _, _ =
                gmsh.model.mesh.getElementProperties(et)

            # kvadratúra + Lagrange
            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order + 1))
            numIntPoints = length(intWeights)

            _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h_ref = reshape(fun, :, numIntPoints)

            dof_elem = pdim * numNodes

            # H mátrix (elem-típus-specifikus, minden elemhez ugyanaz)
            H = zeros(rowsOfH * numIntPoints, dof_elem)
            @inbounds for k in 1:numIntPoints, l in 1:numNodes
                val = h_ref[l, k]
                for c in 1:pdim
                    H[k*pdim-(pdim-c), l*pdim-(pdim-c)] = val
                end
            end

            # Jacobianok előre
            elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)
            jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)
            nip = length(det_all) ÷ length(elemTags_global)

            idxmap = Dict{Int,Int}(); sizehint!(idxmap, length(elemTags_global))
            @inbounds for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            # --- a csoport elemei ---
            for (edim, etag) in dimTags
                edim == dim || continue
                elemTags_local, elemNodeTags_local =
                    gmsh.model.mesh.getElementsByType(et, etag)
                isempty(elemTags_local) && continue

                nloc = length(elemTags_local)

                # szálankénti gyűjtők
                nt = Threads.nthreads()
                Is = [Int[] for _ in 1:nt]
                Js = [Int[] for _ in 1:nt]
                Vs = [Float64[] for _ in 1:nt]
                for t in 1:nt
                    sizehint!(Is[t], dof_elem^2 * cld(nloc, nt))
                    sizehint!(Js[t], dof_elem^2 * cld(nloc, nt))
                    sizehint!(Vs[t], dof_elem^2 * cld(nloc, nt))
                end

                # minták
                Iidx = repeat(1:dof_elem, inner=dof_elem)
                Jidx = repeat(1:dof_elem, outer=dof_elem)

                # per-thread workspace-ek
                M1s  = [zeros(dof_elem, dof_elem) for _ in 1:nt]
                nn2s = [zeros(Int, dof_elem)      for _ in 1:nt]
                rs   = (problem.type == :AxiSymmetric) ? [zeros(numIntPoints) for _ in 1:nt] : Vector{Vector{Float64}}()
                nnet_s = [zeros(Int, numNodes) for _ in 1:nt]

                @batch per=thread for j in 1:nloc
                    tid = Threads.threadid()
                    Iᵗ, Jᵗ, Vᵗ = Is[tid], Js[tid], Vs[tid]
                    M1 = M1s[tid]; nn2 = nn2s[tid]; nnet = nnet_s[tid]
                    (problem.type == :AxiSymmetric) && fill!(rs[tid], 0.0)

                    # node címkék
                    @inbounds for k in 1:numNodes
                        nnet[k] = elemNodeTags_local[(j-1)*numNodes+k]
                    end

                    # jac szeletek
                    gidx      = idxmap[elemTags_local[j]]
                    det_slice = @view det_all[gidx*nip + 1 : (gidx+1)*nip]

                    fill!(M1, 0.0)

                    if problem.type != :AxiSymmetric
                        @inbounds for k in 1:numIntPoints
                            @views H1 = H[k*pdim-(pdim-1):k*pdim, 1:dof_elem]
                            M1 .+= (H1' * H1) * (det_slice[k] * intWeights[k])
                        end
                        M1 .*= (ρ * b)
                    else
                        rloc = rs[tid]
                        @inbounds for k in 1:numIntPoints
                            rloc[k] = dot(@view(h_ref[:,k]), ncoord2[nnet .* 3 .- 2])
                        end
                        @inbounds for k in 1:numIntPoints
                            @views H1 = H[k*pdim-(pdim-1):k*pdim, 1:dof_elem]
                            M1 .+= (H1' * H1) * (rloc[k] * det_slice[k] * intWeights[k])
                        end
                        M1 .*= (2π * ρ * b)
                    end

                    # globális indexek
                    @inbounds for c in 1:pdim
                        nn2[c:pdim:dof_elem] = pdim .* nnet .- (pdim - c)
                    end

                    append!(Iᵗ, nn2[Iidx]); append!(Jᵗ, nn2[Jidx]); append!(Vᵗ, M1[:])
                end

                append!(I, vcat(Is...)); append!(J, vcat(Js...)); append!(V, vcat(Vs...))
            end
        end
    end

    dof = problem.pdim * problem.non
    M = sparse(I, J, V, dof, dof)
    if lumped
        M = spdiagm(vec(sum(M, dims = 2)))
    end
    dropzeros!(M)

    BLAS.set_num_threads(old_blas_threads)
    return SystemMatrix(M, problem)
end

function massMatrixSolidParallel0(problem; elements=[], lumped=true)
    # --- BLAS szálak mentése + átállítása ---
    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    gmsh.model.setCurrent(problem.name)

    # felső becslés I,J,V hosszára
    elemTypes_all, elemTags_all, elemNodeTags_all =
        gmsh.model.mesh.getElements(problem.dim, -1)
    lengthOfIJV = sum([(div(length(elemNodeTags_all[i]), length(elemTags_all[i])) * problem.dim)^2 *
                       length(elemTags_all[i]) for i in 1:length(elemTags_all)])

    I  = Int[]
    J  = Int[]
    V  = Float64[]
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)

    # Globális csomópont-koordináták (axi-hoz kell az R=x komponens)
    nodeTags_all, ncoord_all, _ = gmsh.model.mesh.getNodes(problem.dim, -1, true, false)
    ncoord2 = zeros(3 * problem.non)
    ncoord2[nodeTags_all .* 3 .- 2] = ncoord_all[1:3:end]   # R (x)
    ncoord2[nodeTags_all .* 3 .- 1] = ncoord_all[2:3:end]   # Z/Y
    ncoord2[nodeTags_all .* 3 .- 0] = ncoord_all[3:3:end]   # (z)

    # anyagcsoportonként
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        ρ      = problem.material[ipg].ρ
        dim    = problem.dim
        pdim   = problem.pdim

        rowsOfH, b, is_axi = if problem.dim == 3 && problem.type == :Solid
            (3, 1.0, false)
        elseif problem.dim == 2 && problem.type == :PlaneStress
            (2, problem.thickness, false)
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            (2, 1.0, false)
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
            (2, 1.0, true)
        else
            error("massMatrixSolidParallel: dimension=$(problem.dim), type=$(problem.type)")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # érintett elemtípusok
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x->push!(types_in_group, x), t)
        end

        for et in types_in_group
            elementName, edim, order, numNodes::Int, _, _ =
                gmsh.model.mesh.getElementProperties(et)

            # kvadratúra + alakfüggvények (csak H-hoz kell a "Lagrange")
            intPoints, intWeights =
                gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2*order + 1))
            numIntPoints = length(intWeights)

            _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
            h_ref = reshape(fun, :, numIntPoints)   # (numNodes × numIntPoints)

            # globális, típusonként: elemek és jacobianok bulkban
            elemTags_global, elemNodeTags_global = gmsh.model.mesh.getElementsByType(et)
            jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)
            nip = length(det_all) ÷ length(elemTags_global)

            idxmap = Dict{Int,Int}()             # elemTag -> 0-based index a szeletekhez
            for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            for (edim2, etag) in dimTags
                if edim2 != problem.dim; continue; end
                elemTags_local, elemNodeTags_local =
                    gmsh.model.mesh.getElementsByType(et, etag)
                if isempty(elemTags_local); continue; end

                nloc     = length(elemTags_local)
                dof_elem = pdim * numNodes

                # szálankénti triplet pufferek
                nt  = Threads.nthreads()
                Is  = [Int[] for _ in 1:nt]
                Js  = [Int[] for _ in 1:nt]
                Vs  = [Float64[] for _ in 1:nt]
                for t in 1:nt
                    sizehint!(Is[t], dof_elem^2 * cld(nloc, nt))
                    sizehint!(Js[t], dof_elem^2 * cld(nloc, nt))
                    sizehint!(Vs[t], dof_elem^2 * cld(nloc, nt))
                end

                # scatter indexek (állandó ebben a blokkban)
                Iidx = [l for k in 1:dof_elem, l in 1:dof_elem]
                Jidx = [k for k in 1:dof_elem, l in 1:dof_elem]

                # H mátrix előállítása (numIntPoints*pdim) × dof_elem
                H = zeros(pdim*numIntPoints, dof_elem)
                for k in 1:numIntPoints, l in 1:numNodes
                    @inbounds for kk in 1:pdim
                        H[(k-1)*pdim + kk, l*pdim - (pdim-kk)] = h_ref[l,k]
                    end
                end

                # workspace szálanként
                function make_workspace()
                    (
                        nnet = zeros(Int, numNodes),
                        M1   = MMatrix{dof_elem,dof_elem,Float64}(undef),
                        nn2  = MVector{dof_elem,Int}(undef),
                        rbuf = is_axi ? zeros(Float64, numIntPoints) : Float64[],  # axi: IP-enkénti r
                    )
                end
                W = [make_workspace() for _ in 1:nt]

                @batch per=thread for j in 1:nloc
                    tid = Threads.threadid()
                    w   = W[tid]
                    Iᵗ, Jᵗ, Vᵗ = Is[tid], Js[tid], Vs[tid]

                    fill!(w.M1, 0.0)

                    # elem csomópontok
                    @inbounds for k in 1:numNodes
                        w.nnet[k] = elemNodeTags_local[(j-1)*numNodes + k]
                    end

                    # jacobian szeletek
                    gidx      = idxmap[elemTags_local[j]]
                    det_slice = @view det_all[gidx*nip + 1 : (gidx+1)*nip]

                    if !is_axi
                        # nem axi
                        @inbounds for k in 1:numIntPoints
                            H1 = SMatrix{pdim,dof_elem,Float64}(H[(k-1)*pdim+1 : k*pdim, 1:dof_elem])
                            G  = H1' * H1 * det_slice[k] * intWeights[k]
                            # in-place
                            @inbounds for p in 1:dof_elem
                                @simd for q in 1:dof_elem
                                    w.M1[p,q] += G[p,q]
                                end
                            end
                        end
                        w.M1 .*= ρ * b
                    else
                        # axi: r(IP) = Σ_l h_l(IP) * R_l
                        @inbounds for k in 1:numIntPoints
                            acc = 0.0
                            @inbounds for l in 1:numNodes
                                acc += h_ref[l,k] * ncoord2[w.nnet[l]*3 - 2]  # R_l
                            end
                            w.rbuf[k] = acc
                        end
                        @inbounds for k in 1:numIntPoints
                            H1 = SMatrix{pdim,dof_elem,Float64}(H[(k-1)*pdim+1 : k*pdim, 1:dof_elem])
                            G  = H1' * H1 * det_slice[k] * intWeights[k] * w.rbuf[k]
                            @inbounds for p in 1:dof_elem
                                @simd for q in 1:dof_elem
                                    w.M1[p,q] += G[p,q]
                                end
                            end
                        end
                        w.M1 .*= 2π * ρ * b
                    end

                    # globális dof-indexek
                    @inbounds for c in 1:pdim
                        w.nn2[c:pdim:dof_elem] = pdim .* w.nnet .- (pdim - c)
                    end

                    append!(Iᵗ, w.nn2[Iidx[:]])
                    append!(Jᵗ, w.nn2[Jidx[:]])
                    append!(Vᵗ, w.M1[:])
                end # @batch

                # szálpufferek összefésülése
                append!(I, vcat(Is...))
                append!(J, vcat(Js...))
                append!(V, vcat(Vs...))
            end
        end
    end

    dof = problem.pdim * problem.non
    M = sparse(I, J, V, dof, dof)
    if lumped
        M = spdiagm(vec(sum(M, dims=2)))
    end
    dropzeros!(M)

    # --- BLAS szálak visszaállítása ---
    BLAS.set_num_threads(old_blas_threads)

    return SystemMatrix(M, problem)
end

function massMatrixSolidOld(problem; elements=[], lumped=true)
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
        ρ = problem.material[ipg].ρ
        if problem.dim == 3 && problem.type == :Solid
            rowsOfH = 3
            b = 1
        elseif problem.dim == 2 && problem.type == :PlaneStress
            rowsOfH = 2
            b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            rowsOfH = 2
            b = 1
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
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
                if problem.type != :AxiSymmetric
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
                        M1 *= ρ * b
                        for k in 1:pdim
                            nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                        end
                        append!(I, nn2[Iidx[:]])
                        append!(J, nn2[Jidx[:]])
                        append!(V, M1[:])
                    end
                elseif problem.type == :AxiSymmetric
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
                        M1 *= 2π * ρ * b
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
    return SystemMatrix(M, problem)
end

function massMatrixTruss(problem; lumped=false)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(1, -1)
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
        ρ = problem.material[ipg].ρ
        A = problem.material[ipg].A
        dim = problem.dim
        pdim = problem.pdim
        
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            if edim ≠ 1
                error("stiffnessMatrixTruss: not 1D elements.")
            end
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(1, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                if numNodes ≠ 2
                    error("stiffnessMatrixTruss: truss element must have exactly two nodes.")
                end
                #intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                #numIntPoints = length(intWeights)
                #comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
                #∇h = reshape(dfun, :, numIntPoints)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                #invJac = zeros(3, 3numIntPoints)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                #∂h = zeros(dim, numNodes * numIntPoints)
                #B = zeros(rowsOfB * numIntPoints, pdim * numNodes)
                M1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                #appendlock = ReentrantLock()
                #Threads.@threads for j in 1:length(elemTags[i])
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    #jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    #Jac = reshape(jac, 3, :)
                    #for k in 1:numIntPoints
                    #    invJac[1:3, 3*k-2:3*k] = @inline inv(Jac[1:3, 3*k-2:3*k])'
                    #end
                    #∂h .*= 0
                    #for k in 1:numIntPoints, l in 1:numNodes
                    #    ∂h[1:dim, (k-1)*numNodes+l] .= invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k]
                    #end
                    #B .*= 0
                    #if dim == 2 && rowsOfB == 3
                    #    for k in 1:numIntPoints, l in 1:numNodes
                    #        B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                    #    end
                    #elseif dim == 3 && rowsOfB == 6
                    #    for k in 1:numIntPoints, l in 1:numNodes
                    #        B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                    #    end
                    #else
                    #    error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    #end
                    #K1 .*= 0
                    #for k in 1:numIntPoints
                    #    B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:pdim*numNodes]
                    #    K1 += B1' * D * B1 * b * jacDet[k] * intWeights[k]
                    #end
                    M0 = [2 1; 1 2] / 6
                    x1 = ncoord2[nnet[j, 1] * 3 .- 2]
                    y1 = ncoord2[nnet[j, 1] * 3 .- 1]
                    z1 = ncoord2[nnet[j, 1] * 3 .- 0]
                    x2 = ncoord2[nnet[j, 2] * 3 .- 2]
                    y2 = ncoord2[nnet[j, 2] * 3 .- 1]
                    z2 = ncoord2[nnet[j, 2] * 3 .- 0]
                    L = √((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
                    c1 = (x2 - x1) / L
                    c2 = (y2 - y1) / L
                    c3 = (z2 - z1) / L
                    T = [c1 c2 c3 0 0 0; 0 0 0 c1 c2 c3]
                    m1 = A * ρ * L
                    M1 = m1 * T' * M0 * T
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    #Threads.lock(appendlock)
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, M1[:])
                    #Threads.unlock(appendlock)
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
    return SystemMatrix(M, problem)
end

"""
    dampingMatrix(K, M, ωₘₐₓ; α=0.0, ξ=..., β=...)

Generates the damping matrix for proportional damping. C = αM + βK, or
C = αM + β₁K + β₂KM⁻¹K + β₃KM⁻¹KM⁻¹K + ⋯. The latter corresponds to a damping characteristic
given by a power series in the natural frequencies with odd exponents. ξᵢ (`ξ` in the
arguments) are the values of the individual terms of the series at ωₘₐₓ. βᵢ (`β` in the
arguments) are the coefficients of the series. Either `ξ` or `β` must be specified; each may
be a scalar or a vector. `K` is the stiffness matrix, `M` is the mass matrix, and `ωₘₐₓ` is the
largest natural frequency.

Returns: `dampingMatrix`

Types:
- `K`: SystemMatrix
- `M`: SystemMatrix
- `ωₘₐₓ`: Float64
- `α`: Float64
- `ξ`: Float64 or Vector{Float64}
- `β`: Float64 or Vector{Float64}
- `dampingMatrix`: SystemMatrix

# Examples

```julia
K = stiffnessMatrix(problem)
M = massMatrix(problem; lumped=true)
ωmax = 2π * 1000
C = dampingMatrix(K, M, ωmax; α=0.0, ξ=[0.02, 0.02])
```
"""
function dampingMatrix(K, M, ωₘₐₓ; α=0.0, ξ=0.01, β=[2ξ[i]/(ωₘₐₓ)^(2i-1) for i in 1:length(ξ)])
    if K.model != M.model
        error("dampingMatrix: problem of K and M are not the same.")
    end
    dof, dof = size(M)
    dof2, dof2 = size(K)
    if dof != nnz(M)
        error("dampingMatrix: M is not lumped!")
    end
    if dof != dof2
        error("dampingMatrix: sizes of M and K are not match: $dof <--> $dof2!")
    end
    invM = spdiagm(1 ./ diag(M))
    C = spzeros(dof, dof)
    MK = copy(K)
    iMK = invM * K
    C += α * M
    C += β[1] * MK
    for i in 2:length(β)
        MK *= iMK
        C += β[i] * MK
    end
    dropzeros!(C)
    return SystemMatrix(C, K.model)
end

"""
    elasticSupportMatrix(problem, elSupp)

Solves the elastic support matrix of the `problem`. `elSupp` is a vector of elastic
supports defined in function `elasticSupport`. With the displacementent vector `q` in hand the
reaction force vector `fR` arising from the elastic support can be solved.
(`fR = heatConvMat * q`)

Return: `elSuppMat`

Types:
- `problem`: Problem
- `elSupp`: Vector{Tuple{String, Float64, Float64, Float64}}
- `elSuppMat`: SystemMatrix
"""
function elasticSupportMatrix(problem, elSupports)
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

    if !isa(elSupports, Vector)
        error("elasticSupportMatrix: elastic supports are not arranged in a vector. Put them in [...]")
    end
    pdim = problem.pdim
    DIM = problem.dim
    b = problem.thickness
    non = problem.non
    dof = pdim * non
    ncoord2 = zeros(3 * problem.non)
    for n in 1:length(elSupports)
        name, kx, ky, kz = elSupports[n]
        if problem.pdim == 3
            f = [0, 0, 0]
        elseif problem.pdim == 2
            f = [0, 0]
        elseif problem.pdim == 1
            f = [0]
        else
            error("elasticSupportMatrix: dimension of the problem is $(problem.pdim).")
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
                nnet = zeros(Int, length(elementTags[ii]), numNodes)
                Iidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                Jidx = zeros(Int, numNodes * pdim, numNodes * pdim)
                for k in 1:numNodes*pdim, l in 1:numNodes*pdim
                    Iidx[k, l] = l
                    Jidx[k, l] = k
                end
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                        end
                    end
                end
                C1 = zeros(pdim * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    C1 .*= 0
                    for j in 1:numIntPoints
                        x = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 2]
                        y = 0
                        z = 0
                        if isa(kx, Function) || isa(ky, Function) || isa(kz, Function)
                            y = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 1]
                            z = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 0]
                        end
                        f[1] = isa(kx, Function) ? kx(x, y, z) : kx
                        if problem.pdim > 1
                            f[2] = isa(ky, Function) ? ky(x, y, z) : ky
                        end
                        if problem.pdim == 3
                            f[3] = isa(kz, Function) ? kz(x, y, z) : kz
                        end
                        r = x
                        if pdim == 3
                            kk = [f[1] 0 0; 0 f[2] 0; 0 0 f[3]]
                        elseif pdim == 2
                            kk = [f[1] 0; 0 f[2]]
                        elseif pdim == 1
                            kk = [f[1]]
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
                            Ja = 2π * jacDet[j] * r
                        elseif DIM == 2 && dim == 1 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * r
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
                        else
                            error("elasticSupportMatrix: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                        end
                        C1 += H1' * kk * H1 * Ja * intWeights[j]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                    end
                    append!(I, nn2[Iidx[:]])
                    append!(J, nn2[Jidx[:]])
                    append!(V, C1[:])
                end
            end
        end
    end

    C = sparse(I, J, V, dof, dof)
    dropzeros!(C)
    return SystemMatrix(C, problem)
end

"""
    loadVector(problem, loads)

Solves a load vector of `problem`. `loads` is a tuple of name of physical group 
`name`, coordinates `fx`, `fy` and `fz` of the intensity of distributed force.
It can solve traction or body force depending on the problem.
- In case of 2D problems and Point physical group means concentrated force.
- In case of 2D problems and Line physical group means surface force.
- In case of 2D problems and Surface physical group means body force.
- In case of 3D problems and Point physical group means concentrated force.
- In case of 3D problems and Line physical group means edge force.
- In case of 3D problems and Surface physical group means surface force.
- In case of 3D problems and Volume physical group means body force.

Return: `loadVec`

Types:
- `problem`: Problem
- `loads`: Vector{Tuple{String, Float64, Float64, Float64}}
- `loadVec`: VectorField
"""
function loadVector(problem, loads)
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
                nnet = zeros(Int, length(elementTags[ii]), numNodes)
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                        end
                    end
                end
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    f1 .*= 0
                    for j in 1:numIntPoints
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
                        r = x
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
                            Ja = 2π * jacDet[j] * r
                        elseif DIM == 2 && dim == 1 && problem.type != :AxiSymmetric && problem.type != :AxiSymmetricHeatConduction
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * b
                        elseif DIM == 2 && dim == 1 && (problem.type == :AxiSymmetric || problem.type == :AxiSymmetricHeatConduction)
                            Ja = 2π * √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2) * r
                        elseif DIM == 2 && dim == 0
                            Ja = 1
                        ############ 1D #######################################################
                        else
                            error("loadVector: dimension of the problem is $(problem.dim), dimension of load is $dim.")
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
    type = :null
    if pdim == 3
        type = :v3D
    elseif pdim == 2
        type = :v2D
    elseif pdim == 1
        type = :scalar
    else
        error("loadVector: wrong pdim ($pdim).")
    end
    if type == :v3D || type == :v2D
        return VectorField([], reshape(fp, :, 1), [0.0], [], 1, type, problem)
    elseif type == :scalar
        return ScalarField([], reshape(fp, :, 1), [0.0], [], 1, type, problem)
    end
end

"""
    applyBoundaryConditions!(stiffMat, loadVec, supports)

Applies displacement boundary conditions `supports` on a stiffness matrix
`stiffMat` and load vector `loadVec`. Mesh details are in `problem`. `supports`
is a tuple of `name` of physical group and prescribed displacements `ux`, `uy`
and `uz`.

Returns: nothing

Types:
- `stiffMat`: SystemMatrix 
- `loadVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions!(stiffMat::SystemMatrix, loadVec, supports)
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    dof, dof = size(stiffMat.A)
    massMat = SystemMatrix(spzeros(dof, dof), stiffMat.model)
    dampMat = SystemMatrix(spzeros(dof, dof), stiffMat.model)
    applyBoundaryConditions!(stiffMat, massMat, dampMat, loadVec, supports)
    #massMat.A *= 0
    #dampMat.A *= 0
    return
end

"""
    applyBoundaryConditions(stiffMat, loadVec, supports)

Applies displacement boundary conditions `supports` on a stiffness matrix
`stiffMat` and load vector `loadVec`. Mesh details are in `problem`. `supports`
is a tuple of `name` of physical group and prescribed displacements `ux`, `uy`
and `uz`. Creates a new stiffness matrix and load vector.

Returns: `stiffMat`, `loadVec`

Types:
- `stiffMat`: SystemMatrix 
- `loadVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions(stiffMat0, loadVec0, supports)
    problem = stiffMat.model
    if !isa(supports, Vector)
        error("applyBoundaryConditions: supports are not arranged in a vector. Put them in [...]")
    end
    dof, dof = size(stiffMat0.A)
    massMat = SystemMatrix(spzeros(dof, dof), stiffMat0.model)
    dampMat = SystemMatrix(spzeros(dof, dof), stiffMat0.model)
    stiffMat = copy(stiffMat0)
    loadVec = copy(loadVec0)
    applyBoundaryConditions!(stiffMat, massMat, dampMat, loadVec, supports)
    #massMat.A = []
    #dampMat.A = []
    return stiffMat, loadVec
end

"""
    applyBoundaryConditions!(heatCondMat, heatCapMat, heatFluxVec, supports)

Applies boundary conditions `supports` on a heat conduction matrix
`heatCondMat`, heat capacity matrix `heatCapMat` and heat flux vector `heatFluxVec`. Mesh details are in `problem`. `supports`
is a tuple of `name` of physical group and prescribed temperature `T`.

Returns: nothing

Types:
- `stiffMat`: SystemMatrix 
- `loadVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions!(heatCondMat::SystemMatrix, heatCapMat::SystemMatrix, heatFluxVec, supports; fix=1)
    if !isa(supports, Vector)
        error("applyBoundaryConditions: supports are not arranged in a vector. Put them in [...]")
    end
    dof, dof = size(heatCondMat.A)
    dampMat = SystemMatrix(spzeros(dof, dof), heatCondMat.model)
    applyBoundaryConditions!(heatCondMat, heatCapMat, dampMat, heatFluxVec, supports, fix=fix)
    #dampMat.A = []
    return
end

"""
    getTagForPhysicalName(name)

Returns `tags` of elements of physical group `name`.

Returns: `tags`

Types:
- `name`: String
- `tags`: Vector{Integer}
"""
function getTagForPhysicalName(name)
    dimTags = gmsh.model.getPhysicalGroups(-1)
    i = 1
    while gmsh.model.getPhysicalName(dimTags[i][1], dimTags[i][2]) != name
        i += 1
        if i > length(dimTags)
            error("Physical name '$name' does not exist.")
        end
    end
    return dimTags[i][2]
end

"""
    applyBoundaryConditions!(stiffMat, massMat, dampMat, loadVec, supports)

Applies displacement boundary conditions `supports` on a stiffness matrix
`stiffMat`, mass matrix `massMat`, damping matrix `dampMat` and load vector `loadVec`.
Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Returns: nothing

Types:
- `stiffMat`: SystemMatrix 
- `massMat`: SystemMatrix 
- `dampMat`: SystemMatrix 
- `loadVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions!(stiffMat::SystemMatrix, massMat::SystemMatrix, dampMat::SystemMatrix, loadVec, supports; fix=1)
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    if stiffMat.model != massMat.model || massMat.model != dampMat.model
        error("applyBoundaryConditions!: K, M or C does not belong to the same problem.")
    end
    gmsh.model.setCurrent(stiffMat.model.name)
    dof, dof = size(stiffMat.A)
    pdim = stiffMat.model.pdim

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
                f0 = stiffMat.A[:, nodeTagsX] * uux
            else
                f0 = stiffMat.A[:, nodeTagsX] * ux
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
                f0 = stiffMat.A[:, nodeTagsY] * uuy
            else
                f0 = stiffMat.A[:, nodeTagsY] * uy
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
                f0 = stiffMat.A[:, nodeTagsZ] * uuz
            else
                f0 = stiffMat.A[:, nodeTagsZ] * uz
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
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = fix
                massMat.A[j, :] .= 0
                massMat.A[:, j] .= 0
                massMat.A[j, j] = 1
                dampMat.A[j, :] .= 0
                dampMat.A[:, j] .= 0
                dampMat.A[j, j] = fix
                if isa(ux, Function)
                    loadVec.a[j] = uux[jj]
                else
                    loadVec.a[j] = ux
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
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = fix
                massMat.A[j, :] .= 0
                massMat.A[:, j] .= 0
                massMat.A[j, j] = 1
                dampMat.A[j, :] .= 0
                dampMat.A[:, j] .= 0
                dampMat.A[j, j] = fix
                if isa(uy, Function)
                    loadVec.a[j] = uuy[jj]
                else
                    loadVec.a[j] = uy
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
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = fix
                massMat.A[j, :] .= 0
                massMat.A[:, j] .= 0
                massMat.A[j, j] = 1
                dampMat.A[j, :] .= 0
                dampMat.A[:, j] .= 0
                dampMat.A[j, j] = fix
                if isa(uz, Function)
                    loadVec.a[j] = uuz[jj]
                else
                    loadVec.a[j] = uz
                end
            end
        end
    end

    dropzeros!(stiffMat.A)
    dropzeros!(massMat.A)
    dropzeros!(dampMat.A)
end

"""
    applyBoundaryConditions!(dispVec, supports)

Applies displacement boundary conditions `supports` on a displacement vector
`dispVec`. Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Returns: nothing

Types:
- `problem`: Problem
- `dispVec`: VectorField
- `supports`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyBoundaryConditions!(dispVec::Matrix, supports)
    problem = dispVec.model
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
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
                dispVec[nodeTagsX,:] .= uux
            else
                dispVec[nodeTagsX,:] .= ux
            end
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim - 2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
                dispVec[nodeTagsY,:] .= uuy
            else
                dispVec[nodeTagsY,:] .= uy
            end
        end
        if pdim == 3 && uz != 1im
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
                dispVec[nodeTagsZ,:] .= uuz
            else
                dispVec[nodeTagsZ,:] .= uz
            end
        end
    end
end

function applyBoundaryConditions!(dispVec::VectorField, supports)
    return applyBoundaryConditions!(dispVec.model, dispVec.a, supports)
end

"""
    applyElasticSupport!(stiffMat, elastSupp)

Applies elastic support boundary conditions `elastSupp` on a stiffness matrix
`stiffMat`. Mesh details are in `problem`. `elastSupp` is a tuple of `name`
of physical group and prescribed `kx`, `ky` and `kz` stiffnesses.

Returns: nothing

Types:
- `stiffMat`: SystemMatrix 
- `elastSupp`: Vector{Tuple{String, Float64, Float64, Float64}}
"""
function applyElasticSupport!(stiffMat::SystemMatrix, elastSupp)
    problem = stiffMat.model
    if !isa(elastSupp, Vector)
        error("applyElasticSupport!: elastic supports are not arranged in a vector. Put them in [...]")
    end
    C0 = elasticSupportMatrix(stiffMat.model, elastSupp)
    stiffMat.A .+= C0.A
end

"""
    solveDisplacement(K, q)

Solves the equation K*q=f for the displacement vector `q`. `K` is the stiffness Matrix,
`q` is the load vector.

Return: `q`

Types:
- `K`: SystemMatrix 
- `f`: VectorField 
- `q`: VectorField
"""
function solveDisplacement(K, f)
    type = :null
    if f.type == :v3D
        type = :v3D
    elseif f.type == :v2D
        type = :v2D
    else
        error("solveDisplacement: wrong type of 'f': ($(f.type))")
    end
    return VectorField([], reshape(K.A \ f.a, :, 1), [0.0], [], 1, type, K.model)
end

"""
    solveDisplacement(problem, load, supp)

Solves the displacement vector `q` of `problem` with loads `load` and
supports `supp`.

Return: `q`

Types:
- `problem`: Problem 
- `load`: Vector{Tuple} 
- `supp`: Vector{Tuple}
- `q`: VectorField
"""
function solveDisplacement(problem, load, supp)
    K = stiffnessMatrix(problem)
    f = loadVector(problem, load)
    applyBoundaryConditions!(K, f, supp)
    type = :null
    if f.type == :v3D
        type = :v3D
    elseif f.type == :v2D
        type = :v2D
    else
        error("solveDisplacement: wrong type of 'f': ($(f.type))")
    end
    return VectorField([], reshape(K.A \ f.a, :, 1), [0.0], [], 1, type, problem)
end

"""
    solveDisplacement(problem, load, supp, elasticSupp)

Solves the displacement vector `q` of `problem` with loads `load`, 
supports `supp` and elastic supports `elasticSupp`.

Return: `q`

Types:
- `problem`: Problem 
- `load`: Vector{Tuple} 
- `supp`: Vector{Tuple}
- `q`: VectorField
"""
function solveDisplacement(problem, load, supp, elsupp)
    K = stiffnessMatrix(problem)
    f = loadVector(problem, load)
    applyElasticSupport!(K, elsupp)
    applyBoundaryConditions!(K, f, supp)
    type = :null
    if f.type == :v3D
        type = :v3D
    elseif f.type == :v2D
        type = :v2D
    else
        error("solveDisplacement: wrong type of 'f': ($(f.type))")
    end
    return VectorField([], reshape(K.A \ f.a, :, 1), [0.0], [], 1, type, problem)
end

"""
    solveStrain(q; DoFResults=false)

Solves the strain field `E` from displacement vector `q`. Strain field is given
per elements, so it usually contains jumps at the boundaries of elements. Details
of mesh is available in `problem`. If `DoFResults` is true, `E` is a matrix with
nodal results. In this case `showDoFResults` can be used to show the results 
(otherwise `showStrainResults` or `showElementResults`).

Return: `E`

Types:
- `q`: VectorField
- `E`: TensorField
"""
function solveStrain(q; DoFResults::Bool=false)
    problem = q.model
    gmsh.model.setCurrent(problem.name)

    # --- ellenőrzések ---
    if !(q isa VectorField)
        error("solveStrain: argument must be a VectorField. Now it is '$(typeof(q))'.")
    end
    if q.type != :v3D && q.type != :v2D
        error("solveStrain: argument must be a displacement vector (v2D/v3D). Now it is '$(q.type)'.")
    end
    if q.A != []
        error("solveStrain: q.A != []")
    end

    nsteps = q.nsteps
    dim_global = problem.dim
    pdim = problem.pdim
    non  = problem.non

    # kimenet gyűjtők
    ε = Vector{Matrix{Float64}}()   # elemenként (9*numNodes × nsteps)
    numElem = Int[]                 # elem címkék
    type = :e

    # DoF-átlagoláshoz (csomóponti strain):
    if DoFResults
        E1  = zeros(non * 9, nsteps) # 9 komponens / csomópont
        pcs = zeros(Int, non)         # hány elem érinti a csomópontot
    end

    # --- anyagcsoportonként (physical group) ---
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        ν = problem.material[ipg].ν   # PlaneStress-hez kell a ε_zz képlethez

        # helyi dim/rowsOfB/b a feladat típusától függően
        local_dim = 0
        rowsOfB = 0
        b = 1.0
        if problem.dim == 3 && problem.type == :Solid
            local_dim = 3; rowsOfB = 6; b = 1.0
        elseif problem.dim == 2 && problem.type == :PlaneStress
            local_dim = 2; rowsOfB = 3; b = problem.thickness
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            local_dim = 2; rowsOfB = 3; b = 1.0
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
            local_dim = 2; rowsOfB = 4; b = 1.0
        else
            error("solveStrain: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        # AXI-hoz kell a (R= X) koordináta
        ncoord2 = zeros(3 * non)
        if problem.type == :AxiSymmetric
            nodeTags, ncoord, _ = gmsh.model.mesh.getNodes(local_dim, -1, true, false)
            ncoord2[nodeTags .* 3 .- 2] .= ncoord[1:3:length(ncoord)]  # R (X)
            ncoord2[nodeTags .* 3 .- 1] .= ncoord[2:3:length(ncoord)]  # Z (Y)
            ncoord2[nodeTags .* 3 .- 0] .= ncoord[3:3:length(ncoord)]
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        # --- érintett elemtípusok összegyűjtése ---
        types_in_group = Set{Int}()
        for (edim, etag) in dimTags
            t, _, _ = gmsh.model.mesh.getElements(edim, etag)
            foreach(x->push!(types_in_group, x), t)
        end

        # --- típusonként ---
        for et in types_in_group
            elementName, edim, order, numNodes::Int, localNodeCoord, _ =
                gmsh.model.mesh.getElementProperties(et)

            # "pontok": itt az ELEMENT CSOMÓPONTOK lokális koordinátái (nodeCoord)
            nodeCoord = zeros(numNodes * 3)
            @inbounds for k in 1:local_dim, j in 1:numNodes
                nodeCoord[k + (j-1)*3] = localNodeCoord[k + (j-1)*local_dim]
            end

            # bázisfüggvények node-helyeken
            _, dfun, _ = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
            ∇h_ref = reshape(dfun, :, numNodes)  # (3*numNodes)×numNodes → rendezve
            _, fun, _ = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
            h_ref = reshape(fun, :, numNodes)    # numNodes × numNodes

            # típus összes elemére Jacobianok előre
            elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)
            jac_all, _, _ = gmsh.model.mesh.getJacobians(et, nodeCoord)  # det nem kell a strainhez
            nip = numNodes  # nodeCoord hosszából következik

            # elemTag → globális index a jac_all-hoz
            idxmap = Dict{Int,Int}(); sizehint!(idxmap, length(elemTags_global))
            @inbounds for (gi, gtag) in enumerate(elemTags_global)
                idxmap[gtag] = gi - 1
            end

            # --- a csoport adott entitásain belüli elemek ---
            for (edim, etag) in dimTags
                edim == local_dim || continue
                elemTags_local, elemNodeTags_local = gmsh.model.mesh.getElementsByType(et, etag)
                isempty(elemTags_local) && continue

                nloc = length(elemTags_local)
                dof_elem = pdim * numNodes

                # munkaterek
                invJac = zeros(3, 3*numNodes)                    # 3×(3*numNodes)
                dH     = zeros(local_dim, numNodes * numNodes)   # ∂h (∂/x, ∂/y[,∂/z]) minden csomópontra
                B      = zeros(rowsOfB * numNodes, dof_elem)     # B-mátrix csomópontokra
                nn2    = zeros(Int, dof_elem)
                r      = (problem.type == :AxiSymmetric) ? zeros(numNodes) : Float64[]

                # gyors csomóponti index-minta a B1 kivágásához
                # csomópont k-hoz tartozó rows: (k-1)*rowsOfB+1 : k*rowsOfB
                # dof-oszlopok: 1:dof_elem

                # feldolgozás elemenként
                @inbounds for j in 1:nloc
                    elem = elemTags_local[j]

                    # node-tag-ek
                    nnet_j = similar(nn2, Int, numNodes)
                    for k in 1:numNodes
                        nnet_j[k] = elemNodeTags_local[(j-1)*numNodes + k]
                    end

                    # Jacobian szelet és inverzei node-pontonként
                    gidx = idxmap[elem]
                    jac_slice = @view jac_all[gidx*9*nip + 1 : (gidx+1)*9*nip]

                    @inbounds for k in 1:numNodes
                        Jk = @view jac_slice[(k-1)*9+1 : k*9]
                        # 3×3 mátrix
                        J = reshape(copy(Jk), 3, 3)  # kis mátrix → a copy olcsó; kerülhet SMatrix is, de így kisebb a compile-time
                        IJT = permutedims(inv(J))    # inv(J)'  (3×3)
                        invJac[:, 3k-2:3k] .= IJT
                    end

                    # sugár (AXI): r[k] = Σ_l h_ref[l,k] * R_l
                    if problem.type == :AxiSymmetric
                        @inbounds for k in 1:numNodes
                            r[k] = dot(@view(h_ref[:,k]), ncoord2[nnet_j .* 3 .- 2])
                        end
                    end

                    # ∂h mező (lokális → globális deriváltak)
                    fill!(dH, 0.0)
                    @inbounds for k in 1:numNodes, l in 1:numNodes
                        # invJac[1:local_dim, k*3-2 : k*3-(3-local_dim)] * ∇h_ref[l*3-2 : l*3-(3-local_dim), k]
                        dH[1:local_dim, (k-1)*numNodes + l] =
                            @view(invJac[1:local_dim, (k-1)*3+1 : (k-1)*3+local_dim]) *
                            @view(∇h_ref[l*3-2 : l*3-(3-local_dim), k])
                    end

                    # B összeállítás (csomópontokra)
                    fill!(B, 0.0)
                    if local_dim == 2 && rowsOfB == 3 && problem.type == :PlaneStress
                        @inbounds for k in 1:numNodes, l in 1:numNodes
                            # εxx, εyy, γxy
                            B[k*3-2, l*2-1] = dH[1, (k-1)*numNodes + l]
                            B[k*3-1, l*2-0] = dH[2, (k-1)*numNodes + l]
                            B[k*3-0, l*2-1] = dH[2, (k-1)*numNodes + l]
                            B[k*3-0, l*2-0] = dH[1, (k-1)*numNodes + l]
                        end
                    elseif local_dim == 2 && rowsOfB == 3 && problem.type == :PlaneStrain
                        @inbounds for k in 1:numNodes, l in 1:numNodes
                            # εxx, εyy, γxy
                            B[k*3-2, l*2-1] = dH[1, (k-1)*numNodes + l]
                            B[k*3-1, l*2-0] = dH[2, (k-1)*numNodes + l]
                            B[k*3-0, l*2-1] = dH[2, (k-1)*numNodes + l]
                            B[k*3-0, l*2-0] = dH[1, (k-1)*numNodes + l]
                        end
                    elseif local_dim == 3 && rowsOfB == 6 && problem.type == :Solid
                        @inbounds for k in 1:numNodes, l in 1:numNodes
                            # εxx, εyy, εzz, γyz, γxz, γxy (Voigt 6)
                            B[k*6-5, l*3-2] = dH[1, (k-1)*numNodes + l]   # εxx ← ∂u/∂x
                            B[k*6-4, l*3-1] = dH[2, (k-1)*numNodes + l]   # εyy ← ∂v/∂y
                            B[k*6-3, l*3-0] = dH[3, (k-1)*numNodes + l]   # εzz ← ∂w/∂z

                            B[k*6-2, l*3-1] = dH[3, (k-1)*numNodes + l]   # γyz ← ∂w/∂y
                            B[k*6-2, l*3-0] += dH[2, (k-1)*numNodes + l]  #       + ∂v/∂z

                            B[k*6-1, l*3-2] = dH[3, (k-1)*numNodes + l]   # γxz ← ∂w/∂x
                            B[k*6-1, l*3-0] += dH[1, (k-1)*numNodes + l]  #       + ∂u/∂z

                            B[k*6-0, l*3-2] = dH[2, (k-1)*numNodes + l]   # γxy ← ∂v/∂x
                            B[k*6-0, l*3-1] += dH[1, (k-1)*numNodes + l]  #       + ∂u/∂y
                        end
                    elseif local_dim == 2 && rowsOfB == 4 && problem.type == :AxiSymmetric
                        @inbounds for k in 1:numNodes, l in 1:numNodes
                            # εrr, εzz, εθθ, γrz
                            B[k*4-3, l*2-1] = dH[1, (k-1)*numNodes + l]               # εrr ← ∂u_r/∂r
                            B[k*4-2, l*2-0] = dH[2, (k-1)*numNodes + l]               # εzz ← ∂u_z/∂z
                            B[k*4-1, l*2-1] = (r[k] < 1e-12) ? 0.0 : h_ref[l,k] / r[k] # εθθ ← u_r / r
                            B[k*4-0, l*2-1] = dH[2, (k-1)*numNodes + l]               # γrz ← ∂u_r/∂z + ∂u_z/∂r
                            B[k*4-0, l*2-0] += dH[1, (k-1)*numNodes + l]
                        end
                    else
                        error("solveStrain: rowsOfB=$rowsOfB, dim=$local_dim, type=$(problem.type) not supported.")
                    end

                    # globális dof indexek (nn2)
                    @inbounds for c in 1:pdim
                        nn2[c:pdim:dof_elem] = pdim .* nnet_j .- (pdim - c)
                    end

                    # elemi elmozdulás vektor(ok) (összes lépésre)
                    # disp_sub: dof_elem × nsteps
                    disp_sub = @view q.a[nn2, :]

                    # elemi eredmény: 9*numNodes × nsteps
                    e_elem = zeros(9*numNodes, nsteps)

                    # csomópontonként kivágjuk a B-ből a blokkot, és minden lépésre szorozzuk
                    @inbounds for k in 1:numNodes
                        row1 = (k-1)*rowsOfB + 1
                        row2 = k*rowsOfB
                        B1 = @view B[row1:row2, :]

                        # e0: Voigt vektor (3D:6×, 2D:3× vagy AXI:4×) minden lépésre
                        # számoljuk lépésenként, hogy ne hozzunk létre hatalmas tempókat
                        @inbounds for kk in 1:nsteps
                            e0 = B1 * @view(disp_sub[:, kk])
                            # e-tensor (9) kitöltése a típusnak megfelelően
                            if rowsOfB == 6 && local_dim == 3 && problem.type == :Solid
                                # [εxx, εyy, εzz, γyz, γxz, γxy]
                                e9 = (Float64[
                                    e0[1], e0[6]/2, e0[5]/2,
                                    e0[6]/2, e0[2], e0[4]/2,
                                    e0[5]/2, e0[4]/2, e0[3]
                                ])
                                e_elem[(k-1)*9+1 : k*9, kk] = e9
                                if DoFResults
                                    base = 9*(nnet_j[k]-1)
                                    E1[base+1:base+9, kk] .+= e9
                                end
                            elseif rowsOfB == 3 && local_dim == 2 && problem.type == :PlaneStress
                                # [εxx, εyy, γxy] → 3D tenzor komponensek, εzz = ν/(ν-1)*(εxx+εyy)
                                ezz = ν/(ν-1) * (e0[1] + e0[2])
                                e9 = (Float64[
                                    e0[1], e0[3]/2, 0.0,
                                    e0[3]/2, e0[2], 0.0,
                                    0.0,     0.0,   ezz
                                ])
                                e_elem[(k-1)*9+1 : k*9, kk] = e9
                                if DoFResults
                                    base = 9*(nnet_j[k]-1)
                                    # megjegyzés: eredetiben a DoFResults esetén γxy-t nem felezted; ott "e0[3]" ment
                                    # itt követjük az elem-eredmény konvenciót: γxy/2 → e12
                                    E1[base+1:base+9, kk] .+= e9
                                end
                            elseif rowsOfB == 3 && local_dim == 2 && problem.type == :PlaneStrain
                                # [εxx, εyy, γxy], εzz = 0
                                e9 = (Float64[
                                    e0[1], e0[3]/2, 0.0,
                                    e0[3]/2, e0[2], 0.0,
                                    0.0,     0.0,   0.0
                                ])
                                e_elem[(k-1)*9+1 : k*9, kk] = e9
                                if DoFResults
                                    base = 9*(nnet_j[k]-1)
                                    E1[base+1:base+9, kk] .+= e9
                                end
                            elseif rowsOfB == 4 && local_dim == 2 && problem.type == :AxiSymmetric
                                # [εrr, εzz, εθθ, γrz]
                                e9 = (Float64[
                                    e0[1], e0[4]/2, 0.0,
                                    e0[4]/2, e0[3], 0.0,
                                    0.0,     0.0,   e0[2]
                                ])
                                e_elem[(k-1)*9+1 : k*9, kk] = e9
                                if DoFResults
                                    base = 9*(nnet_j[k]-1)
                                    E1[base+1:base+9, kk] .+= e9
                                end
                            else
                                error("solveStrain: unexpected branch (rowsOfB=$rowsOfB, local_dim=$local_dim, type=$(problem.type)).")
                            end
                        end
                    end

                    # DoFResults-hez számlálás (hányszor érintett a csomópont)
                    if DoFResults
                        pcs[nnet_j] .+= 1
                    else
                        push!(ε, e_elem)
                        push!(numElem, elem)
                    end
                end
            end
        end
    end

    if DoFResults
        # csomóponti átlagolás
        @inbounds for l in 1:non
            cnt = pcs[l]
            if cnt > 0
                E1[9*(l-1)+1 : 9*l, :] ./= cnt
            end
        end
        return TensorField([], E1, q.t, [], nsteps, type, problem)
    else
        return TensorField(ε, [;;], q.t, numElem, nsteps, type, problem)
    end
end

function solveStrainOld(q; DoFResults=false)
    problem = q.model
    gmsh.model.setCurrent(problem.name)

    if !isa(q, VectorField) 
        error("solveStrain:argument must be a VectorField. Now it is '$(typeof(q))'.")
    end
    if q.type != :v3D && q.type != :v2D
        error("solveStrain: argument must be a displacement vector. Now it is '$(q.type)'.")
    end
    if q.A != []
        error("solveStrain: q.A != []")
    end
    nsteps = q.nsteps
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
            rowsOfB = 6
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
                    if dim == 2 && rowsOfB == 3
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*3-0, l*2-0] = B[k*3-2, l*2-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*3-0, l*2-1] = B[k*3-1, l*2-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif dim == 3 && rowsOfB == 6
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*6-5, l*3-2] = B[k*6-2, l*3-1] = B[k*6-0, l*3-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*6-4, l*3-1] = B[k*6-2, l*3-2] = B[k*6-1, l*3-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*6-3, l*3-0] = B[k*6-1, l*3-1] = B[k*6-0, l*3-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif dim == 2 && rowsOfB == 4
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*4-2, l*2-1] = r[k] < 1e-10 ? 0 : h[l, k] / r[k]
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
                        if rowsOfB == 6 && dim == 3 && problem.type == :Solid
                            B1 = B[k*6-5:k*6, 1:3*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * q.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*9+1:k*9, kk] = [e0[1], e0[4]/2, e0[6]/2,
                                        e0[4]/2, e0[2], e0[5]/2,
                                        e0[6]/2, e0[5]/2, e0[3]]
                                end
                                if DoFResults == true
                                    E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[4]/2, e0[6]/2, e0[4]/2, e0[2], e0[5]/2, e0[6]/2, e0[5]/2, e0[3]]
                                end
                            end
                        elseif rowsOfB == 3 && dim == 2 && problem.type == :PlaneStress
                            B1 = B[k*3-2:k*3, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * q.a[nn2, kk]
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
                                e0 = B1 * q.a[nn2, kk]
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
                                e0 = B1 * q.a[nn2, kk]
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
        epsilon = TensorField([], E1, q.t, [], nsteps, type, problem)
        return epsilon
    else
        epsilon = TensorField(ε, [;;], q.t, numElem, nsteps, type, problem)
        return epsilon
    end
end

"""
    solveStress(q; T=..., T₀=..., DoFResults=false)

Solves the stress field `S` from displacement vector `q`. Stress field is given
per elements, so it usually contains jumps at the boundary of elements. Details
of mesh is available in `problem`. If `DoFResults` is true, `S` is a matrix with
nodal results. In this case `showDoFResults` can be used to show the results 
(otherwise `showStressResults` or `showElementResults`).
If the `T` temperature field (and `T₀` initial temperature field if it differs from zero) is given, the
function solves also the thermal stresses.

Return: `S`

Types:
- `q`: VectorField
- `T`: ScalarField
- `T₀`: ScalarField
- `S`: TensorField
"""
function solveStress(q; T=ScalarField([],[;;],[0.0],[],0,:null,q.model), T₀=ScalarField([],reshape(zeros(q.model.non),:,1),[0],[],1,:scalar,q.model), DoFResults=false)
    problem = q.model
    gmsh.model.setCurrent(problem.name)
    
    if !isa(q, VectorField) 
        error("solveStress:argument must be a VectorField. Now it is '$(typeof(q))'.")
    end
    if q.type != :v3D && q.type != :v2D
        error("solveStress: argument must be a displacement vector. Now it is '$(q.type)'.")
    end
    if q.A != []
        error("solveStress: q.A != []")
    end

    type = :s
    nsteps = q.nsteps
    σ = []
    numElem = Int[]
    ncoord2 = zeros(3 * problem.non)
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    #if T₀ == 1im
    #    T₀ = zeros(problem.non)
    #end
    if DoFResults == true
        S1 = zeros(non * 9, nsteps)
        pcs = zeros(Int64, non * dim)
    end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        α = problem.material[ipg].α
        ακ = α * E / ν / (1 - 2ν)
        dim = 0
        if problem.dim == 3 && problem.type == :Solid
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                ν 1-ν ν 0 0 0;
                ν ν 1-ν 0 0 0;
                0 0 0 (1-2ν)/2 0 0;
                0 0 0 0 (1-2ν)/2 0;
                0 0 0 0 0 (1-2ν)/2]

            dim = 3
            rowsOfB = 6
            b = 1
            E0 = [1,1,1,0,0,0]
        elseif problem.dim == 2 && problem.type == :PlaneStress
            D = E / (1 - ν^2) * [1 ν 0;
                ν 1 0;
                0 0 (1-ν)/2]
            dim = 2
            rowsOfB = 3
            b = problem.thickness
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == :PlaneStrain
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                ν 1-ν 0;
                0 0 (1-2ν)/2]
            dim = 2
            rowsOfB = 3
            b = 1
            E0 = [1,1,0]
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                ν 1-ν ν 0;
                ν ν 1-ν 0;
                0 0 0 (1-2ν)/2]
            dim = 2
            rowsOfB = 4
            b = 1
            E0 = [1,1,1,0]
        else
            error("solveStress: dimension is $(problem.dim), problem type is $(problem.type).")
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
                #s0 = zeros(rowsOfB * numNodes)
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
                nn1 = zeros(Int, pdimT * numNodes)
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
                    if pdim == 2 && rowsOfB == 3
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                        end
                    elseif pdim == 3 && rowsOfB == 6
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                            B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif pdim == 2 && rowsOfB == 4
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*4-2, l*2-1] = r[k] < 1e-10 ? 0 : h[l, k] / r[k]
                        end
                    else
                        error("solveStress: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    for k in 1:pdimT
                        nn1[k:pdimT:pdimT*numNodes] = pdimT * nnet[j, 1:numNodes] .- (pdimT - k)
                    end
                    s = zeros(9numNodes, nsteps) # tensors have nine elements
                    for k in 1:numNodes
                        if rowsOfB == 6 && pdim == 3 && problem.type == :Solid
                            H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                            B1 = B[k*rowsOfB-5:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q.a[nn2, kk]
                                if T.type != :null
                                    s0 -= D * E0 * H1 * (T.a[nn1, kk] - T₀.a[nn1]) * α
                                end
                                if DoFResults == false
                                    s[(k-1)*9+1:k*9, kk] = [s0[1], s0[4], s0[6],
                                        s0[4], s0[2], s0[5],
                                        s0[6], s0[5], s0[3]]
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[4], s0[6], s0[4], s0[2], s0[5], s0[6], s0[5], s0[3]]
                                end
                            end
                        elseif rowsOfB == 3 && pdim == 2 && problem.type == :PlaneStress
                            H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                            B1 = B[k*rowsOfB-2:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q.a[nn2, kk]
                                if T.type != :null
                                    s0 -= D * E0 * H1 * (T.a[nn1, kk] - T₀.a[nn1]) * α
                                end
                                if DoFResults == false
                                    s[(k-1)*9+1:k*9, kk] = [s0[1], s0[3], 0,
                                        s0[3], s0[2], 0,
                                        0, 0, 0]
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[3], 0, s0[3], s0[2], 0, 0, 0, 0]
                                end
                            end
                        elseif rowsOfB == 3 && dim == 2 && problem.type == :PlaneStrain
                            H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                            B1 = B[k*rowsOfB-2:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q.a[nn2, kk]
                                if T.type != :null
                                    s0 -= D * E0 * H1 * (T.a[nn1, kk] - T₀.a[nn1]) * α
                                end
                                if DoFResults == false
                                    s[(k-1)*9+1:k*9, kk] = [s0[1], s0[3], 0,
                                        s0[3], s0[2], 0,
                                        0, 0, ν*(s0[1]+s0[2])]
                                    # PlaneStrain: σz ≠ 0
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[3], 0, s0[3], s0[2], 0, 0, 0, ν*(s0[1]+s0[2])]
                                end
                            end
                        elseif rowsOfB == 4 && dim == 2 && problem.type == :AxiSymmetric
                            H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                            B1 = B[k*4-3:k*4, 1:2*numNodes]
                            for kk in 1:nsteps
                                s0 = D * B1 * q.a[nn2, kk]
                                if T.type != :null
                                    s0 -= D * E0 * H1 * (T.a[nn1, kk] - T₀.a[nn1]) * α
                                end
                                if DoFResults == false
                                    s[(k-1)*9+1:k*9, kk] = [s0[1], s0[4], 0,
                                        s0[4], s0[3], 0,
                                        0, 0, s0[2]]
                                end
                                if DoFResults == true
                                    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[4], 0, s0[4], s0[3], 0, 0, 0, s0[2]]
                                end
                            end
                        else
                            error("solveStress: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
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
        for k in 1:9
            for l in 1:non
                S1[k + 9 * l - 9, :] ./= pcs[l]
            end
        end
    end
    if DoFResults == true
        sigma = TensorField([], S1, q.t, [], nsteps, type, problem)
        return sigma
    else
        sigma = TensorField(σ, [;;], q.t, numElem, nsteps, type, problem)
        return sigma
    end
end

function solveAxialForce(q::VectorField)
    problem = q.model
    gmsh.model.setCurrent(problem.name)
    
    if q.type != :v3D && q.type != :v2D
        error("solveStress: argument must be a displacement vector. Now it is '$(q.type)'.")
    end
    if q.A != []
        error("solveStress: q.A != []")
    end

    type = :scalar
    nsteps = q.nsteps
    σ = []
    numElem = Int[]
    ncoord2 = zeros(3 * problem.non)
    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    #if T₀ == 1im
    #    T₀ = zeros(problem.non)
    #end
    #if DoFResults == true
    #    S1 = zeros(non * 9, nsteps)
    #    pcs = zeros(Int64, non * dim)
    #end

    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        A = problem.material[ipg].A
        #α = problem.material[ipg].α
        #ακ = α * E / ν / (1 - 2ν)
        dim = 0

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            if edim ≠ 1
                error("stiffnessMatrixTruss: not 1D elements.")
            end
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(edim, -1, true, false)
            ncoord2[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
            ncoord2[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
            for i in 1:length(elemTypes)
                et = elemTypes[i]
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                #s0 = zeros(rowsOfB * numNodes)
                nodeCoord = zeros(numNodes * 3)
                for k in 1:dim, j = 1:numNodes
                    nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
                end
                #comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
                #∇h = reshape(dfun, :, numNodes)
                #comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
                #h = reshape(fun, :, numNodes)
                nnet = zeros(Int, length(elemTags[i]), numNodes)
                #pdimT = 1
                #H = zeros(pdimT * numNodes, pdimT * numNodes)
                #for j in 1:numNodes
                #    for k in 1:numNodes
                #        for l in 1:pdimT
                #            H[j*pdimT-(pdimT-l), k*pdimT-(pdimT-l)] = h[k, j]
                #        end
                #    end
                #end
                #invJac = zeros(3, 3numNodes)
                #∂h = zeros(3, numNodes * numNodes)
                #B = zeros(rowsOfB * numNodes, pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                #nn1 = zeros(Int, pdimT * numNodes)
                #r = zeros(numNodes)
                B1 = [-1 1]
                for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]
                    for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
                    end
                    #jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    #Jac = reshape(jac, 3, :)
                    #for k in 1:numNodes
                    #    invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
                    #    r[k] = h[:, k]' * ncoord2[nnet[j, :] * 3 .- 2]
                    #end
                    #∂h .*= 0
                    #for k in 1:numNodes, l in 1:numNodes
                    #    ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
                    #end
                    #B .*= 0
                    #if pdim == 2 && rowsOfB == 3
                    #    for k in 1:numNodes, l in 1:numNodes
                    #        B[k*rowsOfB-0, l*pdim-0] = B[k*rowsOfB-2, l*pdim-1] = ∂h[1, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-0, l*pdim-1] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                    #    end
                    #elseif pdim == 3 && rowsOfB == 6
                    #    for k in 1:numNodes, l in 1:numNodes
                    #        B[k*rowsOfB-5, l*pdim-2] = B[k*rowsOfB-2, l*pdim-1] = B[k*rowsOfB-0, l*pdim-0] = ∂h[1, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-4, l*pdim-1] = B[k*rowsOfB-2, l*pdim-2] = B[k*rowsOfB-1, l*pdim-0] = ∂h[2, (k-1)*numNodes+l]
                    #        B[k*rowsOfB-3, l*pdim-0] = B[k*rowsOfB-1, l*pdim-1] = B[k*rowsOfB-0, l*pdim-2] = ∂h[3, (k-1)*numNodes+l]
                    #    end
                    #elseif pdim == 2 && rowsOfB == 4
                    #    for k in 1:numNodes, l in 1:numNodes
                    #        B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes+l]
                    #        B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes+l]
                    #        B[k*4-2, l*2-1] = r[k] < 1e-10 ? 0 : h[l, k] / r[k]
                    #    end
                    #else
                    #    error("solveStress: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    #end
                    push!(numElem, elem)
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    #for k in 1:pdimT
                    #    nn1[k:pdimT:pdimT*numNodes] = pdimT * nnet[j, 1:numNodes] .- (pdimT - k)
                    #end
                    s = zeros(numNodes, nsteps) # tensors have nine elements
                    x1 = ncoord2[nnet[j, 1] * 3 .- 2]
                    y1 = ncoord2[nnet[j, 1] * 3 .- 1]
                    z1 = ncoord2[nnet[j, 1] * 3 .- 0]
                    x2 = ncoord2[nnet[j, 2] * 3 .- 2]
                    y2 = ncoord2[nnet[j, 2] * 3 .- 1]
                    z2 = ncoord2[nnet[j, 2] * 3 .- 0]
                    L = √((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
                    c1 = (x2 - x1) / L
                    c2 = (y2 - y1) / L
                    c3 = (z2 - z1) / L
                    T = [c1 c2 c3 0 0 0; 0 0 0 c1 c2 c3]
                    k1 = A * E / L
                    for k in 1:numNodes
                        #if rowsOfB == 6 && pdim == 3 && problem.type == :Solid
                            #H1 = H[k*pdimT-(pdimT-1):k*pdimT, 1:pdimT*numNodes]
                            #B1 = B[k*rowsOfB-5:k*rowsOfB, 1:pdim*numNodes]
                            for kk in 1:nsteps
                                s0 = k1 * B1 * T * q.a[nn2, kk]
                                #if T.type != :null
                                #    s0 -= D * E0 * H1 * (T.a[nn1, kk] - T₀.a[nn1]) * α
                                #end
                                #if DoFResults == false
                                    s[k, kk] = s0[1]
                                #end
                                #if DoFResults == true
                                #    S1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [s0[1], s0[4], s0[6], s0[4], s0[2], s0[5], s0[6], s0[5], s0[3]]
                                #end
                            end
                    end
                    #if DoFResults == true
                    #    pcs[nnet[j,1:numNodes]] .+= 1
                    #end
                    #if DoFResults == false
                        push!(σ, s)
                    #end
                end
            end
        end
    end
    #if DoFResults == true
    #    for k in 1:9
    #        for l in 1:non
    #            S1[k + 9 * l - 9, :] ./= pcs[l]
    #        end
    #    end
    #end
    #if DoFResults == true
    #    sigma = TensorField([], S1, q.t, [], nsteps, type, problem)
    #    return sigma
    #else
        sigma = ScalarField(σ, [;;], q.t, numElem, nsteps, type, problem)
        return sigma
    #end
end

"""
    solveEigenModes(K, M; n=6, fₘᵢₙ=1.01)

Solves the eigen frequencies and mode shapes of a problem given by stiffness
matrix `K` and the mass matrix `M`. `n` is the number of eigenfrequencies to solve,
and solves the eigenfrequencies greater than `fₘᵢₙ`. Returns the struct of eigenfrequencies
and eigen modes. Results can be presented by `showModalResults` function.

Return: `modes`

Types:
- `K`: SystemMatrix
- `M`: SystemMatrix
- `n`: Int64
- `fₘᵢₙ`: Float64
- `modes`: Eigen 
"""
function solveEigenModes(K, M; n=6, fₘᵢₙ=0.01)
    if K.model != M.model
        error("solveEigenModes: K and M does not belong to the same model.")
    end
    problem = K.model
    ωₘᵢₙ² = (2π * fₘᵢₙ)^2
    ω², ϕ = Arpack.eigs(K.A, M.A, nev=n, which=:LR, sigma=ωₘᵢₙ², maxiter=10000)
    #if real(ω²[1]) > 0.999 && real(ω²[1]) < 1.001
    #    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    #end
    #err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1]) / norm(K * ϕ[:,1])
    #if err > 1e-3 # || true
    #    @warn("The error in the calculation of the smallest eigenvalue is too large: $err")
    #end
    f = sqrt.(abs.(real(ω²))) / 2π
    ϕ1 = real(ϕ)
    return Eigen(f, ϕ1, problem)
end

"""
    solveBucklingModes(K, Knl; n=6)

Solves the critical force multipliers and buckling mode shapes of a problem given by stiffness
matrix `K` and the nonlinear stiffness matrix `Knl`. `n` is the number of buckling modes to solve.
Returns the struct of critical forces and buckling modes. Results can be presented by `showBucklingResults` function.

Return: `modes`

Types:
- `K`: SystemMatrix
- `Knl`: SystemMatrix
- `n`: Int64
- `modes`: Eigen 
"""
function solveBucklingModes(K, Knl; n=6)
    if K.model != Knl.model
        error("solveBucklingModes: K and Knl does not belong to the same model.")
    end
    problem = K.model
    λ, ϕ = Arpack.eigs(K.A, -Knl.A, nev=n, which=:LR, sigma=0, maxiter=1000)
    #if real(ω²[1]) > 0.999 && real(ω²[1]) < 1.001
    #    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    #end
    #err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1]) / norm(K * ϕ[:,1])
    #if err > 1e-3 # || true
    #    error("The error in the calculation of the smallest eigenvalue is too large: $err")
    #end
    f = abs.(real(λ))
    ϕ1 = real(ϕ)
    return Eigen(f, ϕ1, problem)
end

"""
    solveModalAnalysis(problem; constraints=[]; loads=[], n=6)

Solves the first `n` eigenfrequencies and the corresponding 
mode shapes for the `problem`, when `loads` and 
`constraints` are applied. `loads` and `contraints` are optional. 
Result can be presented by `showModalResults` function. 
`loads` and `constraints` can be defined by `load` and `displacementConstraint` functions,
respectively. If `loads` are given, it solves the eigenfrequencies of a prestressed structure.

Return: `modes`

Types:
- `problem`: Problem
- `loads`: Vector{tuples}
- `constraints`: Vector{tuples}
- `n`: Int64
- `modes`: Eigen
"""
function solveModalAnalysis(problem; constraints=[], loads=[], n=6, fₘᵢₙ=0.1)
    if !isa(loads, Vector)
        error("solveModalAnalysis: loads are not arranged in a vector. Put them in [...]")
    end
    if !isa(constraints, Vector)
        error("solveModalAnalysis: constraints are not arranged in a vector. Put them in [...]")
    end
    dof = problem.pdim * problem.non
    K = stiffnessMatrix(problem)
    M = massMatrix(problem)
    if length(constraints) == 0
        return solveEigenModes(K, M, n=n, fₘᵢₙ=fₘᵢₙ)
    elseif length(loads) == 0
        fdof = freeDoFs(problem, constraints)
        cdof = constrainedDoFs(problem, constraints)
        K = SystemMatrix(K.A[fdof, fdof], K.model)
        M = SystemMatrix(M.A[fdof, fdof], M.model)
        #f = zeros(dof)
        #applyBoundaryConditions!(problem, K, M, f, constraints, fix=fₘᵢₙ<1 ? fₘᵢₙ/10 : 0.1)
        mod = solveEigenModes(K, M, n=n, fₘᵢₙ=fₘᵢₙ)
        nn = length(mod.f)
        ϕ1 = zeros(dof, nn)
        ϕ1[fdof,:] .= mod.ϕ
        applyBoundaryConditions!(ϕ1, constraints)
        return Eigen(mod.f, ϕ1, problem)
    else
        fdof = freeDoFs(problem, constraints)
        cdof = constrainedDoFs(problem, constraints)
        ϕ1 = zeros(dof, n)
        f = loadVector(problem, loads)
        applyBoundaryConditions!(K, f, constraints)
        q = solveDisplacement(K, f)

        err = 1
        count = 0
        while err > 1e-3 && count < 10
            count += 1
            q0 = copy(q)
            Knl = nonLinearStiffnessMatrix(q)
            applyBoundaryConditions!(Knl, f, constraints)
            q = solveDisplacement(K + Knl, f)
            err = sum(abs, q.a - q0.a) / (sum(abs, q0.a) == 0 ? 1 : sum(abs, q0.a))
        end
        if count == 10
            @warn("solveModalAnalysis: number of iterations is $count.")
        end
        Knl = nonLinearStiffnessMatrix(q)

        #applyBoundaryConditions!(problem, K, M, Knl, f, constraints, fix=fₘᵢₙ<1 ? fₘᵢₙ/10 : 0.1)
        mod = solveEigenModes(SystemMatrix((K.A + Knl.A)[fdof,fdof], K.model), SystemMatrix(M.A[fdof,fdof], M.model), n=n, fₘᵢₙ=fₘᵢₙ)
        ϕ1[fdof,:] .= mod.ϕ
        applyBoundaryConditions!(VectorField([], ϕ1, [0.0], [], 1, q.type, q.model), constraints)
        return Eigen(mod.f, ϕ1, problem)
    end
end

"""
    solveBuckling(problem, loads, constraints; n=6)

Solves the multipliers for the first `n` critical forces and the corresponding 
buckling shapes for the instability of the `problem`, when `loads` and 
`constraints` are applied. Result can be presented by `showBucklingResults`
function. `loads` and `constraints` can be defined by `load` and `displacementConstraint` functions,
respectively.

Return: `buckling`

Types:
- `problem`: Problem
- `loads`: Vector{tuples}
- `constraints`: Vector{tuples}
- `n`: Int64
- `buckling`: Eigen 
"""
function solveBuckling(problem, loads, constraints; n=6)
    f = loadVector(problem, loads)
    K = stiffnessMatrix(problem)
    applyBoundaryConditions!(K, f, constraints)
    q = solveDisplacement(K, f)

    err = 1
    count = 0
    while err > 1e-3 && count < 10
        count += 1
        q0 = copy(q)
        Knl = nonLinearStiffnessMatrix(q)
        applyBoundaryConditions!(Knl, f, constraints)
        q = solveDisplacement(K + Knl, f)
        err = sum(abs, q.a - q0.a) / (sum(abs, q0.a) == 0 ? 1 : sum(abs, q0.a))
    end
    if count == 10
        @warn("solveBuckling: number of iterations is $count.")
    end
    Knl = nonLinearStiffnessMatrix(q)

    applyBoundaryConditions!(Knl, f, constraints)
    return solveBucklingModes(K, Knl, n=n)
end

"""
    initialDisplacement(problem, name; ux=..., uy=..., uz=...)

Sets the displacement values `ux`, `uy` and `uz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Returns the initial
displacement vector `u0`.

Return: u0

Types:
- `problem`: Problem
- `name`: String 
- `u0`: VectorField
- `ux`: Float64 
- `uy`: Float64 
- `uz`: Float64 
"""
function initialDisplacement(problem, name; ux=1im, uy=1im, uz=1im)
    pdim = problem.pdim
    dim = problem.dim
    u0 = zeros(problem.non * problem.dim)
    phg = getTagForPhysicalName(name)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
    if ux != 1im
        for i in 1:length(nodeTags)
            u0[nodeTags[i]*dim-(dim-1)] = ux
        end
    end
    if uy != 1im
        for i in 1:length(nodeTags)
            u0[nodeTags[i]*dim-(dim-2)] = uy
        end
    end
    if dim == 3 && uz != 1im
        for i in 1:length(nodeTags)
            u0[nodeTags[i]*dim] = uz
        end
    end
    if pdim == 3
        type = :v3D #Symbol(String(type) * "3D")
    elseif pdim == 2
        type = :v2D #Symbol(String(type) * "2D")
    else
        error("initialDisplacement: wrong pdim=$pdim")
    end
    return VectorField([], reshape(u0, :,1), [0], [], 1, type, problem)
end

"""
    initialDisplacement!(name, u0; ux=..., uy=..., uz=...)

Changes the displacement values to `ux`, `uy` and `uz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Original values are in
displacement vector `u0`.

Return: u0

Types:
- `name`: String 
- `ux`: Float64 
- `uy`: Float64 
- `uz`: Float64 
- `u0`: VectorField
"""
function initialDisplacement!(name, u0; ux=1im, uy=1im, uz=1im, type=:u)
    problem = u0.model
    pdim = problem.pdim
    dim = problem.dim
    phg = getTagForPhysicalName(name)
    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
    if ux != 1im
        for i in 1:length(nodeTags)
            u0.a[nodeTags[i]*dim-(dim-1)] = ux
        end
    end
    if uy != 1im
        for i in 1:length(nodeTags)
            u0.a[nodeTags[i]*dim-(dim-2)] = uy
        end
    end
    if dim == 3 && uz != 1im
        for i in 1:length(nodeTags)
            u0.a[nodeTags[i]*dim] = uz
        end
    end
end

"""
    initialVelocity(problem, name; vx=..., vy=..., vz=...)

Sets the velocity values `vx`, `vy` and `vz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Returns the initial
velocity vector `v0`.

Return: v0

Types:
- `problem`: Problem
- `name`: String 
- `vx`: Float64 
- `vy`: Float64 
- `vz`: Float64 
- `v0`: VectorField
"""
function initialVelocity(problem, name; vx=1im, vy=1im, vz=1im)
    return initialDisplacement(problem, name, ux=vx, uy=vy, uz=vz, type=:v)
end

"""
    initialVelocity!(name, v0; vx=..., vy=..., vz=...)

Changes the velocity values `vx`, `vy` and `vz` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Original values are in
velocity vector `v0`.

Returns: nothing

Types:
- `name`: String 
- `v0`: VectorField
- `vx`: Float64 
- `vy`: Float64 
- `vz`: Float64 
"""
function initialVelocity!(name, v0; vx=1im, vy=1im, vz=1im)
    initialDisplacement!(name, v0, ux=vx, uy=vy, uz=vz)
end

"""
    nodalForce!(name, f0; fx=..., fy=..., fz=...)

Changes the force values `fx`, `fy` and `fz` (depending on the dimension of
the problem) at nodes belonging to physical group `name`. Original values are in
load vector `f0`.

Returns: nothing

Types:
- `name`: String 
- `f0`: VectorField
- `fx`: Float64 
- `fy`: Float64 
- `fz`: Float64 
"""
function nodalForce!(f0; fx=1im, fy=1im, fz=1im)
    initialDisplacement!(f0, ux=fx, uy=fy, uz=fz)
end

"""
    nodalAcceleration!(name, a0; ax=..., ay=..., az=...)

Changes the acceleration values `ax`, `ay` and `az` (depending on the dimension of
the `problem`) at nodes belonging to physical group `name`. Original values are in
acceleration vector `a0`.

Returns: nothing

Types:
- `name`: String 
- `a0`: VectorField
- `ax`: Float64
- `ay`: Float64
- `az`: Float64
"""
function nodalAcceleration!(a0; ax=1im, ay=1im, az=1im)
    initialDisplacement!(a0, ux=ax, uy=ay, uz=az)
end

"""
    largestPeriodTime(K, M)

Solves the largest period of time for a dynamic problem given by stiffness
matrix `K` and the mass matrix `M`.

Return: `Δt`

Types:
- `K`: SystemMatrix
- `M`: SystemMatrix
- `Δt`: Float64 
"""
function largestPeriodTime(K, M)
    ω², ϕ = Arpack.eigs(K.A, M.A, nev=1, which=:LR, sigma=0.01, maxiter=10000)
    if real(ω²[1]) > 0.999 && real(ω²[1]) < 1.001
        ω², ϕ = Arpack.eigs(K.A, M.A, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    end
    err = norm(K.A * ϕ[:,1] - ω²[1] * M.A * ϕ[:,1]) / norm(K.A * ϕ[:,1])
    if err > 1e-3 # || true
        @warn("The error in the calculation of the smallest eigenvalue is too large: $err")
    end
    Δt = 2π / √(abs(real(ω²[1])))
    return Δt
end

"""
    smallestPeriodTime(K, M)

Solves the smallest period of time for a dynamic problem given by stiffness
matrix `K` and the mass matrix `M`.

Return: `Δt`

Types:
- `K`: SystemMatrix
- `M`: SystemMatrix
- `Δt`: Float64 
"""
function smallestPeriodTime(K, M)
    ω², ϕ = Arpack.eigs(K.A, M.A, nev=1, which=:LM, maxiter=100)
    
    err = norm(K.A * ϕ[:,1] - ω²[1] * M.A * ϕ[:,1]) / norm(K.A * ϕ[:,1])
    if err > 1e-3 # || true
        @warn("The error in the calculation of the largest eigenvalue is too large: $err")
    end
    Δt = 2π / √(abs(real(ω²[1])))
    return Δt
end

"""
    smallestEigenValue(K, M)

Solves the largest eigenvalue for a transient problem given by stiffness (heat conduction)
matrix `K` and the mass (heat capacity) matrix `M` (`C`).

Return: `λₘₐₓ`

Types:
- `K`: SystemMatrix
- `M`: SystemMatrix
- `λₘₐₓ`: Float64 
"""
function smallestEigenValue(K, C)
    λ, ϕ = Arpack.eigs(K.A, C.A, nev=1, which=:LR, sigma=0.0001, maxiter=10000)
    if real(λ[1]) > 0.999 && real(λ[1]) < 1.001
        λ, ϕ = Arpack.eigs(K.A, C.A, nev=1, which=:LR, sigma=1.01, maxiter=10000)
    end
    err = norm(K.A * ϕ[:,1] - λ[1] * C.A * ϕ[:,1]) / norm(K.A * ϕ[:,1])
    if err > 1e-3 # || true
        @warn("The error in the calculation of the largest eigenvalue is too large: $err")
    end
    λₘₐₓ = abs(real(λ[1]))
    return λₘₐₓ
end

"""
    largestEigenValue(K, M)

Solves the smallest eigenvalue for a transient problem given by stiffness (heat conduction)
matrix `K` and the mass (heat capacity) matrix `M` (`C`).

Return: `λₘᵢₙ`

Types:
- `K`: SystemMatrix
- `M`: SystemMatrix
- `λₘᵢₙ`: Float64 
"""
function largestEigenValue(K, C)
    λ, ϕ = Arpack.eigs(K.A, C.A, nev=1, which=:LM)

    err = norm(K.A * ϕ[:,1] - λ[1] * C.A * ϕ[:,1]) / norm(K.A * ϕ[:,1])
    if err > 1e-3 # || true
        @warn("The error in the calculation of the smallest eigenvalue is too large: $err")
    end
    λₘᵢₙ = abs(real(λ[1]))
    return λₘᵢₙ
end

"""
    CDM(K, M, C, f, u0, v0, T, Δt)

Solves a transient dynamic problem using central difference method (CDM) (explicit).
`K` is the stiffness Matrix, `M` is the mass matrix, `C` is the damping matrix,
`f` is the load vector, `u0` is the initial displacement, `v0` is the initial
velocity, `T` is the upper bound of the time intervall (lower bound is zero)
and `Δt` is the time step size. Returns the displacement vectors and velocity
vectors in each time step arranged in the columns of the two matrices `u` and `v`
and a vector `t` of the time instants used.

The critical (largest allowed) time step is `Δtₘₐₓ = Tₘᵢₙ / π * (√(1 + ξₘₐₓ^2) - ξₘₐₓ)`
where `Tₘᵢₙ` is the time period of the largest eigenfrequency and `ξₘₐₓ` is the largest
modal damping.

Return: `u`, `v`

Types:
- `K`: SystemMatrix
- `M`: SystemMatrix
- `C`: SystemMatrix
- `f`: VectorField
- `u0`: VectorField
- `v0`: VectorField
- `T`: Float64
- `Δt`: Float64 
- `u`: VectorField
- `v`: VectorField
"""
function CDM(K, M, C, f, uu0, vv0, T, Δt)
    if K.model != M.model || M.model != C.model || C.model != f.model
        error("CDM: K, M, C and f does not belong to the same model.")
    end
    invM = spdiagm(1 ./ diag(M.A))
    nsteps = ceil(Int64, T / Δt)
    dof, dof = size(K.A)

    u = zeros(dof, nsteps)
    v = zeros(dof, nsteps)
    t = zeros(nsteps)
    #kene = zeros(nsteps)
    #sene = zeros(nsteps)
    #diss = zeros(nsteps)

    a0 = M.A \ (f.a - K.A * uu0.a - C.A * vv0.a)
    u00 = uu0.a - vv0.a * Δt + a0 * Δt^2 / 2

    u[:, 1] = uu0.a
    v[:, 1] = vv0.a
    t[1] = 0
    #kene[1] = dot(v0' * M, v0) / 2
    #sene[1] = dot(u0' * K, u0) / 2
    u0 = uu0.a[:,1]

    for i in 2:nsteps
        u1 = 2.0 * u0 - u00 + Δt * Δt * invM * (f.a - K.A * u0) - Δt * invM * (C.A * (u0 - u00))
        u[:, i] = u1
        v1 = (u1 - u0) / Δt
        v[:, i] = v1
        t[i] = t[i-1] + Δt
        #kene[i] = dot(v1' * M, v1) / 2
        #sene[i] = dot(u1' * K, u1) / 2
        #diss[i] = dot(v1' * C, v1)
        u00 = u0
        u0 = u1
    end
    return VectorField([], u, t, [], length(t), uu0.type, f.model), VectorField([], v, t, [], length(t), vv0.type, f.model)
end

function CDM(K, M, f, u0, v0, T, Δt)
    C = K * 0.0
    dropzeros!(C.A)
    return CDM(K, M, C, f, u0, v0, T, Δt)
end

"""
    HHT(K, M, f, u0, v0, T, Δt; α=..., δ=..., γ=..., β=...)

Solves a transient dynamic problem using HHT-α method[^1] (implicit).
`K` is the stiffness Matrix, `M` is the mass matrix, `f` is the load vector, 
`u0` is the initial displacement, `v0` is the initial velocity, `T` is the 
upper bound of the time intervall (lower bound is zero) and `Δt` is the time 
step size. Returns the displacement vectors and velocity vectors in each time 
step arranged in the columns of the two matrices `u` and `v` and a vector `t` 
of the time instants used. For the meaning of `α`, `β` and `γ` see [^1]. If
`δ` is given, γ=0.5+δ and β=0.25⋅(0.5+γ)².

[^1]: Hilber, Hans M., Thomas JR Hughes, and Robert L. Taylor. *Improved 
    numerical dissipation for time integration algorithms in structural 
    dynamics*. Earthquake Engineering & Structural Dynamics 5.3 (1977): 283-292.

Return: `u`, `v`

Types:
- `K`: SystemMatrix
- `M`: SystemMatrix
- `f`: VectorField
- `u0`: VectorField
- `v0`: VectorField
- `T`: Float64
- `Δt`: Float64 
- `α`: Float64
- `β`: Float64
- `γ`: Float64
- `δ`: Float64
- `u`: VectorField
- `v`: VectorField
"""
function HHT(K, M, f, uu0, vv0, T, Δt; α=0.0, δ=0.0, γ=0.5 + δ, β=0.25 * (0.5 + γ)^2)
    nsteps = ceil(Int64, T / Δt)
    dof, dof = size(K.A)

    u = zeros(dof, nsteps)
    v = zeros(dof, nsteps)
    t = zeros(nsteps)
    kene = zeros(nsteps)
    sene = zeros(nsteps)
    diss = zeros(nsteps)

    dt = Δt
    dtdt = dt * dt

    c0 = 1.0 / (β * dtdt)
    c1 = γ / (β * dt)
    c2 = 1.0 / (β * dt)
    c3 = 0.5 / β - 1.0
    c4 = γ / β - 1.0
    c5 = dt / 2.0 * (γ / β - 2.0)
    c6 = dt * (1.0 - γ)
    c7 = γ * dt

    a0 = M.A \ (f.a - K.A * uu0.a)

    u[:, 1] = uu0.a
    v[:, 1] = vv0.a
    t[1] = 0
    #kene[1] = dot(v0' * M, v0) / 2
    #sene[1] = dot(u0' * K, u0) / 2
    
    A = (α + 1) * K.A + M.A * c0
    AA = lu(A)

    u0 = uu0.a[:,1]
    v0 = vv0.a[:,1]
    for i in 2:nsteps
        b = f.a + M.A * (u0 * c0 + v0 * c2 + a0 * c3) + α * K.A * u0
        u1 = AA \ b
        u[:, i] = u1
        a1 = (u1 - u0) * c0 - v0 * c2 - a0 * c3
        v1 = v0 + a0 * c6 + a1 * c7
        v[:, i] = v1
        t[i] = t[i-1] + Δt
        #kene[i] = dot(v1' * M, v1) / 2
        #sene[i] = dot(u1' * K, u1) / 2
        #diss[i] = dot(v1' * C, v1)
        u0 = u1
        v0 = v1
        a0 = a1
    end
    return VectorField([], u, t, [], length(t), uu0.type, f.model), VectorField([], v, t, [], length(t), vv0.type, f.model)
end

"""
    CDMaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=..., ξ=..., β=..., show_β=..., show_ξ=...)

Gives some functions (graphs) for accuracy analysis of the CDM method. 
`ωₘᵢₙ` and `ωₘₐₓ` are the square root of smallest and largest eigenvalues of the
**Kϕ**=ω²**Mϕ** eigenvalue problem, `Δt` is the time step size. `type` is one of the
following values:
- `:SR`: spectral radius
- `:PDR`: physical damping ratio
- `:ADR`: algorithmic damping ratio
- `:PE`: period error
For details see [^3]. 
`n` is the number of points in the graph. The damping matrix is assembled in the 
following ways: **C**=α**M**+β**K** or **C**=α**M**+β₁**K**+β₂**KM⁻¹K**+β₃**KM⁻¹KM⁻¹K**+⋅⋅⋅. 
The latter corresponds to the damping characteristic characterized by a power series consisting of powers
of the natural frequencies with odd exponents. ξᵢ (`ξ` in the argument list) are the values ​​of the 
individual members of the series corresponding to the ωₘₐₓ value. βᵢ (`β` in the argument list) are the 
coefficients of the series. (see [^4]) Either `ξ` or `β` must be specified. `ξ` or `β` are scalars or 
vectors. If `show_β` or `show_ξ` is `true`, the corresponding `β` or `ξ` values will be 
sent to the output.
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

[^4]: Serfőző, D., Pere, B.: *An effective reduction method with Caughey damping for 
    spurious oscillations in dynamic problems*, Meccanica, <https://doi.org/10.1007/s11012-025-02036-9>

Return: `xy`

Types:
- `ωₘᵢₙ`: Float64
- `ωₘₐₓ`: Float64
- `Δt`: Float64 
- `n`: Int64
- `α`: Float64
- `β`: Float64 of Vector{Float64}
- `ξ`: Float64 of Vector{Float64}
- `show_β`: Boolean
- `show_ξ`: Boolean
- `xy`: Tuple{Vector{Float64},Vector{Float64}}
"""
function CDMaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, ξ=0.01, β=[2ξ[i]/(ωₘₐₓ)^(2i-1) for i in 1:length(ξ)], show_β=false, show_ξ=false)
    if show_β == true
        println("β = $β")
    end
    if show_ξ == true
        ξ = [β[i] / 2 * ωₘₐₓ^((2i-1)) for i in 1:length(β)]
        println("ξ = $ξ")
    end
    #Tₘᵢₙ /= √(1-ξₘₐₓ^2)
    x = zeros(n)
    y = zeros(n)
    ω = range(ωₘᵢₙ, length=n, stop=ωₘₐₓ)
    for i ∈ 1:n
        ξ = α/ω[i]
        for j in 1:length(β)
            ξ += β[j] / 2 * ω[i]^(2j-1)
        end
        Ω = Δt * ω[i]
        A = [2-2ξ*Ω-Ω^2 2ξ*Ω-1
            1 0]

        eig = eigen(A)
        ρ, idx = findmax(abs, eig.values)
        λ = eig.values[idx]
        σ = real(λ)
        ε = imag(λ)
        if type == :SR
            x[i] = log((ω[i] / 2π) * Δt)
            y[i] = ρ
        elseif type == :ADR
            x[i] = (ω[i] / 2π) * Δt
            Ω0 = √((log(ρ))^2 + (atan(ε,σ))^2 / 4)
            y[i] = -log(ρ) / 2Ω0
        elseif type == :PDR
            x[i] = (ω[i] / 2π) * Δt
            for j in 1:length(β)
                y[i] += β[j] / 2 * (2π * x[i] / Δt) ^ (2j-1)
            end
        elseif type == :PE
            x[i] = (ω[i] / 2π) * Δt
            Ω0 = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = 1 - Ω0/(Δt*ω[i])
        else
            str1 = "CDMaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\n:SR - spectral radius\n"
            str4 = ":PDR - physical damping ratio\n"
            str5 = ":ADR - algorithmic damping ratio\n"
            str6 = ":PE - period error\n"
            str7 = "\nFor details see Serfőző, D., Pere, B.: A method to accurately define arbitrary\n"
            str8 = "algorithmic damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).\n"
            str9 = "https://doi.org/10.1007/s00419-023-02454-9\n"
            error(str1*str2*str3*str4*str5*str6*str7*str8*str9)
        end
    end
    return x, y
end

"""
    HHTaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, δ=0.0, γ=0.5 + δ, β=0.25 * (0.5 + γ)^2)

Gives some functions (graphs) for accuracy analysis of the HHT-α method[^1]. 
`ωₘᵢₙ` and `ωₘₐₓ` are the square root of smallest and largest eigenvalues of the
**Kϕ**=ω²**Mϕ** eigenvalue problem, `Δt` is the time step size. `type` is one of the
following values:
- `:SR`: spectral radius
- `:ADR`: algorithmic damping ratio
- `:PE`: period error
For details see [^2] and [^3]. 
`n` is the number of points in the graph. For the meaning of `α`, `β` and `γ`
see [^1]. If `δ` is given, γ=0.5+δ and β=0.25⋅(0.5+γ)².
Returns a tuple of x and y values of the graph. (Can be plotted with `plot(xy)`)

[^2]: Belytschko, Ted, and Thomas JR, Hughes: *Computational methods for 
    transient analysis*, North-Holland, (1983).

[^3]: Serfőző, D., Pere, B.: *A method to accurately define arbitrary algorithmic
    damping character as viscous damping*. Arch Appl Mech 93, 3581–3595 (2023).
    <https://doi.org/10.1007/s00419-023-02454-9>

Return: `xy`

Types:
- `ωₘᵢₙ`: Float64
- `ωₘₐₓ`: Float64
- `Δt`: Float64 
- `n`: Int64
- `α`: Float64
- `β`: Float64
- `γ`: Float64
- `δ`: Float64
- `xy`: Tuple{Vector{Float64},Vector{Float64}}
"""
function HHTaccuracyAnalysis(ωₘᵢₙ, ωₘₐₓ, Δt, type; n=100, α=0.0, δ=0.0, γ=0.5 + δ, β=0.25 * (0.5 + γ)^2)
    x = zeros(n)
    y = similar(x)
    invT = range(ωₘᵢₙ/2π, length=n, stop=ωₘₐₓ/2π)
    for i ∈ 1:n
        ω = 2π * invT[i]
        A1 = [1 0 -Δt^2*β
            0 1 -Δt*γ
            (1+α)*ω^2 0 1]
        A2 = [1 Δt Δt^2*(0.5-β)
            0 1 Δt*(1-γ)
            α*ω^2 0 0]

        A = A1 \ A2

        eig = eigen(A)
        ρ, idx = findmax(abs, eig.values)
        λ = eig.values[idx]
        σ = real(λ)
        ε = imag(λ)
        if type == :SR
            x[i] = log(invT[i] * Δt)
            y[i] = ρ
        elseif type == :ADR
            x[i] = invT[i] * Δt
            Ω = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = -log(ρ) / 2Ω
            #y[i] = -log(ρ) / atan(ε, σ)
        elseif type == :PE
            x[i] = invT[i] * Δt
            Ω = √(log(ρ)^2 / 4 +atan(ε,σ)^2)
            y[i] = 1 - Ω/(2π*Δt*invT[i])
        else
            str1 = "HHTaccuracyAnalysis: wrong analysis type: $type\n"
            str2 = "Possibilities:\n"
            str3 = "\n:SR - spectral radius\n"
            str5 = ":ADR - algorithmic damping ratio\n"
            str6 = ":PE - period error\n"
            str7 = "\nFor details see Serfőző, D., Pere, B.: A method to accurately define arbitrary\n"
            str8 = "algorithmic damping character as viscous damping. Arch Appl Mech 93, 3581–3595 (2023).\n"
            str9 = "https://doi.org/10.1007/s00419-023-02454-9\n"
            error(str1*str2*str3*str5*str6*str7*str8*str9)
        end
    end
    return x, y
end
