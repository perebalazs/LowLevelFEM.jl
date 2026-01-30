export nodePositionVector, ∇, curl, rot, div, grad
export tangentMatrixConstitutive, tangentMatrixInitialStress
export equivalentNodalForce, nonFollowerLoadVector
export applyDeformationBoundaryConditions!, suppressDeformationAtBoundaries!, suppressDeformationAtBoundaries
export solveDeformation, showDeformationResults
export grad_xy
export materialTangentMatrix, initialStressMatrix, externalTangentFollower, internalForceVector
export IIPiolaKirchhoff

"""
    nodePositionVector(problem)

Returns the position vectors of all mesh nodes as a `VectorField` (initial configuration).
Returning vector is always a 3D vector.

Returns: `R`

Types:
- `problem`: Problem
- `R`: VectorField
"""
function nodePositionVector(problem)
    nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes()
    r = zeros(length(nodeTags) * 3)
    r[nodeTags*3 .- 2] = ncoord[1:3:length(ncoord)]
    r[nodeTags*3 .- 1] = ncoord[2:3:length(ncoord)]
    r[nodeTags*3 .- 0] = ncoord[3:3:length(ncoord)]
    return VectorField([], reshape(r, :,1), [0], [], 1, :v3D, problem)
end

#=
function ∇_old(rr::Union{VectorField, ScalarField, TensorField}; nabla=:grad)
    problem = rr.model
    gmsh.model.setCurrent(problem.name)
    if rr.a == [;;]
        r = elementsToNodes(rr)
    else
        r = rr
    end
    DoFResults=false
    
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

    for ipg in 0:length(problem.material)
        phName = ""
        if ipg == 0
            if rr.model.geometry.nameVolume ≠ ""
                phName = rr.model.geometry.nameVolume
            else
                continue
            end
        else
            phName = problem.material[ipg].phName
        end
        #ν = problem.material[ipg].ν
        dim = 0
        if problem.dim == 3 && r isa VectorField && nabla == :grad
            dim = 3
            pdim = 3
            rowsOfB = 9
            b = 1
            sz = 9
        elseif problem.dim == 3 && r isa VectorField && nabla == :div
            dim = 3
            pdim = 3
            rowsOfB = 9
            b = 1
            sz = 1
        elseif problem.dim == 3 && r isa VectorField && nabla == :curl
            dim = 3
            pdim = 3
            rowsOfB = 9
            b = 1
            sz = 3
        elseif (problem.dim == 1 || problem.dim == 2 || problem.dim == 3) && r isa ScalarField
            dim = 3
            pdim = 1
            rowsOfB = 3
            b = 1
            sz = 3
        elseif r isa TensorField && nabla == :div
            dim = 3
            pdim = 9
            rowsOfB = 3
            b = 1
            sz = 3
        elseif problem.dim == 2 && problem.type == :AxiSymmetric
            dim = 2
            rowsOfB = 4
            b = 1
        else
            error("∇: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
            #nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(dim, -1, true, false)
            nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes()
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
                B = zeros(rowsOfB * numNodes, pdim * numNodes) # spzeros????
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
                    if dim == 3 && r isa VectorField
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*9-(9-1), l*3-(3-1)] = B[k*9-(9-2), l*3-(3-2)] = B[k*9-(9-3), l*3-(3-3)] = ∂h[1, (k-1)*numNodes+l]
                            B[k*9-(9-4), l*3-(3-1)] = B[k*9-(9-5), l*3-(3-2)] = B[k*9-(9-6), l*3-(3-3)] = ∂h[2, (k-1)*numNodes+l]
                            B[k*9-(9-7), l*3-(3-1)] = B[k*9-(9-8), l*3-(3-2)] = B[k*9-(9-9), l*3-(3-3)] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif (dim == 1 || dim == 2 || dim == 3) && rowsOfB == 3 && r isa ScalarField
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*3-(3-1), l] = ∂h[1, (k-1)*numNodes+l]
                            B[k*3-(3-2), l] = ∂h[2, (k-1)*numNodes+l]
                            B[k*3-(3-3), l] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif r isa TensorField && nabla == :div
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*3-(3-1), l*9-(9-1)] = B[k*3-(3-2), l*9-(9-2)] = B[k*3-(3-3), l*9-(9-3)] = ∂h[1, (k-1)*numNodes+l]
                            B[k*3-(3-1), l*9-(9-4)] = B[k*3-(3-2), l*9-(9-5)] = B[k*3-(3-3), l*9-(9-6)] = ∂h[2, (k-1)*numNodes+l]
                            B[k*3-(3-1), l*9-(9-7)] = B[k*3-(3-2), l*9-(9-8)] = B[k*3-(3-3), l*9-(9-9)] = ∂h[3, (k-1)*numNodes+l]
                        end
                    elseif dim == 2 && rowsOfB == 4
                        for k in 1:numNodes, l in 1:numNodes
                            B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes+l]
                            B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes+l]
                            B[k*4-2, l*2-1] = rr[k] < 1e-10 ? 0 : h[l, k] / rr[k]
                        end
                    else
                        error("∇: rows of B is $rowsOfB, dimension of the problem is $dim.")
                    end
                    push!(numElem, elem)
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnet[j, 1:numNodes] .- (pdim - k)
                    end
                    e = zeros(sz * numNodes, nsteps) # tensors have nine elements
                    for k in 1:numNodes
                        if rowsOfB == 9 && dim == 3 && r isa VectorField
                            B1 = B[k*9-8:k*9, 1:3*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    if nabla == :grad
                                        e[(k-1)*9+1:k*9, kk] = [e0[1], e0[2], e0[3],
                                            e0[4], e0[5], e0[6],
                                            e0[7], e0[8], e0[9]]
                                    elseif nabla == :div
                                        e[k, kk] = e0[1] + e0[5] + e0[9]
                                    elseif nabla == :curl
                                        e[(k-1)*3+1:k*3, kk] = [e0[6] - e0[8],
                                            e0[7] - e0[3],
                                            e0[2] - e0[4]]
                                    end
                                end
                                if DoFResults == true
                                    if nabla == :grad
                                        E1[9*nnet[j, k]-8:9*nnet[j,k], kk] .+= [e0[1], e0[2], e0[3], e0[4], e0[5], e0[6], e0[7], e0[8], e0[9]]
                                    elseif nabla == :div
                                        E1[nnet[j, k], kk] .+= e0[1] + e0[5] + e0[9]
                                    elseif nabla == :curl
                                        E1[3*nnet[j, k]-2:3*nnet[j,k], kk] .+= [e0[6] - e0[8], e0[7] - e0[3], e0[2] - e0[4]]
                                    end
                                end
                            end
                        elseif rowsOfB == 3 && (dim == 1 || dim == 2 || dim == 3) && r isa ScalarField && nabla == :grad
                            B1 = B[k*3-2:k*3, 1:numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*3+1:k*3, kk] = [e0[1], e0[2], e0[3]]
                                end
                                if DoFResults == true
                                    E1[3*nnet[j, k]-2:3*nnet[j,k], kk] .+= [e0[1], e0[2], e0[3]]
                                end
                            end
                        elseif r isa TensorField && nabla == :div
                            B1 = B[k*3-2:k*3, 1:9*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * r.a[nn2, kk]
                                if DoFResults == false
                                    e[(k-1)*3+1:k*3, kk] = [e0[1], e0[2], e0[3]]
                                end
                                if DoFResults == true
                                    E1[3*nnet[j, k]-2:3*nnet[j,k], kk] .+= [e0[1], e0[2], e0[3]]
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
                            error("∇: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
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
        for k in 1:rowsOfB
            for l in 1:non
                E1[k + rowsOfB * l - rowsOfB, :] ./= pcs[l]
            end
        end
    end
    if DoFResults == true
        if r isa VectorField && nabla == :grad
            return TensorField([], E1, r.t, [], nsteps, :e, problem)
        elseif r isa VectorField && nabla == :div
            return ScalarField([], E1, r.t, [], nsteps, :scalar, problem)
        elseif r isa VectorField && nabla == :curl
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        elseif r isa ScalarField
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        elseif r isa TensorField && nabla == :div
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        end
    else
        if r isa VectorField && nabla == :grad
            type = :tensor
            return TensorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        elseif r isa VectorField && nabla == :div
            type = :scalar
            return ScalarField(ε, [;;], r.t, numElem, nsteps, type, problem)
        elseif r isa VectorField && nabla == :curl
            type = :v3D
            return VectorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        elseif r isa ScalarField
            type = :v3D
            return VectorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        elseif r isa TensorField && nabla == :div
            type = :v3D
            return VectorField(ε, [;;], r.t, numElem, nsteps, type, problem)
        end
    end
end
=#

# --- gyors cache-ek elemtípusra (et) ---
const _prop_cache = Dict{Int, Tuple{Int,Int,Vector{Float64}}}()   # et => (dim, numNodes, nodeCoord)
const _basis_cache = Dict{Int, Tuple{Matrix{Float64}, Matrix{Float64}}}()  # et => (∇h, h)

# et tulajdonságok + nodeCoord előállítása és cache-elése
function _get_props_cached(et::Int)
    if haskey(_prop_cache, et)
        return _prop_cache[et]
    end
    elementName, dim, order, numNodes::Int, localNodeCoord, numPrimaryNodes =
        gmsh.model.mesh.getElementProperties(et)
    # nodeCoord: a Gmsh a lokális csomópont-koordinátákat dim-enként adja,
    # mi 3-as blokkokba "csomagoljuk", ahogy az eredeti kód is teszi
    nodeCoord = Vector{Float64}(undef, 3*numNodes)
    @inbounds for j in 1:numNodes, k in 1:dim
        nodeCoord[k + (j-1)*3] = localNodeCoord[k + (j-1)*dim]
    end
    # a 3-dim közötti maradék komponenseket (ha dim<3) automatikusan 0.0-ként hagyjuk (undef -> később nem olvassuk ki)
    # de biztonság kedvéért nullázzuk:
    if dim < 3
        @inbounds for j in 1:numNodes, k in (dim+1):3
            nodeCoord[k + (j-1)*3] = 0.0
        end
    end
    _prop_cache[et] = (dim, numNodes, nodeCoord)
    return _prop_cache[et]
end

# bázisfüggvények cache-elése (∇h és h már a helyes alakra "reshape-elve")
function _get_basis_cached(et::Int, nodeCoord::AbstractVector{<:Real}, numNodes::Int)
    if haskey(_basis_cache, et)
        return _basis_cache[et]
    end
    comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
    ∇h = reshape(dfun, :, numNodes)
    comp, fun,  ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "Lagrange")
    h  = reshape(fun,  :, numNodes)
    _basis_cache[et] = (Array(∇h), Array(h))  # saját példány, hogy véletlen se alias-oljuk a Gmsh bufferét
    return _basis_cache[et]
end

"""
    ∇(r::Union{VectorField, ScalarField, TensorField}; nabla=:grad)

Computes derivatives of `r`.
- If `r` is a `ScalarField` and `nabla == :grad`, returns the gradient (a `VectorField`).
- If `r` is a `VectorField` and `nabla == :grad`, returns the gradient (a `TensorField`).
- If `r` is a `VectorField` and `nabla == :curl`, returns the curl (a `VectorField`).
- If `r` is a `VectorField` and `nabla == :div`, returns the divergence (a `ScalarField`).
- If `r` is a `TensorField` and `nabla == :div`, returns the divergence (a `VectorField`).

Returns: `ScalarField`, `VectorField`, or `TensorField`

Types:
- `r`: `ScalarField`, `VectorField`, or `TensorField`
- `nabla`: Symbol

# 3D Examples (assumes `problem` is set as in the ∇ doc setup)

```julia
# One-time 3D setup (assumes examples/Fields/cube.geo exists with physical group "body")
using LowLevelFEM
gmsh.initialize()
gmsh.open("examples/Fields/cube.geo")
mat = material("body", E=210e3, ν=0.3, ρ=7.85e-9)
problem = Problem([mat], type=:Solid)

# 1) Gradient of a 3D scalar field: ∇f → VectorField
f(X,Y,Z) = X^2 + Y*Z
S = scalarField(problem, [field("body", f=f)])
G = ∇(S)  # VectorField with 3 components

# 2) Curl of a 3D vector field: ∇ × v → VectorField
vx(X,Y,Z) = 0
vy(X,Y,Z) = X
vz(X,Y,Z) = 0
V = vectorField(problem, [field("body", fx=vx, fy=vy, fz=vz)])
C = ∇(V, nabla=:curl)  # approx (0, 0, 1) everywhere

# 3) Divergence of a 3D vector field: ∇ ⋅ v → ScalarField
v1(X,Y,Z) = X
v2(X,Y,Z) = Y
v3(X,Y,Z) = Z
V2 = vectorField(problem, [field("body", fx=v1, fy=v2, fz=v3)])
D = ∇(V2, nabla=:div)  # ≈ 3

# 4) Divergence of a 3D tensor field: ∇ · T → VectorField (if T is TensorField)
# For example, a diagonal tensor T with only Tzz = g(Z): div(T) = (0, 0, ∂g/∂Z)
g(Z) = 10 - Z
T = tensorField(problem, [field("body", fz=g)])
DV = ∇(T, nabla=:div)  # VectorField

# Symmetric displacement gradient via operators
# A = (u ∘ ∇ + ∇ ∘ u) / 2
gmsh.finalize()
```
"""
function ∇(rr::Union{VectorField, ScalarField, TensorField}; nabla=:grad)
    problem = rr.model
    gmsh.model.setCurrent(problem.name)

    # nodális -> elemi térre, ha kell
    r = rr.a == [;;] ? elementsToNodes(rr) : rr

    DoFResults = false
    nsteps = r.nsteps
    ε = Vector{Matrix{Float64}}()   # elemenkénti eredmények gyűjtője
    numElem = Int[]

    # --- getNodes egyszer ---
    nodeTags, ncoord, _ = gmsh.model.mesh.getNodes()
    ncoord2 = zeros(3 * problem.non)
    @inbounds begin
        ncoord2[nodeTags .* 3 .- 2] = ncoord[1:3:length(ncoord)]
        ncoord2[nodeTags .* 3 .- 1] = ncoord[2:3:length(ncoord)]
        ncoord2[nodeTags .* 3 .- 0] = ncoord[3:3:length(ncoord)]
    end

    dim_problem = problem.dim
    pdim = problem.pdim
    non  = problem.non
    type = :e

    # DoFResults ágon: előkészítés (megtartjuk a jelen viselkedést)
    if DoFResults
        E1  = zeros(non * 9, nsteps)
        pcs = zeros(Int, non * dim_problem)
    end

    # anyagcsoportok (0: teljes térfogat, különben material.phName)
    for ipg in 0:length(problem.material)
        phName = ""
        if ipg == 0
            if rr.model.geometry.nameVolume ≠ ""
                phName = rr.model.geometry.nameVolume
            else
                continue
            end
        else
            phName = problem.material[ipg].phName
        end

        # kimenet dimenzió beállítás a bemenet típusa + operátor alapján
        local_dim = 0
        local_pdim = pdim
        rowsOfB = 0
        sz = 0

        #if dim_problem == 3 && r isa VectorField && nabla == :grad
        if r.type == :v3D && r isa VectorField && nabla == :grad
            local_dim = 3; local_pdim = 3; rowsOfB = 9; sz = 9
        #elseif dim_problem == 3 && r isa VectorField && nabla == :div
        elseif r.type == :v3D && r isa VectorField && nabla == :div
            local_dim = 3; local_pdim = 3; rowsOfB = 9; sz = 1
        #elseif dim_problem == 3 && r isa VectorField && nabla == :curl
        elseif r.type == :v3D && r isa VectorField && nabla == :curl
            local_dim = 3; local_pdim = 3; rowsOfB = 9; sz = 3
        elseif (dim_problem == 1 || dim_problem == 2 || dim_problem == 3) && r isa ScalarField
            local_dim = 3; local_pdim = 1; rowsOfB = 3; sz = 3
        elseif r isa TensorField && nabla == :div
            local_dim = 3; local_pdim = 9; rowsOfB = 3; sz = 3
        elseif dim_problem == 2 && problem.type == :AxiSymmetric
            local_dim = 2; rowsOfB = 4
        else
            error("∇: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        @inbounds for idm in 1:length(dimTags)
            edim, etag = dimTags[idm]
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

            @inbounds for i in 1:length(elemTypes)
                et = elemTypes[i]
                et = Int64(et)

                # et-hez cache-elt props + bázis
                dim_et, numNodes, nodeCoord = _get_props_cached(et)
                #@disp dim_et
                ∇h_all, h_all = _get_basis_cached(et, nodeCoord, numNodes)

                # et-hez igazított munkatömbök
                invJac = Matrix{Float64}(undef, 3, 3*numNodes)
                ∂h    = Matrix{Float64}(undef, 3, numNodes*numNodes)
                B     = Matrix{Float64}(undef, rowsOfB * numNodes, local_pdim * numNodes)
                nn2   = Vector{Int}(undef, local_pdim * numNodes)
                rrvec = Vector{Float64}(undef, numNodes)
                nnet  = Matrix{Int}(undef, length(elemTags[i]), numNodes)

                # elem-ciklus
                @inbounds for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]

                    # csomóponttagok betöltése gyors hozzáféréshez
                    @inbounds for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes + k]
                    end

                    # Jacobian (elemszint)
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)

                    # inv(Jac_k)' minden lokális csomópontra
                    @inbounds for k in 1:numNodes
                        # 3x3 blokk kivágása és invertálása
                        #Jk = @view Jac[1:3, 3*k-2:3*k]
                        #invJac[1:3, 3*k-2:3*k] = inv(Matrix(Jk))'
                        ###################################################################
                        # lokális dimenzió (1 a vonal, 2 a felület, 3 a térfogat)
                        dim_et = 3
                        dimJ = dim_et

                        # Jk = 3 × dimJ mátrix
                        #Jk = @view Jac[1:3, (k-1)*dimJ+1 : k*dimJ]
                        Jk = @view Jac[1:dimJ, (k-1)*3+1 : k*3]

                        # inv(Jk) = dimJ × 3 (transzponált mapping)
                        invJk = pinv(Matrix(Jk'))        # stabil Moore–Penrose inverz

                        # 3×3 blokkba másolás (extra oszlopok nullák)
                        fill!(@view(invJac[1:3, 3*k-2:3*k]), 0.0)
                        #@views invJac[1:3, 3*k-2 : 3*k-3+dimJ] .= invJk'
                        #@disp (size(invJac[1:3, 3*k-2 : 3*k-3+dimJ]))
                        #@disp (size(invJk))
                        @views invJac[1:3, 3*k-2 : 3*k-3+dimJ] .= invJk
                        ####################################################################

                        # sugár/rr (axi résznél kellhet), az eredeti kód mintájára az x-koordináták lineáris kombinációja
                        rrvec[k] = (@view h_all[:, k])' * (@view ncoord2[nnet[j, :] .* 3 .- 2])
                    end

                    # ∂h = invJac * ∇h (csomópontpárokra)
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numNodes, l in 1:numNodes
                        # dim_et (1..3) komponenseket vesszük
                        #∂h[1:local_dim, (k-1)*numNodes + l] =
                        #    @view(invJac[1:local_dim, 3*k-2:3*k - (3-local_dim)]) *
                        #    @view(∇h_all[l*3-2:l*3 - (3-local_dim), k])
                        ########################################################################
                            # dim_et = element topological dimension (1,2,3)

                        # local derivatives: (dim_et × 1)
                        g = @view ∇h_all[(l-1)*3+1 : (l-1)*3+dim_et, k]

                        # spatial derivatives: (3 × 1)
                        Jslice = @view invJac[:, 3*k-2 : 3*k-2+dim_et-1]
                        dφdx = Jslice * g

                        # store in ∂h
                        @views ∂h[:, (k-1)*numNodes + l] .= dφdx
                        ########################################################################
                        
                    end
                    
                    # B feltöltése
                    fill!(B, 0.0)
                    if local_dim == 3 && r isa VectorField
                        @inbounds for k in 1:numNodes, l in 1:numNodes
                            B[k*9-8, l*3-2] = B[k*9-7, l*3-1] = B[k*9-6, l*3-0] = ∂h[1, (k-1)*numNodes + l]
                            B[k*9-5, l*3-2] = B[k*9-4, l*3-1] = B[k*9-3, l*3-0] = ∂h[2, (k-1)*numNodes + l]
                            B[k*9-2, l*3-2] = B[k*9-1, l*3-1] = B[k*9-0, l*3-0] = ∂h[3, (k-1)*numNodes + l]
                        end
                    elseif (dim_problem == 1 || dim_problem == 2 || dim_problem == 3) && rowsOfB == 3 && r isa ScalarField
                        @inbounds for k in 1:numNodes, l in 1:numNodes
                            B[k*3-2, l] = ∂h[1, (k-1)*numNodes + l]
                            B[k*3-1, l] = ∂h[2, (k-1)*numNodes + l]
                            B[k*3-0, l] = ∂h[3, (k-1)*numNodes + l]
                        end
                    elseif r isa TensorField && nabla == :div
                        @inbounds for k in 1:numNodes, l in 1:numNodes
                            B[k*3-2, l*9-8] = B[k*3-1, l*9-7] = B[k*3-0, l*9-6] = ∂h[1, (k-1)*numNodes + l]
                            B[k*3-2, l*9-5] = B[k*3-1, l*9-4] = B[k*3-0, l*9-3] = ∂h[2, (k-1)*numNodes + l]
                            B[k*3-2, l*9-2] = B[k*3-1, l*9-1] = B[k*3-0, l*9-0] = ∂h[3, (k-1)*numNodes + l]
                        end
                    elseif local_dim == 2 && rowsOfB == 4
                        @inbounds for k in 1:numNodes, l in 1:numNodes
                            B[k*4-3, l*2-1] = B[k*4-0, l*2-0] = ∂h[1, (k-1)*numNodes + l]
                            B[k*4-1, l*2-0] = B[k*4-0, l*2-1] = ∂h[2, (k-1)*numNodes + l]
                            B[k*4-2, l*2-1] = rrvec[k] < 1e-10 ? 0.0 : h_all[l, k] / rrvec[k]
                        end
                    else
                        error("∇: rowsOfB is $rowsOfB, dimension of the problem is $local_dim.")
                    end

                    push!(numElem, elem)

                    # globális szabadsági vektor-indexek (nn2)
                    @inbounds for k in 1:local_pdim
                        vnn = @view nnet[j, 1:numNodes]
                        nn2[k:local_pdim:local_pdim*numNodes] = local_pdim .* vnn .- (local_pdim - k)
                    end

                    # elemi eredmény mátrix (nem DoF-be gyűjtünk)
                    e = Matrix{Float64}(undef, sz * numNodes, nsteps)

                    @inbounds for k in 1:numNodes
                        if rowsOfB == 9 && local_dim == 3 && r isa VectorField
                            B1 = @view B[k*9-8:k*9, 1:3*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * @view r.a[nn2, kk]
                                if nabla == :grad
                                    @views e[(k-1)*9+1:k*9, kk] = e0
                                elseif nabla == :div
                                    e[k, kk] = e0[1] + e0[5] + e0[9]
                                elseif nabla == :curl
                                    @views e[(k-1)*3+1:k*3, kk] = [e0[6] - e0[8], e0[7] - e0[3], e0[2] - e0[4]]
                                end
                            end
                        elseif rowsOfB == 3 && (dim_problem == 1 || dim_problem == 2 || dim_problem == 3) && r isa ScalarField && nabla == :grad
                            B1 = @view B[k*3-2:k*3, 1:numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * @view r.a[nn2, kk]
                                @views e[(k-1)*3+1:k*3, kk] = e0
                            end
                        elseif r isa TensorField && nabla == :div
                            B1 = @view B[k*3-2:k*3, 1:9*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * @view r.a[nn2, kk]
                                @views e[(k-1)*3+1:k*3, kk] = e0
                            end
                        elseif rowsOfB == 4 && local_dim == 2 && problem.type == :AxiSymmetric
                            B1 = @view B[k*4-3:k*4, 1:2*numNodes]
                            for kk in 1:nsteps
                                e0 = B1 * @view r.a[nn2, kk]
                                # 3D tensor „beágyazás” a 2D axihoz – változatlanul hagyva a mintát
                                @views e[(k-1)*9+1:k*9, kk] = [e0[1], e0[4]/2, 0.0,
                                                                e0[4]/2, e0[3],   0.0,
                                                                0.0,    0.0,     e0[2]]
                            end
                        else
                            error("∇: rowsOfB is $rowsOfB, dimension of the problem is $local_dim, problem type is $(problem.type).")
                        end
                    end

                    push!(ε, e)
                end # elem loop
            end # et loop
        end # idm loop
    end # ipg loop

    # DoFResults ág változatlanul (most is false)
    if DoFResults
        for k in 1:rowsOfB, l in 1:non
            E1[k + rowsOfB*l - rowsOfB, :] ./= pcs[l]
        end
        if r isa VectorField && nabla == :grad
            return TensorField([], E1, r.t, [], nsteps, :e, problem)
        elseif r isa VectorField && nabla == :div
            return ScalarField([], E1, r.t, [], nsteps, :scalar, problem)
        elseif r isa VectorField && nabla == :curl
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        elseif r isa ScalarField
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        elseif r isa TensorField && nabla == :div
            return VectorField([], E1, r.t, [], nsteps, :v3D, problem)
        end
    else
        if r isa VectorField && nabla == :grad
            return TensorField(ε, [;;], r.t, numElem, nsteps, :tensor, problem)
        elseif r isa VectorField && nabla == :div
            return ScalarField(ε, [;;], r.t, numElem, nsteps, :scalar, problem)
        elseif r isa VectorField && nabla == :curl
            return VectorField(ε, [;;], r.t, numElem, nsteps, :v3D, problem)
        elseif r isa ScalarField
            return VectorField(ε, [;;], r.t, numElem, nsteps, :v3D, problem)
        elseif r isa TensorField && nabla == :div
            return VectorField(ε, [;;], r.t, numElem, nsteps, :v3D, problem)
        end
    end
end

"""
    curl(r::VectorField)

Solves the rotation of the vector field `r`.
An alternative way to solve `curl` is to use `∇` as a differencial operator.

Return: VectorField

Types:
- `r`: VectorField

# 3D Example (assumes `problem` is set as in the ∇ doc setup)

```julia
# Assumes a 3D mesh with physical group "body".
vx(X, Y, Z) = 0
vy(X, Y, Z) = X
vz(X, Y, Z) = 0
v = vectorField(problem, [field("body", fx=vx, fy=vy, fz=vz)])
D1 = curl(v)
D2 = ∇ × v
```
"""
function curl(r::VectorField)
    return ∇(r, nabla=:curl)
end

"""
    rot(r::VectorField)

Solves the rotation of the vector field `r`. In some countries "rot" denotes the English "curl".
(See the `curl` function.)

Return: VectorField

Types:
- `r`: VectorField
"""
function rot(r::VectorField)
    return ∇(r, nabla=:curl)
end

import Base.div

"""
    div(r::Union{VectorField,TensorField})

Solves the divergence of the vector field or tensor field `r`.
An alternative way to solve `div` is to use `∇` as a differencial operator.

Return: ScalarField or VectorField

Types:

- `r`: VectorField or TensorField

# 3D Examples (assumes `problem` is set as in the ∇ doc setup)

```julia
# Assumes a 3D mesh with physical group "body".

# 1) Divergence of a 3D vector field → ScalarField
v1(X,Y,Z) = X
v2(X,Y,Z) = Y
v3(X,Y,Z) = Z
v = vectorField(problem, [field("body", fx=v1, fy=v2, fz=v3)])
D1 = div(v)
D2 = ∇ ⋅ v

# 2) Divergence of a 3D tensor field → VectorField
fsz(X, Y, Z) = 10 - Z
S = tensorField(problem, [field("body", fz=fsz)])
b1 = -div(S)
b2 = -S ⋅ ∇
```
"""
function div(r::Union{VectorField,TensorField})
    return ∇(r, nabla=:div)
end

"""
    grad(r::Union{ScalarField,VectorField})

Solves the gradient of the scalar field or vector field `r`.
An alternative way to solve `grad` is to use `∇` as a differencial operator.

Return: VectorField or TensorField

Types:

- `r`: ScalarField or VectorField

# 3D Examples

```julia
# Assumes a 3D mesh with physical group "body".

# 1) Gradient of a 3D scalar field → VectorField
f(X,Y,Z) = X^2 + Y*Z
S = scalarField(problem, [field("body", f=f)])
G1 = grad(S)
G2 = ∇(S)

# 2) Gradient of a 3D vector field → TensorField
vx(X,Y,Z) = X
vy(X,Y,Z) = Y
vz(X,Y,Z) = Z
V = vectorField(problem, [field("body", fx=vx, fy=vy, fz=vz)])
T1 = grad(V)
T2 = V ∘ ∇
```
"""
function grad(r::Union{VectorField,ScalarField})
    return ∇(r)
end

"""
    grad_xy(rr::ScalarField)

Számolja a 2D (xy-síkban fekvő) hálón adott skalármező x és y szerinti deriváltját
az elemi csomópontokban. Eredmény: elemenként eltárolt, 2 komponensű vektormeő
(∂r/∂x, ∂r/∂y).

Megjegyzés: ha a keretrendszered csak :v3D-t kezel, a kimenet típusát cseréld
: v3D-re, és egészítsd ki z=0 komponenssel.
"""

function grad_xy(rr::ScalarField)
    problem = rr.model
    gmsh.model.setCurrent(problem.name)

    # nodális -> elemi térre, ha kell (az eredeti mintát követve)
    r = isElementwise(rr) ? elementsToNodes(rr) : rr

    ε = Vector{Matrix{Float64}}()   # elemenkénti eredmények
    numElem = Int[]
    nsteps = r.nsteps

    # --- getNodes egyszer ---
    nodeTags, ncoord, _ = gmsh.model.mesh.getNodes()
    ncoord2 = zeros(3 * problem.non)
    @inbounds begin
        ncoord2[nodeTags .* 3 .- 2] = ncoord[1:3:length(ncoord)]
        ncoord2[nodeTags .* 3 .- 1] = ncoord[2:3:length(ncoord)]
        ncoord2[nodeTags .* 3 .- 0] = ncoord[3:3:length(ncoord)]
    end

    # anyagcsoportok / fizikai csoportok
    for ipg in 0:length(problem.material)
        phName = ""
        if ipg == 0
            if rr.model.geometry.nameVolume ≠ ""
                phName = rr.model.geometry.nameVolume
            else
                continue
            end
        else
            phName = problem.material[ipg].phName
        end

        # csak 2D entitások (felületek) – az xy síkban fekvő hálóra
        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        @inbounds for idm in 1:length(dimTags)
            edim, etag = dimTags[idm]
            if edim != 2
                continue
            end

            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)

            @inbounds for i in 1:length(elemTypes)
                et = Int(elemTypes[i])

                # cache-elt geometriapropok + bázis
                dim_et, numNodes, nodeCoord = _get_props_cached(et)
                ∇h_all, h_all = _get_basis_cached(et, nodeCoord, numNodes)
                # ∇h_all: (3 × numNodes) × numNodes – az első két sora kell (x,y)

                # munkatömbök
                ∂h    = Matrix{Float64}(undef, 2, numNodes*numNodes)                  # dN/dx, dN/dy páronként
                B     = Matrix{Float64}(undef, 2 * numNodes, 1 * numNodes)            # 2 komponensű grad * 1 skálár
                nn2   = Vector{Int}(undef, 1 * numNodes)                               # DoF-indexek
                nnet  = Matrix{Int}(undef, length(elemTags[i]), numNodes)

                # elem-ciklus
                @inbounds for j in 1:length(elemTags[i])
                    elem = elemTags[i][j]

                    # csomóponttagok betöltése gyors hozzáféréshez
                    @inbounds for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[i][(j-1)*numNodes + k]
                    end

                    # Jacobian a Gmsh-től (lokális pontonként 3×3 blokk), de csak a 2×2 (x,y) részt használjuk
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)

                    # ∂h = inv(J₂×₂) * ∇h_ref(1:2,:) csomópontpárokra
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numNodes
                        Jk = @view Jac[1:2, 3*k-2:3*k-1]   # 2×2 rész (xy)
                        invJ2 = inv(Matrix(Jk))
                        for l in 1:numNodes
                            # ∇h_all szűkítése az első két komponensre (ref. deriváltak)
                            dN_ref = @view ∇h_all[l*3-2:l*3-1, k]   # 2×1
                            @views ∂h[1:2, (k-1)*numNodes + l] = invJ2 * dN_ref
                        end
                    end

                    # B feltöltése (2 komponens: ∂/∂x, ∂/∂y)
                    fill!(B, 0.0)
                    @inbounds for k in 1:numNodes, l in 1:numNodes
                        B[k*2-1, l] = ∂h[1, (k-1)*numNodes + l]   # dN_l/dx a k. lokális ponton
                        B[k*2-0, l] = ∂h[2, (k-1)*numNodes + l]   # dN_l/dy
                    end

                    push!(numElem, elem)

                    # globális szabadsági indexek (skálármező -> 1 DoF/csomópont)
                    @inbounds begin
                        nn2[1:1:1*numNodes] = @view(nnet[j, 1:numNodes])
                    end

                    # elemi eredménymátrix: (2*numNodes) × nsteps
                    e = Matrix{Float64}(undef, 2 * numNodes, nsteps)

                    # node-on-element kiértékelés
                    @inbounds for k in 1:numNodes
                        B1 = @view B[k*2-1:k*2, 1:numNodes]   # 2×numNodes
                        for kk in 1:nsteps
                            e0 = B1 * @view r.a[nn2, kk]      # 2×1
                            @views e[(k-1)*2+1:k*2, kk] = e0
                        end
                    end

                    push!(ε, e)
                end # elem
            end # elemType
        end # idm
    end # ipg

    # 2 komponensű vektormező visszaadása (∂x r, ∂y r)
    return VectorField(ε, [;;], r.t, numElem, nsteps, :v2D, problem)
end

function dx(rr::ScalarField)
    problem = rr.model
    gmsh.model.setCurrent(problem.name)

    # ha nodális tér, konvertáljuk elemisre
    r = rr.a == [;;] ? elementsToNodes(rr) : rr
    nsteps = r.nsteps

    # globális node koordináták gyors elérése
    nodeTags, ncoord, _ = gmsh.model.mesh.getNodes()
    non = problem.non

    ncoord2 = zeros(3 * non)
    @inbounds begin
        ncoord2[nodeTags .* 3 .- 2] = ncoord[1:3:end]
        ncoord2[nodeTags .* 3 .- 1] = ncoord[2:3:end]
        ncoord2[nodeTags .* 3 .- 0] = ncoord[3:3:end]
    end

    ε = Vector{Matrix{Float64}}()    # elemi eredmények
    numElem = Int[]                  # elem azonosítók

    # anyagcsoportok végigjárása
    for ipg in 0:length(problem.material)
        phName = ""
        if ipg == 0
            if rr.model.geometry.nameVolume ≠ ""
                phName = rr.model.geometry.nameVolume
            else
                continue
            end
        else
            phName = problem.material[ipg].phName
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        @inbounds for (edim, etag) in dimTags
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            @inbounds for it in 1:length(elemTypes)
                et = Int(elemTypes[it])

                # elem adatok cache-ből
                dim_et, numNodes, nodeCoord = _get_props_cached(et)
                ∇h_all, h_all = _get_basis_cached(et, nodeCoord, numNodes)

                # munkaváltozók
                invJac = Matrix{Float64}(undef, 3, 3*numNodes)
                ∂h    = Matrix{Float64}(undef, 1, numNodes*numNodes)  # csak ∂/∂x
                nn2   = Vector{Int}(undef, numNodes)                   # scalar → 1 DoF/node
                nnet  = Matrix{Int}(undef, length(elemTags[it]), numNodes)

                # elem ciklus
                for j in 1:length(elemTags[it])
                    elem = elemTags[it][j]

                    # node tags → gyors index
                    @inbounds for k in 1:numNodes
                        nnet[j, k] = elemNodeTags[it][(j-1)*numNodes + k]
                    end

                    # Jacobian
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)

                    # invJac minden lokális csomópontra
                    for k in 1:numNodes
                        Jk = @view Jac[1:dim_et, (k-1)*3+1 : k*3]
                        invJk = pinv(Matrix(Jk'))        # stabil

                        fill!(@view(invJac[1:3, 3*k-2:3*k]), 0.0)
                        @views invJac[1:dim_et, 3*k-2:3*k-3+dim_et] .= invJk
                    end

                    # ∂h = ∂/∂x shape-deriváltak
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numNodes, l in 1:numNodes
                        g = @view ∇h_all[(l-1)*3+1 : (l-1)*3+dim_et, k]
                        Jslice = @view invJac[:, 3*k-2 : 3*k-2+dim_et-1]

                        dphidx = (Jslice * g)[1]   # CSAK ∂/∂x
                        ∂h[1, (k-1)*numNodes + l] = dphidx
                    end

                    # globális scalar DoF indexek
                    @inbounds for k in 1:numNodes
                        nn2[k] = nnet[j, k]        # scalar → csak 1 szabadságfok
                    end

                    # elem eredmény (numNodes × nsteps)
                    e = Matrix{Float64}(undef, numNodes, nsteps)

                    @inbounds for k in 1:numNodes
                        idx = (k-1)*numNodes
                        for kk in 1:nsteps
                            s = 0.0
                            for l in 1:numNodes
                                s += ∂h[1, idx + l] * r.a[nn2[l], kk]
                            end
                            e[k, kk] = s
                        end
                    end

                    push!(ε, e)
                    push!(numElem, elem)
                end
            end
        end
    end

    return ScalarField(ε, [;;], r.t, numElem, nsteps, :scalar, problem)
end

"""
    tangentMatrixConstitutive(r::VectorField)

Solves the constitutive part of the tangent matrix (when solving large deformation problems).
(See [^6]) `r` is the position vector field in the current configuration.

Return: SystemMatrix

Types:
    - `r`: VectorField

[^6]: Javier Bonet, Richard D. Wood: *Nonlinear Continuum Mechanics for Finite Element Analysis*, 
    Cambridge University Press, 2008, <https://doi.org/10.1017/CBO9780511755446>
"""
function tangentMatrixConstitutive(r::VectorField)
    problem = r.model
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
    return SystemMatrix(K, problem)
end

"""
    tangentMatrixInitialStress(r::VectorField)

Solves the initial stress part of the tangent matrix (when solving large deformation problems).
(See [^6]) `r` is the position vector field in the current configuration.

Return: SystemMatrix

Types:

- `r`: VectorField
"""
function tangentMatrixInitialStress(r::VectorField)
    problem = r.model
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
    return SystemMatrix(K, problem)
end

"""
    equivalentNodalForce(r::VectorField)

Solves the equivalent nodal force (when solving large deformation problems).
(See [^6]) `r` is the position vector field in the current configuration.

Return: VectorField

Types:

- `r`: VectorField
"""
function equivalentNodalForce(r::VectorField)
    problem = r.model
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
    return VectorField([], reshape(f, :,1), [0], [], 1, :v3D, problem)
end

"""
    nonFollowerLoadVector(r::VectorField, load)

Solves the non-follower load vector (when solving large deformation problems).
`r` is the position vector field in the current configuration.

Return: VectorField

Types:

- `r`: VectorField
"""
function nonFollowerLoadVector(r::VectorField, loads)
    problem = r.model
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
    F0 = ∇(r)
    F0 = elementsToNodes(F0)
    for n in 1:length(loads)
        #name, fx, fy, fz = loads[n]
        name = loads[n].phName
        fx = loads[n].fx
        fy = loads[n].fy
        fz = loads[n].fz
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
                        y = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 1]
                        z = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 0]
                        if fx !== nothing
                            f[1] = _loadvec_helper(fx, h, x, y, z, nnet, j, l)
                        end
                        if fy !== nothing
                            f[2] = _loadvec_helper(fy, h, x, y, z, nnet, j, l)
                        end
                        if fz !== nothing
                            f[3] = _loadvec_helper(fz, h, x, y, z, nnet, j, l)
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
    return VectorField([], reshape(fp, :,1), [0], [], 1, :v3D, problem)
end

"""
    applyDeformationBoundaryConditions!(deformVec, supports)

Applies displacement boundary conditions `supports` on deformation vector `deformVec`.
Mesh details are in `problem`. `supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Returns: nothing

Types:
- `deformVec`: VectorField
- `supports`: Vector{BoundaryCondition}
"""
function applyDeformationBoundaryConditions!(deformVec, supports; fact=1.0)
    problem = deformVec.model
    if !isa(supports, Vector)
        error("applyDeformationBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    dof, st = size(deformVec.a)
    pdim = problem.pdim

    for i in 1:length(supports)
        #name, ux, uy, uz = supports[i]
        name = supports[i].phName
        ux = supports[i].ux
        uy = supports[i].uy
        uz = supports[i].uz
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux !== nothing
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
            elseif ux isa ScalarField
                uux = ux.a[nodeTags]
            end
            jj = 0
            for j ∈ nodeTagsX
                jj += 1
                if ux isa Function || ux isa ScalarField
                    deformVec.a[j] += uux[jj] * fact
                else
                    deformVec.a[j] += ux * fact
                end
            end
        end
        if uy !== nothing
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
            elseif uy isa ScalarField
                uuy = uy.a[nodeTags]
            end
            jj = 0
            for j ∈ nodeTagsY
                jj += 1
                if uy isa Function || uy isa ScalarField
                    deformVec.a[j] += uuy[jj] * fact
                else
                    deformVec.a[j] += uy * fact
                end
            end
        end
        if pdim == 3 && uz !== nothing
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= pdim
            nodeTagsZ .-= (pdim - 3)
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
            elseif uz isa ScalarField
                uuz = uz.a[nodeTags]
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

"""
    suppressDeformationAtBoundaries!(stiffMat, loadVec, supports)

Suppresses the displacements given in `support` in `stiffMat` and `loadVec` 
so that it is only necessary to consider them once during iteration.
`stiffMat` is the stiffness matrix, `loadVec` is the load vector.
`supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Returns: nothing

Types:
- `stiffMat`: SystemMatrix 
- `loadVec`: VectorField
- `supports`: Vector{BoundaryCondition}
"""
function suppressDeformationAtBoundaries!(stiffMat, loadVec, supports)
#function applyDeformationBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports; fix=1)
    if stiffMat.model != loadVec.model
        error("suppressDeformationAtBoundaries!: system matrix and load vector does not belong to the same problem.")
    end
    problem = stiffMat.model
    if !isa(supports, Vector)
        error("applyBoundaryConditions!: supports are not arranged in a vector. Put them in [...]")
    end
    gmsh.model.setCurrent(problem.name)
    dof, dof = size(stiffMat.A)
    pdim = problem.pdim

    for i in 1:length(supports)
        #name, ux, uy, uz = supports[i]
        name = supports[i].phName
        ux = supports[i].ux
        uy = supports[i].uy
        uz = supports[i].uz
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux !== nothing
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim - 1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
                f0 = stiffMat.A[:, nodeTagsX] * uux * 0
            elseif ux isa ScalarField
                uux = ux.a[nodeTags]
                f0 = stiffMat.A[:, nodeTagsX] * uux * 0
            else
                f0 = stiffMat.A[:, nodeTagsX] * ux * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
        if uy !== nothing
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim - 2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
                f0 = stiffMat.A[:, nodeTagsY] * uuy * 0
            elseif uy isa ScalarField
                uuy = uy.a[nodeTags]
                f0 = stiffMat.A[:, nodeTagsY] * uuy * 0
            else
                f0 = stiffMat.A[:, nodeTagsY] * uy * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
        if pdim == 3 && uz !== nothing
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= 3
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
                f0 = stiffMat.A[:, nodeTagsZ] * uuz * 0
            elseif uz isa ScalarField
                uuz = uz.a[nodeTags]
                f0 = stiffMat.A[:, nodeTagsZ] * uuz * 0
            else
                f0 = stiffMat.A[:, nodeTagsZ] * uz * 0
                f0 = sum(f0, dims=2)
            end
            loadVec.a .-= f0
        end
    end

    for i in 1:length(supports)
        #name, ux, uy, uz = supports[i]
        name = supports[i].phName
        ux = supports[i].ux
        uy = supports[i].uy
        uz = supports[i].uz
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if isa(ux, Function) || isa(uy, Function) || isa(uz, Function)
            xx = coord[1:3:length(coord)]
            yy = coord[2:3:length(coord)]
            zz = coord[3:3:length(coord)]
        end
        if ux !== nothing
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
            if isa(ux, Function)
                uux = ux.(xx, yy, zz)
            elseif ux isa ScalarField
                uux = ux.a[nodeTags]
            end
            jj = 0
            for j ∈ nodeTagsX
                jj += 1
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(ux, Function)
                    loadVec.a[j] = uux[jj] * 0
                elseif ux isa ScalarField
                    loadVec.a[j] = uux[jj] * 0
                else
                    loadVec.a[j] = ux * 0
                end
            end
        end
        if uy !== nothing
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= pdim
            nodeTagsY .-= (pdim-2)
            if isa(uy, Function)
                uuy = uy.(xx, yy, zz)
            elseif uy isa ScalarField
                uuy = uy.a[nodeTags]
            end
            jj = 0
            for j ∈ nodeTagsY
                jj += 1
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(uy, Function)
                    loadVec.a[j] = uuy[jj] * 0
                elseif uy isa ScalarField
                    loadVec.a[j] = uuy[jj] * 0
                else
                    loadVec.a[j] = uy * 0
                end
            end
        end
        if pdim == 3 && uz !== nothing
            nodeTagsZ = copy(nodeTags)
            nodeTagsZ *= pdim
            nodeTagsZ .-= (pdim - 3)
            if isa(uz, Function)
                uuz = uz.(xx, yy, zz)
            elseif uz isa ScalarField
                uuz = uz.a[nodeTags]
            end
            jj = 0
            for j ∈ nodeTagsZ
                jj += 1
                stiffMat.A[j, :] .= 0
                stiffMat.A[:, j] .= 0
                stiffMat.A[j, j] = 1
                #massMat[j, :] .= 0
                #massMat[:, j] .= 0
                #massMat[j, j] = 1
                #dampMat[j, :] .= 0
                #dampMat[:, j] .= 0
                #dampMat[j, j] = 1
                if isa(uz, Function)
                    loadVec.a[j] = uuz[jj] * 0
                elseif uz isa ScalarField
                    loadVec.a[j] = uuz[jj] * 0
                else
                    loadVec.a[j] = uz * 0
                end
            end
        end
    end

    dropzeros!(stiffMat.A)
    #dropzeros!(massMat)
    #dropzeros!(dampMat)
end

"""
    suppressDeformationAtBoundaries(stiffMat, loadVec, supports)

Suppresses the displacements given in `support` in `stiffMat` and `loadVec` 
so that it is only necessary to consider them once during iteration.
`stiffMat` is the stiffness matrix, `loadVec` is the load vector.
`supports` is a tuple of `name` of physical group and
prescribed displacements `ux`, `uy` and `uz`.

Return: stiffMat1, loadVec1

Types:
- `stiffMat`: SystemMatrix 
- `loadVec`: VectorField
- `supports`: Vector{BoundaryCondition}
- `stiffMat1`: SystemMatrix 
- `loadVec1`: VectorField
"""
function suppressDeformationAtBoundaries(stiffMat, loadVec, supports)
    if stiffMat.model != loadVec.model
        error("suppressDeformationAtBoundaries!: system matrix and load vector does not belong to the same problem.")
    end
    problem = stiffMat.model
    K1 = copy(stiffMat)
    f1 = copy(loadVec)
    suppressDeformationAtBoundaries!(K1, f1, supports)
    return K1, f1
end

"""
    solveDeformation(problem::Problem, load, supp;
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

Solves the deformed shape of a non-linearly elastic body...
"""
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
            applyDeformationBoundaryConditions!(r0, supp, fact=1 / loadSteps)
        elseif rampedSupport == false && j == 1
            applyDeformationBoundaryConditions!(r0, supp)
        end
        err = 1
        i = 0
        while err > relativeError && i < maxIteration
            i += 1
       
            Kl = tangentMatrixConstitutive(r0)
            Knl = tangentMatrixInitialStress(r0)
            if followerLoad == false
                f = nonFollowerLoadVector(r0, load)
            end
            fnl = equivalentNodalForce(r0)
            K1, f1 = suppressDeformationAtBoundaries(Kl + Knl, fact * f - fnl, supp)
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
    r1 = VectorField([], r1, 1:size(r1, 2), [], size(r1, 2), :v3D, problem)
    if plotConvergence == true
        return r1, e
    else
        return r1
    end
end

"""
    showDeformationResults(r::VectorField, comp; name=String, visible=Boolean)

Shows deformation result, where `r` contains the position vectors of nodes 
in the *current configuration*.
"""
function showDeformationResults(r::VectorField, comp; name=comp, visible=false)
    problem = r.model
    if r.a == [;;]
        r = elementsToNodes(r)
    end
    r0 = nodePositionVector(problem)
    u = copy(r)
    for i in 1:size(r.a, 2)
        u.a[:, i] = r.a[:, i] - r0.a
    end
    return showDoFResults(u, comp, name=name, visible=visible, factor=1)
end

function showDeformationResults(r::VectorField; name="vector", visible=false)
    return showDeformationResults(r, :vector, name=name, visible=visible)
end

#############################################################
# New approach ##############################################
#############################################################

"""
    materialTangentMatrix(
        problem::Problem;
        F::TensorField,
        C::AbstractMatrix,
        energy::Function,
        params
    )

Assembles the **material (constitutive) tangent stiffness matrix**
for a 3D solid under **finite deformation (Total Lagrange formulation)**,
using a user-provided deformation gradient field `F` and a material
tangent matrix `C` given in **6×6 Mandel/Voigt form**. Alternatively,
it can compute the material tangent from a **free-energy density**
`energy` (hyperelastic formulation) via `Tensors.jl`.

This function builds the matrix
```

K_mat = ∫_Ω B(F)ᵀ · C · B(F) dΩ

```
where `B(F)` is the nonlinear strain–displacement matrix associated with
the Green–Lagrange strain, and `C` is the material tangent
∂S/∂E expressed in Mandel notation.

The routine is **material-agnostic**: it does not assume any specific
constitutive law. The material behavior is entirely defined by the
supplied `C`, which may be spatially varying.

---

## Arguments

- `problem::Problem`  
  The finite element problem definition.  
  Must satisfy:
  - `problem.dim == 3`
  - `problem.pdim == 3`

---

## Keyword arguments

- `F::TensorField`  
  Nodal deformation gradient field.
  Each node stores the full 3×3 deformation gradient `F`,
  which is interpolated to Gauss points during integration.

- `C::AbstractMatrix` (size `6×6`)  
  Material tangent matrix in Mandel/Voigt notation.

  Each entry `C[i,j]` may be either:
  - a `Number` (constant material tangent), or
  - a `ScalarField` (spatially varying material tangent component).

  All entries must satisfy:
```

C[i,j] isa Number || C[i,j] isa ScalarField

```

- `energy::Function`  
  Free-energy density `ψ(C, params)` of a hyperelastic material, evaluated
  at the **right Cauchy–Green tensor** `C = FᵀF`. If provided, the tangent
  is computed at Gauss points via `Tensors.jl` and `C` must be omitted.

- `params`  
  Optional parameter container passed through to `energy`.
---

## Mathematical formulation

- Strain measure: **Green–Lagrange strain**
```

E = 1/2 (FᵀF − I)

```
- Variation:
```

δE = sym(Fᵀ δF)

```
- Tangent contribution:
```

δS = C : δE

```
- Element stiffness contribution:
```

K_e = ∫ B(F)ᵀ · C · B(F) dΩ

```
The strain–displacement matrix `B(F)` is constructed consistently with
Mandel notation (shear components scaled by 2).

---

## Notes

- This function assembles **only the material part** of the tangent
stiffness matrix.  
The **geometric stiffness** (stress-dependent part) must be assembled
separately.

- The formulation is suitable for:
- hyperelastic materials,
- finite rotations,
- finite strains,
provided that `C` is consistent with the stress measure used elsewhere.

- If `energy` is supplied, the routine constructs `C = FᵀF` at Gauss points
  and obtains the constitutive tangent from the free-energy function using
  `Tensors.jl`. This enables arbitrary hyperelastic laws without providing
  a closed-form `C`.

- The material tangent `C` is evaluated at Gauss points by interpolating
nodal `ScalarField` values when necessary.

---

## Returns

- `SystemMatrix`  
Sparse global material tangent stiffness matrix of size `(ndof × ndof)`.

---

## Typical usage

```julia
Kmat = materialTangentMatrix(
  problem;
  F = Ffield,
  C = Cvoigt   # 6×6 matrix of Number / ScalarField
)

K = Kint + Kmat
```

---

## See also

- `initialStressMatrix`

- `internalForceVector`

- `externalTangentFollower`

- `TensorField`

- `ScalarField`  
"""
function materialTangentMatrix(
    problem::Problem;
    F::TensorField,
    C::Union{AbstractMatrix,Nothing}=nothing,
    energy::Union{Nothing,Function}=nothing,
    params=nothing)

    if problem.type == :dummy
        return nothing
    end
    @assert xor(C !== nothing, energy !== nothing)
    @assert problem.dim == 3 && problem.pdim == 3

    gmsh.model.setCurrent(problem.name)

    pdim = 3
    dof = pdim * problem.non

    # --- elementwise nodal F
    Fe = nodesToElements(F)
    Fmap = Dict(zip(Fe.numElem, Fe.A))

    # --- preprocess constant / field C (old path)
    Centry = nothing
    if C !== nothing
        Centry = Matrix{Any}(undef, 6, 6)
        for Icart in CartesianIndices(C)
            cij = C[Icart]
            if cij isa ScalarField
                ce = nodesToElements(cij)
                Centry[Icart] = Dict(zip(ce.numElem, ce.A)) # elem => nodal array
            else
                Centry[Icart] = Float64(cij)
            end
        end
    end

    # --- prealloc sparse triplets (pos-style, no push!)
    len = LowLevelFEM.estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, len)
    J = Vector{Int}(undef, len)
    V = Vector{Float64}(undef, len)
    pos = 1

    for mat in problem.material
        for (_, etag) in gmsh.model.getEntitiesForPhysicalName(mat.phName)

            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(problem.dim, etag)

            for it in eachindex(elemTypes)
                et = elemTypes[it]
                _, _, order, numNodes, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                ip, wip =
                    gmsh.model.mesh.getIntegrationPoints(et, "Gauss$(2order+1)")

                # shape functions
                _, hfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, ip, "Lagrange")
                H = reshape(hfun, numNodes, :)

                # gradients
                _, dfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, ip, "GradLagrange")
                ∇h = reshape(dfun, :, length(wip))

                ndofe = 3 * numNodes

                # --- reusable buffers (element-level)
                B = zeros(6, ndofe)
                Ke = zeros(ndofe, ndofe)
                Cgp = zeros(6, 6)
                Fgp = zeros(3, 3)

                # buffers for Ke update: tmp6 = Cgp*B
                tmp6 = zeros(6, ndofe)

                for (e, elem) in enumerate(elemTags[it])
                    nodeTags = elemNodeTags[it][(e-1)*numNodes+1:e*numNodes]
                    Fnode = Fmap[elem]

                    jac, jacDet, _ =
                        gmsh.model.mesh.getJacobian(elem, ip)
                    Jac = reshape(jac, 3, :)

                    fill!(Ke, 0.0)

                    for k in eachindex(wip)
                        # --- invJ via StaticArrays (fast, no heap)
                        Jm = @SMatrix [
                            Jac[1, 3k-2] Jac[1, 3k-1] Jac[1, 3k];
                            Jac[2, 3k-2] Jac[2, 3k-1] Jac[2, 3k];
                            Jac[3, 3k-2] Jac[3, 3k-1] Jac[3, 3k]
                        ]
                        invJ = transpose(inv(Jm))
                        w = jacDet[k] * wip[k]

                        # --- F at Gauss (avoid reshape alloc by manual fill)
                        fill!(Fgp, 0.0)
                        @inbounds for a in 1:numNodes
                            Na = H[a, k]
                            base = 9a - 8

                            Fgp[1, 1] += Na * Fnode[base+0, 1]
                            Fgp[2, 1] += Na * Fnode[base+1, 1]
                            Fgp[3, 1] += Na * Fnode[base+2, 1]

                            Fgp[1, 2] += Na * Fnode[base+3, 1]
                            Fgp[2, 2] += Na * Fnode[base+4, 1]
                            Fgp[3, 2] += Na * Fnode[base+5, 1]

                            Fgp[1, 3] += Na * Fnode[base+6, 1]
                            Fgp[2, 3] += Na * Fnode[base+7, 1]
                            Fgp[3, 3] += Na * Fnode[base+8, 1]
                        end

                        # --- C at Gauss
                        if energy !== nothing
                            Fm = SMatrix{3,3,Float64}(Fgp)
                            Cmat = transpose(Fm) * Fm
                            Ctan = LowLevelFEM.tangent_from_energy(energy, Cmat, params)  # expected 6×6
                            Cgp .= Ctan
                        else
                            @inbounds for i in 1:6, j in 1:6
                                cij = Centry[i, j]
                                if cij isa Float64
                                    Cgp[i, j] = cij
                                else
                                    nod = cij[elem][:, 1]
                                    Cgp[i, j] = dot(nod, @view H[:, k])
                                end
                            end
                        end

                        # --- build B (no allocations in ∇Na, v)
                        fill!(B, 0.0)
                        @inbounds for a in 1:numNodes
                            gh = @view ∇h[3a-2:3a, k]
                            ∇Na1 = invJ[1, 1] * gh[1] + invJ[1, 2] * gh[2] + invJ[1, 3] * gh[3]
                            ∇Na2 = invJ[2, 1] * gh[1] + invJ[2, 2] * gh[2] + invJ[2, 3] * gh[3]
                            ∇Na3 = invJ[3, 1] * gh[1] + invJ[3, 2] * gh[2] + invJ[3, 3] * gh[3]

                            for j in 1:3
                                v1 = Fgp[j, 1]
                                v2 = Fgp[j, 2]
                                v3 = Fgp[j, 3]
                                col = (a - 1) * 3 + j

                                B[1, col] = v1 * ∇Na1
                                B[2, col] = v2 * ∇Na2
                                B[3, col] = v3 * ∇Na3

                                # engineering shear rows (E23,E13,E12) already with factor 2 in B definition
                                B[4, col] = v2 * ∇Na3 + v3 * ∇Na2
                                B[5, col] = v1 * ∇Na3 + v3 * ∇Na1
                                B[6, col] = v1 * ∇Na2 + v2 * ∇Na1
                            end
                        end

                        # --- Ke += w * (B' * (Cgp*B)) using mul!
                        mul!(tmp6, Cgp, B)                    # tmp6 = Cgp*B
                        mul!(Ke, transpose(B), tmp6, w, 1.0)  # Ke = w*B'*tmp6 + 1*Ke
                    end

                    # --- scatter (pos-style)
                    @inbounds for a in 1:ndofe, b in 1:ndofe
                        I[pos] = (nodeTags[div(a - 1, 3)+1] - 1) * 3 + mod(a - 1, 3) + 1
                        J[pos] = (nodeTags[div(b - 1, 3)+1] - 1) * 3 + mod(b - 1, 3) + 1
                        V[pos] = Ke[a, b]
                        pos += 1
                    end
                end
            end
        end
    end

    resize!(I, pos - 1)
    resize!(J, pos - 1)
    resize!(V, pos - 1)

    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end

"""
    initialStressMatrix(
        problem::Problem;
        stress::TensorField,
        energy::Function,
        C::TensorField,
        params
    )

Assembles the **geometric (initial stress) stiffness matrix**
associated with a given stress field, for finite deformation analysis
in a **Total Lagrange–type formulation**.

The function computes the matrix

```

K_geo = ∫_Ω G(S) dΩ

```

where `G(S)` represents the geometric stiffness contribution induced by
the supplied stress tensor field `stress`. The stress measure itself
(`P`, `S`, or `σ`) is **not interpreted** by the routine; it is treated
purely as a second-order tensor field whose contraction with displacement
gradients yields the geometric stiffness.

---

## Arguments

- `problem::Problem`  
  Finite element problem definition.
  Must satisfy:
  - `problem.dim == problem.pdim`

---

## Keyword arguments

- `stress::TensorField`  
  Nodal stress tensor field.
  
  Each node stores a full `dim×dim` tensor, which is interpolated to Gauss
  points using standard Lagrange shape functions. The physical meaning
  of the tensor is left to the caller; typical choices include:
  - First Piola–Kirchhoff stress `P`,
  - Second Piola–Kirchhoff stress `S`,
  - Cauchy stress `σ`,
  
  provided that the chosen stress measure is **consistent with the
  kinematic description used elsewhere** in the formulation.

- `energy::Function`  
  Free-energy density `ψ(C, params)` of a hyperelastic material. If
  provided, the stress field is computed at Gauss points via `Tensors.jl`
  from the supplied `C` and `stress` must be omitted.

- `C::TensorField`  
  Nodal field of the **right Cauchy–Green tensor** `C = FᵀF`. Required
  when using the energy-based path.

- `params`  
  Optional parameter container passed through to `energy`.

---

## Mathematical formulation

At each Gauss point, the stress tensor `S_gp` is obtained by interpolation
from nodal values. The geometric stiffness contribution is computed as

```

g_ab = (∇N_a)ᵀ · S_gp · (∇N_b)

```

and distributed to the displacement components as

```

(K_geo)_ai,bi = g_ab

```

with no coupling between different displacement directions (`i = i` only).

The element contribution reads

```

K_e = ∫ (∇N)ᵀ · S · (∇N) dΩ

```

and is assembled into the global matrix.

---

## Notes

- This function assembles **only the geometric (stress-dependent) part**
  of the tangent stiffness matrix.
  The **material (constitutive) tangent** must be assembled separately,
  e.g. via `materialTangentMatrix`.

- The routine is **stress-measure agnostic**:
  it does not enforce symmetry or specific push-forward/pull-back
  operations. Ensuring energetic and kinematic consistency between the
  stress field and the rest of the formulation is the responsibility of
  the caller.

- If `energy` is supplied, `C` is interpolated to Gauss points and the
  stress is obtained from the free-energy function using `Tensors.jl`.
  This enables arbitrary hyperelastic laws without explicitly providing
  a stress field.

- The formulation corresponds to the classical **initial stress stiffness**
  encountered in large-deformation and stability (buckling) analyses.

---

## Returns

- `SystemMatrix`  
  Sparse global geometric stiffness matrix of size `(ndof × ndof)`.

---

## Typical usage

```julia
Kgeo = initialStressMatrix(
    problem;
    stress = Sfield   # TensorField (P, S, or σ)
)

K = Kmat + Kgeo
```

---

## See also

* `materialTangentMatrix`
* `internalForceVector`
* `TensorField`
"""
function initialStressMatrix(
    problem::Problem;
    S::Union{TensorField,Nothing}=nothing,
    energy::Union{Nothing,Function}=nothing,
    params=nothing,
    C::Union{Nothing,TensorField}=nothing)

    if problem.type == :dummy
        return nothing
    end
    stress= S
    @assert xor(stress !== nothing, energy !== nothing)
    gmsh.model.setCurrent(problem.name)

    dim = problem.dim
    pdim = dim
    dof = pdim * problem.non

    # --- elementwise maps
    Smap = nothing
    Cmap = nothing
    if stress !== nothing
        Se = nodesToElements(stress)
        Smap = Dict(zip(Se.numElem, Se.A))
    else
        Ce = nodesToElements(C)
        Cmap = Dict(zip(Ce.numElem, Ce.A))
    end

    # --- prealloc sparse triplets (pos-style, no push!)
    len = LowLevelFEM.estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, len)
    J = Vector{Int}(undef, len)
    V = Vector{Float64}(undef, len)
    pos = 1

    for mat in problem.material
        for (_, etag) in gmsh.model.getEntitiesForPhysicalName(mat.phName)

            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(dim, etag)

            for it in eachindex(elemTypes)
                et = elemTypes[it]
                _, _, order, numNodes, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                ip, wip =
                    gmsh.model.mesh.getIntegrationPoints(et, "Gauss$(2order+1)")

                # shape funcs
                _, hfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, ip, "Lagrange")
                H = reshape(hfun, numNodes, :)

                # grad shape funcs
                _, dfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, ip, "GradLagrange")
                ∇h = reshape(dfun, :, length(wip))

                # reusable buffers (element-level)
                Ke = zeros(pdim * numNodes, pdim * numNodes)
                Sgp = zeros(dim, dim)

                # energy path buffers (Gauss-level, reused)
                CgpM = @MMatrix zeros(3, 3)  # mutable accumulator

                for (e, elem) in enumerate(elemTags[it])
                    nodeTags = elemNodeTags[it][(e-1)*numNodes+1:e*numNodes]

                    jac, jacDet, _ =
                        gmsh.model.mesh.getJacobian(elem, ip)
                    Jac = reshape(jac, 3, :)

                    fill!(Ke, 0.0)

                    # fetch element nodal arrays once
                    Selem = stress !== nothing ? Smap[elem] : nothing
                    Celem = energy !== nothing ? Cmap[elem] : nothing

                    for k in eachindex(wip)
                        # --- invJ via StaticArrays (fast, no heap)
                        Jm = @SMatrix [
                            Jac[1, 3k-2] Jac[1, 3k-1] Jac[1, 3k];
                            Jac[2, 3k-2] Jac[2, 3k-1] Jac[2, 3k];
                            Jac[3, 3k-2] Jac[3, 3k-1] Jac[3, 3k]
                        ]
                        invJ = transpose(inv(Jm))
                        w = jacDet[k] * wip[k]

                        if stress !== nothing
                            # --- interpolate stress to Gauss
                            fill!(Sgp, 0.0)
                            @inbounds for a in 1:numNodes
                                Na = H[a, k]
                                base = 9a - 8
                                # manual reshape(Selem[...],3,3) without alloc
                                Sgp[1, 1] += Na * Selem[base+0, 1]
                                Sgp[2, 1] += Na * Selem[base+1, 1]
                                Sgp[3, 1] += Na * Selem[base+2, 1]
                                Sgp[1, 2] += Na * Selem[base+3, 1]
                                Sgp[2, 2] += Na * Selem[base+4, 1]
                                Sgp[3, 2] += Na * Selem[base+5, 1]
                                Sgp[1, 3] += Na * Selem[base+6, 1]
                                Sgp[2, 3] += Na * Selem[base+7, 1]
                                Sgp[3, 3] += Na * Selem[base+8, 1]
                            end
                        else
                            # --- interpolate C to Gauss (MMatrix accumulator)
                            fill!(CgpM, 0.0)
                            @inbounds for a in 1:numNodes
                                Na = H[a, k]
                                base = 9a - 8
                                CgpM[1, 1] += Na * Celem[base+0, 1]
                                CgpM[2, 1] += Na * Celem[base+1, 1]
                                CgpM[3, 1] += Na * Celem[base+2, 1]
                                CgpM[1, 2] += Na * Celem[base+3, 1]
                                CgpM[2, 2] += Na * Celem[base+4, 1]
                                CgpM[3, 2] += Na * Celem[base+5, 1]
                                CgpM[1, 3] += Na * Celem[base+6, 1]
                                CgpM[2, 3] += Na * Celem[base+7, 1]
                                CgpM[3, 3] += Na * Celem[base+8, 1]
                            end
                            # convert once, call stress_from_energy (expects SMatrix{3,3})
                            Cgp = SMatrix{3,3,Float64}(CgpM)
                            S = LowLevelFEM.stress_from_energy(energy, Cgp, params)  # should return 3×3 (SMatrix or Matrix)
                            Sgp .= S
                        end

                        # --- geometric stiffness contribution
                        @inbounds for a in 1:numNodes, b in 1:numNodes
                            # ∇N = invJ * grad_hat (manual for speed)
                            ghA = @view ∇h[3a-2:3a, k]
                            ghB = @view ∇h[3b-2:3b, k]

                            ∇Na1 = invJ[1, 1] * ghA[1] + invJ[1, 2] * ghA[2] + invJ[1, 3] * ghA[3]
                            ∇Na2 = invJ[2, 1] * ghA[1] + invJ[2, 2] * ghA[2] + invJ[2, 3] * ghA[3]
                            ∇Na3 = invJ[3, 1] * ghA[1] + invJ[3, 2] * ghA[2] + invJ[3, 3] * ghA[3]

                            ∇Nb1 = invJ[1, 1] * ghB[1] + invJ[1, 2] * ghB[2] + invJ[1, 3] * ghB[3]
                            ∇Nb2 = invJ[2, 1] * ghB[1] + invJ[2, 2] * ghB[2] + invJ[2, 3] * ghB[3]
                            ∇Nb3 = invJ[3, 1] * ghB[1] + invJ[3, 2] * ghB[2] + invJ[3, 3] * ghB[3]

                            # t = Sgp * ∇Nb
                            t1 = Sgp[1, 1] * ∇Nb1 + Sgp[1, 2] * ∇Nb2 + Sgp[1, 3] * ∇Nb3
                            t2 = Sgp[2, 1] * ∇Nb1 + Sgp[2, 2] * ∇Nb2 + Sgp[2, 3] * ∇Nb3
                            t3 = Sgp[3, 1] * ∇Nb1 + Sgp[3, 2] * ∇Nb2 + Sgp[3, 3] * ∇Nb3

                            g = ∇Na1 * t1 + ∇Na2 * t2 + ∇Na3 * t3

                            # distribute to displacement components (i=i only)
                            for i in 1:dim
                                ia = (a - 1) * pdim + i
                                ib = (b - 1) * pdim + i
                                Ke[ia, ib] += g * w
                            end
                        end
                    end

                    # --- scatter (pos-style)
                    @inbounds for a in 1:pdim*numNodes, b in 1:pdim*numNodes
                        I[pos] = (nodeTags[div(a - 1, pdim)+1] - 1) * pdim + mod(a - 1, pdim) + 1
                        J[pos] = (nodeTags[div(b - 1, pdim)+1] - 1) * pdim + mod(b - 1, pdim) + 1
                        V[pos] = Ke[a, b]
                        pos += 1
                    end
                end
            end
        end
    end

    resize!(I, pos - 1)
    resize!(J, pos - 1)
    resize!(V, pos - 1)

    K = sparse(I, J, V, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end

"""
    internalForceVector(
        problem::Problem;
        P::TensorField,
        F::TensorField,
        energy::Function,
        params
    )

Assembles the **internal force vector** associated with a given
first-order stress tensor field `P`, using a finite deformation
formulation. Alternatively, `P` can be computed from a **free-energy
density** `energy` (hyperelastic formulation) via `Tensors.jl`.

The function computes

```

f_int = ∫_Ω B_Pᵀ · P dΩ

```

where `P` is interpolated to Gauss points and contracted with the
gradients of the shape functions. The routine is **stress-measure
agnostic**: it treats `P` purely as a second-order tensor field and does
not enforce any specific constitutive interpretation.

---

## Arguments

- `problem::Problem`  
  Finite element problem definition.
  The routine currently supports only:
```

problem.dim  == 3
problem.pdim == 3

```

- `P::TensorField`  
Nodal stress tensor field.
Each node stores a full `dim×dim` tensor, which is interpolated to Gauss
points using Lagrange shape functions.

Typical choices for `P` include:
- First Piola–Kirchhoff stress,
- Second Piola–Kirchhoff stress (with compatible kinematics),
- Cauchy stress (if used consistently).

The function does **not** distinguish between these cases; ensuring
consistency with the rest of the formulation is the responsibility of
the caller.

- `F::TensorField`  
  Nodal deformation gradient field. Required when using the energy-based
  path to compute `P` internally.

- `energy::Function`  
  Free-energy density `ψ(C, params)` of a hyperelastic material. If
  provided, the routine constructs `C = FᵀF`, computes `S` via `Tensors.jl`,
  and forms `P = F * S` at Gauss points.

- `params`  
  Optional parameter container passed through to `energy`.

---

## Mathematical formulation

At each Gauss point, the stress tensor is obtained as

```

P_gp = Σ_a N_a P_a

```

The element internal force contribution is computed as

```

(f_e)*{ai} = ∫ (P_gp)*{iJ} (∇N_a)_J dΩ

```

where:
- `a` denotes the node index,
- `i` denotes the displacement component,
- `∇N_a` is the gradient of the shape function in the reference
  configuration.

The resulting element vector is assembled into the global internal force
vector.

---

## Notes

- This function assembles **only the internal force vector**.
  The associated tangent contributions must be assembled separately via:
  - `materialTangentMatrix` (constitutive tangent),
  - `initialStressMatrix` (geometric tangent).

- The stress field `P` is interpolated using standard Lagrange shape
  functions. No symmetry or push-forward/pull-back operations are
  applied internally.

- If `energy` is supplied, `P` is computed at Gauss points from the
  free-energy function using `Tensors.jl`, enabling arbitrary
  hyperelastic laws without explicitly providing a stress field.

- The formulation corresponds to a Total Lagrange–type internal force
  expression when `P` is interpreted as the first Piola–Kirchhoff stress.

---

## Returns

- `VectorField`  
  Global internal force vector of size `(ndof)`, returned as a
  `VectorField` with one time step.

---

## Typical usage

```julia
f_int = internalForceVector(problem, Pfield)

R = f_int - f_ext
```

---

## See also

* `materialTangentMatrix`
* `initialStressMatrix`
* `externalTangentFollower`
* `TensorField`
"""
function internalForceVector(
    problem::Problem;
    P::Union{Nothing,TensorField}=nothing,
    F::Union{Nothing,TensorField}=nothing,
    energy::Union{Nothing,Function}=nothing,
    params=nothing)

    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    @assert xor(P !== nothing, energy !== nothing) """
    Either provide `P` (TensorField),
    or provide (`energy`, `params`, `F`) to compute it at Gauss points.
    """

    dim = problem.dim
    pdim = problem.pdim
    non = problem.non
    dof = pdim * non
    @assert dim == 3 && pdim == 3  # jelen implementáció

    # --- elementwise maps
    Pmap = nothing
    Fmap = nothing
    if P !== nothing
        Pe = nodesToElements(P)
        Pmap = Dict(zip(Pe.numElem, Pe.A))
    else
        @assert F !== nothing
        Fe = nodesToElements(F)
        Fmap = Dict(zip(Fe.numElem, Fe.A))
    end

    f = zeros(dof)

    for mat in problem.material
        for (edim, etag) in gmsh.model.getEntitiesForPhysicalName(mat.phName)

            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for it in eachindex(elemTypes)
                et = elemTypes[it]
                _, _, order, numNodes, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                ip, wip =
                    gmsh.model.mesh.getIntegrationPoints(et, "Gauss$(2order + 1)")

                # --- shape functions
                _, hfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, ip, "Lagrange")
                H = reshape(hfun, numNodes, :)

                # --- grad shape
                _, dfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, ip, "GradLagrange")
                ∇h = reshape(dfun, :, length(wip))

                # --- reusable buffers
                fe = zeros(pdim * numNodes)
                Pgp = zeros(3, 3)  # Float64 buffer
                Fgp = zeros(3, 3)  # only used in energy-path

                for (e, elem) in enumerate(elemTags[it])
                    nodeTags = elemNodeTags[it][(e-1)*numNodes+1:e*numNodes]

                    jac, jacDet, _ =
                        gmsh.model.mesh.getJacobian(elem, ip)
                    Jac = reshape(jac, 3, :)

                    Pelem = (P !== nothing) ? Pmap[elem] : nothing
                    Felem = (energy !== nothing) ? Fmap[elem] : nothing

                    fill!(fe, 0.0)

                    for k in eachindex(wip)
                        # --- invJ via StaticArrays
                        Jm = @SMatrix [
                            Jac[1, 3k-2] Jac[1, 3k-1] Jac[1, 3k];
                            Jac[2, 3k-2] Jac[2, 3k-1] Jac[2, 3k];
                            Jac[3, 3k-2] Jac[3, 3k-1] Jac[3, 3k]
                        ]
                        invJ = transpose(inv(Jm))
                        w = jacDet[k] * wip[k]

                        if P !== nothing
                            # --- interpolate P to Gauss into Pgp buffer (no alloc)
                            fill!(Pgp, 0.0)
                            @inbounds for a in 1:numNodes
                                Na = H[a, k]
                                base = 9a - 8
                                # column-major (ugyanaz, mint reshape(x,3,3))
                                Pgp[1, 1] += Na * Pelem[base+0, 1]
                                Pgp[2, 1] += Na * Pelem[base+1, 1]
                                Pgp[3, 1] += Na * Pelem[base+2, 1]

                                Pgp[1, 2] += Na * Pelem[base+3, 1]
                                Pgp[2, 2] += Na * Pelem[base+4, 1]
                                Pgp[3, 2] += Na * Pelem[base+5, 1]

                                Pgp[1, 3] += Na * Pelem[base+6, 1]
                                Pgp[2, 3] += Na * Pelem[base+7, 1]
                                Pgp[3, 3] += Na * Pelem[base+8, 1]
                            end
                        else
                            # --- energy path: Fgp -> C -> S -> P = F*S
                            fill!(Fgp, 0.0)
                            @inbounds for a in 1:numNodes
                                Na = H[a, k]
                                base = 9a - 8
                                Fgp[1, 1] += Na * Felem[base+0, 1]
                                Fgp[2, 1] += Na * Felem[base+1, 1]
                                Fgp[3, 1] += Na * Felem[base+2, 1]
                                Fgp[1, 2] += Na * Felem[base+3, 1]
                                Fgp[2, 2] += Na * Felem[base+4, 1]
                                Fgp[3, 2] += Na * Felem[base+5, 1]
                                Fgp[1, 3] += Na * Felem[base+6, 1]
                                Fgp[2, 3] += Na * Felem[base+7, 1]
                                Fgp[3, 3] += Na * Felem[base+8, 1]
                            end

                            Fm = SMatrix{3,3,Float64}(Fgp)
                            Cmat = transpose(Fm) * Fm
                            S = LowLevelFEM.stress_from_energy(energy, Cmat, params)  # 3×3 (SMatrix ok)
                            Pm = Fm * S

                            # copy back into Pgp buffer (avoid Matrix(Pm))
                            @inbounds begin
                                Pgp[1, 1] = Pm[1, 1]
                                Pgp[1, 2] = Pm[1, 2]
                                Pgp[1, 3] = Pm[1, 3]
                                Pgp[2, 1] = Pm[2, 1]
                                Pgp[2, 2] = Pm[2, 2]
                                Pgp[2, 3] = Pm[2, 3]
                                Pgp[3, 1] = Pm[3, 1]
                                Pgp[3, 2] = Pm[3, 2]
                                Pgp[3, 3] = Pm[3, 3]
                            end
                        end

                        # --- fe accumulation
                        @inbounds for a in 1:numNodes
                            gh = @view ∇h[3a-2:3a, k]
                            ∇Na1 = invJ[1, 1] * gh[1] + invJ[1, 2] * gh[2] + invJ[1, 3] * gh[3]
                            ∇Na2 = invJ[2, 1] * gh[1] + invJ[2, 2] * gh[2] + invJ[2, 3] * gh[3]
                            ∇Na3 = invJ[3, 1] * gh[1] + invJ[3, 2] * gh[2] + invJ[3, 3] * gh[3]

                            base = (a - 1) * pdim
                            # dot(Pgp[i,:], ∇Na) kézzel, hogy ne allokáljon view/vec
                            fe[base+1] += (Pgp[1, 1] * ∇Na1 + Pgp[1, 2] * ∇Na2 + Pgp[1, 3] * ∇Na3) * w
                            fe[base+2] += (Pgp[2, 1] * ∇Na1 + Pgp[2, 2] * ∇Na2 + Pgp[2, 3] * ∇Na3) * w
                            fe[base+3] += (Pgp[3, 1] * ∇Na1 + Pgp[3, 2] * ∇Na2 + Pgp[3, 3] * ∇Na3) * w
                        end
                    end

                    # scatter fe into global f
                    @inbounds for a in 1:numNodes
                        nt = nodeTags[a]
                        base = (a - 1) * pdim
                        f[(nt-1)*pdim+1] += fe[base+1]
                        f[(nt-1)*pdim+2] += fe[base+2]
                        f[(nt-1)*pdim+3] += fe[base+3]
                    end
                end
            end
        end
    end

    return VectorField([], reshape(f, :, 1), [0.0], [], 1, :v3D, problem)
end

@inline function _tangent_load_helper(f, h, x, y, z, nnet, j, l, nsteps)
    if f isa Number
        return fill(f, nsteps)
    elseif f isa Function
        return fill(f(x, y, z), nsteps)
    elseif f isa ScalarField && isNodal(f)
        v = h[:, j]' * f.a[nnet[l, :], :]
        return vec(v)
    else
        error("externalTangentFollowerTL: invalid load type.")
    end
end

"""
    externalTangentFollower(
        problem::Problem,
        loads::Vector{BoundaryCondition};
        F::TensorField
    )

Assembles the **external tangent stiffness matrix associated with
follower (configuration-dependent) surface loads** for a finite
deformation analysis.

The function computes the linearization of the external force vector
with respect to the displacement field when the applied tractions rotate
and/or stretch together with the deforming body.

---

## Arguments

- `problem::Problem`  
  Finite element problem definition.

- `loads::Vector{BoundaryCondition}`  
  Boundary conditions defining surface tractions or pressures acting on
  the current configuration. Supported load types:
  - surface traction components (`fx`, `fy`, `fz`),
  - surface pressure (`p`, only for 3D solids).

  Each boundary condition is interpreted as a **follower load**, i.e.
  its direction and/or magnitude may depend on the deformation.

---

## Keyword arguments

- `F::TensorField`  
  Nodal deformation gradient field.
  Each node stores the full deformation gradient `F`, which is
  interpolated to Gauss points and used to:
  - rotate the applied tractions,
  - compute the Jacobian and inverse deformation gradient,
  - form the consistent external tangent matrix.

---

## Mathematical formulation

The external force vector for follower tractions may be written as

```

f_ext = ∫_Γ Nᵀ · t(F) dΓ

```

where the traction vector depends on the deformation gradient, typically
through a relation of the form

```

t = J F⁻ᵀ t₀    or    t = F · t_loc

```

The present implementation linearizes this expression consistently,
leading to an external tangent contribution

```

K_ext = ∂f_ext / ∂u

```

which includes:
- the variation of the Jacobian and inverse deformation gradient,
- the variation of the traction due to rotation/stretching with `F`.

Both pressure loads and vector-valued surface tractions are supported.

---

## Notes

- This function assembles **only the tangent contribution** of follower
  loads.  
  The corresponding external force vector must be assembled separately,
  e.g. using `loadVector` with deformation-dependent tractions.

- Dead loads (loads fixed in the reference configuration) **must not**
  be passed to this function. For dead loads, the external tangent is
  identically zero.

- The formulation assumes a **Total Lagrange–type description** and is
  intended for large-rotation and large-deformation problems.

- For pressure loads, the reference normal vector is obtained from the
  undeformed geometry and appropriately transformed.

---

## Returns

- `SystemMatrix`  
  Sparse global external tangent stiffness matrix of size `(ndof × ndof)`.

---

## Typical usage

```julia
Kext = externalTangentFollower(
    problem,
    loads;
    F = Ffield
)

K = Kmat + Kgeo - Kext
```

---

## See also

* `loadVector`
* `internalForceVector`
* `materialTangentMatrix`
* `initialStressVector`
"""
function externalTangentFollower(
    problem::Problem,
    loads::Vector{BoundaryCondition};
    F::TensorField)

    if problem.type == :dummy
        return nothing
    end
    gmsh.model.setCurrent(problem.name)

    dim  = problem.dim
    pdim = problem.pdim
    non  = problem.non
    dof  = pdim * non
    @assert dim == 3 && pdim == 3

    # --- elementwise F
    Fe   = nodesToElements(F)
    Fmap = Dict(zip(Fe.numElem, Fe.A))

    #Id = Matrix{Float64}(I, 3, 3)
    Id = [1.0 0 0; 0 1 0; 0 0 1]

    # --- pos-style sparse
    len = LowLevelFEM.estimateLengthOfIJV(problem)
    I = Vector{Int}(undef, len)
    J = Vector{Int}(undef, len)
    V = Vector{Float64}(undef, len)
    pos = 1

    # --- global node coordinates cache
    ncoord2 = zeros(3 * non)

    for bc in loads
        name = bc.phName

        fx = bc.fx
        fy = bc.fy
        fz = bc.fz
        p  = bc.p

        fx = (fx isa ScalarField && isElementwise(fx)) ? elementsToNodes(fx) : fx
        fy = (fy isa ScalarField && isElementwise(fy)) ? elementsToNodes(fy) : fy
        fz = (fz isa ScalarField && isElementwise(fz)) ? elementsToNodes(fz) : fz
        p  = (p  isa ScalarField && isElementwise(p )) ? elementsToNodes(p ) : p

        if p !== nothing && (fx !== nothing || fy !== nothing || fz !== nothing)
            error("externalTangentFollowerTL: p and fx/fy/fz cannot be defined together.")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(name)

        for (edim, etag) in dimTags
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            nodeTagsAll, ncoord, _ =
                gmsh.model.mesh.getNodes(edim, etag, true, false)

            @inbounds begin
                ncoord2[nodeTagsAll*3 .- 2] = ncoord[1:3:end]
                ncoord2[nodeTagsAll*3 .- 1] = ncoord[2:3:end]
                ncoord2[nodeTagsAll*3 .- 0] = ncoord[3:3:end]
            end

            nv = (p !== nothing) ? -normalVector(problem, name) : nothing

            for it in eachindex(elemTypes)
                et = elemTypes[it]
                _, _, order, numNodes, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                ip, wip =
                    gmsh.model.mesh.getIntegrationPoints(et, "Gauss$(2order + 1)")

                _, hfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, ip, "Lagrange")
                h = reshape(hfun, numNodes, :)

                _, dfun, _ =
                    gmsh.model.mesh.getBasisFunctions(et, ip, "GradLagrange")
                ∇h = reshape(dfun, :, length(wip))

                ndofe = pdim * numNodes

                # --- reusable buffers
                Ke   = zeros(ndofe, ndofe)
                Fgp  = zeros(3, 3)
                tgp  = zeros(3)
                dt0  = zeros(3)
                tloc = zeros(3)

                for (l, elem) in enumerate(elemTags[it])
                    nodeTags =
                        @view elemNodeTags[it][(l-1)*numNodes+1:l*numNodes]

                    jac, jacDet, _ =
                        gmsh.model.mesh.getJacobian(elem, ip)
                    Jac = reshape(jac, 3, :)

                    Fnode = Fmap[elem]
                    fill!(Ke, 0.0)

                    for j in eachindex(wip)
                        Jm = @SMatrix [
                            Jac[1,3j-2] Jac[1,3j-1] Jac[1,3j];
                            Jac[2,3j-2] Jac[2,3j-1] Jac[2,3j];
                            Jac[3,3j-2] Jac[3,3j-1] Jac[3,3j]
                        ]
                        invJ = transpose(inv(Jm))
                        w = jacDet[j] * wip[j]

                        # --- F at Gauss (reshape-safe)
                        fill!(Fgp, 0.0)
                        @inbounds for a in 1:numNodes
                            Na = h[a, j]
                            base = 9a - 8
                            Fgp[1,1] += Na * Fnode[base+0,1]
                            Fgp[2,1] += Na * Fnode[base+1,1]
                            Fgp[3,1] += Na * Fnode[base+2,1]
                            Fgp[1,2] += Na * Fnode[base+3,1]
                            Fgp[2,2] += Na * Fnode[base+4,1]
                            Fgp[3,2] += Na * Fnode[base+5,1]
                            Fgp[1,3] += Na * Fnode[base+6,1]
                            Fgp[2,3] += Na * Fnode[base+7,1]
                            Fgp[3,3] += Na * Fnode[base+8,1]
                        end

                        Fm = SMatrix{3,3,Float64}(Fgp)
                        Finv = inv(Fm)
                        FinvT = transpose(Finv)
                        Jgp = det(Fm)

                        # --- Gauss coordinates
                        x = 0.0; y = 0.0; z = 0.0
                        @inbounds for a in 1:numNodes
                            Na = h[a,j]
                            nt = nodeTags[a]
                            x += Na * ncoord2[3nt-2]
                            y += Na * ncoord2[3nt-1]
                            z += Na * ncoord2[3nt-0]
                        end

                        # --- traction
                        fill!(tgp, 0.0)
                        if p !== nothing
                            pval = _tangent_load_helper(p, h, x, y, z, nodeTags, j, 1, 1)[1]
                            for d in 1:3
                                tgp[d] = -pval * (h[:,j]' * nv[d].a[nodeTags,1])
                            end
                        else
                            if fx !== nothing
                                tgp[1] = _tangent_load_helper(fx, h, x, y, z, nodeTags, j, 1, 1)[1]
                            end
                            if fy !== nothing
                                tgp[2] = _tangent_load_helper(fy, h, x, y, z, nodeTags, j, 1, 1)[1]
                            end
                            if fz !== nothing
                                tgp[3] = _tangent_load_helper(fz, h, x, y, z, nodeTags, j, 1, 1)[1]
                            end
                        end

                        tgp .= Fgp * tgp
                        tloc .= Finv * tgp

                        for a in 1:numNodes, b in 1:numNodes
                            ghb = @view ∇h[3b-2:3b, j]
                            ∇Nb = invJ * ghb

                            for jb in 1:3
                                # --- unchanged physics
                                fcol = Finv[:,jb]
                                α = dot(fcol, ∇Nb)

                                dJFmT = Jgp * FinvT * (α*Id - (∇Nb*fcol'))
                                dt0 .= dJFmT * tgp

                                dFt = dot(∇Nb, tloc)
                                @inbounds for ia in 1:3
                                    dt0[ia] += Jgp * FinvT[ia,jb] * dFt
                                end

                                @inbounds for ia in 1:3
                                    Ke[(a-1)*3+ia,(b-1)*3+jb] +=
                                        h[a,j] * dt0[ia] * w
                                end
                            end
                        end
                    end

                    @inbounds for a in 1:ndofe, b in 1:ndofe
                        I[pos] = (nodeTags[div(a-1,3)+1]-1)*3 + mod(a-1,3) + 1
                        J[pos] = (nodeTags[div(b-1,3)+1]-1)*3 + mod(b-1,3) + 1
                        V[pos] = Ke[a,b]
                        pos += 1
                    end
                end
            end
        end
    end

    resize!(I,pos-1); resize!(J,pos-1); resize!(V,pos-1)
    K = sparse(I,J,V,dof,dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end

"""
    IIPiolaKirchhoff(energy, C::TensorField, params)

Computes the nodal Second Piola–Kirchhoff (II PK) stress tensor field
from a given nodal Right Cauchy–Green deformation tensor field.

For each time step and each node, the II PK stress is evaluated as

    S = ∂ψ(C, params) / ∂C

using the provided free energy function `energy` and material parameters
`params`. The stress evaluation is performed directly at the nodal
values of `C`, without Gauss-point integration or spatial averaging.

This function is intended for post-processing, visualization, and
diagnostic purposes. The resulting stress field is energy-consistent
with the constitutive model but is not suitable for assembling internal
forces or tangents, which require Gauss-point evaluation.

Arguments:
- `energy::Function`:
    Free energy density ψ(C, params).
- `C::TensorField`:
    Nodal Right Cauchy–Green deformation tensor field.
- `params::NamedTuple`:
    Material parameters passed to `energy`.

Requirements:
- The problem must be three-dimensional (`dim = pdim = 3`).
- The input field `C` must be defined nodally.

Returns:
- `TensorField`:
    Nodal Second Piola–Kirchhoff stress field with the same time steps
    and model as `C`.

Notes:
- The stress is computed independently at each node and time step.
- No interpolation, integration, or smoothing is performed.
- For use in weak forms or Newton iterations, Gauss-point evaluation
  should be used instead.
"""
function IIPiolaKirchhoff(
    energy::Function,
    C::TensorField,
    params::NamedTuple)

    problem = C.model
    @assert problem.dim == 3 && problem.pdim == 3

    C0 = elementsToNodes(C)
    non = problem.non
    nsteps = C.nsteps
    S = zeros(9 * non, nsteps)

    @inbounds for j in 1:nsteps
        @inbounds for a in 1:non
            base = 9*(a-1)

            # --- C at node (column-major safe)
            Cn = @SMatrix [
                C0.a[base+1,j] C0.a[base+4,j] C0.a[base+7,j];
                C0.a[base+2,j] C0.a[base+5,j] C0.a[base+8,j];
                C0.a[base+3,j] C0.a[base+6,j] C0.a[base+9,j]
            ]

            # --- II PK stress
            Sn = LowLevelFEM.stress_from_energy(energy, Cn, params)

            # --- store back (same layout as TensorField)
            S[base+1,j] = Sn[1,1]
            S[base+2,j] = Sn[2,1]
            S[base+3,j] = Sn[3,1]
            S[base+4,j] = Sn[1,2]
            S[base+5,j] = Sn[2,2]
            S[base+6,j] = Sn[3,2]
            S[base+7,j] = Sn[1,3]
            S[base+8,j] = Sn[2,3]
            S[base+9,j] = Sn[3,3]
        end
    end

    return TensorField([], S, C.t, [], nsteps, :tensor, C.model)
end

# ========= PUBLIC API =========

function stress_from_energy(ψ, C::SMatrix, p)
    try
        return _stress_from_energy(ψ, C, p)
    catch err
        if err isa MethodError && err.f === _stress_from_energy
            error("""
            Energy-based stress computation requires Tensors.jl.

            Please load it first:
                using Tensors
            """)
        else
            rethrow()
        end
    end
end

function tangent_from_energy(ψ, C::SMatrix, p)
    try
        return _tangent_from_energy(ψ, C, p)
    catch err
        if err isa MethodError && err.f === _tangent_from_energy
            error("""
            Energy-based tangent computation requires Tensors.jl.

            Please load it first:
                using Tensors
            """)
        else
            rethrow()
        end
    end
end


# ========= EXTENSION HOOKS (NO METHODS!) =========

function _stress_from_energy end
function _tangent_from_energy end

# ===========================================================================



#=
@inline function _rotation_from_F(F::AbstractMatrix{<:Real})
    U, _, Vt = svd(F)
    R = U * Vt
    if det(R) < 0
        U[:, end] .*= -1
        R = U * Vt
    end
    return R
end
=#

#=
@inline function _loadvec_helper(f, h, x, y, z, nnet, j, l, nsteps)
    if f isa Number
        return fill(f, nsteps)
    elseif f isa Function
        return fill(f(x, y, z), nsteps)
    elseif f isa ScalarField && isNodal(f)
        v = h[:, j]' * f.a[nnet[l, :], :]
        return vec(v)
    else
        error("loadVector: internal error.")
    end
end

function loadVectorNew(problem, loads; F=nothing)
    gmsh.model.setCurrent(problem.name)
    if !isa(loads, Vector)
        error("loadVector: loads are not arranged in a vector. Put them in [...]")
    end

    Fe = nothing
    Fmap = nothing
    if F !== nothing
        Fe = nodesToElements(F)
        Fmap = Dict(zip(Fe.numElem, Fe.A))
    end

    pdim = problem.pdim
    DIM = problem.dim
    b = problem.thickness
    non = problem.non
    dof = pdim * non
    ncoord2 = zeros(3 * problem.non)
    f = nothing
    fp = zeros(dof,1)
    nsteps = 1
    for n in 1:length(loads)
        name = loads[n].phName
        fx = loads[n].fx
        fy = loads[n].fy
        fz = loads[n].fz
        T = loads[n].T
        p = loads[n].p
        hc = loads[n].h
        T0 = loads[n].T
        qn = loads[n].qn
        hs = loads[n].h

        q = loads[n].q

        fxy = loads[n].fxy
        fyz = loads[n].fyz
        fzx = loads[n].fzx
        
        fyx = loads[n].fyx
        fzy = loads[n].fzy
        fxz = loads[n].fxz

        fx  = (fx  isa ScalarField && isElementwise(fx))  ? elementsToNodes(fx)  : fx
        fy  = (fy  isa ScalarField && isElementwise(fy))  ? elementsToNodes(fy)  : fy
        fz  = (fz  isa ScalarField && isElementwise(fz))  ? elementsToNodes(fz)  : fz
        T   = (T   isa ScalarField && isElementwise(T))   ? elementsToNodes(T)   : T
        p   = (p   isa ScalarField && isElementwise(p))   ? elementsToNodes(p)   : p
        hc  = (hc  isa ScalarField && isElementwise(hc))  ? elementsToNodes(hc)  : hc
        T0  = (T0  isa ScalarField && isElementwise(T0))  ? elementsToNodes(T0)  : T0
        qn  = (qn  isa ScalarField && isElementwise(qn))  ? elementsToNodes(qn)  : qn
        hs  = (hs  isa ScalarField && isElementwise(hs))  ? elementsToNodes(hs)  : hs
        q   = (q   isa ScalarField && isElementwise(q))   ? elementsToNodes(q)   : q
        fxy = (fxy isa ScalarField && isElementwise(fxy)) ? elementsToNodes(fxy) : fxy
        fyz = (fyz isa ScalarField && isElementwise(fyz)) ? elementsToNodes(fyz) : fyz
        fzx = (fzx isa ScalarField && isElementwise(fzx)) ? elementsToNodes(fzx) : fzx
        fyx = (fyx isa ScalarField && isElementwise(fyx)) ? elementsToNodes(fyx) : fyx
        fzy = (fzy isa ScalarField && isElementwise(fzy)) ? elementsToNodes(fzy) : fzy
        fxz = (fxz isa ScalarField && isElementwise(fxz)) ? elementsToNodes(fxz) : fxz
    
        nsteps = 1
        for i in [fx, fy , fz, T,  p, hc, T0, qn, hs, q, fxy, fyz, fzx, fyx, fzy, fxz]
            if i isa ScalarField
                nsteps = max(nsteps, i.nsteps)
            end
        end
        if nsteps > size(fp, 2)
            fp = hcat(fp, zeros(dof, nsteps - size(fp,2)))
        end

        (qn !== nothing || hc !== nothing || hs !== nothing || T !== nothing) && (fx !== nothing || fy !== nothing || fz !== nothing) &&
            error("loadVector: qn/h/T∞ and fx/fy/fz cannot be defined in the same BC.")
        if pdim == 1 || pdim == 2 || pdim == 3 || pdim == 9
            f = zeros(pdim, nsteps)
        else
            error("loadVector: dimension of the problem is $(problem.dim).")
        end
        if p !== nothing && DIM == 3 && pdim == 3
            nv = -normalVector(problem, name)
            ex = VectorField(problem, [field(name, fx=1, fy=0, fz=0)])
            ey = VectorField(problem, [field(name, fx=0, fy=1, fz=0)])
            ez = VectorField(problem, [field(name, fx=0, fy=0, fz=1)])
            if p isa Number || p isa ScalarField
                fy = elementsToNodes((nv ⋅ ey) * p)
                fz = elementsToNodes((nv ⋅ ez) * p)
                fx = elementsToNodes((nv ⋅ ex) * p)
            elseif p isa Function
                pp = scalarField(problem, [field(name, f=p)])
                fy = elementsToNodes((nv ⋅ ey) * pp)
                fz = elementsToNodes((nv ⋅ ez) * pp)
                fx = elementsToNodes((nv ⋅ ex) * pp)
            end
        end
        if p !== nothing && DIM ≠ 3
            error("loadVector: pressure can be given on a surface of a 3D solid.")
        end
        fx = fx !== nothing ? fx : 0.0
        fy = fy !== nothing ? fy : 0.0
        fz = fz !== nothing ? fz : 0.0
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
                f1 = zeros(pdim * numNodes, nsteps)
                nn2 = zeros(Int, pdim * numNodes)
                @inbounds for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    for k in 1:numNodes
                        nnet[l, k] = elemNodeTags[ii][(l-1)*numNodes+k]
                    end
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    fill!(f1, 0.0)
                    @inbounds for j in 1:numIntPoints
                        x = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 2]
                        y = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 1]
                        z = h[:, j]' * ncoord2[nnet[l, :] * 3 .- 0]
                        if hc !== nothing && T0 !== nothing
                            if hc isa Function
                                error("heatConvectionVector: h cannot be a function.")
                            end
                            f[1,:] .= _loadvec_helper(T0, h, x, y, z, nnet, j, l, nsteps) * hc
                        elseif qn !== nothing
                            f[1,:] .= _loadvec_helper(qn, h, x, y, z, nnet, j, l, nsteps)
                        elseif hs !== nothing
                            f[1,:] .= _loadvec_helper(hs, h, x, y, z, nnet, j, l, nsteps)
                        elseif q !== nothing
                            f[1,:] .= _loadvec_helper(q, h, x, y, z, nnet, j, l, nsteps)
                        elseif fx !== nothing && pdim <= 3
                            f[1,:] .= _loadvec_helper(fx, h, x, y, z, nnet, j, l, nsteps)
                        end
                        if pdim > 1 && pdim <= 3
                            f[2,:] .= _loadvec_helper(fy, h, x, y, z, nnet, j, l, nsteps)
                        end
                        if pdim == 3
                            f[3,:] .= _loadvec_helper(fz, h, x, y, z, nnet, j, l, nsteps)
                        end
                        if pdim == 9
                            fill!(f, 0.0)
                            if fx !== nothing
                                f[1,:] .= _loadvec_helper(fx, h, x, y, z, nnet, j, l, nsteps)
                            end
                            if fy !== nothing
                                f[5,:] .= _loadvec_helper(fy, h, x, y, z, nnet, j, l, nsteps)
                            end
                            if fz !== nothing
                                f[9,:] .= _loadvec_helper(fz, h, x, y, z, nnet, j, l, nsteps)
                            end
                            if fxy !== nothing
                                f[4,:] .= _loadvec_helper(fxy, h, x, y, z, nnet, j, l, nsteps)
                            end
                            if fyz !== nothing
                                f[8,:] .= _loadvec_helper(fyz, h, x, y, z, nnet, j, l, nsteps)
                            end
                            if fzx !== nothing
                                f[3,:] .= _loadvec_helper(fzx, h, x, y, z, nnet, j, l, nsteps)
                            end
                            if fyx !== nothing
                                f[2,:] .= _loadvec_helper(fyx, h, x, y, z, nnet, j, l, nsteps)
                            end
                            if fzy !== nothing
                                f[6,:] .= _loadvec_helper(fzy, h, x, y, z, nnet, j, l, nsteps)
                            end
                            if fxz !== nothing
                                f[7,:] .= _loadvec_helper(fxz, h, x, y, z, nnet, j, l, nsteps)
                            end
                        end

                        # -------- FOLLOWER LOAD (Total Lagrange, Piola) ----------
                        if F !== nothing
                            # F elemértékek
                            Fnode = Fmap[elem]   # 9*numNodes × 1

                            # F interpoláció Gauss-pontra (pont úgy, ahogy Kext-ben)
                            Fgp = zeros(DIM, DIM)
                            for a in 1:numNodes
                                Na = h[a, j]
                                Fgp .+= Na * reshape(Fnode[9a-8:9a], DIM, DIM)
                            end

                            Rgp = _rotation_from_F(Fgp)

                            Jgp = det(Fgp)
                            FinvT = inv(Fgp)'

                            if !(Jgp > 0)
                                error("Element inversion detected: det(Fgp) = $Jgp at elem=$elem, gp=$j")
                            end

                            # reference surface normal at GP from reference Jacobian
                            t1 = Jac[:, 3*j-2]          # first tangent in reference
                            t2 = Jac[:, 3*j-1]          # second tangent in reference
                            nref = cross(t1, t2)
                            Ja_ref = norm(nref)         # ugyanaz, mint amit lent Ja-ként használsz
                            N0 = nref / Ja_ref          # unit normal in reference

                            # Nanson scalar: da = |J F^{-T} N0| dA
                            scale = norm(Jgp * FinvT * N0)

                            # Piola transzformáció a traction-re
                            for s in 1:nsteps
                                ##f[:,s] .= Jgp * FinvT * f[:,s]
                                #f[:, s] .= Rgp * f[:, s]          # t = R * t_loc  (itt f a traction, lokálisnak tekint
                                #f[:, s] .= Jgp * FinvT * f[:, s]  # Piola: t0 = J F^{-T} t
                                #f[:, s] .= Rgp * f[:, s]        # true follower: rotate traction with body
                                #f[:, s] .*= scale               # ONLY area scaling (no FinvT acting on t)
                                f[:, s] .= Fgp * f[:, s]        # true follower: rotate traction with body
                            end
                        end
                        # ---------------------------------------------------------

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
                        elseif DIM == 1 && dim == 1
                            Ja = Jac[1, 3*j-2] * b
                        elseif DIM == 1 && dim == 0
                            Ja = 1
                        else
                            error("loadVector: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                        end
                        f1 += H1' * f * Ja * intWeights[j]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                    end
                    fp[nn2,1:nsteps] .+= f1
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
    elseif pdim == 9
        type = :tensor
    else
        error("loadVector: wrong pdim ($pdim).")
    end
    ts = [i for i in 0:nsteps-1]
    if type == :v3D || type == :v2D
        return VectorField([], reshape(fp, :, nsteps), ts, [], nsteps, type, problem)

    elseif type == :scalar
        return ScalarField([], reshape(fp, :, nsteps), ts, [], nsteps, type, problem)

    elseif type == :tensor
        return TensorField([], reshape(fp, :, nsteps), ts, [], nsteps, type, problem)
    end
end
=#
