```julia
using LinearAlgebra, SparseArrays, StaticArrays, Polyester
using Gmsh

# Ez a segédfüggvény egyetlen véges elem merevségi mátrixát és a
# globális indexeit számítja ki. A fő függvény a "stiffnessMatrixSolid7" ezt
# hívja meg párhuzamosan.
# A függvény minden szükséges paramétert argumentumként kap, így nem függ a külső
# hatókörben lévő változóktól, ami a párhuzamos futtatás stabilitását garantálja.
function calculate_element_stiffness(elem_tag, pdim, numNodes_elem, D, rowsOfB, b, intPoints_raw, intWeights, ∇h_local_grads, elemNodeTags_phase, i, j, dim)
    
    # A szabadságfokok (DOFs) száma elemenként
    num_dofs_per_element = pdim * numNodes_elem
    
    # Az aktuális elemhez tartozó csomópontok globális azonosítóinak lekérdezése.
    # A vector használata lehetővé teszi a dinamikus méretezést.
    element_node_tags = Vector{Int}(undef, numNodes_elem)
    for k in 1:numNodes_elem
        element_node_tags[k] = elemNodeTags_phase[i][(j-1)*numNodes_elem+k]
    end

    # Az elemhez tartozó szabadságfokok (DOFs) globális azonosítóinak generálása.
    # Pl. egy 2D-s elemnél (pdim=2), ha a csomópont ID 10, akkor a DOFs 19 és 20.
    global_dofs_for_element = Vector{Int}(undef, num_dofs_per_element)
    for k_node in 1:numNodes_elem
        for l_dof in 1:pdim
            global_dofs_for_element[(k_node-1)*pdim + l_dof] = pdim * element_node_tags[k_node] - (pdim - l_dof)
        end
    end
    
    # A Jacobian mátrixok és determinánsok lekérdezése a Gmsh API-ból.
    jac_raw, jacDet_raw, coord_raw = gmsh.model.mesh.getJacobian(elem_tag, intPoints_raw)
    
    # Az inverz Jacobi mátrixok tárolása SMatrix formátumban, ami optimalizált a
    # kisebb, statikus méretű mátrixokhoz.
    invJac_matrices = Vector{SMatrix{dim, dim, Float64}}(undef, length(intWeights))
    for k_int in 1:length(intWeights)
        Jac_k_full = SMatrix{3, 3, Float64}(jac_raw[9*(k_int-1)+1 : 9*k_int])
        invJac_matrices[k_int] = inv(Jac_k_full[1:dim, 1:dim])'
    end
    
    # Elemi merevségi mátrix inicializálása (K_e)
    # Ez a mátrix tárolja az aktuális elem hozzájárulását a globális merevségi mátrixhoz.
    K_element_local = zeros(num_dofs_per_element, num_dofs_per_element)
    
    # B mátrix (alakváltozás-elmozdulás mátrix) inicializálása
    B_int_point = zeros(rowsOfB, num_dofs_per_element)

    # Integrációs pontokon végigfutó ciklus
    for k_int in 1:length(intWeights)
        invJac_k = invJac_matrices[k_int]
        jacDet_k = jacDet_raw[k_int]
        intWeight_k = intWeights[k_int]
        
        # A B-mátrix nullázása minden integrációs ponthoz.
        fill!(B_int_point, 0.0)
        
        # A csomópontokon végigfutó ciklus a B-mátrix feltöltéséhez.
        for l_node in 1:numNodes_elem
            # A globális gradiens (alakfüggvény deriváltja) számítása
            dN_global_vec = invJac_k * ∇h_local_grads[k_int][l_node]

            # A B-mátrix feltöltése a probléma dimenziójától (2D/3D) függően.
            if dim == 2 && rowsOfB == 3
                B_int_point[1, (l_node-1)*pdim + 1] = dN_global_vec[1]
                B_int_point[2, (l_node-1)*pdim + 2] = dN_global_vec[2]
                B_int_point[3, (l_node-1)*pdim + 1] = dN_global_vec[2]
                B_int_point[3, (l_node-1)*pdim + 2] = dN_global_vec[1]
            elseif dim == 3 && rowsOfB == 6
                B_int_point[1, (l_node-1)*pdim+1] = dN_global_vec[1]
                B_int_point[2, (l_node-1)*pdim+2] = dN_global_vec[2]
                B_int_point[3, (l_node-1)*pdim+3] = dN_global_vec[3]
                B_int_point[4, (l_node-1)*pdim+1] = dN_global_vec[2]
                B_int_point[4, (l_node-1)*pdim+2] = dN_global_vec[1]
                B_int_point[5, (l_node-1)*pdim+2] = dN_global_vec[3]
                B_int_point[5, (l_node-1)*pdim+3] = dN_global_vec[2]
                B_int_point[6, (l_node-1)*pdim+1] = dN_global_vec[3]
                B_int_point[6, (l_node-1)*pdim+3] = dN_global_vec[1]
            else
                error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
            end
        end
        # Az elemi merevségi mátrix frissítése az integrációs pontok alapján.
        # K_e = B_t * D * B * det(J) * w
        K_element_local += B_int_point' * D * B_int_point * b * jacDet_k * intWeight_k
    end
    
    # A globális I és J indexvektorok előkészítése a ritka mátrixhoz.
    global_row_dofs = repeat(global_dofs_for_element, 1, num_dofs_per_element)[:]
    global_col_dofs = repeat(global_dofs_for_element', num_dofs_per_element, 1)[:]

    return global_row_dofs, global_col_dofs, K_element_local[:]
end

# Ez a fő függvény felelős a globális merevségi mátrix (K) összeszereléséért.
# Párhuzamosítja az elemi mátrixok számítását a Polyester csomag segítségével.
function stiffnessMatrixSolid7(problem; elements=Int[])
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(problem.dim, -1)
    
    # A globális I, J, V tömbök méretének előzetes becslése a memóriafoglalás optimalizálására.
    lengthOfIJV = 0
    for i in 1:length(elemTags)
        if !isempty(elemTags[i])
            nodes_per_elem = div(length(elemNodeTags[i]), length(elemTags[i]))
            dofs_per_elem = nodes_per_elem * problem.dim
            lengthOfIJV += dofs_per_elem^2 * length(elemTags[i])
        end
    end

    I_global = Vector{Int}(undef, 0)
    J_global = Vector{Int}(undef, 0)
    V_global = Vector{Float64}(undef, 0)
    sizehint!(I_global, lengthOfIJV)
    sizehint!(J_global, lengthOfIJV)
    sizehint!(V_global, lengthOfIJV)

    # A fizikai anyagcsoportokon végigfutó ciklus.
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        E = problem.material[ipg].E
        ν = problem.material[ipg].ν
        dim = problem.dim
        pdim = problem.dim

        # Az anyagjellemző D mátrix (stress-strain) definiálása a probléma típusától függően.
        local D::SMatrix
        local rowsOfB::Int
        local b::Float64

        if dim == 3 && problem.type == :Solid
            D_base = @SMatrix [1-ν ν ν 0 0 0; ν 1-ν ν 0 0 0; ν ν 1-ν 0 0 0; 0 0 0 (1-2ν)/2 0 0; 0 0 0 0 (1-2ν)/2 0; 0 0 0 0 0 (1-2ν)/2]
            D = E / ((1 + ν) * (1 - 2ν)) * D_base
            rowsOfB = 6
            b = 1.0
        elseif dim == 2 && problem.type == :PlaneStress
            D_base = @SMatrix [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
            D = E / (1 - ν^2) * D_base
            rowsOfB = 3
            b = problem.thickness
        elseif dim == 2 && problem.type == :PlaneStrain
            D_base = @SMatrix [1-ν ν 0; ν 1-ν 0; 0 0 (1-2ν)/2]
            D = E / ((1 + ν) * (1 - 2ν)) * D_base
            rowsOfB = 3
            b = 1.0
        else
            error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for idm in 1:length(dimTags)
            dimTag = dimTags[idm]
            edim = dimTag[1]
            etag = dimTag[2]

            elemTypes_phase, elemTags_phase, elemNodeTags_phase = gmsh.model.mesh.getElements(edim, etag)

            for i in 1:length(elemTypes_phase)
                et = elemTypes_phase[i]
                elementName, dim_elem, order, numNodes_elem::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
                intPoints_raw, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order + 1))
                numIntPoints = length(intWeights)
                comp, dfun_raw, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints_raw, "GradLagrange")
                
                ∇h_local_grads = Vector{Vector{SVector{dim, Float64}}}(undef, numIntPoints)
                for k_int in 1:numIntPoints
                    ∇h_local_grads[k_int] = Vector{SVector{dim, Float64}}(undef, numNodes_elem)
                    for l_node in 1:numNodes_elem
                        idx_start = (k_int - 1) * (numNodes_elem * dim) + (l_node - 1) * dim + 1
                        ∇h_local_grads[k_int][l_node] = SVector{dim, Float64}(dfun_raw[idx_start : idx_start + dim - 1])
                    end
                end

                # Thread-local tömbök inicializálása a párhuzamos adatgyűjtéshez.
                all_I_threads = [Vector{Int}() for _ in 1:Threads.nthreads()]
                all_J_threads = [Vector{Int}() for _ in 1:Threads.nthreads()]
                all_V_threads = [Vector{Float64}() for _ in 1:Threads.nthreads()]

                # Párhuzamos ciklus a Polyester.@batch makróval.
                # A makró szálakra osztja a számításokat.
                @batch for j in 1:length(elemTags_phase[i])
                    tid = Threads.threadid()
                    # A segédfüggvény hívása az elemi merevségi mátrix kiszámítására.
                    I_local, J_local, V_local = calculate_element_stiffness(elemTags_phase[i][j], pdim, numNodes_elem, D, rowsOfB, b, intPoints_raw, intWeights, ∇h_local_grads, elemNodeTags_phase, i, j, dim)
                    # Hozzáadás a szálhoz tartozó lokális tömbökhöz.
                    append!(all_I_threads[tid], I_local)
                    append!(all_J_threads[tid], J_local)
                    append!(all_V_threads[tid], V_local)
                end

                # A szálak által generált eredmények egyesítése a globális tömbökbe.
                append!(I_global, vcat(all_I_threads...))
                append!(J_global, vcat(all_J_threads...))
                append!(V_global, vcat(all_V_threads...))
            end
        end
    end
    
    # A végső globális ritka mátrix (K) összeszerelése az összegyűjtött adatokból.
    dof = problem.dim * problem.non
    K = sparse(I_global, J_global, V_global, dof, dof)
    dropzeros!(K)
    return SystemMatrix(K, problem)
end
```