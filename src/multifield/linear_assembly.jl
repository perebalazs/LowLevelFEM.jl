"""
    matmul_sf(A, B)

Matrix–matrix or matrix–field multiplication used internally
when collapsing operator chains.

Supports combinations of:
- matrices
- scalar/vector/tensor fields
- numeric arrays

Returns the resulting matrix or field.
"""
function matmul_sf(A, B)

    m,k = size(A)
    k2,n = size(B)

    @assert k==k2

    C = Matrix{Any}(undef,m,n)

    for i in 1:m, j in 1:n

        s = nothing

        for p in 1:k

            a = A[i,p]
            b = B[p,j]

            term = a*b

            if s === nothing
                s = term
            else
                s += term
            end

        end

        C[i,j] = s

    end

    return C
end

function _field_pdim(f)
    if f isa ScalarField
        return 1
    elseif f isa VectorField
        if f.type == :v2D
            return 2
        elseif f.type == :v3D
            return 3
        else
            error("Unknown VectorField type $(f.type)")
        end
    elseif f isa TensorField
        return 9
    elseif f isa Number
        return 1
    elseif f isa AbstractVector
        return length(f)
    elseif f isa AbstractMatrix
        return length(f)
    else
        error("Unsupported coefficient type $(typeof(f)) in assemble_linear.")
    end
end

mutable struct LinearWorkspace
    invJac::Matrix{Float64}
    ∂h::Matrix{Float64}
    B::Matrix{Float64}
    fe::Vector{Float64}
end

"""
    LinearWorkspace(dim, num_nodes, num_integration_points, outdim, ndofs_loc)

Allocate reusable element-level work arrays for linear-form assembly.
"""
function LinearWorkspace(
    dim::Int,
    num_nodes::Integer,
    num_integration_points::Integer,
    outdim::Integer,
    ndofs_loc::Integer
)
    return LinearWorkspace(
        zeros(3, 3 * num_integration_points),
        zeros(dim, num_nodes * num_integration_points),
        zeros(outdim, ndofs_loc),
        zeros(ndofs_loc)
    )
end

"""
    assemble_element_vector!(ws, elem, gidx, ...)

Compute the local element vector and store it in `ws.fe`.
"""
function assemble_element_vector!(
    ws::LinearWorkspace,
    elem::Integer,
    gidx::Integer,
    jac_all,
    det_all,
    nip::Integer,
    numIntPoints::Integer,
    dim::Integer,
    numNodes::Integer,
    ∇h,
    h,
    intWeights,
    Gprep,
    Wprep,
    op::AbstractOp,
    P::Problem
)
    invJac = ws.invJac
    ∂h = ws.∂h
    B = ws.B
    fe = ws.fe

    jac_slice =
        @view jac_all[gidx * 9 * nip + 1 : (gidx + 1) * 9 * nip]

    jacDet =
        @view det_all[gidx * nip + 1 : (gidx + 1) * nip]

    Jac = reshape(jac_slice, 3, :)

    @inbounds for k in 1:numIntPoints
        #@views invJac[1:3, 3k-2:3k] .=
        #    inv(Jac[1:3, 3k-2:3k])'
        Jk = SMatrix{3,3,Float64}(Jac[1:3, 3k-2:3k])
        invJac[1:3, 3k-2:3k] .= transpose(inv3x3_fast(Jk))
    end

    fill!(∂h, 0.0)

    @inbounds for k in 1:numIntPoints, a in 1:numNodes
        @views begin
            invJk = invJac[1:dim, 3k-2:3k-(3-dim)]
            gha = ∇h[a*3-2:a*3-(3-dim), k]
            dha = ∂h[1:dim, (k-1)*numNodes+a]
            mul!(dha, invJk, gha)
        end
    end

    fill!(fe, 0.0)

    @inbounds for k in 1:numIntPoints
        wcoef = 1.0

        if Wprep !== nothing
            wcoef = _eval_coefficient_at_gp(
                Wprep, elem, view(h, :, k)
            )
            wcoef isa Number ||
                error("assemble_linear: weight must evaluate to a scalar.")
        end

        w = jacDet[k] * intWeights[k] * wcoef

        build_B!(B, op, P, k, h, ∂h, numNodes)

        g_gp = _eval_coefficient_at_gp(
            Gprep, elem, view(h, :, k)
        )

        if g_gp isa Number
            @assert size(B, 1) == 1
            for i in eachindex(fe)
                fe[i] += B[1, i] * g_gp * w
            end
        else
            mul!(fe, transpose(B), g_gp, w, 1.0)
        end
    end

    return nothing
end

"""
    assemble_linear(P::Problem, op, rhs; weight=nothing, domain, threads=:auto)

Assemble a linear finite element operator of the form

    ∫ v ⋅ op ⋅ rhs

where `v` is the test field associated with problem `P`.

Arguments
---------
- `P::Problem` : finite element problem definition
- `op` : operator or matrix chain acting on the right-hand side
- `rhs` : scalar, vector, tensor field, or numeric vector
- `weight` : optional quadrature weight or coefficient
- `domain` : domain specification (`Ω` or `Γ`)
- `gauss` : Number of integration points.
        `:full`       2order+1
        `:reduced`    2order-1
        Int           2order+1 + Int
- `threads` : number of assembly workers; `:auto` uses all Julia threads and
  `1` selects serial assembly

Returns
-------
Global right-hand side vector.
"""
function assemble_linear(
    P::Problem,
    op::AbstractOp,
    g;
    weight = nothing,
    domain = nothing,
    gauss = :full,
    threads = :auto)

    num_threads = resolve_num_threads(threads)
    use_threads = num_threads > 1

    old_blas_threads = LinearAlgebra.BLAS.get_num_threads()
    if use_threads && old_blas_threads != 1
        LinearAlgebra.BLAS.set_num_threads(1)
    end

    try

    gmsh.model.setCurrent(P.name)

    outdim = op_outdim(op, P)

    # allow matrix coefficient for tensor operators
    if g isa AbstractMatrix
        if length(g) != outdim
            error("matrix coefficient size $(size(g)) incompatible with operator output dimension $outdim")
        end
        g = vec(g)
    end
    if g isa AbstractVector
        if length(g) != outdim
            error("assemble_linear: vector coefficient length ($(length(g))) does not match operator output dimension ($outdim).")
        end
    else
        pdim_g = _field_pdim(g)
        if pdim_g != outdim
            error("assemble_linear: coefficient field dimension ($pdim_g) does not match operator output dimension ($outdim).")
        end
    end

    nd = ndofs(P)
    I = Int[]
    V = Float64[]

    phNames = domain === nothing ?
        [mat.phName for mat in P.material] :
        [domain.name]
    
    for phName in phNames
        Gprep = domain === nothing ?
            _prepare_coefficient(g, DomainSpec(:Ω, phName)) :
            _prepare_coefficient(g, domain)
        Wprep = weight === nothing ? nothing :
            domain === nothing ?
            _prepare_coefficient(weight, DomainSpec(:Ω, phName)) :
            _prepare_coefficient(weight, domain)

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        for (edim, etag) in dimTags

            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            for itype in eachindex(elemTypes)

                et = elemTypes[itype]

                _, _, order, numNodes, _, _ =
                    gmsh.model.mesh.getElementProperties(et)

                gorder = gauss == :reduced ? max(1, 2order - 1) :
                    gauss == :full    ? (2order + 1) :
                    gauss isa Int     ? max(1,2order+1 + gauss) :
                    (2order + 1)

                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(gorder))
                #intPoints, intWeights =
                #    gmsh.model.mesh.getIntegrationPoints(
                #        et, "Gauss" * string(2order + 1)
                #    )

                numIntPoints = length(intWeights)

                _, fun, _ =
                    gmsh.model.mesh.getBasisFunctions(
                        et, intPoints, "Lagrange"
                    )

                h = reshape(fun, :, numIntPoints)

                _, dfun, _ =
                    gmsh.model.mesh.getBasisFunctions(
                        et, intPoints, "GradLagrange"
                    )

                ∇h = reshape(dfun, :, numIntPoints)

                # Get Jacobians for all elements of this type
                elemTags_global, _ = gmsh.model.mesh.getElementsByType(et)

                jac_all, det_all, _ = gmsh.model.mesh.getJacobians(et, intPoints)

                numElementsGlobal = length(elemTags_global)
                nip = length(det_all) ÷ numElementsGlobal

                @assert nip == numIntPoints

                # Map element tag to its position in jac_all and det_all
                idxmap = Dict{eltype(elemTags_global),Int}()
                sizehint!(idxmap, numElementsGlobal)

                @inbounds for (gi, gtag) in enumerate(elemTags_global)
                    idxmap[gtag] = gi - 1
                end

                nel = length(elemTags[itype])
                ndofs_loc = P.pdim * numNodes

                block_start = length(V) + 1
                block_entries = nel * ndofs_loc

                resize!(I, length(I) + block_entries)
                resize!(V, length(V) + block_entries)

                workspaces = [
                    LinearWorkspace(
                        P.dim,
                        numNodes,
                        numIntPoints,
                        outdim,
                        ndofs_loc
                    )
                    for _ in 1:num_threads
                ]

                let Ithread = I, Vthread = V
                    if use_threads
                        _run_workers(num_threads) do worker
                            ws = workspaces[worker]
                            for e in _worker_range(nel, num_threads, worker)
                            elem = elemTags[itype][e]
                            gidx = idxmap[elem]

                            assemble_element_vector!(
                                ws,
                                elem,
                                gidx,
                                jac_all,
                                det_all,
                                nip,
                                numIntPoints,
                                P.dim,
                                numNodes,
                                ∇h,
                                h,
                                intWeights,
                                Gprep,
                                Wprep,
                                op,
                                P
                            )

                            element_start =
                                block_start + (e - 1) * ndofs_loc

                            @inbounds for a_loc in 1:ndofs_loc
                                node = div(a_loc - 1, P.pdim) + 1
                                comp = mod(a_loc - 1, P.pdim) + 1

                                gnode =
                                    elemNodeTags[itype][
                                        (e - 1) * numNodes + node
                                    ]

                                p = element_start + a_loc - 1
                                Ithread[p] =
                                    (gnode - 1) * P.pdim + comp
                                Vthread[p] = ws.fe[a_loc]
                            end
                            end
                        end
                    else
                        @inbounds for e in 1:nel
                            ws = workspaces[1]
                            elem = elemTags[itype][e]
                            gidx = idxmap[elem]

                            assemble_element_vector!(
                                ws,
                                elem,
                                gidx,
                                jac_all,
                                det_all,
                                nip,
                                numIntPoints,
                                P.dim,
                                numNodes,
                                ∇h,
                                h,
                                intWeights,
                                Gprep,
                                Wprep,
                                op,
                                P
                            )

                            element_start =
                                block_start + (e - 1) * ndofs_loc

                            for a_loc in 1:ndofs_loc
                                node = div(a_loc - 1, P.pdim) + 1
                                comp = mod(a_loc - 1, P.pdim) + 1

                                gnode =
                                    elemNodeTags[itype][
                                        (e - 1) * numNodes + node
                                    ]

                                p = element_start + a_loc - 1
                                Ithread[p] =
                                    (gnode - 1) * P.pdim + comp
                                Vthread[p] = ws.fe[a_loc]
                            end
                        end
                    end
                end
            end
        end
    end

    rhs = collect(sparsevec(I, V, nd))

    if P.pdim == 1
        return ScalarField([], reshape(rhs, :, 1), [0], [], 1, :scalar, P)

    elseif P.pdim == 2 || P.pdim == 3
        type = P.pdim == 2 ? :v2D : :v3D
        return VectorField([], reshape(rhs, :, 1), [0], [], 1, type, P)

    elseif P.pdim == 9
        return TensorField([], reshape(rhs, :, 1), [0], [], 1, :tensor, P)

    else
        error("assemble_linear: unsupported pdim $(P.pdim).")
    end
    finally
        if use_threads && old_blas_threads != 1
            LinearAlgebra.BLAS.set_num_threads(old_blas_threads)
        end
    end
end

"""
    compliance9_iso(E, nu; penalty=1e8)

Return a 9×9 compliance-like matrix acting on vec(σ) with ordering
(11,12,13,21,22,23,31,32,33). Symmetric part follows isotropic linear elasticity,
antisymmetric part is penalized with `penalty`.
"""
function compliance9_iso(E, nu; penalty=1e8)
    G = E / (2 * (1 + nu))

    # Voigt compliance (engineering shear):
    # [ε11, ε22, ε33, γ12, γ13, γ23] = S6 * [σ11, σ22, σ33, σ12, σ13, σ23]
    S6 = zeros(6, 6)
    S6[1, 1] = 1 / E
    S6[2, 2] = 1 / E
    S6[3, 3] = 1 / E
    S6[1, 2] = -nu / E
    S6[1, 3] = -nu / E
    S6[2, 1] = -nu / E
    S6[2, 3] = -nu / E
    S6[3, 1] = -nu / E
    S6[3, 2] = -nu / E
    S6[4, 4] = 1 / G
    S6[5, 5] = 1 / G
    S6[6, 6] = 1 / G

    # Map 9 -> 6 (take symmetric + engineering shear)
    # v9 = [σ11 σ12 σ13 σ21 σ22 σ23 σ31 σ32 σ33]
    # v6 = [σ11, σ22, σ33, σ12, σ13, σ23] with σ12 := (σ12+σ21)/2 etc.
    P = zeros(6, 9)
    # normals
    P[1, 1] = 1.0       # σ11
    P[2, 5] = 1.0       # σ22
    P[3, 9] = 1.0       # σ33
    # shear sym averages
    P[4, 2] = 0.5
    P[4, 4] = 0.5   # σ12
    P[5, 3] = 0.5
    P[5, 7] = 0.5   # σ13
    P[6, 6] = 0.5
    P[6, 8] = 0.5   # σ23

    # Sym part compliance in 9-space: Ssym = P' * S6 * P
    Ssym = P' * S6 * P

    # Antisym penalty: penalize (σ12-σ21), (σ13-σ31), (σ23-σ32)
    # Add penalty * Q'Q where Q*v9 = [σ12-σ21, σ13-σ31, σ23-σ32]
    Q = zeros(3, 9)
    Q[1, 2] = 1.0
    Q[1, 4] = -1.0
    Q[2, 3] = 1.0
    Q[2, 7] = -1.0
    Q[3, 6] = 1.0
    Q[3, 8] = -1.0
    Santi = penalty * (Q' * Q)

    return Ssym + Santi
end

