function copy(A::AbstractField)
    T = typeof(A)
    
    a = [copy(block) for block in A.A]
    b = copy(A.a)
    c = copy(A.t)
    d = copy(A.numElem)
    e = copy(A.nsteps)
    f = A.type
    g = A.model
    return T(a, b, c, d, e, f, g, A.kind)
end

import Base:show
"""
    Base.show(io::IO, ::MIME"text/plain", M::Union{ScalarField,VectorField,TensorField})

Plain-text display in REPL and notebooks.
"""
function Base.show(io::IO, ::MIME"text/plain", M::Union{ScalarField,VectorField,TensorField})
    # fejléc külön sorban (szebb notebookban)
    if isNodal(M)
        type = "nodal "
    else
        type = "elementwise "
    end
    
    println(io, type, nameof(typeof(M)))
    
    io2 = IOContext(io, :compact => false)
    
    if isNodal(M)
        Base.show(io2, M.a)
    else
        Base.show(io2, M.A)
    end
end

import Base:display
function display(M::Union{ScalarField,VectorField,TensorField})
    if isNodal(M)
        display("$(nameof(typeof(M))) $(M.a)")
    else
        display("$(nameof(typeof(M))) $(M.A)")
    end
end

"""
    ∂(r::ScalarField, dir::Int) -> ScalarField

Compute the spatial derivative of a scalar field in the global `x`, `y`, or `z`
direction. The input field may be nodal or elementwise; nodal fields are
automatically converted to elementwise form using `nodesToElements(r)`.

# Arguments
- `r::ScalarField`: scalar field to differentiate.
- `dir::Int`: spatial direction:
    * `1` → ∂/∂x  
    * `2` → ∂/∂y  
    * `3` → ∂/∂z

# Returns
A new `ScalarField` in **elementwise representation**, containing  
the gradient component ∂r/∂(dir) at all nodes of each element.

# Notes
- Works for 1D, 2D, and 3D meshes.
- Uses the element Jacobian and shape-function derivatives in global coordinates.
- The global problem mesh is accessed via `r.model`.
"""
function ∂(r::ScalarField, dir::Int)
    #@info "∂: nodal"
    problem = r.model
    gmsh.model.setCurrent(problem.name)

    @assert dir ∈ 1:3 "d(r,dir): dir must be 1 (x), 2 (y), or 3 (z)."

    # ha nodális tér, konvertáljuk elemisre
    #r = rr.a == [;;] ? elementsToNodes(rr) : rr
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
            if r.model.geometry.nameVolume ≠ ""
                phName = r.model.geometry.nameVolume
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

                # munkatömbök
                invJac = Matrix{Float64}(undef, 3, 3*numNodes)
                ∂h    = Matrix{Float64}(undef, 1, numNodes*numNodes)  # ∂/∂dir
                nn2   = Vector{Int}(undef, numNodes)
                nnet  = Matrix{Int}(undef, length(elemTags[it]), numNodes)

                dim_et = 3
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
                        invJk = pinv(Matrix(Jk'))       # stabil

                        fill!(@view(invJac[1:3, 3*k-2:3*k]), 0.0)
                        @views invJac[1:dim_et, 3*k-2:3*k-3+dim_et] .= invJk
                    end

                    # --- ∂h = ∂/∂dir (shape-derivált) ---
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numNodes, l in 1:numNodes
                        g = @view ∇h_all[(l-1)*3+1 : (l-1)*3+dim_et, k]
                        Jslice = @view invJac[:, 3*k-2 : 3*k-2+dim_et-1]

                        dphi = (Jslice * g)[dir]   # dir=1,2,3 → x,y,z irány
                        ∂h[1, (k-1)*numNodes + l] = dphi
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
    ∂e(r::ScalarField, dir::Int) -> ScalarField

Compute the spatial derivative of an **elementwise** scalar field purely at the
element level, without nodal conversion. This avoids mixing contributions from
neighboring elements and is suitable for discontinuous fields or post-processing
of elementwise quantities.

# Arguments
- `r::ScalarField`: must already be elementwise (`isElementwise(r) == true`).
- `dir::Int`: spatial direction (1 → x, 2 → y, 3 → z).

# Returns
An elementwise `ScalarField` containing the derivative ∂r/∂(dir) at the local
nodes of each element.

# Notes
- Uses the element-local Jacobian per integration point.
- Safer for interfaces or discontinuities than nodal differentiation.
- Only elementwise → elementwise transformations are performed.
"""
function ∂e(rr::ScalarField, dir::Int)
    #@info "∂e: elementwise scalar derivative"
    problem = rr.model
    gmsh.model.setCurrent(problem.name)

    dim = problem.dim              # 1, 2 vagy 3
    @assert 1 ≤ dir ≤ dim "∂e(rr,dir): dir must be between 1 and problem.dim."

    # elementwise scalar field kell
    @assert rr.A != [] "∂e: rr must be elementwise ScalarField (A ≠ [])."

    nsteps = rr.nsteps

    # elemTag → index in rr.A / rr.numElem
    elem_index = Dict{Int,Int}()
    @inbounds for (i, e) in enumerate(rr.numElem)
        elem_index[e] = i
    end

    ε       = Vector{Matrix{Float64}}()  # elementwise results
    numElem = Int[]                      # elemTags az eredményhez

    # --- iterate physical groups (ugyanúgy, mint a többi függvényben) ---
    for ipg in 0:length(problem.material)
        phName =
            ipg == 0 ?
                (rr.model.geometry.nameVolume ≠ "" ? rr.model.geometry.nameVolume : "") :
                problem.material[ipg].phName

        if phName == ""
            continue
        end

        dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

        @inbounds for (edim, etag) in dimTags
            elemTypes, elemTags, elemNodeTags =
                gmsh.model.mesh.getElements(edim, etag)

            @inbounds for it in 1:length(elemTypes)
                et = Int(elemTypes[it])

                # cached FE data
                dim_et, numNodes, nodeCoord = _get_props_cached(et)
                ∇h_all, h_all = _get_basis_cached(et, nodeCoord, numNodes)

                # csak olyan elemtípussal foglalkozzunk, ahol az rr-ben van elem
                # (elem_index úgyis szűr, itt nem muszáj, de nem árt)

                # invJac: dim × (dim_et * numNodes) – minden csomóponthoz egy dim×dim_et blokk
                invJac = Matrix{Float64}(undef, dim, dim_et * numNodes)
                ∂h     = Matrix{Float64}(undef, 1, numNodes * numNodes)   # csak a kiválasztott dir komponens

                # elem ciklus
                @inbounds for j in 1:length(elemTags[it])
                    elem = elemTags[it][j]

                    # csak azok az elemek, amelyekhez rr-ben tényleg van adat
                    eidx = get(elem_index, elem, 0)
                    if eidx == 0
                        continue
                    end

                    # Jacobianok az elem csomópontjain
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
                    Jac = reshape(jac, 3, :)   # gmsh mindig 3 sorral jön, de csak 1..dim kell

                    # invJac blokkok minden lokális csomópontra
                    @inbounds for k in 1:numNodes
                        # Jk: dim × dim_et
                        Jk = @view Jac[1:dim, (k-1)*3+1 : (k-1)*3+dim_et]
                        invJk = pinv(Matrix(Jk'))   # dim_et × dim  →  dim × dim_et

                        @views invJac[:, (k-1)*dim_et+1 : k*dim_et] .= invJk
                    end

                    # ∂h = ∂φ_l/∂x_dir globális koordinátákban
                    fill!(∂h, 0.0)
                    @inbounds for k in 1:numNodes, l in 1:numNodes
                        # lokális grad (dim_et komponens)
                        g = @view ∇h_all[(l-1)*3+1 : (l-1)*3+dim_et, k]

                        # dim × dim_et blokk a k-adik csomópontnál
                        Jslice = @view invJac[:, (k-1)*dim_et+1 : k*dim_et]

                        # globális grad (dim komponens) → ebből választjuk a dir-ediket
                        dphi_vec = Jslice * g
                        dphi = dphi_vec[dir]

                        ∂h[1, (k-1)*numNodes + l] = dphi
                    end

                    # csomóponti értékek ezen az elemen (numNodes × nsteps)
                    vals = rr.A[eidx]
                    @assert size(vals, 1) == numNodes "∂e: numNodes mismatch with rr.A."

                    # elem eredmény: numNodes × nsteps
                    e = Matrix{Float64}(undef, numNodes, nsteps)

                    # ∑_l ∂φ_l * value_l
                    @inbounds for k in 1:numNodes
                        idx = (k-1)*numNodes
                        for kk in 1:nsteps
                            s = 0.0
                            for l in 1:numNodes
                                s += ∂h[1, idx + l] * vals[l, kk]
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

    return ScalarField(ε, [;;], rr.t, numElem, nsteps, :scalar, problem)
end

"""
    ∂x(r::ScalarField) -> ScalarField

Differentiate a scalar field with respect to global Cartesian `x`.
Equivalent to `∇(r)[1]`; see [`∇`](@ref) for supported geometries.

Nodal input is evaluated through element interpolation. Element-wise input
preserves its support and discontinuities. The output is element-wise and
contains all input columns. The derivative is that of the finite-element
interpolant, so an analytic function must be representable by the element
space before its exact derivative can be expected.

```julia
s = ScalarField(P, "body", (x,y,z) -> 2x + 3y)
dx = ∂x(s)  # constant 2 on a mesh reproducing linear functions
```
"""
∂x(r::ScalarField) = ∇(r)[1]
"""
    ∂y(r::ScalarField) -> ScalarField

Differentiate with respect to global Cartesian `y`. Equivalent to
`∇(r)[2]`; input support, columns, and discontinuities are handled as
for [`∂x`](@ref). The result is element-wise.

On an xy-planar mesh, `∂z` returns zero. On embedded curves or surfaces,
Cartesian components follow the tangential physical gradient from [`∇`](@ref).
"""
∂y(r::ScalarField) = ∇(r)[2]
"""
    ∂z(r::ScalarField) -> ScalarField

Differentiate with respect to global Cartesian `z`. Equivalent to
`∇(r)[3]`; input support, columns, and discontinuities are handled as
for [`∂x`](@ref). The result is element-wise.

On an xy-planar mesh, `∂z` returns zero. On embedded curves or surfaces,
Cartesian components follow the tangential physical gradient from [`∇`](@ref).
"""
∂z(r::ScalarField) = ∇(r)[3]
#∂x(r::ScalarField) = isNodal(r) ? ∂(r, 1) : ∂e(r, 1)
#∂y(r::ScalarField) = isNodal(r) ? ∂(r, 2) : ∂e(r, 2)
#∂z(r::ScalarField) = isNodal(r) ? ∂(r, 3) : ∂e(r, 3)

"""
    ∂t(s::Union{ScalarField,VectorField,TensorField})

Compute the time derivative of a time-dependent field using second-order
finite differences in time.

The field values are assumed to be stored at discrete time instants `s.t`
and are differentiated independently at each spatial degree of freedom.
The result has the same type and layout (nodal or elementwise) as the input
field.

### Time discretization
- For interior time steps, a second-order central difference is used:
```

∂s/∂t(tᵢ) ≈ (s(tᵢ₊₁) − s(tᵢ₋₁)) / (tᵢ₊₁ − tᵢ₋₁)

```
- At the first and last time steps, the derivative is obtained by
second-order extrapolation to preserve overall accuracy.
- If only two time steps are present, a first-order difference is used
and assigned to both time levels.
- If `s.nsteps == 1`, a zero field is returned.

The time grid `s.t` is not required to be uniform.

### Notes
- The operation is purely algebraic in time and does not involve any
spatial operators or time-integration schemes.
- The derivative is computed in-place on newly allocated storage and
does not modify the original field.
- For elementwise fields, the time derivative is computed independently
for each element.

### Returns
A new field of the same concrete type as `s` (`ScalarField`, `VectorField`,
or `TensorField`) containing the time derivative.
"""
function ∂t(s::Union{ScalarField,VectorField,TensorField})
    if s.nsteps == 1
        return 0s
    end
    T = typeof(s)
    if isNodal(s)
        a = similar(s.a)
        time_derivative!(a, s.t, s.a)
        return T([], a, s.t, [], s.nsteps, s.type, s.model)
    elseif isElementwise(s)
        A = [similar(s.A[e]) for e in eachindex(s.A)]
        for e in eachindex(s.A)
            time_derivative!(A[e], s.t, s.A[e])
        end
        return T(A, [;;], s.t, s.numElem, s.nsteps, s.type, s.model)
    else
        error("∂t: internal error")
    end
end

function time_derivative!(a, t, src)
    nsteps = length(t)
    if nsteps == 2
        @inbounds begin
            a[:,1] .= (src[:,2] - src[:,1]) / (t[2] - t[1])
            a[:,2] .= a[:,1]
        end
    else
        @inbounds begin
            for i in 2:nsteps-1
                a[:,i] .= (src[:,i+1] - src[:,i-1]) / (t[i+1] - t[i-1])
            end
            a[:,1]   .= 2*(src[:,2]-src[:,1])/(t[2]-t[1]) - a[:,2]
            a[:,end] .= 2*(src[:,end]-src[:,end-1])/(t[end]-t[end-1]) - a[:,end-1]
        end
    end
end

"""
    Base.getindex(v::VectorField, k::Int)

    v[k] -> ScalarField

Extract the k-th component (1,2,3) of a `VectorField` as a `ScalarField`.

# Arguments
- `k::Int`: component index (1 → x, 2 → y, 3 → z).

# Returns
A `ScalarField` in the same representation as the parent field (elementwise or nodal).

# Errors
Raises an error for any non-scalar indexing:
- `v[1:2]`    → ERROR  
- `v[[1,3]]`  → ERROR  
- `v[:]`      → ERROR

Only single-component extraction `v[k]` is allowed.
The extracted component is marked `kind=:field`, since it is not a
complete eigenvector of the original vector problem.
"""
function Base.getindex(v::VectorField, k::Int)
    pdim =
        v.type == :v2D ? 2 :
        v.type == :v3D ? 3 :
        error("Unknown VectorField type $(v.type)")

    @assert 1 ≤ k ≤ pdim "VectorField has exactly $pdim components."

    if isElementwise(v)
        Aout = Vector{Matrix{Float64}}(undef, length(v.A))

        @inbounds for e in eachindex(v.A)
            Ae = v.A[e]

            size(Ae, 1) % pdim == 0 ||
                error(
                    "VectorField: incompatible local data size " *
                    "for element $(v.numElem[e])."
                )

            Aout[e] = Ae[k:pdim:end, :]
        end

        return ScalarField(
            Aout,
            [;;],
            v.t,
            v.numElem,
            v.nsteps,
            :scalar,
            v.model,
            :field
        )
    end

    if isNodal(v)
        return ScalarField(
            [],
            v.a[k:pdim:end, :],
            v.t,
            v.numElem,
            v.nsteps,
            :scalar,
            v.model,
            :field
        )
    end

    error("VectorField: no data found in A or a.")
end

#=
function Base.getindex(v::VectorField, k::Int)
    pdim = v.model.pdim
    pdim =
        v.type == :v2D ? 2 :
        v.type == :v3D ? 3 :
        error("Unknown VectorField type $(v.type)")
    @assert 1 ≤ k ≤ pdim "VectorField has exactly $pdim components."

    # ha elementwise mező
    if v.A != []
        numElem = length(v.A)
        numNodes = size(v.A[1], 1) ÷ pdim   # 3 komponens miatt
        nsteps = v.nsteps

        Aout = Vector{Matrix{Float64}}(undef, numElem)

        @inbounds for i in 1:numElem
            Ai = v.A[i]
            rows = k:pdim:pdim*numNodes
            Aout[i] = Ai[rows, :]
        end

        return ScalarField(Aout, [;;], v.t, v.numElem, v.nsteps, :scalar, v.model)
    end

    # ha nodal mező
    if v.a != [;;]
        numNodes = size(v.a, 1) ÷ pdim
        rows = k:pdim:pdim*numNodes
        aout = v.a[rows, :]

        return ScalarField([], aout, v.t, v.numElem, v.nsteps, :scalar, v.model)
    end

    error("VectorField: no data found in A or a.")
end
=#

"""
    Base.getindex(T::TensorField, i::Int, j::Int)

    T[i,j] -> ScalarField

Extract a single tensor component (i,j) from a 3×3 `TensorField`.

# Arguments
- `i::Int`, `j::Int`: tensor indices in 1:3.

# Returns
A `ScalarField` corresponding to component Tᵢⱼ.

# Notes
Block extraction is done elementwise or nodally depending on the representation
of `T`.

# Errors
Raises an error if i or j is outside 1:3.
"""
function Base.getindex(T::TensorField, i::Int, j::Int)
    @assert 1 ≤ i ≤ 3 && 1 ≤ j ≤ 3 "TensorField indexing must be T[i,j], with i,j ∈ 1:3"

    block = (j-1)*3 + i   # which tensor component
    # --- elementwise ---
    if T.A != []
        numElem  = length(T.A)
        #numNodes = size(T.A[1], 1) ÷ 9     # 9 tensor components
        #Aout = Vector{Matrix{Float64}}(undef, numElem)

        #@inbounds for e in 1:numElem
        #    Ae = T.A[e]
        #    rows = block:9:9*numNodes
        #    Aout[e] = Ae[rows, :]
        #end
        Aout = Vector{Matrix{Float64}}(undef, length(T.A))

        @inbounds for e in eachindex(T.A)
            Ae = T.A[e]

            size(Ae, 1) % 9 == 0 ||
                error(
                    "TensorField: incompatible local data size " *
                    "for element $(T.numElem[e])."
                )

            Aout[e] = Ae[block:9:end, :]
        end

        return ScalarField(Aout, [;;], T.t, T.numElem, T.nsteps, :scalar, T.model, T.kind)
    end

    # --- nodal ---
    if T.a != [;;]
        numNodes = size(T.a, 1) ÷ 9
        rows = block:9:9*numNodes
        aout = T.a[rows, :]

        return ScalarField([], aout, T.t, T.numElem, T.nsteps, :scalar, T.model, T.kind)
    end

    error("TensorField: no data in A or a.")
end

"""
    Base.getindex(T::TensorField, i::Int, ::Colon)

    T[i,:] -> VectorField

Extract the i-th row of a 3×3 `TensorField` as a 3-component `VectorField`.

Equivalent to: [T[i,1], T[i,2], T[i,3]]

# Arguments
- `i::Int`: row index (1 ≤ i ≤ 3).

# Returns
A `VectorField` of size 3.

# Notes
Colon indexing is required; `T[i]` is not allowed.
"""
function Base.getindex(T::TensorField, i::Int, ::Colon)
    @assert 1 ≤ i ≤ 3 "TensorField row index must be 1..3"

    comps = ScalarField[]
    for j in 1:3
        push!(comps, T[i,j])   # already returns ScalarField
    end
    return VectorField(comps)
end

"""
    Base.getindex(T::TensorField, ::Colon, j::Int)

    T[:,j] -> VectorField

Extract the j-th column of a 3×3 `TensorField` as a 3-component `VectorField`.

Equivalent to: [T[1,j], T[2,j], T[3,j]]

# Arguments
- `j::Int`: column index (1 ≤ j ≤ 3).

# Returns
A `VectorField` of size 3.

# Notes
Colon indexing is required; `T[j]` is not allowed.
"""
function Base.getindex(T::TensorField, ::Colon, j::Int)
    @assert 1 ≤ j ≤ 3 "TensorField column index must be 1..3"

    comps = ScalarField[]
    for i in 1:3
        push!(comps, T[i,j])
    end
    return VectorField(comps)
end

"""
    Base.getindex(T::TensorField, k::Int)

T[k] -> ScalarField

Return the k-th component of the tensor in column-major order:

    [T11, T21, T31, T12, T22, T32, T13, T23, T33]

# Arguments
- `k::Int`: component index (1 ≤ k ≤ 9)

# Returns
ScalarField
"""
function Base.getindex(T::TensorField, k::Int)
    @assert 1 ≤ k ≤ 9 "TensorField index must be 1..9"

    i = mod1(k, 3)
    j = (k-1) ÷ 3 + 1

    return T[i,j]
end

"""
    Base.getindex(T::TensorField, I...)

Invalid index pattern for `TensorField`.

Allowed:
- `T[i,j]`    → ScalarField
- `T[i,:]`    → VectorField (row)
- `T[:,j]`    → VectorField (column)

Not allowed:
- `T[i]`
- `T[:]`
- `T[1:2,1]`
- `T[:, :]`
- `T[[1,3],2]`
"""
function Base.getindex(T::TensorField, I...)
    error("""
    Invalid TensorField indexing.

    Allowed:
      T[i,j]    → ScalarField
      T[i,:]    → VectorField (row)
      T[:,j]    → VectorField (column)

    Everything else is invalid.
    """)
end

"""
    Base.setindex!(v::VectorField, s::ScalarField, k::Int)

    v[k] = s

Assign a scalar component into a `VectorField`.

This operation overwrites the k-th component (k = 1, 2, 3) of the vector field.

## Requirements
- `v` must be a `VectorField`
- `s` must be a `ScalarField`
- both fields must belong to the **same problem**
- both must have the **same number of time steps**
- storage layout must match:
    - If `v` is **elementwise** (`v.A ≠ []`), then `s` must also be elementwise
      and must have the same number of elements (`numElem`).
    - If `v` is **nodal** (`v.a ≠ [;;]`), then `s` must also be nodal and must have
      the same nodal size.

Only the rows belonging to component `k` are modified.  
Other components remain unchanged.

## Example
```julia
v = VectorField([s1, s2, s3])

v[1] = s4       # overwrite first component
v[3] = 2s1      # scale and assign into third component
```

## Errors

An informative error is thrown if:

* k ∉ 1:3
* the fields belong to different problems
* storage layouts do not match (nodal vs. elementwise)
* nodal or element counts differ
"""
function Base.setindex!(v::VectorField, s::ScalarField, k::Int)
    @assert 1 ≤ k ≤ 3 "VectorField has exactly 3 components."

    # check problem consistency
    @assert s.model === v.model "ScalarField and VectorField belong to different problems."
    @assert s.nsteps == v.nsteps "Time step mismatch in VectorField assignment."

    # elementwise case
    if v.A != []
        @assert s.A != [] "Cannot assign nodal field into elementwise VectorField."
        @assert length(s.A) == length(v.A) "numElem mismatch in VectorField assignment."

        numNodes = size(v.A[1], 1) ÷ 3
        @assert size(s.A[1], 1) == numNodes "Node count mismatch in VectorField assignment."

        @inbounds for e in 1:length(v.A)
            # overwrite only the block belonging to component k
            rows = k:3:3*numNodes
            v.A[e][rows, :] .= s.A[e]
        end
        return v
    end

    # nodal case
    if v.a != [;;]
        @assert s.a != [;;] "Cannot assign elementwise field into nodal VectorField."

        numNodes = size(v.a, 1) ÷ 3
        @assert size(s.a, 1) == numNodes "Node count mismatch in VectorField assignment."

        rows = k:3:3*numNodes
        v.a[rows, :] .= s.a
        return v
    end

    error("VectorField: no data to assign into (empty A and a).")
end

"""
    Base.setindex!(T::TensorField, s::ScalarField, i::Int, j::Int)

    T[i, j] = s
    T[i, :] = v
    T[:, j] = v

Assign a scalar or vector component into a `TensorField`.

## Supported forms

### 1. Assign scalar component
`T[i, j] = s`

- Inserts the scalar field `s` into the tensor component `(i, j)`  
  (with i,j ∈ 1:3).
- Only the rows of the `(i,j)` block are overwritten.

### 2. Assign row vector
`T[i, :] = v`

- `v` must be a `VectorField` with 3 components.
- Replaces the entire i-th tensor row with the components of `v`.

### 3. Assign column vector
`T[:, j] = v`

- `v` must be a `VectorField`.
- Replaces the entire j-th tensor column with the components of `v`.

## Requirements
- Assigned field(s) must belong to the **same problem** as `T`.
- The number of time steps must match.
- Storage mode must match:
    - Elementwise tensor → only elementwise scalar/vector may be assigned.
    - Nodal tensor → only nodal scalar/vector may be assigned.
- Element counts or nodal sizes must be identical.

## Examples
```julia
# scalar insertion
T[1,1] = s1
T[2,3] = s2

# row insertion
T[1,:] = v1    # v1 is VectorField([s11, s12, s13])

# column insertion
T[:,3] = v2
```

## Errors

Errors are thrown when:

* indices are out of bounds
* shapes do not match (nodal vs elementwise mismatch)
* component sizes differ
* attempting unsupported forms like `T[1] = ...` or `T[:, :] = ...`
"""
function Base.setindex!(T::TensorField, s::ScalarField, i::Int, j::Int)
    @assert 1 ≤ i ≤ 3 && 1 ≤ j ≤ 3 "TensorField indexing must be T[i,j] with i,j ∈ 1:3."

    @assert s.model === T.model "TensorField and ScalarField belong to different problems."
    @assert s.nsteps == T.nsteps "Time step mismatch in TensorField assignment."

    block = (j-1)*3 + i  # component index

    # elementwise
    if T.A != []
        @assert s.A != [] "Cannot assign nodal field into elementwise TensorField."
        @assert length(s.A) == length(T.A) "numElem mismatch in TensorField assignment."

        numNodes = size(T.A[1], 1) ÷ 9
        @assert size(s.A[1], 1) == numNodes "Node count mismatch in TensorField assignment."

        rows = block : 9 : 9*numNodes

        @inbounds for e in 1:length(T.A)
            T.A[e][rows, :] .= s.A[e]
        end
        return T
    end

    # nodal
    if T.a != [;;]
        @assert s.a != [;;] "Cannot assign elementwise field into nodal TensorField."

        numNodes = size(T.a, 1) ÷ 9
        @assert size(s.a, 1) == numNodes

        rows = block : 9 : 9*numNodes
        T.a[rows, :] .= s.a
        return T
    end

    error("TensorField: empty storage (A and a both empty).")
end

function Base.setindex!(T::TensorField, v::VectorField, ::Colon, j::Int)
    @assert 1 ≤ j ≤ 3 "TensorField column index must be 1..3."
    @assert v.model === T.model "Model mismatch in TensorField assignment."
    @assert v.nsteps == T.nsteps "Time step mismatch."

    # expand column assignment
    for i in 1:3
        T[i, j] = v[i]   # reuses the ScalarField assignment above
    end

    return T
end

function Base.setindex!(T::TensorField, v::VectorField, i::Int, ::Colon)
    @assert 1 ≤ i ≤ 3 "TensorField row index must be 1..3."
    @assert v.model === T.model "Model mismatch in TensorField assignment."
    @assert v.nsteps == T.nsteps "Time step mismatch."

    # expand row assignment
    for j in 1:3
        T[i, j] = v[j]   # reuses ScalarField setter
    end

    return T
end

"""
    Base.setindex!(T::TensorField, s::ScalarField, k::Int)

T[k] = s

Assign the k-th tensor component in column-major order:

    [T11, T21, T31, T12, T22, T32, T13, T23, T33]
"""
function Base.setindex!(T::TensorField, s::ScalarField, k::Int)
    @assert 1 ≤ k ≤ 9 "TensorField index must be 1..9."

    i = mod1(k, 3)
    j = (k-1) ÷ 3 + 1

    T[i,j] = s   # reuse existing setter

    return T
end

"""
    mergeFields(uu::Vector{<:AbstractField})

Concatenates multiple Fields along a common pseudo-time axis.

Interpretation of `nsteps`:
- If `nsteps == 1`, the Field is treated as an iteration state
  (pseudo-time index: 1, 2, 3, ...).
- If `nsteps > 1`, the Field is treated as a time history, and its
  internal time vector `t` is preserved.

When concatenating multiple time histories, continuity of time is
enforced by shifting each subsequent segment such that it starts
after the previous one. The time increment between segments is
estimated from the last two time points of the previous Field.
If this information is unavailable (`nsteps == 1`), a unit increment
(`Δt = 1`) is assumed.

This function is intended for post-processing and visualization
(e.g. merging iteration histories or successive load/time segments),
not for time-accurate integration.

All input Fields must:
- have `kind=:field`; eigenmode fields cannot be merged on a time axis,
- have the same number of degrees of freedom,
- belong to the same model,
- share the same Field type (ScalarField, VectorField, or TensorField).

Returns a Field of the same concrete type as the first element of `uu`.

# Exaples

u1 = solveDisplacement(problem, load=[load1])   # nsteps = 1
u2 = solveDisplacement(problem, load=[load2])   # nsteps = 1
u3 = solveDisplacement(problem, load=[load3])   # nsteps = 1

u_steps = mergeFields([u1, u2, u3])

# Warning

The time increment between concatenated segments is estimated heuristically.

When a Field has `nsteps == 1`, it is interpreted as an iteration state
and a unit pseudo-time increment (`Δt = 1`) is assumed when concatenating
with subsequent segments.

This heuristic ensures continuity of the time axis but does NOT represent
physical time. The resulting time vector should therefore be used only for
post-processing and visualization, not for time-accurate analysis.
"""
function mergeFields(uu::Vector{<:LowLevelFEM.AbstractField})
    nn = length(uu)
    nn ≥ 1 || error("mergeFields: empty input vector.")

    all(u -> u.kind === :field, uu) ||
        throw(ArgumentError(
            "mergeFields: eigenmode fields cannot be merged as time histories " *
            "because their t values are eigenvalues."
        ))

    u0 = uu[1]
    T = typeof(u0)

    nodal = isNodal(u0)
    elemwise = isElementwise(u0)

    nodal || elemwise || error("mergeFields: invalid field storage.")

    for u in uu
        typeof(u) == T ||
            error("mergeFields: all fields must have the same concrete type.")
        u.model === u0.model ||
            error("mergeFields: all fields must belong to the same model.")
        u.type == u0.type ||
            error("mergeFields: all fields must have the same field type.")
        isNodal(u) == nodal ||
            error("mergeFields: cannot mix nodal and elementwise fields.")
    end

    n = sum(u.nsteps for u in uu)
    t = zeros(n)

    function shifted_time!(t, c, ui, prev)
        ns = ui.nsteps

        if c == 0
            if ns == 1
                t[1] = 1.0
            else
                t[1:ns] .= ui.t
            end
            return
        end

        dt = prev.nsteps > 1 ? prev.t[end] - prev.t[end-1] : 1.0

        if ns == 1
            t[c+1] = t[c] + dt
        else
            t[c+1:c+ns] .= (ui.t .- ui.t[1]) .+ (t[c] + dt)
        end
    end

    if nodal
        ndofs = size(u0.a, 1)
        a = zeros(ndofs, n)

        c = 0
        for i in 1:nn
            ui = uu[i]
            size(ui.a, 1) == ndofs ||
                error("mergeFields: nodal dof mismatch.")

            ns = ui.nsteps
            a[:, c+1:c+ns] .= ui.a
            shifted_time!(t, c, ui, i == 1 ? ui : uu[i-1])
            c += ns
        end

        return T([], a, t, [], n, u0.type, u0.model)

    else
        numElem = copy(u0.numElem)
        ne = length(numElem)
        A = Vector{Matrix{Float64}}(undef, ne)

        index0 = Dict(e => k for (k, e) in enumerate(numElem))

        for u in uu
            Set(u.numElem) == Set(numElem) ||
                error("mergeFields: element sets must be identical.")
        end

        for (ie, elem) in enumerate(numElem)
            rows = size(u0.A[ie], 1)
            A[ie] = zeros(rows, n)
        end

        c = 0
        for i in 1:nn
            ui = uu[i]
            idx = Dict(e => k for (k, e) in enumerate(ui.numElem))

            ns = ui.nsteps
            for elem in numElem
                ie0 = index0[elem]
                iei = idx[elem]

                size(ui.A[iei], 1) == size(A[ie0], 1) ||
                    error("mergeFields: element $elem has incompatible local size.")

                A[ie0][:, c+1:c+ns] .= ui.A[iei]
            end

            shifted_time!(t, c, ui, i == 1 ? ui : uu[i-1])
            c += ns
        end

        return T(A, [;;], t, numElem, n, u0.type, u0.model)
    end
end
#=
function mergeFields(uu::Vector{<:AbstractField})
    nn = length(uu)
    @assert nn ≥ 1

    isElementwise(uu[1]) && error("mergeFields: vector of fields must contain nodal fields.")

    # --- összes időlépés
    n = 0
    for u in uu
        n += u.nsteps
    end

    ndofs = size(uu[1].a, 1)
    a = zeros(ndofs, n)
    t = zeros(n)

    # --- első blokk
    ns = uu[1].nsteps
    a[:, 1:ns] .= uu[1].a

    if ns == 1
        t[1] = 1.0                  # iterációs lépés
    else
        t[1:ns] .= uu[1].t          # valódi idő
    end
    c = ns

    # --- további blokkok
    for i in 2:nn
        ui = uu[i]
        ns = ui.nsteps

        a[:, c+1:c+ns] .= ui.a

        # időlépés becslés az előző szegmensből
        if uu[i-1].nsteps > 1
            dt = uu[i-1].t[end] - uu[i-1].t[end-1]
        else
            dt = 1.0
        end

        t[c+1:c+ns] .= ui.t .+ (t[c] + dt)
        c += ns
    end

    # --- visszaépítés: ugyanaz a Field típus, mint az első
    T = typeof(uu[1])
    return T([], a, t, [], n, uu[1].type, uu[1].model)
end
=#

