###############################################################################
#                                                                             #
#                      Scalar fields operations                               #
#                                                                             #
###############################################################################

import Base.*
import Base./
import Base.+
import Base.-
import Base.log
import Base.sqrt
import Base.cbrt
import Base.abs
export mapScalarField
using StaticArrays

"""
    *(A::ScalarField, B::ScalarField)

Performs element-wise multiplication of two `ScalarField` objects on the same set of elements.

Returns: `ScalarField`

# Examples
```julia
C = A * B
```
"""
function *(A::ScalarField, B::ScalarField)
    A = isNodal(A) ? nodesToElements(A) : A
    B = isNodal(B) ? nodesToElements(B) : B

    A.nsteps == B.nsteps || error("*(ScalarField, ScalarField): nsteps mismatch")
    nsteps = A.nsteps

    # elem → index gyorsleképezés
    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    # közös elemek
    sec = intersect(A.numElem, B.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    @inbounds @views for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        # kisebbik csomópontszám a biztonság kedvéért
        n = min(size(AAi, 1), size(BBi, 1))
        @views C[ii] = AAi[1:n, :] .* BBi[1:n, :]
        num[ii] = e
    end

    return ScalarField(C, [;;], A.t, num, nsteps, :scalar, A.model)
end
#=
function *(AA::ScalarField, BB::ScalarField)
    if isNodal(AA)
        A = nodesToElements(AA)
    else
        A = AA
    end
    if isNodal(BB)
        B = nodesToElements(BB)
    else
        B = BB
    end
    sz = 0
    if A.nsteps != B.nsteps
        error("*(ScalarField, ScalarField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
    end
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    @inbounds for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    @inbounds for i in eachindex(sec)
        n = length(B.A[i])
        D = zeros(n, nsteps)
        if n != sz
            #D = zeros(n, nsteps)
            sz = n
        end
        @inbounds for j in 1:n
            @inbounds for k in 1:nsteps
                D[j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][j, k]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, num, A.nsteps, :scalar, A.model)
end
=#

"""
    /(A::ScalarField, B::ScalarField)

Performs element-wise division of two `ScalarField` objects on the same set of elements.

Returns: `ScalarField`

# Examples
```julia
C = A / B
```
"""
function /(AA::ScalarField, BB::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    A.nsteps == B.nsteps || error("/(ScalarField, ScalarField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps))")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(A.numElem, B.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        n = min(size(AAi, 1), size(BBi, 1))
        @views C[ii] = AAi[1:n, :] ./ BBi[1:n, :]
        num[ii] = e
    end

    return ScalarField(C, [;;], A.t, num, A.nsteps, :scalar, A.model)
end

"""
    +(A::ScalarField, B::ScalarField)

Performs element-wise addition of two `ScalarField` objects on the same set of elements.

Returns: `ScalarField`

# Examples
```julia
C = A + B
```
"""
function +(A::ScalarField, B::ScalarField)
    if isElementwise(A) && isElementwise(B)
        if A.type == B.type
            A.nsteps == B.nsteps || error("+(A::ScalarField, B::ScalarField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")

            a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
            b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

            sec = intersect(A.numElem, B.numElem)
            dif1 = setdiff(A.numElem, B.numElem)
            dif2 = setdiff(B.numElem, A.numElem)

            total = length(sec) + length(dif1) + length(dif2)
            C = Vector{Matrix{Float64}}(undef, total)
            num = Vector{Int}(undef, total)

            pos = 1

            @inbounds for e in sec
                ia = a_index[e]
                ib = b_index[e]
                C[pos] = A.A[ia] .+ B.A[ib]
                num[pos] = e
                pos += 1
            end

            @inbounds for e in dif1
                ia = a_index[e]
                C[pos] = A.A[ia]
                num[pos] = e
                pos += 1
            end

            @inbounds for e in dif2
                ib = b_index[e]
                C[pos] = B.A[ib]
                num[pos] = e
                pos += 1
            end

            return ScalarField(C, [;;], A.t, num, A.nsteps, A.type, A.model)
        else
            error("+(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if A.type == B.type
            return ScalarField([], A.a + B.a, A.t, [], A.nsteps, A.type, A.model)
        else
            error("+(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("+(ScalarField, ScalarField): internal error")
    end
end

"""
    -(A::ScalarField, B::ScalarField)

Performs element-wise subtraction of two `ScalarField` objects on the same set of elements.

Returns: `ScalarField`

# Examples
```julia
C = A - B
```
"""
function -(A::ScalarField, B::ScalarField)
    if isElementwise(A) && isElementwise(B)
        if A.type == B.type
            A.nsteps == B.nsteps || error("-(A::ScalarField, B::ScalarField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")

            a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
            b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

            sec = intersect(A.numElem, B.numElem)
            dif1 = setdiff(A.numElem, B.numElem)
            dif2 = setdiff(B.numElem, A.numElem)

            total = length(sec) + length(dif1) + length(dif2)
            C = Vector{Matrix{Float64}}(undef, total)
            num = Vector{Int}(undef, total)

            pos = 1

            @inbounds for e in sec
                ia = a_index[e]
                ib = b_index[e]
                C[pos] = A.A[ia] .- B.A[ib]
                num[pos] = e
                pos += 1
            end

            @inbounds for e in dif1
                ia = a_index[e]
                C[pos] = A.A[ia]
                num[pos] = e
                pos += 1
            end

            @inbounds for e in dif2
                ib = b_index[e]
                C[pos] = -B.A[ib]
                num[pos] = e
                pos += 1
            end

            return ScalarField(C, [;;], A.t, num, A.nsteps, A.type, A.model)
        else
            error("-(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if A.type == B.type
            return ScalarField([], A.a - B.a, A.t, [], A.nsteps, A.type, A.model)
        else
            error("-(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("-(ScalarField, ScalarField): internal error")
    end
end

function -(A::ScalarField)
    return A * (-1)
end

"""
    *(A::ScalarField, b::Number)

Performs multiplication of a ScalarField objects and a Number.

Return: ScalarField

# Examples
```julia
C = A * 2.0
```
"""
function *(A::ScalarField, b::Number)
    if isElementwise(A)
        n = length(A.A)
        C = Vector{Matrix{Float64}}(undef, n)
        @inbounds for i in 1:n
            C[i] = A.A[i] .* b
        end
        return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
    else
        return ScalarField(A.A, A.a * b, A.t, A.numElem, A.nsteps, A.type, A.model)
    end
end

"""
    *(b::Number, A::ScalarField)

Performs multiplication of a ScalarField objects and a Number.

Return: ScalarField

# Examples
```julia
C = 2.0 * A
```
"""
function *(b::Number, A::ScalarField)
    return A * b
end

"""
    /(A::ScalarField, b::Number)

Elementwise division of a scalar field by a constant.

Each element-wise matrix of the scalar field is divided by the scalar `b`.
If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the elementwise divided values.
"""
function /(AA::ScalarField, b::Number)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = A.A[i] ./ b
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    /(b::Number, A::ScalarField)

Elementwise division of a constant by a scalar field.

Each element-wise matrix of the scalar field is used as the divisor of the
constant `b`, i.e. `b ./ A`.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the elementwise divided values.
"""
function /(b::Number, AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = b ./ A.A[i]
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    log(A::ScalarField)

Elementwise natural logarithm of a scalar field.

Applies the natural logarithm to each entry of every element-wise matrix
of the scalar field.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the elementwise logarithmic values.
"""
function log(AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = log.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    sqrt(A::ScalarField)

Elementwise square root of a scalar field.

Applies the square root to each entry of every element-wise matrix
of the scalar field.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the elementwise square-rooted values.
"""
function sqrt(AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = sqrt.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    cbrt(A::ScalarField)

Elementwise cubic root of a scalar field.

Applies the cubic root to each entry of every element-wise matrix
of the scalar field.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the elementwise cubic-rooted values.
"""
function cbrt(AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = cbrt.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    abs(A::ScalarField)

Elementwise absolute value of a scalar field.

Applies the absolute value to each entry of every element-wise matrix
of the scalar field.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the elementwise absolute values.
"""
function abs(AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = abs.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    +(A::ScalarField, b::Number)

Add a constant offset to a scalar field.

The scalar `b` is added elementwise to each entry of every element-wise
matrix of the scalar field.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the shifted values.
"""
function +(AA::ScalarField, b::Number)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = A.A[i] .+ b
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    -(A::ScalarField, b::Number)

Subtract a constant offset from a scalar field.

The scalar `b` is subtracted elementwise from each entry of every element-wise
matrix of the scalar field.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the shifted values.
"""
function -(AA::ScalarField, b::Number)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = A.A[i] .- b
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    +(b::Number, A::ScalarField)

Add a constant offset to a scalar field.

The scalar `b` is added elementwise to each entry of every element-wise
matrix of the scalar field.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the shifted values.
"""
function +(b::Number, AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = A.A[i] .+ b
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    -(b::Number, A::ScalarField)

Subtract a constant offset from a scalar field.

The scalar `b` is subtracted elementwise from each entry of every element-wise
matrix of the scalar field.

If the field is nodal, it is first converted to elementwise form.

# Returns
- A new `ScalarField` containing the shifted values.
"""
function -(b::Number, AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = A.A[i] .- b
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

"""
    mapScalarField(f, A::ScalarField)

Apply a function elementwise to a scalar field.

The function `f` is applied to each element-wise matrix of the scalar field.
If the field is nodal, it is first converted to elementwise form.

This is a low-level helper used to implement elementwise scalar-field
operations such as `abs`, `+`, `-`, `log`, `sqrt`, etc.

# Returns
- A new `ScalarField` containing the transformed values.
"""
@inline function mapScalarField(f::F, AA::ScalarField) where {F}
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = f.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

#=
function /(AA::ScalarField, b::Number)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = A.A[i] ./ b
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

function /(b::Number, AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = b ./ A.A[i]
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

function log(AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = log.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

function sqrt(AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = sqrt.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

function cbrt(AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = cbrt.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

function abs(AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = abs.(A.A[i])
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

function +(AA::ScalarField, b::Number)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = A.A[i] .+ b
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

function -(AA::ScalarField, b::Number)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = A.A[i] .- b
    end
    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
end

@inline function mapScalarField(f::F, AA::ScalarField) where {F}
    A = isNodal(AA) ? nodesToElements(AA) : AA
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    @inbounds for i in 1:n
        C[i] = f(A.A[i])
    end
    return ScalarField(C, [:,], A.t, A.numElem, A.nsteps, A.type, A.model)
end
=#

###############################################################################
#                                                                             #
#                      Vector fields operations                               #
#                                                                             #
###############################################################################

#import LinearAlgebra: dot
import LinearAlgebra.dot
import LinearAlgebra.⋅
import LinearAlgebra.×
import Base.∘
import LinearAlgebra.norm
import LinearAlgebra.diagm

"""
    *(A::ScalarField, B::VectorField)

Scales a `VectorField` by a `ScalarField` element-wise on matching elements.

Returns: `VectorField`

# Examples
```julia
v2 = s .* v  # equivalent to s * v
```
"""
function *(AA::ScalarField, BB::VectorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    B.type == :v3D || error("*(AA::ScalarField, BB::VectorField): Vector field must be 3 dimansional.")
    if A.nsteps != B.nsteps || A.nsteps != 1
        error("*(AA::ScalarField, BB::VectorField): Number of nsteps in AA must be one or equal to BB.nsteps.")
    end

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(B.numElem, A.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = B.nsteps
    a_has_single_step = (A.nsteps == 1)

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        n_nodes = min(size(AAi, 1), size(BBi, 1) ÷ 3)
        rows = 3n_nodes

        result = Matrix{Float64}(undef, rows, nsteps)

        if n_nodes > 0
            @inbounds @views for j in 1:n_nodes
                dest = result[3j-2:3j, :]
                src = BBi[3j-2:3j, :]
                if a_has_single_step
                    dest .= src .* AAi[j, 1]
                else
                    dest .= src .* AAi[j, :]
                end
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return VectorField(C, [;;], B.t, num, B.nsteps, B.type, B.model)
end

"""
    *(B::VectorField, A::ScalarField)

Scales a `VectorField` by a `ScalarField` element-wise on matching elements.

Returns: `VectorField`
"""
function *(BB::VectorField, AA::ScalarField)
    return AA * BB
end

"""
    /(B::VectorField, A::ScalarField)

Divides a `VectorField` by a `ScalarField` element-wise on matching elements.

Returns: `VectorField`
"""
function /(BB::VectorField, AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    B.type == :v3D || error("/(BB::VectorField, AA::ScalarField): Vector field must be 3 dimansional.")
    if A.nsteps != B.nsteps || A.nsteps != 1
        error("/(BB::VectorField, AA::ScalarField): Number of nsteps in AA must be one or equal to BB.nsteps.")
    end

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(B.numElem, A.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = B.nsteps
    a_has_single_step = (A.nsteps == 1)

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        n_nodes = min(size(AAi, 1), size(BBi, 1) ÷ 3)
        rows = 3n_nodes

        result = Matrix{Float64}(undef, rows, nsteps)

        if n_nodes > 0
            @inbounds @views for j in 1:n_nodes
                dest = result[3j-2:3j, :]
                src = BBi[3j-2:3j, :]
                if a_has_single_step
                    dest .= src ./ AAi[j, 1]
                else
                    dest .= src ./ AAi[j, :]
                end
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return VectorField(C, [;;], B.t, num, B.nsteps, BB.type, B.model)
end

function +(A::VectorField, B::VectorField)
    if isElementwise(A) && isElementwise(B)
        if A.type == B.type
            #if length(A.A) != length(B.A)
            #    error("+(A::VectorField, B::VectorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            #end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("+(A::VectorField, B::VectorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
            b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

            sec = intersect(A.numElem, B.numElem)
            dif1 = setdiff(A.numElem, B.numElem)
            dif2 = setdiff(B.numElem, A.numElem)

            total = length(sec) + length(dif1) + length(dif2)
            C = Vector{Matrix{Float64}}(undef, total)
            num = Vector{Int}(undef, total)

            pos = 1

            @inbounds for e in sec
                ia = a_index[e]
                ib = b_index[e]
                C[pos] = A.A[ia] .+ B.A[ib]
                num[pos] = e
                pos += 1
            end

            @inbounds for e in dif1
                ia = a_index[e]
                C[pos] = A.A[ia]
                num[pos] = e
                pos += 1
            end

            @inbounds for e in dif2
                ib = b_index[e]
                C[pos] = B.A[ib]
                num[pos] = e
                pos += 1
            end

            return VectorField(C, [;;], A.t, num, A.nsteps, A.type, A.model)
        else
            error("+(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if A.type == B.type
            return VectorField([], A.a + B.a, A.t, [], A.nsteps, A.type, A.model)
        else
            error("+(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("+(VectorField, VectorField): internal error")
    end
end

function -(A::VectorField, B::VectorField)
    if isElementwise(A) && isElementwise(B)
        if A.type == B.type
            #if length(A.A) != length(B.A)
            #    error("+(A::VectorField, B::VectorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            #end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("-(A::VectorField, B::VectorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
            b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

            sec = intersect(A.numElem, B.numElem)
            dif1 = setdiff(A.numElem, B.numElem)
            dif2 = setdiff(B.numElem, A.numElem)

            total = length(sec) + length(dif1) + length(dif2)
            C = Vector{Matrix{Float64}}(undef, total)
            num = Vector{Int}(undef, total)

            pos = 1

            @inbounds for e in sec
                ia = a_index[e]
                ib = b_index[e]
                C[pos] = A.A[ia] .- B.A[ib]
                num[pos] = e
                pos += 1
            end

            @inbounds for e in dif1
                ia = a_index[e]
                C[pos] = A.A[ia]
                num[pos] = e
                pos += 1
            end

            @inbounds for e in dif2
                ib = b_index[e]
                C[pos] = -B.A[ib]
                num[pos] = e
                pos += 1
            end

            return VectorField(C, [;;], A.t, num, A.nsteps, A.type, A.model)
        else
            error("-(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if A.type == B.type
            return VectorField([], A.a - B.a, A.t, [], A.nsteps, A.type, A.model)
        else
            error("-(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("-(VectorField, VectorField): internal error")
    end
end

function -(A::VectorField)
    return A * (-1)
end

function *(A::VectorField, b::Number)
    if isElementwise(A)
        n = length(A.A)
        C = Vector{Matrix{Float64}}(undef, n)
        @inbounds for i in 1:n
            C[i] = A.A[i] .* b
        end
        return VectorField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
    elseif A.a != [;;]
        return VectorField([], A.a .* b, A.t, [], A.nsteps, A.type, A.model)
    else
        error("*(VectorField, b): internal error")
    end
end

function *(b::Number, A::VectorField)
    return A * b
end

function /(A::VectorField, b::Number)
    if isElementwise(A)
        n = length(A.A)
        C = Vector{Matrix{Float64}}(undef, n)
        @inbounds for i in 1:n
            C[i] = A.A[i] ./ b
        end
        return VectorField(C, [;;], A.t, A.numElem, A.nsteps, A.type, A.model)
    elseif A.a != [;;]
        return VectorField([], A.a ./ b, A.t, [], A.nsteps, A.type, A.model)
    else
        error("/(VectorField, b): internal error")
    end
end

function dot(AA::VectorField, BB::VectorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    if A.nsteps != B.nsteps
        error("*(VectorField, VectorField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
    end
    if A.type != :v3D && B.type != :v3D
        error("*(AA::VectorField, BB::VectorField): AA and BB must be 3D vectors.")
    end

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(B.numElem, A.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = B.nsteps

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        rows = min(size(AAi, 1), size(BBi, 1))
        n_nodes = fld(rows, 3)
        used_rows = 3n_nodes

        result = Matrix{Float64}(undef, n_nodes, nsteps)

        if n_nodes > 0
            Ab = view(AAi, 1:used_rows, 1:nsteps)
            Bb = view(BBi, 1:used_rows, 1:nsteps)

            @inbounds @views for j in 1:n_nodes
                a_block = Ab[3j-2:3j, :]
                b_block = Bb[3j-2:3j, :]
                row = result[j, :]
                row .= a_block[1, :] .* b_block[1, :] .+
                        a_block[2, :] .* b_block[2, :] .+
                        a_block[3, :] .* b_block[3, :]
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return ScalarField(C, [;;], A.t, num, A.nsteps, :scalar, A.model)
end

#=
function ⋅(AA::VectorField, BB::VectorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if BB.A == []
        B = nodesToElements(BB)
    else
        B = BB
    end
    if true #(A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
        if length(A.A) != length(B.A)
            error("*(A::VectorField, B::VectorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
        end
        if A.numElem != B.numElem
            error("*(A::VectorField, B::VectorField): vector fields are not compatible.")
        end
        nsteps = A.nsteps
        nsteps2 = B.nsteps
        if nsteps != nsteps2
            error("*(A::VectorField, B::VectorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
        end
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 3
            m = length(B.A[i]) ÷ 3
            if n != m
                error("*(A::VectorField, B::VectorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
            end
            D = zeros(n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[j, k] = reshape(A.A[i][3j-2:3j, k], 3, 1) ⋅ reshape(B.A[i][3j-2:3j, k], 3, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
    else
        error("*(A::VectorField, B::VectorField): VectorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end
=#

function *(A::VectorField, B::VectorField)
    return dot(A, B)
end

function ∘(AA::VectorField, BB::VectorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    if A.nsteps != B.nsteps
        error("∘(VectorField, VectorField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
    end
    if A.type != :v3D && B.type != :v3D
        error("∘(AA::VectorField, BB::VectorField): AA and BB must be 3D vectors.")
    end

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(B.numElem, A.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = B.nsteps

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        rows = min(size(AAi, 1), size(BBi, 1))
        n_nodes = fld(rows, 3)
        used_rows = 3n_nodes

        result = Matrix{Float64}(undef, 9n_nodes, nsteps)

        if n_nodes > 0
            Ab = view(AAi, 1:used_rows, 1:nsteps)
            Bb = view(BBi, 1:used_rows, 1:nsteps)

            @inbounds @views for j in 1:n_nodes
                dest = result[9j-8:9j, :]
                a_block = Ab[3j-2:3j, :]
                b_block = Bb[3j-2:3j, :]

                dest[1, :] .= a_block[1, :] .* b_block[1, :]
                dest[2, :] .= a_block[2, :] .* b_block[1, :]
                dest[3, :] .= a_block[3, :] .* b_block[1, :]

                dest[4, :] .= a_block[1, :] .* b_block[2, :]
                dest[5, :] .= a_block[2, :] .* b_block[2, :]
                dest[6, :] .= a_block[3, :] .* b_block[2, :]

                dest[7, :] .= a_block[1, :] .* b_block[3, :]
                dest[8, :] .= a_block[2, :] .* b_block[3, :]
                dest[9, :] .= a_block[3, :] .* b_block[3, :]
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return TensorField(C, [;;], A.t, num, A.nsteps, :e, A.model)
end

#=
function ∘(AA::VectorField, BB::VectorField)
    if AA.A == []
        a = nodesToElements(AA)
    else
        a = AA
    end
    if BB.A == []
        b = nodesToElements(BB)
    else
        b = BB
    end
    G = MMatrix{3,3}([0.0 0 0; 0 0 0; 0 0 0])
    if length(a.A) != 0 && length(b.A) != 0
        if a isa VectorField && b isa VectorField
            nsteps = a.nsteps
            C = []
            for i in 1:length(a.A)
                n = length(a.A[i]) ÷ 3
                H = zeros(9n, nsteps)
                for j in 1:n
                    for k in 1:nsteps
                        e = reshape(a.A[i][3j-2:3j, k], 3, 1)
                        f = reshape(b.A[i][3j-2:3j, k], 3, 1)
                        for p in 1:3, q in 1:3
                            G[p, q] = e[p, 1] * f[q, 1]
                        end
                        H[9j-8:9j, k] = G[:]
                    end
                end
                push!(C, H)
            end
            aa = [;;]
            return TensorField(C, aa, a.t, a.numElem, a.nsteps, :e, a.model)
        else
            error("∘(a::VectorField, b::VectorField): a and b are not VectorField(s).")
        end
    else
        error("∘(a::VectorField, b::VectorField): data at nodes is not yet implemented.")
    end
end
=#

"""
    ×(a::VectorField, b::VectorField)

Element-wise 3D vector cross product on matching elements.

Returns: `VectorField`

# Examples
```julia
w = u × v
```
"""
function ×(AA::VectorField, BB::VectorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    if A.nsteps != B.nsteps
        error("×(VectorField, VectorField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
    end
    if A.type != :v3D && B.type != :v3D
        error("×(AA::VectorField, BB::VectorField): AA and BB must be 3D vectors.")
    end

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(B.numElem, A.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = B.nsteps

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        rows = min(size(AAi, 1), size(BBi, 1))
        n_nodes = fld(rows, 3)
        used_rows = 3n_nodes

        result = Matrix{Float64}(undef, used_rows, nsteps)

        if n_nodes > 0
            Ab = view(AAi, 1:used_rows, 1:nsteps)
            Bb = view(BBi, 1:used_rows, 1:nsteps)

            @inbounds @views for j in 1:n_nodes
                dest = result[3j-2:3j, :]
                a_block = Ab[3j-2:3j, :]
                b_block = Bb[3j-2:3j, :]
                dest[1, :] .= a_block[2, :] .* b_block[3, :] .- a_block[3, :] .* b_block[2, :]
                dest[2, :] .= a_block[3, :] .* b_block[1, :] .- a_block[1, :] .* b_block[3, :]
                dest[3, :] .= a_block[1, :] .* b_block[2, :] .- a_block[2, :] .* b_block[1, :]
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return VectorField(C, [;;], A.t, num, A.nsteps, :v3D, A.model)
end

"""
    norm(A::VectorField)

Element-wise Euclidean norm of a `VectorField`.

Returns: `ScalarField`
"""
function norm(AA::VectorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    nsteps = A.nsteps

    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds for i in 1:n
        Ai = A.A[i]
        rows = size(Ai, 1)
        n_nodes = fld(rows, 3)
        used_rows = 3n_nodes

        result = Matrix{Float64}(undef, n_nodes, nsteps)

        if n_nodes > 0
            Ab = view(Ai, 1:used_rows, 1:nsteps)
            @inbounds @views for j in 1:n_nodes
                block = Ab[3j-2:3j, :]
                result[j, :] .= sqrt.(block[1, :].^2 .+ block[2, :].^2 .+ block[3, :].^2)
            end
        end

        C[i] = result
    end

    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, :scalar, A.model)
end

"""
    diagm(A::VectorField)

Creates a diagonal `TensorField` from a `VectorField` (dim=3), i.e., places vector
components on the tensor diagonal for each node/element.

Returns: `TensorField`
"""
function diagm(AA::VectorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    if A isa VectorField && A.type == :v3D
        nsteps = A.nsteps
        n = length(A.A)
        C = Vector{Matrix{Float64}}(undef, n)

        @inbounds for i in 1:n
            Ai = A.A[i]
            rows = size(Ai, 1)
            n_nodes = fld(rows, 3)
            used_rows = 3n_nodes

            result = Matrix{Float64}(undef, 9n_nodes, nsteps)

            if n_nodes > 0
                Ab = view(Ai, 1:used_rows, 1:nsteps)
                @inbounds @views for j in 1:n_nodes
                    dest = result[9j-8:9j, :]
                    block = Ab[3j-2:3j, :]
                    dest .= 0.0
                    dest[1, :] .= block[1, :]
                    dest[5, :] .= block[2, :]
                    dest[9, :] .= block[3, :]
                end
            end

            C[i] = result
        end

        return TensorField(C, [;;], A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("diagm(A::VectorField): A is not a VectorField or dim != 3.")
    end
end

export ⋅

###############################################################################
#                                                                             #
#                      Tensor fields operations                               #
#                                                                             #
###############################################################################

import Base.inv
import Base.transpose
import Base.adjoint
import LinearAlgebra.det
import LinearAlgebra.eigen

export det
export unitTensor
export trace

"""
    *(A::TensorField, B::TensorField)

Tensor contraction (matrix multiplication) for each element/node: reshapes 9×1 blocks
into 3×3, multiplies, then flattens back.

Returns: `TensorField`
"""
function *(AA::TensorField, BB::TensorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    A.nsteps == B.nsteps || error("*(A::TensorField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(A.numElem, B.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = A.nsteps

    @inbounds @views for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        rows = min(size(AAi, 1), size(BBi, 1))
        nblocks = fld(rows, 9)
        used_rows = 9nblocks

        result = Matrix{Float64}(undef, used_rows, nsteps)

        if nblocks > 0
            @inbounds for j in 1:nblocks
                Ab = AAi[9j-8:9j, :]
                Bb = BBi[9j-8:9j, :]
                for k in 1:nsteps
                    MA = SMatrix{3,3}(reshape(Ab[:, k], 3, 3))
                    MB = SMatrix{3,3}(reshape(Bb[:, k], 3, 3))
                    dest = view(result, 9j-8:9j, k)
                    dest .= (MA * MB)[:]
                end
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return TensorField(C, [;;], A.t, num, A.nsteps, A.type, A.model)
end

"""
    ⋅(A::TensorField, B::TensorField)

Element-wise (Hadamard) product followed by summation of all components, yielding a
scalar per tensor (i.e., Frobenius inner product).

Returns: `ScalarField`
"""
function ⋅(AA::TensorField, BB::TensorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    A.nsteps == B.nsteps || error("*(A::TensorField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(A.numElem, B.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = A.nsteps

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        rows = min(size(AAi, 1), size(BBi, 1))
        nblocks = fld(rows, 9)
        used_rows = 9nblocks

        result = Matrix{Float64}(undef, nblocks, nsteps)

        if nblocks > 0
            Ab = AAi[1:used_rows, 1:nsteps]
            Bb = BBi[1:used_rows, 1:nsteps]
            @inbounds for j in 1:nblocks
                blockA = Ab[9j-8:9j, :]
                blockB = Bb[9j-8:9j, :]
                row = view(result, j, :)
                row .= 0.0
                @inbounds for r in 1:9
                    row .+= blockA[r, :] .* blockB[r, :]
                end
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return ScalarField(C, [;;], A.t, num, A.nsteps, :scalar, A.model)
end

function +(AA::TensorField, BB::TensorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    A.nsteps == B.nsteps || error("+(A::TensorField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(A.numElem, B.numElem)
    dif1 = setdiff(A.numElem, B.numElem)
    dif2 = setdiff(B.numElem, A.numElem)

    total = length(sec) + length(dif1) + length(dif2)
    C = Vector{Matrix{Float64}}(undef, total)
    num = Vector{Int}(undef, total)

    pos = 1

    @inbounds for e in sec
        ia = a_index[e]
        ib = b_index[e]
        C[pos] = A.A[ia] .+ B.A[ib]
        num[pos] = e
        pos += 1
    end

    @inbounds for e in dif1
        ia = a_index[e]
        C[pos] = A.A[ia]
        num[pos] = e
        pos += 1
    end

    @inbounds for e in dif2
        ib = b_index[e]
        C[pos] = B.A[ib]
        num[pos] = e
        pos += 1
    end

    return TensorField(C, [;;], A.t, num, A.nsteps, A.type, A.model)
end

function -(AA::TensorField, BB::TensorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    A.nsteps == B.nsteps || error("-(A::TensorField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(A.numElem, B.numElem)
    dif1 = setdiff(A.numElem, B.numElem)
    dif2 = setdiff(B.numElem, A.numElem)

    total = length(sec) + length(dif1) + length(dif2)
    C = Vector{Matrix{Float64}}(undef, total)
    num = Vector{Int}(undef, total)

    pos = 1

    @inbounds for e in sec
        ia = a_index[e]
        ib = b_index[e]
        C[pos] = A.A[ia] .- B.A[ib]
        num[pos] = e
        pos += 1
    end

    @inbounds for e in dif1
        ia = a_index[e]
        C[pos] = A.A[ia]
        num[pos] = e
        pos += 1
    end

    @inbounds for e in dif2
        ib = b_index[e]
        C[pos] = -B.A[ib]
        num[pos] = e
        pos += 1
    end

    return TensorField(C, [;;], A.t, num, A.nsteps, A.type, A.model)
end

function -(A::TensorField)
    return A * (-1)
end

function *(AA::TensorField, BB::VectorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    b = isNodal(BB) ? nodesToElements(BB) : BB

    b.type == :v3D || error("*(A::TensorField, b::VectorField): Multiply by 2D vector is not yet implemented.")
    A.nsteps == b.nsteps || error("*(A::TensorField, b::VectorField): nsteps of A=$(A.nsteps) != nsteps of b=$(b.nsteps)")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(b.numElem))

    sec = intersect(A.numElem, b.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = A.nsteps

    @inbounds @views for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = b.A[ib]

        rowsA = size(AAi, 1)
        rowsB = size(BBi, 1)
        nblocks = min(fld(rowsA, 9), fld(rowsB, 3))
        usedA = 9nblocks
        usedB = 3nblocks

        result = Matrix{Float64}(undef, usedB, nsteps)

        if nblocks > 0
            Ab = AAi[1:usedA, 1:nsteps]
            Bb = BBi[1:usedB, 1:nsteps]
            @inbounds for j in 1:nblocks
                blockA = Ab[9j-8:9j, :]
                blockB = Bb[3j-2:3j, :]
                dest = view(result, 3j-2:3j, :)
                for k in 1:nsteps
                    MA = SMatrix{3,3}(reshape(blockA[:, k], 3, 3))
                    v = SVector{3}(blockB[:, k])
                    dest[:, k] = MA * v
                end
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return VectorField(C, [;;], A.t, num, A.nsteps, b.type, A.model)
end

function *(BB::VectorField, AA::TensorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    b = isNodal(BB) ? nodesToElements(BB) : BB

    b.type == :v3D || error("*(A::TensorField, b::VectorField): Multiply by 2D vector is not yet implemented.")
    A.nsteps == b.nsteps || error("*(A::TensorField, b::VectorField): nsteps of A=$(A.nsteps) != nsteps of b=$(b.nsteps)")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(b.numElem))

    sec = intersect(A.numElem, b.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = A.nsteps

    @inbounds @views for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = b.A[ib]

        rowsA = size(AAi, 1)
        rowsB = size(BBi, 1)
        nblocks = min(fld(rowsA, 9), fld(rowsB, 3))
        usedA = 9nblocks
        usedB = 3nblocks

        result = Matrix{Float64}(undef, usedB, nsteps)

        if nblocks > 0
            Ab = AAi[1:usedA, 1:nsteps]
            Bb = BBi[1:usedB, 1:nsteps]
            @inbounds for j in 1:nblocks
                blockA = Ab[9j-8:9j, :]
                blockB = Bb[3j-2:3j, :]
                dest = view(result, 3j-2:3j, :)
                for k in 1:nsteps
                    M = SMatrix{3,3}(reshape(blockA[:, k], 3, 3))
                    v = SVector{3}(blockB[:, k])
                    dest[:, k] = M' * v
                end
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return VectorField(C, [;;], A.t, num, A.nsteps, b.type, A.model)
end

function *(AA::TensorField, b::Number)
    if isElementwise(AA)
        n = length(AA.A)
        C = Vector{Matrix{Float64}}(undef, n)
        @inbounds for i in 1:n
            C[i] = AA.A[i] .* b
        end
        return TensorField(C, [;;], AA.t, AA.numElem, AA.nsteps, AA.type, AA.model)
    elseif AA.a != [;;]
        return TensorField(AA.A, AA.a .* b, AA.t, AA.numElem, AA.nsteps, AA.type, AA.model)
    else
        error("*(A::TensorField, b): TensorField type ($(AA.type)) is not yet implemented.")
    end
end

function *(b::Number, AA::TensorField)
    return AA * b
end

function /(AA::TensorField, b::Number)
    if isElementwise(AA)
        n = length(AA.A)
        C = Vector{Matrix{Float64}}(undef, n)
        @inbounds for i in 1:n
            C[i] = AA.A[i] ./ b
        end
        return TensorField(C, [;;], AA.t, AA.numElem, AA.nsteps, AA.type, AA.model)
    elseif AA.a != [;;]
        return TensorField(AA.A, AA.a ./ b, AA.t, AA.numElem, AA.nsteps, AA.type, AA.model)
    else
        error("/(A::TensorField, b): TensorField type ($(AA.type)) is not yet implemented.")
    end
end

function *(AA::ScalarField, BB::TensorField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    (A.nsteps == 1 || A.nsteps == B.nsteps) || error("*(ScalarField, TensorField): nsteps of ScalarField ($(A.nsteps)) must be one or equal to TensorField ($(B.nsteps))")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(B.numElem, A.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = B.nsteps
    a_has_single_step = (A.nsteps == 1)

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        rows = size(BBi, 1)
        nblocks = fld(rows, 9)
        used_rows = 9nblocks

        result = Matrix{Float64}(undef, used_rows, nsteps)

        if nblocks > 0
            Bb = view(BBi, 1:used_rows, 1:nsteps)
            @inbounds for j in 1:nblocks
                dest = view(result, 9j-8:9j, :)
                blockB = view(Bb, 9j-8:9j, :)
                if a_has_single_step
                    dest .= blockB .* AAi[j, 1]
                else
                    dest .= blockB .* AAi[j, :]
                end
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return TensorField(C, [;;], B.t, num, B.nsteps, B.type, B.model)
end

function *(BB::TensorField, AA::ScalarField)
    return AA * BB
end

function /(BB::TensorField, AA::ScalarField)
    A = isNodal(AA) ? nodesToElements(AA) : AA
    B = isNodal(BB) ? nodesToElements(BB) : BB

    (A.nsteps == 1 || A.nsteps == B.nsteps) || error("/(TensorField, ScalarField): nsteps of ScalarField ($(A.nsteps)) must be one or equal to TensorField ($(B.nsteps))")

    a_index = Dict(e => i for (i, e) in enumerate(A.numElem))
    b_index = Dict(e => i for (i, e) in enumerate(B.numElem))

    sec = intersect(B.numElem, A.numElem)
    nsec = length(sec)

    C = Vector{Matrix{Float64}}(undef, nsec)
    num = Vector{Int}(undef, nsec)

    nsteps = B.nsteps
    a_has_single_step = (A.nsteps == 1)

    @inbounds for (ii, e) in enumerate(sec)
        ia = a_index[e]
        ib = b_index[e]

        AAi = A.A[ia]
        BBi = B.A[ib]

        rows = size(BBi, 1)
        nblocks = fld(rows, 9)
        used_rows = 9nblocks

        result = Matrix{Float64}(undef, used_rows, nsteps)

        if nblocks > 0
            Bb = view(BBi, 1:used_rows, 1:nsteps)
            @inbounds for j in 1:nblocks
                dest = view(result, 9j-8:9j, :)
                blockB = view(Bb, 9j-8:9j, :)
                if a_has_single_step
                    dest .= blockB ./ AAi[j, 1]
                else
                    dest .= blockB ./ AAi[j, :]
                end
            end
        end

        C[ii] = result
        num[ii] = e
    end

    return TensorField(C, [;;], B.t, num, B.nsteps, B.type, B.model)
end

"""
    transpose(A::TensorField)

Transposes each 3×3 tensor block.

Returns: `TensorField`
"""
function transpose(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds @views for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, 9nblocks, nsteps)

        @inbounds for j in 1:nblocks
            block = Ai[9j-8:9j, :]
            dest = result[9j-8:9j, :]
            for k in 1:nsteps
                M = reshape(block[:, k], 3, 3)
                dest[:, k] .= transpose(M)[:]
            end
        end

        C[i] = result
    end

    return TensorField(C, [;;], A.t, A.numElem, A.nsteps, :e, A.model)
end

"""
    adjoint(A::TensorField)

Adjoint (conjugate transpose) of each 3×3 tensor block.

Returns: `TensorField`
"""
function adjoint(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds @views for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, 9nblocks, nsteps)

        @inbounds for j in 1:nblocks
            block = Ai[9j-8:9j, :]
            dest = result[9j-8:9j, :]
            for k in 1:nsteps
                M = reshape(block[:, k], 3, 3)
                dest[:, k] .= adjoint(M)[:]
            end
        end

        C[i] = result
    end

    return TensorField(C, [;;], A.t, A.numElem, A.nsteps, :e, A.model)
end

"""
    unitTensor(A::TensorField)

Creates an identity tensor field (I) with the same element structure and time steps as `A`.

Returns: `TensorField`
"""
function unitTensor(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)
    ident = SMatrix{3,3}(I)

    @inbounds for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, 9nblocks, nsteps)

        @inbounds for j in 1:nblocks
            dest = view(result, 9j-8:9j, :)
            for k in 1:nsteps
                dest[:, k] .= ident[:]
            end
        end

        C[i] = result
    end

    return TensorField(C, [;;], A.t, A.numElem, A.nsteps, :e, A.model)
end

"""
    trace(A::TensorField)

Computes the trace of each 3×3 tensor block.

Returns: `ScalarField`
"""
function trace(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds @views for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, nblocks, nsteps)

        @inbounds for j in 1:nblocks
            block = Ai[9j-8:9j, :]
            row = result[j, :]
            for k in 1:nsteps
                row[k] = block[1, k] + block[5, k] + block[9, k]
            end
        end

        C[i] = result
    end

    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, :scalar, A.model)
end

"""
    det(A::TensorField)

Computes the determinant of each 3×3 tensor block.

Returns: `ScalarField`
"""
function det(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds @views for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, nblocks, nsteps)

        @inbounds for j in 1:nblocks
            block = Ai[9j-8:9j, :]
            row = result[j, :]
            for k in 1:nsteps
                M = SMatrix{3,3}(reshape(block[:, k], 3, 3))
                row[k] = det(M)
            end
        end

        C[i] = result
    end

    return ScalarField(C, [;;], A.t, A.numElem, A.nsteps, :scalar, A.model)
end

"""
    inv(A::TensorField)

Matrix inverse of each 3×3 tensor block.

Returns: `TensorField`
"""
function inv(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds @views for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, 9nblocks, nsteps)

        @inbounds for j in 1:nblocks
            block = Ai[9j-8:9j, :]
            dest = result[9j-8:9j, :]
            for k in 1:nsteps
                M = SMatrix{3,3}(reshape(block[:, k], 3, 3))
                dest[:, k] .= inv(M)[:]
            end
        end

        C[i] = result
    end

    return TensorField(C, [;;], A.t, A.numElem, A.nsteps, :e, A.model)
end

function eigen(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A isa TensorField
        nsteps = A.nsteps
        C = Vector{Matrix{Float64}}(undef, length(A.A))
        c = Vector{Matrix{Float64}}(undef, length(A.A))
        @inbounds @views for (idx, Ai) in enumerate(A.A)
            nblocks = fld(size(Ai, 1), 9)
            eigenvectors = Matrix{Float64}(undef, 9nblocks, nsteps)
            eigenvalues = Matrix{Float64}(undef, 3nblocks, nsteps)
            @inbounds for j in 1:nblocks
                block = Ai[9j-8:9j, :]
                vec_dest = eigenvectors[9j-8:9j, :]
                val_dest = eigenvalues[3j-2:3j, :]
                for k in 1:nsteps
                    E = reshape(block[:, k], 3, 3)
                    vals, vecs = eigen(E, sortby=-)
                    for col in 1:3
                        v = view(vecs, :, col)
                        nrm = norm(v)
                        vecs[:, col] .= v ./ (nrm > 1e-8 ? nrm : 1.0)
                    end
                    vec_dest[:, k] .= vecs[:]
                    val_dest[:, k] .= vals
                end
            end
            C[idx] = eigenvectors
            c[idx] = eigenvalues
        end
        return VectorField(c, [;;], A.t, A.numElem, A.nsteps, :v3D, A.model), TensorField(C, [;;], A.t, A.numElem, A.nsteps, :tensor, A.model)
    else
        error("eigen(A::TensorField): A is not a TensorField.")
    end
end

function sqrt(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds @views for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, 9nblocks, nsteps)

        @inbounds for j in 1:nblocks
            block = Ai[9j-8:9j, :]
            dest = result[9j-8:9j, :]
            for k in 1:nsteps
                M = reshape(block[:, k], 3, 3)
                dest[:, k] .= sqrt(M)[:]
            end
        end

        C[i] = result
    end

    return TensorField(C, [;;], A.t, A.numElem, A.nsteps, :e, A.model)
end

function cbrt(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds @views for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, 9nblocks, nsteps)

        @inbounds for j in 1:nblocks
            block = Ai[9j-8:9j, :]
            dest = result[9j-8:9j, :]
            for k in 1:nsteps
                M = reshape(block[:, k], 3, 3)
                dest[:, k] .= cbrt(M)[:]
            end
        end

        C[i] = result
    end

    return TensorField(C, [;;], A.t, A.numElem, A.nsteps, :e, A.model)
end

function log(AA::TensorField)
    A = AA.A == [] ? nodesToElements(AA) : AA

    nsteps = A.nsteps
    n = length(A.A)
    C = Vector{Matrix{Float64}}(undef, n)

    @inbounds @views for i in 1:n
        Ai = A.A[i]
        nblocks = fld(size(Ai, 1), 9)
        result = Matrix{Float64}(undef, 9nblocks, nsteps)

        @inbounds for j in 1:nblocks
            block = Ai[9j-8:9j, :]
            dest = result[9j-8:9j, :]
            for k in 1:nsteps
                M = reshape(block[:, k], 3, 3)
                dest[:, k] .= log(M)[:]
            end
        end

        C[i] = result
    end

    return TensorField(C, [;;], A.t, A.numElem, A.nsteps, :e, A.model)
end

###############################################################################
#                                                                             #
#                      SystemMatrix operations                                #
#                                                                             #
###############################################################################

"""
    *(A::Union{SystemMatrix,Matrix}, B::Union{ScalarField,VectorField,TensorField})

Matrix–vector multiplication between a system matrix and a nodal vector field.

If the vector field is defined elementwise, it is automatically converted to
nodal representation before multiplication.

# Returns
- `VectorField` containing the nodal result with one time step.
"""
function *(A::Union{SystemMatrix,Matrix}, BB::Union{ScalarField,VectorField,TensorField})
    B = isNodal(BB) ? elementsToNodes(BB) : BB
    T = typeof(BB)

    #if B.a != [;;] && B.nsteps == 1
    type = B.type
    C = A isa Matrix ? A * B.a : A.A * B.a
    return T([], reshape(C, :, BB.nsteps), BB.t, [], BB.nsteps, type, B.model)
    #else
    #    error("*(A, B::Union{ScalarField,VectorField,TensorField}): vector field must be nodal with a single time step.")
    #end
end

#=
"""
    *(A::Union{SystemMatrix,Matrix}, B::ScalarField)

Matrix–vector multiplication between a system matrix and a nodal scalar field.

If the scalar field is defined elementwise, it is automatically converted to
nodal representation before multiplication.

# Returns
- `ScalarField` containing the nodal result with one time step.
"""
function *(A::Union{SystemMatrix,Matrix}, BB::ScalarField)
    B = BB.A != [] ? elementsToNodes(BB) : BB

    if B.a != [;;] && B.nsteps == 1
        C = A isa Matrix ? A * B.a : A.A * B.a
        return ScalarField([], reshape(C, :, 1), [0.0], [], 1, :scalar, B.model)
    else
        error("*(A, B::ScalarField): scalar field must be nodal with a single time step.")
    end
end
=#


###############################################################################
# Linear solves
###############################################################################

import Base:\
"""
    \\(A::Union{SystemMatrix,Matrix}, b::Union{ScalarField,VectorField,TensorField})

Solves the linear system `A * x = b` for a nodal scalar, vectoror tensor field right-hand side.

If the field is defined elementwise, it will be converted to nodal form
before solving.

# Returns
- `ScalarField` or `VectorField` or `TensorField` containing the solution.
"""
function \(A::Union{SystemMatrix,Matrix}, BB::Union{ScalarField,VectorField,TensorField})
    B = isNodal(BB) ? elementsToNodes(BB) : BB
    T = typeof(BB)

    #if B.a != [;;] && B.nsteps == 1
        C = A isa Matrix ? A \ B.a : A.A \ B.a
        return T([], reshape(C, :, 1), BB.t, [], BB.nsteps, BB.type, BB.model)
    #else
    #    error("\\(A, b::Union{ScalarField,VectorField,TensorField}): scalar field must be nodal with a single time step.")
    #end
end

#=
"""
    \\(A::Union{SystemMatrix,Matrix}, b::VectorField)

Solves the linear system `A * x = b` for a nodal vector field right-hand side.

If the vector field is defined elementwise, it is converted to nodal form
before solving.

# Returns
- `VectorField` containing the solution.
"""

function \(A::Union{SystemMatrix,Matrix}, BB::VectorField)
    B = BB.A != [] ? elementsToNodes(BB) : BB

    if B.a != [;;] && B.nsteps == 1
        C = A isa Matrix ? A \ B.a : A.A \ B.a
        return VectorField([], reshape(C, :, 1), [0.0], [], 1, B.type, B.model)
    else
        error("\\(A, b::VectorField): vector field must be nodal with a single time step.")
    end
end
=#

###############################################################################
# Sparse RHS solves
###############################################################################

"""
    ldiv_sparse!(X, K, F)

Solves the sparse linear system `K * X = F` column-by-column, where `F` is a
sparse matrix representing multiple right-hand sides.

# Returns
- Sparse matrix `X` containing the solution.
"""
function ldiv_sparse!(X::SparseMatrixCSC, K::Union{SystemMatrix,SparseMatrixCSC}, F::SparseMatrixCSC)
    Ffac = K isa SystemMatrix ? lu(K.A) : lu(K)
    n, m = size(F)
    x = zeros(n)
    b = zeros(n)

    I, J, V = Int[], Int[], Float64[]

    for j in 1:m
        idx, val = findnz(view(F, :, j))
        fill!(b, 0.0)
        b[idx] .= val
        ldiv!(x, Ffac, b)
        nz = findall(!iszero, x)
        append!(I, nz)
        append!(J, fill(j, length(nz)))
        append!(V, x[nz])
    end

    return sparse(I, J, V, n, m)
end


"""
    \\(K::Union{SystemMatrix,SparseMatrixCSC}, F::SparseMatrixCSC)

Solves a sparse linear system with multiple right-hand sides.

# Returns
- Sparse matrix containing the solution.
"""
function \(K::Union{SystemMatrix,SparseMatrixCSC}, F::SparseMatrixCSC)
    n = K isa SystemMatrix ? size(K.A, 1) : size(K, 1)
    X = spzeros(n, size(F, 2))
    return ldiv_sparse!(X, K, F)
end


###############################################################################
# Algebraic operations
###############################################################################

"""
    *(A::SystemMatrix, c::Number)
    *(c::Number, A::SystemMatrix)

Scalar multiplication of a system matrix.
"""
function *(A::SystemMatrix, b::Number)
    SystemMatrix(A.A * b, A.model)
end
function *(b::Number, A::SystemMatrix)
    SystemMatrix(A.A * b, A.model)
end


"""
    +(A::SystemMatrix, B::SystemMatrix)
    -(A::SystemMatrix, B::SystemMatrix)

Addition and subtraction of system matrices.
"""
function +(A::SystemMatrix, B::SystemMatrix)
    SystemMatrix(A.A + B.A, A.model)
end
function -(A::SystemMatrix, B::SystemMatrix)
    SystemMatrix(A.A - B.A, A.model)
end
function -(A::SystemMatrix)
    SystemMatrix(-A.A, A.model)
end


###############################################################################
# Optional but recommended utilities
###############################################################################

import Base:copy
"""
    copy(K::SystemMatrix)

Returns a deep copy of the system matrix.
"""
copy(K::SystemMatrix) = SystemMatrix(copy(K.A), K.model)


import Base:transpose,adjoint
import SparseArrays: sparse
"""
    transpose(K::SystemMatrix)
    adjoint(K::SystemMatrix)

Transpose / adjoint of a system matrix.
"""
transpose(K::SystemMatrix) = SystemMatrix(sparse(transpose(K.A)), K.test_model, K.model)
adjoint(K::SystemMatrix)   = SystemMatrix(sparse(adjoint(K.A)), K.test_model, K.model)


import LinearAlgebra: issymmetric
"""
    issymmetric(K::SystemMatrix)

Checks whether the system matrix is symmetric.
"""
issymmetric(K::SystemMatrix) = issymmetric(K.A)

import Base:getindex,setindex!,size,axes,eltype

"""
    getindex(K::SystemMatrix, I...)

Indexing operation for `SystemMatrix`.

Forwards all indexing operations to the underlying sparse matrix `K.A`,
allowing a `SystemMatrix` to be indexed in the same way as a
`SparseMatrixCSC`.

Examples include:
- `K[i, j]`
- `K[:, j]`, `K[i, :]`
- `K[a:b, c:d]`
- `K[v1, v2]` where `v1` and `v2` are index vectors.
"""
@inline getindex(K::SystemMatrix, I...) = getindex(K.A, I...)


"""
    setindex!(K::SystemMatrix, v, I...)

In-place assignment for `SystemMatrix`.

Forwards indexed assignment to the underlying sparse matrix `K.A`,
enabling modifications such as:
- `K[i, j] = v`
- `K[a:b, c:d] .= v`
- `K[v1, v2] .= submatrix`

Note that assignment follows the semantics and performance characteristics
of `SparseMatrixCSC`.
"""
@inline setindex!(K::SystemMatrix, v, I...) = setindex!(K.A, v, I...)


"""
    size(K::SystemMatrix)

Return the size of the system matrix.

Equivalent to `size(K.A)`.
"""
size(K::SystemMatrix) = size(K.A)


"""
    axes(K::SystemMatrix)

Return the valid index ranges for the system matrix.

Equivalent to `axes(K.A)`.
"""
axes(K::SystemMatrix) = axes(K.A)


"""
    eltype(K::SystemMatrix)

Return the element type of the system matrix.

Equivalent to `eltype(K.A)`.
"""
eltype(K::SystemMatrix) = eltype(K.A)

###############################################################################
#                                                                             #
#                      Transformation operations                              #
#                                                                             #
###############################################################################

function *(A::Transformation, B::Transformation)
    if A.non != B.non || A.dim != B.dim
        error("*(A::Transformation, B::Transformation): size missmatch non = $(A.non) ≠ $(B.non), dim = $(A.dim) ≠ $(B.dim).")
    end
    return Transformation(dropzeros(A.T * B.T), A.non, A.dim)
end

function *(A::Transformation, BB::VectorField)
    if BB.A != []
        B = elementsToNodes(BB)
    else
        B = BB
    end
    n = size(B.a, 1)
    non = A.non
    dim = A.dim
    if dim * non == n
        if B.A == []
            v = A.T * B.a
            return VectorField([], v, B.t, [], length(B.t), B.type, B.model)
        else
            error("*(A::Transformation, B::VectorField): B contains element data instead of nodal data.")
        end
    else
        error("*(A::Transformation, B::VectorField): size missmatch A.dim * A.non = $dim * $non ≠ $n = size(B.a, 1).")
    end
end

function *(BB::VectorField, A::Transformation)
    if BB.A != []
        B = elementsToNodes(BB)
    else
        B = BB
    end
    n = size(B.a, 1)
    non = A.non
    dim = A.dim
    if dim * non == n
        if B.A == []
            v = (B.a' * A.T)'
            return VectorField([], v, B.t, [], length(B.t), B.type, B.model)
        else
            error("*(B::VectorField, A::Transformation): B contains element data instead of nodal data.")
        end
    else
        error("*(B::VectorField, A::Transformation): size missmatch A.dim * A.non = $dim * $non ≠ $n = size(B.a, 1).")
    end
end

function *(A::Transformation, BB::TensorField)
    if BB.A != []
        B = elementsToNodes(BB)
    else
        B = BB
    end
    n = size(B.a, 1)
    m = size(B.a, 2)
    non = A.non
    dim = A.dim
    if B.A == []
        C = zeros(3non, 3)
        D = zeros(3non, 3)
        E = zeros(n, m)
        T = []
        I = []
        J = []
        V = Float64[]
        T1 = zeros(9)
        I0 = [1, 2, 3, 1, 2, 3, 1, 2, 3]
        J0 = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        if dim == 2
            for i in 1:non
                T1 = [A.T[2i-1, 2i-1], A.T[2i, 2i-1], 0, A.T[2i-1, 2i], A.T[2i, 2i], 0, 0, 0, 1]
                Idx = I0 .+ (3i - 3)
                Jdx = J0 .+ (3i - 3)
                append!(I, Idx)
                append!(J, Jdx)
                append!(V, T1)
            end
            fn(x, y) = y
            T = sparse(I, J, V, 3non, 3non, fn)
            dropzeros!(T)
        elseif dim == 3
            T = A.T
        else
            error("*(A::Transformation, B::TensorField): dim of A is $dim")
        end
        for k in 1:m
            for i in 1:non
                for j = 1:3
                    C[3i-2:3i, j] = B.a[9i-9+3j-2:9i-9+3j, k]
                end
            end
            D = T * C
            for i in 1:non
                for j = 1:3
                    E[9i-9+3j-2:9i-9+3j, k] = D[3i-2:3i, j]
                end
            end
        end
        return TensorField([], E, B.t, [], length(B.t), B.type, B.model)
    else
        error("*(A::Transformation, B::TensorField): B contains element data instead of nodal data.")
    end
end

function *(BB::TensorField, A::Transformation)
    if BB.A != []
        B = elementsToNodes(BB)
    else
        B = BB
    end
    n = size(B.a, 1)
    m = size(B.a, 2)
    non = A.non
    dim = A.dim
    if B.A == []
        C = zeros(3, 3non)
        D = zeros(3, 3non)
        E = zeros(n, m)
        T = []
        I = []
        J = []
        V = Float64[]
        T1 = zeros(9)
        I0 = [1, 2, 3, 1, 2, 3, 1, 2, 3]
        J0 = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        if dim == 2
            for i in 1:non
                T1 = [A.T[2i-1, 2i-1], A.T[2i, 2i-1], 0, A.T[2i-1, 2i], A.T[2i, 2i], 0, 0, 0, 1]
                Idx = I0 .+ (3i - 3)
                Jdx = J0 .+ (3i - 3)
                append!(I, Idx)
                append!(J, Jdx)
                append!(V, T1)
            end
            fn(x, y) = y
            T = sparse(I, J, V, 3non, 3non, fn)
            dropzeros!(T)
        elseif dim == 3
            T = A.T
        else
            error("*(A::Transformation, B::TensorField): dim of A is $dim")
        end
        for k in 1:m
            for i in 1:non
                for j = 1:3
                    C[1:3, 3i-3+j] = B.a[9i-9+3j-2:9i-9+3j, k]
                end
            end
            D = C * T
            for i in 1:non
                for j = 1:3
                    E[9i-9+3j-2:9i-9+3j, k] = D[1:3, 3i-3+j]
                end
            end
        end
        return TensorField([], E, B.t, [], length(B.t), B.type, B.model)
    else
        error("*(B::TensorField, A::Transformation): B contains element data instead of nodal data.")
    end
end

function transpose(A::Transformation)
    return Transformation(transpose(A.T), A.non, A.dim)
end

function adjoint(A::Transformation)
    return Transformation(adjoint(A.T), A.non, A.dim)
end

function *(A::Transformation, B::SystemMatrix)
    n = size(B.A, 1)
    non = A.non
    dim = A.dim
    if dim * non == n
        return SystemMatrix(dropzeros(A.T * B.A), B.model)
    else
        error("*(A::Transformation, B::SystemMatrix): size missmatch dim * non = $dim * $non ≠ $n.")
    end
end

function *(B::SystemMatrix, A::Transformation)
    n = size(B.A, 1)
    non = A.non
    dim = A.dim
    if dim * non == n
        return SystemMatrix(dropzeros(B.A * A.T), B.model)
    else
        error("*(A::Transformation, B::SystemMatrix): size missmatch dim * non = $dim * $non ≠ $n.")
    end
end

###############################################################################
#                                                                             #
#                  Differential operators operations                          #
#                                                                             #
###############################################################################

"""
    ∘(A::Union{ScalarField,VectorField}, D::Function)

Right application of differential operator `D` to field `A`.
- If `D == ∇` and `A` is `ScalarField`: returns `grad(A)`.
- If `D == ∇` and `A` is `VectorField`: returns `grad(A)`.

Returns: `VectorField` or `TensorField`

# Examples
```julia
# 3D (assumes `problem` and a "body" physical group are defined)
S = scalarField(problem, [field("body", f=(x,y,z)->x*y)])
G = S ∘ ∇      # grad of scalar field
V = vectorField(problem, [field("body", fx=x->x, fy=y->y, fz=z->z)])
H = V ∘ ∇      # grad of vector field (tensor)
```
"""
function ∘(A::Union{ScalarField,VectorField}, D::Function)
    if D == ∇
        return grad(A)
    else
        error("⋅(A::Union{ScalarField,VectorField}, D::Function): D is not a differential operator.")
    end
end

"""
    ∘(D::Function, A::Union{ScalarField,VectorField})

Left application of differential operator `D` to field `A`.
- If `D == ∇` and `A` is `ScalarField`: returns `grad(A)`.
- If `D == ∇` and `A` is `VectorField`: returns `grad(A)'` (transpose).

Returns: `VectorField` or `TensorField`

# Examples
```julia
# 3D (assumes `problem` and a "body" physical group are defined)
V = vectorField(problem, [field("body", fx=x->x, fy=y->y, fz=z->z)])
T = ∇ ∘ V      # equals grad(V)'
```
"""
function ∘(D::Function, A::Union{ScalarField,VectorField})
    if D == ∇
        if A isa ScalarField
            return grad(A)
        elseif A isa VectorField
            return grad(A)'
        end
    else
        error("⋅(A::Union{ScalarField,VectorField}, D::Function): D is not a differential operator.")
    end
end

"""
    ⋅(A::Union{VectorField,TensorField}, D::Function)

Right contraction with the differential operator. With `D == ∇`:
- If `A` is `VectorField`: returns `div(A)` (scalar field).
- If `A` is `TensorField`: returns `div(A)` (vector field).

Returns: `ScalarField` or `VectorField`

# Examples
```julia
# 3D (assumes `problem` and a "body" physical group are defined)
V = vectorField(problem, [field("body", fx=x->x, fy=y->y, fz=z->z)])
divV = V ⋅ ∇   # ScalarField
```
"""
function ⋅(A::Union{VectorField,TensorField}, D::Function)
    if D == ∇
        return div(A)
    else
        error("⋅(A::Union{VectorField,TensorField}, D::Function): D is not a differential operator.")
    end
end

"""
    ⋅(D::Function, A::Union{VectorField,TensorField})

Left contraction with the differential operator. With `D == ∇`:
- If `A` is `VectorField`: returns `div(A)`.
- If `A` is `TensorField`: returns `div(A')`.

Returns: `ScalarField` or `VectorField`

# Examples
```julia
# 3D (assumes `problem` and a "body" physical group are defined)
T = tensorField(problem, [field("body", fz=z->z)])
DV = ∇ ⋅ T     # VectorField (divergence of tensor)
```
"""
function ⋅(D::Function, A::Union{VectorField,TensorField})
    if D == ∇
        if A isa VectorField
            return div(A)
        elseif A isa TensorField
            return div(A')
        end
    else
        error("⋅(D::Function, A::Union{VectorField,TensorField}): D is not a differential operator.")
    end
end

"""
    ×(D::Function, A::VectorField)

Left curl. With `D == ∇`, returns `curl(A)`.

Returns: `VectorField`

# Examples
```julia
# 3D (assumes `problem` and a "body" physical group are defined)
V = vectorField(problem, [field("body", fx=x->0, fy=x->x, fz=z->0)])
C = ∇ × V
```
"""
function ×(D::Function, A::VectorField)
    if D == ∇
        return curl(A)
    else
        error("×(D::Function, A::VectorField): D is not a differential operator.")
    end
end

"""
    ×(A::VectorField, D::Function)

Right curl with sign convention. With `D == ∇`, returns `-curl(A)`.

Returns: `VectorField`

# Examples
```julia
# 3D (assumes `problem` and a "body" physical group are defined)
V = vectorField(problem, [field("body", fx=x->0, fy=x->x, fz=z->0)])
Cneg = V × ∇   # -curl(V)
```
"""
function ×(A::VectorField, D::Function)
    if D == ∇
        return curl(A) * (-1)
    else
        error("×(A::VectorField, D::Function): D is not a differential operator.")
    end
end
