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

"""
    *(A::ScalarField, B::ScalarField)

Performs element-wise multiplication of two `ScalarField` objects on the same set of elements.

Returns: `ScalarField`

# Examples
```julia
C = A * B
```
"""
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
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i])
        D = zeros(n, nsteps)
        if n != sz
            #D = zeros(n, nsteps)
            sz = n
        end
        for j in 1:n
            for k in 1:nsteps
                D[j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][j, k]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, num, A.nsteps, :scalar, A.model)
end

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
        error("/(ScalarField, ScalarField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
    end
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i])
        D = zeros(n, nsteps)
        if n != sz
            #D = zeros(n, nsteps)
            sz = n
        end
        for j in 1:n
            for k in 1:nsteps
                D[j, k] = A.A[indS[i]][j, k] / B.A[indT[i]][j, k]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, num, A.nsteps, :scalar, A.model)
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
            #if length(A.A) != length(B.A)
            #    error("+(A::ScalarField, B::ScalarField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            #end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("+(A::ScalarField, B::ScalarField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            idx01 = findall(j -> j in sec, A.numElem)
            idx02 = findall(j -> j in sec, B.numElem)
            for i in eachindex(idx01)
                D = A.A[idx01[i]] + B.A[idx02[i]]
                append!(num, A.numElem[idx01[i]] #=sec[i]=#)
                push!(C, D)
            end
            idx1 = findall(j -> j in dif1, A.numElem)
            for i in idx1 #eachindex(dif1)
                D = A.A[i]
                append!(num, A.numElem[i] #=dif1[i]=#)
                push!(C, D)
            end
            idx2 = findall(j -> j in dif2, B.numElem)
            for i in idx2 #eachindex(dif2)
                D = B.A[i]
                append!(num, B.numElem[i] #=dif2[i]=#)
                push!(C, D)
            end
            a = [;;]
            return ScalarField(C, a, A.t, num, A.nsteps, A.type, A.model)
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
            #if length(A.A) != length(B.A)
            #    error("+(A::ScalarField, B::ScalarField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            #end
            nsteps = A.nsteps
            nsteps2 = B.nsteps
            if nsteps != nsteps2
                error("-(A::ScalarField, B::ScalarField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
            end
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            idx01 = findall(j -> j in sec, A.numElem)
            idx02 = findall(j -> j in sec, B.numElem)
            for i in eachindex(idx01)
                D = A.A[idx01[i]] - B.A[idx02[i]]
                append!(num, A.numElem[idx01[i]] #=sec[i]=#)
                push!(C, D)
            end
            idx1 = findall(j -> j in dif1, A.numElem)
            for i in idx1 #eachindex(dif1)
                #display("dif1 = $dif1")
                D = A.A[i]
                append!(num, A.numElem[i] #=dif1[i]=#)
                push!(C, D)
            end
            idx2 = findall(j -> j in dif2, B.numElem)
            for i in idx2 #eachindex(dif2)
                D = -B.A[i]
                append!(num, B.numElem[i] #=dif2[i]=#)
                push!(C, D)
            end
            a = [;;]
            return ScalarField(C, a, A.t, num, A.nsteps, A.type, A.model)
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
function *(AA::ScalarField, b::Number)
    if isNodal(AA)
        A = nodesToElements(AA)
    else
        A = AA
    end
    C = []
    for i in 1:length(A.A)
        D = A.A[i] * b
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
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
function *(b::Number, AA::ScalarField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    C = []
    for i in 1:length(A.A)
        D = A.A[i] * b
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
end

function /(AA::ScalarField, b::Number)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    C = []
    for i in 1:length(A.A)
        D = A.A[i] / b
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
end

function /(b::Number, AA::ScalarField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    C = []
    for i in 1:length(A.A)
        D = b ./ A.A[i]
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
end

function log(AA::ScalarField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    nsteps = A.nsteps
    C = []
    for i in 1:length(A.A)
        n = length(A.A[i])
        D = zeros(n, nsteps)
        for j in 1:n
            for k in 1:nsteps
                D[j, k] = log(A.A[i][j, k])
            end
        end
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
end

function sqrt(AA::ScalarField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    nsteps = A.nsteps
    C = []
    for i in 1:length(A.A)
        n = length(A.A[i])
        D = zeros(n, nsteps)
        for j in 1:n
            for k in 1:nsteps
                D[j, k] = sqrt(A.A[i][j, k])
            end
        end
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, A.numElem, A.nsteps, A.type, A.model)
end

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
    if BB.type != :v3D
        error("*(AA::ScalarField, BB::VectorField): Vector field must be 3 dimansional.")
    end
    if A.nsteps != B.nsteps || A.nsteps != 1
        error("*(AA::ScalarField, BB::VectorField): Number of nsteps in AA must be one or equal to BB.nsteps.")
    end
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i]) ÷ 3
        D = zeros(3n, nsteps)
        if n != sz
            sz = n
        end
        for j in 1:n
            for k in 1:nsteps
                kk = AA.nsteps > 1 ? k : 1
                D[3j-2:3j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][3j-2:3j, kk]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return VectorField(C, a, B.t, num, B.nsteps, BB.type, B.model)
end

"""
    *(B::VectorField, A::ScalarField)

Scales a `VectorField` by a `ScalarField` element-wise on matching elements.

Returns: `VectorField`
"""
function *(BB::VectorField, AA::ScalarField)
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
    if BB.type != :v3D
        error("*(BB::VectorField, AA::ScalarField): Vector field must be 3 dimansional.")
    end
    if A.nsteps != B.nsteps || A.nsteps != 1
        error("*(BB::VectorField, AA::ScalarField): Number of nsteps in AA must be one or equal to BB.nsteps.")
    end
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i]) ÷ 3
        D = zeros(3n, nsteps)
        if n != sz
            sz = n
        end
        for j in 1:n
            for k in 1:nsteps
                kk = AA.nsteps > 1 ? k : 1
                D[3j-2:3j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][3j-2:3j, kk]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return VectorField(C, a, B.t, num, B.nsteps, BB.type, B.model)
end

"""
    /(B::VectorField, A::ScalarField)

Divides a `VectorField` by a `ScalarField` element-wise on matching elements.

Returns: `VectorField`
"""
function /(BB::VectorField, AA::ScalarField)
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
    if BB.type != :v3D
        error("/(BB::VectorField, AA::ScalarField): Vector field must be 3 dimansional.")
    end
    if A.nsteps != B.nsteps || A.nsteps != 1
        error("/(BB::VectorField, AA::ScalarField): Number of nsteps in AA must be one or equal to BB.nsteps.")
    end
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i]) ÷ 3
        D = zeros(3n, nsteps)
        if n != sz
            sz = n
        end
        for j in 1:n
            for k in 1:nsteps
                kk = AA.nsteps > 1 ? k : 1
                D[3j-2:3j, k] = A.A[indS[i]][j, k] / B.A[indT[i]][3j-2:3j, kk]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return VectorField(C, a, B.t, num, B.nsteps, BB.type, B.model)
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
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            idx01 = findall(j -> j in sec, A.numElem)
            idx02 = findall(j -> j in sec, B.numElem)
            for i in eachindex(idx01)
                D = A.A[idx01[i]] + B.A[idx02[i]]
                append!(num, A.numElem[idx01[i]] #=sec[i]=#)
                push!(C, D)
            end
            idx1 = findall(j -> j in dif1, A.numElem)
            for i in idx1 #eachindex(dif1)
                D = A.A[i]
                append!(num, A.numElem[i] #=dif1[i]=#)
                push!(C, D)
            end
            idx2 = findall(j -> j in dif2, B.numElem)
            for i in idx2 #eachindex(dif2)
                D = B.A[i]
                append!(num, B.numElem[i] #=dif2[i]=#)
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, num, A.nsteps, A.type, A.model)
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
            sec = intersect(A.numElem, B.numElem)
            ind1 = []
            ind2 = []
            sizehint!(ind1, length(sec))
            sizehint!(ind2, length(sec))
            for i in sec
                append!(ind1, findall(j -> j == i, A.numElem))
                append!(ind2, findall(j -> j == i, B.numElem))
            end
            dif1 = setdiff(A.numElem, B.numElem)
            ind3 = []
            sizehint!(ind3, length(dif1))
            for i in dif1
                append!(ind3, findall(j -> j == i, A.numElem))
            end
            dif2 = setdiff(B.numElem, A.numElem)
            ind4 = []
            sizehint!(ind4, length(dif2))
            for i in dif2
                append!(ind4, findall(j -> j == i, B.numElem))
            end
            C = []
            num = []
            sizehint!(C, length(sec) + length(dif1) + length(dif2))
            sizehint!(num, length(sec) + length(dif1) + length(dif2))
            idx01 = findall(j -> j in sec, A.numElem)
            idx02 = findall(j -> j in sec, B.numElem)
            for i in eachindex(idx01)
                D = A.A[idx01[i]] - B.A[idx02[i]]
                append!(num, A.numElem[idx01[i]] #=sec[i]=#)
                push!(C, D)
            end
            idx1 = findall(j -> j in dif1, A.numElem)
            for i in idx1 #eachindex(dif1)
                D = A.A[i]
                append!(num, A.numElem[i] #=dif1[i]=#)
                push!(C, D)
            end
            idx2 = findall(j -> j in dif2, B.numElem)
            for i in idx2 #eachindex(dif2)
                D = -B.A[i]
                append!(num, B.numElem[i] #=dif2[i]=#)
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, num, A.nsteps, A.type, A.model)
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
    if length(A.A) != 0
        if A.type == :v3D || A.type == :v2D || A.type == :v3D || A.type == :v2D || true
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                D = A.A[i] * b
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, A.numElem, A.nsteps, A.type, A.model)
        else
            error("*(A::VectorField, b): VectorField type ($(A.type) or $(b.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0
        return VectorField([], A.a .* b, A.t, [], A.nsteps, A.type, A.model)
    else
        error("*(VectorField, b): internal error")
    end
end

function *(b::Number, A::VectorField)
    if length(A.A) != 0
        if A.type == :u3D || A.type == :u2D || A.type == :f3D || A.type == :f2D || true
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                D = A.A[i] * b
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, A.numElem, A.nsteps, A.type, A.model)
        else
            error("*(A::VectorField, b): VectorField type ($(A.type) or $(b.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0
        return VectorField([], A.a .* b, A.t, [], A.nsteps, A.type, A.model)
    else
        error("*(b, VectorField): internal error")
    end
end

function /(A::VectorField, b::Number)
    if length(A.A) != 0
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            D = A.A[i] / b
            push!(C, D)
        end
        a = [;;]
        return VectorField(C, a, A.t, A.numElem, A.nsteps, A.type, A.model)
    elseif length(A.a) != 0
        return VectorField([], A.a / b, A.t, [], A.nsteps, A.type, A.model)
    else
        error("/(VectorField, b): internal error")
    end
end

function dot(AA::VectorField, BB::VectorField)
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
        error("*(VectorField, VectorField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
    end
    if A.type != :v3D && B.type != :v3D
        error("*(AA::VectorField, BB::VectorField): AA and BB must be 3D vectors.")
    end
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i])
        D = zeros(n ÷ 3, nsteps)
        if n != sz
            #D = zeros(n, nsteps)
            sz = n
        end
        for j in 1:n ÷ 3
            for k in 1:nsteps
                D[j, k] = A.A[indS[i]][3j-2, k] * B.A[indT[i]][3j-2, k] +
                        A.A[indS[i]][3j-1, k] * B.A[indT[i]][3j-1, k] +
                        A.A[indS[i]][3j,   k] * B.A[indT[i]][3j,   k]
                 #A.A[indS[i]][j, k] * B.A[indT[i]][j, k]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, num, A.nsteps, :scalar, A.model)
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
        error("∘(VectorField, VectorField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
    end
    if A.type != :v3D && B.type != :v3D
        error("∘(AA::VectorField, BB::VectorField): AA and BB must be 3D vectors.")
    end
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i])
        D = zeros(3n, nsteps)
        if n != sz
            #D = zeros(n, nsteps)
            sz = n
        end
        for j in 1:(n ÷ 3)
            for k in 1:nsteps
                ax = A.A[indS[i]][3j-2, k]
                ay = A.A[indS[i]][3j-1, k]
                az = A.A[indS[i]][3j,   k]

                bx = B.A[indT[i]][3j-2, k]
                by = B.A[indT[i]][3j-1, k]
                bz = B.A[indT[i]][3j,   k]

                D[9j-8,k] = ax * bx
                D[9j-7,k] = ay * bx
                D[9j-6,k] = az * bx

                D[9j-5,k] = ax * by
                D[9j-4,k] = ay * by
                D[9j-3,k] = az * by

                D[9j-2,k] = ax * bz
                D[9j-1,k] = ay * bz
                D[9j  ,k] = az * bz
                #D[j, k] = A.A[indS[i]][3j-2, k] * B.A[indT[i]][3j-2, k] +
                #        A.A[indS[i]][3j-1, k] * B.A[indT[i]][3j-1, k] +
                #        A.A[indS[i]][3j,   k] * B.A[indT[i]][3j,   k]
                 #A.A[indS[i]][j, k] * B.A[indT[i]][j, k]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return TensorField(C, a, A.t, num, A.nsteps, :e, A.model)
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
        error("×(VectorField, VectorField): nsteps od A and B are not equal ($(A.nsteps) != $(B.nsteps)")
    end
    if A.type != :v3D && B.type != :v3D
        error("×(AA::VectorField, BB::VectorField): AA and BB must be 3D vectors.")
    end
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i])
        D = zeros(n, nsteps)
        if n != sz
            #D = zeros(n, nsteps)
            sz = n
        end
        for j in 1:(n ÷ 3)
            for k in 1:nsteps
                e = SVector{3}(A.A[indS[i]][3j-2:3j, k])
                f = SVector{3}(B.A[indT[i]][3j-2:3j, k])
                g = e × f

                D[3j-2:3j,k] .= g[1:3]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return VectorField(C, a, A.t, num, A.nsteps, :v3D, A.model)
end

"""
    norm(A::VectorField)

Element-wise Euclidean norm of a `VectorField`.

Returns: `ScalarField`
"""
function norm(AA::VectorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    nsteps = A.nsteps
    C = []
    for i in 1:length(A.A)
        n = length(A.A[i]) ÷ 3
        D = zeros(n, nsteps)
        for j in 1:n
            for k in 1:nsteps
                D[j, k] = norm(reshape(A.A[i][3j-2:3j, k], 3, 1))
            end
        end
        push!(C, D)
    end
    a = [;;]
    return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
end

"""
    diagm(A::VectorField)

Creates a diagonal `TensorField` from a `VectorField` (dim=3), i.e., places vector
components on the tensor diagonal for each node/element.

Returns: `TensorField`
"""
function diagm(AA::VectorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A isa VectorField && A.model.dim == 3
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 3
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    e = reshape(A.A[i][3j-2:3j, k], 3, 1)
                    F = [e[1] 0 0; 0 e[2] 0; 0 0 e[3]]
                    D[9j-8:9j, k] = reshape(F, 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
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
            error("*(A::TensorField, B::TensorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
        end
        if A.numElem != B.numElem
            error("*(A::TensorField, B::TensorField): tensor fields are not compatible.")
        end
        nsteps = A.nsteps
        nsteps2 = B.nsteps
        if nsteps != nsteps2
            error("*(A::TensorField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
        end
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            m = length(B.A[i]) ÷ 9
            if n != m
                error("*(A::TensorField, B::TensorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
            end
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = reshape(reshape(A.A[i][9j-8:9j, k], 3, 3) * reshape(B.A[i][9j-8:9j, k], 3, 3), 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("*(A::TensorField, B::TensorField): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

"""
    ⋅(A::TensorField, B::TensorField)

Element-wise (Hadamard) product followed by summation of all components, yielding a
scalar per tensor (i.e., Frobenius inner product).

Returns: `ScalarField`
"""
function ⋅(AA::TensorField, BB::TensorField)
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
            error("*(A::TensorField, B::TensorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
        end
        if A.numElem != B.numElem
            error("*(A::TensorField, B::TensorField): tensor fields are not compatible.")
        end
        nsteps = A.nsteps
        nsteps2 = B.nsteps
        if nsteps != nsteps2
            error("*(A::TensorField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
        end
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            m = length(B.A[i]) ÷ 9
            if n != m
                error("*(A::TensorField, B::TensorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
            end
            D = zeros(n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    E = reshape(A.A[i][9j-8:9j, k], 3, 3) .* reshape(B.A[i][9j-8:9j, k], 3, 3)
                    D[j, k] = sum(E)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
    else
        error("*(A::TensorField, B::TensorField): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function +(AA::TensorField, BB::TensorField)
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
    if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F) || true
        if length(A.A) != length(B.A)
            error("+(A::TensorField, B::TensorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
        end
        nsteps = A.nsteps
        nsteps2 = B.nsteps
        if nsteps != nsteps2
            error("+(A::TensorField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
        end
        sec = intersect(A.numElem, B.numElem)
        ind1 = []
        ind2 = []
        sizehint!(ind1, length(sec))
        sizehint!(ind2, length(sec))
        for i in sec
            append!(ind1, findall(j -> j == i, A.numElem))
            append!(ind2, findall(j -> j == i, B.numElem))
        end
        dif1 = setdiff(A.numElem, B.numElem)
        ind3 = []
        sizehint!(ind3, length(dif1))
        for i in dif1
            append!(ind3, findall(j -> j == i, A.numElem))
        end
        dif2 = setdiff(B.numElem, A.numElem)
        ind4 = []
        sizehint!(ind4, length(dif2))
        for i in dif2
            append!(ind4, findall(j -> j == i, B.numElem))
        end
        C = []
        num = []
        sizehint!(C, length(sec) + length(dif1) + length(dif2))
        sizehint!(num, length(sec) + length(dif1) + length(dif2))
        for i in eachindex(sec)
            #n = length(A.A[i]) ÷ 9
            #m = length(B.A[i]) ÷ 9
            #if n != m
            #    error("+(A::TensorField, B::TensorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
            #end
            D = A.A[i] + B.A[i]
            append!(num, sec[i])
            push!(C, D)
        end
        for i in eachindex(dif1)
            D = A.A[i]
            append!(num, dif1[i])
            push!(C, D)
        end
        for i in eachindex(dif2)
            D = B.A[i]
            append!(num, dif2[i])
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, num, A.nsteps, :e, A.model)
    else
        error("+(A::TensorField, B::TensorField): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function -(AA::TensorField, BB::TensorField)
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
    if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F) || true
        if length(A.A) != length(B.A)
            error("-(A::TensorField, B::TensorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
        end
        nsteps = A.nsteps
        nsteps2 = B.nsteps
        if nsteps != nsteps2
            error("-(A::TensorField, B::TensorField): nsteps of A=$(A.nsteps) != nsteps of B=$(B.nsteps)")
        end
        sec = intersect(A.numElem, B.numElem)
        ind1 = []
        ind2 = []
        sizehint!(ind1, length(sec))
        sizehint!(ind2, length(sec))
        for i in sec
            append!(ind1, findall(j -> j == i, A.numElem))
            append!(ind2, findall(j -> j == i, B.numElem))
        end
        dif1 = setdiff(A.numElem, B.numElem)
        ind3 = []
        sizehint!(ind3, length(dif1))
        for i in dif1
            append!(ind3, findall(j -> j == i, A.numElem))
        end
        dif2 = setdiff(B.numElem, A.numElem)
        ind4 = []
        sizehint!(ind4, length(dif2))
        for i in dif2
            append!(ind4, findall(j -> j == i, B.numElem))
        end
        C = []
        num = []
        sizehint!(C, length(sec) + length(dif1) + length(dif2))
        sizehint!(num, length(sec) + length(dif1) + length(dif2))
        for i in eachindex(sec)
            D = A.A[i] - B.A[i]
            append!(num, sec[i])
            push!(C, D)
        end
        for i in eachindex(dif1)
            D = A.A[i]
            append!(num, dif1[i])
            push!(C, D)
        end
        for i in eachindex(dif2)
            D = -B.A[i]
            append!(num, dif2[i])
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, num, A.nsteps, :e, A.model)
    else
        error("-(A::TensorField, B::TensorField): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function -(A::TensorField)
    return A * (-1)
end

function *(AA::TensorField, BB::VectorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if BB.A == []
        b = nodesToElements(BB)
    else
        b = BB
    end
    #if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
    #if b.model.dim == 3
    if b.type == :v3D
        if length(A.A) != length(b.A)
            error("*(A::TensorField, b::VectorField): size of A=$(length(A.A)) != size of b=$(length(b.A))")
        end
        if A.numElem != b.numElem
            error("*(A::TensorField, b::VectorField): tensor fields are not compatible.")
        end
        nsteps = A.nsteps
        nsteps2 = b.nsteps
        if nsteps != nsteps2
            error("*(A::TensorField, b::VectorField): nsteps of A=$(A.nsteps) != nsteps of b=$(b.nsteps)")
        end
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            m = length(b.A[i]) ÷ 3
            if n != m
                error("*(A::TensorField, b::VectorField): size of A.A[$i]=$(9n) != size of 3*b.A[$j]=$(3m)")
            end
            D = zeros(3n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[3j-2:3j, k] = reshape(reshape(A.A[i][9j-8:9j, k], 3, 3) * reshape(b.A[i][3j-2:3j, k], 3, 1), 3, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return VectorField(C, a, A.t, A.numElem, A.nsteps, b.type, A.model)
    else
        error("*(A::TensorField, b::VectorField): Multiply by 2D vector is not yet implemented.")
    end
end

function *(BB::VectorField, AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if BB.A == []
        b = nodesToElements(BB)
    else
        b = BB
    end
    #if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
    #if b.model.dim == 3
    if b.type == :v3D
        if length(A.A) != length(b.A)
            error("*(A::TensorField, b::VectorField): size of A=$(length(A.A)) != size of b=$(length(b.A))")
        end
        if A.numElem != b.numElem
            error("*(A::TensorField, b::VectorField): tensor fields are not compatible.")
        end
        nsteps = A.nsteps
        nsteps2 = b.nsteps
        if nsteps != nsteps2
            error("*(A::TensorField, b::VectorField): nsteps of A=$(A.nsteps) != nsteps of b=$(b.nsteps)")
        end
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            m = length(b.A[i]) ÷ 3
            if n != m
                error("*(A::TensorField, b::VectorField): size of A.A[$i]=$(9n) != size of 3*b.A[$j]=$(3m)")
            end
            D = zeros(3n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[3j-2:3j, k] = reshape(reshape(b.A[i][3j-2:3j, k], 3, 1)' * reshape(A.A[i][9j-8:9j, k], 3, 3), 3, 1)'
                end
            end
            push!(C, D)
        end
        a = [;;]
        return VectorField(C, a, A.t, A.numElem, A.nsteps, b.type, A.model)
    else
        error("*(A::TensorField, b::VectorField): Multiply by 2D vector is not yet implemented.")
    end
end

function *(AA::TensorField, b::Number)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A.type == :s || A.type == :e || true
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            D = A.A[i] * b
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("*(A::TensorField, b): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function *(b::Number, AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A.type == :s || A.type == :e || true
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            D = A.A[i] * b
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("*(A::TensorField, b): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function /(AA::TensorField, b::Number)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A.type == :s || A.type == :e || true
        C = []
        for i in 1:length(A.A)
            D = A.A[i] / b
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("/(A::TensorField, b): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function *(AA::ScalarField, BB::TensorField)
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
    sz = 0
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i]) ÷ 9
        if n != sz
            D = zeros(9n, nsteps)
            sz = n
        end
        for j in 1:n
            for k in 1:nsteps
                D[9j-8:9j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][9j-8:9j, k]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return TensorField(C, a, B.t, num, B.nsteps, :e, B.model)
end

function *(BB::TensorField, AA::ScalarField)
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
    sz = 0
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i]) ÷ 9
        if n != sz
            D = zeros(9n, nsteps)
            sz = n
        end
        for j in 1:n
            for k in 1:nsteps
                D[9j-8:9j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][9j-8:9j, k]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return TensorField(C, a, B.t, num, B.nsteps, :e, B.model)
end

function /(BB::TensorField, AA::ScalarField)
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
    sz = 0
    nsteps = B.nsteps
    sec = intersect(B.numElem, A.numElem)
    indS = []
    indT = []
    sizehint!(indS, length(sec))
    sizehint!(indT, length(sec))
    for i in sec
        append!(indS, findall(j -> j == i, A.numElem))
        append!(indT, findall(j -> j == i, B.numElem))
    end
    C = []
    num = []
    sizehint!(C, length(sec))
    sizehint!(num, length(sec))
    D = []
    for i in eachindex(sec)
        n = length(B.A[i]) ÷ 9
        if n != sz
            D = zeros(9n, nsteps)
            sz = n
        end
        for j in 1:n
            for k in 1:nsteps
                D[9j-8:9j, k] = B.A[indT[i]][9j-8:9j, k] / A.A[indS[i]][j, k]
            end
        end
        append!(num, sec[i])
        push!(C, D)
    end
    a = [;;]
    return TensorField(C, a, B.t, num, B.nsteps, :e, B.model)
end

"""
    transpose(A::TensorField)

Transposes each 3×3 tensor block.

Returns: `TensorField`
"""
function transpose(A::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    #if A.type == :s || A.type == :e || A.type == :F
    if true
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = reshape(transpose(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("transpose(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
    end
end

"""
    adjoint(A::TensorField)

Adjoint (conjugate transpose) of each 3×3 tensor block.

Returns: `TensorField`
"""
function adjoint(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if true
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = reshape(adjoint(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("adjoint(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
    end
end

"""
    unitTensor(A::TensorField)

Creates an identity tensor field (I) with the same element structure and time steps as `A`.

Returns: `TensorField`
"""
function unitTensor(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if true #A.type == :s || A.type == :e || A.type == :F
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = reshape([1 0 0; 0 1 0; 0 0 1], 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("unit(A::TensorField): TensorField type ($(A.type) is not yet implemented.")
    end
end

"""
    trace(A::TensorField)

Computes the trace of each 3×3 tensor block.

Returns: `ScalarField`
"""
function trace(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    sz = 0
    if true #A.type == :s || A.type == :e || A.type == :F
        nsteps = A.nsteps
        C = []
        D = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            if sz != n
                D = zeros(n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    trace = A.A[i][9j-8, k] + A.A[i][9j-4, k] + A.A[i][9j, k]
                    D[j, k] = trace
                end
            end
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
    else
        error("trace(A::TensorField): TensorField type ($(A.type) is not yet implemented.")
    end
end

"""
    det(A::TensorField)

Computes the determinant of each 3×3 tensor block.

Returns: `ScalarField`
"""
function det(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    sz = 0
    if true #A.type == :s || A.type == :e || A.type == :F
        nsteps = A.nsteps
        C = []
        D = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            if sz != n
                D = zeros(n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    d = LinearAlgebra.det(reshape(A.A[i][9j-8:9j, k], 3, 3))
                    D[j, k] = d
                end
            end
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :scalar, A.model)
    else
        error("det(A::TensorField): TensorField type ($(A.type) is not yet implemented.")
    end
end

"""
    inv(A::TensorField)

Matrix inverse of each 3×3 tensor block.

Returns: `TensorField`
"""
function inv(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A isa TensorField
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = reshape(inv(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("inv(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
    end
end

function eigen(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A isa TensorField
        nsteps = A.nsteps
        C = []
        c = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            D = zeros(9n, nsteps)
            h = zeros(3n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    E = reshape(A.A[i][9j-8:9j, k], 3, 3)
                    f, G = eigen(E, sortby=-)
                    #GG = mapslices(v -> v ./ norm(v), G, dims=1)
                    G[:,1] ./= norm(G[:,1]) > 1e-8 ? norm(G[:,1]) : 1
                    G[:,2] ./= norm(G[:,2]) > 1e-8 ? norm(G[:,2]) : 1
                    G[:,3] ./= norm(G[:,3]) > 1e-8 ? norm(G[:,3]) : 1
                    D[9j-8:9j, k] = reshape(G, 9, 1)
                    h[3j-2:3j, k] = reshape(f, 3, 1)
                end
            end
            push!(C, D)
            push!(c, h)
        end
        a = [;;]
        return VectorField(c, a, A.t, A.numElem, A.nsteps, :v3D, A.model), TensorField(C, a, A.t, A.numElem, A.nsteps, :tensor, A.model)
    else
        error("eigen(A::TensorField): A is not a TensorField.")
    end
end

function sqrt(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A.type == :s || A.type == :e || A.type == :F
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = reshape(sqrt(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("sqrt(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
    end
end

function log(AA::TensorField)
    if AA.A == []
        A = nodesToElements(AA)
    else
        A = AA
    end
    if A.type == :s || A.type == :e
        nsteps = A.nsteps
        C = []
        for i in 1:length(A.A)
            n = length(A.A[i]) ÷ 9
            D = zeros(9n, nsteps)
            for j in 1:n
                for k in 1:nsteps
                    D[9j-8:9j, k] = reshape(log(reshape(A.A[i][9j-8:9j, k], 3, 3)), 9, 1)
                end
            end
            push!(C, D)
        end
        a = [;;]
        return TensorField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("log(A::TensorField): TensorField type ($(A.type)) is not yet implemented.")
    end
end

###############################################################################
#                                                                             #
#                      System matrix operations                               #
#                                                                             #
###############################################################################

function *(A::Union{SystemMatrix,Matrix}, BB::VectorField)
    if BB.A != []
        B = elementsToNodes(BB)
    else
        B = BB
    end
    type = :none
    if B.a != [;;] && B.nsteps == 1
        if B.type == :v2D
            type = :v2D
        elseif B.type == :v3D
            type = :v3D
        elseif B.type == :other
            type = :other
        else
            error("*(A::Union{SystemMatrix,Matrix}, B::VectorField): neither 2D nor 3D (type=$(B.type))?")
        end
        C = A isa Matrix ? A * B.a : A.A * B.a
        return VectorField([], reshape(C, :, 1), [0.0], [], 1, type, B.model)
    else
        error("*(A::Union{SystemMatrix,Matrix}, B::VectorField): Type of data is not nodal or more than one time steps ($(B.nsteps)).")
    end
end

function *(A::Union{SystemMatrix,Matrix}, BB::ScalarField)
    if BB.A != []
        B = elementsToNodes(BB)
    else
        B = BB
    end
    type = :none
    if B.a != [;;] && B.nsteps == 1
        type = :scalar
        C = A isa Matrix ? A * B.a : A.A * B.a
        return ScalarField([], reshape(C, :, 1), [0.0], [], 1, type, B.model)
    else
        error("*(A::SystemMatrix, B::ScalarField): Type of data is not nodal or more than one time steps ($(B.nsteps)).")
    end
end

import Base.\
function \(A::Union{SystemMatrix,Matrix}, BB::ScalarField)
    if BB.A != []
        B = elementsToNodes(BB)
    else
        B = BB
    end
    type = :none
    if B.a != [;;] && B.nsteps == 1
        type = :scalar
        C = A isa Matrix ? A \ B.a : A.A \ B.a
        return ScalarField([], reshape(C, :, 1), [0.0], [], 1, type, B.model)
    else
        error("\\(A::SystemMatrix, B::ScalarField): Type of data is not nodal or more than one time steps ($(B.nsteps)).")
    end
end

function \(A::Union{SystemMatrix,Matrix}, BB::VectorField)
    if BB.A != []
        B = elementsToNodes(BB)
    else
        B = BB
    end
    type = :none
    if B.a != [;;] && B.nsteps == 1
        if B.type == :v2D
            type = :v2D
        elseif B.type == :v3D
            type = :v3D
        elseif B.type == :other
            type = :other
        else
            error("\\(A::SystemMatrix, B::VectorField): ")
        end
        C = A isa Matrix ? A \ B.a : A.A \ B.a
        return VectorField([], reshape(C, :, 1), [0.0], [], 1, type, B.model)
    else
        error("\\(A::SystemMatrix, B::VectorField): Type of data is not nodal or more than one time steps ($(B.nsteps)).")
    end
end

function \(A::Union{SystemMatrix,Matrix,SparseMatrixCSC}, b::SparseMatrixCSC)
    m, n = size(b)
    c = zeros(m)
    d = zeros(m, n)
    AA = A isa SystemMatrix ? lu(A.A) : lu(A)
    for i in 1:n
        c .= b[:, i]
        d[:, i] = AA \ c
    end
    return d
end

function ldiv_sparse!(X::SparseMatrixCSC, K::SystemMatrix, F::SparseMatrixCSC)
    n, m = size(K.A, 1), size(F, 2)
    Ffac = lu(K.A)
    x = zeros(n)
    b = zeros(n)
    I, V, J = Int[], Float64[], Int[]

    for j in 1:m
        idx, val = findnz(view(F, :, j))
        b .= 0
        b[idx] .= val
        ldiv!(x, Ffac, b)
        nz = findall(!iszero, x)
        append!(I, nz)
        append!(J, fill(j, length(nz)))
        append!(V, x[nz])
    end

    X = sparse(I, J, V, n, m)
    return X
end

function \(K::SystemMatrix, F::SparseMatrixCSC)
    X = spzeros(size(K.A, 1), size(F, 2))
    return ldiv_sparse!(X, K, F)
end

function *(A::SystemMatrix, b::Number)
    return SystemMatrix(A.A * b, A.model)
end

function *(b::Number, A::SystemMatrix)
    return SystemMatrix(A.A * b, A.model)
end

function +(A::SystemMatrix, B::SystemMatrix)
    return SystemMatrix(A.A + B.A, A.model)
end

function -(A::SystemMatrix, B::SystemMatrix)
    return SystemMatrix(A.A - B.A, A.model)
end

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
