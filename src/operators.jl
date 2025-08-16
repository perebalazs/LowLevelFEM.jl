###############################################################################
#                                                                             #
#                      Scalar fields operations                               #
#                                                                             #
###############################################################################

import Base.*
import Base./
import Base.+
import Base.-

function *(A::ScalarField, B::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
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
            if n != sz
                D = zeros(n, nsteps)
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
        return ScalarField(C, a, A.t, num, A.nsteps, :e, A.model)
    else
        error("*(ScalarField, ScalarField): data at nodes is not yet implemented.")
    end
end

function /(A::ScalarField, B::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
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
            if n != sz
                D = zeros(n, nsteps)
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
        return ScalarField(C, a, A.t, num, A.nsteps, :e, A.model)
    else
        error("/(ScalarField, ScalarField): data at nodes is not yet implemented.")
    end
end

function +(A::ScalarField, B::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
        if A.type == B.type
            if length(A.A) != length(B.A)
                error("+(A::ScalarField, B::ScalarField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
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
            for i in eachindex(sec)
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

import Base.-
function -(A::ScalarField, B::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
        if A.type == B.type
            if length(A.A) != length(B.A)
                error("-(A::ScalarField, B::ScalarField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
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
            return ScalarField(C, a, A.t, num, A.nsteps, A.type, A.model)
        else
            error("+(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if A.type == B.type
            return ScalarField([], A.a - B.a, A.t, [], A.nsteps, A.type, A.model)
        else
            error("+(A::ScalarField, B::ScalarField): ScalarField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("+(ScalarField, ScalarField): internal error")
    end
end

function *(A::ScalarField, b::Number)
    if length(A.A) != 0
        C = []
        for i in 1:length(A.A)
            D = A.A[i] * b
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :sc, A.model)
    else
        error("*(ScalarField, Any): data at nodes is not yet implemented.")
    end
end

function *(b::Number, A::ScalarField)
    if length(A.A) != 0
        C = []
        for i in 1:length(A.A)
            D = A.A[i] * b
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :sc, A.model)
    else
        error("*(Any, ScalarField): data at nodes is not yet implemented.")
    end
end

function /(A::ScalarField, b::Number)
    if length(A.A) != 0
        C = []
        for i in 1:length(A.A)
            D = A.A[i] / b
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :sc, A.model)
    else
        error("/(ScalarField, Any): data at nodes is not yet implemented.")
    end
end

function /(b::Number, A::ScalarField)
    if length(A.A) != 0
        C = []
        for i in 1:length(A.A)
            D = b ./ A.A[i]
            push!(C, D)
        end
        a = [;;]
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :sc, A.model)
    else
        error("/(Any, ScalarField): data at nodes is not yet implemented.")
    end
end
