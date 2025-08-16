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

function log(A::ScalarField)
    if length(A.A) != 0
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
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("log(ScalarField): data at nodes is not yet implemented.")
    end
end

function sqrt(A::ScalarField)
    if length(A.A) != 0
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
    else
        error("sqrt(ScalarField): data at nodes is not yet implemented.")
    end
end

###############################################################################
#                                                                             #
#                      Vector fields operations                               #
#                                                                             #
###############################################################################

import LinearAlgebra.⋅
import Base.∘
import LinearAlgebra.norm
import LinearAlgebra.diagm

function *(A::ScalarField, B::VectorField)
    if length(A.A) != 0 && length(B.A) != 0
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
            n = length(B.A[i]) ÷ 3
            if n != sz
                D = zeros(3n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    D[3j-2:3j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][3j-2:3j, k]
                end
            end
            append!(num, sec[i])
            push!(C, D)
        end
        a = [;;]
        return VectorField(C, a, B.t, num, B.nsteps, :vector, B.model)
    else
        error("*(ScalarField, VectorField): data at nodes is not yet implemented.")
    end
end

function *(B::VectorField, A::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
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
            n = length(B.A[i]) ÷ 3
            if n != sz
                D = zeros(3n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    D[3j-2:3j, k] = A.A[indS[i]][j, k] * B.A[indT[i]][3j-2:3j, k]
                end
            end
            append!(num, sec[i])
            push!(C, D)
        end
        a = [;;]
        return VectorField(C, a, B.t, num, B.nsteps, :vector, B.model)
    else
        error("*(VectorField, ScalarField): data at nodes is not yet implemented.")
    end
end

function /(B::VectorField, A::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
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
            n = length(B.A[i]) ÷ 3
            if n != sz
                D = zeros(3n, nsteps)
                sz = n
            end
            for j in 1:n
                for k in 1:nsteps
                    D[3j-2:3j, k] = B.A[indT[i]][3j-2:3j, k] / A.A[indS[i]][j, k]
                end
            end
            append!(num, sec[i])
            push!(C, D)
        end
        a = [;;]
        return VectorField(C, a, B.t, num, B.nsteps, :vector, B.model)
    else
        error("/(VectorField, ScalarField): data at nodes is not yet implemented.")
    end
end

function +(A::VectorField, B::VectorField)
    if length(A.A) != 0 && length(B.A) != 0
        if (A.type == :u3D && B.type == :u3D) || (A.type == :u2D && B.type == :u2D) || (A.type == :f3D && B.type == :f3D) || (A.type == :f2D && B.type == :f2D)
            if length(A.A) != length(B.A)
                error("+(A::VectorField, B::VectorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
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
            for i in eachindex(sec)
                #n = length(A.A[i]) ÷ 9
                #m = length(B.A[i]) ÷ 9
                #if n != m
                #    error("+(A::VectorField, B::VectorField): size of A.A[$i]=$(9n) != size of B.A[$j]=$(9m)")
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
            return VectorField(C, a, A.t, num, A.nsteps, A.type, A.model)
        else
            error("+(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if (A.type == :u3D && B.type == :u3D) || (A.type == :u2D && B.type == :u2D) || (A.type == :f3D && B.type == :f3D) || (A.type == :f2D && B.type == :f2D)
            return VectorField([], A.a + B.a, A.t, [], A.nsteps, A.type, A.model)
        else
            error("+(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("+(VectorField, VectorField): internal error")
    end
end

function -(A::VectorField, B::VectorField)
    if length(A.A) != 0 && length(B.A) != 0
        #if (A.type == :u3D && B.type == :u3D) || (A.type == :u2D && B.type == :u2D) || (A.type == :f3D && B.type == :f3D) || (A.type == :f2D && B.type == :f2D)
        if A.type == B.type
            if length(A.A) != length(B.A)
                error("-(A::VectorField, B::VectorField): size of A=$(length(A.A)) != size of B=$(length(B.A))")
            end
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
            return VectorField(C, a, A.t, num, A.nsteps, A.type, A.model)
        else
            error("-(A::VectorField, B::VectorField): Operation with type ($(A.type) and $(B.type)) is not supported.")
        end
    elseif length(A.a) != 0 && length(B.a) != 0
        if (A.type == :u3D && B.type == :u3D) || (A.type == :u2D && B.type == :u2D) || (A.type == :f3D && B.type == :f3D) || (A.type == :f2D && B.type == :f2D)
            return VectorField([], A.a - B.a, A.t, [], A.nsteps, A.type, A.model)
        else
            error("-(A::VectorField, B::VectorField): VectorField type ($(A.type) and $(B.type)) is not yet implemented.")
        end
    else
        error("-(VectorField, VectorField): internal error")
    end
end

function *(A::VectorField, b::Number)
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
        if A.type == :u3D || A.type == :u2D || A.type == :f3D || A.type == :f2D
            nsteps = A.nsteps
            C = []
            for i in 1:length(A.A)
                D = A.A[i] / b
                push!(C, D)
            end
            a = [;;]
            return VectorField(C, a, A.t, A.numElem, A.nsteps, A.type, A.type)
        else
            error("/(A::VectorField, b): VectorField type ($(A.type) or $(B.type)) is not yet implemented.")
        end
    elseif length(A.a) != 0
        return VectorField([], A.a / b, A.t, [], A.nsteps, A.type, A.model)
    else
        error("/(VectorField, b): internal error")
    end
end

function ⋅(A::VectorField, B::VectorField)
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
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("*(A::VectorField, B::VectorField): VectorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function *(A::VectorField, B::VectorField)
    return A ⋅ B
end

function ∘(a::VectorField, b::VectorField)
    G = MMatrix{3,3}([0.0 0 0;0 0 0;0 0 0])
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
                            G[p, q] = e[p,1] * f[q,1]
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

function norm(A::VectorField)
    if A.A != [] #(A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
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
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("norm(A::VectorField): data at nodes is not yet implemented.")
    end
end

function diagm(A::VectorField)
    if length(A.A) != 0
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
    else
        error("diagm(VectorField): data at nodes is not yet implemented.")
    end
end

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

function *(A::TensorField, B::TensorField)
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

function ⋅(A::TensorField, B::TensorField)
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
        return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
    else
        error("*(A::TensorField, B::TensorField): TensorField type ($(A.type) or $(B.type)) is not yet implemented.")
    end
end

function +(A::TensorField, B::TensorField)
    if length(A.A) != 0 && length(B.A) != 0
        if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
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
    elseif A.a != [;;] && B.a != [;;]
        return TensorField([], A.a + B.a, A.t, [], A.steps, A.type, A.model)
    else
        error("+(TensorField, TensorField): internal error.")
    end
end

function -(A::TensorField, B::TensorField)
    if length(A.A) != 0 && length(B.A) != 0
        if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
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
    elseif A.a != [;;] && B.a != [;;]
        return TensorField([], A.a - B.a, A.t, [], A.nsteps, A.type, A.model)
    else
        error("-(TensorField, TensorField): internal error.")
    end
end

function *(A::TensorField, b::VectorField)
    #if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
    if b.model.dim == 3
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

function *(b::VectorField, A::TensorField)
    #if (A.type == :s || A.type == :e || A.type == :F) && (B.type == :s || B.type == :e || B.type == :F)
    if b.model.dim == 3
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

function *(A::TensorField, b::Number)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
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
    elseif A.a != [;;]
        return TensorField([], A.a * b, A.t, [], A.nsteps, A.type, A.model)
    else
        error("*(TensorField, Any): internal error.")
    end
end

function *(b::Number, A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
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
    elseif A.a != [;;]
        return TensorField([], A.a * b, A.t, [], A.nsteps, A.type, A.model)
    else
        error("*(Any, TensorField): internal error.")
    end
end

function /(A::TensorField, b::Number)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
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
    else
        error("/(TensorField, Any): data at nodes is not yet implemented.")
    end
end

function *(A::ScalarField, B::TensorField)
    if length(A.A) != 0 && length(B.A) != 0
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
    else
        error("*(ScalarField, TensorField): data at nodes is not yet implemented.")
    end
end

function *(B::TensorField, A::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
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
    else
        error("*(TensorField, ScalarField): data at nodes is not yet implemented.")
    end
end

function /(B::TensorField, A::ScalarField)
    if length(A.A) != 0 && length(B.A) != 0
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
    else
        error("/(TensorField, ScalarField): data at nodes is not yet implemented.")
    end
end

function transpose(A::TensorField)
    if length(A.A) != 0
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
    else
        error("transpose(TensorField): data at nodes is not yet implemented.")
    end
end

function adjoint(A::TensorField)
    if length(A.A) != 0
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
    else
        error("adjoint(TensorField): data at nodes is not yet implemented.")
    end
end

export unitTensor
function unitTensor(A::TensorField)
    if length(A.A) != 0
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
    else
        error("transpose(TensorField): data at nodes is not yet implemented.")
    end
end

export trace
function trace(A::TensorField)
    if length(A.A) != 0
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
            return ScalarField(C, a, A.t, A.numElem, A.nsteps, :e, A.model)
        else
            error("trace(A::TensorField): TensorField type ($(A.type) is not yet implemented.")
        end
    else
        error("trace(TensorField): data at nodes is not yet implemented.")
    end
end

function det(A::TensorField)
    if length(A.A) != 0
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
            return ScalarField(C, a, A.t, A.numElem, A.nsteps, :sc, A.model)
        else
            error("det(A::TensorField): TensorField type ($(A.type) is not yet implemented.")
        end
    else
        error("det(TensorField): data at nodes is not yet implemented.")
    end
end

function inv(A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
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
    else
        error("inv(TensorField): data at nodes is not yet implemented.")
    end
end

function eigen(A::TensorField)
    if length(A.A) != 0
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
                        GG = mapslices(v -> v ./ norm(v), G, dims=1)
                        D[9j-8:9j, k] = reshape(GG, 9, 1)
                        h[3j-2:3j, k] = reshape(f, 3, 1)
                    end
                end
                push!(C, D)
                push!(c, h)
            end
            a = [;;]
            return VectorField(c, a, A.t, A.numElem, A.nsteps, :e3D, A.model), TensorField(C, a, A.t, A.numElem, A.nsteps, :Q, A.model)
        else
            error("eigen(A::TensorField): A is not a TensorField.")
        end
    else
        error("eigen(TensorField): data at nodes is not yet implemented.")
    end
end

function sqrt(A::TensorField)
    if length(A.A) != 0
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
    else
        error("sqrt(TensorField): data at nodes is not yet implemented.")
    end
end

function log(A::TensorField)
    if length(A.A) != 0
        if A.type == :s || A.type == :e || A.type == :F
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
    else
        error("log(TensorField): data at nodes is not yet implemented.")
    end
end

###############################################################################
#                                                                             #
#                      System matrix operations                               #
#                                                                             #
###############################################################################

function *(A::Union{SystemMatrix,Matrix}, B::VectorField)
    type = :none
    if B.a != [;;] && B.nsteps == 1
        if B.type == :u2D
            type = :f2D
        elseif B.type == :u3D
            type = :f3D
        elseif B.type == :other
            type = :other
        else
            error("*(A::Union{SystemMatrix,Matrix}, B::VectorField): neither 2D nor 3D (type=$(B.type))?")
        end
        C = A isa Matrix ? A * B.a : A.A * B.a
        return VectorField([], reshape(C, :,1), [0.0], [], 1, type, B.model)
    else
        error("*(A::Union{SystemMatrix,Matrix}, B::VectorField): Type of data is not nodal or more than one time steps ($(B.nsteps)).")
    end
end

function *(A::Union{SystemMatrix,Matrix}, B::ScalarField)
    type = :none
    if B.a != [;;] && B.nsteps == 1
        if B.type == :T
            type = :qn
        elseif B.type == :other
            type = :other
        else
            error("*(A::SystemMatrix, B::ScalarField): ")
        end
        C = A isa Matrix ? A * B.a : A.A * B.a
        return ScalarField([], reshape(C, :,1), [0.0], [], 1, type, B.model)
    else
        error("*(A::SystemMatrix, B::ScalarField): Type of data is not nodal or more than one time steps ($(B.nsteps)).")
    end
end

import Base.\
function \(A::Union{SystemMatrix,Matrix}, B::ScalarField)
    type = :none
    if B.a != [;;] && B.nsteps == 1
        if B.type == :qn
            type = :T
        elseif B.type == :other
            type = :other
        else
            error("\\(A::SystemMatrix, B::ScalarField): type = $(B.type)")
        end
        C = A isa Matrix ? A \ B.a : A.A \ B.a
        return ScalarField([], reshape(C, :,1), [0.0], [], 1, type, B.model)
    else
        error("\\(A::SystemMatrix, B::ScalarField): Type of data is not nodal or more than one time steps ($(B.nsteps)).")
    end
end

function \(A::Union{SystemMatrix,Matrix}, B::VectorField)
    type = :none
    if B.a != [;;] && B.nsteps == 1
        if B.type == :f2D
            type = :u2D
        elseif B.type == :f3D
            type = :u3D
        elseif B.type == :other
            type = :other
        else
            error("\\(A::SystemMatrix, B::VectorField): ")
        end
        C = A isa Matrix ? A \ B.a : A.A \ B.a
        return VectorField([], reshape(C, :,1), [0.0], [], 1, type, B.model)
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
        c .= b[:,i]
        d[:,i] = AA \ c
    end
    return d
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

function *(A::Transformation, B::VectorField)
    n = size(B.a, 1)
    non = A.non
    dim = A.dim
    if dim * non == n
        if B.A == []
            v =  A.T * B.a
            return VectorField([], v, B.t, [], length(B.t), B.type, B.model)
        else
            error("*(A::Transformation, B::VectorField): B contains element data instead of nodal data.")
        end
    else
        error("*(A::Transformation, B::VectorField): size missmatch A.dim * A.non = $dim * $non ≠ $n = size(B.a, 1).")
    end
end

function *(B::VectorField, A::Transformation)
    n = size(B.a, 1)
    non = A.non
    dim = A.dim
    if dim * non == n
        if B.A == []
            v =  (B.a' * A.T)'
            return VectorField([], v, B.t, [], length(B.t), B.type, B.model)
        else
            error("*(B::VectorField, A::Transformation): B contains element data instead of nodal data.")
        end
    else
        error("*(B::VectorField, A::Transformation): size missmatch A.dim * A.non = $dim * $non ≠ $n = size(B.a, 1).")
    end
end

function *(A::Transformation, B::TensorField)
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

function *(B::TensorField, A::Transformation)
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
