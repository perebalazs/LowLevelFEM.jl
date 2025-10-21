module DebugTools

export @showfields, @showstruct, @showdef, @showtype, @showmem, @showmethods, @disp, @showsize

"Listázza a struct mezőit és azok aktuális értékét."
macro showfields(expr)
    return quote
        x = $(esc(expr))
        T = typeof(x)
        println("$(T):")
        for n in fieldnames(T)
            println("  ", n, " = ", getfield(x, n))
        end
        nothing
    end
end

"Listázza a struct mezőit típusokkal és aktuális értékekkel."
macro showstruct(expr)
    return quote
        x = $(esc(expr))
        T = typeof(x)
        println("$(T):")
        for (n, t) in zip(fieldnames(T), fieldtypes(T))
            println("  ", n, "::", t, " = ", getfield(x, n))
        end
        nothing
    end
end

"Csak a struct definícióját (mezőnevek és típusok) írja ki, értékek nélkül."
macro showdef(expr)
    return quote
        T = typeof($(esc(expr)))
        println("struct $(T):")
        for (n, t) in zip(fieldnames(T), fieldtypes(T))
            println("  ", n, "::", t)
        end
        nothing
    end
end

"Kiírja az adott kifejezés típusát és méretét (ha ismert)."
macro showtype(expr)
    return quote
        x = $(esc(expr))
        T = typeof(x)
        println("Type: ", T)
        try
            println("Size: ", Base.summarysize(x), " bytes")
        catch
            # summarysize may fail for some types
        end
        nothing
    end
end

"Megméri a kifejezés futási idejét és memória-allokációját (BenchmarkTools nélkül)."
macro showmem(expr)
    return quote
        local result = @timed $(esc(expr))
        println("Time: ", round(result.time * 1000, digits=3), " ms")
        println("Allocations: ", result.bytes, " bytes")
        result.value
    end
end

"Kiírja a megadott függvény összes metódusát és azok aláírásait."
macro showmethods(expr)
    return quote
        f = $(esc(expr))
        println("Methods for $(f):")
        for m in methods(f)
            println("  ", m)
        end
        nothing
    end
end

macro disp(expr)
    s = string(expr)
    return :(display($s * " = " * string($(esc(expr)))))
end

macro showsize(expr)
    return quote
        local x = $(esc(expr))
        local s = Base.summarysize(x)
        local units = ["B", "kB", "MB", "GB"]
        local val = s
        local i = 1
        while val > 1024 && i < length(units)
            val /= 1024
            i += 1
        end
        println("Size of ", $(string(expr)), ": ",
                round(val, digits=3), " ", units[i],
                "  (", s, " bytes)")
        nothing
    end
end

end # module DebugTools

