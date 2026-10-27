module DebugTools

export @showfields, @showstruct, @showdef, @showtype, @showmem, @showmethods, @disp, @showsize

"""
    @showfields expr

Evaluate `expr` once and print its type, field names and current field values.
Return `nothing`; the value must support structural field access.
"""
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

"""
    @showstruct expr

Evaluate `expr` once and print its type plus field names, declared types and
current values. Return `nothing`.
"""
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

"""
    @showdef expr

Evaluate `expr` once and print a structural summary of its type: field names
and declared field types, without current values. Return `nothing`.
"""
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

"""
    @showtype expr

Evaluate `expr` once and print its type and total reachable object size from
`Base.summarysize` when available. Return `nothing`; this is not array shape.
"""
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

"""
    @showmem expr

Evaluate `expr` with `@timed`, print elapsed milliseconds and allocated bytes,
and return the expression result. First-call compilation is included; this
single sample is not a statistical benchmark.
"""
macro showmem(expr)
    return quote
        local result = @timed $(esc(expr))
        println("Time: ", round(result.time * 1000, digits=3), " ms")
        println("Allocations: ", result.bytes, " bytes")
        result.value
    end
end

"""
    @showmethods expr

Evaluate a callable expression once and print its available methods and
signatures. Return `nothing`.
"""
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

"""
    @disp expr

Evaluate `expr` once and display its source text followed by its string value.
The macro returns the result of `display` (normally `nothing`).
"""
macro disp(expr)
    s = string(expr)
    return :(display($s * " = " * string($(esc(expr)))))
end

"""
    @showsize expr

Evaluate `expr` once and print its total reachable memory from
`Base.summarysize`, using B/kB/MB/GB with factors of 1024. Return `nothing`.
This reports object storage, not allocations during evaluation.
"""
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

