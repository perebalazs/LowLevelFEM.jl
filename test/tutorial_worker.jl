# One fresh Julia process per notebook. Python extracts the original code cells
# using only its standard library; Julia executes them with notebook soft scope.
using REPL, Random, Test

const _tutorial_workdir, _tutorial_name, _tutorial_checks = ARGS
const _tutorial_skipped_gui = String[]
Random.seed!(1234)

function _tutorial_skip_gui(name)
    push!(_tutorial_skipped_gui, name)
    println("HEADLESS_SKIP ", name)
    return nothing
end

function _tutorial_headless(ex)
    ex isa Expr || return ex
    if ex.head === :call
        name = string(ex.args[1])
        if name in ("openPreProcessor", "openPostProcessor",
                    "LowLevelFEM.openPreProcessor", "LowLevelFEM.openPostProcessor",
                    "gmsh.fltk.run", "LowLevelFEM.gmsh.fltk.run")
            return :(_tutorial_skip_gui($name))
        end
        if name in ("Pkg.add", "Pkg.develop", "Pkg.rm", "Pkg.update", "Pkg.resolve", "Pkg.instantiate")
            return :(error("Notebook package/environment changes are disabled: " * $name))
        end
    end
    return Expr(ex.head, map(_tutorial_headless, ex.args)...)
end

try
    cd(_tutorial_workdir) do
        cells = sort(filter(name -> startswith(name, "notebook_cell_") && endswith(name, ".jl"), readdir()))
        for cell in cells
            println("CELL ", _tutorial_name, " ", cell)
            flush(stdout)
            transform = ex -> REPL.softscope(_tutorial_headless(ex))
            Base.include_string(transform, Main, read(cell, String), joinpath(_tutorial_workdir, cell))
        end
        println("CELLS_COMPLETED ", length(cells))
        include(_tutorial_checks)
    end
finally
    if isdefined(Main, :LowLevelFEM)
        g = Main.LowLevelFEM.gmsh
        if g.isInitialized() != 0
            g.clear()
            g.finalize()
        end
    end
end
