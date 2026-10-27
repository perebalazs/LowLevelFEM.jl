#!/usr/bin/env python3
"""Build and check the docs locally, without running Pkg operations or deployment."""
import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

BUILD = r'''
using Documenter
root, output = ARGS
pushfirst!(LOAD_PATH,root)
using LowLevelFEM

# Read only the existing navigation constant; do not execute docs/make.jl.
parsed = Meta.parseall(read(joinpath(root,"docs","make.jl"),String))
nav = only(filter(parsed.args) do expr
    expr isa Expr && expr.head === :const &&
    expr.args[1] isa Expr && expr.args[1].head === :(=) &&
    expr.args[1].args[1] === :DOC_PAGES
end)
pages = Core.eval(Main,nav.args[1].args[2])
common = (
    root=joinpath(root,"docs"),
    sitename="LowLevelFEM",
    pages=pages,
    remotes=nothing,
    doctest=false,
)
html() = Documenter.HTML(prettyurls=false, edit_link=nothing,
    repolink="https://github.com/perebalazs/LowLevelFEM.jl",
    size_threshold_warn=300_000, size_threshold=600_000,
    assets=["assets/plausible.js","assets/custom.css"])

println("Checking package-owned exported docstring inclusion...")
mktempdir() do strict_output
    makedocs(; common..., build=strict_output, format=html(),
        modules=[LowLevelFEM], checkdocs=:exports,
        checkdocs_ignored_modules=[LowLevelFEM.gmsh])
end

# docs/make.jl has no modules filter: include the same external method docs so
# references inherited from Base/LinearAlgebra are checked in this pass too.
println("Building with the unfiltered docstrings used by docs/make.jl...")
makedocs(; common..., build=output, format=html())
println("Built documentation: ",joinpath(output,"index.html"))
'''

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--julia',default='julia')
    parser.add_argument('--output',type=Path,help='new or empty directory; default: temporary build directory')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    if args.output:
        output=args.output.resolve()
        if output.exists() and (not output.is_dir() or any(output.iterdir())):
            parser.error('--output must be a new or empty directory')
    else:
        output=Path(tempfile.mkdtemp(prefix='llfem-doc-build-'))/'html'
    with tempfile.TemporaryDirectory(prefix='llfem-doc-builder-') as folder:
        script=Path(folder)/'build.jl'
        script.write_text(BUILD)
        result=subprocess.run([args.julia,'--startup-file=no','--compiled-modules=existing',
            '--project='+str(root/'docs'),str(script),str(root),str(output)],cwd=root,check=False)
    return result.returncode

if __name__=='__main__':
    sys.exit(main())
