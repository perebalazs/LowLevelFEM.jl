# Installation

Install the package from the Julia package registry:

```julia
using Pkg
Pkg.add("LowLevelFEM")
```

Load the package:

```julia
using LowLevelFEM
```

LowLevelFEM re-exports `gmsh` via `gmsh_jll`, so no separate `Gmsh.jl` package is required for standard workflows.

## Verify Setup

```julia
using LowLevelFEM

gmsh.initialize()
gmsh.finalize()
```

If this runs without error, your setup is ready.
