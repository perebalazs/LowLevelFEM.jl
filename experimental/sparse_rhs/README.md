# Sparse right-hand-side experiment

These files were moved unchanged from `src/`. They are standalone calibration
and benchmarking experiments, not part of the LowLevelFEM runtime or test suite.
The production linear backend is in `src/solve.jl`.

To try the experiment, use a separate Julia environment with `Polyester` and
`BenchmarkTools` installed, then include these files in order:

```julia
include("sparseRHS.jl")
include("sparseRHS-test.jl")
```

The calibration stores its selected threshold in
`~/.julia/umfpack_threshold.toml`. Experimental dependencies are intentionally
not added to LowLevelFEM's runtime dependencies.
