# Performance and Parallel Assembly

## Warm-up and timing

Separate package loading, first-call compilation and repeated execution when
measuring a model. Keep mesh, interpolation order, quadrature, optional packages
and thread count identical between comparisons. The package precompile workload
covers representative DSL paths; it cannot cover every user-defined coefficient
or specialization. See `test/precompile_latency.jl` for separate first/repeated
call measurements.

## Threads

Start Julia with the desired worker count, for example `julia --threads=4`.
Bilinear and linear DSL assembly defaults to `threads=:auto`; set `threads=1`
for a serial comparison or pass a positive worker count for a particular integral.
Small meshes may be faster serially because scheduling and reduction have a cost.
Do not mutate fields or the Gmsh model concurrently with assembly.

```julia
K1 = ∫(Grad(P) ⋅ Grad(P); threads=1)
Kp = ∫(Grad(P) ⋅ Grad(P); threads=:auto)
```

## Repeated assembly and memory

Direct CSC assembly is the default. Reuse a compatible structural pattern when
the mesh, integration domain and trial/test spaces are unchanged:

```julia
pattern = build_csc_pattern(P, P; Ω="body")
fill!(pattern.nzval, 0.0)
K = ∫(Grad(P) ⋅ Grad(P); Ω="body", csc_matrix=pattern)
```

Assembly adds to the supplied values. Reset before an independent assembly,
copy results that must survive reuse, and never remove structural zeros with
`dropzeros!`. Parallel CSC assembly uses private value buffers per worker,
so worker count also affects memory use. `element_chunk_size` controls work
partitioning on supported CSC bilinear/compound paths. `assembly=:ijv` retains
the triplet route for comparison; keyword support varies by expression method.
See [assembly details](../reference/multifield.md).

## Solver cost

Use a solver appropriate to the matrix and constraints. For mixed indefinite
systems, an unsuitable CG/Cholesky choice cannot be repaired by adding threads.
Optional iterative backends and preconditioners are described under
[solveField](../reference/solution.md). Changing quadrature or approximation
space changes the numerical model and must be validated separately.
