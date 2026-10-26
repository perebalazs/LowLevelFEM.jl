# Test suite

Run from the package root with the existing project environment:

```sh
julia --project=. --threads=2 test/runtests.jl
```

`test/runtests.jl` remains the entry point for `Pkg.test()` too. The test
dependencies are already declared by the package; this suite adds none.

The runner always executes two separately named groups:

- **Modern architecture (DSL):** field data and operators, weak-form assembly,
  multifield systems, boundary conditions/MPC, `solveField`, `solveEigen`,
  contact and invalid inputs. FEM matrices and loads use the `∫` DSL.
- **Legacy:** old elasticity assembly/solves and compatibility comparisons.
  These files can eventually be removed together with their runner includes.

Each modern file owns one theme. `test_multifield.jl` covers block algebra and
coupled solving; `test_boundary_conditions.jl` owns the combined coordinate/MPC
regressions. Eigen constraints are also exercised through `solveEigen` because
modal reconstruction has a different contract from a static solve.

`helpers.jl` provides isolated Gmsh fixtures and numerical assertions. Models use
explicit small mesh sizes. Parallel checks use `threads=:auto` against
`threads=1`; a process with at least two Julia threads exercises both paths.
Single-threaded execution remains supported. Eigenvectors are checked by
residuals, normalization and constraints, without assuming signs or a particular
basis in a repeated eigenspace. ARPACK cases supply deterministic start vectors.

Independent references include exact P1 element matrices, polynomial static
solutions, constant-strain energies, the two-element eigenvalue λ=12, and flat
contact gaps, interface matrices and a finite-penalty compression solution.

The old plane-stress reference and tolerance, integration values 4/4.5/1 and
stored-column values 4/8, unit load resultants, solver/transformation regressions,
and development eigen spectra are retained. The old cube smoke model now uses
two elements per direction. Backward-compatible solver keywords and the
legacy/DSL stiffness comparison live only in `legacy/`.

Checks promoted from `examples/dev/eigen/` are adapted to field-valued default
output and explicit `output=:raw`. Eigen field copy/conversion/persistence checks
belong to `test_fields.jl`; invalid inputs belong to `test_validation.jl`.
Contact tests target `ContactGap` weak forms, not obsolete notebook members.

No transient analysis or time integrator is tested. Multiple static RHS columns
and stored-column selection exercise field representation only. Temporary
persistence files use `mktempdir`; tests do not write into the package sources.

## Operator regression conventions

`TensorDiv` uses the full nine-component, column-major `TensorField` storage in
both 2D and 3D. Nonsymmetric affine tensors test the component mapping through
CSC/IJV assembly and serial/automatic threading. The original reference
energies 5 and 155 are retained, with an additional 2D reference of 26.

`SymGrad` follows the documented engineering strain order
`[xx, yy, zz, xy, yz, zx]`, consistent with `VOIGT3D`. Unequal shear coefficients
give the reference energy 112. The three formerly broken expectations are
ordinary passing tests after the operator fixes in `src/multifield.jl`.

Custom 3D coefficient matrices prepared for the former `SymGrad` order
`[xx, yy, zz, xy, xz, yz]` need the corresponding row and column permutation:
`p = [1, 2, 3, 4, 6, 5]; D_new = D_old[p, p]`. Isotropic `D(:Solid, mat)` is
unchanged by this permutation.

## Tutorial notebook integration checks

The core suite checks small numerical models. Notebook execution is a separate,
opt-in integration check and does not increase normal `runtests.jl` runtime.
It uses the existing Julia environment and Python 3's standard library, without
IJulia, nbconvert, pip packages, or package installation. Run from the package root:

```sh
python3 test/run_tutorials.py --list
python3 test/run_tutorials.py
python3 test/run_tutorials.py buckling_and_prestressed_modes.ipynb
python3 test/run_tutorials.py --all --timeout 900
```

The catalog comes from notebook links in `docs/src/tutorials/`. The default run
covers the plane-stress and coupled-beam tutorials, including postprocessing.
It adds analytical plane-stress checks and beam equilibrium/constraint checks.
Other selected notebooks are checked for successful cell execution; additional
numerical checks can be registered in `tutorial_checks.jl`. `--all` can be
expensive, especially for 3D contact and curved high-order meshes. Transient
elasticity and space-time dynamics are explicitly excluded from this scope.

Each notebook runs in a fresh Julia process with notebook-style soft scope,
fixed RNG seed, copied input assets and a temporary working directory. The
original notebooks and geometry files are not edited. Cells execute in stored
order, with original cell numbers in the logs. Only explicit preprocessor,
postprocessor and `gmsh.fltk.run()` GUI-opening calls are replaced with logged
headless no-ops; view creation and numerical postprocessing still execute.
This does not verify the GUI's appearance or behavior.

Failures/timeouts produce a nonzero exit status. Logs, copied inputs and a JSON
report are retained in the printed temporary report directory, with excluded
and unselected notebooks listed explicitly. Missing optional packages or stale
notebook APIs are failures; the runner does not install packages or silently
skip numerical cells. Use `--julia /path/to/julia` to select an existing runtime.
