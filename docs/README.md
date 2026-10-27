# Documentation maintenance

The existing `make.jl` remains the publishing entry point; its unused RSS
pipeline has been removed and its navigation constant is unchanged. Reference
topics distinguish the retained engineering API from direct DSL workflows and
legacy Poisson/nonlinear implementations. Engineering names such as
`stiffnessMatrix` and `solveDisplacement` are intended to remain as the backend
migrates to DSL-backed wrappers. Existing URLs remain available.
The detailed contact/operator explanations and notebook links are preserved.

## Public API docstrings

From the repository root, using the already configured Julia environment:

```sh
python3 docs/check_api_docs.py
python3 docs/check_api_docs.py --strict
```

The normal check fails for missing docstrings on defined exports and for new
undefined exports. It warns about four known undefined names (`getEigenVectors`,
`nodalForce!`, `nodalAcceleration!`, `probe_field`); strict mode fails for those
as well. No callable stubs are fabricated to hide the export inconsistencies.
The published API audit also maps all exports to their primary reference topics.
The external bundled `gmsh` API is not duplicated.

The check covers Julia's loaded bindings, including aliases and imported macros.
Presence alone is not a semantic completeness test: review new overloads,
keywords, defaults, field layouts, return values and optional extension behavior.
Base/LinearAlgebra field operations are also documented by method and summarized
under field algebra/calculus; not all of these are LowLevelFEM exports.

## Executable documentation examples

```sh
python3 docs/check_examples.py
python3 docs/check_examples.py --with-optional
```

Only self-contained `julia` fences whose first line is
`# docs-check: unique-name` are executed. A marker ending with `; optional`
requires `--with-optional` and already installed optional packages. Each runs in
a fresh Julia process, with two threads by default, in a temporary directory.
`--julia`, `--threads` and `--timeout` configure the run. No packages are installed.
GUI launch is omitted in checked examples; Gmsh views can still be constructed.

The checked examples cover 3D elasticity/postprocessing, the exact linear scalar
diffusion solution, an exact two-field mass system, a discrete eigenpair with
mass normalization, and optional anisotropic energy/tangent conventions using
Tensors. Existing illustrative fragments assume their documented surrounding
model and are not silently treated as standalone scripts.

The independent notebook runner under `test/run_tutorials.py` continues to test
selected complete tutorials. Documentation examples do not replace notebook or
package regression tests, mesh convergence studies, or validation of a new
physical formulation.

## Strict local HTML build without publishing

```sh
python3 docs/build_local.py
```

This uses the existing docs environment and installed dependencies. It reads only
`DOC_PAGES` from `make.jl`. A first build checks package-owned exported-doc
inclusion, and a second build uses the same unfiltered docstring selection as
`make.jl`, so inherited Base/LinearAlgebra cross-references are checked too.
It does not execute `Pkg.develop` or `deploydocs`. HTML is written
to a printed temporary path. `--output /absolute/new-or-empty-directory` selects
another location; a nonempty directory is rejected. Remote source/edit links are
disabled for this local preview.

The external `gmsh` module is excluded from Documenter's recursive inclusion
check. Package-owned exported docstrings must be included; the separate API
checker also catches undocumented bindings and undefined exports. Doctests
remain disabled in `make.jl`; marked numerical examples are checked separately.

All three helpers use Python's standard library and launch Julia with
`--compiled-modules=existing` to use existing caches without building new ones.
They do not edit project/dependency files. An existing installed Documenter is
required for the HTML build. The local helpers were validated with Julia 1.13.3.

## Scope and known findings

FEM implementation code is unchanged: `src/*.jl` edits are confined to
docstrings. In `docs/`, the unused RSS build blocks, generator, tracked backup
and generated feed have been removed; the News pages remain. The audit records the four undefined exports, the incomplete
`Adv(P)` convenience constructor, and the Tensors/DSL tangent-order difference.
The material guide documents and checks the two-axis tangent permutation.
Dedicated nonlinear and Poisson assembly remain explicitly marked as legacy;
shared field calculus and DSL material-response functions remain modern.

No dependencies, CI workflows, Codecov configuration, branches or commits are
changed. Publishing settings are unchanged apart from removing the unused RSS
steps. Existing transient documentation is retained without adding new transient
examples.
