# Source layout

All production definitions belong to the `LowLevelFEM` module. The entry files
`general.jl` and `multifield.jl` load their parts in the original evaluation
order; the subdirectories do not introduce nested modules or change the API.

| Directory | Contents |
| --- | --- |
| `general/` | Materials and problems, system matrices, fields, boundary conditions, field calculus and builders, mapping, integration, postprocessing, I/O and mesh helpers |
| `multifield/` | Differential operators, coefficient preparation, bilinear and linear assembly, reductions, basic and compound weak-form expressions, block boundary data and constitutive matrices |

The other physics and solver entry files remain alongside these loaders.
Optional package integrations live in `ext/`; constitutive tensor conversion
uses the common `VOIGT3D` order `[xx, yy, zz, xy, yz, zx]` in `nonlinear.jl`.

The standalone sparse-RHS experiment lives in `experimental/sparse_rhs/` and
is not loaded by the package or the test suite.
