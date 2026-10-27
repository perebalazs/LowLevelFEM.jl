# API Reference

```@meta
CurrentModule = LowLevelFEM
```

Find functions by the task you are performing. Use the problem-oriented
engineering API for standard workflows, or the weak-form DSL directly for
custom formulations. The engineering names remain part of the planned API;
their internal migration to DSL-backed wrappers is described below.

| Task | Reference | Practical guide |
|:--|:--|:--|
| Build a mesh, choose materials and unknowns | [Preprocessing](preprocessing.md) | [First problem](../getting-started/first-problem.md) |
| Construct, convert or save fields | [Fields and storage](fields.md) | [Nodal vs element fields](../explanations/nodal-vs-element-fields.md) |
| Combine field values and tensors | [Field algebra](operators.md) | [Operators and fields](../manual/operators-and-fields.md) |
| Differentiate a computed field | [Field calculus](field-calculus.md) | [Field calculus and DSL operators](../manual/operators-and-fields.md) |
| Assemble weak forms and mixed blocks | [Weak forms and assembly](multifield.md) | [Assembly and solvers](../manual/assembly-and-solvers.md) |
| Prescribe constraints and coordinate bases | [Boundary conditions and MPC](boundary-conditions.md) | [Boundary conditions and loads](../manual/boundary-and-loads.md) |
| Solve a linear system or eigenproblem | [Solution](solution.md) | [Matrix-level workflow](../explanations/matrix-level-workflow.md) |
| Evaluate an energy-based material | [Constitutive models](constitutive.md) | [Nonlinear formulations](../tutorials/nonlinear-large-deformation.md) |
| Model contact | [Contact](contact.md) | [Contact tutorials](../tutorials/contact-mechanics.md) |
| Inspect, integrate or display results | [Postprocessing](postprocessing.md) | [Postprocessing guide](../manual/postprocessing-and-visualization.md) |

## Engineering workflow and legacy implementations

- [API layers: engineering wrappers, DSL and legacy assembly](legacy.md)
- [Engineering linear mechanics and eigen helpers](linear.md)
- [Engineering heat assembly and solution](heat.md)
- [Legacy nonlinear assembly](nonlinear.md)
- [Legacy Poisson and dedicated differential matrices](poisson.md)

The sidebar retains its existing page names for link compatibility; this topic
map distinguishes retained engineering interfaces from legacy implementations. A function has one
primary docstring entry, with links from related topics instead of copied
reference text.

For a complete alphabetical inventory, aliases, and known API inconsistencies,
see the [public API audit](api-audit.md). The bundled `gmsh` module is re-exported;
its own API is documented by Gmsh.
