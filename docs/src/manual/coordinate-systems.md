# Coordinate Systems

LowLevelFEM supports local coordinate systems for nodal degrees of freedom. Coordinate transformations can be defined globally or by user-defined functions, making it possible to model curvilinear or spatially varying local axes.

Typical applications include:

- beam and shell elements with local axes
- cylindrical and curvilinear coordinate systems
- boundary conditions defined in local coordinates
- loads acting in local directions

Coordinate systems can be assigned before assembly and are automatically taken into account during the solution process.

For implementation details and available functions, see the API reference below.

## Related API

- [Reference: Preprocessing](../reference/preprocessing.md)
