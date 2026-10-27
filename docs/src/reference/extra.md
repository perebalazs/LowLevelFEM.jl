# Projection Utilities and Diagnostics

```@meta
CurrentModule = LowLevelFEM
```

Specialized field projection and volume operations.

```@docs
projectScalarField
fieldsToVolume
```

## Interactive diagnostics

```@docs
LowLevelFEM.@showfields
LowLevelFEM.@showstruct
LowLevelFEM.@showdef
LowLevelFEM.@showtype
LowLevelFEM.@showmem
LowLevelFEM.@showmethods
LowLevelFEM.@disp
LowLevelFEM.@showsize
```

## Bundled Gmsh module

The re-exported module provides access to the bundled Gmsh API. The full
external API is maintained in the Gmsh documentation, rather than duplicated
in this package's reference.

See the [Gmsh reference manual](https://gmsh.info/doc/texinfo/) for the external API.
