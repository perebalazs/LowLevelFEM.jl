# Weak-Form DSL Design

LowLevelFEM provides a DSL for expressing weak forms directly using operator objects.

## Core Idea

Expressions such as

```julia
K = ∫(Grad(Pu) ⋅ Grad(Pu); Ω="domain")
```

build finite element operators explicitly from composable pieces.

## Status

Placeholder page for a detailed design walkthrough.
