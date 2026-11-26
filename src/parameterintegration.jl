# Integration over an independent scalar parameter (e.g. shell thickness).
const _ParameterField = Union{ScalarField,VectorField,TensorField}

# An n-point Gauss-Legendre rule integrates polynomials through degree 2n-1.
function _parameter_quadrature(order)
    order isa Integer && !(order isa Bool) && 1 <= order <= typemax(Int) ||
        throw(ArgumentError("Parameter integration: order must be a positive integer."))
    n = Int(order) ÷ 2 + 1
    # Golub-Welsch: eigenvalues of the Legendre Jacobi matrix, with the
    # squared first components of the normalized eigenvectors as weights.
    offdiag = [k / sqrt(4.0*k*k - 1) for k in 1:n-1]
    rule = eigen(SymTridiagonal(zeros(n), offdiag))
    return rule.values, 2 .* abs2.(rule.vectors[1, :])
end

function _parameter_element_field(f::_ParameterField, group)
    g = nodesToElements(f; onPhysicalGroup=group)
    isElementwise(g) || throw(ArgumentError("Parameter integration: field has no element values."))
    length(g.numElem) == length(g.A) == length(unique(g.numElem)) ||
        throw(ArgumentError("Parameter integration: invalid or duplicate element tags."))
    length(g.t) == g.nsteps && g.nsteps >= 1 && all(A -> size(A, 2) == g.nsteps, g.A) ||
        throw(ArgumentError("Parameter integration: inconsistent field time columns."))
    isempty(group) && return g

    # Restrict by element tag. Do not extrapolate or average values across
    # adjacent elements, even if the field covers only part of this group.
    gmsh.model.setCurrent(g.model.name)
    entities = gmsh.model.getEntitiesForPhysicalName(group)
    isempty(entities) && throw(ArgumentError("Parameter integration: physical group '$group' not found."))
    tags = Int[]
    for (dim, tag) in entities
        _, elements, _ = gmsh.model.mesh.getElements(dim, tag)
        for block in elements
            append!(tags, block)
        end
    end
    unique!(tags)
    isempty(tags) && throw(ArgumentError("Parameter integration: physical group '$group' has no elements."))
    index = Dict(e => i for (i, e) in enumerate(g.numElem))
    selected = map(tags) do e
        haskey(index, e) || throw(ArgumentError("Parameter integration: field does not cover physical group '$group' (missing element $e)."))
        index[e]
    end
    return typeof(g)(g.A[selected], [;;], g.t, tags, g.nsteps, g.type, g.model, g.kind)
end

# Only the accumulator owns new arrays. Input field blocks are read-only.
_parameter_weighted(x::Number, w, group) = w*x
function _parameter_weighted(x::_ParameterField, w, group)
    f = _parameter_element_field(x, group)
    return typeof(f)([w .* A for A in f.A], [;;], copy(f.t), copy(f.numElem),
                     f.nsteps, f.type, f.model, f.kind)
end
function _parameter_array_check(x)
    all(v -> v isa Union{Number,ScalarField}, x) ||
        throw(ArgumentError("Parameter integration: array entries must be numbers or ScalarFields."))
end
function _parameter_weighted(x::Union{AbstractVector,AbstractMatrix}, w, group)
    _parameter_array_check(x)
    return map(v -> _parameter_weighted(v, w, group), x)
end
_parameter_weighted(x, w, group) =
    throw(ArgumentError("Parameter integration: unsupported integrand result $(typeof(x))."))

_parameter_accumulate(a::Number, b::Number, w, group) = a + w*b

function _parameter_add_constant!(a::ScalarField, b::Number)
    for A in a.A
        A .+= b
    end
    return a
end
_parameter_accumulate(a::ScalarField, b::Number, w, group) =
    _parameter_add_constant!(a, w*b)
_parameter_accumulate(a::Number, b::ScalarField, w, group) =
    _parameter_add_constant!(_parameter_weighted(b, w, group), a)

function _parameter_accumulate(a::_ParameterField, b::_ParameterField, w, group)
    f = _parameter_element_field(b, group)
    typeof(a) === typeof(f) && a.model === f.model && a.type === f.type && a.kind === f.kind ||
        throw(ArgumentError("Parameter integration: incompatible field types or Problems."))
    a.nsteps == f.nsteps && a.t == f.t ||
        throw(ArgumentError("Parameter integration: field time grids must agree."))
    length(a.numElem) == length(f.numElem) ||
        throw(ArgumentError("Parameter integration: field element supports must agree."))
    # Element ordering may differ, but the support and local layouts must agree.
    index = a.numElem == f.numElem ? nothing : Dict(e => i for (i, e) in enumerate(f.numElem))
    for (i, e) in enumerate(a.numElem)
        j = index === nothing ? i : get(index, e, 0)
        j != 0 || throw(ArgumentError("Parameter integration: field element supports must agree (missing element $e)."))
        size(a.A[i]) == size(f.A[j]) ||
            throw(ArgumentError("Parameter integration: incompatible local field sizes on element $e."))
        a.A[i] .+= w .* f.A[j]
    end
    return a
end

function _parameter_accumulate(a::Union{AbstractVector,AbstractMatrix},
                               b::Union{AbstractVector,AbstractMatrix}, w, group)
    axes(a) == axes(b) || throw(ArgumentError("Parameter integration: array dimensions must agree."))
    _parameter_array_check(b)
    return map((x, y) -> _parameter_accumulate(x, y, w, group), a, b)
end
_parameter_accumulate(a, b, w, group) =
    throw(ArgumentError("Parameter integration: incompatible integrand results $(typeof(a)) and $(typeof(b))."))

"""
    integrate(f, a::Real, b::Real; order=5, onPhysicalGroup="")
    ∫(f, a::Real, b::Real; order=5, onPhysicalGroup="")

Integrate a callable `f(z)` over an independent scalar parameter `z` from
`a` to `b`, using Gauss-Legendre quadrature. Spatial coordinates and stored
field time columns are unchanged. Typical applications are thickness or
cross-section integrals of material or stress fields.

`order` is a positive integer specifying the polynomial degree to integrate
exactly (up to floating-point error). The rule uses `order ÷ 2 + 1` points.
The default `order=5` uses three points. Use separate integrals for layers or
other discontinuities of `f`. Reversed bounds give a signed integral; equal
bounds give a zero of the integrand's shape (evaluated once at `a`).

The integrand may return a number, `ScalarField`, `VectorField`, `TensorField`,
or a vector/matrix with mixed numeric and `ScalarField` entries. The result
has the same scalar/field type or array shape. Field entries are returned in
**elementwise** representation; numbers stay numbers unless that entry is a
field at another quadrature point. Input values are never modified.

`onPhysicalGroup` optionally selects the spatial element support (e.g. a shell
surface); nodal data are copied to that group's element nodes. Every
participating elementwise field must cover the group. Without this keyword,
nodal fields use their Problem's material groups, and elementwise fields keep
their support. Each field-valued component must have compatible Problem,
element tags, local sizes and time labels at all quadrature points. Missing
element values are errors, not implicit zero values or nodal averages.

# Examples

```julia
A = sum(∫(ζ -> Q[k],     z[k], z[k+1]) for k in 1:length(Q))
B = sum(∫(ζ -> ζ*Q[k],   z[k], z[k+1]) for k in 1:length(Q))
D_b = sum(∫(ζ -> ζ^2*Q[k], z[k], z[k+1]) for k in 1:length(Q))
D_mb = [A B; B D_b]
K = ∫(B_mb' ⋅ D_mb ⋅ B_mb; Γ="shell")
```

`Q[k]` may contain numeric constants and spatially varying scalar fields.
See the [parameter integration guide](@ref parameter-integration).
"""
function integrate(f, a::Real, b::Real; order=5, onPhysicalGroup="")
    x, y = Float64(a), Float64(b)
    isfinite(x) && isfinite(y) || throw(ArgumentError("Parameter integration: bounds must be finite."))
    onPhysicalGroup isa AbstractString || throw(ArgumentError("Parameter integration: onPhysicalGroup must be a string."))
    group = String(onPhysicalGroup)
    points, weights = _parameter_quadrature(order)
    x == y && return _parameter_weighted(f(x), 0.0, group)
    half = y/2 - x/2
    center = x/2 + y/2
    result = _parameter_weighted(f(center + half*points[1]), half*weights[1], group)
    for k in 2:length(points)
        result = _parameter_accumulate(result, f(center + half*points[k]), half*weights[k], group)
    end
    return result
end

∫(f, a::Real, b::Real; kwargs...) = integrate(f, a, b; kwargs...)
