export ∇

@inline function _physical_gradient_new(J::AbstractVector{<:Real}, offset::Int,
                                        dim::Int, g1::Float64,
                                        g2::Float64, g3::Float64)
    a1 = J[offset + 1]
    a2 = J[offset + 2]
    a3 = J[offset + 3]

    if dim == 1
        d = a1*a1 + a2*a2 + a3*a3
        q = g1 / d
        return a1*q, a2*q, a3*q
    end

    b1 = J[offset + 4]
    b2 = J[offset + 5]
    b3 = J[offset + 6]

    if dim == 2
        aa = a1*a1 + a2*a2 + a3*a3
        ab = a1*b1 + a2*b2 + a3*b3
        bb = b1*b1 + b2*b2 + b3*b3
        d = aa*bb - ab*ab
        q1 = (bb*g1 - ab*g2) / d
        q2 = (aa*g2 - ab*g1) / d
        return a1*q1 + b1*q2,
               a2*q1 + b2*q2,
               a3*q1 + b3*q2
    end

    c1 = J[offset + 7]
    c2 = J[offset + 8]
    c3 = J[offset + 9]

    bc1 = b2*c3 - b3*c2
    bc2 = b3*c1 - b1*c3
    bc3 = b1*c2 - b2*c1
    ca1 = c2*a3 - c3*a2
    ca2 = c3*a1 - c1*a3
    ca3 = c1*a2 - c2*a1
    ab1 = a2*b3 - a3*b2
    ab2 = a3*b1 - a1*b3
    ab3 = a1*b2 - a2*b1
    d = a1*bc1 + a2*bc2 + a3*bc3

    return (g1*bc1 + g2*ca1 + g3*ab1) / d,
           (g1*bc2 + g2*ca2 + g3*ab2) / d,
           (g1*bc3 + g2*ca3 + g3*ab3) / d
end

"""
    ∇(r::Union{ScalarField,VectorField,TensorField}; nabla=:grad)

Compute spatial derivatives at the element nodes and return an element-wise
field. This experimental implementation retrieves the Jacobians of all
elements of a given type and entity in one Gmsh call and evaluates the field
derivatives directly, without assembling element derivative matrices.

Nodal input is converted to element-wise storage with `nodesToElements`, which
preserves its values at the element nodes. Element-wise input is used directly,
so discontinuities between adjacent elements are preserved. The result is
always element-wise because derivatives can also be discontinuous between
adjacent elements. For element-wise input, `r.numElem` defines the domain of
the operation; consequently, fields restricted to a surface, curve, or another
physical group remain restricted to that group.

Supported operations:

- `ScalarField`, `nabla=:grad` -> `VectorField`
- three-component `VectorField`, `nabla=:grad` -> `TensorField`
- three-component `VectorField`, `nabla=:div` -> `ScalarField`
- three-component `VectorField`, `nabla=:curl` -> `VectorField`
- `TensorField`, `nabla=:div` -> `VectorField`
- axisymmetric two-component `VectorField`, `nabla=:grad` -> `TensorField`
"""
function ∇(rr::Union{ScalarField,VectorField,TensorField}; nabla=:grad)
    problem = rr.model
    gmsh.model.setCurrent(problem.name)

    r = isElementwise(rr) ? rr : nodesToElements(rr)
    nsteps = r.nsteps

    mode, output_size = if r isa ScalarField && nabla == :grad
        (:scalar_grad, 3)
    elseif r isa VectorField && r.type == :v3D && nabla == :grad
        (:vector_grad, 9)
    elseif r isa VectorField && r.type == :v3D && nabla == :div
        (:vector_div, 1)
    elseif r isa VectorField && r.type == :v3D && nabla == :curl
        (:vector_curl, 3)
    elseif r isa TensorField && nabla == :div
        (:tensor_div, 3)
    elseif r isa VectorField && r.type == :v2D &&
           problem.type == :AxiSymmetric && nabla == :grad
        (:axisymmetric_grad, 9)
    else
        error("∇: unsupported field type, field layout, or operator: " *
              "$(typeof(r)), $(r.type), nabla=$nabla, problem type=$(problem.type).")
    end

    number_of_elements = length(r.numElem)
    number_of_elements > 0 || error("∇: the input field contains no elements.")
    values = Vector{Matrix{Float64}}(undef, number_of_elements)
    numElem = copy(r.numElem)
    field_index = Dict{Int,Int}()
    sizehint!(field_index, number_of_elements)
    for (i, elem) in pairs(numElem)
        haskey(field_index, elem) &&
            error("∇: duplicate element tag $elem in the input field.")
        field_index[elem] = i
    end
    processed = falses(number_of_elements)
    number_processed = 0

    # The element list stored by the field defines the differentiation domain.
    # This is essential for fields restricted to a surface or another physical
    # group: the Problem's material groups need not describe the field domain.
    _, _, field_dim, _ = gmsh.model.mesh.getElement(first(numElem))

    for (edim, etag) in gmsh.model.getEntities(field_dim)
        elemTypes, elemTags, _ = gmsh.model.mesh.getElements(edim, etag)

        for itype in eachindex(elemTypes)
            tags = elemTags[itype]
            any_selected = false
            @inbounds for elem0 in tags
                elem = Int(elem0)
                if haskey(field_index, elem)
                    any_selected = true
                    break
                end
            end
            any_selected || continue

            et = Int(elemTypes[itype])
            dim_et, numNodes, nodeCoord = _get_props_cached(et)
            grad_h, h = _get_basis_cached(et, nodeCoord, numNodes)
            nel = length(tags)

            jac_all, _, coord_all =
                gmsh.model.mesh.getJacobians(et, nodeCoord, etag)

            expected_jacobians = 9 * numNodes * nel
            length(jac_all) == expected_jacobians ||
                error("∇: unexpected number of Jacobian entries for element type $et " *
                      "on entity $etag: $(length(jac_all)) instead of $expected_jacobians.")

            @inbounds for ie in 1:nel
                elem = Int(tags[ie])
                output_index = get(field_index, elem, 0)
                output_index == 0 && continue
                processed[output_index] &&
                    error("∇: mesh element $elem was encountered more than once.")
                processed[output_index] = true
                number_processed += 1

                re = r.A[output_index]
                e = zeros(Float64, output_size * numNodes, nsteps)
                jac_element_offset = (ie - 1) * 9 * numNodes
                coord_element_offset = (ie - 1) * 3 * numNodes

            if mode === :scalar_grad
                for k in 1:numNodes, l in 1:numNodes
                    gbase = 3*(l - 1)
                    jbase = jac_element_offset + 9*(k - 1)
                    gx, gy, gz = _physical_gradient_new(
                        jac_all, jbase, dim_et,
                        grad_h[gbase + 1, k],
                        grad_h[gbase + 2, k],
                        grad_h[gbase + 3, k])
                    for step in 1:nsteps
                        v = re[l, step]
                        e[3k - 2, step] += gx*v
                        e[3k - 1, step] += gy*v
                        e[3k,     step] += gz*v
                    end
                end

            elseif mode === :vector_grad
                for k in 1:numNodes, l in 1:numNodes
                    gbase = 3*(l - 1)
                    jbase = jac_element_offset + 9*(k - 1)
                    gx, gy, gz = _physical_gradient_new(
                        jac_all, jbase, dim_et,
                        grad_h[gbase + 1, k],
                        grad_h[gbase + 2, k],
                        grad_h[gbase + 3, k])
                    nodebase = 3*(l - 1)
                    for step in 1:nsteps
                        vx = re[nodebase + 1, step]
                        vy = re[nodebase + 2, step]
                        vz = re[nodebase + 3, step]
                        obase = 9*(k - 1)
                        e[obase + 1, step] += gx*vx
                        e[obase + 2, step] += gx*vy
                        e[obase + 3, step] += gx*vz
                        e[obase + 4, step] += gy*vx
                        e[obase + 5, step] += gy*vy
                        e[obase + 6, step] += gy*vz
                        e[obase + 7, step] += gz*vx
                        e[obase + 8, step] += gz*vy
                        e[obase + 9, step] += gz*vz
                    end
                end

            elseif mode === :vector_div
                for k in 1:numNodes, l in 1:numNodes
                    gbase = 3*(l - 1)
                    jbase = jac_element_offset + 9*(k - 1)
                    gx, gy, gz = _physical_gradient_new(
                        jac_all, jbase, dim_et,
                        grad_h[gbase + 1, k],
                        grad_h[gbase + 2, k],
                        grad_h[gbase + 3, k])
                    nodebase = 3*(l - 1)
                    for step in 1:nsteps
                        e[k, step] += gx*re[nodebase + 1, step] +
                                      gy*re[nodebase + 2, step] +
                                      gz*re[nodebase + 3, step]
                    end
                end

            elseif mode === :vector_curl
                for k in 1:numNodes, l in 1:numNodes
                    gbase = 3*(l - 1)
                    jbase = jac_element_offset + 9*(k - 1)
                    gx, gy, gz = _physical_gradient_new(
                        jac_all, jbase, dim_et,
                        grad_h[gbase + 1, k],
                        grad_h[gbase + 2, k],
                        grad_h[gbase + 3, k])
                    nodebase = 3*(l - 1)
                    for step in 1:nsteps
                        vx = re[nodebase + 1, step]
                        vy = re[nodebase + 2, step]
                        vz = re[nodebase + 3, step]
                        obase = 3*(k - 1)
                        e[obase + 1, step] += gy*vz - gz*vy
                        e[obase + 2, step] += gz*vx - gx*vz
                        e[obase + 3, step] += gx*vy - gy*vx
                    end
                end

            elseif mode === :tensor_div
                for k in 1:numNodes, l in 1:numNodes
                    gbase = 3*(l - 1)
                    jbase = jac_element_offset + 9*(k - 1)
                    gx, gy, gz = _physical_gradient_new(
                        jac_all, jbase, dim_et,
                        grad_h[gbase + 1, k],
                        grad_h[gbase + 2, k],
                        grad_h[gbase + 3, k])
                    nodebase = 9*(l - 1)
                    for step in 1:nsteps
                        obase = 3*(k - 1)
                        e[obase + 1, step] += gx*re[nodebase + 1, step] +
                                                  gy*re[nodebase + 4, step] +
                                                  gz*re[nodebase + 7, step]
                        e[obase + 2, step] += gx*re[nodebase + 2, step] +
                                                  gy*re[nodebase + 5, step] +
                                                  gz*re[nodebase + 8, step]
                        e[obase + 3, step] += gx*re[nodebase + 3, step] +
                                                  gy*re[nodebase + 6, step] +
                                                  gz*re[nodebase + 9, step]
                    end
                end

            else # :axisymmetric_grad
                for k in 1:numNodes, l in 1:numNodes
                    gbase = 3*(l - 1)
                    jbase = jac_element_offset + 9*(k - 1)
                    gx, gy, _ = _physical_gradient_new(
                        jac_all, jbase, dim_et,
                        grad_h[gbase + 1, k],
                        grad_h[gbase + 2, k],
                        grad_h[gbase + 3, k])
                    nodebase = 2*(l - 1)
                    radius = coord_all[coord_element_offset + 3k - 2]
                    hoop = abs(radius) < 1e-10 ? 0.0 : h[l, k] / radius
                    for step in 1:nsteps
                        ur = re[nodebase + 1, step]
                        uz = re[nodebase + 2, step]
                        obase = 9*(k - 1)
                        e[obase + 1, step] += gx*ur
                        shear = 0.5*(gy*ur + gx*uz)
                        e[obase + 2, step] += shear
                        e[obase + 4, step] += shear
                        e[obase + 5, step] += gy*uz
                        e[obase + 9, step] += hoop*ur
                    end
                end
            end

                values[output_index] = e
            end
        end
    end

    number_processed == number_of_elements ||
        error("∇: processed $number_processed mesh elements, but the input field " *
              "contains $number_of_elements elements.")

    if mode === :vector_grad || mode === :axisymmetric_grad
        return TensorField(values, [;;], r.t, numElem, nsteps, :tensor, problem)
    elseif mode === :vector_div
        return ScalarField(values, [;;], r.t, numElem, nsteps, :scalar, problem)
    else
        return VectorField(values, [;;], r.t, numElem, nsteps, :v3D, problem)
    end
end
