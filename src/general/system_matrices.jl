"""
    SystemMatrix(A::SparseMatrixCSC{Float64}, trial::Problem, test::Problem=trial)
    SystemMatrix(blocks::Matrix{SystemMatrix})
    SystemMatrix(K::SystemMatrix, trial::Problem, test::Problem=trial)

Sparse finite-element matrix with trial/test field metadata.

The sparse-array constructor wraps an existing matrix. Rebinding an existing
`SystemMatrix` changes metadata while sharing its `.A` array; it does not copy
entries. Matrix multiplication by a compatible field returns a nodal field.
`Symmetric(K)` declares a symmetric system for solver dispatch; it does not
verify the physical definiteness of the problem.

The block constructor builds a global matrix from `SystemMatrix` blocks.
Literal zero blocks are accepted where the row/column Problems can be inferred
from other blocks. `.problems` and `.offsets` record the block order.

This function assembles a global sparse matrix from a block structure where each
block represents the discretization of a bilinear form between two fields
(test and trial problems). The blocks are assumed to correspond to different
fields in a multifield formulation.

The function performs the following steps:

1. Determine the ordering of **trial problems** from the columns of the block matrix.

2. Determine the ordering of **test problems** from the rows of the block matrix.

3. Verify that the block system is **square in the field sense**, i.e.
   
       trial_problems[i] === test_problems[i]
   
   for all fields.

4. Compute global DOF offsets for each field.

5. Assemble the global sparse matrix by inserting the entries of each block
   at the appropriate offset position.

Each block `blocks[i,j]` must satisfy

    blocks[i,j].test_model  → test Problem (row space)
    blocks[i,j].model       → trial Problem (column space)
    blocks[i,j].A           → sparse block matrix

The resulting global matrix has the structure

    [ A₁₁  A₁₂  ⋯ ]
    [ A₂₁  A₂₂  ⋯ ]
    [ ⋮    ⋮    ⋱ ]

where `Aᵢⱼ` couples the `j`-th trial field with the `i`-th test field.

# Arguments

- `blocks::Matrix{SystemMatrix}`: Matrix of block matrices representing the
  bilinear operators between fields.

# Returns

- `SystemMatrix`: A global system matrix containing the assembled sparse matrix
  and metadata about the involved `Problem`s and their DOF offsets.

# Errors

An error is thrown if:

- no trial or test problems are found in the block matrix,
- the block system is not square in the field sense,
- a block contains inconsistent matrix dimensions.

# Example

```julia
Auu = ∫(ε(Pu) ⋅ C ⋅ ε(Pu))
Aup = ∫(Div(Pu) ⋅ Id(Pp))
Apu = ∫(Div(Pu) ⋅ Id(Pp))'
App = ∫(Grad(Pp) ⋅ Grad(Pp) * δ)

K = SystemMatrix([
    Auu  Aup
    Apu  App
])
```

This assembles the global block system

```
[ Kuu  Kup ]
[ Kpu  Kpp ]
```

For replacing field and test_field, use:

```julia
M = SystemMatrix(M, new_field)
```

or

```julia
M = SystemMatrix(M, new_field, new_test_field)
```

corresponding to a two-field problem.  
"""
struct SystemMatrix
    A::SparseMatrixCSC
    model::Union{Problem,Nothing}
    test_model::Union{Problem,Nothing}
    problems::Union{Vector{Problem},Nothing}      # multifield metadata
    offsets::Union{Vector{Int},Nothing}           # global starting indices

    # ----------------------------------------------------------
    # Single-field constructors
    # ----------------------------------------------------------

    SystemMatrix(
        A::SparseMatrixCSC{Float64},
        model::Problem,
        test_model::Problem
    ) = new(A, model, test_model, nothing, nothing)

    SystemMatrix(
        A::SystemMatrix,
        model::Problem
    ) = new(A.A, model, model, nothing, nothing)

    SystemMatrix(
        A::SystemMatrix,
        model::Problem,
        test_model::Problem
    ) = new(A.A, model, test_model, nothing, nothing)

    SystemMatrix(
        A::SparseMatrixCSC{Float64},
        model::Problem
    ) = new(A, model, model, nothing, nothing)

    SystemMatrix(
        A::SparseMatrixCSC{Float64},
        model::Union{Problem,Nothing},
        test_model::Union{Problem,Nothing},
        problems::Union{Vector{Problem},Nothing},
        offsets::Union{Vector{Int},Nothing}
    ) = new(A, model, test_model, problems, offsets)


    # ----------------------------------------------------------
    # Multifield block constructor
    #
    # The global sparse matrix is assembled directly in CSC
    # format from the already assembled CSC blocks.
    #
    # No IJV representation and no call to sparse(I,J,V) are
    # used here.
    # ----------------------------------------------------------

    function SystemMatrix(blocks::Matrix{SystemMatrix})

        nrows, ncols = size(blocks)

        # ------------------------------------------------------
        # 1) Collect trial problems (column-based order)
        # ------------------------------------------------------
        trial_problems = Problem[]

        function push_unique!(vec, P)
            P === nothing && return
            if all(q -> q !== P, vec)
                push!(vec, P)
            end
        end

        for j in 1:ncols
            for i in 1:nrows
                blk = blocks[i, j]
                if blk.model !== nothing
                    push_unique!(trial_problems, blk.model)
                    break
                end
            end
        end

        isempty(trial_problems) &&
            error("No trial Problems found in block matrix.")

        # ------------------------------------------------------
        # 2) Collect test problems (row-based order)
        # ------------------------------------------------------
        test_problems = Problem[]

        for i in 1:nrows
            for j in 1:ncols
                blk = blocks[i, j]
                if blk.test_model !== nothing
                    push_unique!(test_problems, blk.test_model)
                    break
                end
            end
        end

        isempty(test_problems) &&
            error("No test Problems found in block matrix.")

        # ------------------------------------------------------
        # 3) Field-level square check
        # ------------------------------------------------------
        if length(trial_problems) != length(test_problems) ||
           any(trial_problems[i] !== test_problems[i]
               for i in eachindex(trial_problems))
            error(
                "Block system is not square in field sense. " *
                "Trial and test spaces differ."
            )
        end

        problems = trial_problems
        nproblems = length(problems)

        # ------------------------------------------------------
        # 4) Compute global offsets
        # ------------------------------------------------------
        offsets = Vector{Int}(undef, nproblems)
        offsets[1] = 0

        for i in 2:nproblems
            offsets[i] =
                offsets[i - 1] + ndofs(problems[i - 1])
        end

        total_dofs =
            offsets[end] + ndofs(problems[end])

        # ------------------------------------------------------
        # Helper: find Problem index by object identity
        # ------------------------------------------------------
        function problem_index(P)
            for (k, q) in enumerate(problems)
                if q === P
                    return k
                end
            end

            error("Problem not found in metadata.")
        end

        # ------------------------------------------------------
        # 5) Map input blocks to their actual field positions
        #
        # blockmap[i,j] corresponds to:
        #
        #   test  = problems[i]
        #   trial = problems[j]
        #
        # This makes the CSC assembly independent of the
        # physical ordering of the supplied block matrix.
        # ------------------------------------------------------
        blockmap =
            Matrix{Union{SystemMatrix,Nothing}}(
                nothing,
                nproblems,
                nproblems
            )

        for bi in 1:nrows
            for bj in 1:ncols

                blk = blocks[bi, bj]

                rowP = blk.test_model
                colP = blk.model

                rowP === nothing &&
                    error("Block ($bi,$bj) is missing test_model.")

                colP === nothing &&
                    error("Block ($bi,$bj) is missing model.")

                iP = problem_index(rowP)
                jP = problem_index(colP)

                # Each field pair must already have been reduced
                # to a single SystemMatrix before this constructor.
                blockmap[iP, jP] === nothing ||
                    error(
                        "Multiple blocks found for the same " *
                        "(test_model, model) pair."
                    )

                A = blk.A

                expected_size =
                    (ndofs(rowP), ndofs(colP))

                size(A) == expected_size ||
                    error(
                        "Block size mismatch in block ($bi,$bj): " *
                        "got $(size(A)), expected $expected_size."
                    )

                blockmap[iP, jP] = blk
            end
        end

        # ------------------------------------------------------
        # 6) Determine exact number of stored nonzeros
        # ------------------------------------------------------
        total_nnz = 0

        for iP in 1:nproblems
            for jP in 1:nproblems
                blk = blockmap[iP, jP]
                blk === nothing && continue

                total_nnz += nnz(blk.A)
            end
        end

        # ------------------------------------------------------
        # 7) Allocate final CSC storage
        # ------------------------------------------------------
        colptr = Vector{Int}(undef, total_dofs + 1)
        rowval = Vector{Int}(undef, total_nnz)
        nzval  = Vector{Float64}(undef, total_nnz)

        pos = 1
        colptr[1] = 1

        # ------------------------------------------------------
        # 8) Assemble the global matrix column by column
        #
        # CSC requirement:
        #
        #   - columns are written consecutively
        #   - row indices inside each column remain sorted
        #
        # Since row blocks are traversed in global row order
        # and each input CSC block already has sorted row indices,
        # the resulting global row indices are sorted as well.
        # ------------------------------------------------------
        @inbounds for jP in 1:nproblems

            ncol_local = ndofs(problems[jP])
            col_offset = offsets[jP]

            for local_col in 1:ncol_local

                global_col = col_offset + local_col

                for iP in 1:nproblems

                    blk = blockmap[iP, jP]
                    blk === nothing && continue

                    A = blk.A
                    row_offset = offsets[iP]

                    p1 = A.colptr[local_col]
                    p2 = A.colptr[local_col + 1] - 1

                    for p in p1:p2
                        rowval[pos] =
                            row_offset + A.rowval[p]

                        nzval[pos] =
                            A.nzval[p]

                        pos += 1
                    end
                end

                colptr[global_col + 1] = pos
            end
        end

        # Defensive consistency check
        pos == total_nnz + 1 ||
            error(
                "Internal CSC assembly error: " *
                "expected $total_nnz stored entries, " *
                "but wrote $(pos - 1)."
            )

        # ------------------------------------------------------
        # 9) Construct SparseMatrixCSC directly
        # ------------------------------------------------------
        A_big = SparseMatrixCSC(
            total_dofs,
            total_dofs,
            colptr,
            rowval,
            nzval
        )

        return new(
            A_big,
            nothing,
            nothing,
            problems,
            offsets
        )
    end
end

"""
    _nodesOnPhysicalGroup(problem::Problem, physicalName::String)

Return the node tags belonging to a Gmsh physical group.

The node tags are collected from all finite elements of the physical group,
duplicates are removed, and the result is sorted in ascending node-tag order.

This ordering is intentionally identical to the node ordering used by the
contact formulation.
"""
function _nodesOnPhysicalGroup(
    problem::Problem,
    physicalName::String
    )

    gmsh.model.setCurrent(problem.name)

    dimTags =
        gmsh.model.getEntitiesForPhysicalName(physicalName)

    isempty(dimTags) &&
        error(
            "_nodesOnPhysicalGroup: physical group " *
            "'$physicalName' was not found."
        )

    nodeTags = Int[]

    for (dim, tag) in dimTags

        _, _, nodesByType =
            gmsh.model.mesh.getElements(dim, tag)

        for nodes in nodesByType
            append!(nodeTags, nodes)
        end
    end

    isempty(nodeTags) &&
        error(
            "_nodesOnPhysicalGroup: no mesh nodes were found " *
            "in physical group '$physicalName'."
        )

    sort!(unique!(nodeTags))

    return nodeTags
end

"""
    _nodalDoFs(problem::Problem, nodeTags::AbstractVector{<:Integer})

Return the global degrees of freedom associated with the supplied node tags.

The node order is preserved. For a vector field, all components belonging
to one node remain consecutive.

For example, for a three-dimensional vector field and

    nodeTags = [4, 7]

the returned DoF ordering is

    [4x, 4y, 4z, 7x, 7y, 7z]

in terms of the global finite-element numbering.
"""
function _nodalDoFs(
    problem::Problem,
    nodeTags::AbstractVector{<:Integer}
    )

    pdim = problem.pdim
    dofs = Vector{Int}(undef, pdim * length(nodeTags))

    k = 1

    @inbounds for nodeTag0 in nodeTags

        nodeTag = Int(nodeTag0)

        1 <= nodeTag <= problem.non ||
            error(
                "_nodalDoFs: node tag $nodeTag is outside the " *
                "valid range 1:$(problem.non)."
            )

        offset = (nodeTag - 1) * pdim

        for component in 1:pdim
            dofs[k] = offset + component
            k += 1
        end
    end

    return dofs
end

"""
    subSystemMatrix(K::SystemMatrix, nodeTags::AbstractVector{<:Integer})

Restrict a field-level system matrix to the supplied mesh nodes.

Rows are selected from the test field and columns from the trial field.
The order of `nodeTags` is preserved, therefore the function can also be
used as a permutation of the original nodal system.

The returned matrix acts on an anonymous reduced algebraic space and
therefore has `model === nothing` and `test_model === nothing`.

This operation is intended for field-level `SystemMatrix` objects and is
not defined for already assembled multifield systems.
"""
function subSystemMatrix(
    K::SystemMatrix,
    nodeTags::AbstractVector{<:Integer}
    )

    K.problems === nothing ||
        error(
            "subSystemMatrix: reduction of an assembled multifield " *
            "SystemMatrix is not supported."
        )

    K.model === nothing &&
        error(
            "subSystemMatrix: the SystemMatrix has no trial Problem."
        )

    K.test_model === nothing &&
        error(
            "subSystemMatrix: the SystemMatrix has no test Problem."
        )

    expectedSize =
        (ndofs(K.test_model), ndofs(K.model))

    size(K.A) == expectedSize ||
        error(
            "subSystemMatrix: matrix size $(size(K.A)) is incompatible " *
            "with its test/trial Problems; expected $expectedSize."
        )

    rowDoFs = _nodalDoFs(K.test_model, nodeTags)
    colDoFs = _nodalDoFs(K.model, nodeTags)

    A = K.A[rowDoFs, colDoFs]

    dropzeros!(A)

    return SystemMatrix(
        A,
        nothing,
        nothing,
        nothing,
        nothing
    )
end

"""
    subSystemMatrix(K::SystemMatrix, physicalName::String)

Restrict a field-level system matrix to the nodes of a Gmsh physical group.

Nodes are ordered by ascending Gmsh node tag. This is the same ordering used
by the contact formulation, so the returned matrix is directly compatible
with contact-space operators constructed on the same physical group.
"""
function subSystemMatrix(
    K::SystemMatrix,
    physicalName::String
    )

    K.model === nothing &&
        error(
            "subSystemMatrix: the SystemMatrix has no trial Problem."
        )

    K.test_model === nothing &&
        error(
            "subSystemMatrix: the SystemMatrix has no test Problem."
        )

    K.model.name == K.test_model.name ||
        error(
            "subSystemMatrix: trial and test Problems must belong " *
            "to the same Gmsh model when selecting a physical group."
        )

    nodeTags =
        _nodesOnPhysicalGroup(
            K.model,
            physicalName
        )

    return subSystemMatrix(
        K,
        nodeTags
    )
end

"""
    subSystemMatrix(K::SystemMatrix; onPhysicalGroup::String)

Keyword form of [`subSystemMatrix`](@ref).
"""
function subSystemMatrix(
    K::SystemMatrix;
    onPhysicalGroup::String
    )

    return subSystemMatrix(
        K,
        onPhysicalGroup
    )
end

"""
    Base.show(io::IO, M::SystemMatrix)

Internal function to display `SystemMatrix` as a `SparseMatrixCSC{Float64}`.
"""
function Base.show(io::IO, M::SystemMatrix)
    show(io, M.A)
end

function Base.show(io::IO, ::MIME"text/plain", M::SystemMatrix)
    show(io, M.A)
    #println(io, "SystemMatrix ($(size(M.A,1)) × $(size(M.A,2))), nnz = $(nnz(M.A))")
end

import Base:display
function display(M::SystemMatrix)
    display(M.A)
end

