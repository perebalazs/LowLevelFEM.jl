using LowLevelFEM
using LinearAlgebra
using SparseArrays
using Test
using Arpack
using IterativeSolvers, Preconditioners, IncompleteLU

include("helpers.jl")

@testset "LowLevelFEM.jl" begin
    @testset "Modern architecture (DSL)" begin
        include("test_fields.jl")
        include("test_operators.jl")
        include("test_weakforms.jl")
        include("test_multifield.jl")
        include("test_boundary_conditions.jl")
        include("test_solvefield.jl")
        include("test_solveeigen.jl")
        include("test_transient.jl")
        include("test_transient_remote.jl")
        include("test_transient_modify.jl")
        include("test_review_regressions.jl")
        include("test_contact.jl")
        include("test_validation.jl")
        include("test_directional_rect_mesh.jl")
    end
    @testset "Legacy" begin
        include("legacy/test_legacy_elasticity.jl")
        include("legacy/test_legacy_compatibility.jl")
    end
end
