using LowLevelFEM
using LinearAlgebra
using SparseArrays
using Test

@testset "LowLevelFEM.jl" begin
    include("test_fields.jl")
    include("test_operators.jl")
    include("test_assembly.jl")
    include("test_validation.jl")
    include("test_cube.jl")
end