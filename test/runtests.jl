using LowLevelFEM
using Test

@testset "LowLevelFEM.jl" begin
    # Write your tests here.
    @test LowLevelFEM.greet_your_package_name() == "Hello LowLevelFEM!"
    @test LowLevelFEM.greet_your_package_name() != "Hello world!"
end
