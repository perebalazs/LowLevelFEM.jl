include("../src/LowLevelFEM.jl")
import LowLevelFEM as FEM
using Test

@testset "LowLevelFEM.jl" begin
    # Write your tests here.
    #@test LowLevelFEM.greet_your_package_name() == "Hello LowLevelFEM!"
    #@test LowLevelFEM.greet_your_package_name() != "Hello world!"
    returntest = "test", 0, 0, 0
    @test FEM.displacementConstraint("test", ux=0, uy=0, uz=0) == returntest
    display(FEM.displacementConstraint("test", ux=0, uy=0, uz=0))
end
