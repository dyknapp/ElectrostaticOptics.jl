using ElectrostaticOptics
using Test
using Aqua
using JET

@testset "ElectrostaticOptics.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ElectrostaticOptics)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(ElectrostaticOptics; target_defined_modules = true)
    end
    # Write your tests here.
end
