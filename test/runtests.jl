using ElectrostaticOptics
using Test
using Aqua
using JET

@testset "ElectrostaticOptics.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ElectrostaticOptics, deps_compat=(check_extras=false,))
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(ElectrostaticOptics)
    end
    # Write your tests here.
end
