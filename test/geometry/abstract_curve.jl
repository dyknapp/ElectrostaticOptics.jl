using Test
using StaticArrays
using ElectrostaticOptics
using ForwardDiff
using LinearAlgebra
using JET

@testset "LineSeg geometry" begin
    p₁ = SVector(0.1, 1.1)
    p₂ = SVector(50.0, 5.0*π)
    line = LineSeg(p₁, p₂)

    @test line(0) ≈ p₁
    @test line(1) ≈ p₂

    @test_throws AssertionError line(-0.1)
    @test_throws AssertionError line(10.0)
    @test_throws AssertionError line(2)
end

@testset "tangent" begin
    p₁ = SVector(983.98, 1840.16)
    p₂ = SVector(-2090.0, -390.546778)
    line = LineSeg(p₁, p₂)

    tangent_vec = SVector((p₂[1] - p₁[1]), (p₂[2] - p₁[2]))
    tangent_vec = tangent_vec / sqrt((p₂[1] - p₁[1])^2 + (p₂[2] - p₁[2])^2)
    @test tangent(line, 0.5) ≈ tangent_vec
    @test tangent(line, 0) ≈ tangent_vec
    @test tangent(line, 1) ≈ tangent_vec
    @test tangent(line, 0.3475472309) ≈ tangent_vec

    f(t) = SVector(cosh(t), sinh(t))
    hyperbola = FuncCurve(f, -1.0, 1.0)
    @test tangent(hyperbola, 0.5) ≈ [0.0, 1.0]
    JET.@test_call tangent(line, 0.5)
end