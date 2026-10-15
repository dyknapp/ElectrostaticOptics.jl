using Test
using StaticArrays
using ElectrostaticOptics
using ForwardDiff
using LinearAlgebra
using JET
using Supposition

ϵ_test = 10 * eps()
min_len = 1.0e-6
max_len = 1.0e+6

coord_gen = Supposition.Data.Floats{Float64}(minimum = -max_len, maximum = max_len, nans = false, infs = false)
vec_gen = map(SVector{2, Float64}, Supposition.Data.Vectors(coord_gen, min_size = 2, max_size = 2))
t_gen = Supposition.Data.Floats{Float64}(minimum = 0.0, maximum = 1.0, nans = false, infs = false)
len_gen = Supposition.Data.Floats{Float64}(minimum = min_len, maximum = max_len, nans = false, infs = false)
angle_gen = Supposition.Data.Floats{Float64}(minimum = 0.0, maximum = 2.0 * π, nans = false, infs = false)

@testset "LineSeg geometry" begin
    p₁ = SVector(0.1, 1.1)
    p₂ = SVector(50.0, 5.0 * π)
    line = LineSeg(p₁, p₂)

    @test line(0) ≈ p₁
    @test line(1) ≈ p₂

    @test_throws AssertionError line(-0.1)
    @test_throws AssertionError line(10.0)
    @test_throws AssertionError line(2)

    @check function line_is_linear(p₁ = vec_gen, p₂ = vec_gen, t = t_gen)
        if norm(p₂ - p₁) ≤ ϵ_test
            reject!()
        end

        scale = max(norm(p₁), norm(p₂))
        line = LineSeg(p₁, p₂)

        actual = (line(t) - p₁)
        expected = t * (p₂ - p₁)

        residual = norm(actual - expected)
        tolerance = scale * ϵ_test
        event!("residual", residual)
        event!("tolerance", tolerance)

        return residual ≤ tolerance
    end
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

    @check function line_tangent_constant(p₁ = vec_gen, p₂ = vec_gen, t = t_gen)
        if norm(p₂ - p₁) ≤ 1.0e-6
            reject!()
        end

        tangent_vec = SVector((p₂[1] - p₁[1]), (p₂[2] - p₁[2]))
        tangent_vec = tangent_vec / sqrt((p₂[1] - p₁[1])^2 + (p₂[2] - p₁[2])^2)

        scale = max(norm(p₁), norm(p₂))
        line = LineSeg(p₁, p₂)

        residual = norm(tangent(line, t) - tangent_vec)
        tolerance = max(scale * ϵ_test, ϵ_test)
        event!("residual", residual)
        event!("tolerance", tolerance)

        return residual ≤ tolerance
    end

    @check function tangent_length(p₁ = vec_gen, p₂ = vec_gen, t = t_gen)
        if norm(p₂ - p₁) ≤ 1.0e-6
            reject!()
        end
        line = LineSeg(p₁, p₂)
        tang = tangent(line, t)
        return abs(1 - norm(tang)) < ϵ_test
    end

    @check function tangent_normal_angle(p₁ = vec_gen, p₂ = vec_gen, t = t_gen)
        if norm(p₂ - p₁) ≤ 1.0e-6
            reject!()
        end

        line = LineSeg(p₁, p₂)

        result = dot(tangent(line, t), normal(line, t))
        event!("dot product", result)

        return result ≤ ϵ_test
    end

    f(t) = SVector(cosh(t), sinh(t))
    hyperbola = FuncCurve(f, -1.0, 1.0)
    @test tangent(hyperbola, 0.5) ≈ [0.0, 1.0]
    JET.@test_call tangent(line, 0.5)
end

@testset "normal" begin
    @check function normal_length(p₁ = vec_gen, p₂ = vec_gen, t = t_gen)
        if norm(p₂ - p₁) ≤ 1.0e-6
            reject!()
        end
        line = LineSeg(p₁, p₂)
        nor = normal(line, t)
        return abs(1 - norm(nor)) < ϵ_test
    end
end

@testset "circle" begin
    o = SVector(0.273738, 28347.83784)
    r = 3854.23982789
    θ₁ = 0.0
    θ₂ = 1.0 * π

    arc = CircArc(o, r, θ₁, θ₂)
    @test norm(arc(0.94839) - o) ≈ r

    @check function circle_radius_constant(o = vec_gen, r = len_gen, θ₁ = angle_gen, θ₂ = angle_gen, t = t_gen)
        arc = CircArc(o, r, θ₁, θ₂)

        if abs(θ₂ - θ₁) < ϵ_test
            reject!()
        end

        circ_rad = norm(arc(t) - o)
        residual = abs(circ_rad - r)
        tolerance = max(r * ϵ_test, norm(o * ϵ_test), ϵ_test)
        event!("actual radius", circ_rad)
        event!("residual", residual)
        event!("tolerance", tolerance)

        return residual ≤ tolerance
    end
end

@testset "measure" begin
    @check function line_measure(p₁ = vec_gen, p₂ = vec_gen, t = t_gen)
        if norm(p₂ - p₁) ≤ 1.0e-6
            reject!()
        end
        line = LineSeg(p₁, p₂)
        expected = norm(p₂ - p₁)
        meas = measure(line, t)
        residual = abs(expected - meas)
        tolerance = max(expected * ϵ_test, ϵ_test)

        event!("measure", meas)
        event!("residual", residual)
        event!("tolerance", tolerance)

        return residual ≤ tolerance
    end

    @check function circle_measure(o = vec_gen, r = len_gen, θ₁ = angle_gen, θ₂ = angle_gen, t = t_gen)
        arc = CircArc(o, r, θ₁, θ₂)

        if abs(θ₂ - θ₁) < ϵ_test
            reject!()
        end

        expected = r * abs(θ₂ - θ₁)
        meas = measure(arc, t)

        residual = abs(expected - meas)
        tolerance = max(expected * ϵ_test, ϵ_test)

        event!("measure", meas)
        event!("expected", expected)
        event!("tolerance", tolerance)

        return residual ≤ tolerance
    end
end
