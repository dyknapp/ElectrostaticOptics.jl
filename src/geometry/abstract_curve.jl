# 2026 10 14, Daniel Knapp, Amsterdam : initial version

abstract type AbstractCurve2D end


val(ξ::Real) = ξ
val(ξ::ForwardDiff.Dual) = ForwardDiff.value(ξ)


function tangent(c::AbstractCurve2D, t::Real)::SVector{2}
    return SVector(normalize(ForwardDiff.derivative(c, t)))
end

function normal(c::AbstractCurve2D, t::Real)::SVector{2}
    τ = tangent(c, t)
    return SVector(τ[2], -τ[1])
end

function measure(c::AbstractCurve2D, t::Real)::Real
    return norm(ForwardDiff.derivative(c, t))
end


struct LineSeg{T <: Real} <: AbstractCurve2D
    p₁::SVector{2, T}
    p₂::SVector{2, T}
end
function (L::LineSeg)(ξ::Real)::SVector{2}
    @assert 0 ≤ val(ξ) ≤ 1
    return (1 - ξ) * L.p₁ + ξ * L.p₂
end

struct CircArc{T <: Real} <: AbstractCurve2D
    o::SVector{2, T}
    r::Real
    θ₁::Real
    θ₂::Real

    function CircArc(o::SVector{2}, r::Real, θ₁::Real, θ₂::Real)
        @assert abs(θ₂ - θ₁) ≤ 2π
    end
end
function (C::CircArc)(ξ::T)::SVector{2} where {T}
    @assert 0 ≤ ξ ≤ 1
    θ = (1 - ξ) * C.θ₁ + ξ * C.θ₂
    return SVector(C.o[1] + C.r * cos(θ), T.o[2] + C.r * sin(θ))
end

struct FuncCurve{F, T <: Real} <: AbstractCurve2D
    f::F # Must be a function t::Real → SVector{2,T} where t∈[t₁,t₂]
    # These get remapped so that you can call f(ξ∈[0,1])
    t₁::T
    t₂::T
end
function (C::FuncCurve)(ξ::T)::SVector{2} where {T}
    @assert 0 ≤ val(ξ) ≤ 1
    t = C.t₁ + ξ*(C.t₂ - C.t₁)
    return C.f(t)
end


@recipe function f(C::AbstractCurve2D)
    seriestype --> :path
    linecolor  --> :black
    arrow      --> true

    t = LinRange(0.0, 1.0, 128)
    pts = C.(t)

    return [pt[1] for pt in pts], [pt[2] for pt in pts]
end