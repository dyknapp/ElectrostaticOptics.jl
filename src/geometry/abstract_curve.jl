# 2026 10 14, Daniel Knapp, Amsterdam : initial version
# 2026 10 15, Daniel Knapp, Amsterdam : Attach to AbstractMapping, add affine_transform and translation
# 2026 10 18, Daniel Knapp, Amsterdam : Add EllipticArc

abstract type AbstractCurve2D <: AbstractMapping{2, 1} end

# Adapt reference coordinates to the existing scalar curve interface.
point(c::AbstractCurve2D, ξ::SVector{1}) = c(only(ξ))


val(ξ::Real) = ξ
val(ξ::ForwardDiff.Dual) = val(ForwardDiff.value(ξ))


function tangent(c::AbstractCurve2D, t::Real)::SVector{2}
    return normalize(jacobian(c, t)[:, 1])
end

# Derivative of an existing curve, intended for continuity checks.
# TODO: ForwardDiff is the baseline, curves should be given hand-written implementations
function derivative_curve(c::AbstractCurve2D)::AbstractCurve2D
    return FuncCurve(t -> jacobian(c, SVector(t))[:, 1], 0, 1)
end

function normal(c::AbstractCurve2D, t::Real)::SVector{2}
    τ = tangent(c, t)
    return SVector(τ[2], -τ[1])
end

# generic fallback to use FuncCurve's remapping
function trim(c::AbstractCurve2D, interval::Tuple{Real, Real})
    a, b = interval
    @assert 0 ≤ a < b ≤ 1
    return FuncCurve(c, a, b)
end

struct LineSeg{T <: Real} <: AbstractCurve2D
    p₁::SVector{2, T}
    p₂::SVector{2, T}
end
function (L::LineSeg)(ξ::Real)::SVector{2}
    @assert 0 ≤ val(ξ) ≤ 1
    return (1 - ξ) * L.p₁ + ξ * L.p₂
end
function affine_transform(L::LineSeg, mat::SMatrix{2, 2, T}) where {T}
    new_p₁ = mat * L.p₁
    new_p₂ = mat * L.p₂
    return LineSeg(new_p₁, new_p₂)
end
function translation(L::LineSeg, trans::SVector{2})
    new_p₁ = L.p₁ + trans
    new_p₂ = L.p₂ + trans
    return LineSeg(new_p₁, new_p₂)
end

struct CircArc{T <: Real} <: AbstractCurve2D
    o::SVector{2, T}
    r::Real
    θ₁::Real
    θ₂::Real

    function CircArc(o::SVector{2}, r::T, θ₁::T, θ₂::T) where {T}
        @assert abs(θ₂ - θ₁) ≤ 2π
        return new{T}(o, r, θ₁, θ₂)
    end
end
function (C::CircArc)(ξ::Real)::SVector{2}
    @assert 0 ≤ val(ξ) ≤ 1
    θ = (1 - ξ) * C.θ₁ + ξ * C.θ₂
    return SVector(C.o[1] + C.r * cos(θ), C.o[2] + C.r * sin(θ))
end

function affine_transform(C::CircArc, mat::SMatrix{2, 2})
    return EllipticArc(mat * C.o, C.r * mat, C.θ₁, C.θ₂)
end

function translation(C::CircArc, trans::SVector{2})
    o = C.o + trans
    return CircArc(o, C.r, C.θ₁, C.θ₂)
end

struct EllipticArc{T <: Real} <: AbstractCurve2D
    o::SVector{2, T}
    B::SMatrix{2, 2, T, 4}
    θ₁::T
    θ₂::T

    function EllipticArc(o::SVector{2, <:Real}, B::SMatrix{2, 2, <:Real}, θ₁::Real, θ₂::Real)
        Δθ = θ₂ - θ₁
        @assert abs(Δθ) ≤ 2 * oftype(float(Δθ), π)
        T = promote_type(eltype(o), eltype(B), typeof(θ₁), typeof(θ₂))
        return new{T}(o, B, θ₁, θ₂)
    end
end

function (C::EllipticArc)(ξ::Real)::SVector{2}
    @assert 0 ≤ val(ξ) ≤ 1
    θ = (1 - ξ) * C.θ₁ + ξ * C.θ₂
    return C.o + C.B * SVector(cos(θ), sin(θ))
end

function jacobian(C::EllipticArc, ξ::SVector{1})
    t = only(ξ)
    @assert 0 ≤ val(t) ≤ 1
    θ = (1 - t) * C.θ₁ + t * C.θ₂
    velocity = (C.θ₂ - C.θ₁) * (C.B * SVector(-sin(θ), cos(θ)))
    return SMatrix{2, 1}(velocity)
end

function affine_transform(C::EllipticArc, mat::SMatrix{2, 2})
    return EllipticArc(mat * C.o, mat * C.B, C.θ₁, C.θ₂)
end

function translation(C::EllipticArc, trans::SVector{2})
    return EllipticArc(C.o + trans, C.B, C.θ₁, C.θ₂)
end

function affine_transform(C::AbstractCurve2D, mat::SMatrix{2, 2}, trans::SVector{2})
    return translation(affine_transform(C, mat), trans)
end

struct FuncCurve{F, T <: Real} <: AbstractCurve2D
    f::F # Must be a function t::Real → SVector{2,T} where t∈[t₁,t₂]
    # These get remapped so that you can call f(ξ∈[0,1])
    t₁::T
    t₂::T
end
function (C::FuncCurve)(ξ::Real)::SVector{2}
    @assert 0 ≤ val(ξ) ≤ 1
    t = C.t₁ + ξ * (C.t₂ - C.t₁)
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
