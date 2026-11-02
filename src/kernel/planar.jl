# 2026 11 2, Daniel Knapp, Amsterdam : initial version

# General rule:
# Kernel should own the quadrature rules.

# Avoid division by caching log(r₀)
struct PlanarKernel{T} <: AbstractKernel{T}
    kₑ::T
    log_r₀::T
end

function PlanarKernel(::Type{T}, kₑ::T, r₀::T) where {T <: Real}
    return PlanarKernel(
        kₑ,
        log(r₀),
    )
end

function (kernel::PlanarKernel{T})(pt₁::SVector{2, T}, pt₂::SVector{2, T}) where {T}
    r = norm(pt₁ - pt₂)
    return 2 * kernel.kₑ * (kernel.log_r₀ - log(r))
end
