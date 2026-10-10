# 2026 10 04, Daniel Knapp, Amsterdam : initial version (mostly writing out conventions)
# 2026 10 10, Daniel Knapp, Amsterdam : Build AbstractMapping (delete initial version)

using ForwardDiff
using StaticArrays

# =====================================================================================
# AbstractMapping
#
# Hold onto the basic differential geometry stuff in this abstract supertype
#
# A mapping ℝᴰ → ℝᴺ
abstract type AbstractMapping{N, D} end

# Evaluation of the actual value of the mapping
function point(m::AbstractMapping{D}, ξ::SVector{D}) where {D} end
function (m::AbstractMapping{N, D})(ξ::SVector{D})::SMatrix{N} where {N, D}
    return point(m, ξ)
end

# Jacobian
# This is the dumb default implementation that uses automatic differentiation.
# You should implement the actual analytical jacobian if you care about performacne!
function jacobian(m::AbstractMapping{N, D}, ξ)::SMatrix{N, D} where {N, D}
    J = ForwardDiff.jacobian(ξ → m(ξ), ξ)
    return SMatrix{N, D}
end

# Gram matrix
gram(J::AbstractMatrix)::AbstractMatrix = J' * J
function gram(m::AbstractMapping{N, D}, ξ)::AbstractMatrix where {N, D}
    return gram(jacobian(m, ξ))
end

# measure for integration
measure(J::AbstractMatrix) = T(sqrt(det(gram(J))))
function measure(m::AbstractMapping{N, D}, ξ)::AbstractMatrix where {N, D}
    return measure(gram(jacobian(m, ξ)))
end


# Composition
#
# Compose a mapping ℝᴰ → ℝᴹ with a mapping ℝᴹ → ℝᴺ
# Gives overall a mapping ℝᴰ → ℝᴺ
#
# Based on Julia ComposedFunction implementation: https://github.com/JuliaLang/julia/blob/d1c37793dd2ab0de6bca636e1d7f2ceb43150a9c/base/operators.jl#L1069-L1099
struct ComposedMapping{N, M, D, F <: AbstractMapping{N, M}, G <: AbstractMapping{M, D}} <: AbstractMapping{N, D}
    outer::F
    inner::G
end

function compose(f::AbstractMapping{N, M}, g::AbstractMapping{M, D}) where {N, M, D}
    return ComposedMapping{N, M, D, typeof(F), typeof(G)}(f, g)
end

Base.:∘(f::AbstractMapping, g::AbstractMapping) = compose(f, g)


# =====================================================================================
# AbstractElemet
#
# Mapping of unit cell to actual geometry
abstract type AbstractElement{N, D, T} <: AbstractMapping{N, D} end

# AbstractKernel
#
# A kernel defines the interaction between two points in space, e.g. the Coulomb interaction for electrostatics.
# Since this package has a focus on 2D BEM, even for electrostatics (or other possible future things), the kernels will look a bit funny
#   because we're going to be integrating out various degrees of freedom.
#
# Abstract interface for kernels:
#
# - patch(pt₁, pt₂) → Real (evaluate the kernel between two points)
#       Just evaluate the formula between two points.
#       WARNING: this will almost certainly be infinite if pt₁ = pt₂.
#           You probably want to avoid that happening in your code (singular quadrature), or handle it explicitly (limiting value for point collocation)
abstract type AbstractKernel{S <: AbstractSpace, T <: Real} end
(kernel::AbstractKernel{S, T})(pt₁::SVector{2, T}, pt₂::SVector{2, T}) where {S, T} = error("Not implemented for $(typeof(kernel))")
# TODO: is there a nice general way to ask for diagonal/adjacent element/far separated element blocks?
