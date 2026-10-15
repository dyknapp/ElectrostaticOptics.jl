# 2026 10 15, Daniel Knapp, Amsterdam : initial version

# Universal filleting between two curves using a circle.
# The key idea is that the center of the circle should lie along the normal
#   of the trimmed edge of both curves.
# We need to search for fillet radius ±R

# TODO: handle no-fillet case
# TODO: handling for no fillet found case

const ITS_MAX = 16

# Assume curve A has the shared vertex at ξ=1 and curve B at ξ=0
function solve_fillet(fillet_radius::Real, curves::Tuple{AbstractCurve2D, AbstractCurve2D})
    trims = SVector(0.9, 0.1) # Ansatz to slightly favor the correct orientation
    candidate_trims = SVector(0.999, 0.001)
    cent = SVector(0.0, 0.0)

    # Newton's method
    for it in 1:ITS_MAX # cut off if too many iterations
        # fcc = "fillet center candidate", p/m = plus/minus, curve 1/2
        fcc_p_1(ξ) = curves[1](ξ[1]) + fillet_radius * normal(curves[1], ξ[1])
        fcc_m_1(ξ) = curves[1](ξ[1]) - fillet_radius * normal(curves[1], ξ[1])
        fcc_p_2(ξ) = curves[2](ξ[2]) + fillet_radius * normal(curves[2], ξ[2])
        fcc_m_2(ξ) = curves[2](ξ[2]) - fillet_radius * normal(curves[2], ξ[2])

        # Which side should we be looking at?
        # (Whichever one is closer to being a fillet)
        if norm(fcc_p_1(trims) - fcc_p_2(trims)) < norm(fcc_m_1(trims) - fcc_m_2(trims))
            residual = ξ -> (fcc_p_1(ξ) - fcc_p_2(ξ))
            center = ξ -> (fcc_p_1(ξ) + fcc_p_2(ξ)) / 2
        else
            residual = ξ -> (fcc_m_1(ξ) - fcc_m_2(ξ))
            center = ξ -> (fcc_m_1(ξ) + fcc_m_2(ξ)) / 2
        end

        println("Iteration ", it, " ", residual(trims))

        # Newton iteration step, damped
        γ = 1.0
        J = ForwardDiff.jacobian(residual, trims)
        for _ in 1:10
            candidate_trims = trims - γ * (J \ residual(trims))
            println(γ, ", ", candidate_trims)
            if all([(0 ≤ val(ξ) ≤ 1) for ξ in candidate_trims])
                break
            end
            γ /= 2.0
        end
        trims = candidate_trims

        # De Morgan wrt the damping loop
        if any([!(0 ≤ val(ξ) ≤ 1) for ξ in trims])
            println(trims)
            error("Invalid ξ detected.  Fillet solve terminated.")
            # Converged to an invalid solution
        end
        if it == ITS_MAX
            error("Fillet solve did not converge by iteration ITS_MAX.")
        end
        if norm(residual(trims)) < 1.0e-12
            cent = center(trims)
            return trims, cent
        end
    end
    error("Fillet solver should never reach this error!!")
end

function solve_fillet_arc(fillet_radius::Real, curves::Tuple{AbstractCurve2D, AbstractCurve2D})
    trims, center = solve_fillet(fillet_radius, curves)

    trim_pt_A = curves[1](trims[1])
    trim_pt_B = curves[2](trims[2])

    vec_A = trim_pt_A - center
    vec_B = trim_pt_B - center

    θ₁ = atan(vec_A[2], vec_A[1])
    θ₂ = atan(vec_B[2], vec_B[1])

    τ_A = tangent(curves[1], trims[1])
    # 2D cross-product gives clockwise vs counterclockwise
    dir = sign(vec_A[1] * τ_A[2] - vec_A[2] * τ_A[1])

    Δθ = dir * mod(dir * (θ₂ - θ₁), 2π)
    θ₂ = θ₁ + Δθ

    fillet_arc = CircArc(center, fillet_radius, θ₁, θ₂)

    # Check that the fillet lines up nicely
    @assert(norm(tangent(curves[1], trims[1]) - tangent(fillet_arc, 0.0)) < 1.0e-11)

    return trims, fillet_arc
end
