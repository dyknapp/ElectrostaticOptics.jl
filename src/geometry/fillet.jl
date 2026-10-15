# 2026 10 15, Daniel Knapp, Amsterdam : initial version

# Universal filleting between two curves using a circle.
# The key idea is that the center of the circle should lie along the normal
#   of the trimmed edge of both curves.
# We need to search for fillet radius ±R

# Assume curve A has the shared vertex at ξ=1 and curve B at ξ=0
function solve_fillet(fillet_radius::Real, curves::Tuple{AbstractCurve2D, AbstractCurve2D})
    trims = SVector(0.9, 0.1) # Ansatz to slightly favor the correct orientation
    cent = SVector(0.0, 0.0)

    # Newton's method
    for it in 1:16 # cut off if too many iterations
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

        # Newton iteration step
        J = ForwardDiff.jacobian(residual, trims)
        trims = trims - J \ residual(trims)

        if any([!(0 ≤ val(ξ) ≤ 1) for ξ in trims])
            println("Invalid ξ detected.  Fillet solve terminated.")
            break # Converged to an invalid solution
        end
        if norm(residual(trims)) < 1.0e-12
            cent = center(trims)
            println("Fillet solve converged, center ", cent)
            break
        end
    end

    return trims, cent
end
