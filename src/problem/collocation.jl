# 2026 10 29, Daniel Knapp, Amsterdam : initial version
# 2026 11 04, Daniel Knapp, Amsterdam : :ArcLength collocation ponit assignment

struct Row
    panel_idx::Int64
    ξ::Float64
    pt::SVector{2, Float64}
end

struct CollocationProblem{T}
    mesh::Mesh
    space::FunctionSpace{T}
    points::Vector{Vector{Float64}} # panels x DOFs
    row_defs::Vector{Row} # panel number, ξ, panel(ξ)

    kernel::AbstractKernel{T}
end

# TODO: set up caching of quadrature rules?  Revisit when writing quadrature logic.
# TODO: move to dofs.jl so that global DOF assignment can share the same logic?
function assign_points(n::Int, method::Symbol)::Vector{Float64}
    # In general, important to exclude endpoints (although consistent endpoint inclusion is fine)
    # Collocation points at ξ=0,1 will cause BEM matrix to be singular if they overlap
    if method == :Uniform
        return [Float64(i / (n + 1)) for i in 1:n]
    elseif method == :UniformWithEndpoint
        # includes ξ=0
        return [Float64((i - 1) / n) for i in 1:n]
    elseif method == :GaussLegendre
        x, _ = gausslegendre(n)
        # Map (-1, 1) to (0, 1) for ξ
        return @. (x + 1) / 2
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end

# Construct collocation problem, including collocation point assignment
# offset_fac is only relevant when method == :ArcLength
function CollocationProblem(mesh::Mesh, space::FunctionSpace, kernel::AbstractKernel, method::Symbol; offset_fac = 0.0)
    points = Vector{Vector{Float64}}(undef, length(mesh.panels))

    if method == :ArcLength
        for cc in connected_components(mesh.topology)
            # First, measure the total arc length and the number of collocation points needed
            al = 0.0
            pts = 0
            for panel_idx in cc
                al += arc_length(mesh.panels[panel_idx])
                pts += space.assigned_dofs[panel_idx]
            end
            pt_spacing = al / pts

            # Iterate through the curves and assign points with equal arc length spacing
            @assert 0 ≤ offset_fac < 1.0
            traversed_length = pt_spacing * offset_fac
            available_length = 0.0
            previous_panel_lengths = 0.0

            for panel_idx in cc
                this_panel_length = arc_length(mesh.panels[panel_idx])
                available_length += this_panel_length

                ξs_panel = Float64[]

                while traversed_length < available_length
                    current_panel_length = traversed_length - previous_panel_lengths
                    ξ_pt = length_ξ(mesh.panels[panel_idx], current_panel_length)
                    push!(ξs_panel, ξ_pt)

                    traversed_length += pt_spacing
                end

                points[panel_idx] = ξs_panel

                previous_panel_lengths += this_panel_length
            end
        end
    else
        for idx in eachindex(mesh.panels)
            points[idx] = assign_points(space.assigned_dofs[idx], method)
        end
    end

    assigned_pts = sum(length(ξs) for ξs in points)
    @assert assigned_pts == space.n_dofs

    row_defs = Vector{Row}()
    for (idx, ξs) in enumerate(points)
        panel = mesh.panels[idx]
        for ξ in ξs
            push!(row_defs, Row(idx, ξ, panel(ξ)))
        end
    end

    return CollocationProblem(
        mesh,
        space,
        points,
        row_defs,
        kernel,
    )
end


# Completely unoptimized and as-simple-as-possible reference
function reference_matrix(problem::CollocationProblem{T}) where {T <: Real}
    space = problem.space

    # BEM matrix
    N = problem.space.n_dofs
    mat = zeros(N, N)

    # blank vector to use as unit vector
    u = zeros(N)
    for io in 1:N # this is the surce global DOF
        # prep the unit vector representing global DOF
        fill!(u, 0)
        u[io] = 1

        for (ii, bem_panel) in enumerate(space.bem_panels)
            # Does the unit vector u have support in bem_panel??
            coeffs = local_coefficients(space.maps[ii], u)
            if iszero(coeffs)
                continue # No overlap!
            end
            # There is support!

            # Reconstruct the basis function on bem_panel for the global DOF
            density(ξ) = sum(
                coeffs[l] * bem_panel.basis(l, ξ) for l in eachindex(coeffs)
            )

            # Let's scatter into the BEM matrix:
            # bem_panel is the source, enumerated row_idx/row is the target collocation point's row
            for (row_idx, row) in enumerate(problem.row_defs)
                # This is the integrand:
                f(ξ) = problem.kernel(bem_panel.panel(ξ), row.pt) *
                    density(ξ) * measure(bem_panel.panel, ξ)

                # If the target collocation point is on a source panel, we have a divergence
                # (Basis functions can have support on two panels!)
                bounds = if (ii == row.panel_idx) && (0 < row.ξ < 1)
                    (0.0, row.ξ, 1.0)
                else
                    (0.0, 1.0)
                end

                # Too high error bounds will cause domain error.
                value, err = quadgk(f, bounds...; rtol = 1.0e-9, atol = 1.0e-9)
                mat[row_idx, io] += value
            end
        end
    end

    return mat
end
