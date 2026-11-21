# 2026 11 16, Daniel Knapp, Amsterdam : initial version

# We don't need a function space; it is defined already by the collocation points (e.g. Gauss-Legendre)
struct NRow
    panel_idx::Int64
    ξ::Float64
    pt::SVector{2, Float64}
    weight::Float64
end

struct NystromProblem{T}
    mesh::Mesh
    # quadrature_type some kind of enum for picking different quadrature methods maybe TODO?
    row_defs::Vector{NRow} # panel number, ξ, panel(ξ), node weight

    kernel::AbstractKernel{T}
end

# For now, hard-code Gauss-Legendre and same number of pts.  Exapnd later if needed1
function NystromProblemProblem(mesh::Mesh, colloc_order::Int, kernel::AbstractKernel; offset_fac = 0.0)
    points = Vector{Vector{Float64}}(undef, length(mesh.panels))
    weights = Vector{Vector{Float64}}(undef, length(mesh.panels))
    for idx in eachindex(mesh.panels)
        x, w = gausslegendre(colloc_order)
        points[idx] = @. (x + 1) / 2
        weights[idx] = @. w / 2
    end

    assigned_pts = sum(length(ξs) for ξs in points)
    @assert assigned_pts == space.n_dofs

    row_defs = Vector{Row}()
    for (idx, (ξs, ws)) in enumerate(zip(points, weights))
        panel = mesh.panels[idx]
        for (ξ, w) in zip(ξs, ws)
            push!(row_defs, NRow(idx, ξ, panel(ξ), w))
        end
    end

    return NystromProblem(
        mesh,
        row_defs,
        kernel,
    )
end

# # Completely unoptimized and as-simple-as-possible reference
# function reference_matrix(problem::NystromProblem{T}) where {T <: Real}

# end
