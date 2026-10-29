# 2026 10 29, Daniel Knapp, Amsterdam : initial version

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

    # kernel::TODO
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
function CollocationProblem(mesh::Mesh, space::FunctionSpace, method::Symbol)
    points = Vector{Vector{Float64}}()
    for idx in eachindex(mesh.panels)
        push!(points, assign_points(space.assigned_dofs[idx], method))
    end

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
    )
end
