# 2026 10 15, Daniel Knapp, Amsterdam : initial version

const SNAP_TOL = 1.0e-9

# TODO: handling for different types of incorrect mesh initializations
struct Mesh
    panels::Vector{AbstractCurve2D}
    topology::SimpleDiGraph

    function Mesh(panels::Vector{AbstractCurve2D}, topology::SimpleDiGraph)
        @assert length(panels) == nv(topology) # Initial check, is topology generally consistent with panels?
        # Next, go through and actually check connectivity
        for edge in edges(topology)
            println("Verifying edge ", edge)
            curve_A = panels[src(edge)]
            curve_B = panels[dst(edge)]
            @assert(norm(curve_A(1.0) - curve_B(0.0)) ≤ SNAP_TOL)
        end
        return new(panels, topology)
    end
end

function empty_mesh()
    return Mesh(AbstractCurve2D[], SimpleDiGraph())
end

function concat_meshes(mesh1::Mesh, mesh2::Mesh)
    return Mesh(vcat(mesh1.panels, mesh2.panels), blockdiag(mesh1.topology, mesh2.topology))
end

function circle(o::SVector{2}, r::T)::Mesh where {T}
    # Split the counterclockwise boundary into two arcs.
    arcs = AbstractCurve2D[CircArc(o, r, 0.0, 1.0 * π), CircArc(o, r, 1.0 * π, 2.0 * π)]
    topo = SimpleDiGraph([false true; true false])
    return Mesh(arcs, topo)
end


function rectangle(center::SVector{2}, width::Real, height::Real, θ::Real)
    # Draw a counterclockwise rectangle, rotate, and then move it.
    rectsides = [
        LineSeg(SVector(-width / 2, height / 2), SVector(-width / 2, -height / 2)),
        LineSeg(SVector(-width / 2, -height / 2), SVector(width / 2, -height / 2)),
        LineSeg(SVector(width / 2, -height / 2), SVector(width / 2, height / 2)),
        LineSeg(SVector(width / 2, height / 2), SVector(-width / 2, height / 2)),
    ]
    rot_mat = SMatrix{2, 2}([cos(θ) sin(θ); -sin(θ) cos(θ)])
    final_sides = AbstractCurve2D[]
    for side in rectsides
        side = affine_transform(side, rot_mat)
        push!(final_sides, translation(side, center))
    end

    topology = SimpleDiGraph(4)
    add_edge!(topology, 1, 2)
    add_edge!(topology, 2, 3)
    add_edge!(topology, 3, 4)
    add_edge!(topology, 4, 1)

    return Mesh(final_sides, topology)
end

@recipe function f(M::Mesh)
    seriestype   --> :path
    arrow        --> true
    aspect_ratio --> :equal
    labels       --> false

    for panel in M.panels
        @series begin
            t = LinRange(0.0, 1.0, 128)
            pts = panel.(t)
            [pt[1] for pt in pts], [pt[2] for pt in pts]
        end
    end
end
