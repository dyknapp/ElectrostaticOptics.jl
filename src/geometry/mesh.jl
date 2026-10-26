# 2026 10 15, Daniel Knapp, Amsterdam : initial version
# 2026 10 26, Daniel Knapp, Amsterdam : smoothness & filleting

const SNAP_TOL = 1.0e-9

@enum Smoothness begin
    DISCONTINUOUS = -1
    G0 = 0
    G1 = 1
    G2 = 2
    G∞ = -1000
end

const MAX_SMOOTHNESS_CHECK = 1

# TODO: handling for different types of incorrect mesh initializations
# TODO: far in the future: generalize to allow 3D too.
struct Mesh
    panels::Vector{AbstractCurve2D}
    topology::SimpleDiGraph
    smoothness_map::Dict{Tuple{Int, Int}, Smoothness}

    function Mesh(panels::Vector{AbstractCurve2D}, topology::SimpleDiGraph)
        @assert length(panels) == nv(topology) # Initial check, is topology generally consistent with panels?
        # Next, go through and actually check connectivity
        # TODO: flagging and healing discontinuities
        for edge in edges(topology)
            println("Verifying edge ", edge)
            curve_A = panels[src(edge)]
            curve_B = panels[dst(edge)]
            @assert(norm(curve_A(1.0) - curve_B(0.0)) ≤ SNAP_TOL)
        end
        # Go through and assign smoothness
        smoothness_map = Dict{Tuple{Int, Int}, Smoothness}()
        for edge in edges(topology)
            curve_A = panels[src(edge)] # vertex at ξ=1
            curve_B = panels[dst(edge)] # vertex at ξ=0
            smoothness = 0
            dcA, dcB = curve_A, curve_B
            while smoothness ≤ MAX_SMOOTHNESS_CHECK
                dcA = derivative_curve(dcA)
                dcB = derivative_curve(dcB)
                if !all((dcA(1.0) / norm(dcA(1.0))) ≈ (dcB(0.0) / norm(dcB(0.0))))
                    break
                else
                    smoothness += 1
                end
            end
            smoothness_map[(src(edge), dst(edge))] = Smoothness(smoothness)
            if smoothness_map[(src(edge), dst(edge))] == Smoothness(MAX_SMOOTHNESS_CHECK + 1)
                println("Assigning edge ", edge, " smoothness: ", smoothness_map[(src(edge), dst(edge))], " (at least)")
            else
                println("Assigning edge ", edge, " smoothness: ", smoothness_map[(src(edge), dst(edge))])
            end
        end
        return new(panels, topology, smoothness_map)
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

# helper functionf or inserting a fillet
# For vertices connected v1→v3, return a vertex v2 such that the graph is connected v1→v2→v3
function insert_vertex!(g::SimpleDiGraph, v1, v3)
    @assert has_edge(g, v1, v3) # we can't insert into a non-existent edge.
    @assert add_vertex!(g) # did we succesfully add a vertex?
    v2 = nv(g) # grab the new vertex we added

    # update the graph with the new vertex
    rem_edge!(g, v1, v3)
    add_edge!(g, v1, v2)
    add_edge!(g, v2, v3)
    return v2
end

# remember, vertices in our graphs are panels, edges are joints between panels
function smoothness(m::Mesh, e::Graphs.SimpleEdge{Int64})
    return m.smoothness_map[(src(e), dst(e))]
end


# TODO: handling for overlapping fillets.
function fillet_corners(m::Mesh, radius::Real)
    panels = copy(m.panels)
    topology = copy(m.topology)
    # for edge in edges(topology)
    #     println(edge)
    # end

    # Bounds in each ORIGINAL panel's parameter space
    n = length(m.panels)
    # keep track of the trimming bounds:
    lo, hi = zeros(n), ones(n)

    # Read only from m; modify only the output panels and topology
    for edge in edges(m.topology)
        if !(smoothness(m, edge) == G0)
            continue
        end

        v1, v3 = src(edge), dst(edge)
        (a, b), arc = solve_fillet_arc(radius, (m.panels[v1], m.panels[v3]))

        hi[v1] = a
        lo[v3] = b

        println("Inserting edge: ", edge)
        # Update topology with new fillet
        v2 = insert_vertex!(topology, v1, v3)
        # add the fillet panel
        push!(panels, arc)
        # check consistency
        @assert length(panels) == v2
    end

    # Combine endpoint cuts before trimming each panel
    for i in eachindex(m.panels)
        @assert (0 ≤ lo[i] < hi[i] ≤ 1)

        if lo[i] != 0 || hi[i] != 1
            panels[i] = trim(m.panels[i], (lo[i], hi[i]))
        end
    end

    return Mesh(panels, topology)
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
