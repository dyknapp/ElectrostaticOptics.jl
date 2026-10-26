# 2026 10 26, Daniel Knapp, Amsterdam : initial version
#
# The idea here is to eventually allow global DOFs to couple to multiple local DOFs.
# For instance, the charge distribution (G¹ geometry) should be continuous (C⁰) across panels, that already adds two constraints.
# When the geometry is smoother (e.g. Gᴺ), the local basis should follow (as Cᴺ⁻¹).
# As a starting point though, just have an identity mapping, so discontinuous charge distribution.

# Reminder: N-dimensional space, D-dimensional boundary
# B is the number of basis functions
# TODO generalize: abstract type BoundaryPanel{N, D, B} end

# BoundaryPanel folds together the BEM DOFs and the geometry panels
# basis is stored individually to allow mixed polynomial orders in a larger BEM problem
# L is number of basis functions, G is number of global DOFs they create
struct BoundaryPanel{L}
    panel::AbstractCurve2D
    basis::AbstractPanelBasis{L} # L basis functions
end

struct PanelDofMap{T}
    global_idcs::Vector{Int}
    weights::Matrix{T}  # L rows × G columns
end

function local_coefficients(map::PanelDofMap, u)
    return map.weights * view(u, map.global_idcs)
end

struct FunctionSpace{T}
    bem_panels::Vector{BoundaryPanel}
    maps::Vector{PanelDofMap{T}}
    n_dofs::Int # number of total DOFs

    # These get inherited from Mesh.jl
    topology::SimpleDiGraph
end

# I think this is a first pass at this, and it should be superseded by a more general function
# For now, just the same polynomial basis functions on each panel
function build_discontinuous_function_space(mesh::Mesh, basis::AbstractPanelBasis{L}) where {L}
    g = mesh.topology
    for v in vertices(g)
        # Check that each node has one incoming and one outgoing edge
        @assert ((indegree(g, v) == 1) && (outdegree(g, v) == 1))
    end

    # Go through the closed loops
    dof_counter = 0 # global DOF counter
    maps = PanelDofMap{Float64}[]
    bps = BoundaryPanel[]

    cycles = simplecycles(g)
    for cycle in cycles
        for panel_idx in cycle
            map = PanelDofMap(
                collect((dof_counter + 1):(dof_counter + L)),
                Matrix{Float64}(I, L, L),
            )
            dof_counter += L
            push!(maps, map)
            println(map)

            bp = BoundaryPanel(mesh.panels[panel_idx], basis)
            push!(bps, bp)
            println(bp)
        end
    end

    return FunctionSpace(
        bps,
        maps,
        dof_counter,
        mesh.topology,
    )
end
