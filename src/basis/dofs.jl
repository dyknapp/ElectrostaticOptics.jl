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
    n_dofs::Int # number of total DOFs

    # These get inherited from Mesh.jl
    topology::SimpleDiGraph
end

@enum BasisConstraint begin
    DISCONTINUOUS = -2
    VALUE = -1
    C0 = 0
    C1 = 1
end

function power_vector(n::Int, x::Float64)
    return transpose((x) .^ collect(0:(n - 1)))
end

function derivative_matrix(n::Int)
    D = zeros(Float64, n, n)
    for idx in 2:n
        D[idx - 1, idx] = Float64(idx - 1)
    end

    return D
end

function basis_matrix(input_basis::PolynomialBasis{n})
    L = zeros(Float64, n, n)
    for (col, poly) in enumerate(input_basis.mappings)
        L[:, col] = poly.p.coeffs
    end

    return L
end

# This should be generalized, for instance, if rational bases are used in the future
function generate_constraint_rows(
        ξ::Float64,
        panel::AbstractCurve2D,
        constraint::BasisConstraint,
        input_basis::PolynomialBasis{L, M}
    )
    constraint_rows = Matrix{Float64}[]
    pv = power_vector(L, ξ)
    BM = basis_matrix(input_basis)

    if (constraint == C0) || (constraint == C1) || (constraint == VALUE)
        push!(constraint_rows, pv * BM)

        if constraint == C1
            D = derivative_matrix(L)
            speed = measure(panel, ξ) # For true C1 rather than G1 continuity, adjust for the arc lengths
            push!(constraint_rows, pv * (D * BM) ./ speed)
        end
    end

    return constraint_rows
end

function build_function_space(mesh::Mesh, input_basis::AbstractBasis{L}, constraints_map::Dict{Int, BasisConstraint}) where {L}
    g = mesh.topology
    for v in vertices(g)
        # Check that each node has one incoming and one outgoing edge
        @assert ((indegree(g, v) == 1) && (outdegree(g, v) == 1))
    end

    # Go through the closed loops
    dof_counter = 0 # global DOF counter
    bps = BoundaryPanel[]
    global_dofs = []

    for edge in edges(g)
        src_constraints = generate_constraint_rows(
            0.0,
            constraints_map[src(edge)],
            input_basis,
        )
        dst_constraints = generate_constraint_rows(
            1.0,
            constraints_map[dst(edge)],
            input_basis,
        )
        constraints = [src_constraints; dst_constraints]

        if length(constraints) > L
            throw(ErrorException("Specified basis is of insufficient dimension to satisfy constraints."))
        end

        remaining_dofs = L - length(constraints)
        ξ_additional = [i / (remaining_dofs + 1) for i in 1:remaining_dofs]

        if remaining_dofs > 0
            for ξ in ξ_additional
                append!(
                    constraints, generate_constraint_rows(
                        ξ,
                        VALUE,
                        input_basis,
                    )
                )
            end
        end

        # We should have matched the number of global DOFs to the dimension of this basis
        @assert length(constraints) == L

        # TODO: avoid redundancy by passing this to the constraint generator?
        BM = basis_matrix(input_basis)

        M = reduce(vcat, constraints)

        # Rebuild new DOF basis
        dof_basis_matrix = BM * inv(M)
        mappings = SVector([Polynomial(dof_basis_matrix[row, :] for row in 1:L)])
        dof_basis = PolynomialBasis(mappings)

        push!(
            bps,
            BoundaryPanel(
                panel,
                dof_basis,
            )
        )

        # avoid overcounting; for shared DOFs, only count the src side of each edge.
        global_dofs += length(src_constraints) + remaining_dofs
    end
    return
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
