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

# TODO: move to static arrays
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
    assigned_dofs::Vector{Int64} # For collocation
    # supported_dofs::Vector{Int64} TODO: for galerkin
    n_dofs::Int # number of total DOFs

    # These get inherited from Mesh.jl
    # TODO: I don't think FunctionSpace needs topology
    topology::SimpleDiGraph
end

@enum BasisConstraint begin
    UNCONSTRAINED = -2
    VALUE = -1
    C0 = 0
    C1 = 1
end

# The example in matrices_polynomials.jl uses the "coefficient space" representation of polynomials.
# This is ok, but it ends up incurring a lot of linear transformations back and forth.
# We can get around this by working in the "basis function space" (I'm totally making up words here).
#
# VALUE CONSTRAINT
# With basis functions bᵢ(ξ) the charge density σ(ξ) = ∑ lᵢ bᵢ(ξ) where lᵢ is the "basis function space" vector
# The summation is just the 1xL matrix bᵢ acting on lᵢ, giving a constraint row for the constraint matrix.
# We can get around the basis change back to "coefficient space" and functions like power_vector by setting up the
#   calculation directly in "basis function space" like this.
#
# SLOPE CONSTRAINT
# To constraint the slope, we just differentiate σ:
# dσ/dξ = ∑ lᵢ bᵢ'(ξ)
# So if we just set up a derivative row vector, we're good to go.

function generate_constraint_rows(
        ξ::Float64,
        panel::AbstractCurve2D,
        constraint::BasisConstraint,
        input_basis::AbstractPanelBasis{L},
    ) where {L}
    constraint_rows = Matrix{Float64}[]

    if (constraint == C0) || (constraint == C1) || (constraint == VALUE)
        # Row of basis values evaluated at ξ (See VALUE CONSTRAINT comment, above)
        val_row = reshape([input_basis.mappings[j](ξ) for j in 1:L], 1, L)
        push!(constraint_rows, val_row)

        # Nested ifs work because a joint is C0 if it is C1 (idk if this is technically correct, but you know what I mean)
        if (constraint == C1)
            # See SLOPE CONSTRAINT comment, above
            speed = measure(panel, ξ)
            # basis functions are mappings R¹→R¹ so the jacobian is literally just the derivative, hence only() to turn it into a float
            deriv_row = reshape(
                [only(jacobian(input_basis.mappings[j], ξ)) for j in 1:L],
                1,
                L,
            )
            deriv_row ./= speed # account for arc length speed difference
            push!(constraint_rows, deriv_row)
        end
    end

    return constraint_rows
end

# Now organized in two passes:
# PASS 1: assign shared global DOFs on joints
# PASS 2: Build PanelDofMap for each panel
function build_function_space(
        mesh::Mesh,
        input_basis::AbstractPanelBasis{L},
        constraints_map::Dict{Tuple{Int, Int}, BasisConstraint},
    ) where {L}
    g = mesh.topology

    for v in vertices(g)
        @assert (indegree(g, v) == 1) && (outdegree(g, v) == 1)
    end

    dof_counter = 0

    # PASS 1
    edge_dofs = Dict{Tuple{Int, Int}, Vector{Int}}()
    for edge in edges(g)
        joint = (src(edge), dst(edge))
        c = constraints_map[joint]

        if (c == C0) || (c == VALUE)
            # Only one constraint: add 1 Dof
            edge_dofs[joint] = [dof_counter + 1]
            dof_counter += 1
        elseif (c == C1)
            # two constraints
            edge_dofs[joint] = [dof_counter + 1, dof_counter + 2]
            dof_counter += 2
        else
            # UNCONSTRAINED is the only case that comes to mind for this case, at the moment
            edge_dofs[joint] = Int[]
        end
    end

    # PASS 2
    bps = BoundaryPanel[]
    maps = PanelDofMap{Float64}[]
    assigned_dofs = zeros(Int64, nv(g))

    for v in vertices(g)
        panel = mesh.panels[v]

        # The below only works because we asserted that each vertex has one edge in and one edge out at the beginning
        joint_in = (only(inneighbors(g, v)), v)
        joint_out = (v, only(outneighbors(g, v)))

        src_constraints = generate_constraint_rows(
            0.0,
            panel,
            constraints_map[joint_in],
            input_basis,
        )
        dst_constraints = generate_constraint_rows(
            1.0,
            panel,
            constraints_map[joint_out],
            input_basis,
        )

        n_endpoint_constraints = length(src_constraints) + length(dst_constraints)
        if n_endpoint_constraints > L
            throw(Base.ErrorException("Basis is overconstrained!"))
        end

        # Soak up the remaining DOFs.
        remaining_dofs = L - n_endpoint_constraints
        ξ_additional = [i / (remaining_dofs + 1) for i in 1:remaining_dofs]

        interior_constraints = Matrix{Float64}[]
        for ξ in ξ_additional
            append!(
                interior_constraints,
                generate_constraint_rows(ξ, panel, VALUE, input_basis),
            )
        end

        interior_dofs = collect((dof_counter + 1):(dof_counter + remaining_dofs))
        dof_counter += remaining_dofs
        assigned_dofs[v] = remaining_dofs + length(edge_dofs[joint_out])

        # Assemble the LxL constraint matrix
        # Reordered nicely to flow along the graph
        constraints = [src_constraints; interior_constraints; dst_constraints]
        @assert length(constraints) == L
        M = reduce(vcat, constraints)

        panel_global_idcs = [edge_dofs[joint_in]; interior_dofs; edge_dofs[joint_out]]
        @assert length(panel_global_idcs) == L

        push!(bps, BoundaryPanel(panel, input_basis))
        push!(maps, PanelDofMap(panel_global_idcs, inv(M)))
    end

    # Check that we assigned the right number of DOFs to panels
    @assert sum(assigned_dofs) == dof_counter

    return FunctionSpace(
        bps,
        maps,
        assigned_dofs,
        dof_counter,
        mesh.topology,
    )
end
