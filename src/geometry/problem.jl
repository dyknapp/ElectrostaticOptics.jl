# 2026 10 18, Daniel Knapp, Amsterdam : initial version


struct BEMproblem{B, K, F}
    # The geometry + topology info
    mesh::Mesh

    # # The charge density basis function(s') definition
    # basis::B

    # # The kernel for the BEM probelm
    # kernel::K

    # The kernel for this BEM problem
    kernel::K

    # panel_dofs is a mapping from geometry panel index
    #   to a vector of global indices for the panel's corresponding DOFs.
    # By "DOF" (degree of freedom) we mean a single basis function's coefficient
    # In general, geometry panels can share DOFs because
    #   we (sometimes) want the charge density to vary smoothly from panel to panel.
    #   (it wouldn't make sense to forcefully cut off the basis functions between panels)
    panel_dofs::Vector{Vector{Int}}

    # The inverse mapping from DOF to geometry.
    # So that we can quickly retrieve the source geometry for a specific DOF
    dof_panels::Vector{Vector{Int}}

    # TODO: reassess in the future if we need both mappings

    # The global array actually holding the specific DOFs.
    # The DOFs store the specific coordinates on their panels.
    dofs::Vector{F}

    # number of dofs.  N for NxN BEM matrix
    n_dofs::Int
end

# Assign a value to an index in a vector, padding with zeros if it doesn't already exist
function setindex_grow!(v::Vector{T}, val, i::Int) where {T}
    old_len = length(v)
    if i > old_len
        resize!(v, i)
        # Zero-fill the newly created gap between the old end and i - 1
        v[(old_len + 1):(i - 1)] .= zero(T)
    end
    v[i] = val
    return v
end

function pushindex_grow!(v::Vector{Vector{T}}, val, i::Int) where {T}
    old_len = length(v)
    if i > old_len
        resize!(v, i)
        # Fill the new gaps with empty inner vectors
        for j in (old_len + 1):i
            v[j] = T[]
        end
    end
    push!(v[i], val)
    return v
end

function assign_dofs(mesh::Mesh, dofs_per_panel::Int)
    g = mesh.topology
    for v in vertices(g)
        # Check that each node has one incoming and one outgoing edge
        @assert ((indegree(g, v) == 1) && (outdegree(g, v) == 1))
    end

    # Go through the closed loops
    dof_counter = 1 # global DOF counter
    panel_dofs = Vector{Vector{Int}}()

    cycles = simplecycles(g)
    for cycle in cycles
        for panel in cycle
            for _ in 1:dofs_per_panel
                # Assign DOF
                pushindex_grow!(panel_dofs, dof_counter, panel)
                dof_counter += 1
            end
        end
    end

    return panel_dofs, dof_counter
end

# TODO: something like this to check the Mesh geometry/topology
function check_mapping_agreement(p::BEMproblem)::Bool
    (; panel_dofs, dof_panels, dofs, n_dofs) = p

    if length(dofs) != n_dofs || length(dof_panels) != n_dofs
        return false
    end

    n_panels = length(panel_dofs)

    # Forward check: every (dof -> panel) link exists in panel_dofs
    for (dof, panels) in enumerate(dof_panels)
        # Surjection check: every DOF must belong to at least one panel
        # Sanity checks, map can't be empty or have repeat panels
        if isempty(panels)
            return false
        end
        if !allunique(panels)
            return false
        end
        # check that the mapping checks out
        for panel in panels
            # are we in bounds?
            if !checkbounds(Bool, panel_dofs, panel)
                return false
            end
            # do the check
            if !(dof in panel_dofs[panel])
                return false
            end
        end
    end

    # Backward check: every (panel -> dof) link exists in dof_panels
    for (panel, p_dofs) in enumerate(panel_dofs)
        if isempty(p_dofs)
            return false
        end
        if !allunique(p_dofs)
            return false
        end

        for dof in p_dofs
            if !checkbounds(Bool, dof_panels, dof)
                return false
            end
            if !(panel in dof_panels[dof])
                return false
            end
        end
    end

    return true
end
