# 2026 10 18, Daniel Knapp, Amsterdam : initial version


struct BEMproblem{B, F}
    # The geometry + topology info
    mesh::Mesh

    # The charge density basis function(s') definition
    basis::B

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
