function to_inti_panel(
        mesh::Mesh;
        subdivisions = ones(Int, length(mesh.panels)),
        quad_order::Int = 8,
    )
    inti_panels = Inti.EntityKey[]
    for (idx, panel) in enumerate(mesh.panels)
        push!(
            inti_panels,
            Inti.parametric_curve(
                panel, 0.0, 1.0;
                labels = ["panel_$idx"],
            )
        )
    end

    domain = Inti.Domain(inti_panels)

    # number of elements on each original geometry panel
    counts = Dict(
        inti_panels[i] => subdivisions[i]
            for i in eachindex(inti_panels)
    )
    imesh = Inti.meshgen(domain, counts)

    # Construct Gauss-Legendre points
    Q = Inti.Quadrature(imesh, Inti.GaussLegendre(quad_order))
    panel_nodes = [
        Inti.dom2qtags(Q, Inti.Domain(panel))
            for panel in inti_panels
    ]

    return (
        mesh = mesh,
        entities = inti_panels,
        domain = domain,
        imesh = imesh,
        Q = Q,
        panel_nodes = panel_nodes,
    )
end
