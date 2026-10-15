# 2026 10 10, Daniel Knapp, Egmond aan Zee: Start writing out basic implementation ideas (collocation derivation and mesh)
# 2026 10 14, Daniel Knapp, Amsterdam     : Add Literate.jl

using ElectrostaticOptics
using Documenter
using DocumenterCitations

if haskey(ENV, "GITHUB_ACTIONS")
    ENV["JULIA_DEBUG"] = "Documenter"
end

deployconfig = Documenter.auto_detect_deploy_system()
Documenter.post_status(deployconfig; type = "pending", repo = "github.com/dyknapp/ElectrostaticOptics.jl.git")

using Literate
using Plots # to not capture precompilation output
using StaticArrays
using Graphs


Literate.markdown(joinpath(@__DIR__, "literate", "example.jl"), joinpath(@__DIR__, "src", "generated"))
Literate.markdown(joinpath(@__DIR__, "literate", "abstract_curve.jl"), joinpath(@__DIR__, "src", "generated"))
Literate.markdown(joinpath(@__DIR__, "literate", "mesh.jl"), joinpath(@__DIR__, "src", "generated"))

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style = :numeric,
)

DocMeta.setdocmeta!(ElectrostaticOptics, :DocTestSetup, :(using ElectrostaticOptics); recursive = true)

makedocs(;
    plugins = [bib],
    modules = [ElectrostaticOptics],
    authors = "Daniel Knapp <daniel.y.knapp@gmail.com> and contributors",
    sitename = "ElectrostaticOptics.jl",
    format = Documenter.HTML(;
        canonical = "https://dyknapp.github.io/ElectrostaticOptics.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Mesh" => "abstract_mesh.md",
        "Collocation" => "point_collocation.md",
        "Step 1: Collocation" => "constant_charge_collocation.md",
        "AbstractCurve2D" => "generated/abstract_curve.md",
        "Mesh" => "generated/mesh.md",
        "Example" => "generated/example.md",
    ],
)

deploydocs(;
    repo = "github.com/dyknapp/ElectrostaticOptics.jl",
    devbranch = "main",
)
