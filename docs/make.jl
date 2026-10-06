using ElectrostaticOptics
using Documenter
using DocumenterCitations

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
        "Collocation" => "point_collocation.md",
    ],
)

deploydocs(;
    repo = "github.com/dyknapp/ElectrostaticOptics.jl",
    devbranch = "main",
)
