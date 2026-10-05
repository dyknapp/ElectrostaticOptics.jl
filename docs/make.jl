using ElectrostaticOptics
using Documenter

DocMeta.setdocmeta!(ElectrostaticOptics, :DocTestSetup, :(using ElectrostaticOptics); recursive = true)

makedocs(;
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
    ],
)

deploydocs(;
    repo = "github.com/dyknapp/ElectrostaticOptics.jl",
    devbranch = "main",
)
