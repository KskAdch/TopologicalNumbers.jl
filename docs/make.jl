using TopologicalNumbers
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "docs.bib"); style=:numeric)

DocMeta.setdocmeta!(TopologicalNumbers, :DocTestSetup, :(using TopologicalNumbers); recursive=true)

makedocs(bib;
    modules=[TopologicalNumbers],
    authors="Keisuke Adachi <18s2002x@gmail.com>, Minoru Kanega <phys_chibaraki@yahoo.co.jp>",
    repo="https://github.com/KskAdch/TopologicalNumbers.jl/blob/{commit}{path}#{line}",
    sitename="TopologicalNumbers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://KskAdch.github.io/TopologicalNumbers.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API" => "internal.md",
        "References" => "references.md"
    ]
)

deploydocs(;
    repo="github.com/KskAdch/TopologicalNumbers.jl",
    devbranch="main"
)
