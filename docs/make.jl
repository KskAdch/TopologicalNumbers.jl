using TopologicalNumbers
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "docs.bib"); style=:numeric)

DocMeta.setdocmeta!(
    TopologicalNumbers, :DocTestSetup, :(using TopologicalNumbers); recursive=true
)

makedocs(;
    plugins=[bib],
    modules=[TopologicalNumbers],
    authors="Keisuke Adachi <18s2002x@gmail.com>, Minoru Kanega <phys_chibaraki@yahoo.co.jp>",
    repo="https://github.com/KskAdch/TopologicalNumbers.jl/blob/{commit}{path}#{line}",
    sitename="TopologicalNumbers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://KskAdch.github.io/TopologicalNumbers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "One-dimensional case" =>
                ["SSH model" => "1D/SSH.md", "Kitaev Chain" => "1D/Kitaev-Chain.md"],
            "Two-dimensional case (Chern)" => [
                "Square Lattice w/ Flux" => "2D/flux.md",
                "Haldane model" => "2D/Haldane.md",
                "Kitaev honeycomb model" => "2D/Kitaev-Honeycomb.md",
            ],
            "Two-dimensional case (Z2)" => [
                "Thouless pumping" => "2D/Thouless.md",
                "Kane-Mele model" => "2D/Kane-Mele.md",
                "BHZ model" => "2D/BHZ.md",
            ],
            "Three-dimensional case (Z)" => ["Weyl semimetal" => "3D/Weyl.md"],
            "Four-dimensional case (Z)" =>
                ["Lattice Dirac model" => "4D/LatticeDirac.md"],
        ],
        "Library" => ["Public" => "lib/public.md", "Internal" => "lib/internal.md"],
        "References" => "references.md",
    ],
)

deploydocs(; repo="github.com/KskAdch/TopologicalNumbers.jl", devbranch="main")
