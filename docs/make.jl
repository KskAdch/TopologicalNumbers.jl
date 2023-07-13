using TopologicalNumbers
using Documenter

DocMeta.setdocmeta!(TopologicalNumbers, :DocTestSetup, :(using TopologicalNumbers); recursive=true)

makedocs(;
    modules=[TopologicalNumbers],
    authors="Keisuke Adachi, Minoru Kanega",
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
    ],
)

deploydocs(;
    repo="github.com/KskAdch/TopologicalNumbers.jl",
    devbranch="main",
)
