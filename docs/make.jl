using Documenter, KM3io

makedocs(;
    modules=[KM3io],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://git.km3net.de/jschumann/Neurthino.jl/blob/{commit}{path}#L{line}",
    sitename="Neurthino.jl",
    authors="Johannes Schumann, Tamas Gal",
    assets=String[],
)
