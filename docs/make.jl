using Documenter, Neurthino

makedocs(;
    modules = [Neurthino],
    authors = "Johannes Schumann, Tamas Gal",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/logo.ico"],
    ),
    pages=[
        "Introduction" => "index.md",
        "API" => "api.md",
    ],
    repo="https://github.com/KM3NeT/Neurthino.jl/blob/{commit}{path}#L{line}",
    sitename="Neurthino.jl",
)

deploydocs(;
    repo="github.com/KM3NeT/Neurthino.jl",
)
