using BiweightStats
using Documenter

setup = quote
    using BiweightStats
    using StableRNGs
    rng = StableRNG(1123)
end

DocMeta.setdocmeta!(BiweightStats, :DocTestSetup, setup; recursive=true)

makedocs(;
    modules=[BiweightStats],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/mileslucas/BiweightStats.jl/blob/{commit}{path}#{line}",
    sitename="BiweightStats.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mileslucas.github.io/BiweightStats.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mileslucas/BiweightStats.jl",
    devbranch="main",
)
