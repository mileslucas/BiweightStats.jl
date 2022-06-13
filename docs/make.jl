using Biweight
using Documenter

setup = quote
    using Biweight
    using StableRNGs
    rng = StableRNG(1123)
end

DocMeta.setdocmeta!(Biweight, :DocTestSetup, setup; recursive=true)

makedocs(;
    modules=[Biweight],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/mileslucas/Biweight.jl/blob/{commit}{path}#{line}",
    sitename="Biweight.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mileslucas.github.io/Biweight.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mileslucas/Biweight.jl",
    devbranch="main",
)
