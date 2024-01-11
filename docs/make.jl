using ParametricGrobnerBases
using Documenter

DocMeta.setdocmeta!(ParametricGrobnerBases, :DocTestSetup, :(using ParametricGrobnerBases); recursive=true)

makedocs(;
    modules=[ParametricGrobnerBases],
    authors="Andreas Bøgh Poulsen",
    sitename="ParametricGrobnerBases.jl",
    format=Documenter.HTML(;
        canonical="https://0708andreas.github.io/ParametricGrobnerBases.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/0708andreas/ParametricGrobnerBases.jl",
    devbranch="main",
)
