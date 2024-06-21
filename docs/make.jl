using ParametricGroebnerBases
using Documenter

DocMeta.setdocmeta!(ParametricGroebnerBases, :DocTestSetup, :(using ParametricGroebnerBases); recursive=true)

makedocs(;
    modules=[ParametricGroebnerBases],
    authors="Andreas Bøgh Poulsen",
    sitename="ParametricGroebnerBases.jl",
    format=Documenter.HTML(;
        canonical="https://0708andreas.github.io/ParametricGroebnerBases.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
	"Reference" => "reference.md",
	"Example" => "example.md"
    ],
)

deploydocs(;
    repo="github.com/0708andreas/ParametricGroebnerBases.jl",
    devbranch="main",
)
