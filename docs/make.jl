using LambertsProblem
using Documenter

DocMeta.setdocmeta!(LambertsProblem, :DocTestSetup, :(using LambertsProblem); recursive=true)

makedocs(;
    modules=[LambertsProblem],
    authors="Burton Yale <burtonyale@gmail.com>",
    sitename="LambertsProblem.jl",
    format=Documenter.HTML(;
        canonical="https://burtony3.github.io/LambertsProblem.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/burtony3/LambertsProblem.jl",
    devbranch="main",
)
