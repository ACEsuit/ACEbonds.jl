using ACEbonds
using Documenter

DocMeta.setdocmeta!(ACEbonds, :DocTestSetup, :(using ACEbonds); recursive=true)

makedocs(;
    modules=[ACEbonds],
    authors="Christoph Ortner <christophortner0@gmail.com> and contributors",
    repo="https://github.com/ACEsuit/ACEbonds.jl/blob/{commit}{path}#{line}",
    sitename="ACEbonds.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ACEsuit.github.io/ACEbonds.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ACEsuit/ACEbonds.jl",
    devbranch="main",
)
