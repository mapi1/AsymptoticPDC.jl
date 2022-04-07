using AsymptoticPDC
using Documenter

makedocs(;
    modules=[AsymptoticPDC],
    authors="Marius Pille <marius.pille@hu-berlin.de> and contributors",
    repo="https://github.com/mapi1/AsymptoticPDC/blob/{commit}{path}#L{line}",
    sitename="Cardio.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mapi1.github.io/AsymptoticPDC.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "Examples.md"
    ],
)

deploydocs(;
    repo="github.com/mapi1/AsymptoticPDC",
)