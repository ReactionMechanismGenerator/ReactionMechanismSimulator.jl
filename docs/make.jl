using Documenter, ReactionMechanismSimulator

makedocs(
    format = Documenter.HTML(),
    sitename = "ReactionMechanismSimulator.jl",
    modules = [ReactionMechanismSimulator],
    pages = [
        "Home" => "index.md",
        "Input.md",
        "Simulating.md",
        "Analysis.md",
    ]
)

deploydocs(
    repo = "github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl.git",
    target = "build",
)
