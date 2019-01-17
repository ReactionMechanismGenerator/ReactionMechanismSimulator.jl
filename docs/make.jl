using Documenter, ReactionMechanismSimulator

makedocs(
    format = :html,
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
    repo = "github.com/mjohnson541/ReactionMechanismSimulator.jl.git"
)
