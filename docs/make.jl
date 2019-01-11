using Documenter, ReactionMechanismSimulator

makedocs(
    format = :html,
    sitename = "ReactionMechanismSimulator.jl",
    modules = [ReactionMechanismSimulator],
    pages = [
        "Input.md",
        "Simulating.md",
        "Analysis.md",
    ]
)

deploydocs(
    repo = "github.com/mjohnson541/ReactionMechanismSimulator.jl.git",
    target = "build",
    julia  = "1.0",
    deps = nothing,
    make = nothing
)
