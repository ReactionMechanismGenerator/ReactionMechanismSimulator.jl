using Documenter, RMS

makedocs(
    format = :html,
    sitename = "ReactionMechanismSimulator.jl",
    modules = [RMS],
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
