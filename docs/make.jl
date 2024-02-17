import Logging
Logging.disable_logging(Logging.Warn)

using Documenter, ReactionMechanismSimulator

makedocs(
    format = Documenter.HTML(),
    sitename = "ReactionMechanismSimulator.jl",
    modules = [ReactionMechanismSimulator],
    pages = [
        "Home" => "index.md",
        "Installation.md",
        "Input.md",
        "Simulating.md",
        "Analysis.md",
        "AutomaticMechanismAnalysis.md"
    ],
    warnonly=true,
)

deploydocs(
    repo = "github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl.git",
    target = "build",
)
