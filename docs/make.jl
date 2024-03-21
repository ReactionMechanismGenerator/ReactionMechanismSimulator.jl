import Logging
Logging.disable_logging(Logging.Warn)

using CondaPkg
using PythonCall
packages = keys(CondaPkg.current_packages())
if !("rmg" in packages) && !("rmgmolecule" in packages)
    CondaPkg.add("rmgmolecule"; version=">=0.3.0", channel="mjohnson541")
end
if !(v"3.7" <= PythonCall.C.python_version() && PythonCall.C.python_version() <= v"3.9")
    CondaPkg.add("python"; version="3.9")
end
if !("matplotlib" in packages)
    CondaPkg.add("matplotlib")
end
if !("rdkit" in packages)
    CondaPkg.add("rdkit")
end
if !("pydot" in packages)
    CondaPkg.add("pydot")
end
const Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
Pkg.build("PythonCall")

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
