import Logging
Logging.disable_logging(Logging.Warn)

using PythonCall
using CondaPkg
const Chem = PythonCall.pynew()
const molecule = PythonCall.pynew()
const fragment = PythonCall.pynew()
const pydot = PythonCall.pynew()

packages = keys(CondaPkg.current_packages())

if !("rmg" in packages) && !("rmgmolecule" in packages)
    @info "missing rmg and rmgmolecule installing rmgmolecule..."
    if !(v"3.7" <= PythonCall.C.python_version() && PythonCall.C.python_version() <= v"3.9")
        @info "python version was not in 3.7-3.9 changing python version"
        CondaPkg.add("python"; version="3.9")
    end
    CondaPkg.add("rmgmolecule"; version=">=0.3.0", channel="mjohnson541")
    CondaPkg.add("matplotlib", channel="conda-forge")
    CondaPkg.add("rdkit", channel="conda-forge")
    CondaPkg.add("pydot", channel="conda-forge")

    Pkgc = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    Pkgc.build("PythonCall")
end

PythonCall.pycopy!(Chem, pyimport("rdkit.Chem"))
try
    PythonCall.pycopy!(molecule, pyimport("rmgpy.molecule"))
    PythonCall.pycopy!(fragment, pyimport("rmgpy.molecule.fragment"))
catch e
    PythonCall.pycopy!(molecule, pyimport("molecule.molecule"))
    PythonCall.pycopy!(fragment, pyimport("molecule.molecule.fragment"))
end
PythonCall.pycopy!(pydot, pyimport("pydot"))

include("Constants.jl")
include("Tools.jl")
include("Calculators/RateUncertainty.jl")
include("Calculators/ThermoUncertainty.jl")
include("Calculators/Thermo.jl")
include("Calculators/Diffusion.jl")
include("Calculators/Rate.jl")
include("Calculators/Viscosity.jl")
include("Calculators/kLAkH.jl")
include("Species.jl")
include("Solvent.jl")
include("Reaction.jl")
include("Phase.jl")
include("PhaseState.jl")
include("Interface.jl")
include("Domain.jl")
include("Parse.jl")
include("ModelReduction.jl")
include("Reactor.jl")
include("ThreadedSensitivities.jl")
include("Simulation.jl")
include("TransitorySensitivities.jl")
include("AutomaticMechanismAnalysis.jl")
include("Debugging.jl")
include("EdgeAnalysis.jl")
include("Plotting.jl")
include("fluxdiagrams.jl")
