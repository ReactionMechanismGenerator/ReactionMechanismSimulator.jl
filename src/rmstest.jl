import Logging
Logging.disable_logging(Logging.Warn)
ENV["JULIA_CONDAPKG_BACKEND"] = "MicroMamba"
using CondaPkg 

packages = keys(CondaPkg.current_packages())

if !("rmg" in packages) && !("rmgmolecule" in packages)
    @info "missing rmg and rmgmolecule installing rmgmolecule..."
    if "python" in packages
        py_version = VersionNumber(CondaPkg.current_packages()["python"][:version])
    else
        py_version = nothing
    end
    if py_version === nothing || !(v"3.7" <= py_version && py_version <= v"3.9")
        @info "python version was not in 3.7-3.9 changing python version"
        CondaPkg.add("python"; version="3.9",resolve=false)
        CondaPkg.add("rmgmolecule"; version=">=0.3.0", channel="mjohnson541",resolve=false)
        CondaPkg.add("matplotlib", channel="conda-forge",resolve=false)
        CondaPkg.add("rdkit", channel="conda-forge",resolve=false)
        CondaPkg.add("pydot", channel="conda-forge",resolve=false,version=">=2.0")
        CondaPkg.resolve()
    else
        CondaPkg.add("rmgmolecule"; version=">=0.3.0", channel="mjohnson541",resolve=false)
        CondaPkg.add("matplotlib", channel="conda-forge",resolve=false)
        CondaPkg.add("rdkit", channel="conda-forge",resolve=false)
        CondaPkg.add("pydot", channel="conda-forge",resolve=false,version=">=2.0")
        CondaPkg.resolve()
    end
end

using PythonCall

const Chem = PythonCall.pynew()
const molecule = PythonCall.pynew()
const fragment = PythonCall.pynew()
const pydot = PythonCall.pynew()

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
