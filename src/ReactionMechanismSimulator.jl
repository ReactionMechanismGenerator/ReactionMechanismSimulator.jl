module ReactionMechanismSimulator
using CondaPkg
using Logging
using PythonCall

has_rmgpy = try; PythonCall.pyimport("rmgpy"); true; catch e; false; end
has_rmgmolecule = try; PythonCall.pyimport("molecule"); true; catch e; false; end

if !has_rmgpy && !has_rmgmolecule
    @info "missing rmg and rmgmolecule installing rmgmolecule..."
    if !(v"3.7" <= PythonCall.C.python_version() && PythonCall.C.python_version() <= v"3.9")
        @info "python version was not in 3.7-3.9 changing python version"
        CondaPkg.add("python"; version=">=3.9",resolve=false)
    end
    
    CondaPkg.add("rmgmolecule"; version=">=0.3.0", channel="mjohnson541",resolve=false)
    CondaPkg.add("matplotlib", channel="conda-forge",resolve=false)
    CondaPkg.add("rdkit", channel="conda-forge",resolve=false)
    CondaPkg.add("pydot", channel="conda-forge",resolve=false,version=">=2.0")
    CondaPkg.resolve()
end

using PythonCall
const Chem = PythonCall.pynew()
const Desc = PythonCall.pynew()
const molecule = PythonCall.pynew()
const chemkin = PythonCall.pynew()
const species = PythonCall.pynew()
const reaction = PythonCall.pynew()
const nasa = PythonCall.pynew()
const wilhoit = PythonCall.pynew()
const arrhenius = PythonCall.pynew()
const falloff = PythonCall.pynew()
const chebyshev = PythonCall.pynew()
const solvation = PythonCall.pynew()
const fragment = PythonCall.pynew()
const pydot = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(Chem, pyimport("rdkit.Chem"))
    PythonCall.pycopy!(Desc, pyimport("rdkit.Chem.Descriptors"))
    try
        PythonCall.pycopy!(molecule, pyimport("rmgpy.molecule"))
        PythonCall.pycopy!(chemkin, pyimport("rmgpy.chemkin"))
        PythonCall.pycopy!(species, pyimport("rmgpy.species"))
        PythonCall.pycopy!(reaction, pyimport("rmgpy.reaction"))
        PythonCall.pycopy!(nasa, pyimport("rmgpy.thermo.nasa"))
        PythonCall.pycopy!(wilhoit, pyimport("rmgpy.thermo.wilhoit"))
        PythonCall.pycopy!(arrhenius, pyimport("rmgpy.kinetics.arrhenius"))
        PythonCall.pycopy!(falloff, pyimport("rmgpy.kinetics.falloff"))
        PythonCall.pycopy!(chebyshev, pyimport("rmgpy.kinetics.chebyshev"))
        PythonCall.pycopy!(solvation, pyimport("rmgpy.data.solvation"))
        PythonCall.pycopy!(fragment, pyimport("rmgpy.molecule.fragment"))
    catch e
        PythonCall.pycopy!(molecule, pyimport("molecule.molecule"))
        PythonCall.pycopy!(chemkin, pyimport("molecule.chemkin"))
        PythonCall.pycopy!(species, pyimport("molecule.species"))
        PythonCall.pycopy!(reaction, pyimport("molecule.reaction"))
        PythonCall.pycopy!(nasa, pyimport("molecule.thermo.nasa"))
        PythonCall.pycopy!(wilhoit, pyimport("molecule.thermo.wilhoit"))
        PythonCall.pycopy!(arrhenius, pyimport("molecule.kinetics.arrhenius"))
        PythonCall.pycopy!(falloff, pyimport("molecule.kinetics.falloff"))
        PythonCall.pycopy!(chebyshev, pyimport("molecule.kinetics.chebyshev"))
        PythonCall.pycopy!(solvation, pyimport("molecule.data.solvation"))
        PythonCall.pycopy!(fragment, pyimport("molecule.molecule.fragment"))
    end
    PythonCall.pycopy!(pydot, pyimport("pydot"))
end
include("Constants.jl")
include("Tools.jl")
include("Calculators/RateUncertainty.jl")
include("Calculators/ThermoUncertainty.jl")
include("Calculators/Thermo.jl")
include("Calculators/Diffusion.jl")
include("Calculators/Rate.jl")
include("Calculators/Viscosity.jl")
include("Calculators/Thermovec.jl")
include("Calculators/Ratevec.jl")
include("Calculators/kLAkH.jl")
include("Species.jl")
include("Solvent.jl")
include("Reaction.jl")
include("Phase.jl")
include("PhaseState.jl")
include("Interface.jl")
include("Domain.jl")
include("yml.jl")
include("Parse.jl")
include("ModelReduction.jl")
include("Reactor.jl")
include("ThreadedSensitivities.jl")
include("Simulation.jl")
include("TransitorySensitivities.jl")
include("AutomaticMechanismAnalysis.jl")
include("EdgeAnalysis.jl")
include("Debugging.jl")
include("Plotting.jl")
include("fluxdiagrams.jl")
end
