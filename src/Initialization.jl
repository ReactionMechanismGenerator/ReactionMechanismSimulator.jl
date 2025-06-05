using PythonCall
using CondaPkg

function checkmolecule()
    has_rmgpy = try; PythonCall.pyimport("rmgpy"); true; catch e; false; end
    has_rmgmolecule = try; PythonCall.pyimport("molecule"); true; catch e; false; end
    return has_rmgpy,has_rmgmolecule
end

function installmolecule()
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

function loadsubmodules()
    include("Constants.jl")
    include("Tools.jl")
    include("Calculators/RateUncertainty.jl")
    include("Calculators/ThermoUncertainty.jl")
    include("Calculators/Thermo.jl")
    include("Calculators/Diffusion.jl")
    include("Calculators/Rate.jl")
    include("Calculators/ThermoCoverageDependence")
    include("Calculators/RateCoverageDependence.jl")
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