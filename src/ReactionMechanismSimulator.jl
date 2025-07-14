module ReactionMechanismSimulator
using CondaPkg
using Logging
using PythonCall

include("Initialization.jl")

has_rmgpy,has_rmgmolecule = checkmolecule()

if !has_rmgpy && !has_rmgmolecule
    @info "missing rmg and rmgmolecule installing rmgmolecule..."
    installmolecule()
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
loadsubmodules()
end
