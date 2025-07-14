import Logging
Logging.disable_logging(Logging.Warn)
using CondaPkg 
using PythonCall

include("Initialization.jl")

has_rmgpy,has_rmgmolecule = checkmolecule()

if !has_rmgpy && !has_rmgmolecule
    @info "missing rmg and rmgmolecule installing rmgmolecule..."
    installmolecule()
end

const Chem = PythonCall.pynew()
const molecule = PythonCall.pynew()
const fragment = PythonCall.pynew()
const pydot = PythonCall.pynew()

PythonCall.pycopy!(Chem, PythonCall.pyimport("rdkit.Chem"))
try
    PythonCall.pycopy!(molecule, PythonCall.pyimport("rmgpy.molecule"))
    PythonCall.pycopy!(fragment, PythonCall.pyimport("rmgpy.molecule.fragment"))
catch e
    PythonCall.pycopy!(molecule, PythonCall.pyimport("molecule.molecule"))
    PythonCall.pycopy!(fragment, PythonCall.pyimport("molecule.molecule.fragment"))
end
PythonCall.pycopy!(pydot, PythonCall.pyimport("pydot"))

loadsubmodules()
