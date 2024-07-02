using CondaPkg
using PythonCall
packages = keys(CondaPkg.current_packages())
if !("rmg" in packages) && !("rmgmolecule" in packages)

    if !(v"3.7" <= PythonCall.C.python_version() && PythonCall.C.python_version() <= v"3.9")
        CondaPkg.add("python"; version=">=3.9")
    end
    CondaPkg.add("rmgmolecule"; version=">=0.3.0", channel="mjohnson541")
    CondaPkg.add("matplotlib", channel="conda-forge")
    CondaPkg.add("rdkit", channel="conda-forge")
    CondaPkg.add("pydot", channel="conda-forge")
end

const Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
Pkg.build("PythonCall")