using CondaPkg
using PythonCall
packages = keys(CondaPkg.current_packages())
if !("rmg" in packages) && !("rmgmolecule" in packages) || !(v"3.7" <= PythonCall.C.python_version() && PythonCall.C.python_version() <= v"3.9")
    const Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    CondaPkg.add("python"; version="3.9")
    CondaPkg.add("matplotlib")
    CondaPkg.add("rdkit")
    CondaPkg.add("pydot")
    CondaPkg.add("rmgmolecule"; version=">=0.3.0", channel="mjohnson541")
    Pkg.build("PythonCall")
end