using PyCall
using Conda
packages = Conda._installed_packages()
if !("rmg" in packages) && !("rmgmolecule" in packages) && (PyCall.pyversion.major != 3 || PyCall.pyversion.minor != 7)
    const Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    Conda.rm("mamba")
    Conda.add("conda=4")
    Conda.add("mamba")
    Conda.update()
    Conda.add("python=3.7")
    try
        Conda.rm("numpy") #get around MKL problem
    catch e 
    end
    Conda.add("nomkl")
    Conda.add("numpy")
    Conda.add_channel("mjohnson541")
    Conda.add_channel("hwpang")
    Conda.add("rmgmolecule")
    Pkg.build("PyCall")
end