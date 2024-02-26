using PyCall
if PyCall.pyversion.major != 3 || PyCall.pyversion.minor != 7
    using Conda
    const Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    try
        Conda.add("python=3.7")
    catch e
        Conda.add("conda=22")
        Conda.add("python=3.7")
    end
    
    try
        Conda.rm("numpy") #get around MKL problem
    catch e 
    end
    Conda.add("nomkl")
    Conda.add("numpy")
    try 
        pyimport("rmgpy")
    catch e
        Conda.add_channel("hwpang")
        Conda.add("rmgmolecule")
    end
    Pkg.build("PyCall")
end