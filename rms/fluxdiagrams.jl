using PyCall
using SparseArrays
using Images
@pyimport rmgpy.molecule as molecule
@pyimport pydot

function draw(spc::Species,path::String=".")
    name = spc.name
    fname = string(name,".png")
    if !(path in readdir("."))
        mkdir(path)
    else
        if fname in readdir(path)
            return
        end
    end
    if spc.inchi != ""
        mol = molecule.Molecule()[:fromInChI](spc.inchi)
    elseif spc.smiles != ""
        mol = molecule.Molecule()[:fromSMILES](spc.smiles)
    else
        throw(error("no smiles or inchi for molecule $name"))
    end
    mol[:draw](joinpath(path,fname))
end

function drawspecies(phase::T) where {T<:AbstractPhase}
    for spc in phase.species
        draw(spc,"species")
    end
end
