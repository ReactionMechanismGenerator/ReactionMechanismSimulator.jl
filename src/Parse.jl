using Unitful
using YAML
using PyCall
using StaticArrays

module Calc
    include("Calculators/RateUncertainty.jl")
    include("Calculators/ThermoUncertainty.jl")
    include("Calculators/Thermo.jl")
    include("Calculators/Diffusion.jl")
    include("Calculators/Rate.jl")
    include("Calculators/Viscosity.jl")
end
module Spc
    include("Calculators/ThermoUncertainty.jl")
    include("Calculators/Thermo.jl")
    include("Calculators/Diffusion.jl")
    include("Species.jl")
end
module Rxn
    include("Calculators/RateUncertainty.jl")
    include("Calculators/ThermoUncertainty.jl")
    include("Calculators/Thermo.jl")
    include("Calculators/Diffusion.jl")
    include("Calculators/Rate.jl")
    include("Species.jl")
    include("Reaction.jl")
end

const unitsdict = Dict()
const elementdict = Dict([1=>"H",6=>"C",8=>"O",7=>"N",17=>"Cl",16=>"S",18=>"Ar",10=>"Ne",2=>"He",
        15=>"P",9=>"F",35=>"Br",53=>"I",289=>"Fl"])

const allowedfcnlist = vcat(names(Calc),names(Spc),names(Rxn))

"""
Zhao et al 2003
"""
const mcgowanvolumes = Dict(["H"=> upreferred(8.71u"mL/mol").val, "He"=> upreferred(6.75u"mL/mol").val,
        "C"=> upreferred(16.35u"mL/mol").val, "N"=> upreferred(14.39u"mL/mol").val,
         "O"=> upreferred(12.43u"mL/mol").val, "F"=> upreferred(10.47u"mL/mol").val,
         "Ne"=> upreferred(8.51u"mL/mol").val,"Si"=> upreferred(26.83u"mL/mol").val, "P"=> upreferred(24.87u"mL/mol").val,
        "S"=> upreferred(22.91u"mL/mol").val, "Cl"=> upreferred(20.95u"mL/mol").val,
         "Ar"=> upreferred(18.99u"mL/mol").val,"Br"=> upreferred(26.21u"mL/mol").val,])

function parseentry(q,units=nothing)
    """
    parses single YAML lines into strings, numbers and arrays of Numbers in SI units
    """
    if isa(q,String)
        q = parsestring(q)
    end

    if isa(q,String)
        return q
    elseif isa(q,Quantity)
        return upreferred(q).val
    elseif isa(q,Number)
        if units == nothing
            return q
        else
            return upreferred(Quantity(q,units)).val
        end
    elseif isa(q,AbstractArray)
        if units == nothing
            if typeof(q[1]) <: AbstractArray
                return convert(Array,hcat(q...)')
            else
                return q
            end
        else
            if typeof(q[1]) <: AbstractArray
                return convert(Array,hcat(upreferred(Quantity(q,units)).val...)')
            else
                return upreferred(Quantity(q,units)).val
            end
        end
    elseif isa(q,Dict)
        return q
    else
        tq = typeof(q)
        throw(error("unable to parse $q of type $tq"))
    end

end
export parseentry

function parsestring(s)
    """
    Detects if a given string is actually a quantity input using metaprogramming
    the string is parsed by julia's parser and then analyzed to check if it is in
    the format of a quantity input
    """
    try
        ex = Meta.parse(s)
        if !isa(ex,Symbol) && !isa(ex,String) && length(ex.args) == 3 && ex.args[1] == :* &&
            (isa(ex.args[2],Number) || (ex.args[2].head == :vec && all([isa(x,Number) for x in ex.args[2].args]))) &&
            ex.args[3].args[1] == Symbol("@u_str") &&
            isa(ex.args[3].args[2],LineNumberNode) && ex.args[3].args[2].line == 1 &&
            ex.args[3].args[2].file == :none && isa(ex.args[3].args[3],String)
            return eval(ex)
        else
            return s
        end
    catch
        return s
    end
end
export parsestring

function getatomdictfromrdkit(mol)
    """
    retrives the number of each type of atom and the number of bonds of an rdkit molecule
    """
    atmD = Dict{String,Int64}()
    for atm in mol.GetAtoms()
        v = elementdict[atm.GetAtomicNum()]
        if v in keys(atmD)
            atmD[v] += 1
        else
            atmD[v] = 1
        end
    end
    nbonds = length(mol.GetBonds())
    return atmD,nbonds
end
export getatomdictfromrdkit

getatomdictsmiles(smiles) = getatomdictfromrdkit(Chem.AddHs(Chem.MolFromSmiles(smiles)))
export getatomdictsmiles
getatomdictinchi(inchi) = getatomdictfromrdkit(Chem.AddHs(Chem.MolFromInchi(inchi)))
export getatomdictinchi

function getspeciesradius(atomdict::Dict{String,Int64},nbonds::Int64)
    """
    estimates the McGowan radius by calculating the McGowan Volume
    from the number of each type of atom and the number of bonds
    as described in Zhao et al. 2003
    """
    Vtot = 0.0
    for (atm,n) in atomdict
        Vtot += n*mcgowanvolumes[atm]
    end
    Vtot -= nbonds*upreferred(6.56u"mL/mol").val
    V = Vtot/Na
    r = (0.75/Base.pi*V)^(1.0/3.0)
    return r
end
export getspeciesradius

function fcndict2obj(d::T,ymlunitsdict::Q) where {T,Q<:Any}
    """
    constructs an object from a dictionary by recursively constructing
    the objects in its fields
    the dictionary entry for "type" denotes the object being constructed
    which must be in the allowedfcnlist
    """
    strfcn = d["type"]
    fcn = Symbol(strfcn)
    @assert fcn in allowedfcnlist fcn
    kwexprs = Array{Expr,1}()
    for (key,val) in d
        if key == "type"
            continue
        end
        if isa(val,AbstractDict) && "type" in keys(val)
            val = fcndict2obj(val,ymlunitsdict)
        elseif isa(val,AbstractArray) && typeof(val).parameters[1] <: AbstractDict
            val = [fcndict2obj(x,ymlunitsdict) for x in val]
        else
            if (strfcn,key) in keys(unitsdict)
                if unitsdict[(strfcn,key)] in keys(ymlunitsdict)
                    val = parseentry(val,ymlunitsdict[unitsdict[(strfcn,key)]])
                else
                    val = parseentry(val)
                end
            else
                #@warn("assuming values for $strfcn,$key are in SI units")
                val = parseentry(val)
            end
        end
        kwexpr = Expr(:kw,Symbol(key),val)
        push!(kwexprs,kwexpr)
    end
    ex = Expr(:call,fcn,Expr(:parameters,kwexprs...))
    return eval(ex)
end
export fcndict2obj

"""
examines the input file and parses as appropriate
for chemkin (.inp files) with or without species dictionaries
calls rmgpy to convert it to yml and then parses the yml file
for rms (.rms, .yml, .rmg) it parses it directly
"""
function readinput(fname::String;spcdict::String="",output::String="chem.rms")
    extension = split(fname,".")[end]
    if extension in ["rms","yml","rmg"]
        return readinputyml(fname)
    elseif extension in ["inp"]
        if spcdict == ""
            convertchemkin2yml(fname,output=output)
        else
            convertchemkin2yml(fname,spcdictpath=spcdict,output=output)
        end
        return readinputyml(output)
    else
        throw(error("file extension $extension not understood, use .inp for chemkin files or .rms for rms yaml files"))
    end
end
export readinput

"""
parses a YAML input file into a dictionary containing
partitions of Species and Reaction objects
"""
function readinputyml(fname::String)
    D = YAML.load(open(fname))
    outdict = Dict()

    #Units
    ymlunitsdict=D["Units"] #for now don't use

    #Solvents
    if "Solvents" in keys(D)
        s = D["Solvents"]
        for sol in s
            sol["type"] = "Solvent"
        end
        solvents = map(fcndict2obj,s)
        outDict["Solvents"] = solvents
    end

    #phases
    spcindex = 1
    spcdict = Dict()
    for p in D["Phases"]
        name = p["name"]
        spclist = Array{Species,1}()
        for d in p["Species"]
            d["index"] = spcindex
            spcindex += 1
            spcname = d["name"]
            #attempt to generate molecular information from rdkit if possible
            if !("atomnums" in keys(d)) || !("bondnum" in keys(d))
                if "smiles" in keys(d)
                    try
                        d["atomnums"],d["bondnum"] = getatomdictsmiles(d["smiles"])
                    catch
                        @warn("failed to generate molecular information from smiles for species $spcname")
                    end
                elseif "inchi" in keys(d)
                    try
                        d["atomnums"],d["bondnum"] = getatomdictinchi(d["inchi"])
                    catch
                        @warn("failed to generate molecular information from inchi for species $spcname")
                    end
                end
            end
            #if diffusion model is not present and we can calculate molecular radius and stokesdiffusivity
            if !("diffusion" in keys(d)) && "atomnums" in keys(d) && "bondnum" in keys(d)
                try
                    diffD = Dict()
                    diffD["type"] = "StokesDiffusivity"
                    diffD["r"] = getspeciesradius(d["atomnums"],d["bondnum"])
                    d["diffusion"] = diffD
                    if !("radius" in keys(d))
                        d["radius"] = diffD["r"]
                    end
                catch
                    @warn("failed to generate StokesDiffusivity model for species $spcname")
                end
            end

            spc = fcndict2obj(d,ymlunitsdict)
            push!(spclist,spc)
            spcdict[spc.name] = (spc,p["name"])
        end
        outdict[name] = Dict()
        outdict[name]["Species"] = spclist
        outdict[name]["Reactions"] = Array{ElementaryReaction,1}()
    end

    #reactions
    rxnindex = 1
    for rxn in D["Reactions"]
        phs = Array{String,1}()
        rxn["reactants"] = convert(Array{Any},rxn["reactants"])
        rxn["index"] = rxnindex
        rxnindex += 1
        for (i,item) in enumerate(rxn["reactants"])
            spc,ph = spcdict[item]
            push!(phs,ph)
            rxn["reactants"][i] = spc
        end
        rxn["products"] = convert(Array{Any},rxn["products"])
        for (i,item) in enumerate(rxn["products"])
            spc,ph = spcdict[item]
            push!(phs,ph)
            rxn["products"][i] = spc
        end
        rxn["reactants"] = [r for r in rxn["reactants"]]
        rxn["reactants"] = SVector(rxn["reactants"]...)
        rxn["products"] = [r for r in rxn["products"]]
        rxn["products"] = SVector(rxn["products"]...)
        rxn["reactantinds"] =  [r.index for r in rxn["reactants"]]
        rxn["reactantinds"] = SVector(rxn["reactantinds"]...)
        rxn["productinds"] =  [r.index for r in rxn["products"]]
        rxn["productinds"] = SVector(rxn["productinds"]...)
        if "efficiencies" in keys(rxn["kinetics"])
            rxn["kinetics"]["efficiencies"] = Dict{Int64,Float64}([spcdict[name][1].index=>val-1.0 for (name,val) in rxn["kinetics"]["efficiencies"]]) #in RMS we correct [M] rather than calculate it so we subtract 1
        end
        r = fcndict2obj(rxn,ymlunitsdict)
        unique!(phs)
        if length(phs) == 1
            push!(outdict[phs[1]]["Reactions"],r)
        else
            if Set(phs) in keys(outDict)
                push!(outdict[Set(phs)],r)
            else
                outdict[Set(phs)]=[r]
            end
        end
    end

    return outdict
end

export readinputyml
