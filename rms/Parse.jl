using Unitful
using YAML
using PyCall

@pyimport rdkit.Chem as Chem

include("Tools.jl")
include("Calculators.jl")
include("Species.jl")
include("Reaction.jl")

module Calc
    include("Calculators.jl")
end
module Spc
    include("Species.jl")
end
module Rxn
    include("Reaction.jl")
end

const unitsdict = Dict()
const allowedfcnlist = vcat(names(Calc),names(Spc),names(Rxn))

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
    elseif isa(q,Number) || isa(q,AbstractArray)
        if units == nothing
            return q
        else
            return upreferred(Quantity(q,units)).val
        end
    elseif isa(q,Dict)
        return q
    else
        tq = typeof(q)
        throw(error("unable to parse $q of type $tq"))
    end

end

function parsestring(s)
    """
    Detects if a given string is actually a quantity input using metaprogramming
    the string is parsed by julia's parser and then analyzed to check if it is in
    the format of a quantity input
    """
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
end

function fcndict2obj(d::T,ymlunitsdict::Q) where {T,Q<:Any}
    """
    constructs an object from a dictionary by recursively constructing
    the objects in its fields
    the dictionary entry for "type" denotes the object being constructed
    which must be in the allowedfcnlist
    """
    strfcn = d["type"]
    fcn = Symbol(strfcn)
    @assert fcn in allowedfcnlist
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
                @warn("assuming values for $strfcn,$key are in SI units")
                val = parseentry(val)
            end
        end
        kwexpr = Expr(:kw,Symbol(key),val)
        push!(kwexprs,kwexpr)
    end
    ex = Expr(:call,fcn,Expr(:parameters,kwexprs...))
    return eval(ex)
end

function readinput(fname::String)
    """
    parses a YAML input file into a dictionary containing
    partitions of Species and Reaction objects
    """
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
        rxn["products"] = [r for r in rxn["products"]]
        rxn["reactantinds"] = [r.index for r in rxn["reactants"]]
        rxn["productinds"] = [r.index for r in rxn["products"]]
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
