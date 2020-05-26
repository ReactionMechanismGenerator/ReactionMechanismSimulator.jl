using Parameters
import Base: length

abstract type AbstractPhase end
export AbstractPhase

abstract type IdealPhase <: AbstractPhase end
export IdealPhase

struct EmptyPhase <: AbstractPhase end
export EmptyPhase

@with_kw struct IdealGas{Q<:AbstractReaction} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{Q,1}
    spcdict::Dict{String,Int64}
    diffusionlimited::Bool = false
end
IdealGas(species,reactions; name="",diffusionlimited=false) = IdealGas(species=species,reactions=reactions,name=name,
diffusionlimited=diffusionlimited,spcdict=Dict([sp.name=>sp.index for sp in species]))
export IdealGas

@with_kw struct IdealDiluteSolution{Q<:AbstractReaction} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{Q,1}
    solvent::Solvent
    spcdict::Dict{String,Int64}
    diffusionlimited::Bool = true
end
IdealDiluteSolution(species,reactions,solvent; name="",diffusionlimited=true) = IdealDiluteSolution(species=species,reactions=reactions,
solvent=solvent,name=name,diffusionlimited=diffusionlimited,spcdict=Dict([sp.name=>sp.index for sp in species]))
export IdealDiluteSolution

@with_kw struct HomogeneousCatalyst{Q<:AbstractReaction} <: AbstractPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{Q,1}
    spcdict::Dict{String,Int64}
end
export HomogeneousCatalyst

"""
Split the reactions into groups with the same kinetics
"""
function splitreactionsbykinetics(rxns)
    tps = []
    rxnlists = []
    for (i,rxn) in enumerate(rxns)
        typ = getkineticstype(rxn.kinetics)
        tind = findfirst(x->x==typ,tps)
        if tind == nothing
            push!(rxnlists,[])
            push!(tps,typ)
            push!(rxnlists[end],rxn)
        else
            push!(rxnlists[tind],rxn)
        end
    end
    return (tps,rxnlists)
end
export splitreactionsbykinetics

"""
create vectorized kinetics calculators for the reactions
"""
function getveckinetics(rxns)
    tps,rxnlists = splitreactionsbykinetics(rxns)
    posinds = Array{Int64,1}()
    fs = []
    otherrxns = Array{ElementaryReaction,1}()
    vecinds = Array{Int64,1}()
    otherrxninds = Array{Int64,1}()
    for (i,tp) in enumerate(tps)
        if typeof(tp)<:Tuple
            typ = tp[1]
        else
            typ = tp
        end
        rinds = [findfirst(isequal(rxn),rxns) for rxn in rxnlists[i]]
        fcn = Symbol(typ * "vec")
        if !(fcn in allowedfcnlist) || occursin("Troe",typ) #no vectorized kinetics or for the moment stuff with efficiencies
            append!(otherrxns,rxnlists[i])
            append!(otherrxninds,rinds)
        else
            append!(vecinds,rinds)
            fcn = eval(fcn)
            x = fcn([x.kinetics for x in rxnlists[i]])
            push!(fs,x)
            if posinds == Array{Int64,1}()
                push!(posinds,length(rinds))
            else 
                push!(posinds,length(rinds)+posinds[end])
            end
        end
    end
    vectuple = tuple(fs...)
    return (vectuple,vecinds,otherrxns,otherrxninds,posinds)
end
export getveckinetics

"""
create vectorized thermo calculator for species
"""
function getvecthermo(spcs)
    thermo = [sp.thermo for sp in spcs]
    typeassert.(thermo,NASA)
    return NASAvec([sp.thermo for sp in spcs])
end
export getvecthermo

length(p::T) where {T<:AbstractPhase} = 1
export length

iterate(p::T) where {T<:AbstractPhase} = p
export iterate

Broadcast.broadcastable(p::T) where {T<:AbstractPhase} = Ref(p)
export broadcastable
