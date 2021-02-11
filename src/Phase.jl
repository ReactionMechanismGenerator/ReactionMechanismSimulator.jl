using Parameters
import Base: length

using SparseArrays

abstract type AbstractPhase end
export AbstractPhase

abstract type IdealPhase <: AbstractPhase end
export IdealPhase

struct EmptyPhase <: AbstractPhase end
export EmptyPhase

include("Calculators/Ratevec.jl")
include("Calculators/Thermovec.jl")
include("Reaction.jl")

@with_kw struct IdealGas{W<:Tuple,W2,W3} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{ElementaryReaction,1}
    spcdict::Dict{String,Int64}
    stoichmatrix::W2
    Nrp::Array{Float64,1}
    veckinetics::W
    veckineticsinds::Array{Int64,1}
    vecthermo::NASAvec
    otherreactions::Array{ElementaryReaction,1}
    electronchange::W3
    reversibility::Array{Bool,1}
    diffusionlimited::Bool = false
end
IdealGas(species,reactions; name="",diffusionlimited=false) = IdealGas(species=species,reactions=reactions,name=name,
diffusionlimited=diffusionlimited,spcdict=Dict([sp.name=>sp.index for sp in species]))

function IdealGas(species,reactions; name="",diffusionlimited=false)
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,products=rxn.products,
        productinds=rxn.productinds,kinetics=rxn.kinetics,radicalchange=rxn.radicalchange,reversible=rxn.reversible,pairs=rxn.pairs) for (i,rxn) in enumerate(rxns)]
    therm = getvecthermo(species)
    M,Nrp = getstoichmatrix(species,rxns)
    echangevec = getfield.(rxns,:electronchange)
    if all(echangevec .== 0)
        electronchange = nothing
    else 
        electronchange = convert(echangevec,Array{Float64,1})
    end
    reversibility = getfield.(rxns,:reversible)
    return IdealGas(species=species,reactions=rxns,name=name,
        spcdict=Dict([sp.name=>sp.index for sp in species]),stoichmatrix=M,Nrp=Nrp,veckinetics=vectuple, 
        veckineticsinds=posinds, vecthermo=therm, otherreactions=otherrxns, electronchange=electronchange, 
        reversibility=reversibility,diffusionlimited=diffusionlimited,)
end
export IdealGas

@with_kw struct IdealDiluteSolution{W<:Tuple,W2,W3} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{ElementaryReaction,1}
    solvent::Solvent
    stoichmatrix::W2
    Nrp::Array{Float64,1}
    veckinetics::W
    veckineticsinds::Array{Int64,1}
    vecthermo::NASAvec
    otherreactions::Array{ElementaryReaction,1}
    electronchange::W3
    spcdict::Dict{String,Int64}
    reversibility::Array{Bool,1}
    diffusionlimited::Bool = true
end
IdealDiluteSolution(species,reactions,solvent; name="",diffusionlimited=true) = IdealDiluteSolution(species=species,reactions=reactions,
solvent=solvent,name=name,diffusionlimited=diffusionlimited,spcdict=Dict([sp.name=>sp.index for sp in species]))
function IdealDiluteSolution(species,reactions,solvent; name="",diffusionlimited=true)
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,products=rxn.products,
        productinds=rxn.productinds,kinetics=rxn.kinetics,radicalchange=rxn.radicalchange,reversible=rxn.reversible,pairs=rxn.pairs) for (i,rxn) in enumerate(rxns)]
    therm = getvecthermo(species)
    M,Nrp = getstoichmatrix(species,rxns)
    echangevec = getfield.(rxns,:electronchange)
    if all(echangevec .== 0)
        electronchange = nothing
    else 
        electronchange = convert(echangevec,Array{Float64,1})
    end
    reversibility = getfield.(rxns,:reversible)
    return IdealDiluteSolution(species=species,reactions=rxns,solvent=solvent,name=name,
        spcdict=Dict([sp.name=>sp.index for sp in species]),stoichmatrix=M,Nrp=Nrp,veckinetics=vectuple,
        veckineticsinds=posinds,vecthermo=therm,otherreactions=otherrxns,electronchange=electronchange,
        reversibility=reversibility,diffusionlimited=diffusionlimited)
end
export IdealDiluteSolution

@with_kw struct IdealSurface{W<:Tuple,W2,W3,W4,W5} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{ElementaryReaction,1}
    stoichmatrix::W2
    Nrp::W5
    veckinetics::W
    veckineticsinds::Array{Int64,1}
    vecthermo::W4
    otherreactions::Array{ElementaryReaction,1}
    electronchange::W3
    spcdict::Dict{String,Int64}
    reversibility::Array{Bool,1}
    sitedensity::Float64
    diffusionlimited::Bool = false
end
function IdealSurface(species,reactions,sitedensity;name="",diffusionlimited=false)
    @assert diffusionlimited==false "diffusionlimited=true not supported for IdealSurface"
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,products=rxn.products,
        productinds=rxn.productinds,kinetics=rxn.kinetics,radicalchange=rxn.radicalchange,electronchange=rxn.electronchange,reversible=rxn.reversible,pairs=rxn.pairs) for (i,rxn) in enumerate(rxns)]
    therm = getvecthermo(species)
    M,Nrp = getstoichmatrix(species,rxns)
    echangevec = getfield.(rxns,:electronchange).*F
    if all(echangevec .== 0)
        electronchange = nothing
    else 
        electronchange = convert(typeof(Nrp),echangevec)
    end
    reversibility = getfield.(rxns,:reversible)
    return IdealSurface(species=species,reactions=rxns,name=name,
        spcdict=Dict([sp.name=>sp.index for sp in species]),stoichmatrix=M,Nrp=Nrp,veckinetics=vectuple,
        veckineticsinds=posinds,vecthermo=therm,otherreactions=otherrxns,electronchange=electronchange,
        reversibility=reversibility,sitedensity=sitedensity,diffusionlimited=diffusionlimited)
end
export IdealSurface

"""
construct the stochiometric matrix for the reactions and the reaction molecule # change
"""
function getstoichmatrix(spcs,rxns)
    M = spzeros(length(rxns),length(spcs))
    Nrp = zeros(length(rxns))
    for (i,rxn) in enumerate(rxns)
        Nrp[i] = Float64(length(rxn.productinds) - length(rxn.reactantinds))
        for ind in rxn.reactantinds
            M[i,ind] += 1.0
        end
        for ind in rxn.productinds
            M[i,ind] -= 1.0
        end
    end
    return M,Nrp
end

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
            typ = split(tp[1],".")[end] #this split needs done because in pyrms the type names come as ReactionMechanismSimulator.Arrhenius instead of Arrhenius in RMS proper
        else
            typ = split(tp,".")[end]
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
    if isa(thermo[1],NASA)
        typeassert.(thermo,NASA)
        return NASAvec([sp.thermo for sp in spcs])
    elseif isa(thermo[1],ConstantG)
        typeassert.(thermo,ConstantG)
        return ConstantGvec([th.G for th in thermo],thermo[1].T)
    else
        t = typeof(thermo[1])
        error("Thermo type $t unsupported!")
    end
end
export getvecthermo

function getC0(ph::X,T) where {X<:Union{IdealDiluteSolution,IdealGas}}
    return 1.0e5/(R*T)
end

function getC0(ph::X,T) where {X<:IdealSurface}
    return ph.sitedensity
end

export getC0

length(p::T) where {T<:AbstractPhase} = 1
export length

iterate(p::T) where {T<:AbstractPhase} = p
export iterate

Broadcast.broadcastable(p::T) where {T<:AbstractPhase} = Ref(p)
export broadcastable
