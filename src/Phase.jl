using Parameters
import Base: length

using SparseArrays

abstract type AbstractPhase end
export AbstractPhase

abstract type IdealPhase <: AbstractPhase end
export IdealPhase

struct EmptyPhase <: AbstractPhase end
export EmptyPhase

abstract type FragmentBasedIdealPhase <: IdealPhase end
export FragmentBasedIdealPhase

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
    rxnarray::Array{Int64,2}
    veckinetics::W
    veckineticsinds::Array{Int64,1}
    vecthermo::NASAvec
    otherreactions::Array{ElementaryReaction,1}
    electronchange::W3
    reversibility::Array{Bool,1}
    forwardability::Array{Bool,1}
    diffusionlimited::Bool = false
end

function IdealGas(species,reactions; name="",diffusionlimited=false)
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,
        products=rxn.products,productinds=rxn.productinds,kinetics=rxn.kinetics,electronchange=rxn.electronchange,
        radicalchange=rxn.radicalchange,reversible=rxn.reversible,forwardable=rxn.forwardable,pairs=rxn.pairs,comment=rxn.comment) for (i,rxn) in enumerate(rxns)]
    therm = getvecthermo(species)
    rxnarray = getreactionindices(species,rxns)
    M,Nrp = getstoichmatrix(rxnarray,species)
    echangevec = getfield.(rxns,:electronchange)
    electronchange = convert(Array{Float64,1},echangevec)
    reversibility = getfield.(rxns,:reversible)
    forwardability = getfield.(rxns,:forwardable)
    return IdealGas(species=species,reactions=rxns,name=name,
        spcdict=Dict([sp.name=>i for (i,sp) in enumerate(species)]),stoichmatrix=M,Nrp=Nrp,rxnarray=rxnarray,veckinetics=vectuple, 
        veckineticsinds=posinds, vecthermo=therm, otherreactions=otherrxns, electronchange=electronchange, 
        reversibility=reversibility,forwardability=forwardability,diffusionlimited=diffusionlimited,)
end
export IdealGas

@with_kw struct IdealDiluteSolution{W<:Tuple,W2,W3} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{ElementaryReaction,1}
    solvent::Solvent
    stoichmatrix::W2
    Nrp::Array{Float64,1}
    rxnarray::Array{Int64,2}
    veckinetics::W
    veckineticsinds::Array{Int64,1}
    vecthermo::NASAvec
    otherreactions::Array{ElementaryReaction,1}
    electronchange::W3
    spcdict::Dict{String,Int64}
    reversibility::Array{Bool,1}
    forwardability::Array{Bool,1}
    diffusionlimited::Bool = true
end

function IdealDiluteSolution(species,reactions,solvent; name="",diffusionlimited=true)
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,
        products=rxn.products,productinds=rxn.productinds,kinetics=rxn.kinetics,electronchange=rxn.electronchange,
        radicalchange=rxn.radicalchange,reversible=rxn.reversible,forwardable=rxn.forwardable,pairs=rxn.pairs,comment=rxn.comment) for (i,rxn) in enumerate(rxns)]
    therm = getvecthermo(species)
    rxnarray = getreactionindices(species,rxns)
    M,Nrp = getstoichmatrix(rxnarray,species)
    echangevec = getfield.(rxns,:electronchange)
    electronchange = convert(Array{Float64,1},echangevec)
    reversibility = getfield.(rxns,:reversible)
    forwardability = getfield.(rxns,:forwardable)

    return IdealDiluteSolution(species=species,reactions=rxns,solvent=solvent,name=name,
        spcdict=Dict([sp.name=>i for (i,sp) in enumerate(species)]),stoichmatrix=M,Nrp=Nrp,rxnarray=rxnarray,veckinetics=vectuple,
        veckineticsinds=posinds,vecthermo=therm,otherreactions=otherrxns,electronchange=electronchange,
        reversibility=reversibility,forwardability=forwardability,diffusionlimited=diffusionlimited)
end
export IdealDiluteSolution

@with_kw struct IdealSurface{W<:Tuple,W2,W3,W4,W5} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{ElementaryReaction,1}
    stoichmatrix::W2
    Nrp::W5
    rxnarray::Array{Int64,2}
    veckinetics::W
    veckineticsinds::Array{Int64,1}
    vecthermo::W4
    otherreactions::Array{ElementaryReaction,1}
    electronchange::W3
    spcdict::Dict{String,Int64}
    reversibility::Array{Bool,1}
    forwardability::Array{Bool,1}
    sitedensity::Float64
    diffusionlimited::Bool = false
end
function IdealSurface(species,reactions,sitedensity;name="",diffusionlimited=false)
    @assert diffusionlimited==false "diffusionlimited=true not supported for IdealSurface"
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,
        products=rxn.products,productinds=rxn.productinds,kinetics=rxn.kinetics,electronchange=rxn.electronchange,
        radicalchange=rxn.radicalchange,reversible=rxn.reversible,forwardable=rxn.forwardable,pairs=rxn.pairs,comment=rxn.comment) for (i,rxn) in enumerate(rxns)]
    therm = getvecthermo(species)
    rxnarray = getreactionindices(species,rxns)
    M,Nrp = getstoichmatrix(rxnarray,species)
    echangevec = getfield.(rxns,:electronchange).*F
    electronchange = convert(Array{Float64,1},echangevec)
    reversibility = getfield.(rxns,:reversible)
    forwardability = getfield.(rxns,:forwardable)
    return IdealSurface(species=species,reactions=rxns,name=name,
        spcdict=Dict([sp.name=>i for (i,sp) in enumerate(species)]),stoichmatrix=M,Nrp=Nrp,rxnarray=rxnarray,veckinetics=vectuple,
        veckineticsinds=posinds,vecthermo=therm,otherreactions=otherrxns,electronchange=electronchange,
        reversibility=reversibility,forwardability=forwardability,sitedensity=sitedensity,diffusionlimited=diffusionlimited)
end
export IdealSurface

@with_kw struct FragmentBasedIdealFilm{W<:Tuple,W2,W3} <: FragmentBasedIdealPhase
    name::String = ""
    species::Array{Species,1}
    fragments::Array{Species,1}
    fragmentintermediates::Array{Species,1}
    reactions::Array{ElementaryReaction,1}
    stoichmatrix::W2
    Nrp::Array{Float64,1}
    rxnarray::Array{Int64,2}
    veckinetics::W
    veckineticsinds::Array{Int64,1}
    vecthermo::NASAvec
    otherreactions::Array{ElementaryReaction,1}
    electronchange::W3
    spcdict::Dict{String,Int64}
    reversibility::Array{Bool,1}
    forwardability::Array{Bool,1}
    diffusionlimited::Bool = false
end

function FragmentBasedIdealFilm(species, reactions; name="", diffusionlimited=false)
    @assert diffusionlimited==false "diffusionlimited=true not supported for FragmentBasedIdealFilm"
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,
        products=rxn.products,productinds=rxn.productinds,kinetics=rxn.kinetics,electronchange=rxn.electronchange,
        radicalchange=rxn.radicalchange,reversible=rxn.reversible,forwardable=rxn.forwardable,pairs=rxn.pairs,
        fragmentbasedreactants=rxn.fragmentbasedreactants,fragmentbasedreactantinds=rxn.fragmentbasedreactantinds,
        fragmentbasedproducts=rxn.fragmentbasedproducts,fragmentbasedproductinds=rxn.fragmentbasedproductinds,
        comment=rxn.comment) for (i,rxn) in enumerate(rxns)]
    therm = getvecthermo(species)
    rxnarray = getreactionindices(species,rxns)
    M,Nrp = getstoichmatrix(rxnarray,species)
    echangevec = getfield.(rxns,:electronchange)
    if all(echangevec .== 0)
        electronchange = nothing
    else 
        electronchange = convert(echangevec,Array{Float64,1})
    end
    reversibility = getfield.(rxns,:reversible)
    forwardability = getfield.(rxns,:forwardable)

    fragments = [spc for spc in species if spc.isfragment]
    fragmentintermediates = [spc for spc in species if spc.isfragmentintermediate]

    return FragmentBasedIdealFilm(species=species,fragments=fragments,fragmentintermediates=fragmentintermediates,
        reactions=rxns,name=name,
        spcdict=Dict([sp.name=>i for (i,sp) in enumerate(species)]),stoichmatrix=M,Nrp=Nrp,rxnarray=rxnarray,veckinetics=vectuple,
        veckineticsinds=posinds,vecthermo=therm,otherreactions=otherrxns,electronchange=electronchange,
        reversibility=reversibility,forwardability=forwardability,diffusionlimited=diffusionlimited)
end
export FragmentBasedIdealFilm

function getphasespecies(phase::FragmentBasedIdealPhase)
    return phase.fragments
end

function getphasespecies(phase::AbstractPhase)
    return phase.species
end

"""
construct the stochiometric matrix for the reactions and the reaction molecule # change
"""
function getstoichmatrix(rxnarray,spcs)
    M = spzeros(size(rxnarray)[2],length(spcs))
    Nrp = zeros(size(rxnarray)[2])
    for i in 1:size(rxnarray)[2]
        n = 0.0
        for (j,ind) in enumerate(rxnarray[:,i])
            if ind == 0
                continue
            elseif j > 3
                M[i,ind] -= 1.0
                n += 1.0
            else
                M[i,ind] += 1.0
                n -= 1.0
            end
        end
        Nrp[i] = n
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

function getC0(ph::X,T) where {X<:Union{FragmentBasedIdealFilm,IdealDiluteSolution,IdealGas}}
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

function getreactionindices(spcs,rxns) where {Q<:AbstractPhase}
    arr = zeros(Int64,(8,length(rxns)))
    names = [spc.name for spc in spcs]
    for (i,rxn) in enumerate(rxns)
        inds = [findfirst(isequal(spc),spcs) for spc in rxn.reactants]
        for (j,spc) in enumerate(rxn.reactants)
            ind = findfirst(isequal(spc),spcs)
            arr[j,i] = ind
            rxn.reactantinds[j] = ind
        end
        for (j,spc) in enumerate(rxn.products)
            ind = findfirst(isequal(spc),spcs)
            arr[j+4,i] = ind
            rxn.productinds[j] = ind
        end
        if hasproperty(rxn.kinetics,:efficiencies) && length(rxn.kinetics.nameefficiencies) > 0
            while length(rxn.kinetics.efficiencies) > 0
                pop!(rxn.kinetics.efficiencies)
            end
            for (key,val) in rxn.kinetics.nameefficiencies
                ind = findfirst(isequal(key),names)
                if !(ind === nothing)
                    rxn.kinetics.efficiencies[ind] = val
                end
            end
        end
    end
    return arr
end
export getreactionindices
