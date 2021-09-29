using Parameters
using Statistics

@with_kw struct DebugSpecies
    name::String
    G298::Float64 = 0.0
    ratio::Float64 = 0.0
    dy::Float64 = 0.0
    tol::Float64 = 0.0
    index::Int64 = 0
end

@with_kw struct DebugReaction
    reactants::Array{DebugSpecies,1}
    products::Array{DebugSpecies,1}
    rxnstring::String
    tol::Float64 = 0.0
    kf::Float64 = 0.0
    krev::Float64 = 0.0
    rt::Float64 = 0.0
    ratio::Float64 = 0.0
    T::Float64 = 0.0
    P::Float64 = 0.0
    index::Int64 = 0
end

struct DebugMech
    spcs::Array{DebugSpecies,1}
    rxns::Array{DebugReaction,1}
end

function getdebugreaction(rxn::ElementaryReaction;tol=0.0,kf=0.0,krev=0.0,rt=0.0,ratio=0.0,T=0.0,P=0.0,index=0)
    return DebugReaction([getdebugspecies(spc) for spc in rxn.reactants],
            [getdebugspecies(spc) for spc in rxn.products],
            getrxnstr(rxn),tol,kf,krev,rt,ratio,T,P,index)
end

export getdebugreaction

function getdebugspecies(spc::Species;dy=0.0,ratio=0.0,tol=0.0,index=0)
    G298 = getGibbs(spc.thermo,298.0)
    return DebugSpecies(spc.name,G298,ratio,dy,tol,index)
end

export getdebugspecies
