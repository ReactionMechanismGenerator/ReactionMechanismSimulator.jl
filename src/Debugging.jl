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

"""
This function prints a report analyzing rates and dydt to determine reactions
and thermochemistry worth looking at
sim should be a simulation object created from a solution object corresponding to a crash
tol is the ratio relative to the median absolute rate or dn_i/dt value above which reactions
and species will be reported
"""
function analyzecrash(sim::Simulation;tol=1e6)
    rxns = Array{DebugReaction,1}()
    spcs = Array{DebugSpecies,1}()
    t = sim.sol.t[end]
    rts = rates(sim,t)
    rmedian = median(abs.([rt for rt in rts if !isnan(rt) && rt != 0.0]))
    dydt = zeros(length(sim.sol.u[end]))
    dydtreactor!(dydt,sim.sol.u[end],0.0,sim.domain,sim.interfaces;p=sim.p)
    dymedian = median(abs.([dy for dy in dydt if !isnan(dy) && dy != 0.0]))
    y = sim.sol(t)
    p = sim.domain.p
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(sim.domain,y,t,p)
    for (i,rt) in enumerate(rts)
        if isnan(rt)
            push!(rxns,getdebugreaction(sim.reactions[i];tol=tol,kf=kfs[i],krev=krevs[i],rt=rt,ratio=NaN,index=i))
        elseif abs(rt/rmedian) > tol
            ratio = abs(rt/rmedian)
            push!(rxns,getdebugreaction(sim.reactions[i];tol=tol,kf=kfs[i],krev=krev[i],rt=rt,ratio=ratio,index=i))
        end
    end
    for (i,dy) in enumerate(dydt)
        if i > length(sim.species)
            continue
        end
        if isnan(dy)
            push!(spcs,getdebugspecies(sim.species[i];dy=dy,ratio=NaN,tol=tol,index=i))
        elseif abs(dy/dymedian) > tol
            ratio = abs(dy/dymedian)
            push!(spcs,getdebugspecies(sim.species[i];dy=dy,ratio=ratio,tol=tol,index=i))
        end
    end
    return DebugMech(spcs,rxns)
end

"""
This function prints a report analyzing rates and dydt to determine reactions
and thermochemistry worth looking at
sim should be a simulation object created from a solution object corresponding to a crash
tol is the ratio relative to the median absolute rate or dn_i/dt value above which reactions
and species will be reported
"""
function analyzecrash(sim::SystemSimulation;tol=1e6)
    rxns = Array{DebugReaction,1}()
    spcs = Array{DebugSpecies,1}()
    t = sim.sol.t[end]
    rts = rates(sim,t)
    rmedian = median(abs.([rt for rt in rts if !isnan(rt) && rt != 0.0]))
    dydt = zeros(length(sim.sol.u[end]))
    dydtreactor!(dydt,sim.sol.u[end],0.0,tuple([s.domain for s in sim.sims]),sim.interfaces;p=sim.p)
    dymedian = median(abs.([dy for dy in dydt if !isnan(dy) && dy != 0.0]))
    y = sim.sol(t)
    p = sim.domain.p
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(sim.domain,y,t,p)
    for (i,rt) in enumerate(rts)
        if isnan(rt)
            push!(rxns,getdebugreaction(sim.reactions[i];tol=tol,kf=kfs[i],krev=krevs[i],rt=rt,ratio=NaN,index=i))
        elseif abs(rt/rmedian) > tol
            ratio = abs(rt/rmedian)
            push!(rxns,getdebugreaction(sim.reactions[i];tol=tol,kf=kfs[i],krev=krev[i],rt=rt,ratio=ratio,index=i))
        end
    end
    for (i,dy) in enumerate(dydt)
        if i > length(sim.species)
            continue
        end
        if isnan(dy)
            push!(spcs,getdebugspecies(sim.species[i];dy=dy,ratio=NaN,tol=tol,index=i))
        elseif abs(dy/dymedian) > tol
            ratio = abs(dy/dymedian)
            push!(spcs,getdebugspecies(sim.species[i];dy=dy,ratio=ratio,tol=tol,index=i))
        end
    end
    return DebugMech(spcs,rxns)
end
export analyzecrash
