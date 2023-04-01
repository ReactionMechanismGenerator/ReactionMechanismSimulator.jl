using Parameters
using Statistics

"""
Object for storing useful debugging information about potentially problematic species
tol is the tolerance analyzed at, G is the Gibbs free energy, dy is
the flux to that species and ratio is the ratio between dy and the flux scale
index is the index of the species
"""
@with_kw struct DebugSpecies
    name::String
    G::Float64 = 0.0
    ratio::Float64 = 0.0
    dy::Float64 = 0.0
    tol::Float64 = 0.0
    index::Int64 = 0
end

"""
Object for storing useful debugging information about potentially problematic reactions
tol is the tolerance analyzed at, kf and krev are the forward and reverse rate coefficients
rt is the rate of the reaction, ratio is the ratio between rt and the rate scale
T and P are the temperature and pressure and index is the index of the reaction
"""
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
    return DebugReaction([getdebugspecies(spc,T) for spc in rxn.reactants],
            [getdebugspecies(spc,T) for spc in rxn.products],
            getrxnstr(rxn),tol,kf,krev,rt,ratio,T,P,index)
end

export getdebugreaction

function getdebugspecies(spc::Species,T::Float64;dy=0.0,ratio=0.0,tol=0.0,index=0)
    G = getGibbs(spc.thermo,T)
    return DebugSpecies(spc.name,G,ratio,dy,tol,index)
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
    T = getT(sim,t)
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
            push!(rxns,getdebugreaction(sim.reactions[i];tol=tol,kf=kfs[i],krev=krevs[i],rt=rt,T=T,ratio=NaN,index=i))
        elseif abs(rt/rmedian) > tol
            ratio = abs(rt/rmedian)
            push!(rxns,getdebugreaction(sim.reactions[i];tol=tol,kf=kfs[i],krev=krevs[i],rt=rt,T=T,ratio=ratio,index=i))
        end
    end
    for (i,dy) in enumerate(dydt)
        if i > length(sim.species)
            continue
        end
        if isnan(dy)
            push!(spcs,getdebugspecies(sim.species[i],T;dy=dy,ratio=NaN,tol=tol,index=i))
        elseif abs(dy/dymedian) > tol
            ratio = abs(dy/dymedian)
            push!(spcs,getdebugspecies(sim.species[i],T;dy=dy,ratio=ratio,tol=tol,index=i))
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
function analyzecrash(ssys::SystemSimulation;tol=1e6)
    rxns = Array{DebugReaction,1}()
    spcs = Array{DebugSpecies,1}()
    t = ssys.sol.t[end]
    rts = rates(ssys,t)
    rmedian = median(abs.([rt for rt in rts if !isnan(rt) && rt != 0.0]))
    dydt = zeros(length(ssys.sol.u[end]))
    dydtreactor!(dydt,ssys.sol.u[end],0.0,tuple([s.domain for s in ssys.sims]...),ssys.interfaces;p=ssys.p)
    dymedian = median(abs.([dy for dy in dydt if !isnan(dy) && dy != 0.0]))
    y = ssys.sol(t)
    p = ssys.p
    rdindex = 1
    sdindex = 1
    for sim in ssys.sims
        ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(sim.domain,y,t,p)
        for (i,rt) in enumerate(rts[rdindex:rdindex-1+length(sim.domain.phase.reactions)])
            if isnan(rt)
                push!(rxns,getdebugreaction(sim.reactions[i];tol=tol,kf=kfs[i],krev=krevs[i],rt=rt,ratio=NaN,index=i))
            elseif abs(rt/rmedian) > tol
                ratio = abs(rt/rmedian)
                push!(rxns,getdebugreaction(sim.reactions[i];tol=tol,kf=kfs[i],krev=krevs[i],rt=rt,ratio=ratio,index=i))
            end
        end
        rdindex += length(sim.domain.phase.reactions)
        for (i,dy) in enumerate(dydt[sdindex:sdindex-1+length(getphasespecies(sim.domain.phase))])
            if i > length(sim.species)
                continue
            end
            if isnan(dy)
                push!(spcs,getdebugspecies(sim.species[i],T;dy=dy,ratio=NaN,tol=tol,index=i))
            elseif abs(dy/dymedian) > tol
                ratio = abs(dy/dymedian)
                push!(spcs,getdebugspecies(sim.species[i],T;dy=dy,ratio=ratio,tol=tol,index=i))
            end
        end
        sdindex += length(getphasespecies(sim.domain.phase))
    end
    return DebugMech(spcs,rxns)
end
export analyzecrash

function getdebugmechstring(mech::DebugMech)
    s = "Crash Analysis Report\n"
    for rxn in mech.rxns
        s *= getdebugrxnstring(rxn)
    end
    for spc in mech.spcs
        s *= getdebugspeciesstring(spc)
    end
    return s
end

export getdebugmechstring

function getdebugrxnstring(rxn)
    s = ""
    rt = rxn.rt
    tol = rxn.tol
    ratio = rxn.ratio
    kf = rxn.kf
    krev = rxn.krev
    index = rxn.index
    if isnan(rt)
        s *= "NaN for Reaction:\n"
    elseif ratio > tol
        s *= "rt/rmedian > $tol, rt=$rt, rt/rmedian=$ratio \nReaction:\n"
    end
    s *= "Index: $index\n"
    s *= rxn.rxnstring * "\n"
    if length(rxn.reactants) == 1
        kunitsf = "s^-1"
    elseif length(rxn.reactants) == 2
        kunitsf = "m^3/(mol*s)"
    elseif length(rxn.reactants) == 3
        kunitsf = "m^6/(mol^2*s)"
    else
        error("cannot accomodate reactants>3")
    end
    if length(rxn.products) == 1
        kunitsr = "s^-1"
    elseif length(rxn.products) == 2
        kunitsr = "m^3/(mol*s)"
    elseif length(rxn.products) == 3
        kunitsr = "m^6/(mol^2*s)"
    else
        error("cannot accomodate products>3")
    end
    s *= "kf: $kf $kunitsf\n"
    s *= "krev: $krev $kunitsr\n"
    s *= "Reactants: \n"
    for reactant in rxn.reactants
        s *= getdebugspeciesstring(reactant)
    end
    s *= "Products: \n"
    for product in rxn.products
        s *= getdebugspeciesstring(product)
    end
    return s
end

export getdebugreactionstring

function getdebugspeciesstring(spc)
    s = ""
    dy = spc.dy
    ratio = spc.ratio
    tol = spc.tol
    index = spc.index
    if isnan(dy)
        s *= "NaN for Species net flux:\n"
    elseif ratio> tol
        s *= "dydt/dydtmedian > $tol, dydt=$dy, dydt/dydtmedian=$ratio \nSpecies:\n"
    end
    if index != 0
        s *= "Index: $index\n"
    end
    G = spc.G/4184.0
    name = spc.name
    s *= "\t$name G(T): $G kcal/mol\n"
    return s
end

export getdebugspeciesstring

function printcrashanalysis(mech::DebugMech)
    println(getdebugmechstring(mech))
end

function printcrashanalysis(sim::Simulation;tol=1e6)
    printcrashanalysis(analyzecrash(sim;tol=tol))
end

export printcrashanalysis

"""
This calculates the collision limit of H + H -> H2 as an upper bound
"""
function calccollisionlimit(T)
    mu = 0.00050397
    sigma = 2.05e-10
    epsilon = 1205.6
    Tr = T*kB*Na/epsilon
    collintegral = 1.16145*Tr^(-0.14874)+0.52487*exp(-0.7732 * Tr)+2.16178*exp(-2.437887*Tr)
    kcoll = sqrt(8.0*pi*kB*T*Na/mu)*sigma^2*collintegral*Na
    return kcoll
end

export calccollisionlimit

"""
Compares bimolecular and trimolecular reactions with the collision limit for H+H->H2
to check whether they are physical, a report with only a title and kcoll means no violators
"""
function analyzecolllimit(phase,Tmin,Tmax,Pmin,Pmax)
    println("Collision Limit Report (Comparisons are done with collision limit for H+H->H2)")
    println("kcoll=$kcoll")
    for T in [Tmin,Tmax]
        for P in [Pmin,Pmax]
            kcoll = calccollisionlimit(T)
            cpdivR,hdivRT,sdivR = calcHSCpdless(phase.vecthermo,T)
            Gs = (hdivRT.-sdivR)*(R*T)
            if phase.diffusionlimited
                diffs = getfield.(phase.species,:diffusion)(T=T,mu=0.0,P=P)
            else
                diffs = Array{Float64,1}()
            end
            kfs,krevs = getkfkrevs(phase,T,P,P/(R*T),1.0,Gs,diffs,V=1.0/C)
            for (i,rxn) in enumerate(phase.reactions)
                boo = false
                if len(rxn.reactants) > 2
                    if kfs[i] > kcoll
                        if !boo
                            println(getrxnstr(rxn))
                            boo = true
                        end
                        kf = kfs[i]
                        ratio = kfs[i]/kcoll
                        println("Violated in forward direction kf=$kf ratio=$ratio T=$T P=$P")
                    end
                end
                if len(rxn.products) > 2
                    if krevs[i] > kcoll
                        if !boo
                            println(getrxnstr(rxn))
                            boo = true
                        end
                        krev = krevs[i]
                        ratio = krevs[i]/kcoll
                        println("Violated in forward direction krev=$krev ratio=$ratio T=$T P=$P")
                    end
                end
            end
        end
    end
end

export analyzecolllimit
