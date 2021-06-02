"""
Tools for model edge analysis for automatic mechanism generation
"""

using Logging
using Sundials
using SparseArrays
using DiffEqBase: build_solution

abstract type AbstractTerminationCriterion end

struct TerminationTime <: AbstractTerminationCriterion
    time::Float64
end
struct TerminationConversion <: AbstractTerminationCriterion
    species::Species
    conversion::Float64
end
struct TerminationRateRatio <: AbstractTerminationCriterion
    ratio::Float64
end

export TerminationTime
export TerminationConversion
export TerminationRateRatio

"""
Generate appropriate core Simulation/SystemSimulation object
"""
function tosim(sol,domains,inters,p)
    return Simulation(sol,domains,inters,p)
end

"""
Generate appropriate core Simulation/SystemSimulation object
"""
function tosim(sol,domains::Tuple,inters,p)
   return SystemSimulation(sol,domains,inters,p) 
end

"""
Generate appropriate edge Simulation/SystemSimulation object
"""
function getsim(inte,react,coreedgedomain,inters,p,coretoedgespcmap)
    ycoreedge = getycoreedge(inte.u,coretoedgespcmap,coreedgedomain.indexes[end])
    sol = build_solution(react.ode,inte.alg,[0.0,inte.t],[ycoreedge,ycoreedge])
    sim = Simulation(sol,coreedgedomain,inters,p)
    return sim
end

"""
Generate appropriate edge Simulation/SystemSimulation object
"""
function getsim(inte,react,coreedgedomains::Tuple,inters,p,coretoedgespcmap)
    ycoreedge = getycoreedge(inte.u,coretoedgespcmap,coreedgedomains[end].indexes[end])
    sol = build_solution(react.ode,inte.alg,[0.0,inte.t],[ycoreedge,ycoreedge])
    ssys = SystemSimulation(sol,coreedgedomains,inters,p)
    return ssys
end

"""
Generate a state vector appropriate for both the edge and core from the core state vector
"""
function getycoreedge(y,coretoedgespcmap,numedgespc)
    ycoreedge = zeros(numedgespc)
    for (coreind,edgeind) in coretoedgespcmap
        ycoreedge[edgeind] = y[coreind]
    end
    return ycoreedge
end

"""
Calculate key flux and concentration related quantities for edge analysis
"""
@inline function calcfluxes(sim::Simulation)
    t = sim.sol.t[end]
    dydt = zeros(sim.domain.indexes[end])
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR,phi = calcthermo(sim.domain,sim.sol.u[end],sim.sol.t[end],DiffEqBase.NullParameters())
    rts,frts,rrts = addreactionratecontributionsforwardreverse!(dydt,sim.domain.rxnarray,cs,kfs,krevs,V)
    calcdomainderivatives!(sim.domain,dydt,[];t=t,T=T,P=P,Us=Us,Hs=Hs,V=V,C=C,ns=ns,N=N,Cvave=Cvave)
    return dydt,rts,frts,rrts,cs
end

"""
Calculate key flux and concentration related quantities for edge analysis
"""
@inline function calcfluxes(ssys::SystemSimulation)
    cstot = zeros(sum([length(sim.domain.phase.species) for sim in ssys.sims]))
    dydt = zeros(length(ssys.sol.u[end]))
    y = ssys.sol.u[end]
    t = ssys.sol.t[end]
    p = ssys.p
    domains = [sim.domain for sim in ssys.sims]
    interfaces = ssys.interfaces
    domain = domains[1]
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR,phi = calcthermo(domain,y,t,DiffEqBase.NullParameters())
    vns = Array{Any,1}(undef,length(domains))
    vns[1] = ns
    vcs = Array{Any,1}(undef,length(domains))
    vcs[1] = cs
    cstot[domain.indexes[1]:domain.indexes[2]] = cs
    vT = Array{Any,1}(undef,length(domains))
    vT[1] = T
    vP = Array{Any,1}(undef,length(domains))
    vP[1] = P
    vV = Array{Any,1}(undef,length(domains))
    vV[1] = V
    vC = Array{Any,1}(undef,length(domains))
    vC[1] = C
    vN = Array{Any,1}(undef,length(domains))
    vN[1] = N
    vmu = Array{Any,1}(undef,length(domains))
    vmu[1] = mu
    vkfs = Array{Any,1}(undef,length(domains))
    vkfs[1] = kfs
    vkrevs = Array{Any,1}(undef,length(domains))
    vkrevs[1] = krevs
    vHs = Array{Any,1}(undef,length(domains))
    vHs[1] = Hs
    vUs = Array{Any,1}(undef,length(domains))
    vUs[1] = Us
    vGs = Array{Any,1}(undef,length(domains))
    vGs[1] = Gs
    vdiffs = Array{Any,1}(undef,length(domains))
    vdiffs[1] = diffs
    vCvave = Array{Any,1}(undef,length(domains))
    vCvave[1] = Cvave
    vcpdivR = Array{Any,1}(undef,length(domains))
    vcpdivR[1] = cpdivR
    vphi = Array{Any,1}(undef,length(domains))
    vphi[1] = phi
    rtsall,frtsall,rrtsall = addreactionratecontributionsforwardreverse!(dydt,domain.rxnarray,cstot,kfs,krevs,V)
    rtsall = [rts]
    frtsall = [frts]
    rrtsall = [rrts]
    for (i,domain) in enumerate(@views domains[2:end])
        k = i + 1
        vns[k],vcs[k],vT[k],vP[k],vV[k],vC[k],vN[k],vmu[k],vkfs[k],vkrevs[k],vHs[k],vUs[k],vGs[k],vdiffs[k],vCvave[k],vcpdivR[k],vphi[k] = calcthermo(domain,y,t,DiffEqBase.NullParameters())
        cstot[domain.indexes[1]:domain.indexes[2]] .= vcs[k]
        rts,frts,rrts = addreactionratecontributionsforwardreverse!(dydt,domain.rxnarray,cstot,vkfs[k],vkrevs[k],vV[k])
        push!(rtsall,rts)
        push!(frtsall,frts)
        push!(rrtsall,rrts)
    end
    for (i,inter) in enumerate(interfaces)
        if isa(inter,ReactiveInternalInterface)
            kfs,krevs = getkfskrevs(inter,vT[inter.domaininds[1]],vT[inter.domaininds[2]],vphi[inter.domaininds[1]],vphi[inter.domaininds[2]],vGs[inter.domaininds[1]],vGs[inter.domaininds[2]],cstot)
            rts,frts,rrts = addreactionratecontributionsforwardreverse!(dydt,inter.rxnarray,cstot,kfs.*p[inter.parameterindexes[1]:inter.parameterindexes[2]],krevs.*p[inter.parameterindexes[1]:inter.parameterindexes[2]],inter.A)
            push!(rtsall,rts)
            push!(frtsall,frts)
            push!(rrtsall,rrts)
        elseif isa(inter,ReactiveInternalInterfaceConstantTPhi)
            rts,frts,rrts = addreactionratecontributionsforwardreverse!(dydt,inter.rxnarray,cstot,inter.kfs,inter.krevs,inter.A)
            push!(rtsall,rts)
            push!(frtsall,frts)
            push!(rrtsall,rrts)
        elseif isa(inter,AbstractReactiveInternalInterface)
            typ = typeof(inter)
            error("No handling available for AbstractReactiveInternalInterface: $typ")
        end
    end
    for (i,domain) in enumerate(domains)
        calcdomainderivatives!(domain,dydt,interfaces;t=t,T=vT[i],P=vP[i],Us=vUs[i],Hs=vHs[i],V=vV[i],C=vC[i],ns=vns[i],N=vN[i],Cvave=vCvave[i])
    end
    return dydt,rtsall,frtsall,rrtsall,cstot
end
export calcfluxes