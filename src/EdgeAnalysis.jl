"""
Tools for model edge analysis for automatic mechanism generation
"""

using Logging
using Sundials
using SparseArrays
using SciMLBase: build_solution
using Base.Iterators: flatten

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
function getsim(inte,react,edgereact,coreedgedomain,inters,p,coretoedgespcmap)
    ycoreedge = getycoreedge(inte.u,coretoedgespcmap,coreedgedomain.indexes[end],react,coreedgedomain)
    sol = build_solution(edgereact.ode,inte.alg,[0.0,inte.t],[ycoreedge,ycoreedge])
    sim = Simulation(sol,coreedgedomain,inters,p)
    return sim
end

"""
Generate appropriate edge Simulation/SystemSimulation object
"""
function getsim(inte,react,edgereact,coreedgedomains::Tuple,inters,p,coretoedgespcmap)
    ycoreedge = getycoreedge(inte.u,coretoedgespcmap,coreedgedomains[end].indexes[end],react,coreedgedomains)
    sol = build_solution(edgereact.ode,inte.alg,[0.0,inte.t],[ycoreedge,ycoreedge])
    ssys = SystemSimulation(sol,coreedgedomains,inters,p)
    return ssys
end

"""
Generate a state vector appropriate for both the edge and core from the core state vector
"""
function getycoreedge(y,coretoedgespcmap,edgelen,react,coreedgedomains::Tuple)
    ycoreedge = zeros(edgelen)
    for (coreind,edgeind) in coretoedgespcmap
        ycoreedge[edgeind] = y[coreind]
    end
    for (i,domain) in enumerate(coreedgedomains)
        for j in 3:length(domain.indexes)
            @inbounds ycoreedge[domain.indexes[j]] = y[react.domain[i].indexes[j]]
        end
    end
    return ycoreedge
end

"""
Generate a state vector appropriate for both the edge and core from the core state vector
"""
function getycoreedge(y,coretoedgespcmap,edgelen,react,coreedgedomain)
    ycoreedge = zeros(edgelen)
    for (coreind,edgeind) in coretoedgespcmap
        ycoreedge[edgeind] = y[coreind]
    end
    for j in 3:length(coreedgedomain.indexes)
        @inbounds ycoreedge[coreedgedomain.indexes[j]] = y[react.domain.indexes[j]]
    end
    return ycoreedge
end

"""
Calculate key flux and concentration related quantities for edge analysis
"""
@inline function calcfluxes(sim::Simulation)
    t = sim.sol.t[end]
    dydt = zeros(sim.domain.indexes[end])
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR,phi = calcthermo(sim.domain,sim.sol.u[end],sim.sol.t[end],SciMLBase.NullParameters())
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
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR,phi = calcthermo(domain,y,t,SciMLBase.NullParameters())
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
    rts,frts,rrts = addreactionratecontributionsforwardreverse!(dydt,domain.rxnarray,cstot,kfs,krevs,V)
    rtsall = [rts]
    frtsall = [frts]
    rrtsall = [rrts]
    for (i,domain) in enumerate(@views domains[2:end])
        k = i + 1
        @inbounds vns[k],vcs[k],vT[k],vP[k],vV[k],vC[k],vN[k],vmu[k],vkfs[k],vkrevs[k],vHs[k],vUs[k],vGs[k],vdiffs[k],vCvave[k],vcpdivR[k],vphi[k] = calcthermo(domain,y,t,SciMLBase.NullParameters())
        @inbounds cstot[domain.indexes[1]:domain.indexes[2]] .= vcs[k]
        rts,frts,rrts = addreactionratecontributionsforwardreverse!(dydt,domain.rxnarray,cstot,vkfs[k],vkrevs[k],vV[k])
        push!(rtsall,rts)
        push!(frtsall,frts)
        push!(rrtsall,rrts)
    end
    for (i,inter) in enumerate(interfaces)
        if isa(inter,ReactiveInternalInterface)
            @inbounds kfs,krevs = getkfskrevs(inter,vT[inter.domaininds[1]],vT[inter.domaininds[2]],vphi[inter.domaininds[1]],vphi[inter.domaininds[2]],vGs[inter.domaininds[1]],vGs[inter.domaininds[2]],cstot)
            @inbounds rts,frts,rrts = addreactionratecontributionsforwardreverse!(dydt,inter.rxnarray,cstot,kfs.*p[inter.parameterindexes[1]:inter.parameterindexes[2]],krevs.*p[inter.parameterindexes[1]:inter.parameterindexes[2]],inter.A)
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
        @inbounds calcdomainderivatives!(domain,dydt,interfaces;t=t,T=vT[i],P=vP[i],Us=vUs[i],Hs=vHs[i],V=vV[i],C=vC[i],ns=vns[i],N=vN[i],Cvave=vCvave[i])
    end
    return dydt,collect(flatten(rtsall)),collect(flatten(frtsall)),collect(flatten(rrtsall)),cstot
end
export calcfluxes

"""
Precalculate important indices and maps for use in edge analysis
"""
function getkeyselectioninds(coreedgedomains,coreedgeinters,domains,inters)
    @inbounds corespcsinds = collect(flatten([coreedgedomains[i].indexes[1]:coreedgedomains[i].indexes[1]+domains[i].indexes[2]-domains[i].indexes[1] for i = 1:length(domains)]))
    @inbounds edgespcsinds = collect(flatten([coreedgedomains[i].indexes[1]+domains[i].indexes[2]-domains[i].indexes[1]+1:coreedgedomains[i].indexes[2] for i = 1:length(domains)]))
    corerxninds = Array{Int64,1}()
    edgerxninds = Array{Int64,1}()
    Nrxns = sum(length(d.phase.reactions) for d in coreedgedomains)
    reactantindices = zeros(Int64,(4,Nrxns))
    productindices  = zeros(Int64,(4,Nrxns))
    coretoedgespcmap = Dict{Int64,Int64}()
    coretoedgerxnmap = Dict{Int64,Int64}()
    spcindexcore = 0
    spcindexedge = 0
    rxnindexcore = 0
    rxnindexedge = 0
    ind = 1
    for i = 1:length(domains)
        for (j,spc) in enumerate(domains[i].phase.species)
            edgeind = findfirst(isequal(spc),coreedgedomains[i].phase.species)
            @inbounds coretoedgespcmap[j+spcindexcore] = edgeind+spcindexedge
        end
        for j = 3:length(domains[i].indexes)
            @inbounds coretoedgespcmap[domains[i].indexes[j]] = coreedgedomains[i].indexes[j]
        end
        for (j,rxn) in enumerate(coreedgedomains[i].phase.reactions)
            @inbounds coreind = findfirst(x->rxn.reactants==x.reactants && rxn.products==x.products && rxn.kinetics==x.kinetics,domains[i].phase.reactions)
            if coreind === nothing
                push!(edgerxninds,j+rxnindexedge)
            else
                coretoedgerxnmap[coreind+rxnindexcore] = j+rxnindexedge
                push!(corerxninds,j+rxnindexedge)
            end
        end
        @inbounds spcindexcore += length(domains[i].phase.species)
        @inbounds spcindexedge += length(coreedgedomains[i].phase.species)
        @inbounds rxnindexcore += length(domains[i].phase.reactions)
        @inbounds rxnindexedge += length(coreedgedomains[i].phase.reactions)

        @inbounds indend = length(coreedgedomains[i].phase.reactions)
        @inbounds reactantindices[:,ind:ind+indend-1] = coreedgedomains[i].rxnarray[1:4,:]
        @inbounds productindices[:,ind:ind+indend-1] = coreedgedomains[i].rxnarray[5:8,:]
        ind += indend
    end

    for i = 1:length(inters)
        if @inbounds isa(inters[i],ReactiveInternalInterface)
            @inbounds push!(corerxnrangearray,index:index+length(inters[i].reactions))
            @inbounds push!(edgerxnrangearray,index+length(inters[i].reactions):index+length(coreedgeinters[i].reactions))
            @inbounds index += length(coreedgeinters[i].phase.reactions)
            for (j,rxn) in enumerate(coreedgeinters[i].reactions)
                @inbounds coreind = findfirst(x->rxn.reactants==x.reactants && rxn.products==x.products && rxn.kinetics==x.kinetics,inters[i].phase.reactions)
                if coreind === nothing
                    push!(edgerxninds,j+rxnindexedge)
                else
                    coretoedgerxnmap[coreind+rxnindexcore] = j+rxnindexedge
                    push!(corerxninds,j+rxnindexedge)
                end
            end
            @inbounds rxnindexcore += length(inters[i].reactions)
            @inbounds rxnindexedge += length(coreedgeinters[i].reactions)
            @inbounds indend = length(inters[i].reactions)
            @inbounds reactantindices[:,ind:ind+indend] = coreedgeinters[i].rxnarray[1:4,:]
            @inbounds productindices[:,ind:ind+indend] = coreedgeinters[i].rxnarray[5:8,:]
            ind += indend
        end
    end


    return corespcsinds,collect(flatten(corerxninds)),edgespcsinds,collect(flatten(edgerxninds)),reactantindices,productindices,coretoedgespcmap,coretoedgerxnmap
end

"""
Precalculate important indices and maps for use in edge analysis
"""
function getkeyselectioninds(coreedgedomain::AbstractDomain,coreedgeinters,domain,inters)
    corespcsinds = 1:length(domain.phase.species)
    edgespcsinds = length(domain.phase.species)+1:length(coreedgedomain.phase.species)
    corerxninds = zeros(Int64,length(domain.phase.reactions))
    corerxncount = 1
    edgerxninds = zeros(Int64,length(coreedgedomain.phase.reactions)-length(domain.phase.reactions))
    edgerxncount = 1
    coretoedgerxnmap = Dict{Int64,Int64}()
    for (j,rxn) in enumerate(coreedgedomain.phase.reactions)
        coreind = findfirst(x->rxn.reactants==x.reactants && rxn.products==x.products && rxn.kinetics==x.kinetics,domain.phase.reactions)
        if coreind === nothing
            @inbounds edgerxninds[edgerxncount] = j
            edgerxncount += 1
        else
            @inbounds coretoedgerxnmap[coreind] = j
            @inbounds corerxninds[corerxncount] = j
            corerxncount += 1
        end
    end

    @views @inbounds reactantindices = coreedgedomain.rxnarray[1:4,:]
    @views @inbounds productindices = coreedgedomain.rxnarray[5:8,:]
    coretoedgespcmap = Dict([i=>findfirst(isequal(spc),coreedgedomain.phase.species) for (i,spc) in enumerate(domain.phase.species)])
    @simd for j = 3:length(domain.indexes)
        @inbounds coretoedgespcmap[domain.indexes[j]] = coreedgedomain.indexes[j]
    end
    return corespcsinds,corerxninds,edgespcsinds,edgerxninds,reactantindices,productindices,coretoedgespcmap,coretoedgerxnmap
end

"""
Process flux information into useful quantities for edge analysis
"""
function processfluxes(sim::SystemSimulation,
        corespcsinds,corerxninds,edgespcsinds,edgerxninds)

    dydt,rts,frts,rrts,cs = calcfluxes(sim)

    corespeciesrates = abs.(dydt[corespcsinds])
    charrate = sqrt(dot(corespeciesrates,corespeciesrates))
    edgespeciesrates = abs.(dydt[edgespcsinds])
    edgereactionrates = rts[edgerxninds]
    corespeciesrateratios = corespeciesrates./charrate
    edgespeciesrateratios = edgespeciesrates./charrate
    corereactionrates = rts[corerxninds]
    corespeciesconcentrations = cs[corespcsinds]
    corespeciesconsumptionrates = zeros(length(corespeciesconcentrations))
    corespeciesproductionrates = zeros(length(corespeciesconcentrations))

    #process core species consumption and production rates
    index = 1
    for d in getfield.(sim.sims,:domain)
        for i = 1:size(d.rxnarray)[2]
            if @inbounds any(d.rxnarray[:,i].>length(corespeciesconcentrations))
                continue
            end
            for j = 1:4
                if @inbounds d.rxnarray[j,i] != 0
                    @inbounds corespeciesconsumptionrates[d.rxnarray[j,i]] += frts[i+index]
                    @inbounds corespeciesproductionrates[d.rxnarray[j,i]] += rrts[i+index]
                else
                    break
                end
            end
            for j = 5:8
                if @inbounds d.rxnarray[j,i] != 0
                    @inbounds corespeciesproductionrates[d.rxnarray[j,i]] += frts[i+index]
                    @inbounds corespeciesconsumptionrates[d.rxnarray[j,i]] += rrts[i+index]
                else
                    break
                end
            end
        end
        @inbounds index += size(d.rxnarray)[2]
    end
    for d in sim.interfaces
        if hasproperty(d,:rxnarray)
            @inbounds for i = 1:size(d.rxnarray)[2]
                if @inbounds any(d.rxnarray[:,i].>length(corespeciesconcentrations))
                    continue
                end
                for j = 1:4
                    if @inbounds d.rxnarray[j,i] != 0
                        @inbounds corespeciesconsumptionrates[d.rxnarray[j,i]] += frts[i+index]
                        @inbounds corespeciesproductionrates[d.rxnarray[j,i]] += rrts[i+index]
                    else
                        break
                    end
                end
                for j = 5:8
                    if @inbounds d.rxnarray[j,i] != 0
                        @inbounds corespeciesproductionrates[d.rxnarray[j,i]] += frts[i+index]
                        @inbounds corespeciesconsumptionrates[d.rxnarray[j,i]] += rrts[i+index]
                    else
                        break
                    end
                end
            end
            index += size(d.rxnarray)[2]
        end
    end

    return dydt,rts,frts,rrts,cs,corespeciesrates,charrate,edgespeciesrates,edgereactionrates,corespeciesrateratios,edgespeciesrateratios,corereactionrates,corespeciesconcentrations,corespeciesproductionrates,corespeciesconsumptionrates
end

"""
Process flux information into useful quantities for edge analysis
"""
function processfluxes(sim::Simulation,corespcsinds,corerxninds,edgespcsinds,edgerxninds)

    dydt,rts,frts,rrts,cs = calcfluxes(sim)

    @inbounds corespeciesrates = abs.(dydt[corespcsinds])
    charrate = sqrt(dot(corespeciesrates,corespeciesrates))
    @inbounds edgespeciesrates = abs.(dydt[edgespcsinds])
    @inbounds edgereactionrates = rts[edgerxninds]
    corespeciesrateratios = corespeciesrates./charrate
    edgespeciesrateratios = edgespeciesrates./charrate
    @inbounds corereactionrates = rts[corerxninds]
    @inbounds corespeciesconcentrations = cs[corespcsinds]
    corespeciesconsumptionrates = zeros(length(corespeciesconcentrations))
    corespeciesproductionrates = zeros(length(corespeciesconcentrations))

    #process core species consumption and production rates
    d = sim.domain
    @inbounds for i = 1:size(d.rxnarray)[2]
        if @inbounds  any(d.rxnarray[:,i].>length(corespeciesconcentrations))
            continue
        end
        for j = 1:4
            if @inbounds  d.rxnarray[j,i] != 0
                @inbounds corespeciesconsumptionrates[d.rxnarray[j,i]] += frts[i]
                @inbounds corespeciesproductionrates[d.rxnarray[j,i]] += rrts[i]
            else
                break
            end
        end
        for j = 5:8
            if @inbounds  d.rxnarray[j,i] != 0
                @inbounds corespeciesproductionrates[d.rxnarray[j,i]] += frts[i]
                @inbounds corespeciesconsumptionrates[d.rxnarray[j,i]] += rrts[i]
            else
                break
            end
        end
    end

    return dydt,rts,frts,rrts,cs,corespeciesrates,charrate,edgespeciesrates,edgereactionrates,corespeciesrateratios,edgespeciesrateratios,corereactionrates,corespeciesconcentrations,corespeciesproductionrates,corespeciesconsumptionrates
end

export processfluxes

"""
Calculate branching numbers for appropriate reactions for use in evaluating
the branching criterion: 1.0 < branchfactor * max(branchingratio,branchingratiomax) * rateratio^branchingindex
"""
function calcbranchingnumbers(sim,reactantinds,productinds,corespcsinds,corerxninds,edgerxninds,edgereactionrates,corespeciesrateratios,
        corespeciesconsumptionrates,branchfactor,branchingratiomax,branchingindex)
    branchingnums = zeros(length(edgereactionrates))
    for ind in 1:length(edgereactionrates)
        index = edgerxninds[ind]
        reactionrate = edgereactionrates[ind]
        if reactionrate > 0
            @views reactantside = reactantinds[:,index]
            @views productside = productinds[:,index]
        else
            @views reactantside = productinds[:,index]
            @views productside = reactantinds[:,index]
        end

        @inbounds rade = [sim.species[i].radicalelectrons for i in productside if i != 0]

        if maximum(rade) > 1
            continue
        end

        for spcindex in reactantside
            if spcindex == 0
                continue
            elseif spcindex in corespcsinds
                if @inbounds sim.species[spcindex].radicalelectrons != 1
                    continue
                end
                @inbounds consumption = corespeciesconsumptionrates[spcindex]
                if consumption != 0
                    br = reactionrate / consumption
                    @inbounds rr = corespeciesrateratios[spcindex]
                    if br > branchingratiomax
                        br = branchingratiomax
                    end

                    bnum = branchfactor * br * rr^branchingindex

                    if @inbounds bnum > branchingnums[ind]
                        @inbounds branchingnums[ind] = bnum
                    end
                end
            end
        end
    end
   return branchingnums
end

export calcbranchingnumbers

"""
determine species pairings that are concentrated enough that they should be reacted
"""
function updatefilterthresholds!(sim::Simulation,corespcsinds,corespeciesconcentrations,charrate,
        unimolecularthreshold,bimolecularthreshold,trimolecularthreshold,tolmovetocore,
        filterthreshold)
    unimolecularthresholdrateconstant,bimolecularthresholdrateconstant,trimolecularthresholdrateconstant = getthresholdrateconstants(sim,sim.domain.phase,filterthreshold)

    unimolecularthresholdval = tolmovetocore * charrate / unimolecularthresholdrateconstant
    bimolecularthresholdval = tolmovetocore * charrate / bimolecularthresholdrateconstant
    trimolecularthresholdval = tolmovetocore * charrate / trimolecularthresholdrateconstant

    for i in 1:length(corespcsinds)
        if @inbounds  !unimolecularthreshold[i]
            if @inbounds  corespeciesconcentrations[i] > unimolecularthresholdval
                @inbounds unimolecularthreshold[i] = true
            end
        end
    end
    for i in 1:length(corespcsinds)
        for j in i:length(corespcsinds)
            if @inbounds   !bimolecularthreshold[i,j]
                if @inbounds  corespeciesconcentrations[i] * corespeciesconcentrations[j] > bimolecularthresholdval
                    @inbounds bimolecularthreshold[i,j] = true
                end
            end
        end
    end
    for i in 1:length(corespcsinds)
        for j in i:length(corespcsinds)
            for k in j:length(corespcsinds)
                if @inbounds  !trimolecularthreshold[i,j,k]
                    if @inbounds  corespeciesconcentrations[i] * corespeciesconcentrations[j] * corespeciesconcentrations[k] > trimolecularthresholdval
                        @inbounds trimolecularthreshold[i,j,k] = true
                    end
                end
            end
        end
    end
end

export updatefilterthresholds!

"""
Determine species/reactions that should be added to the model core, react thresholding and
whether the simulation should be interrupted or terminated
"""
function identifyobjects!(sim,corespcsinds,corerxninds,edgespcsinds,
        edgerxninds,reactantinds,productinds,unimolecularthreshold,bimolecularthreshold,
        trimolecularthreshold,maxedgespeciesrateratios,tolmovetocore,tolinterruptsimulation,
        ignoreoverallfluxcriterion,filterreactions,maxnumobjsperiter,branchfactor,branchingratiomax,
        branchingindex,terminateatmaxobjects,termination,y0,invalidobjects,firsttime,
        filterthreshold,transitorydict,checktransitory)

    rxnarray = vcat()
    t = sim.sol.t[end]
    y = sim.sol.u[end]
    breakflag = false
    numcorespc = length(corespcsinds)
    numcorerxns = length(corerxninds)
    invalidobjectsprintboolean = true
    terminated = false
    conversion = 0.0

    (dydt,rts,frts,rrts,cs,corespeciesratse,charrate,edgespeciesrates,
    edgereactionrates,corespeciesrateratios,edgespeciesrateratios,
    corereactionrates,corespeciesconcentrations,corespeciesproductionrates,
    corespeciesconsumptionrates) = processfluxes(sim,corespcsinds,corerxninds,edgespcsinds,edgerxninds)

    for i = 1:length(edgespeciesrateratios)
        if @inbounds  edgespeciesrateratios[i] > maxedgespeciesrateratios[i]
            @inbounds maxedgespeciesrateratios[i] = edgespeciesrateratios[i]
        end
    end

    if charrate == 0 && length(edgereactionrates) > 0
        maxspeciesindex = argmax(edgespeciesrates)
        @inbounds maxspeciesrate = edgespeciesrates[maxspeciesindex]
        @inbounds index = edgespcsinds[maxspeciesindex]
        @inbounds name = sim.names[index]
        @info "at time $t s, species $name was added to model core to avoid singularity"
        @inbounds push!(invalidobjects,sim.species[index])
        return (false,true,0.0)
    end

    if branchfactor != 0.0 && !firsttime
        branchingnums = calcbranchingnumbers(sim,reactantinds,productinds,corespcsinds,corerxninds,edgerxninds,edgereactionrates,
            corespeciesrateratios,corespeciesconsumptionrates,branchfactor,branchingratiomax,branchingindex)
    end

    if filterreactions
        updatefilterthresholds!(sim,corespcsinds,corespeciesconcentrations,charrate,
            unimolecularthreshold,bimolecularthreshold,trimolecularthreshold,tolmovetocore,
            filterthreshold)
    end

    for term in termination
        if isa(term, TerminationTime)
            if t > term.time
                terminated = true
                @info "at time $t sec, reached target termination time"
            end
        elseif isa(term, TerminationRateRatio)
            if maxcharrate != 0.0 && charrate / maxcharrate < term.ratio
                terminated = true
                ratio = term.ratio
                @info "At time $t sec, reached target termination rate ratio $ratio"
            end
        else isa(term, TerminationConversion)
            index = findfirst(isequal(term.species.name),sim.names)
            @inbounds conversion = 1 - (y[index] / y0[index])
            @inbounds name = sim.species[index].name
            if conversion >= term.conversion
                terminated = true
                @info "At time $t sec, reeached target termination conversion $conversion of $name"
            end
        end
    end

    newobjectinds = Array{Int64,1}()
    newobjects = []
    newobjectvals = Array{Float64,1}()
    newobjecttype = []

    tempnewobjects = []
    tempnewobjectinds = Array{Int64,1}()
    tempnewobjectvals = Array{Float64,1}()
    tempnewobjecttype = []

    interrupt = false

    #movement of species to core based on rate ratios

    if !ignoreoverallfluxcriterion
        for (i,ind) in enumerate(edgespcsinds)
            @inbounds rr = edgespeciesrateratios[i]
            @inbounds obj = sim.species[ind]
            name = obj.name
            if rr > tolmovetocore
                if !(obj in newobjects || obj in invalidobjects)
                    push!(tempnewobjects,obj)
                    push!(tempnewobjectinds,ind)
                    push!(tempnewobjectvals,rr)
                    push!(tempnewobjecttype,"rr")
                end
            end
            if rr > tolinterruptsimulation
                @info "at time $t sec, species $name at $rr exceeded the minimum rate for simulation interruption of $tolinterruptsimulation"
                interrupt = true
            end
        end

        sortedinds = reverse(sortperm(tempnewobjectvals))

        for q in sortedinds
            @inbounds push!(newobjects,tempnewobjects[q])
            @inbounds push!(newobjectinds,tempnewobjectinds[q])
            @inbounds push!(newobjectvals,tempnewobjectvals[q])
            @inbounds push!(newobjecttype,tempnewobjecttype[q])
        end

        tempnewobjects = []
        tempnewobjectinds = Array{Int64,1}()
        tempnewobjectvals = Array{Float64,1}()
        tempnewobjecttype = []
    end

    if branchfactor != 0.0 && !firsttime
        for (i,ind) in enumerate(edgerxninds)
            @inbounds bnum = branchingnums[i]
            if bnum > 1
                @inbounds obj = sim.reactions[ind]
                if !(obj in newobjects || obj in invalidobjects)
                    push!(tempnewobjects,obj)
                    push!(tempnewobjectinds,ind)
                    push!(tempnewobjectvals,bnum)
                    push!(tempnewobjecttype,"branching")
                end
            end
        end
        sortedinds = reverse(sortperm(tempnewobjectvals))

        for q in sortedinds
            @inbounds push!(newobjects,tempnewobjects[q])
            @inbounds push!(newobjectinds,tempnewobjectinds[q])
            @inbounds push!(newobjectvals,tempnewobjectvals[q])
            @inbounds push!(newobjecttype,tempnewobjecttype[q])
        end

        tempnewobjects = []
        tempnewobjectinds = Array{Int64,1}()
        tempnewobjectvals = Array{Float64,1}()
        tempnewobjecttype = []
    end


    if (checktransitory || terminated) && length(transitorydict) > 0 && !firsttime
        transitoryoutdict = Dict()
        TS = transitorysensitivitiesfulltrapezoidal(sim,t;normalized=false)
        TS .*= sim.p'
        @views @inbounds TS = TS[:,length(sim.names)+1:end]
        for (spcname,tol) in transitorydict
            spcind = findfirst(isequal(spcname),sim.names)
            @views @inbounds  tsscale = maximum(abs.(TS[spcind,:]))
            for rind in edgerxninds
                @inbounds sens = TS[spcind,rind]/tsscale
                if abs(sens) > tol
                    @inbounds obj = sim.reactions[rind]
                    if !(obj in newobjects || obj in invalidobjects)
                        push!(tempnewobjects,obj)
                        push!(tempnewobjectinds,rind)
                        @inbounds transitoryoutdict[rind] = spcname
                        push!(tempnewobjectvals,sens)
                        push!(tempnewobjecttype,"transitorysensitivity")
                    end
                end
            end
        end

        sortedinds = reverse(sortperm(tempnewobjectvals))

        for q in sortedinds
            @inbounds push!(newobjects,tempnewobjects[q])
            @inbounds push!(newobjectinds,tempnewobjectinds[q])
            @inbounds push!(newobjectvals,tempnewobjectvals[q])
            @inbounds push!(newobjecttype,tempnewobjecttype[q])
        end

        tempnewobjects = []
        tempnewobjectinds = Array{Int64,1}()
        tempnewobjectvals = Array{Float64,1}()
        tempnewobjecttype = []
    end

    if length(invalidobjects) + length(newobjects) > maxnumobjsperiter
        if invalidobjectsprintboolean
            @info "Exceeded max number of objects...removing excess objects"
            invalidobjectsprintboolean = false
        end
        num = maxnumobjsperiter - length(invalidobjects)
        @inbounds newobjects = newobjects[1:num]
        @inbounds newobjectinds = newobjectinds[1:num]
        @inbounds newobjectvals = newobjectvals[1:num]
        @inbounds newobjecttype = newobjecttype[1:num]
    end

    if terminateatmaxobjects && length(invalidobjects) + length(newobjects) >= maxnumobjsperiter
        @info "Reached max number of objects...preparing to terminate"
        interrupt = true
    end

    if length(newobjects) > 0
        for (i,obj) in enumerate(newobjects)
            @inbounds val = newobjectvals[i]
            @inbounds ind = newobjectinds[i]
            if isa(obj, Species)
                name = obj.name
                @info "At time $t sec, species $name at rate ratio $val exceeded the minimum rate for moving to model core of $tolmovetocore"
            elseif isa(obj,ElementaryReaction)
                rstr = getrxnstr(obj)
                if @inbounds newobjecttype[i] == "branching"
                    @info "at time $t sec, reaction $rstr at a branching number of $val exceeded the threshold of 1 for moving to model core"
                elseif @inbounds newobjecttype[i] == "transitorysensitivity"
                    sens = val
                    @inbounds spcname = transitoryoutdict[ind]
                    @inbounds tol = transitorydict[spcname]
                    @info "at time $t sec, reaction $rstr at a normalized transitory sensitivity from $spcname of $sens exceeded the threshold of $tol for moving to model core"
                end
            end
        end

        append!(invalidobjects,newobjects)
    end

    if interrupt
        @info "Terminating simulation due to interrupt"
    end

    if terminated
        for term in termination
            if isa(term, TerminationConversion)
                index = findfirst(isequal(term.species.name),sim.names)
                @inbounds conversion = 1 - (y[index] / y0[index])
                @inbounds name = sim.species[index].name
                @info "$name conversion: $conversion"
            end
        end
    end

    return (terminated,interrupt,conversion)
end

export identifyobjects!

"""
run edge analysis to determine objects (species/reactions) that should be added to model core
"""
function selectobjects(react,edgereact,coreedgedomains,coreedgeinters,domains,inters,
                corep,coreedgep,tolmovetocore,tolinterruptsimulation,ignoreoverallfluxcriterion,filterreactions,
                maxnumobjsperiter,tolbranchrxntocore,branchingratiomax,
                branchingindex,terminateatmaxobjects,termination,
                filterthreshold,transitorydict,transitorystepperiod;
                atol=1e-20,rtol=1e-6,solver=CVODE_BDF())

    (corespcsinds,corerxninds,edgespcsinds,edgerxninds,reactantindices,
                productindices,coretoedgespcmap,coretoedgerxnmap) = getkeyselectioninds(coreedgedomains,coreedgeinters,domains,inters)

    unimolecularthreshold = falses(length(corespcsinds))
    bimolecularthreshold = falses((length(corespcsinds),length(corespcsinds)))
    trimolecularthreshold = falses((length(corespcsinds),length(corespcsinds),length(corespcsinds)))
    maxedgespeciesrateratios = zeros(length(edgespcsinds))
    invalidobjects = []
    terminated = false
    conversion = 0.0
    code = :Success

    if tolbranchrxntocore != 0.0
        branchfactor = 1.0/tolbranchrxntocore
    else
        branchfactor = 0.0
    end

    @inbounds tf = react.ode.tspan[2]
    inte = init(react.ode,solver,abstol=atol,reltol=rtol);

    t = inte.t
    sim = getsim(inte,react,edgereact,coreedgedomains,coreedgeinters,corep,coretoedgespcmap)

    y0 = sim.sol[end]
    spcsaddindices = Array{Int64,1}()
    firsttime = true

    nsteps = 0
    while t < tf && Symbol(code) == :Success
        step!(inte)
        nsteps += 1
        code = check_error(inte)
        t = inte.t
        checktransitory = nsteps % transitorystepperiod == 1
        sim = getsim(inte,react,edgereact,coreedgedomains,coreedgeinters,coreedgep,coretoedgespcmap)
        terminated,interrupt,conversion = identifyobjects!(sim,corespcsinds,corerxninds,edgespcsinds,
            edgerxninds,reactantindices,productindices,unimolecularthreshold,bimolecularthreshold,
                trimolecularthreshold,maxedgespeciesrateratios,tolmovetocore,tolinterruptsimulation,ignoreoverallfluxcriterion,filterreactions,
                maxnumobjsperiter,branchfactor,branchingratiomax,
                branchingindex,terminateatmaxobjects,termination,y0,invalidobjects,firsttime,
                filterthreshold,transitorydict,checktransitory)
        if firsttime
            firsttime = false
        end
        if terminated || interrupt
            break
        end
    end

    if Symbol(code) == :Success
        return (terminated,false,invalidobjects,unimolecularthreshold,
            bimolecularthreshold,trimolecularthreshold,maxedgespeciesrateratios,t,conversion)
    else
        @error "Solver failed with code $code resurrecting job"
        dydt,rts,frts,rrts,cs,corespeciesrates,charrate,edgespeciesrates,edgereactionrates,
        corespeciesrateratios,edgespeciesrateratios,corereactionrates,corespeciesconcentrations,
        corespeciesproductionrates,corespeciesconsumptionrates = processfluxes(sim,corespcsinds,corerxninds,edgespcsinds,edgerxninds)
        @inbounds ind = edgespcsinds[argmax(edgespeciesrates)]
        @inbounds invalidobjects = [sim.species[ind]]
        return (terminated,true,invalidobjects,unimolecularthreshold,
            bimolecularthreshold,trimolecularthreshold,maxedgespeciesrateratios,t,conversion)
    end

end

export selectobjects

"""
calculate threshold rate constants for different numbers of reactants
(to some extent the filter assumes rate constants can't be faster than these thresholds)
"""
function getthresholdrateconstants(sim::Simulation,phase,filterthreshold)
    return 2.08366122e10*getT(sim,sim.sol.t[end]),filterthreshold,filterthreshold/1.0e3
end

"""
calculate threshold rate constants for different numbers of reactants
(to some extent the filter assumes rate constants can't be faster than these thresholds)
"""
function getthresholdrateconstants(sim::Simulation,phase::IdealDiluteSolution,filterthreshold)
    T = getT(sim,sim.sol.t[end])
    mu = phase.solvent.mu(T)
    return 2.08366122e10*T,22.2*T/mu,0.11*T/mu
end

export getthresholdrateconstants
