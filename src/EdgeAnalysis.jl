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

"""
Precalculate important indices and maps for use in edge analysis
"""
function getkeyselectioninds(coreeedgedomains,coreedgeinters,domains,inters)
    corespcsinds = flatten([coreedgedomains[i].indexes[1]:coreedgedomains[i].indexes[1]+domains[i].indexes[2]-domains[i].indexes[1] for i = 1:length(domains)])
    edgespcsinds = flatten([coreedgedomains[i].indexes[1]+domains[i].indexes[2]-domains[i].indexes[1]:coreedgedomains[i].indexes[2] for i = 1:length(domains)])
    corerxninds = Array{Int64,1}()
    edgerxninds = Array{Int64,1}()
    reactantindices = zeros(Int64,(3,length(corerxninds)))
    productindices  = zeros(Int64,(3,length(corerxninds)))
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
            coretoedgespcmap[j+spcindexcore] = edgeind+spcindexedge
        end
        for j = 3:length(domains[i].indexes)
            coretoedgespcmap[domains[i].indexes[j]] = coreedgedomains[i].indexes[j]
        end
        for (j,rxn) in enumerate(coreedgedomains[i].phase.reactions)
            coreind = findfirst(isequal(rxn),domains[i].phase.reactions)
            if coreind === nothing
                push!(edgerxninds,j+indexedge)
            else 
                coretoedgerxnmap[coreind+indexcore] = j+indexedge
                push!(corerxninds,j+indexedge)
            end
        end
        spcindexcore += length(domains[i].phase.species)
        spcindexedge += length(coreedgedomains[i].phase.species)
        rxnindexcore += length(domains[i].phase.reactions)
        rxnindexedge += length(coreedgedomains[i].phase.reactions)
        
        indend = length(domains[i].reactions)
        reactantindices[:,ind:ind+indend-1] = domains[i].rxnarray[1:3,:]
        productindices[:,ind:ind+indend-1] = domains[i].rxnarray[4:6,:]
        ind += indend
    end
        
    for i = 1:length(inters)
        if isa(inters[i],ReactiveInternalInterface)
            push!(corerxnrangearray,index:index+length(inters[i].reactions))
            push!(edgerxnrangearray,index+length(inters[i].reactions):index+length(coreedgeinters[i].reactions))
            index += length(coreedgeinters[i].phase.reactions)
            
            indend = length(inters[i].reactions)
            reactantindices[:,ind:ind+indend] = inters[i].rxnarray[1:3,:]
            productindices[:,ind:ind+indend] = inters[i].rxnarray[4:6,:]
            ind += indend
        end
    end
    corerxninds = flatten(corerxnrangearray)
    edgerxninds = flatten(edgerxnrangearray)
    
    return corespcsinds,corerxninds,edgespcsinds,edgerxninds,reactantindices,productindices,coretoedgespcmap,coretoedgerxnmap
end

"""
Precalculate important indices and maps for use in edge analysis
"""
function getkeyselectioninds(coreedgedomain::AbstractDomain,coreedgeinters,domain,inters)
    corespcsinds = 1:length(domain.phase.species)
    edgespcsinds = length(domain.phase.species)+1:length(coreedgedomain.phase.species)
    corerxninds = 1:length(domain.phase.reactions)
    edgerxninds = length(domain.phase.reactions)+1:length(coreedgedomain.phase.reactions)
    reactantindices = coreedgedomain.rxnarray[1:3,:]
    productindices = coreedgedomain.rxnarray[4:6,:]
    coretoedgespcmap = Dict([i=>findfirst(isequal(spc),coreedgedomain.phase.species) for (i,spc) in enumerate(domain.phase.species)])
    coretoedgerxnmap = Dict([i=>findfirst(isequal(rxn),coreedgedomain.phase.reactions) for (i,rxn) in enumerate(domain.phase.reactions)])
    for j = 3:length(domain.indexes)
        coretoedgespcmap[domain.indexes[j]] = coreedgedomain.indexes[j]
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
        for i = index:index+size(d.rxnarray)[2]
            if any(d.rxnarray[:,i].>length(corespeciesconcentrations))
                continue
            end
            for j = 1:3
                if d.rxnarray[j,i] != 0
                    corespeciesconsumptionrates[d.rxnarray[j,i]] += frts[i]
                    corespeciesproductionrates[d.rxnarray[j,i]] += rrts[i]
                else
                    break
                end
            end
            for j = 4:6
                if d.rxnarray[j,i] != 0
                    corespeciesproductionrates[d.rxnarray[j,i]] += frts[i]
                    corespeciesconsumptionrates[d.rxnarray[j,i]] += rrts[i]
                else
                    break
                end
            end
        end
        index += size(d.rxnarray)[2]
    end
    for d in inters
        for i = index:index+size(d.rxnarray)[2]
            if any(d.rxnarray[:,i].>length(corespeciesconcentrations))
                continue
            end
            for j = 1:3
                if d.rxnarray[j,i] != 0
                    corespeciesconsumptionrates[d.rxnarray[j,i]] += frts[i]
                    corespeciesproductionrates[d.rxnarray[j,i]] += rrts[i]
                else
                    break
                end
            end
            for j = 4:6
                if d.rxnarray[j,i] != 0
                    corespeciesproductionrates[d.rxnarray[j,i]] += frts[i]
                    corespeciesconsumptionrates[d.rxnarray[j,i]] += rrts[i]
                else
                    break
                end
            end
        end
        index += size(d.rxnarray)[2]
    end
    
    return dydt,rts,frts,rrts,cs,corespeciesratse,charrate,edgespeciesrates,edgereactionrates,corespeciesrateratios,edgespeciesrateratios,corereactionrates,corespeciesconcentrations,corespeciesproductionrates,corespeciesconsumptionrates
end

"""
Process flux information into useful quantities for edge analysis
"""
function processfluxes(sim::Simulation,corespcsinds,corerxninds,edgespcsinds,edgerxninds)
    
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
    d = sim.domain
    for i = 1:size(d.rxnarray)[2]
        if any(d.rxnarray[:,i].>length(corespeciesconcentrations))
            continue
        end
        for j = 1:3
            if d.rxnarray[j,i] != 0
                corespeciesconsumptionrates[d.rxnarray[j,i]] += frts[i]
                corespeciesproductionrates[d.rxnarray[j,i]] += rrts[i]
            else
                break
            end
        end
        for j = 4:6
            if d.rxnarray[j,i] != 0
                corespeciesproductionrates[d.rxnarray[j,i]] += frts[i]
                corespeciesconsumptionrates[d.rxnarray[j,i]] += rrts[i]
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
function calcbranchingnumbers(sim,reactantinds,productinds,corespcsinds,corerxninds,edgereactionrates,corespeciesrateratios,
        corespeciesconsumptionrates,branchfactor,branchingratiomax,branchingindex)
    branchingnums = zeros(length(edgereactionrates))
    for index in 1:length(edgereactionrates)
        reactionrate = edgereactionrates[index]
            
        if reactionrate > 0
            reactantside = reactantinds[:,index+length(corerxninds)]
            productside = productinds[:,index+length(corerxninds)]
        else
            reactantside = productinds[:,index+length(corerxninds)]
            productside = reactantinds[:,index+length(corerxninds)]
        end
            
        rade = [sim.species[i].radicalelectrons for i in productside if i != 0]
            
        if maximum(rade) > 1
            continue
        end
            
        for spcindex in reactantside
            if spcindex == 0 
                continue
            elseif spcindex < length(corespcsinds)
                if sim.species[spcindex].radicalelectrons != 1
                    continue
                end
                consumption = corespeciesconsumptionrates[spcindex]
                if consumption != 0
                    br = reactionrate / consumption
                    rr = corespeciesrateratios[spcindex]
                    if br > branchingratiomax
                        br = branchingratiomax
                    end
                        
                    bnum = branchfactor * br * rr^branchingindex
                        
                    if bnum > branchingnums[index]
                        branchingnums[index] = bnum
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
function updatefilterthresholds!(sim,corespcsinds,corespeciesconcentrations,charrate,
        unimolecularthreshold,bimolecularthreshold,trimolecularthreshold,tolmovetocore,
        filterthreshold)
    unimolecularthresholdrateconstant,bimolecularthresholdrateconstant,trimolecularthresholdrateconstant = getthresholdrateconstants(sim,sim.domain.phase,filterthreshold)
    
    unimolecularthresholdval = tolmovetocore * charrate / unimolecularthresholdrateconstant
    bimolecularthresholdval = tolmovetocore * charrate / bimolecularthresholdrateconstant
    trimolecularthresholdval = tolmovetocore * charrate / trimolecularthresholdrateconstant
        
    for i in 1:length(corespcsinds)
        if !unimolecularthreshold[i]
            if corespeciesconcentrations[i] > unimolecularthresholdval
                unimolecularthreshold[i] = true
            end
        end
    end
    for i in 1:length(corespcsinds)
        for j in i:length(corespcsinds)
            if  !bimolecularthreshold[i,j]
                if corespeciesconcentrations[i] * corespeciesconcentrations[j] > bimolecularthresholdval
                    bimolecularthreshold[i,j] = true
                end
            end
        end
    end
    for i in 1:length(corespcsinds)
        for j in i:length(corespcsinds)
            for k in j:length(corespcsinds)
                if !trimolecularthreshold[i,j,k]
                    if corespeciesconcentrations[i] * corespeciesconcentrations[j] * corespeciesconcentrations[k] > trimolecularthresholdval
                        trimolecularthreshold[i,j,k] = true
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
        filterthreshold)
    
    rxnarray = vcat()
    t = sim.sol.t[end]
    y = sim.sol.u[end]
    breakflag = false
    numcorespc = length(corespcsinds)
    numcorerxns = length(corerxninds)
    invalidobjectsprintboolean = true
    terminated = false
    
    (dydt,rts,frts,rrts,cs,corespeciesratse,charrate,edgespeciesrates,
    edgereactionrates,corespeciesrateratios,edgespeciesrateratios,
    corereactionrates,corespeciesconcentrations,corespeciesproductionrates,
    corespeciesconsumptionrates) = processfluxes(sim,corespcsinds,corerxninds,edgespcsinds,edgerxninds)
    
    for i = 1:length(edgespeciesrateratios)
        if edgespeciesrateratios[i] > maxedgespeciesrateratios[i]
            maxedgespeciesrateratios[i] = edgespeciesrateratios[i]
        end
    end
    
    if charrate == 0 && length(edgereactionrates) > 0
        maxspeciesindex = argmax(edgespeciesrates)
        maxspeciesrate = edgespeciesrates[maxspeciesindex]
        name = sim.names[maxspeciesindex]
        @info "at time $t s, species $name was added to model core to avoid singularity"
        push!(invalidobjects,sim.species[maxspeciesindex])
        return (false,true)
    end
    
    if branchfactor != 0.0 && !firsttime
        branchingnums = calcbranchingnumbers(sim,reactantinds,productinds,corespcsinds,corerxninds,edgereactionrates,
            corespeciesrateratios,corespeciesconsumptionrates,branchfactor,branchingratiomax,branchingindex)
    end
    
    if filterreactions
        updatefilterthresholds!(sim,corespcsinds,corespeciesconcentrations,charrate,
            unimolecularthreshold,bimolecularthreshold,trimolecularthreshold,tolmovetocore,
            filterthreshold)
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
            rr = edgespeciesrateratios[i]
            obj = sim.species[ind]
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
            push!(newobjects,tempnewobjects[q])
            push!(newobjectinds,tempnewobjectinds[q])
            push!(newobjectvals,tempnewobjectvals[q])
            push!(newobjecttype,tempnewobjecttype[q])
        end
        
        tempnewobjects = []
        tempnewobjectinds = Array{Int64,1}()
        tempnewobjectvals = Array{Float64,1}()
        tempnewobjecttype = []
    end
    
    if branchfactor != 0.0 && !firsttime
        for (i,ind) in enumerate(edgerxninds)
            bnum = branchingnums[i]
            if bnum > 1
                obj = sim.reactions[ind+numcorerxns]
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
            push!(newobjects,tempnewobjects[q])
            push!(newobjectinds,tempnewobjectinds[q])
            push!(newobjectvals,tempnewobjectvals[q])
            push!(newobjecttype,tempnewobjecttype[q])
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
        newobjects = newobjects[1:num]
        newobjectinds = newobjectinds[1:num]
        newobjectvals = newobjectvals[1:num]
        newobjecttype = newobjecttype[1:num]
    end
    
    if terminateatmaxobjects && length(invalidobjects) + length(newobjects) >= maxnumobjsperiter
        @info "Reached max number of objects...preparing to terminate"
        interrupt = true
    end
    
    if length(newobjects) > 0
        for (i,obj) in enumerate(newobjects)
            val = newobjectvals[i]
            ind = newobjectinds[i]
            if isa(obj, Species)
                name = obj.name
                @info "At time $t sec, species $name at rate ratio $val exceeded the minimum rate for moving to model core of $tolmovetocore"
            elseif isa(obj,ElementaryReaction)
                rstr = getrxnstr(obj)
                @info "at time $t sec, reaction $rstr at a branching number of $val exceeded the threshold of 1 for moving to model core"
            end
        end
        
        append!(invalidobjects,newobjects)
    end
    
    if interrupt
        @info "Terminating simulation due to interrupt"
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
            conversion = 1 - (y[index] / y0[index])
            name = sim.species[index].name
            if conversion >= term.conversion
                terminated = true
                @info "At time $t sec, reeached target termination conversion $conversion of $name"
            end
        end
    end
    
    if terminated
        for term in termination
            if isa(term, TerminationConversion)
                index = findfirst(isequal(term.species.name),sim.names)
                conversion = 1 - (y[index] / y0[index])
                name = sim.species[index].name
                @info "$name conversion: $conversion"
            end
        end
    end
    
    return (terminated,interrupt) 
end

export identifyobjects!

"""
run edge analysis to determine objects (species/reactions) that should be added to model core
"""
function selectobjects(react,coreedgedomains,coreedgeinters,domains,inters,
                p,tolmovetocore,tolinterruptsimulation,ignoreoverallfluxcriterion,filterreactions,
                maxnumobjsperiter,tolbranchrxntocore,branchingratiomax,
                branchingindex,terminateatmaxobjects,termination,
                filterthreshold;
                atol=1e-20,rtol=1e-6,solver=CVODE_BDF())
            
    (corespcsinds,corerxninds,edgespcsinds,edgerxninds,reactantindices,
                productindices,coretoedgespcmap,coretoedgerxnmap) = getkeyselectioninds(coreedgedomains,coreedgeinters,domains,inters)
    
    unimolecularthreshold = falses(length(corespcsinds))
    bimolecularthreshold = falses((length(corespcsinds),length(corespcsinds)))
    trimolecularthreshold = falses((length(corespcsinds),length(corespcsinds),length(corespcsinds)))
    maxedgespeciesrateratios = zeros(length(edgespcsinds))
    invalidobjects = []
    terminated = false
    
    if tolbranchrxntocore != 0.0
        branchfactor = 1.0/tolbranchrxntocore
    else
        branchfactor = 0.0
    end
    
    tf = react.ode.tspan[2]
    inte = init(react.ode,solver,abstol=atol,reltol=rtol);
    
    t = inte.t
    sim = getsim(inte,react,coreedgedomains,inters,p,coretoedgespcmap)
    
    y0 = sim.sol[end]
    spcsaddindices = Array{Int64,1}()
    firsttime = true
    
    n = 1
    while t < tf
        for i = 1:n
            step!(inte)
        end
        t = inte.t
        sim = getsim(inte,react,coreedgedomains,inters,p,coretoedgespcmap)
        terminated,interrupt = identifyobjects!(sim,corespcsinds,corerxninds,edgespcsinds,
            edgerxninds,reactantindices,productindices,unimolecularthreshold,bimolecularthreshold,
                trimolecularthreshold,maxedgespeciesrateratios,tolmovetocore,tolinterruptsimulation,ignoreoverallfluxcriterion,filterreactions,
                maxnumobjsperiter,branchfactor,branchingratiomax,
                branchingindex,terminateatmaxobjects,termination,y0,invalidobjects,firsttime,
                filterthreshold)
        if firsttime
            firsttime = false
        end
        if terminated || interrupt
            break
        end
    end
    return (terminated,invalidobjects,unimolecularthreshold,
        bimolecularthreshold,trimolecularthreshold,maxedgespeciesrateratios)
end

export selectspecies

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