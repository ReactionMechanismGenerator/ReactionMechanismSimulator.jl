struct ReactionPath
    forward::Bool
    spcsinds::Array{Int64,1}
    rxninds::Array{Int64,1}
    spcind::Int64
    branchfracts::Array{Float64,1}
    branchfract::Array{Float64,1}
    branchind::Array{Int64,1}
end
export ReactionPath

"""
Take a single step following the fluxes in the network from a single species
checks whether the flux is moving the forward/reverse direction
based on production/loss < 1/steptol for forward and
production/loss > steptol for reverse direction and stops accordingly
"""
function stepnetwork(sim,statespcind,ropp,ropl,rts;forward=true,i=0,steptol=1e-2)
    ratio = sum(ropp[:,statespcind])/sum(ropl[:,statespcind])
    if forward
        if ratio > 1.0/steptol
            return (0,[0,])
        elseif i == 0
            rxnind = argmax(ropl[:,statespcind])
        else
            rxnind = reverse(sortperm(ropl[:,statespcind]))[i]
        end
        rt = rts[rxnind]
        if rt > 0.0
            spcinds = sim.reactions[rxnind].productinds
        else
            spcinds = sim.reactions[rxnind].reactantinds
        end
        return rxnind,spcinds
    else
        if ratio < steptol
            return (0,[0,])
        elseif i == 0
            rxnind = argmax(ropp[:,statespcind])
        else
            rxnind = reverse(sortperm(ropp[:,statespcind]))[i]
        end
        rt = rts[rxnind]
        if rt > 0.0
            spcinds = sim.reactions[rxnind].reactantinds
        else
            spcinds = sim.reactions[rxnind].productinds
        end
        return rxnind,spcinds
    end
end

"""
Starting from a given species indicated by spcind attempt to follow the fluxes
in the forward or reverse direction to find a reaction indicated by rxnind
it will follow in all directions until the flux tracked reaches branchtol fraction
of the original flux, the flux stops (based on steptol) or it finds the reaction
if it finds the reaction it will then follow the flux past that reaction using
a flux pairs algorithm until the flux stops
"""
function follow(sim,rxnind,spcind,ropp,ropl,rts,forward;steptol=1e-2,branchtol=5e-2)

    rps = [ReactionPath(forward,[spcind],Array{Int64,1}(),spcind,Array{Float64,1}(),
                    Array{Float64,1}([1.0]),Array{Int64,1}())]
    boo = true
    boo2 = true
    i = 0
    rpouts = Array{ReactionPath,1}()
    newrps = Array{ReactionPath,1}()
    while boo #extend off the center species
        newrps = extendpath(rps[1],sim,ropp,ropl,rts,i=i,steptol=steptol)
        if length(newrps) == 0
            return ReactionPath(forward,Array{Int64,1}(),Array{Int64,1}(),spcind,Array{Float64,1}(),
                    Array{Float64,1}(),Array{Int64,1}())
        elseif newrps[1].rxninds[1] == rxnind
            rpouts = newrps
            boo2 = false
            break
        end
        if newrps[1].branchfract[1] < branchtol
            boo = false
        else
            append!(rps,newrps)
            i += 1
        end
    end
    deleteat!(rps,1)

    while boo2
        ind = argmax([rp.branchfract[1] for rp in rps])
        rp = rps[ind]
        if rp.branchfract[1] < branchtol
            return ReactionPath(forward,Array{Int64,1}(),Array{Int64,1}(),spcind,Array{Float64,1}(),
                    Array{Float64,1}(),Array{Int64,1})
        end
        boo3 = true
        i = 0
        while boo3 #extend off the center species
            newrps = extendpath(rp,sim,ropp,ropl,rts,i=i,steptol=steptol)
            if length(newrps) == 0
                deleteat!(rps,ind)
                break
            end
            if newrps[1].rxninds[end] == rxnind
                rpouts = newrps
                boo2 = false
                break
            end
            if newrps[1].branchfract[1] < branchtol
                deleteat!(rps,ind)
                boo3 = false
            else
                deleteat!(rps,ind)
                append!(rps,newrps)
                i += 1
            end
        end
        if length(rps) == 0
            return ReactionPath(forward,Array{Int64,1}(),Array{Int64,1}(),spcind,Array{Float64,1}(),
                    Array{Float64,1}(),Array{Int64,1}())
        end
    end

    rpout = newrps[argmax([getsimilarity(sim.species[nrp.spcsinds[end-1]],sim.species[nrp.spcsinds[end]]) for nrp in rpouts])]
    while true
        newrps = extendpath(rpout,sim,ropp,ropl,rts,i=0,steptol=steptol)
        if length(newrps) == 0 || newrps[1].branchfract[1] < branchtol
            return rpout
        else
            rpout = newrps[argmax([getsimilarity(sim.species[rpout.spcsinds[end]],
                            sim.domain.phase.species[nrp.spcsinds[end]]) for nrp in newrps])]
        end
    end
end
export follow

"""
expand a ReactionPath object out a single step
"""
function extendpath(rp,sim,ropp,ropl,rts;i=0,steptol=1e-2)
    rxnind,spcinds = stepnetwork(sim,rp.spcsinds[end],ropp,ropl,rts,forward=rp.forward,i=i,steptol=steptol)
    if rxnind == 0
        sind = rp.spcsinds[end]
        name = sim.names[sind]
        forward = rp.forward
        return Array{ReactionPath,1}()
    else
        rps = Array{ReactionPath,1}()
        for spcind in spcinds
            rpnew = deepcopy(rp)
            push!(rpnew.rxninds,rxnind)
            push!(rpnew.spcsinds,spcind)
            push!(rpnew.branchind,i)
            if rpnew.forward
                push!(rpnew.branchfracts,ropl[rxnind,rp.spcsinds[end]]/sum(ropl[:,rp.spcsinds[end]]))
            else
                push!(rpnew.branchfracts,ropp[rxnind,rp.spcsinds[end]]/sum(ropp[:,rp.spcsinds[end]]))
            end
            rpnew.branchfract[1] = rp.branchfract[1]*rpnew.branchfracts[end]
            push!(rps,rpnew)
        end
        return rps
    end
end

struct Branching
    spcind::Int64
    rxninds::Array{Int64,1}
    branchingratios::Array{Float64,1}
end
export Branching

"""
Identifies calculates and stores the reactions and
branching ratios associated with the "reactants" for the input reaction
defined by rxnind as defined by the direction the reaction is occuring in
"""
function getbranching(sim,rxnind,ropl,rts)
    rt = rts[rxnind]
    if rt > 0
        spcs = sim.reactions[rxnind].reactants
    else
        spcs = sim.reactions[rxnind].products
    end
    spcsinds = [findfirst(isequal(spc.name),sim.names) for spc in spcs]
    boospcs = falses(length(spcsinds))
    branches = Array{Branching,1}()
    for (i,spc) in enumerate(spcs)
        ind = findfirst(isequal(spc.name),sim.names)
        rinds = reverse(sortperm(ropl[:,ind]))
        rvals = (ropl[:,ind]/sum(ropl[:,ind]))[rinds]
        indend = findfirst(isequal(0.0),rvals)
        rinds = rinds[1:indend]
        rvals = rvals[1:indend]
        push!(branches,Branching(ind,rinds,rvals))
    end
    return branches
end
export getbranching
