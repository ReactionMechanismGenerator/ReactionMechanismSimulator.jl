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
