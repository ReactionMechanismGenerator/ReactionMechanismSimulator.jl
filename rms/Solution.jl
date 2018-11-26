using DifferentialEquations
import DifferentialEquations.DiffEqBase: AbstractODESolution, HermiteInterpolation,AbstractDiffEqInterpolation

abstract type AbstractSolution end

struct BatchSolution{Q<:AbstractODESolution,W<:AbstractDomain,L<:AbstractArray,G<:Function} <: AbstractSolution
    sol::Q
    domain::W
    names::L
    N::G
end

function BatchSolution(sol::Q,domain::W) where {Q<:AbstractODESolution,W<:AbstractDomain}
    names = getfield.(domain.phase.species,:name)
    Ns = sum(hcat(sol.interp.u...)[domain.indexes[1]:domain.indexes[2],:],dims=1)
    Nderivs = sum(hcat(sol.interp.du...)[domain.indexes[1]:domain.indexes[2],:],dims=1)
    N = HermiteInterpolation(sol.interp.t,Ns,Nderivs)
    F(t::T) where {T<:Real} = N(t,nothing,Val{0},sol.prob.p,:left)
    return BatchSolution(sol,domain,names,F)
end

function molefractions(bsol::Q; name::W,t::E) where {Q<:AbstractSolution, W<:String, E<:Real}
    @assert name in bsol.names
    ind = findfirst(isequal(name),bsol.names)
    return bsol(t)[ind]/bsol.N(t)
end

function molefractions(bsol::Q; t::E) where {Q<:AbstractSolution,E<:Real}
    return bsol.sol(t)[bsol.domain.indexes[1]:bsol.domain.indexes[2]]./bsol.N(t)
end

function molefractions(bsol::Q) where {Q<:AbstractSolution}
    return bsol.sol.u./bsol.N.u
end

getT(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:Union{ConstantTPDomain,ConstantTVDomain},K<:Real,Q,G,L} = bsol.domain.T
getT(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:ConstantVDomain,K<:Real,Q,G,L} = bsol.sol(t)[bsol.domain.indexes[3]]

getV(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.N(t)*getT(bsol,t)*R /bsol.domain.P
getV(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain},K<:Real,Q,G,L} = bsol.domain.V

getP(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.domain.P
getP(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:ConstantTVDomain,K<:Real,Q,G,L} = 1.0e6
getP(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:ConstantVDomain,K<:Real,Q,G,L} = bsol.N(t)*R*getT(bsol,t)/getV(bsol,t)

getC(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.domain.P/(R*bsol.domain.T)
getC(bsol::BatchSolution{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain},K<:Real,Q,G,L} = bsol.N(t)/bsol.domain.V

"""
calculates the rates of production/loss at a given time point
this outputs a sparse matrix of  num reactions xnum species containing the production/loss
rate of that species associated with that reaction
"""
function rops(bsol,t)
    ropmat = spzeros(length(bsol.domain.phase.reactions),length(bsol.domain.phase.species))
    xs = molefractions(bsol,t=t)
    T = getT(bsol,t)
    V = getV(bsol,t)
    P = getP(bsol,t)
    if :Gs in fieldnames(typeof(bsol.domain))
        Gs = domain.Gs
    else
        Gs = calcgibbs(bsol.domain.phase,T)
    end
    if :solvent in fieldnames(typeof(bsol.domain.phase)) && typeof(bsol.domain.phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if bsol.domain.phase.diffusionlimited
        diffs = getfield.(phase.species,:diffusion)(T=T,mu=mu,P=P)
    else
        diffs = Array{typeof(T),1}()
    end
    kfs,krevs = getkfkrevs(phase=bsol.domain.phase,T=T,P=P,C=1.0/V,N=1.0,ns=xs,Gs=Gs,diffs=diffs)
    cs = xs./V
    for (i,rxn) in enumerate(bsol.domain.phase.reactions)
        R = getrate(rxn,cs,kfs,krevs)
        for ind in rxn.reactantinds
            ropmat[i,ind] = -R
        end
        for ind in rxn.productinds
            ropmat[i,ind] = R
        end
    end
    return ropmat
end
