using Parameters
using DifferentialEquations
using ForwardDiff
include("Phase.jl")
include("PhaseState.jl")
include("Domain.jl")

abstract type AbstractReactor end
export AbstractReactor

struct BatchSingleDomainReactor{D<:AbstractDomain}
    domain::D
    ode::ODEProblem
end


function BatchReactor(domain::T,y0::Array{W,1},tspan::Tuple) where {T<:AbstractDomain,W<:Real}
    dydt(y::Array{T,1},p::Nothing,t::T) where {T<:Real} = dydtBatchReactor!(y,t,domain)
    ode = ODEProblem(dydt,y0,tspan)
    return BatchSingleDomainReactor(domain,ode)
end
export BatchReactor

@inline function getrate(rxn::T,cs::Array{W,1},kfs::Array{Q,1},krevs::Array{Q,1}) where {T<:AbstractReaction,Q,W<:Real}
    Nreact = length(rxn.reactantinds)
    Nprod = length(rxn.productinds)
    R = 0.0
    if Nreact == 1
        @fastmath @inbounds R += kfs[rxn.index]*cs[rxn.reactantinds[1]]
    elseif Nreact == 2
        @fastmath @inbounds R += kfs[rxn.index]*cs[rxn.reactantinds[1]]*cs[rxn.reactantinds[2]]
    elseif Nreact == 3
        @fastmath @inbounds R += kfs[rxn.index]*cs[rxn.reactantinds[1]]*cs[rxn.reactantinds[2]]*cs[rxn.reactantinds[3]]
    end

    if Nprod == 1
        @fastmath @inbounds R -= krevs[rxn.index]*cs[rxn.productinds[1]]
    elseif Nprod == 2
        @fastmath @inbounds R -= krevs[rxn.index]*cs[rxn.productinds[1]]*cs[rxn.productinds[2]]
    elseif Nprod == 3
        @fastmath @inbounds R -= krevs[rxn.index]*cs[rxn.productinds[1]]*cs[rxn.productinds[2]]*cs[rxn.productinds[3]]
    end

    return R
end
export getrate

@inline function addreactionratecontribution!(dydt::Array{Q,1},rxn::ElementaryReaction,cs::Array{W,1},kfs::Array{Z,1},krevs::Array{Z,1}) where {Q<:Real,Z<:Real,T<:Integer,W<:Real}
    R = getrate(rxn,cs,kfs,krevs)
    for ind in rxn.reactantinds
        @fastmath @inbounds dydt[ind] -= R
    end
    for ind in rxn.productinds
        @fastmath @inbounds dydt[ind] += R
    end
end
export addreactionratecontribution!

function dydtBatchReactor!(y::Array{U,1},t::Z,domain::Q) where {Z<:Real,U<:Real,J<:Integer,Q<:AbstractDomain}
    dydt = zeros(U,length(y))
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs = calcthermo(domain,y,t)
    @simd for rxn in domain.phase.reactions
        addreactionratecontribution!(dydt,rxn,cs,kfs,krevs)
    end
    dydt *= V
    calcdomainderivatives!(domain,dydt;T=T,Us=Us,V=V,C=C,ns=ns,N=N)
    return dydt
end
export dydtBatchReactor!

function jacobianbatch!(J::Q,y::U,p::W,t::Z,domain::V) where {Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
    if domain.t[1] == t && domain.jacuptodate == true
        return domain.jacobian
    else
        f(y::Array{T,1}) where {T<:Real} = dydtBatchReactor!(y,domain.t[1],domain;sensitivity=false)
        ForwardDiff.jacobian!(domain.jacobian,f,y[1:length(domain.phase.species)])
        domain.jacuptodate[1] = true
        return domain.jacobian
    end
end

function ratederivative(;d::W,cs::Q,V::Y,T::Y2,kfs::Z,krevs::X,sparse::Bool=false) where {W<:Union{ConstantTPDomain,ConstantTVDomain},Q<:AbstractArray,Y2<:Real,Y<:Real,Z<:AbstractArray,X<:AbstractArray}
    Nspcs = length(cs)
    rxns = d.phase.reactions
    Nrxns = length(rxns)

    if sparse
        ratederiv = spzeros(Nspcs,Nspcs+Nrxns)
    else
        ratederiv = zeros(Nspcs,Nspcs+Nrxns)
    end
    RTinv = 1.0/(R*T)

    for (j,rxn) in enumerate(rxns)
        Nreact = length(rxn.reactantinds)
        Nprod = length(rxn.productinds)

        if Nreact == 1
            rind1 = rxn.reactantinds[1]
            fderiv = cs[rind1]
        elseif Nreact == 2
            rind1,rind2 = rxn.reactantinds
            fderiv = cs[rind1]*cs[rind2]
        else
            rind1,rind2,rind3 = rxn.reactantinds
            fderiv = cs[rind1]*cs[rind2]*cs[rind3]
        end

        if Nprod == 1
            pind1 = rxn.productinds[1]
            rderiv = krevs[j]/kfs[j]*cs[pind1]
        elseif Nprod == 2
            pind1,pind2 = rxn.productinds
            rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]
        else
            pind1,pind2,pind3 = rxn.productinds
        end

        flux = fderiv-rderiv
        gderiv = rderiv*kfs[j]*RTinv

        deriv = zeros(Nspcs)

        deriv[rind1] += gderiv
        if Nreact > 1
            deriv[rind2] += gderiv
            if Nreact > 2
                deriv[rind3] == gderiv
            end
        end

        deriv[pind1] -= gderiv
        if Nprod > 1
            deriv[pind2] -= gderiv
            if Nprod > 2
                deriv[pind3] -= gderiv
            end
        end

        ratederiv[rind1,j] -= flux
        ratederiv[rind1,Nrxns+1:Nrxns+Nspcs] .-= deriv
        if Nreact > 1
            ratederiv[rind2,j] -= flux
            ratederiv[rind2,Nrxns+1:Nrxns+Nspcs] .-= deriv
            if Nreact > 2
                ratederiv[rind3,j] -= flux
                ratederiv[rind3,Nrxns+1:Nrxns+Nspcs] .-= deriv
            end
        end

        ratederiv[pind1,j] += flux
        ratederiv[pind1,Nrxns+1:Nrxns+Nspcs] .+= deriv
        if Nprod > 1
            ratederiv[pind2,j] += flux
            ratederiv[pind2,Nrxns+1:Nrxns+Nspcs] .+= deriv
            if Nprod > 2
                ratederiv[pind3,j] += flux
                ratederiv[pind3,Nrxns+1:Nrxns+Nspcs] .+= deriv
            end
        end
    end
    return V*ratederiv
end
export ratederivative
