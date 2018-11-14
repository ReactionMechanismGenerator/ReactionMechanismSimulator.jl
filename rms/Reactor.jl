using Parameters
using DifferentialEquations
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
