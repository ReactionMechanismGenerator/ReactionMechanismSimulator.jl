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

function BatchSingleDomainReactor(domain::T,tspan::Tuple) where {T<:AbstractConstantKDomain}
    recalcgibbs!(domain.phase,domain.state)
    kfs,krevs = getkfkrevs(domain.phase,domain.state)
    N = length(domain.phase.species)
    dydt(y::Array{T,1},p::Nothing,t::T) where {T<:AbstractFloat} = dydtBatchReactor!(y,t,domain,kfs,krevs,N)
    ode = ODEProblem(dydt,domain.state.ns,tspan)
    return BatchSingleDomainReactor(domain,ode)
end

function BatchSingleDomainReactor(domain::T,tspan::Tuple) where {T<:AbstractVariableKDomain}
    N = length(domain.phase.species)+length(domain.indexes)-2 #Note that for multi domain reactors in the future indexing will have to be managed more carefully
    y0 = deepcopy(domain.state.ns)
    if length(y0) != N
        push!(y0,domain.state.T)
    end
    dydt(y::Array{T,1},p::Nothing,t::T) where {T<:AbstractFloat} = dydtBatchReactor!(y,t,domain,N)
    ode = ODEProblem(dydt,y0,tspan)
    return BatchSingleDomainReactor(domain,ode)
end
export BatchSingleDomainReactor

@inline function getrate(rxn::T,state::MolarState,kfs::Array{Q,1},krevs::Array{Q,1}) where {T<:AbstractReaction,Q<:AbstractFloat}
    Nreact = length(rxn.reactantinds)
    Nprod = length(rxn.productinds)
    R = 0.0
    if Nreact == 1
        @fastmath @inbounds R += kfs[rxn.index]*state.cs[rxn.reactantinds[1]]
    elseif Nreact == 2
        @fastmath @inbounds R += kfs[rxn.index]*state.cs[rxn.reactantinds[1]]*state.cs[rxn.reactantinds[2]]
    elseif Nreact == 3
        @fastmath @inbounds R += kfs[rxn.index]*state.cs[rxn.reactantinds[1]]*state.cs[rxn.reactantinds[2]]*state.cs[rxn.reactantinds[3]]
    end

    if Nprod == 1
        @fastmath @inbounds R -= krevs[rxn.index]*state.cs[rxn.productinds[1]]
    elseif Nprod == 2
        @fastmath @inbounds R -= krevs[rxn.index]*state.cs[rxn.productinds[1]]*state.cs[rxn.productinds[2]]
    elseif Nprod == 3
        @fastmath @inbounds R -= krevs[rxn.index]*state.cs[rxn.productinds[1]]*state.cs[rxn.productinds[2]]*state.cs[rxn.productinds[3]]
    end

    return R
end
export getrate

@inline function addreactionratecontribution!(dydt::Array{Q,1},rxn::ElementaryReaction,st::MolarState,kfs::Array{Q,1},krevs::Array{Q,1}) where {Q<:Number,T<:Integer}
    R = getrate(rxn,st,kfs,krevs)
    for ind in rxn.reactantinds
        @fastmath @inbounds dydt[ind] -= R
    end
    for ind in rxn.productinds
        @fastmath @inbounds dydt[ind] += R
    end
end
export addreactionratecontribution!

function dydtBatchReactor!(y::Array{T,1},t::T,domain::Q,kfs::Array{T,1},krevs::Array{T,1},N::J) where {T<:AbstractFloat,J<:Integer,Q<:AbstractConstantKDomain}
    dydt = zeros(N)
    calcthermo!(domain,y,t)
    @simd for rxn in domain.phase.reactions
        addreactionratecontribution!(dydt,rxn,domain.state,kfs,krevs)
    end
    dydt *= domain.state.V
    calcdomainderivatives!(domain,dydt)
    return dydt
end

function dydtBatchReactor!(y::Array{T,1},t::T,domain::Q,N::J) where {J<:Integer,T<:Any,Q<:AbstractVariableKDomain}
    dydt = zeros(N)
    calcthermo!(domain,y,t)
    kfs,krevs = getkfkrevs(domain.phase,domain.state)
    @simd for rxn in domain.phase.reactions
        addreactionratecontribution!(dydt,rxn,domain.state,kfs,krevs)
    end
    dydt *= domain.state.V
    calcdomainderivatives!(domain,dydt)
    return dydt
end
export dydtBatchReactor!
