using Parameters
using DifferentialEquations
include("Phase.jl")
include("State.jl")
include("PhaseState.jl")
include("Domain.jl")


abstract type AbstractReactor end
export AbstractReactor

struct BatchSingleDomainReactor{D<:AbstractDomain}
    domain::D
    ode::ODEProblem
end

function BatchSingleDomainReactor(domain::T,tspan::Tuple) where {T<:AbstractConstantKDomain}
    kfs = map(x->getkf(x,domain.phase,domain.state),domain.phase.reactions)
    println(kfs)
    recalcgibbs!(domain.phase,domain.state)
    keqs = map(x->getKc(x,domain.phase,domain.state),domain.phase.reactions)
    println(keqs)
    krevs = kfs./keqs
    println(krevs)
    N = length(domain.phase.species)
    dydt(y::Array{T,1},p::Nothing,t::T) where {T<:AbstractFloat} = dydtBatchReactor!(y,t,domain,kfs,krevs,N)
    println(dydt(domain.state.ns,nothing,0.0))
    ode = ODEProblem(dydt,domain.state.ns,tspan)
    return BatchSingleDomainReactor(domain,ode)
end

function BatchSingleDomainReactor(domain::T,tspan::Tuple) where {T<:AbstractVariableKDomain}
    N = length(domain.phase.species)
    dydt(y::Array{T,1},p::Nothing,t::T) where {T<:AbstractFloat} = dydtBatchReactor!(y,t,domain,N)
    ode = ODEProblem(dydt,domain.state.ns,tspan)
    return BatchSingleDomainReactor(domain,ode)
end
export BatchSingleDomainReactor

function getrate(rxn::T,state::MolarState,kfs::Array{Q,1},krevs::Array{Q,1}) where {T<:AbstractReaction,Q<:AbstractFloat}
    return kfs[rxn.index]*prod(state.cs[rxn.reactantinds])-krevs[rxn.index]*prod(state.cs[rxn.productinds])
end
export getrate

function addreactionratecontribution!(dydt::Array{Q,1},rxn::ElementaryReaction,st::MolarState,kfs::Array{Q,1},krevs::Array{Q,1}) where {Q<:Number,T<:Integer}
    R = getrate(rxn,st,kfs,krevs)
    for ind in rxn.reactantinds
        dydt[ind] -= R*st.V
    end
    for ind in rxn.productinds
        dydt[ind] += R*st.V
    end
end
export addreactionratecontribution!

function dydtBatchReactor!(y::Array{T,1},t::T,domain::Q,kfs::Array{T,1},krevs::Array{T,1},N::J) where {T<:AbstractFloat,J<:Integer,Q<:AbstractConstantKDomain}
    dydt = zeros(N)
    calcthermo!(domain,y,t)
    for rxn in domain.phase.reactions
        addreactionratecontribution!(dydt,rxn,domain.state,kfs,krevs)
    end
    calcdomainderivatives!(domain,dydt)
    return dydt
end

function dydtBatchReactor!(y::Array{T,1},t::T,domain::Q,N::J) where {J<:Integer,T<:Any,Q<:AbstractVariableKDomain}
    dydt = zeros(N)

    calcthermo!(domain)
    for rxn in domain.phase.reactions
        addreactionratecontribution!(dydt,rxn,domain.phase,domain.state)
    end
    calcdomainderivatives!(domain,dydt)
    return dydt
end
export dydtBatchReactor!
