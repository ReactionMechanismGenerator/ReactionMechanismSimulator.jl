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
