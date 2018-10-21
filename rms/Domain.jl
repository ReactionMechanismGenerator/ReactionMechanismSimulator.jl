using Parameters
include("Phase.jl")
include("State.jl")
include("Interface.jl")

abstract type AbstractDomain end
export AbstractDomain

@with_kw mutable struct ConstantTPDomain{T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} <: AbstractDomain
    state::T
    phase::N
    iterfaces::Array{Q,1} = Array{EmptyInterface,1}()
end
export ConstantTPDomain

@with_kw mutable struct ConstantVDomain{T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} <: AbstractDomain
    state::T
    phase::N
    iterfaces::Array{Q,1} = Array{EmptyInterface,1}()
end
export ConstantVDomain

@with_kw mutable struct ConstantTVDomain{T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} <: AbstractDomain
    state::T
    phase::N
    iterfaces::Array{Q,1} = Array{EmptyInterface,1}()
end
export ConstantTVDomain 
