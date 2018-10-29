using Parameters
include("Constants.jl")
include("Phase.jl")
include("State.jl")
include("Interface.jl")

abstract type AbstractDomain end
export AbstractDomain

abstract type AbstractConstantKDomain <: AbstractDomain  end
export AbstractConstantKDomain

abstract type AbstractVariableKDomain <: AbstractDomain end
export AbstractVariableKDomain

@with_kw struct ConstantTPDomain{T<:AbstractState,N<:AbstractPhase,S<:Integer} <: AbstractConstantKDomain
    state::T
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Array{S,1}
end
ConstantTPDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{EmptyInterface,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} = ConstantTPDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index])
export ConstantTPDomain

@with_kw struct ConstantVDomain{T<:AbstractState,N<:AbstractPhase,S<:Integer} <: AbstractVariableKDomain
    state::T
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Array{S,1}
end
ConstantVDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{AbstractInterface,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} = ConstantVDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index,0])
export ConstantVDomain

@with_kw struct ConstantTVDomain{T<:AbstractState,N<:AbstractPhase,S<:Integer} <: AbstractConstantKDomain
    state::T
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Array{S,1}
end
ConstantTVDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{EmptyInterface,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} = ConstantTVDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index])
export ConstantTVDomain
end
export ConstantTVDomain 
