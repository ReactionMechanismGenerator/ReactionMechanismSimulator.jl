include("Phase.jl")
include("Reaction.jl")

abstract type AbstractInterface end
export AbstractInterface

abstract type AbstractBoundaryInterface <: AbstractInterface end
export AbstractBoundaryInterface

abstract type AbstractInternalInterface <: AbstractInterface end
export AbstractInternalInterface

struct EmptyInterface <: AbstractInterface end
export EmptyInterface

@with_kw struct IdealGasCatalystInterface{T<:AbstractPhase,N<:AbstractPhase,Q<:AbstractReaction} <: AbstractInternalInterface
    gas::T
    catalyst::N
    reactions::Array{Q,1}
end
export IdealGasCatalystInterface
