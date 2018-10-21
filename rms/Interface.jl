include("Phase.jl")
include("Reactions.jl")

abstract type AbstractInterface end
export AbstractInterface

abstract type AbstractBoundaryInterface <: AbstractInterface end
export AbstractBoundaryInterface

abstract type AbstractInternalInterface <: AbstractInterface end
export AbstractInternalInterface

struct EmptyInterface <: AbstractInterface end
export EmptyInterface

@with_kw struct IdealGasCatalystInterface{T<:AbstractPhase,N<:AbstractPhase} <: AbstractInternalInterface
    gas::T
    catalyst::N
    reactions::Array{Reactions,1}
end
export IdealGasCatalystInterface
