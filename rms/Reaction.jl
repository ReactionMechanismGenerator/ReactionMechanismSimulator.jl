using Parameters

include("Calculators.jl")
include("Species.jl")

abstract type AbstractReaction end
export AbstractReaction

@with_kw struct Reaction{T<:AbstractRate,N<:AbstractSpecies} <: AbstractReaction
    reactants::Array{N,1}
    products::Array{N,1}
    kinetics::T
end
export Reaction
