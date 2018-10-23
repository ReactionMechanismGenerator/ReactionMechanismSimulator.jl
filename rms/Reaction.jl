using Parameters

include("Calculators.jl")
include("Species.jl")

abstract type AbstractReaction end
export AbstractReaction

@with_kw struct ElementaryReaction{T<:AbstractRate,N<:AbstractSpecies,Q<:Integer} <: AbstractReaction
    index::Q
    reactants::Array{N,1}
    reactantInds::Array{Q,1}
    products::Array{N,1}
    productInds::Array{Q,1}
    kinetics::T
end
export ElementaryReaction
