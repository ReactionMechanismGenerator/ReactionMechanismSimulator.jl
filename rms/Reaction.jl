using Parameters
import Base: length
include("Calculators.jl")
include("Species.jl")

abstract type AbstractReaction end
export AbstractReaction

@with_kw struct ElementaryReaction{T<:AbstractRate,N<:AbstractSpecies,Q<:Integer} <: AbstractReaction
    index::Q
    reactants::Array{N,1}
    reactantinds::Array{Q,1}
    products::Array{N,1}
    productinds::Array{Q,1}
    kinetics::T
end
export ElementaryReaction

length(r::T) where {T<:AbstractReaction}= 1
export length
