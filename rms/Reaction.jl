using Parameters
include("./Calculators.jl")
include(".Species.jl")

abstract type AbstractReaction

@with_kw struct Reaction{T<:AbstractRate} <: AbstractReaction
    reactants::Array{Species,1}
    products::Array{Species,1}
    kinetics::T
end
