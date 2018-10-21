using Parameters
include("Calculators.jl")

abstract type AbstractSpecies end
export AbstractSpecies

@with_kw struct Species{T<:AbstractThermo,N<:AbstractDiffusivity} <: AbstractSpecies
    name::String
    index::Integer
    inchi::String = ""
    smiles::String = ""
    thermo::T
    diffusion::N = EmptyDiffusivity()
end
export Species
