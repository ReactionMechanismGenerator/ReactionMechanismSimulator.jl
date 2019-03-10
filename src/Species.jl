using Parameters
import Base: length

abstract type AbstractSpecies end
export AbstractSpecies

struct EmptySpecies <: AbstractSpecies end

@with_kw struct Species{T<:AbstractThermo,N<:AbstractDiffusivity} <: AbstractSpecies
    name::String
    index::Integer
    inchi::String = ""
    smiles::String = ""
    thermo::T
    atomnums::Dict{String,Int64} = Dict()
    bondnum::Int64=-1
    diffusion::N = EmptyDiffusivity()
    radius::Float64 = 0.0
    radicalelectrons::Int64 = -100
end
export Species

length(r::T) where {T<:AbstractSpecies}= 1
export length
