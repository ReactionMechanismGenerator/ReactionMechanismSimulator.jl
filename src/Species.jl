using Parameters
import Base: length

abstract type AbstractSpecies end
export AbstractSpecies

struct EmptySpecies <: AbstractSpecies end

@with_kw struct Species{T<:AbstractThermo,N<:AbstractDiffusivity,N1<:AbstractHenryLawConstant,N2<:AbstractLiquidVolumetricMassTransferCoefficient} <: AbstractSpecies
    name::String
    index::Integer
    inchi::String = ""
    smiles::String = ""
    adjlist::String = ""
    thermo::T
    atomnums::Dict{String,Int64} = Dict()
    bondnum::Int64=-1
    diffusion::N = EmptyDiffusivity()
    radius::Float64 = 0.0
    radicalelectrons::Int64 = -100
    molecularweight::Float64 = 0.0
    henrylawconstant::N1 = EmptyHenryLawConstant()
    liquidvolumetricmasstransfercoefficient::N2 = EmptyLiquidVolumetricMassTransferCoefficient()
    comment::String = ""
    isfragment::Bool = false
    isfragmentintermediate::Bool = false
end
export Species

length(r::T) where {T<:AbstractSpecies}= 1
export length
