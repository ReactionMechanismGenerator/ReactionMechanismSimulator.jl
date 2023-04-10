import Base: length

abstract type AbstractSpecies end
export AbstractSpecies

struct EmptySpecies <: AbstractSpecies end

struct Species{T<:AbstractThermo,N<:AbstractDiffusivity,N1<:AbstractHenryLawConstant,N2<:AbstractLiquidVolumetricMassTransferCoefficient} <: AbstractSpecies
    name::String
    index::Integer
    inchi::String
    smiles::String
    adjlist::String
    thermo::T
    atomnums::Dict
    bondnum::Int64
    diffusion::N
    radius::Float64
    radicalelectrons::Int64
    molecularweight::Float64
    henrylawconstant::N1
    liquidvolumetricmasstransfercoefficient::N2
    comment::String
    isfragment::Bool
    isfragmentintermediate::Bool
end

function Species(; name::String, index::Integer, inchi::String="", smiles::String="", adjlist::String="", 
                   thermo::T=nothing, atomnums::Dict=Dict(), bondnum::Int64=-1, diffusion::N=EmptyDiffusivity(),
                   radius::Float64=0.0, radicalelectrons::Int64=-100, molecularweight::Float64=0.0,
                   henrylawconstant::N1=EmptyHenryLawConstant(), 
                   liquidvolumetricmasstransfercoefficient::N2=EmptyLiquidVolumetricMassTransferCoefficient(),
                   comment::String="",isfragment::Bool=false,isfragmentintermediate::Bool=false,
                   ) where {T<:AbstractThermo,N<:AbstractDiffusivity,N1<:AbstractHenryLawConstant,N2<:AbstractLiquidVolumetricMassTransferCoefficient}
    return Species(name, index, inchi, smiles, adjlist, thermo, atomnums, bondnum, diffusion, radius, radicalelectrons, molecularweight, henrylawconstant, liquidvolumetricmasstransfercoefficient, comment, isfragment, isfragmentintermediate)
end
export Species

length(r::T) where {T<:AbstractSpecies}= 1
export length
