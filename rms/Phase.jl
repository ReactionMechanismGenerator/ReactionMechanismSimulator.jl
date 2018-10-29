using Parameters
import Base: length
include("Constants.jl")
include("Species.jl")
include("Reaction.jl")
include("Solvent.jl")

abstract type AbstractPhase end
export AbstractPhase

abstract type IdealPhase <: AbstractPhase end
export IdealPhase

struct EmptyPhase <: AbstractPhase end
export EmptyPhase

@with_kw struct IdealGas{Q<:AbstractReaction} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{Q,1}
    spcdict::Dict{String,Int64}
    diffusionlimited::Bool = false
end
IdealGas(species,reactions; name="",diffusionlimited=false) = IdealGas(species=species,reactions=reactions,name=name,
diffusionlimited=diffusionlimited,spcdict=Dict([sp.name=>sp.index for sp in species]))
export IdealGas

@with_kw struct IdealDiluteSolution{Q<:AbstractReaction} <: IdealPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{Q,1}
    solvent::Solvent
    diffusionlimited::Bool = true
end
IdealDiluteSolution(species,reactions; name="",diffusionlimited=true) = IdealDiluteSolution(species=species,reactions=reactions,name=name,
diffusionlimited=diffusionlimited,spcdict=Dict([sp.name=>sp.index for sp in species]))
export IdealDiluteSolution

@with_kw struct HomogeneousCatalyst{Q<:AbstractReaction} <: AbstractPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{Q,1}
end
export HomogeneousCatalyst

length(p::T) where {T<:AbstractPhase} = 1
export length
