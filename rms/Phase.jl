using Parameters
import Base: length

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
    spcdict::Dict{String,Int64}
    diffusionlimited::Bool = true
end
IdealDiluteSolution(species,reactions,solvent; name="",diffusionlimited=true) = IdealDiluteSolution(species=species,reactions=reactions,
solvent=solvent,name=name,diffusionlimited=diffusionlimited,spcdict=Dict([sp.name=>sp.index for sp in species]))
export IdealDiluteSolution

@with_kw struct HomogeneousCatalyst{Q<:AbstractReaction} <: AbstractPhase
    name::String = ""
    species::Array{Species,1}
    reactions::Array{Q,1}
    spcdict::Dict{String,Int64}
end
export HomogeneousCatalyst

length(p::T) where {T<:AbstractPhase} = 1
export length

iterate(p::T) where {T<:AbstractPhase} = p
export iterate

Broadcast.broadcastable(p::T) where {T<:AbstractPhase} = Ref(p)
export broadcastable
