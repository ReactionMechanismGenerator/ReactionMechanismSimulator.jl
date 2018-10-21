using Parameters
include("Species.jl")
include("Reaction.jl")
include("Solvent.jl")

abstract type AbstractPhase end
export AbstractPhase

struct EmptyPhase <: AbstractPhase end
export EmptyPhase

@with_kw struct IdealGas <: AbstractPhase
    name::String = ""
    species::Array{Species,1} = Array{Species,1}()
    reactions::Array{Reaction,1} = Array{Reaction,1}()
end
export IdealGas

@with_kw struct IdealDiluteSolution <: AbstractPhase
    name::String = ""
    species::Array{Species,1} = Array{Species,1}()
    reactions::Array{Reaction,1} = Array{Reaction,1}()
    solvent::Solvent
end
export IdealDiluteSolution

@with_kw struct HomogeneousCatalyst <: AbstractPhase
    name::String = ""
    species::Array{Species,1} = Array{Species,1}()
    reactions::Array{Reaction,1} = Array{Reaction,1}()
end
export HomogeneousCatalyst
