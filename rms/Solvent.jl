using Parameters
include("Calculators.jl")

abstract type AbstractSolvent end
export AbstractSolvent

struct EmptySolvent <: AbstractSolvent end
export EmptySolvent

struct Solvent <: AbstractSolvent
    name::String
    mu::AbstractViscosity
end
export Solvent
