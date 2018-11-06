using Parameters
include("Calculators.jl")

abstract type AbstractSolvent end
export AbstractSolvent

struct Solvent <: AbstractSolvent
    name::String
    mu::AbstractViscosity
end
export Solvent
