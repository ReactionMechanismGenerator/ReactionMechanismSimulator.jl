using Parameters
include("Calculators.jl")
using Calculators

abstract type AbstractSolvent end
export AbstractSolvent

struct Solvent
    name::String
    mu::AbstractViscosity
end
export Solvent 
