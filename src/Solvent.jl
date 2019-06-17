using Parameters

abstract type AbstractSolvent end
export AbstractSolvent

struct EmptySolvent <: AbstractSolvent end
export EmptySolvent

@with_kw struct Solvent <: AbstractSolvent
    name::String
    mu::AbstractViscosity
end
export Solvent
