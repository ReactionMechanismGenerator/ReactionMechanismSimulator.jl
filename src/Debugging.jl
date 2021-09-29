using Parameters
using Statistics

"""
Object for storing useful debugging information about potentially problematic species
tol is the tolerance analyzed at, G is the Gibbs free energy, dy is
the flux to that species and ratio is the ratio between dy and the flux scale
index is the index of the species
"""
@with_kw struct DebugSpecies
    name::String
    G::Float64 = 0.0
    ratio::Float64 = 0.0
    dy::Float64 = 0.0
    tol::Float64 = 0.0
    index::Int64 = 0
end

"""
Object for storing useful debugging information about potentially problematic reactions
tol is the tolerance analyzed at, kf and krev are the forward and reverse rate coefficients
rt is the rate of the reaction, ratio is the ratio between rt and the rate scale
T and P are the temperature and pressure and index is the index of the reaction
"""
@with_kw struct DebugReaction
    reactants::Array{DebugSpecies,1}
    products::Array{DebugSpecies,1}
    rxnstring::String
    tol::Float64 = 0.0
    kf::Float64 = 0.0
    krev::Float64 = 0.0
    rt::Float64 = 0.0
    ratio::Float64 = 0.0
    T::Float64 = 0.0
    P::Float64 = 0.0
    index::Int64 = 0
end

struct DebugMech
    spcs::Array{DebugSpecies,1}
    rxns::Array{DebugReaction,1}
end
