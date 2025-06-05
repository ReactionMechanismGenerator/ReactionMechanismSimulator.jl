using Parameters

abstract type AbstractThermoCoverageDependence end
export AbstractThermoCoverageDependence

struct EmptyThermoCoverageDependence <: AbstractThermoCoverageDependence end
@inline getcovdepenthalpycorrection(covdep::EmptyThermoCoverageDependence,T,coverages) = zero(T)
@inline getcovdepentropycorrection(covdep::EmptyThermoCoverageDependence,T,coverages) = zero(T)
@inline getcovdepheatcapacitycorrection(covdep::EmptyThermoCoverageDependence,T,coverages) = zero(T)
export EmptyThermoCoverageDependence

"""
Polynomial energy coverage dependence
correction is applied to enthalpy
poly is a dictionary mapping species names to polynomial coverage dependence corrections to their energy/enthalpy
polynomial ordering is a_1 + a_2 * theta + a_3 * theta^2 etc.
"""
@with_kw struct PolynomialThermoEnergyCoverageDependence <: AbstractThermoCoverageDependence
    polys::Dict{String,Vector{Float64}} = Dict{String,Vector{Float64}}()
    indpolys::Dict{Int64,Vector{Float64}} = Dict{Int64,Vector{Float64}}()
end
@inline function PolynomialThermoEnergyCoverageDependence(polys)
    return PolynomialThermoEnergyCoverageDependence(polys,Dict{Int64,Vector{Float64}}())
end
@inline function getcovdepenthalpycorrection(covdep::PolynomialThermoEnergyCoverageDependence,T,coverages)
    E = zero(coverages[1])
    for (ind,poly) in covdep.indpolys
        E += evalpoly(coverages[ind],poly)
    end
    return E
end
@inline getcovdepentropycorrection(covdep::PolynomialThermoEnergyCoverageDependence,T,coverages) = zero(coverages[1])
@inline getcovdepheatcapacitycorrection(covdep::PolynomialThermoEnergyCoverageDependence,T,coverages) = zero(coverages[1])
export PolynomialThermoEnergyCoverageDependence