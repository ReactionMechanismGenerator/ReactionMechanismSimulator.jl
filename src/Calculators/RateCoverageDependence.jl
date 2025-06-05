using Parameters

abstract type AbstractRateCoverageDependence end
export AbstractRateCoverageDependence

struct EmptyRateCoverageDependence <: AbstractRateCoverageDependence end
export EmptyRateCoverageDependence

getcovdepactivationbarriercorrection(covdep::EmptyRateCoverageDependence,T,coverages) = 0.0
getcovdepfactorcorrection(covdep::EmptyRateCoverageDependence,T,coverages) = 1.0

export getcovdepactivationbarriercorrection
export getcovdepfactorcorrection

"""
Polynomial rate coefficient coverage dependence
correction form is: 10^(sum(a_i*theta_k^i,i)) * theta_k^m * exp(sum(E_i*theta_k^i,i)/(RT))
dictionaries map the name of the species (k) to polynomials with ordering a_1 + a_2 * theta + a_3 * theta^2 etc.
Epolys corresponds to E_i (array of polynomial coefficients), avals corresponds to a_i (array of polynomial coefficients), 
ms corresponds to m (float)
"""
@with_kw struct PolynomialRateCoverageDependence{N,K,L} <: AbstractRateCoverageDependence
    avals::Dict{String,L} = Dict()
    Epolys::Dict{String,K} = Dict()
    ms::Dict{String,L} = Dict()
    indavals::Dict{N,L} = Dict()
    indEpolys::Dict{N,K} = Dict()
    indms::Dict{N,L} = Dict()
end
@inline function getcovdepactivationbarriercorrection(covdep::PolynomialRateCoverageDependence,T,coverages)
    if length(coverages) == 0
        return zero(typeof(coverages).parameters[1])
    end
    E = 0.0
    for (ind,Epoly) in covdep.indEpolys
        E += evalpoly(coverages[ind],indEpoly)
    end
    return E
end
@inline function getcovdepfactorcorrection(covdep::PolynomialRateCoverageDependence,T,coverages)
    if length(coverages) == 0
        return zero(typeof(coverages).parameters[1])
    end
    av = 0.0
    for (ind,a) in covdep.indavals
       av += coverages[ind]*a
    end

    if av != 0
        A = 10.0^av
    else 
        A = 1.0
    end

    for (ind,m) in covdep.indms 
        A *= coverages[ind]^m
    end
    
    return A
end
export PolynomialRateCoverageDependence