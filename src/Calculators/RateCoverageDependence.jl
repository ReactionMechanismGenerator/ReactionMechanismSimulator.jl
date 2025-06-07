using Parameters

abstract type AbstractRateCoverageDependence end
export AbstractRateCoverageDependence

struct EmptyRateCoverageDependence <: AbstractRateCoverageDependence end
export EmptyRateCoverageDependence

getcovdepactivationbarriercorrection(covdep::EmptyRateCoverageDependence,T,coverages) = 0.0
getcovdepfactorcorrection(covdep::EmptyRateCoverageDependence,T,coverages) = 1.0

"""
Polynomial rate coefficient coverage dependence
correction form is: 10^(sum(a_i*theta_k^i,i)) * theta_k^m * exp(sum(E_i*theta_k^i,i)/(RT))
dictionaries map the name of the species (k) to polynomials with ordering a_1 + a_2 * theta + a_3 * theta^2 etc.
Epolys corresponds to E_i (array of polynomial coefficients), avals corresponds to a_i (array of polynomial coefficients), 
ms corresponds to m (float)
"""
@with_kw struct PolynomialRateCoverageDependence <: AbstractRateCoverageDependence
    avals::Dict{String,Vector{Float64}} = Dict{String,Vector{Float64}}()
    Epolys::Dict{String,Vector{Float64}} = Dict{String,Vector{Float64}}()
    ms::Dict{String,Vector{Float64}} = Dict{String,Float64}()
    indavals::Dict{Int64,Vector{Float64}} = Dict{Int64,Vector{Float64}}()
    indEpolys::Dict{Int64,Vector{Float64}} = Dict{Int64,Vector{Float64}}()
    indms::Dict{Int64,Vector{Float64}} = Dict{Int64,Float64}()
end
@inline function PolynomialRateCoverageDependence(Epolys,ms=Dict{String,Float64}(),avals=Dict{String,Vector{Float64}}())
    return PolynomialRateCoverageDependence(;Epolys=Epolys,
            ms=ms,avals=avals,indavals=Dict{Int64,Vector{Float64}}(),indEpolys=Dict{Int64,Vector{Float64}}(),
            indms=Dict{Int64,Float64}())
end
@inline function getcovdepactivationbarriercorrection(covdep::PolynomialRateCoverageDependence,T,coverages)
    if coverages === nothing
        return 0.0
    elseif length(coverages) == 0
        return zero(typeof(coverages).parameters[1])
    end
    E = 0.0
    
    for (ind,Epoly) in covdep.indEpolys
        E += evalpoly(coverages[ind],Epoly)
    end
    return E
end
export getcovdepactivationbarriercorrection
@inline function getcovdepfactorcorrection(covdep::PolynomialRateCoverageDependence,T,coverages)
    if coverages === nothing
        return 1.0
    elseif length(coverages) == 0
        return one(typeof(coverages).parameters[1])
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

export getcovdepactivationbarriercorrection
export getcovdepfactorcorrection