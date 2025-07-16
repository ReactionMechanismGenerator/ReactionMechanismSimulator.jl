using Parameters

abstract type AbstractThermo end
export AbstractThermo

@inline function getGibbs(th::P,T::N;coverages=nothing) where {N<:Number,P}
    return @fastmath (getEnthalpy(th,T)+getcovdepenthalpycorrection(th.covdep,T,coverages))-T*(getEntropy(th,T)+getcovdepentropycorrection(th.covdep,T,coverages))
end

@with_kw struct NASApolynomial{H<:AbstractThermoCoverageDependence} <: AbstractThermo
    coefs::Array{Float64,1}
    Tmin::Float64
    Tmax::Float64
    covdep::H = EmptyThermoCoverageDependence()
end
export NASApolynomial

@inline function calcHSCpdless(poly::NASApolynomial,T::X;coverages=nothing) where {X<:Real}
    if length(poly.coefs) != 7
        Tpoly0 = T
        Tpoly1 = T*T
        Tpoly2 = Tpoly1*T
        Tpoly3 = Tpoly2*T
        Tpoly4 = 1.0/T
        Tpoly5 = Tpoly4*Tpoly4
        Tpoly6 = log(T)

        ct0 = poly.coefs[1]*Tpoly5
        ct1 = poly.coefs[2]*Tpoly4
        ct2 = poly.coefs[3]
        ct3 = poly.coefs[4]*Tpoly0
        ct4 = poly.coefs[5]*Tpoly1
        ct5 = poly.coefs[6]*Tpoly2
        ct6 = poly.coefs[7]*Tpoly3

        cpdivR = ct0 + ct1 + ct2 + ct3 + ct4 + ct5 + ct6

        hdivRT = -ct0+Tpoly6*ct1+ct2+0.5*ct3+0.33333333333*ct4+0.25*ct5+0.2*ct6+poly.coefs[8]*Tpoly4
        sdivR = -0.5*ct0 - ct1 + Tpoly6*ct2 + ct3 + 0.5*ct4 + 0.33333333333*ct5 + 0.25*ct6 + poly.coefs[9]
    else
        Tpoly0 = T
        Tpoly1 = T*T
        Tpoly2 = Tpoly1*T
        Tpoly3 = Tpoly2*T
        Tpoly4 = 1.0/T
        #Tpoly5 = Tpoly4*Tpoly4
        Tpoly6 = log(T)

        #ct0 = 0.0
        #ct1 = 0.0
        ct2 = poly.coefs[1]
        ct3 = poly.coefs[2]*Tpoly0
        ct4 = poly.coefs[3]*Tpoly1
        ct5 = poly.coefs[4]*Tpoly2
        ct6 = poly.coefs[5]*Tpoly3

        cpdivR = ct2 + ct3 + ct4 + ct5 + ct6

        hdivRT =ct2+0.5*ct3+0.33333333333*ct4+0.25*ct5+0.2*ct6+poly.coefs[6]*Tpoly4
        sdivR = Tpoly6*ct2 + ct3 + 0.5*ct4 + 0.33333333333*ct5 + 0.25*ct6 + poly.coefs[7]
    end
    hdivRT += getcovdepenthalpycorrection(th.covdep,T,coverages)/(R*T)
    sdivR += getcovdepentropycorrection(th.covdep,T,coverages)/R
    cpdivR += getcovdepheatcapacitycorrection(th.covdep,T,coverages)/R
    return (cpdivR,hdivRT,sdivR)::Tuple{Float64,Float64,Float64}
end
export calcHSCpdless

@inline function getHeatCapacity(poly::NASApolynomial,T::N,coverages=nothing) where {N<:Number}
    if length(poly.coefs) == 9
        return @views @inbounds @fastmath evalpoly(T,poly.coefs[1:7])/T^2*R + getcovdepheatcapacitycorrection(poly.covdep,T,coverages)
    elseif length(poly.coefs) == 7
        return @views @inbounds @fastmath evalpoly(T,poly.coefs[1:5])*R + getcovdepheatcapacitycorrection(poly.covdep,T,coverages)
    else
        throw(error("NASA polynomial has a number of coefficients not equal to 9 or 7"))
    end
end

@inline function getEntropy(poly::NASApolynomial,T::N,coverages=nothing) where {N<:Number}
    if length(poly.coefs) == 9
        return @views @inbounds @fastmath ((-poly.coefs[1]/(2*T)-poly.coefs[2])/T+poly.coefs[3]*log(T)+T*evalpoly(T,poly.coefs[4:end-2]./(1:4))+poly.coefs[end])*R + getcovdepentropycorrection(poly.covdep,T,coverages)
    elseif length(poly.coefs) == 7
        return @views @inbounds @fastmath (poly.coefs[1]*log(T)+T*evalpoly(T,poly.coefs[2:end-2]./(1:4))+poly.coefs[end])*R + getcovdepentropycorrection(poly.covdep,T,coverages)
    else
        throw(error("NASA polynomial has a number of coefficients not equal to 9 or 7"))
    end
end

@inline function getEnthalpy(poly::NASApolynomial,T::N,coverages=nothing) where {N<:Number}
    if length(poly.coefs) == 9
        return @views @inbounds @fastmath ((-poly.coefs[1]/T+poly.coefs[2]*log(T))/T+evalpoly(T,poly.coefs[3:end-2]./(1:5)))*R*T+poly.coefs[end-1]*R+th.covdep(T,coverages) + getcovdepenthalpycorrection(poly.covdep,T,coverages)
    elseif length(poly.coefs) == 7
        return @views @inbounds @fastmath evalpoly(T,poly.coefs[1:end-2]./(1:5))*R*T+poly.coefs[end-1]*R + getcovdepenthalpycorrection(poly.covdep,T,coverages)
    else
        throw(error("NASA polynomial has a number of coefficients not equal to 9 or 7"))
    end
end

@with_kw struct NASA{T<:AbstractThermoUncertainty,H<:AbstractThermoCoverageDependence} <: AbstractThermo
    polys::Array{NASApolynomial,1}
    unc::T = EmptyThermoUncertainty()
    covdep::H = EmptyThermoCoverageDependence()
end
function NASA(polys,unc=EmptyThermoUncertainty(),covdep=EmptyThermoCoverageDependence())
    ps = [polys.covdep isa EmptyThermoCoverageDependence ? NASApolynomial(p.coefs,p.Tmin,p.Tmax,covdep) : p for p in polys]
    return NASA(ps,unc,covdep)
end
export NASA

@inline function selectPoly(nasa::NASA,T::N) where {N<:Real}
    """
    retrieve the nasa polynomial corresponding to the T range
    """
    for p in nasa.polys
        if T<=p.Tmax
            return p
        end
    end
    return nasa.polys[end]
end
export selectPoly

@inline getHeatCapacity(nasa::NASA,T::N) where {N<:Number} = getHeatCapacity(selectPoly(nasa,T),T) #we do not add covdep corrections in NASA because it is also added in the NASAPolynomial, so it would double count
@inline getEntropy(nasa::NASA,T::N) where {N<:Number} = getEntropy(selectPoly(nasa,T),T)
@inline getEnthalpy(nasa::NASA,T::N) where {N<:Number} = getEnthalpy(selectPoly(nasa,T),T)
@inline getGibbs(nasa::NASA,T::N) where {N<:Number} = getGibbs(selectPoly(nasa,T),T)
@inline calcHSCpdless(nasa::NASA,T::N) where {N<:Real}= calcHSCpdless(selectPoly(nasa,T),T)

@with_kw struct Wilhoit{N,Q,T,P,U,R<:Number,M<:AbstractThermoUncertainty,H<:AbstractThermoCoverageDependence} <: AbstractThermo
    Cp0::N
    Cpinf::T
    coefs::Array{Q,1}
    H0::P
    S0::U
    B::R
    unc::M = EmptyThermoUncertainty()
    covdep::H = EmptyThermoCoverageDependence()
end
export Wilhoit

@inline function getHeatCapacity(w::Wilhoit,T::N,coverages=nothing) where {N<:Number}
    @fastmath y = T/(T+w.B)
    return @fastmath w.Cp0 + (w.Cpinf-w.Cp0)*y^2*(1+(y-1)*evalpoly(y,w.coefs)) + getcovdepheatcapacitycorrection(th.covdep,T,coverages)
end

@inline function getEnthalpy(w::Wilhoit,T::N,coverages=nothing) where {N<:Number}
    @fastmath y = T/(T+w.B)
    return @views @fastmath @inbounds w.H0 + w.Cp0 * T - (w.Cpinf - w.Cp0) * T * (
            y * y * ((3 * w.coefs[1] + sum(w.coefs[2:end])) / 6. +
                     (4 * w.coefs[2] + sum(w.coefs[3:end])) * y / 12. +
                     (5 * w.coefs[3] + w.coefs[4]) * y^2 / 20. +
                     w.coefs[4] * y^3 / 5.) +
            (2 + sum(w.coefs)) * (y / 2. - 1 + (1.0 / y - 1.) * log(w.B + T))
        ) + getcovdepenthalpycorrection(th.covdep,T,coverages)
end

@inline function getEntropy(w::Wilhoit,T::N,coverages=nothing) where {N<:Number}
    @fastmath y = T/(T+w.B)
    return @fastmath w.S0 + w.Cpinf*log(T)-(w.Cpinf-w.Cp0)*(log(y)+y*(1+y*evalpoly(y,w.coefs./(2:5)))) + getcovdepentropycorrection(th.covdep,T,coverages)
end

@with_kw struct ConstantG{B<:Number,J} <: AbstractThermo
    G::B
    T::J
end
export ConstantG

@inline function getGibbs(cg::ConstantG,T::B,coverages=nothing) where {B<:Number}
    return cg.G
end

export getGibbs, getEntropy, getEnthalpy, getHeatCapacity
