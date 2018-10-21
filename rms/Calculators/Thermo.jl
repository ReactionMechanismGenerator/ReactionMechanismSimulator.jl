using Parameters

include("../Constants.jl")
include("../Tools.jl") #evalpoly comes from here
include("ThermoUncertainty.jl")

abstract type AbstractThermo end

function getGibbs(th::P,T::N) where {N<:Number,P<:AbstractThermo}
    return getEnthalpy(th,T)-T*getEntropy(th,T)
end

@with_kw struct NASApolynomial{N,Q,P<:Number} <: AbstractThermo
    coefs::Array{N,1}
    Tmin::Q
    Tmax::P
end

function getHeatCapacity(poly::NASApolynomial,T::N) where {N<:Number}
    if length(poly.coefs) == 9
        return evalpoly(T,poly.coefs[1:7])/T^2*R
    elseif length(poly.coefs) == 7
        return evalpoly(T,poly.coefs[1:5])*R
    else
        throw(error("NASA polynomial has a number of coefficients not equal to 9 or 7"))
    end
end

function getEntropy(poly::NASApolynomial,T::N) where {N<:Number}
    if length(poly.coefs) == 9
        return ((-poly.coefs[1]/(2*T)-poly.coefs[2])/T+poly.coefs[3]*log(T)+T*evalpoly(T,poly.coefs[4:end-2]./(1:4))+poly.coefs[end])*R
    elseif length(poly.coefs) == 7
        return (poly.coefs[1]*log(T)+T*evalpoly(T,poly.coefs[2:end-2]./(1:4))+poly.coefs[end])*R
    else
        throw(error("NASA polynomial has a number of coefficients not equal to 9 or 7"))
    end
end

function getEnthalpy(poly::NASApolynomial,T::N) where {N<:Number}
    if length(poly.coefs) == 9
        return ((-poly.coefs[1]/T+poly.coefs[2]*log(T))/T+evalpoly(T,poly.coefs[3:end-2]./(1:5)))*R*T+poly.coefs[end-1]*R
    elseif length(poly.coefs) == 7
        return evalpoly(T,poly.coefs[1:end-2]./(1:5))*R*T+poly.coefs[end-1]*R
    else
        throw(error("NASA polynomial has a number of coefficients not equal to 9 or 7"))
    end
end

@with_kw struct NASA{T<:AbstractThermoUncertainty} <: AbstractThermo
    polys::Array{NASApolynomial,1}
    unc::T = EmptyThermoUncertainty()
end

function selectPoly(nasa::NASA,T::N) where {N<:Number}
    """
    retrieve the nasa polynomial corresponding to the T range
    """
    for p in nasa.polys
        if T<p.Tmax
            if T>p.Tmin
                return p
            end
        end
    end
    throw(error(String("No valid NASA polynomial at T=$T")))
end

getHeatCapacity(nasa::NASA,T::N) where {N<:Number} = getHeatCapacity(selectPoly(nasa,T),T)
getEntropy(nasa::NASA,T::N) where {N<:Number} = getEntropy(selectPoly(nasa,T),T)
getEnthalpy(nasa::NASA,T::N) where {N<:Number} = getEnthalpy(selectPoly(nasa,T),T)
getGibbs(nasa::NASA,T::N) where {N<:Number} = getGibbs(selectPoly(nasa,T),T)

@with_kw struct Wilhoit{N,Q,T,P,U,R<:Number,M<:AbstractThermoUncertainty} <: AbstractThermo
    Cp0::N
    Cpinf::T
    coefs::Array{Q,1}
    H0::P
    S0::U
    B::R
    unc::M = EmptyThermoUncertainty()
end

function getHeatCapacity(w::Wilhoit,T::N) where {N<:Number}
    y = T/(T+w.B)
    return w.Cp0 + (w.Cpinf-w.Cp0)*y^2*(1+(y-1)*evalpoly(y,w.coefs))
end

function getEnthalpy(w::Wilhoit,T::N) where {N<:Number}
    y = T/(T+w.B)
    return w.H0 + w.Cp0 * T - (w.Cpinf - w.Cp0) * T * (
            y * y * ((3 * w.coefs[1] + sum(w.coefs[2:end])) / 6. +
                     (4 * w.coefs[2] + sum(w.coefs[3:end])) * y / 12. +
                     (5 * w.coefs[3] + w.coefs[4]) * y^2 / 20. +
                     w.coefs[4] * y^3 / 5.) +
            (2 + sum(w.coefs)) * (y / 2. - 1 + (1.0 / y - 1.) * log(w.B + T))
        )
end

function getEntropy(w::Wilhoit,T::N) where {N<:Number}
    y = T/(T+w.B)
    return w.S0 + w.Cpinf*log(T)-(w.Cpinf-w.Cp0)*(log(y)+y*(1+y*evalpoly(y,w.coefs./(2:5))))
end
