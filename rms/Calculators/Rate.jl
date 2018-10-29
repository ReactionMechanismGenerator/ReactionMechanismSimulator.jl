using Parameters
using IterTools

include("../Tools.jl")
include("./RateUncertainty.jl")
include("../Constants.jl")

abstract type AbstractRate end
export AbstractRate

@with_kw struct Arrhenius{N,K,Q<:Number,P<:AbstractRateUncertainty} <: AbstractRate
        A::N
        n::K
        Ea::Q
        unc::P = EmptyRateUncertainty()
end
(arr::Arrhenius)(;T::Q,P::N=0.0,C::S=0.0) where {Q,N,S<:Number} = arr.A*T^arr.n*exp(-arr.Ea/(R*T))
(arr::Arrhenius)(T::Q;P::N=0.0,C::S=0.0) where {Q,N,S<:Number} = arr.A*T^arr.n*exp(-arr.Ea/(R*T))
export Arrhenius

@with_kw struct PdepArrhenius{T<:Number,Q<:AbstractRateUncertainty} <: AbstractRate
    Ps::Array{T,1}
    arrs::Array{Arrhenius,1}
    unc::Q = EmptyRateUncertainty()
end
PdepArrhenius(Ps::Array{Q,1},arrs::Array{Arrhenius,1}) where {Q<:Number} = PdepArrhenius(sort(Ps),arrs)

function (parr::PdepArrhenius)(;T::Q=nothing,P::R=nothing,C::S=0.0) where {Q,R,S<:Number}
    inds = getBoundingIndsSorted(P,parr.Ps)

    if length(inds) == 1
        return parr.arrs[inds[1]](T=T)
    else
        highk = parr.arrs[inds[2]](T=T)
        lowk = parr.arrs[inds[1]](T=T)
        Plow = parr.Ps[inds[1]]
        Phigh = parr.Ps[inds[2]]
        return lowk*10^(log10(P/Plow)/log10(Phigh/Plow)*log10(highk/lowk))
    end
end
export PdepArrhenius

@with_kw struct MultiArrhenius{Q<:AbstractRateUncertainty} <: AbstractRate
    arrs::Array{Arrhenius,1}
    unc::Q = EmptyRateUncertainty()
end

function (marr::MultiArrhenius)(;T::Q=nothing,P::R=0.0,C::S=0.0) where {Q,R,S<:Number}
    out = 0.0
    for arr in marr.arrs
        out += arr(T=T)
    end
    return out
end

@with_kw struct MultiPdepArrhenius{Q<:AbstractRateUncertainty} <: AbstractRate
    parrs::Array{PdepArrhenius,1}
    unc::Q = EmptyRateUncertainty()
end

function (parr::MultiPdepArrhenius)(;T::Q=nothing,P::R=0.0,C::S=0.0) where {Q,R,S<:Number}
    out = 0.0
    for pdar in parr.parrs
        out += pdar(T=T,P=P)
    end
    return out
end
export MultiPdepArrhenius

@with_kw struct ThirdBody{N<:Integer,K<:AbstractFloat,Q<:AbstractRateUncertainty} <: AbstractRate
    arr::Arrhenius
    efficiencies::Dict{N,K} = Dict()
    unc::Q = EmptyRateUncertainty()
end

(tbarr::ThirdBody)(;T::Q=nothing,P::R=0.0,C::S=nothing) where {Q,R,S<:Number} = C*tbarr.arr(T=T)
export ThirdBody

@with_kw struct Lindemann{N<:Integer,K<:AbstractFloat,Q<:AbstractRateUncertainty} <: AbstractRate
    arrhigh::Arrhenius
    arrlow::Arrhenius
    efficiencies::Dict{N,K} = Dict()
    unc::Q = EmptyRateUncertainty()
end

function (lnd::Lindemann)(;T::Q=nothing,P::R=0.0,C::S=nothing) where {Q,R,S<:Number}
    k0 = lnd.arrlow(T=T)
    kinf = lnd.arrhigh(T=T)
    Pr = k0*C/kinf
    return kinf*Pr/(1.0+Pr)
end
export Lindemann

@with_kw struct Troe{P,Q,F,L<:Number,N<:Integer,K<:AbstractFloat,R<:AbstractRateUncertainty} <: AbstractRate
    arrhigh::Arrhenius
    arrlow::Arrhenius
    a::P
    T3::Q
    T1::F
    T2::L
    efficiencies::Dict{N,K} = Dict()
    unc::R = EmptyRateUncertainty()
end

function (tr::Troe)(;T::Q=nothing,P::R=0.0,C::S=nothing) where {Q,R,S<:Number}
    k0 = tr.arrlow(T=T)
    kinf = tr.arrhigh(T=T)
    Pr = k0*C/kinf

    if tr.T1 == 0 && tr.T3 == 0
        F = 1.0
    else
        Fcent = (1-tr.a)*exp(-T/tr.T3)+tr.a*exp(-T/tr.T1)
        if tr.T2 != 0
            Fcent += exp(-tr.T2/T)
        end
        d = 0.14
        n = 0.75-1.27*log10(Fcent)
        c = -0.4-0.67*log10(Fcent)
        F = 10^(log10(Fcent)/(1+((log10(Pr)+c)/(n-d*(log10(Pr))))^2))
    end
    return kinf*(Pr/(1+Pr))*F
end
export Troe

@with_kw struct Chebyshev{T,Q,S,V,B<:Number}
    coefs::Array{T,2}
    Tmin::Q
    Tmax::S
    Pmin::V
    Pmax::B
end

function evalChebyshevPolynomial(ch::Chebyshev,n::N,x::Q) where {N<:Integer,Q<:Number}
    """
    evaluate the nth order Chebyshev Polynomial at x
    """
    if n==0
        return 1
    elseif n == 1
        return x
    else
        T0 = 1
        T1 = x
        for i in 1:(n-1)
            T = 2*x*T1-T0
            T0 = T1
            T1 = T
        end
        return T
    end
end

function getRedTemp(ch::Chebyshev,T::N) where {N<:Number}
    """
    return a reduced temperature corresponding to the given tempeprature
    for the Chebyshev polynomial, maps the inverse of temperautre onto [-1,1]
    """
    return (2.0/T-1.0/ch.Tmin-1.0/ch.Tmax)/(1.0/ch.Tmax-1.0/ch.Tmin)
end

function getRedPress(ch::Chebyshev,P::N) where {N<:Number}
    """
    return a reduced pressure corresponding to the given temperature
    for the Chebyshev polynomial maps the logarithm of pressure onto [-1,1]
    """
    return (2.0*log10(P)-log10(ch.Pmin)-log10(ch.Pmax))/(log10(ch.Pmax)-log10(ch.Pmin))
end

function (ch::Chebyshev)(;T::N,P::Q=0.0) where {N,Q<:Number}
    k = 0.0
    Tred = getRedTemp(ch,T)
    Pred = getRedPress(ch,P)
    Tlen,Plen = size(ch.coefs)
    for i = 1:Tlen
        for j = 1:Plen
            k += ch.coefs[i,j]*evalChebyshevPolynomial(ch,i-1,Tred)*evalChebyshevPolynomial(ch,j-1,Pred)
        end
    end
    return 10^k
end

export Chebyshev
