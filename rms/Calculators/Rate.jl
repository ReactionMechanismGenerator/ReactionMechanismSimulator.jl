using Parameters
using IterTools

include("../Tools.jl")
include("./RateUncertainty.jl")
include("../Constants.jl")

abstract type AbstractRate end

@with_kw struct Arrhenius{N,K,Q<:Number,P<:AbstractRateUncertainty} <: AbstractRate
        A::N
        n::K
        Ea::Q
        unc::P = EmptyRateUncertainty()
end
(arr::Arrhenius)(;T::Q,P::N=0.0,C::S=0.0) where {Q,N,S<:Number} = arr.A*T^arr.n*exp(-arr.Ea/(R*T))
(arr::Arrhenius)(T::Q;P::N=0.0,C::S=0.0) where {Q,N,S<:Number} = arr.A*T^arr.n*exp(-arr.Ea/(R*T))

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

@with_kw struct ThirdBody{N<:Integer,K<:AbstractFloat,Q<:AbstractRateUncertainty} <: AbstractRate
    arr::Arrhenius
    efficiencies::Dict{N,K} = Dict()
    unc::Q = EmptyRateUncertainty()
end

(tbarr::ThirdBody)(;T::Q=nothing,P::R=0.0,C::S=nothing) where {Q,R,S<:Number} = C*tbarr.arr(T=T)

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
