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
