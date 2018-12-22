using Parameters
using IterTools

abstract type AbstractRate end
export AbstractRate

abstract type AbstractFalloffRate <: AbstractRate end
export AbstractFalloffRate

@with_kw struct Arrhenius{N<:Real,K<:Real,Q<:Real,P<:AbstractRateUncertainty} <: AbstractRate
        A::N
        n::K
        Ea::Q
        unc::P = EmptyRateUncertainty()
end
@inline (arr::Arrhenius)(;T::Q,P::N=0.0,C::S=0.0) where {Q<:Real,N<:Real,S<:Real} = @fastmath arr.A*T^arr.n*exp(-arr.Ea/(R*T))::Q
@inline (arr::Arrhenius)(T::Q;P::N=0.0,C::S=0.0) where {Q<:Real,N<:Real,S<:Real} = @fastmath arr.A*T^arr.n*exp(-arr.Ea/(R*T))::Q
export Arrhenius

@with_kw struct PdepArrhenius{T<:Real,Q<:AbstractRateUncertainty,Z<:AbstractRate} <: AbstractRate
    Ps::Array{T,1}
    arrs::Array{Z,1}
    unc::Q = EmptyRateUncertainty()
end
PdepArrhenius(Ps::Array{Q,1},arrs::Array{Z,1}) where {Q<:Real,Z<:AbstractRate} = PdepArrhenius(sort(Ps),arrs)

@inline function (parr::PdepArrhenius)(;T::Q=nothing,P::V=nothing,C::S=0.0) where {Q<:Real,V<:Real,S<:Real}
    inds = getBoundingIndsSorted(P,parr.Ps)::Tuple{Int64,Int64}

    if inds[2] == -1
        return @inbounds parr.arrs[inds[1]](T=T)::Q
    else
        @inbounds highk = parr.arrs[inds[2]](T=T)::Q
        @inbounds lowk = parr.arrs[inds[1]](T=T)::Q
        @inbounds Plow = parr.Ps[inds[1]]
        @inbounds Phigh = parr.Ps[inds[2]]
        return @inbounds @fastmath lowk*10^(log10(P/Plow)/log10(Phigh/Plow)*log10(highk/lowk))
    end
end
export PdepArrhenius

@with_kw struct MultiArrhenius{Q<:AbstractRateUncertainty} <: AbstractRate
    arrs::Array{Arrhenius,1}
    unc::Q = EmptyRateUncertainty()
end

@inline function (marr::MultiArrhenius)(;T::Q=nothing,P::R=0.0,C::S=0.0) where {Q<:Real,R<:Real,S<:Real}
    out = 0.0
    for arr in marr.arrs
        @fastmath out += arr(T)::Q
    end
    return out
end
export MultiArrhenius

@with_kw struct MultiPdepArrhenius{Q<:AbstractRateUncertainty} <: AbstractRate
    parrs::Array{PdepArrhenius,1}
    unc::Q = EmptyRateUncertainty()
end

@inline function (parr::MultiPdepArrhenius)(;T::Q=nothing,P::R=0.0,C::S=0.0) where {Q<:Real,R<:Real,S<:Real}
    out = 0.0
    for pdar in parr.parrs
        @fastmath out += pdar(T=T,P=P)::Q
    end
    return out
end
export MultiPdepArrhenius

@with_kw struct ThirdBody{N<:Integer,K<:AbstractFloat,Q<:AbstractRateUncertainty} <: AbstractFalloffRate
    arr::Arrhenius
    efficiencies::Dict{N,K} = Dict()
    unc::Q = EmptyRateUncertainty()
end

(tbarr::ThirdBody)(;T::Q=nothing,P::R=0.0,C::S=nothing) where {Q<:Real,R<:Real,S<:Real} = C*(tbarr.arr(T)::Q)
export ThirdBody

@with_kw struct Lindemann{N<:Integer,K<:AbstractFloat,Q<:AbstractRateUncertainty} <: AbstractFalloffRate
    arrhigh::Arrhenius
    arrlow::Arrhenius
    efficiencies::Dict{N,K} = Dict()
    unc::Q = EmptyRateUncertainty()
end

@inline function (lnd::Lindemann)(;T::Q=nothing,P::R=0.0,C::S=nothing) where {Q<:Real,R<:Real,S<:Real}
    k0 = lnd.arrlow(T=T)::Q
    kinf = lnd.arrhigh(T=T)::Q
    @fastmath Pr = k0*C/kinf
    return @fastmath kinf*Pr/(1.0+Pr)
end
export Lindemann

@with_kw struct Troe{P<:Real,Q<:Real,F<:Real,L<:Real,N<:Integer,K<:AbstractFloat,R<:AbstractRateUncertainty} <: AbstractFalloffRate
    arrhigh::Arrhenius
    arrlow::Arrhenius
    a::P
    T3::Q
    T1::F
    T2::L
    efficiencies::Dict{N,K} = Dict()
    unc::R = EmptyRateUncertainty()
end

@inline function (tr::Troe)(;T::Q=nothing,P::R=0.0,C::S=nothing) where {Q<:Real,R<:Real,S<:Real}
    k0 = tr.arrlow(T=T)::Q
    kinf = tr.arrhigh(T=T)::Q
    @fastmath Pr = k0*C/kinf

    if tr.T1 == 0.0 && tr.T3 == 0.0
        F = 1.0
    else
        @fastmath Fcent = (1-tr.a)*exp(-T/tr.T3)+tr.a*exp(-T/tr.T1)
        if tr.T2 != 0.0
            @fastmath Fcent += exp(-tr.T2/T)
        end
        d = 0.14
        @fastmath n = 0.75-1.27*log10(Fcent)
        @fastmath c = -0.4-0.67*log10(Fcent)
        @fastmath F = 10^(log10(Fcent)/(1+((log10(Pr)+c)/(n-d*(log10(Pr))))^2))
    end
    return @fastmath kinf*(Pr/(1+Pr))*F::S
end
export Troe

@with_kw struct Chebyshev{T<:Real,Q<:Real,S<:Real,V<:Real,B<:Real} <: AbstractRate
    coefs::Array{T,2}
    Tmin::Q
    Tmax::S
    Pmin::V
    Pmax::B
end

@inline function evalChebyshevPolynomial(ch::Chebyshev,n::N,x::Q) where {N<:Integer,Q<:Real}
    """
    evaluate the nth order Chebyshev Polynomial at x
    """
    if n==0
        return 1.0
    elseif n == 1
        return x
    else
        T = 0.0
        T0 = 1.0
        T1 = x
        for i in 1:(n-1)
            @fastmath T = 2*x*T1-T0
            T0 = T1
            T1 = T
        end
        return T
    end
end
export evalChebyshevPolynomial

@inline function getredtemp(ch::Chebyshev,T::N) where {N<:Real}
    """
    return a reduced temperature corresponding to the given tempeprature
    for the Chebyshev polynomial, maps the inverse of temperautre onto [-1,1]
    """
    return @fastmath (2.0/T-1.0/ch.Tmin-1.0/ch.Tmax)/(1.0/ch.Tmax-1.0/ch.Tmin)
end
export getredtemp

@inline function getredpress(ch::Chebyshev,P::N) where {N<:Real}
    """
    return a reduced pressure corresponding to the given temperature
    for the Chebyshev polynomial maps the logarithm of pressure onto [-1,1]
    """
    return @fastmath (2.0*log10(P)-log10(ch.Pmin)-log10(ch.Pmax))/(log10(ch.Pmax)-log10(ch.Pmin))
end
export getredpress

@inline function (ch::Chebyshev)(;T::N,P::Q=0.0,C::B=0.0) where {N<:Real,B<:Real,Q<:Real}
    k = 0.0
    Tred = getredtemp(ch,T)
    Pred = getredpress(ch,P)
    Tlen,Plen = size(ch.coefs)
    @simd for i = 1:Tlen
        @simd for j = 1:Plen
            @fastmath @inbounds k += ch.coefs[i,j]*evalChebyshevPolynomial(ch,i-1,Tred)*evalChebyshevPolynomial(ch,j-1,Pred)
        end
    end
    return @fastmath 10.0^k
end


export Chebyshev
