using Parameters

abstract type AbstractRatevec end
export AbstractRatevec

@with_kw struct Arrheniusvec{N<:Array,K<:Array,Q<:Array} <: AbstractRatevec
        A::N
        n::K
        Ea::Q
end
function Arrheniusvec(arrs::T) where {T<:AbstractArray}
    A = zeros(length(arrs)) 
    n = zeros(length(arrs)) 
    Ea = zeros(length(arrs))
    for (i,aval) in enumerate(arrs)
        A[i] = aval.A
        n[i] = aval.n
        Ea[i] = aval.Ea
    end
    return Arrheniusvec(A=A,n=n,Ea=Ea)
end
@inline (arr::Arrheniusvec)(;T::Q,P::N=0.0,C::S=0.0) where {Q<:Real,N<:Real,S<:Real} = @fastmath @inbounds arr.A.*exp.(arr.n.*log(T).-arr.Ea.*(1.0/(R*T)))::Array{Q,1}
@inline (arr::Arrheniusvec)(T::Q;P::N=0.0,C::S=0.0) where {Q<:Real,N<:Real,S<:Real} = @fastmath @inbounds arr.A.*exp.(arr.n.*log(T).-arr.Ea.*(1.0/(R*T)))::Array{Q,1}
export Arrheniusvec

@with_kw struct Chebyshevvec{T<:AbstractArray,Q<:Real,S<:Real,V<:Real,B<:Real} <: AbstractRate
    coefs::T
    Tmin::Q
    Tmax::S
    Pmin::V
    Pmax::B
end

function Chebyshevvec(chevs::T) where {T<:AbstractArray}
    coefs = zeros(length(chevs),size(chevs[1].coefs)[1],size(chevs[1].coefs)[2])
    for k in 1:size(chevs[1].coefs)[2]
        for j in 1:size(chevs[1].coefs)[1]
            for i in 1:length(chevs)
                coefs[i,j,k] = chevs[i].coefs[j,k]
            end
        end
    end
    return Chebyshevvec(coefs=coefs,Tmin=chevs[1].Tmin,Tmax=chevs[1].Tmax,Pmin=chevs[1].Pmin,Pmax=chevs[1].Pmax)
end
export Chebyshevvec

@inline function evalChebyshevPolynomial(ch::Chebyshevvec,n::N,x::Q) where {N<:Integer,Q<:Real}
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

@inline function getredtemp(ch::Chebyshevvec,T::N) where {N<:Real}
    """
    return a reduced temperature corresponding to the given tempeprature
    for the Chebyshev polynomial, maps the inverse of temperautre onto [-1,1]
    """
    return @fastmath (2.0/T-1.0/ch.Tmin-1.0/ch.Tmax)/(1.0/ch.Tmax-1.0/ch.Tmin)
end
export getredtemp

@inline function getredpress(ch::Chebyshevvec,P::N) where {N<:Real}
    """
    return a reduced pressure corresponding to the given temperature
    for the Chebyshev polynomial maps the logarithm of pressure onto [-1,1]
    """
    return @fastmath (2.0*log10(P)-log10(ch.Pmin)-log10(ch.Pmax))/(log10(ch.Pmax)-log10(ch.Pmin))
end
export getredpress

@inline function (ch::Chebyshevvec)(;T::N,P::Q=0.0,C::B=0.0) where {N<:Real,B<:Real,Q<:Real}
    k = zeros(size(ch.coefs)[1])
    Tred = getredtemp(ch,T)
    Pred = getredpress(ch,P)
    klen,Tlen,Plen = size(ch.coefs)
    @simd for j = 1:Plen
        @simd for i = 1:Tlen
            @fastmath @inbounds k .+= ch.coefs[:,i,j].*(evalChebyshevPolynomial(ch,i-1,Tred)*evalChebyshevPolynomial(ch,j-1,Pred))
        end
    end
    return @fastmath 10.0.^k
end
