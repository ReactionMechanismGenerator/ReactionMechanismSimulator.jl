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


