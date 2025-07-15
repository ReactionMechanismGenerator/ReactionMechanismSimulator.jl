using Parameters
using ReverseDiff
using ForwardDiff

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
@inline (arr::Arrheniusvec)(;T::Q,P::N=0.0,C::S=0.0,phi=0.0,dGrxns=0.0) where {Q<:Real,N<:Real,S<:Real} = @fastmath @inbounds arr.A.*exp.(arr.n.*log(T).-arr.Ea.*(1.0/(R*T)))
@inline (arr::Arrheniusvec)(T::Q;P::N=0.0,C::S=0.0,phi=0.0) where {Q<:Real,N<:Real,S<:Real} = @fastmath @inbounds arr.A.*exp.(arr.n.*log(T).-arr.Ea.*(1.0/(R*T)))
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

@inline function (ch::Chebyshevvec)(;T::N,P::Q=0.0,C::B=0.0,phi=0.0,dGrxns=0.0) where {N<:Real,B<:Real,Q<:Real}
    k = zeros(N,size(ch.coefs)[1])
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

@with_kw struct Troevec{P<:AbstractArray,Q<:AbstractArray,F<:AbstractArray,L<:AbstractArray,N<:Integer,K<:AbstractFloat,R<:AbstractRateUncertainty} <: AbstractFalloffRate
    arrhigh::Arrheniusvec
    arrlow::Arrheniusvec
    a::P
    T3::Q
    T1::F
    T2::L
    efficiencies::Array{Dict{N,K},1}
    unc::R = EmptyRateUncertainty()
end

function Troevec(troes::T) where {T<:AbstractArray}
    
    Ahigh = zeros(length(troes))
    nhigh = zeros(length(troes))
    Ehigh = zeros(length(troes))
    Alow = zeros(length(troes))
    nlow = zeros(length(troes))
    Elow = zeros(length(troes))
    a = zeros(length(troes))
    T3 = zeros(length(troes))
    T1 = zeros(length(troes))
    T2 = zeros(length(troes))
    efficiencies = [Dict{Int64,Float64}() for x in troes]
    for (i,falloff) in enumerate(troes)
        if isa(falloff, ThirdBody)
            Alow[i] = falloff.arr.A
            nlow[i] = falloff.arr.n
            Elow[i] = falloff.arr.Ea
            Ahigh[i] = 1e20  #very high limited by energy transfer process
            T3[i] = Inf
            T1[i] = Inf
            if length(efficiencies) > 0
                efficiencies[i] = falloff.efficiencies
            end
        elseif isa(falloff,Lindemann)
            Alow[i] = falloff.arrlow.A
            nlow[i] = falloff.arrlow.n
            Elow[i] = falloff.arrlow.Ea
            Ahigh[i] = falloff.arrhigh.A
            nhigh[i] = falloff.arrhigh.n
            Ehigh[i] = falloff.arrhigh.Ea
            T3[i] = Inf
            T1[i] = Inf
            if length(efficiencies) > 0
                efficiencies[i] = falloff.efficiencies
            end
        elseif isa(falloff,Troe)
            Alow[i] = falloff.arrlow.A
            nlow[i] = falloff.arrlow.n
            Elow[i] = falloff.arrlow.Ea
            Ahigh[i] = falloff.arrhigh.A
            nhigh[i] = falloff.arrhigh.n
            Ehigh[i] = falloff.arrhigh.Ea
            a[i] = falloff.a
            T3[i] = falloff.T3
            T1[i] = falloff.T1
            T2[i] = falloff.T2
            if length(efficiencies) > 0
                efficiencies[i] = falloff.efficiencies
            end
        else
            val = typeof(falloff)
            error("could not process type $val in Troevec creation")
        end
    end
    return Troevec(arrhigh=Arrheniusvec(Ahigh,nhigh,Ehigh),arrlow=Arrheniusvec(Alow,nlow,Elow),
        a=a,T3=T3,T1=T1,T2=T2,efficiencies=efficiencies)
end
export Troevec
    
@inline function (tr::Troevec)(;T::Q=nothing,P::R=0.0,C::S=nothing,phi=0.0,dGrxns=0.0) where {Q<:Real,R<:Real,S<:Real}
    k0 = tr.arrlow(T=T)
    kinf = tr.arrhigh(T=T)
    @fastmath Pr = k0.*C./kinf
    @fastmath log10Pr = log10.(Pr)
    @fastmath log10Fcent = log10.((1.0.-tr.a).*exp.(-T./tr.T3).+tr.a.*exp.(-T./tr.T1).+exp.(-tr.T2./T))
    d = 0.14
    @fastmath n = 0.75.-1.27.*log10Fcent
    @fastmath c = -0.4.-0.67.*log10Fcent
    F = 10.0.^((.!isinf.(tr.T1)) .* @fastmath (log10Fcent./(1.0.+((log10Pr.+c)./(n.-d.*(log10Pr))).^2)))
    return @fastmath ((k0.*C)./(1.0.+Pr)).*F
end

@with_kw struct PdepArrheniusvec{T<:Real,Q<:AbstractRateUncertainty,Z<:Arrheniusvec} <: AbstractRate
    Ps::Array{T,1}
    arrvecs::Array{Z,1}
    unc::Q = EmptyRateUncertainty()
end

function PdepArrheniusvec(pdeparrs::T) where {T<:AbstractArray}
    Ps = pdeparrs[1].Ps
    arrs = [Array{Arrhenius,1}() for i = 1:length(Ps)]
    for ind in 1:length(pdeparrs)
        for i =1:length(Ps)
            push!(arrs[i],pdeparrs[ind].arrs[i])
        end  
    end
    return PdepArrheniusvec(;Ps=Ps,arrvecs=Arrheniusvec.(arrs))
end
export PdepArrheniusvec

@inline function (parr::PdepArrheniusvec)(;T::Q=nothing,P::V=nothing,C::S=0.0,phi=0.0,dGrxns=0.0) where {Q<:Real,V<:Real,S<:Real}
    inds = getBoundingIndsSorted(P,parr.Ps)::Tuple{Int64,Int64}
    if inds[2] == -1
        return @inbounds parr.arrvecs[inds[1]](T=T)
    else
        @inbounds highk = parr.arrvecs[inds[2]](T=T)
        @inbounds lowk = parr.arrvecs[inds[1]](T=T)
        @inbounds Plow = parr.Ps[inds[1]]
        @inbounds Phigh = parr.Ps[inds[2]]
        return @inbounds @fastmath lowk.*10.0.^(log10.(P./Plow)/log10.(Phigh./Plow)*log10.(highk./lowk))
    end
end
export PdepArrheniusvec
