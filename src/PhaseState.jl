using Parameters
using SpecialFunctions
using LinearAlgebra
using Tracker
using ReverseDiff
using RecursiveArrayTools
using Logging

@inline function calcgibbs(ph::U,T::W) where {U<:IdealPhase,W<:Real}
    return getGibbs.(getfield.(ph.species,:thermo),T)
end
export calcgibbs

@inline function calcenthalpyinternalgibbs(ph::U,T::W,P::Z,V::Q) where {U<:IdealPhase,W,Z,Q<:Real}
    Hs = getEnthalpy.(getfield.(ph.species,:thermo),T)
    Us = Hs .- P*V
    Gs = Hs .- T.*getEntropy.(getfield.(ph.species,:thermo),T)
    return Hs,Us,Gs
end


@inline function calcenthalpyinternalgibbs(ph::Union{IdealGas,IdealSurface},T::W,P::Z,V::Q) where {W,Z,Q<:Real}
    Hs = getEnthalpy.(getfield.(ph.species,:thermo),T)
    Us = Hs .- R*T
    Gs = Hs .- T.*getEntropy.(getfield.(ph.species,:thermo),T)
    return Hs,Us,Gs
end

export calcenthalpyinternalgibbs

function makespcsvector(phase,spcdict)
    y0 = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key,val) in spcdict
        if key == "T"
            continue
        elseif key == "P"
            continue
        elseif key == "V"
            continue
        elseif key == "Hin"
            continue
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end
    return y0
end

export makespcsvector

@inline function getkf(rxn::ElementaryReaction,ph,T,P,C,ns,V,phi,dGrxn,d)
    if isdefined(rxn.kinetics,:efficiencies) && length(rxn.kinetics.efficiencies) > 0
        @views @inbounds @fastmath C += sum([ns[i]*val for (i,val) in rxn.kinetics.efficiencies])/V
    end
    return rxn.kinetics(T=T,P=P,C=C,phi=phi,dGrxn=dGrxn,d=d)
end
export getkf

@inline function getkfs(ph::U,T::W1,P::W2,C::W3,ns::Q,V::W4,phi,dGrxns,d) where {U,W1,W2,W3,W4<:Real,Q<:AbstractArray}
    kfs = similar(ns,length(ph.reactions))
    i = 1
    oldind = 1
    ind = 0
    while i <= length(ph.veckineticsinds) #vectorized kinetics
        @inbounds ind = ph.veckineticsinds[i]
        @inbounds kfs[oldind:ind] = ph.veckinetics[i](;T=T,P=P,C=C)
        oldind = ind+1
        i += 1
    end
    @simd for i in ind+1:length(ph.reactions)
        @inbounds kfs[i] = getkf(ph.reactions[i],ph,T,P,C,ns,V,phi,dGrxns[i],d)
    end
    return kfs
end

export getkfs

"""
Calculates the diffusion limited rate coefficient
for 1 spc returns Inf
for 2 spc calculates using the Smolchowski equation
for >2 spc calculates using the Generalized Smolchowski equation
Equations from Flegg 2016
"""
@inline function getDiffusiveRate(spcs::Q,spcsinds::Q2,diffs::Array{W,1}) where {Q<:AbstractArray,W<:Real,Q2<:AbstractArray}
    if length(spcs) == 1
        return Inf
    elseif length(spcs) == 2
        @fastmath @inbounds kf = 4.0*Base.pi*(diffs[spcsinds[1]]+diffs[spcsinds[2]])*(spcs[1].radius+spcs[2].radius)*Na
    else
        @views @inbounds diffusivity = diffs[spcsinds]
        N = length(spcs)
        @fastmath a = (3.0*length(spcs)-5.0)/2.0
        @fastmath Dinv = 1.0./diffusivity
        @fastmath Dbar = 1.0./reverse(cumsum(Dinv))
        @fastmath Dhat = diffusivity .+ Dbar
        @fastmath @inbounds deltaN = sum(Dinv)/sum([sum([1.0/(diffusivity[i]*diffusivity[m]) for m in 1:N-1 if i>m]) for i in 2:N])
        @views @fastmath @inbounds kf = prod(Dhat[2:end].^1.5)*4*Base.pi^(a+1)/gamma(a)*(sum(getfield.(spcs,:radius))/sqrt(deltaN))^(2*a)*Na^(N-1)
    end
    return kf
end
export getDiffusiveRate

@inline function getKc(rxn::ElementaryReaction,ph::U,T::Z,dGrxn::Q,phi::V=0.0) where {U<:AbstractPhase,V,Q,Z<:Real}
    Nreact = length(rxn.reactantinds)
    Nprod = length(rxn.productinds)
    return @inbounds @fastmath exp((dGrxn+rxn.electronchange*(phi*F))/(-R*T))*(getC0(ph,T))^(Nprod-Nreact)
end

@inline function getKcs(ph::U,T::Z,dGrxns::Q) where {U<:AbstractPhase,Q,Z<:Real}
    return @fastmath @inbounds exp.(dGrxns./(-R*T) .+ ph.Nrp.*log(getC0(ph,T)));
end

@inline function getKcs(ph::U,T::Z,dGrxns::Q,phi::V) where {U<:AbstractPhase,Q,Z<:Real,V<:Real}
    return @fastmath @inbounds exp.(dGrxns./(-R*T) .+ ph.Nrp.*log(getC0(ph,T)));
end

@inline function getKcs(ph,T,Gs1,Gs2,dGrxns)
    return @fastmath @inbounds exp.(dGrxns./(-R*T) .+ ph.Nrp1.*log(getC0(ph.domain1.phase,T)) .+ ph.Nrp2.*log(getC0(ph.domain2.phase,T)));
end

@inline function getKcs(ph1,ph2,T,Nrp1,Nrp2,dGrxns)
    return @fastmath @inbounds exp.(dGrxns./(-R*T) .+ Nrp1.*log(getC0(ph1,T)) .+ Nrp2.*log(getC0(ph2,T)));
end

export getKcs

"""
Calculates the forward and reverse rate coefficients for a given reaction, phase and state
Maintains diffusion limitations if the phase has diffusionlimited=true
"""
@inline function getkfkrev(rxn::ElementaryReaction,ph::U,T::W1,P::W2,C::W3,N::W4,ns::Q1,dGrxn::Q2,diffs::Q3,V::W5,phi::W8,d;kf::W6=-1.0,f::W7=-1.0) where {U<:AbstractPhase,W8,W6,W7,W5,W4,W1,W2,W3<:Real,Q1,Q2,Q3<:AbstractArray}
    if signbit(kf)
        if signbit(f)
            kf = getkf(rxn,ph,T,P,C,ns,V,phi,dGrxn,d)
        else
            kf = getkf(rxn,ph,T,P,C,ns,V,phi,dGrxn,d)*f
        end
    end
    Kc = getKc(rxn,ph,T,dGrxn,phi)
    @fastmath krev = kf/Kc
    if ph.diffusionlimited
        if length(rxn.reactants) == 1
            if length(rxn.products) > 1
                krevdiff = getDiffusiveRate(rxn.products,rxn.productinds,diffs)
                @fastmath krev = krev*krevdiff/(krev+krevdiff)
                @fastmath kf = Kc*krev
            end
        elseif length(rxn.products) == 1
            kfdiff = getDiffusiveRate(rxn.reactants,rxn.reactantinds,diffs)
            @fastmath kf = kf*kfdiff/(kf+kfdiff)
            @fastmath krev = kf/Kc
        elseif length(rxn.products) == length(rxn.reactants)
            kfdiff = getDiffusiveRate(rxn.reactants,rxn.reactantinds,diffs)
            krevdiff = getDiffusiveRate(rxn.products,rxn.productinds,diffs)
            @fastmath kff = kf*kfdiff/(kf+kfdiff)
            @fastmath krevr = krev*krevdiff/(krev+krevdiff)
            @fastmath kfr = Kc*krevr
            if kff > kfr
                kf = kfr
                krev = krevr
            else
                kf = kff
                @fastmath krev = kf/Kc
            end
        end
    end
    kf *= rxn.forwardable
    krev *= rxn.reversible
    return (kf,krev)
end
export getkfkrev

@inline function getkfkrevs(phase::U,T::W1,P::W2,C::W3,N::W4,ns::Q1,Gs::Q2,diffs::Q3,V::W5,phi::W7,d;kfs::W6=nothing) where {U,W7,W6,W5<:Real,W1<:Real,W2<:Real,W3,W4,Q1<:AbstractArray,Q2,Q3<:AbstractArray}
    dGrxns = -(phase.stoichmatrix*Gs).+phase.electronchange.*(phi*F)
    if !phase.diffusionlimited && kfs === nothing
        kfs = getkfs(phase,T,P,C,ns,V,phi,dGrxns,d)
        if phi == 0.0
            krev = @fastmath kfs./getKcs(phase,T,dGrxns)
        else 
            krev = @fastmath kfs./getKcs(phase,T,dGrxns,phi)
        end
    elseif !phase.diffusionlimited && !(kfs === nothing)
        if phi == 0.0
            krev = @fastmath kfs./getKcs(phase,T,dGrxns)
        else 
            krev = @fastmath kfs./getKcs(phase,T,dGrxns,phi)
        end
    elseif phase.diffusionlimited && !(kfs === nothing)
        len = length(phase.reactions)
        krev = zeros(typeof(N),len)
        @simd for i = 1:len
           @fastmath @inbounds kfs[i],krev[i] = getkfkrev(phase.reactions[i],phase,T,P,C,N,ns,dGrxns[i],diffs,V,phi,d;kf=kfs[i])
        end
    else
        len = length(phase.reactions)
        kfs = zeros(typeof(N),len)
        krev = zeros(typeof(N),len)
        @simd for i = 1:len
           @fastmath @inbounds kfs[i],krev[i] = getkfkrev(phase.reactions[i],phase,T,P,C,N,ns,dGrxns[i],diffs,V,phi,d)
        end
    end
    kfs .*= phase.forwardability
    krev .*= phase.reversibility
    return kfs,krev
end

@inline function getkfkrevs(phase::U,T::W1,P::W2,C::W3,N::W4,ns::Q1,Gs::Q2,diffs::Q3,V::W5,phi::W7,d;kfs::W6=nothing) where {U,W7,W6,W5<:Real,W1<:Real,W2<:Real,W3,W4,Q1<:AbstractArray,Q2<:Union{ReverseDiff.TrackedArray,Tracker.TrackedArray},Q3<:AbstractArray} #autodiff p
    dGrxns = -(phase.stoichmatrix*Gs).+phase.electronchange.*(phi*F)
    if !phase.diffusionlimited && kfs === nothing
        kfs = getkfs(phase,T,P,C,ns,V,phi,dGrxns,d)
        if phi == 0.0
            krev = @fastmath kfs./getKcs(phase,T,dGrxns)
        else 
            krev = @fastmath kfs./getKcs(phase,T,dGrxns,phi)
        end
    elseif !phase.diffusionlimited && !(kfs === nothing)
        if phi == 0.0
            krev = @fastmath kfs./getKcs(phase,T,dGrxns)
        else 
            krev = @fastmath kfs./getKcs(phase,T,dGrxns,phi)
        end
    elseif phase.diffusionlimited && !(kfs === nothing)
        len = length(phase.reactions)
        krev = similar(kfs)
        kfsdiff = similar(kfs)
        @simd for i = 1:len
           @fastmath @inbounds kfsdiff[i],krev[i] = getkfkrev(phase.reactions[i],phase,T,P,C,N,ns,dGrxns[i],diffs,V,phi,d;kf=kfs[i])
        end
        return kfsdiff, krev
    else
        len = length(phase.reactions)
        kfs = zeros(typeof(Gs[1]),len)ss
        krev = zeros(typeof(Gs[1]),len)
        @simd for i = 1:len
           @fastmath @inbounds kfs[i],krev[i] = getkfkrev(phase.reactions[i],phase,T,P,C,N,ns,dGrxns[i],diffs,V,phi,d)
        end
    end
    return kfs,krev
end

@inline function getkfkrevs(phase::U,T::W1,P::W2,C::W3,N::W4,ns::Q1,Gs::Array{Q2,1},diffs::Q3,V::W5,phi::W7,d;kfs::W6=nothing) where {U,W7,W6,W5<:Real,W1<:Real,W2<:Real,W3,W4,Q1<:AbstractArray,Q2<:ForwardDiff.Dual,Q3<:AbstractArray} #autodiff p
    dGrxns = -(phase.stoichmatrix*Gs).+phase.electronchange.*(phi*F)
    if !phase.diffusionlimited && kfs === nothing
        kfs = getkfs(phase,T,P,C,ns,V,phi,dGrxns,d)
        if phi == 0.0
            krev = @fastmath kfs./getKcs(phase,T,dGrxns)
        else 
            krev = @fastmath kfs./getKcs(phase,T,dGrxns,phi)
        end
    elseif !phase.diffusionlimited && !(kfs === nothing)
        if phi == 0.0
            krev = @fastmath kfs./getKcs(phase,T,dGrxns)
        else 
            krev = @fastmath kfs./getKcs(phase,T,dGrxns,phi)
        end
    elseif phase.diffusionlimited && !(kfs === nothing)
        len = length(phase.reactions)
        krev = similar(kfs)
        @simd for i = 1:len
           @fastmath @inbounds kfs[i],krev[i] = getkfkrev(phase.reactions[i],phase,T,P,C,N,ns,dGrxns[i],diffs,V,phi,d;kf=kfs[i])
        end
    else
        len = length(phase.reactions)
        kfs = zeros(typeof(Gs[1]),len)
        krev = zeros(typeof(Gs[1]),len)
        @simd for i = 1:len
           @fastmath @inbounds kfs[i],krev[i] = getkfkrev(phase.reactions[i],phase,T,P,C,N,ns,dGrxns[i],diffs,V,phi,d)
        end
    end
    kfs .*= phase.forwardability
    krev .*= phase.reversibility
    return kfs,krev
end

export getkfkrevs
