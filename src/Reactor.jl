using Parameters
using DifferentialEquations
using ForwardDiff

abstract type AbstractReactor end
export AbstractReactor

struct Reactor{D<:AbstractDomain} <: AbstractReactor
    domain::D
    ode::ODEProblem
end

function Reactor(domain::T,y0::Array{W,1},tspan::Tuple) where {T<:AbstractDomain,W<:Real}
    dydt(y::Array{T,1},p::Nothing,t::Q) where {T<:Real,Q<:Real} = dydtreactor!(y,t,domain)
    ode = ODEProblem(dydt,y0,tspan)
    return Reactor(domain,ode)
end
export Reactor

@inline function getrate(rxn::T,cs::Array{W,1},kfs::Array{Q,1},krevs::Array{Q,1}) where {T<:AbstractReaction,Q,W<:Real}
    Nreact = length(rxn.reactantinds)
    Nprod = length(rxn.productinds)
    R = 0.0
    if Nreact == 1
        @fastmath @inbounds R += kfs[rxn.index]*cs[rxn.reactantinds[1]]
    elseif Nreact == 2
        @fastmath @inbounds R += kfs[rxn.index]*cs[rxn.reactantinds[1]]*cs[rxn.reactantinds[2]]
    elseif Nreact == 3
        @fastmath @inbounds R += kfs[rxn.index]*cs[rxn.reactantinds[1]]*cs[rxn.reactantinds[2]]*cs[rxn.reactantinds[3]]
    end

    if Nprod == 1
        @fastmath @inbounds R -= krevs[rxn.index]*cs[rxn.productinds[1]]
    elseif Nprod == 2
        @fastmath @inbounds R -= krevs[rxn.index]*cs[rxn.productinds[1]]*cs[rxn.productinds[2]]
    elseif Nprod == 3
        @fastmath @inbounds R -= krevs[rxn.index]*cs[rxn.productinds[1]]*cs[rxn.productinds[2]]*cs[rxn.productinds[3]]
    end

    return R
end
export getrate

@inline function addreactionratecontributions!(dydt::Array{Q,1},rarray::Array{UInt16,2},cs::Array{W,1},kfs::Array{Z,1},krevs::Array{Z,1}) where {Q<:Real,Z<:Real,T<:Integer,W<:Real}
    @inbounds @simd for i = 1:size(rarray)[2]
        if @inbounds rarray[2,i] == 0
            @inbounds @fastmath fR = kfs[i]*cs[rarray[1,i]]
        elseif @inbounds rarray[3,i] == 0
            @inbounds @fastmath fR = kfs[i]*cs[rarray[1,i]]*cs[rarray[2,i]]
        else
            @inbounds @fastmath fR = kfs[i]*cs[rarray[1,i]]*cs[rarray[2,i]]*cs[rarray[3,i]]
        end
        if @inbounds rarray[5,i] == 0
            @inbounds @fastmath rR = krevs[i]*cs[rarray[4,i]]
        elseif @inbounds rarray[6,i] == 0
            @inbounds @fastmath rR = krevs[i]*cs[rarray[4,i]]*cs[rarray[5,i]]
        else
            @inbounds @fastmath rR = krevs[i]*cs[rarray[4,i]]*cs[rarray[5,i]]*cs[rarray[6,i]]
        end
        @fastmath R = fR - rR
        @inbounds @fastmath dydt[rarray[1,i]] -= R
        if @inbounds rarray[2,i] != 0
            @inbounds @fastmath dydt[rarray[2,i]] -= R
            if @inbounds rarray[3,i] != 0
                @inbounds @fastmath dydt[rarray[3,i]] -= R
            end
        end
        @inbounds @fastmath dydt[rarray[4,i]] += R
        if @inbounds rarray[5,i] != 0
            @inbounds @fastmath dydt[rarray[5,i]] += R
            if @inbounds rarray[6,i] != 0
                @inbounds @fastmath dydt[rarray[6,i]] += R
            end
        end
    end
end
export addreactionratecontributions!

@inline function dydtreactor!(y::Array{U,1},t::Z,domain::Q;sensitivity::Bool=true) where {Z<:Real,U<:Real,J<:Integer,Q<:AbstractDomain}
    dydt = zeros(U,length(y))
    if sensitivity #if sensitivity isn't explicitly set to false set it to domain.sensitivity
        sensitivity = domain.sensitivity
    end
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave = calcthermo(domain,y,t)
    addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
    dydt *= V
    if sensitivity && isa(domain,ConstantVDomain)
        wV = copy(dydt[domain.indexes[1]:domain.indexes[2]])
    else
        wV = Array{U,1}()
    end
    calcdomainderivatives!(domain,dydt;t=t,T=T,P=P,Us=Us,V=V,C=C,ns=ns,N=N,Cvave=Cvave)
    if sensitivity
        Nspcs = length(cs)
        Nrxns = length(domain.phase.reactions)
        jacobian!(domain.jacobian,y,nothing,t,domain)
        dgdk = ratederivative(domain;cs=cs,V=V,T=T,kfs=kfs,Us=Us,Cvave=Cvave,N=N,krevs=krevs,wV=wV,sparse=false)
        for j  = 1:Nspcs+Nrxns #kfs and Gfs
            for i = 1:(domain.indexes[end]-domain.indexes[1]+1) #species
                for z in 1:(domain.indexes[end]-domain.indexes[1]+1)
                    dydt[j*Nspcs+i] += domain.jacobian[i,z]*y[j*Nspcs+z]
                end
                dydt[j*Nspcs+i] += dgdk[i,j]
            end
        end
    end
    return dydt
end
export dydtreactor!

function jacobian!(J::Q,y::U,p::W,t::Z,domain::V) where {Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
    if domain.t[1] == t && domain.jacuptodate == true
        return
    else
        f(y::Array{T,1}) where {T<:Real} = dydtreactor!(y,domain.t[1],domain;sensitivity=false)
        ForwardDiff.jacobian!(domain.jacobian,f,y[1:domain.indexes[end]-domain.indexes[1]+1])
        domain.jacuptodate[1] = true
        return
    end
end
export jacobian!

function ratederivative(d::W; cs::Q,V::Y,T::Y2,Us::Z3,Cvave::Z4,N::Z5,kfs::Z,krevs::X,wV::Q2,sparse::Bool=false) where {W<:Union{ConstantTPDomain,ConstantTVDomain},Z4<:Real,Z5<:Real,Z3<:AbstractArray,Q2<:AbstractArray,Q<:AbstractArray,Y2<:Real,Y<:Real,Z<:AbstractArray,X<:AbstractArray}
    Nspcs = length(cs)
    rxns = d.phase.reactions
    Nrxns = length(rxns)

    if sparse
        ratederiv = spzeros(Nspcs,Nspcs+Nrxns)
    else
        ratederiv = zeros(Nspcs,Nspcs+Nrxns)
    end
    RTinv = 1.0/(R*T)

    for (j,rxn) in enumerate(rxns)
        Nreact = length(rxn.reactantinds)
        Nprod = length(rxn.productinds)

        if Nreact == 1
            rind1 = rxn.reactantinds[1]
            fderiv = cs[rind1]
        elseif Nreact == 2
            rind1,rind2 = rxn.reactantinds
            fderiv = cs[rind1]*cs[rind2]
        else
            rind1,rind2,rind3 = rxn.reactantinds
            fderiv = cs[rind1]*cs[rind2]*cs[rind3]
        end

        if Nprod == 1
            pind1 = rxn.productinds[1]
            rderiv = krevs[j]/kfs[j]*cs[pind1]
        elseif Nprod == 2
            pind1,pind2 = rxn.productinds
            rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]
        else
            pind1,pind2,pind3 = rxn.productinds
            rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]*cs[pind3]
        end

        flux = fderiv-rderiv
        gderiv = rderiv*kfs[j]*RTinv

        deriv = zeros(Nspcs)

        deriv[rind1] += gderiv
        if Nreact > 1
            deriv[rind2] += gderiv
            if Nreact > 2
                deriv[rind3] == gderiv
            end
        end

        deriv[pind1] -= gderiv
        if Nprod > 1
            deriv[pind2] -= gderiv
            if Nprod > 2
                deriv[pind3] -= gderiv
            end
        end

        ratederiv[rind1,j] -= flux
        ratederiv[rind1,Nrxns+1:Nrxns+Nspcs] .-= deriv
        if Nreact > 1
            ratederiv[rind2,j] -= flux
            ratederiv[rind2,Nrxns+1:Nrxns+Nspcs] .-= deriv
            if Nreact > 2
                ratederiv[rind3,j] -= flux
                ratederiv[rind3,Nrxns+1:Nrxns+Nspcs] .-= deriv
            end
        end

        ratederiv[pind1,j] += flux
        ratederiv[pind1,Nrxns+1:Nrxns+Nspcs] .+= deriv
        if Nprod > 1
            ratederiv[pind2,j] += flux
            ratederiv[pind2,Nrxns+1:Nrxns+Nspcs] .+= deriv
            if Nprod > 2
                ratederiv[pind3,j] += flux
                ratederiv[pind3,Nrxns+1:Nrxns+Nspcs] .+= deriv
            end
        end
    end
    return V*ratederiv
end

function ratederivative(d::ConstantVDomain; cs::Q,V::Y,T::Y2,Us::Z3,Cvave::Y3,N::Y2,kfs::Z,krevs::X,wV::Q2,sparse::Bool=false) where {Z3<:AbstractArray,Q<:AbstractArray,Q2<:AbstractArray,Y3<:Real,Y2<:Real,Y<:Real,Z<:AbstractArray,X<:AbstractArray}
    Nspcs = length(cs)
    rxns = d.phase.reactions
    Nrxns = length(rxns)

    if sparse
        ratederiv = spzeros(Nspcs+1,Nspcs+Nrxns)
    else
        ratederiv = zeros(Nspcs+1,Nspcs+Nrxns)
    end
    RTinv = 1.0/(R*T)

    for (j,rxn) in enumerate(rxns)
        Nreact = length(rxn.reactantinds)
        Nprod = length(rxn.productinds)

        if Nreact == 1
            rind1 = rxn.reactantinds[1]
            fderiv = cs[rind1]
        elseif Nreact == 2
            rind1,rind2 = rxn.reactantinds
            fderiv = cs[rind1]*cs[rind2]
        else
            rind1,rind2,rind3 = rxn.reactantinds
            fderiv = cs[rind1]*cs[rind2]*cs[rind3]
        end

        if Nprod == 1
            pind1 = rxn.productinds[1]
            rderiv = krevs[j]/kfs[j]*cs[pind1]
        elseif Nprod == 2
            pind1,pind2 = rxn.productinds
            rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]
        else
            pind1,pind2,pind3 = rxn.productinds
            rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]*cs[pind3]
        end

        flux = fderiv-rderiv
        gderiv = rderiv*kfs[j]*RTinv

        deriv = zeros(Nspcs)

        deriv[rind1] += gderiv
        if Nreact > 1
            deriv[rind2] += gderiv
            if Nreact > 2
                deriv[rind3] == gderiv
            end
        end

        deriv[pind1] -= gderiv
        if Nprod > 1
            deriv[pind2] -= gderiv
            if Nprod > 2
                deriv[pind3] -= gderiv
            end
        end

        ratederiv[rind1,j] -= flux
        ratederiv[rind1,Nrxns+1:Nrxns+Nspcs] .-= deriv
        if Nreact > 1
            ratederiv[rind2,j] -= flux
            ratederiv[rind2,Nrxns+1:Nrxns+Nspcs] .-= deriv
            if Nreact > 2
                ratederiv[rind3,j] -= flux
                ratederiv[rind3,Nrxns+1:Nrxns+Nspcs] .-= deriv
            end
        end

        ratederiv[pind1,j] += flux
        ratederiv[pind1,Nrxns+1:Nrxns+Nspcs] .+= deriv
        if Nprod > 1
            ratederiv[pind2,j] += flux
            ratederiv[pind2,Nrxns+1:Nrxns+Nspcs] .+= deriv
            if Nprod > 2
                ratederiv[pind3,j] += flux
                ratederiv[pind3,Nrxns+1:Nrxns+Nspcs] .+= deriv
            end
        end
    end
    ratederiv *= V
    #Temperature stuff
    @views ratederiv[end,:] += (Us'*ratederiv[1:end-1,:])[1,:]
    ratederiv[end,1:Nspcs] += wV
    ratederiv[end,:] /= (N*Cvave)
    return ratederiv
end
export ratederivative
