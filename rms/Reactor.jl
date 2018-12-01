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

@inline function addreactionratecontribution!(dydt::Array{Q,1},rxn::ElementaryReaction,cs::Array{W,1},kfs::Array{Z,1},krevs::Array{Z,1}) where {Q<:Real,Z<:Real,T<:Integer,W<:Real}
    R = getrate(rxn,cs,kfs,krevs)
    for ind in rxn.reactantinds
        @fastmath @inbounds dydt[ind] -= R
    end
    for ind in rxn.productinds
        @fastmath @inbounds dydt[ind] += R
    end
end
export addreactionratecontribution!

@inline function dydtreactor!(y::Array{U,1},t::Z,domain::Q;sensitivity::Bool=true) where {Z<:Real,U<:Real,J<:Integer,Q<:AbstractDomain}
    dydt = zeros(U,length(y))
    if sensitivity #if sensitivity isn't explicitly set to false set it to domain.sensitivity
        sensitivity = domain.sensitivity
    end
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs = calcthermo(domain,y,t)
    for rxn in domain.phase.reactions
        addreactionratecontribution!(dydt,rxn,cs,kfs,krevs)
    end
    dydt *= V
    calcdomainderivatives!(domain,dydt;T=T,Us=Us,V=V,C=C,ns=ns,N=N)
    if sensitivity
        Nspcs = length(cs)
        Nrxns = length(domain.phase.reactions)
        jacobian!(domain.jacobian,y,nothing,t,domain)
        dgdk = ratederivative(d=domain,cs=cs,V=V,T=T,kfs=kfs,krevs=krevs,sparse=false)
        for j  = 1:Nspcs+Nrxns #kfs and Gfs
            for i = 1:Nspcs #species
                for z in 1:Nspcs
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
        return domain.jacobian
    else
        f(y::Array{T,1}) where {T<:Real} = dydtBatchReactor!(y,domain.t[1],domain;sensitivity=false)
        ForwardDiff.jacobian!(domain.jacobian,f,y[1:length(domain.phase.species)])
        domain.jacuptodate[1] = true
        return domain.jacobian
    end
end
export jacobian!

function ratederivative(;d::W,cs::Q,V::Y,T::Y2,kfs::Z,krevs::X,sparse::Bool=false) where {W<:Union{ConstantTPDomain,ConstantTVDomain},Q<:AbstractArray,Y2<:Real,Y<:Real,Z<:AbstractArray,X<:AbstractArray}
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
export ratederivative
# function jacobianbatchreactor(y::Array{T,1},t::T,domain::Q,kfs::Array{T,1},krevs::Array{T,1},N::J,sparse::Bool=false) where {T<:Real,J<:Integer,Q<:AbstractConstantKDomain}
#     if sparse
#         jac = spzeros((N,N))
#     else
#         jac = zeros((N,N))
#     end
#     jacobianbatchreactor!(y,t,domain,kfs,krevs,jac;zero=false)
# end
# export jacobianbatchreactor
#
# @inline function spreadpartials!(jac::S,deriv::T,inds::V,ind::Q,N::Q) where {S<:AbstractArray, T<:Real, V<:AbstractArray, Q<:Integer}
#     if N == 1
#         jac[inds[1],ind] += deriv
#     elseif N == 2
#         jac[inds[1],ind] += deriv
#         jac[inds[2],ind] += deriv
#     elseif N == 3
#         jac[inds[1],ind] += deriv
#         jac[inds[2],ind] += deriv
#         jac[inds[3],ind] += deriv
#     end
# end
#
# function jacobianbatchreactor!(y::Array{T,1},t::T,domain::ConstantTPDomain,kfs::Array{T,1},krevs::Array{T,1},jac::P;zero::Bool=true) where {P<:AbstractArray,T<:Real,J<:Integer}
#     if zero
#         jac .= 0
#     end #need some setup here
#     for (i,rxn) in enumerate(domain.phase.reactions)
#         Nreact = length(rxn.reactantinds)
#         Nprod = length(rxn.productinds)
#         kf = kfs[i]
#         krev = krevs[i]
#         if Nreact == 1
#             ind1 = rxn.reactantinds[1]
#             jac[ind1,ind1] -= kf
#             if Nprod == 1
#                 jac[rxn.productinds[1],ind1] += kf
#             elseif Nprod == 2
#                 jac[rxn.productinds[1],ind1] += kf
#                 jac[rxn.productinds[2],ind1] += kf
#             elseif Nprod == 3
#                 jac[rxn.productinds[1],ind1] += kf
#                 jac[rxn.productinds[2],ind1] += kf
#                 jac[rxn.productinds[3],ind1] += kf
#             end
#         elseif Nreact == 2
#             ind1,ind2 = rxn.reactantinds
#             corr = -kf*state.cs[ind1]*state.cs[ind2]/state.C #correction for the partial of the volume term
#             deriv = kf*state.cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*state.cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind1,Nprod)
#             for i in 1:length(domain.phase.species)
#                 jac[ind1,i] -= corr
#                 jac[ind2,i] -= corr
#                 spreadpartials!(jac,corr,rxn.productinds,i,Nprod)
#             end
#         elseif Nreact == 3
#             ind1,ind2,ind3 = rxn.reactantinds
#             corr = -2.0*kf*state.cs[ind1]*state.cs[ind2]*state.cs[ind3]/state.C
#             deriv = kf*state.cs[ind1]*state.cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind3,Nprod)
#             deriv = kf*state.cs[ind1]*state.cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*state.cs[ind3]*state.cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             jac[ind3,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind1,Nprod)
#             for i in 1:length(domain.phase.species)
#                 jac[ind1,i] -= corr
#                 jac[ind2,i] -= corr
#                 jac[ind3,i] -= corr
#                 spreadpartials!(jac,corr,rxn.productinds,i,Nprod)
#             end
#         end
#         #reverse direction
#         if Nprod == 1
#             ind1 = rxn.productinds[1]
#             jac[ind1,ind1] -= krev
#             if Nprod == 1
#                 jac[rxn.reactantinds[1],ind1] += krev
#             elseif Nprod == 2
#                 jac[rxn.reactantinds[1],ind1] += krev
#                 jac[rxn.reactantinds[2],ind1] += krev
#             elseif Nprod == 3
#                 jac[rxn.reactantinds[1],ind1] += krev
#                 jac[rxn.reactantinds[2],ind1] += krev
#                 jac[rxn.reactantinds[3],ind1] += krev
#             end
#         elseif Nreact == 2
#             ind1,ind2 = rxn.productinds
#             corr = -krev*state.cs[ind1]*state.cs[ind2]/state.C #correction for the partial of the volume term
#             deriv = krev*state.cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*state.cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind1,Nreact)
#             for i in 1:length(domain.phase.species)
#                 jac[ind1,i] -= corr
#                 jac[ind2,i] -= corr
#                 spreadpartials!(jac,corr,rxn.reactantinds,i,Nreact)
#             end
#         elseif Nreact == 3
#             ind1,ind2,ind3 = rxn.productinds
#             corr = -2.0*krev*state.cs[ind1]*state.cs[ind2]*state.cs[ind3]/state.C
#             deriv = krev*state.cs[ind1]*state.cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind3,Nreact)
#             deriv = kf*state.cs[ind1]*state.cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*state.cs[ind3]*state.cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             jac[ind3,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind1,Nreact)
#             for i in 1:length(domain.phase.species)
#                 jac[ind1,i] -= corr
#                 jac[ind2,i] -= corr
#                 jac[ind3,i] -= corr
#                 spreadpartials!(jac,corr,rxn.reactantinds,i,Nreact)
#             end
#         end
#     end
#     for ind in domain.constantspeciesinds
#         jac[ind,:] .= 0
#     end
#     domain.state.jacobian = jac
#     return jac
# end
# function jacobianbatchreactor!(y::Array{T,1},t::T,domain::ConstantTVDomain,kfs::Array{T,1},krevs::Array{T,1},jac::P;zero::Bool=true) where {P<:AbstractArray,T<:Real,J<:Integer}
#     if zero
#         jac .= 0
#     end #need some setup here
#     temp_10_13 = 0.0
#     for (i,rxn) in enumerate(domain.phase.reactions)
#         Nreact = length(rxn.reactantinds)
#         Nprod = length(rxn.productinds)
#         kf = kfs[i]
#         krev = krevs[i]
#         if Nreact == 1
#             ind1 = rxn.reactantinds[1]
#             jac[ind1,ind1] -= kf
#             if Nprod == 1
#                 jac[rxn.productinds[1],ind1] += kf
#             elseif Nprod == 2
#                 jac[rxn.productinds[1],ind1] += kf
#                 jac[rxn.productinds[2],ind1] += kf
#             elseif Nprod == 3
#                 jac[rxn.productinds[1],ind1] += kf
#                 jac[rxn.productinds[2],ind1] += kf
#                 jac[rxn.productinds[3],ind1] += kf
#             end
#         elseif Nreact == 2
#             ind1,ind2 = rxn.reactantinds
#             deriv = kf*state.cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*state.cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind1,Nprod)
#         elseif Nreact == 3
#             ind1,ind2,ind3 = rxn.reactantinds
#             deriv = kf*state.cs[ind1]*state.cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind3,Nprod)
#             deriv = kf*state.cs[ind1]*state.cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*state.cs[ind3]*state.cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             jac[ind3,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind1,Nprod)
#         end
#         if jac[10,13] != temp_10_13
#             temp_10_13 = jac[10,13]
#         end
#         #reverse direction
#         if Nprod == 1
#             ind1 = rxn.productinds[1]
#             jac[ind1,ind1] -= krev
#             if Nprod == 1
#                 jac[rxn.reactantinds[1],ind1] += krev
#             elseif Nprod == 2
#                 jac[rxn.reactantinds[1],ind1] += krev
#                 jac[rxn.reactantinds[2],ind1] += krev
#             elseif Nprod == 3
#                 jac[rxn.reactantinds[1],ind1] += krev
#                 jac[rxn.reactantinds[2],ind1] += krev
#                 jac[rxn.reactantinds[3],ind1] += krev
#             end
#         elseif Nreact == 2
#             ind1,ind2 = rxn.productinds
#             deriv = krev*state.cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*state.cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind1,Nreact)
#         elseif Nreact == 3
#             ind1,ind2,ind3 = rxn.productinds
#             deriv = krev*state.cs[ind1]*state.cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind3,Nreact)
#             deriv = kf*state.cs[ind1]*state.cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*state.cs[ind3]*state.cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             jac[ind3,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind1,Nreact)
#         end
#         if jac[10,13] != temp_10_13
#             temp_10_13 = jac[10,13]
#         end
#     end
#     for ind in domain.constantspeciesinds
#         jac[ind,:] .= 0
#     end
#     domain.state.jacobian = jac
#     return jac
# end
# export jacobianbatchreactor!
