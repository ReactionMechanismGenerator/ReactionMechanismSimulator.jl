using Parameters
using DiffEqBase
using ForwardDiff
using Sundials
using ModelingToolkit
abstract type AbstractReactor end
export AbstractReactor

struct Reactor{D,Q} <: AbstractReactor
    domain::D
    ode::ODEProblem
    recommendedsolver::Q
    forwardsensitivities::Bool
end

function Reactor(domain::T,y0::Array{W,1},tspan::Tuple,interfaces::Z=[];p::X=DiffEqBase.NullParameters(),forwardsensitivities=false,forwarddiff=false,modelingtoolkit=false) where {T<:AbstractDomain,W<:Real,Z<:AbstractArray,X}
    dydt(dy::X,y::T,p::V,t::Q) where {X,T,Q,V} = dydtreactor!(dy,y,t,domain,interfaces,p=p)
    jacy!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q,V} = jacobiany!(J,y,p,t,domain,interfaces,nothing)
    jacyforwarddiff!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q,V} = jacobianyforwarddiff!(J,y,p,t,domain,interfaces,nothing)
    jacp!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q,V} = jacobianp!(J,y,p,t,domain,interfaces,nothing)
    jacpforwarddiff!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q,V} = jacobianpforwarddiff!(J,y,p,t,domain,interfaces,nothing)
    
    if (forwardsensitivities || !forwarddiff) && domain isa Union{ConstantTPDomain,ConstantVDomain,ConstantPDomain,ParametrizedTPDomain,ParametrizedVDomain,ParametrizedPDomain,ConstantTVDomain,ParametrizedTConstantVDomain,ConstantTAPhiDomain}
        if !forwardsensitivities
            odefcn = ODEFunction(dydt;jac=jacy!,paramjac=jacp!)
        else
            odefcn = ODEFunction(dydt;paramjac=jacp!)
        end
    else
        odefcn = ODEFunction(dydt;jac=jacyforwarddiff!,paramjac=jacpforwarddiff!) 
    end
    if forwardsensitivities
        ode = ODEForwardSensitivityProblem(odefcn,y0,tspan,p)
        recsolver = Sundials.CVODE_BDF(linear_solver=:GMRES)
    else
        ode = ODEProblem(odefcn,y0,tspan,p)
        recsolver  = Sundials.CVODE_BDF()
    end
    if modelingtoolkit
        sys = modelingtoolkitize(ode)
        jac = eval(ModelingToolkit.generate_jacobian(sys)[2])
        if (forwardsensitivities || !forwarddiff) && domain isa Union{ConstantTPDomain,ConstantVDomain,ConstantPDomain,ParametrizedTPDomain,ParametrizedVDomain,ParametrizedPDomain,ConstantTVDomain,ParametrizedTConstantVDomain,ConstantTAPhiDomain}
            odefcn = ODEFunction(dydt;jac=jac,paramjac=jacp!)
        else 
            odefcn = ODEFunction(dydt;jac=jac,paramjac=jacpforwarddiff!)
        end
        if forwardsensitivities
            ode = ODEForwardSensitivityProblem(odefcn,y0,tspan,p)
        else
            ode = ODEProblem(odefcn,y0,tspan,p)
        end
    end
    return Reactor(domain,ode,recsolver,forwardsensitivities)
end
function Reactor(domains::T,y0s::W,tspan::W2,interfaces::Z=[],ps::X=DiffEqBase.NullParameters();forwardsensitivities=false,modelingtoolkit=false) where {T<:Tuple,W,Z<:AbstractArray,X,W2}
    #adjust indexing
    y0 = zeros(sum(length(y) for y in y0s))
    Nvars = 0
    for domain in domains
        Nvars += domain.indexes[end]
    end
    n = Nvars
    k = 1
    for (j,domain) in enumerate(domains)
        Nspcs = length(domain.phase.species)
        Ntherm = length(domain.indexes) - 2
        for i = 1:6, j = 1:size(domain.rxnarray)[2]
            if domain.rxnarray[i,j] != 0
                domain.rxnarray[i,j] += k-1
            end
        end
        domain.indexes[1] = k
        k += Nspcs
        domain.indexes[2] = k-1
        y0[domain.indexes[1]:domain.indexes[2]] = y0s[j][1:Nspcs]
        for m = 3:2+Ntherm
            domain.indexes[m] = n
            y0[n] = y0s[j][Nspcs+m-2]
            n -= 1
        end
    end
    
    p = Array{Float64,1}()
    phases = []
    phaseinds = Array{Tuple,1}()
    for (i,domain) in enumerate(domains)
        if domain.phase in phases
            ind = findfirst(phases.==domain.phase)
            domain.parameterindexes[1] = phaseinds[ind][1]
            domain.parameterindexes[2] = phaseinds[ind][2]
            ds = [domains[ind] for ind in findall(x->x.phase==domain.phase,domains)]
            if any(.!isa.(ds,AbstractConstantKDomain)) && any(isa.(ds,AbstractConstantKDomain)) #handle different p formats for the same phase
                for d in ds
                    if !isa(d,AbstractConstantKDomain)
                        p[domain.parameterindexes[1]:domain.parameterindexes[2]] .= d.p
                        break
                    end
                end
                for d in ds
                    if isa(d,AbstractConstantKDomain)
                        if forwardsensitivities == true #error if running into sensitivity bug
                            error("forward sensitivities is not supported for domain combinations that share rate coefficients, but treat them differently ex: ConstantTPDomain and ConstantVDomain")
                        end
                        d.alternativepformat = true
                    end
                end
            end
        else 
            domain.parameterindexes[1] = length(p)+1
            domain.parameterindexes[2] = length(p)+length(ps[i])
            push!(phaseinds,(length(p)+1,length(p)+length(ps[i])))
            push!(phases,domain.phase)
            p = vcat(p,ps[i])
        end
    end
    
    dydt(dy::X,y::T,p::V,t::Q) where {X,T,Q,V} = dydtreactor!(dy,y,t,domains,interfaces,p=p)
    jacy!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q,V} = jacobianyforwarddiff!(J,y,p,t,domains,interfaces,nothing)
    jacp!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q,V} = jacobianpforwarddiff!(J,y,p,t,domains,interfaces,nothing)

    
    if forwardsensitivities
        odefcn = ODEFunction(dydt;paramjac=jacp!)
        ode = ODEForwardSensitivityProblem(odefcn,y0,tspan,p)
        recsolver = Sundials.CVODE_BDF(linear_solver=:GMRES)
        if modelingtoolkit
            sys = modelingtoolkitize(ode)
            jac = eval(ModelingToolkit.generate_jacobian(sys)[2])
            odefcn = ODEFunction(dydt;jac=jac,paramjac=jacp!)
            ode = ODEForwardSensitivityProblem(odefcn,y0,tspan,p)
        end
    else
        odefcn = ODEFunction(dydt;jac=jacy!,paramjac=jacp!)
        ode = ODEProblem(odefcn,y0,tspan,p)
        recsolver  = Sundials.CVODE_BDF()
        if modelingtoolkit
            sys = modelingtoolkitize(ode)
            jac = eval(ModelingToolkit.generate_jacobian(sys)[2])
            odefcn = ODEFunction(dydt;jac=jac,paramjac=jacp!)
            ode = ODEProblem(odefcn,y0,tspan,p)
        end
    end
    return Reactor(domains,ode,recsolver,forwardsensitivities),y0,p
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

@inline function addreactionratecontributions!(dydt::Q,rarray::Array{W2,2},cs::W,kfs::Z,krevs::Y) where {Q,Z,Y,T,W,W2}
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

@inline function dydtreactor!(dydt::RC,y::U,t::Z,domain::Q,interfaces::B;p::RV=DiffEqBase.NullParameters(),sensitivity::Bool=true) where {RC,RV,B<:AbstractArray,Z,U,Q<:AbstractDomain}    
    dydt .= 0.0
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR,phi = calcthermo(domain,y,t,p)
    addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
    dydt .*= V
    calcdomainderivatives!(domain,dydt,interfaces;t=t,T=T,P=P,Us=Us,Hs=Hs,V=V,C=C,ns=ns,N=N,Cvave=Cvave)
    return dydt
end
@inline function dydtreactor!(dydt::RC,y::U,t::Z,domains::Q,interfaces::B;p::RV=DiffEqBase.NullParameters(),sensitivity::Bool=true) where {RC,RV,B<:AbstractArray,Z,U,Q<:Tuple}    
    cstot = zeros(typeof(y).parameters[1],length(y))
    dydt .= 0.0
    domain = domains[1]
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,phi = calcthermo(domain,y,t,p)
    vns = Array{Any,1}(undef,length(domains))
    vns[1] = ns
    vcs = Array{Any,1}(undef,length(domains))
    vcs[1] = cs
    cstot[domain.indexes[1]:domain.indexes[2]] = cs
    vT = Array{Any,1}(undef,length(domains))
    vT[1] = T
    vP = Array{Any,1}(undef,length(domains))
    vP[1] = P
    vV = Array{Any,1}(undef,length(domains))
    vV[1] = V
    vC = Array{Any,1}(undef,length(domains))
    vC[1] = C
    vN = Array{Any,1}(undef,length(domains))
    vN[1] = N
    vmu = Array{Any,1}(undef,length(domains))
    vmu[1] = mu
    vkfs = Array{Any,1}(undef,length(domains))
    vkfs[1] = kfs
    vkrevs = Array{Any,1}(undef,length(domains))
    vkrevs[1] = krevs
    vHs = Array{Any,1}(undef,length(domains))
    vHs[1] = Hs
    vUs = Array{Any,1}(undef,length(domains))
    vUs[1] = Us
    vGs = Array{Any,1}(undef,length(domains))
    vGs[1] = Gs
    vdiffs = Array{Any,1}(undef,length(domains))
    vdiffs[1] = diffs
    vCvave = Array{Any,1}(undef,length(domains))
    vCvave[1] = Cvave
    vphi = Array{Any,1}(undef,length(domains))
    vphi[1] = phi
    addreactionratecontributions!(dydt,domain.rxnarray,cstot,kfs,krevs)
    @views dydt[domain.indexes[1]:domain.indexes[2]] .*= V
    for (i,domain) in enumerate(@views domains[2:end])
        k = i + 1
        vns[k],vcs[k],vT[k],vP[k],vV[k],vC[k],vN[k],vmu[k],vkfs[k],vkrevs[k],vHs[k],vUs[k],vGs[k],vdiffs[k],vCvave[k],vphi[k] = calcthermo(domain,y,t,p)
        cstot[domain.indexes[1]:domain.indexes[2]] .= vcs[k]
        addreactionratecontributions!(dydt,domain.rxnarray,cstot,vkfs[k],vkrevs[k])
        @views dydt[domain.indexes[1]:domain.indexes[2]] .*= vV[k]
    end
    for (i,domain) in enumerate(domains)
        calcdomainderivatives!(domain,dydt,interfaces;t=t,T=vT[i],P=vP[i],Us=vUs[i],Hs=vHs[i],V=vV[i],C=vC[i],ns=vns[i],N=vN[i],Cvave=vCvave[i])
    end
    return dydt
end
export dydtreactor!

function jacobianyforwarddiff!(J::Q,y::U,p::W,t::Z,domain::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
    f(dy::X,y::Array{T,1}) where {T<:Real,X} = dydtreactor!(dy,y,t,domain,interfaces;p=p,sensitivity=false)
    ForwardDiff.jacobian!(J,f,zeros(size(y)),y)
end
export jacobianyforwarddiff!

function jacobianyforwarddiff(y::U,p::W,t::Z,domain::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
    J = zeros(length(y),length(y))
    jacobianyforwarddiff!(J,y,p,t,domain,interfaces,colorvec)
    return J
end
export jacobianyforwarddiff

function jacobianyforwarddiff!(J::Q,y::U,p::W,t::Z,domains::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:Tuple}
    f(dy::X,y::Array{T,1}) where {T<:Real,X} = dydtreactor!(dy,y,t,domains,interfaces;p=p,sensitivity=false)
    ForwardDiff.jacobian!(J,f,zeros(size(y)),y)
end

function jacobianyforwarddiff(y::U,p::W,t::Z,domains::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:Tuple}
    J = zeros(length(y),length(y))
    jacobianyforwarddiff!(J,y,p,t,domains,interfaces,colorvec)
    return J
end
# function jacobiany!(J::Q,y::U,p::W,t::Z,domain::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2<:AbstractArray,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
#     f(y::Array{T,1}) where {T<:Real} = dydtreactor!(y,t,domain,interfaces;p=p,sensitivity=false)
#     forwarddiff_color_jacobian!(J,f,y,colorvec=colorvec)
# end
# function jacobiany!(J::Q,y::U,p::W,t::Z,domain::Q4,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,Q4<:Union{}}
#     ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave = calcthermo(domain,y,t,p)
#     jacobiany!(y,t,domain,kfs,krevs,J;zero=true)
# end
# export jacobiany!

function jacobianpforwarddiff!(J::Q,y::U,p::W,t::Z,domain::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
    function f(dy::X,p::Array{T,1}) where {X,T<:Real} 
        dydtreactor!(dy,y,t,domain,interfaces;p=p,sensitivity=false)
    end
    dy = zeros(length(y))
    ForwardDiff.jacobian!(J,f,dy,p)
end

function jacobianpforwarddiff!(J::Q,y::U,p::W,t::Z,domains::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:Tuple}
    function f(dy::X,p::Array{T,1}) where {X,T<:Real} 
        dydtreactor!(dy,y,t,domains,interfaces;p=p,sensitivity=false)
    end
    dy = zeros(length(y))
    ForwardDiff.jacobian!(J,f,dy,p)
end

export jacobianpforwarddiff!

function jacobianpforwarddiff(y::U,p::W,t::Z,domain::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
    J = zeros(length(y),length(domain.phase.species)+length(domain.phase.reactions))
    jacobianpforwarddiff!(J,y,p,t,domain,interfaces,colorvec)
end

function jacobianpforwarddiff(y::U,p::W,t::Z,domains::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:Tuple}
    J = zeros(length(y),length(domains.phase.species)+length(domains.phase.reactions))
    jacobianpforwarddiff!(J,y,p,t,domains,interfaces,colorvec)
end

export jacobianpforwarddiff

# function jacobianp!(J::Q,y::U,p::W,t::Z,domain::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2<:AbstractArray,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
#     f(p::Array{T,1}) where {T<:Real} = dydtreactor!(y,domain.t[1],domain,interfaces;p=p,sensitivity=false)
#     forwarddiff_color_jacobian!(J,f,p,colorvec=colorvec)
# end
# function jacobianp!(J::Q,y::U,p::W,t::Z,domain::Q4,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,Q4<:Union{ConstantTPDomain,ConstantTVDomain}}
#     ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave = calcthermo(domain,y,t,p)
#     dydt = zeros(length(y))
#     addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
#     @views wV = dydt[domain.indexes[1]:domain.indexes[2]]
#     jacobianp!(domain;cs=cs,V=V,T=T,Us=Us,Cvave=Cvave,N=N,kfs=kfs,krevs=krevs,wV=wV,ratederiv=J)
# end

function jacobiany(y::U,p::W,t::Z,domain::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
    J = zeros(length(y),length(y))
    jacobiany!(J,y,p,t,domain,interfaces,colorvec)
    return J
end
export jacobiany

function jacobianp(y::U,p::W,t::Z,domain::V,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,U<:AbstractArray,W,Z<:Real,V<:AbstractDomain}
    J = zeros(length(y),length(domain.phase.species)+length(domain.phase.reactions))
    jacobianp!(J,y,p,t,domain,interfaces,colorvec)
    return J
end
export jacobianp

@inline function _spreadreactantpartials!(jac::S,deriv::Float64,rxnarray::Array{Int64,2},rxnind::Int64,ind::Int64) where {S<:AbstractArray}
    @inbounds jac[rxnarray[4,rxnind],ind] += deriv
    if @inbounds rxnarray[5,rxnind] !== 0
        @inbounds jac[rxnarray[5,rxnind],ind] += deriv
        if @inbounds rxnarray[6,rxnind] !== 0
            @inbounds jac[rxnarray[6,rxnind],ind] += deriv
        end
    end
end
@inline function _spreadproductpartials!(jac::S,deriv::Float64,rxnarray::Array{Int64,2},rxnind::Int64,ind::Int64) where {S<:AbstractArray}
    @inbounds jac[rxnarray[1,rxnind],ind] += deriv
    if @inbounds rxnarray[2,rxnind] !== 0
        @inbounds jac[rxnarray[2,rxnind],ind] += deriv
        if @inbounds rxnarray[3,rxnind] !== 0
            @inbounds jac[rxnarray[3,rxnind],ind] += deriv
        end
    end
end

@inline function _jacobianynswrtns!(jac::S,rxnarray::Array{Int64,2},rxnind::Int64,cs::Array{Float64,1},kf::Float64,krev::Float64) where {S<:AbstractArray}
    k=kf
    if rxnarray[2,rxnind] == 0
        deriv = k
        @inbounds jac[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= deriv
        @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[1,rxnind])
    elseif rxnarray[3,rxnind] == 0
        if rxnarray[1,rxnind] == rxnarray[2,rxnind]
            @inbounds @fastmath deriv = 2.0*k*cs[rxnarray[1,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= 2.0*deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[1,rxnind])
        else
            @inbounds @fastmath deriv = k*cs[rxnarray[2,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[1,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[1,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[1,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[2,rxnind]] -= deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[2,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[2,rxnind])
        end
    else
        if rxnarray[1,rxnind]==rxnarray[2,rxnind] && rxnarray[1,rxnind]==rxnarray[3,rxnind]
            @inbounds @fastmath deriv = 3.0*k*cs[rxnarray[1,rxnind]]*cs[rxnarray[1,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= 3.0*deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[1,rxnind])
        elseif rxnarray[1,rxnind]==rxnarray[2,rxnind]
            @inbounds @fastmath deriv = 2.0*k*cs[rxnarray[1,rxnind]]*cs[rxnarray[3,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= 2.0*deriv
            @inbounds jac[rxnarray[3,rxnind],rxnarray[1,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[1,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[1,rxnind]]*cs[rxnarray[1,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[3,rxnind]] -= 2.0*deriv
            @inbounds jac[rxnarray[3,rxnind],rxnarray[3,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[3,rxnind])
        elseif rxnarray[2,rxnind]==rxnarray[3,rxnind]
            @inbounds @fastmath deriv = k*cs[rxnarray[2,rxnind]]*cs[rxnarray[2,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[1,rxnind]] -= 2.0*deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[1,rxnind])
            @inbounds @fastmath deriv = 2.0*k*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[2,rxnind]] -= deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[2,rxnind]] -= 2.0*deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[2,rxnind])
        elseif rxnarray[1,rxnind]==rxnarray[3,rxnind]
            @inbounds @fastmath deriv = 2.0*k*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= 2.0*deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[1,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[1,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[1,rxnind]]*cs[rxnarray[1,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[2,rxnind]] -= 2.0*deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[2,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[2,rxnind])
        else
            @inbounds @fastmath deriv = k*cs[rxnarray[2,rxnind]]*cs[rxnarray[3,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[1,rxnind]] -= deriv
            @inbounds jac[rxnarray[3,rxnind],rxnarray[1,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[1,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[1,rxnind]]*cs[rxnarray[3,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[2,rxnind]] -= deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[2,rxnind]] -= deriv
            @inbounds jac[rxnarray[3,rxnind],rxnarray[2,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[2,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]
            @inbounds jac[rxnarray[1,rxnind],rxnarray[3,rxnind]] -= deriv
            @inbounds jac[rxnarray[2,rxnind],rxnarray[3,rxnind]] -= deriv
            @inbounds jac[rxnarray[3,rxnind],rxnarray[3,rxnind]] -= deriv
            @inbounds _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,rxnarray[3,rxnind])
        end
    end
    k=krev
    if rxnarray[5,rxnind] == 0
        deriv = k
        @inbounds jac[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= deriv
        @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[4,rxnind])
    elseif rxnarray[6,rxnind] == 0
        if rxnarray[4,rxnind] == rxnarray[5,rxnind]
            @inbounds @fastmath deriv = 2.0*k*cs[rxnarray[4,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= 2.0*deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[4,rxnind])
        else
            @inbounds @fastmath deriv = k*cs[rxnarray[5,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[4,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[4,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[4,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[5,rxnind]] -= deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[5,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[5,rxnind])
        end
    else
        if rxnarray[4,rxnind]==rxnarray[5,rxnind] && rxnarray[4,rxnind]==rxnarray[6,rxnind]
            @inbounds @fastmath deriv = 3.0*k*cs[rxnarray[4,rxnind]]*cs[rxnarray[4,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= 3.0*deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[4,rxnind])
        elseif rxnarray[4,rxnind]==rxnarray[5,rxnind]
            @inbounds @fastmath deriv = 2.0*k*cs[rxnarray[4,rxnind]]*cs[rxnarray[6,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= 2.0*deriv
            @inbounds jac[rxnarray[6,rxnind],rxnarray[4,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[4,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[4,rxnind]]*cs[rxnarray[4,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[6,rxnind]] -= 2.0*deriv
            @inbounds jac[rxnarray[6,rxnind],rxnarray[6,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[6,rxnind])
        elseif rxnarray[5,rxnind]==rxnarray[6,rxnind]
            @inbounds @fastmath deriv = k*cs[rxnarray[5,rxnind]]*cs[rxnarray[5,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[4,rxnind]] -= 2.0*deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[4,rxnind])
            @inbounds @fastmath deriv = 2.0*k*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[5,rxnind]] -= deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[5,rxnind]] -= 2.0*deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[5,rxnind])
        elseif rxnarray[4,rxnind]==rxnarray[6,rxnind]
            @inbounds @fastmath deriv = 2.0*k*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= 2.0*deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[4,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[4,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[4,rxnind]]*cs[rxnarray[4,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[5,rxnind]] -= 2.0*deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[5,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[5,rxnind])
        else
            @inbounds @fastmath deriv = k*cs[rxnarray[5,rxnind]]*cs[rxnarray[6,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[4,rxnind]] -= deriv
            @inbounds jac[rxnarray[6,rxnind],rxnarray[4,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[4,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[4,rxnind]]*cs[rxnarray[6,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[5,rxnind]] -= deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[5,rxnind]] -= deriv
            @inbounds jac[rxnarray[6,rxnind],rxnarray[5,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[5,rxnind])
            @inbounds @fastmath deriv = k*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]
            @inbounds jac[rxnarray[4,rxnind],rxnarray[6,rxnind]] -= deriv
            @inbounds jac[rxnarray[5,rxnind],rxnarray[6,rxnind]] -= deriv
            @inbounds jac[rxnarray[6,rxnind],rxnarray[6,rxnind]] -= deriv
            @inbounds _spreadproductpartials!(jac,deriv,rxnarray,rxnind,rxnarray[6,rxnind])
        end
    end
end

@inline function _jacobianynswrtV!(jac::S,Vind::Int64,rxnarray::Array{Int64,2},rxnind::Int64,cs::Array{Float64,1},kf::Float64,krev::Float64) where {S<:AbstractArray}
    k=kf
    if rxnarray[2,rxnind]==0
        nothing
    elseif rxnarray[3,rxnind]== 0
        @inbounds @fastmath deriv = -k*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]
        @inbounds jac[rxnarray[1,rxnind],Vind] -= deriv
        @inbounds jac[rxnarray[2,rxnind],Vind] -= deriv
        _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,Vind)
    else
        @inbounds @fastmath deriv = -2.0*k*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]*cs[rxnarray[3,rxnind]]
        @inbounds jac[rxnarray[1,rxnind],Vind] -= deriv
        @inbounds jac[rxnarray[2,rxnind],Vind] -= deriv
        @inbounds jac[rxnarray[3,rxnind],Vind] -= deriv
        _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,Vind)
    end
    k=krev
    if rxnarray[5,rxnind]==0
        nothing
    elseif rxnarray[6,rxnind] == 0
        @inbounds @fastmath deriv = -k*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]
        @inbounds jac[rxnarray[4,rxnind],Vind] -= deriv
        @inbounds jac[rxnarray[5,rxnind],Vind] -= deriv
        _spreadproductpartials!(jac,deriv,rxnarray,rxnind,Vind)
    else
        @inbounds @fastmath deriv = -2.0*k*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]*cs[rxnarray[6,rxnind]]
        @inbounds jac[rxnarray[4,rxnind],Vind] -= deriv
        @inbounds jac[rxnarray[5,rxnind],Vind] -= deriv
        @inbounds jac[rxnarray[6,rxnind],Vind] -= deriv
        _spreadproductpartials!(jac,deriv,rxnarray,rxnind,Vind)
    end
end

"""
This function calculates the ns partials in jacobiany involving k derivatives. dkdx is either dkdni and dkdV. x is either ni or V.
"""
@inline function _jacobianykderiv!(jac::S,xind::Int64,dkfdx::Float64,dkrevdx::Float64,rxnarray::Array{Int64,2},rxnind::Int64,cs::Array{Float64,1},V::Float64) where {S<:AbstractArray}
    dkdx = dkfdx
    if rxnarray[2,rxnind] == 0
        @inbounds @fastmath deriv = dkdx*cs[rxnarray[1,rxnind]]*V
        @inbounds jac[rxnarray[1,rxnind],xind] -= deriv
        _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,xind)
    elseif rxnarray[3,rxnind] == 0
        @inbounds @fastmath deriv = dkdx*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]*V
        @inbounds jac[rxnarray[1,rxnind],xind] -= deriv
        @inbounds jac[rxnarray[2,rxnind],xind] -= deriv
        _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,xind)
    else
        @inbounds @fastmath deriv = dkdx*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]*cs[rxnarray[3,rxnind]]*V
        @inbounds jac[rxnarray[1,rxnind],xind] -= deriv
        @inbounds jac[rxnarray[2,rxnind],xind] -= deriv
        @inbounds jac[rxnarray[3,rxnind],xind] -= deriv
        _spreadreactantpartials!(jac,deriv,rxnarray,rxnind,xind)
    end
    dkdx = dkrevdx
    if rxnarray[5,rxnind] == 0
        @inbounds @fastmath deriv = dkdx*cs[rxnarray[4,rxnind]]*V
        @inbounds jac[rxnarray[4,rxnind],xind] -= deriv
        _spreadproductpartials!(jac,deriv,rxnarray,rxnind,xind)
    elseif rxnarray[6,rxnind] == 0
        @inbounds @fastmath deriv = dkdx*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]*V
        @inbounds jac[rxnarray[4,rxnind],xind] -= deriv
        @inbounds jac[rxnarray[5,rxnind],xind] -= deriv
        _spreadproductpartials!(jac,deriv,rxnarray,rxnind,xind)
    else
        @inbounds @fastmath deriv = dkdx*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]*cs[rxnarray[6,rxnind]]*V
        @inbounds jac[rxnarray[4,rxnind],xind] -= deriv
        @inbounds jac[rxnarray[5,rxnind],xind] -= deriv
        @inbounds jac[rxnarray[6,rxnind],xind] -= deriv
        _spreadproductpartials!(jac,deriv,rxnarray,rxnind,xind)
    end
end

function _dydttherm(dy::X,x::T,y::Q,p::W,t::Z,domain::D,interfaces::Q3,ind::T1) where {X,Q3<:AbstractArray,T<:Real,Q<:AbstractArray,Z<:Real,D<:AbstractDomain,T1<:Integer,W}
    v = [ i != ind ? convert(typeof(x),z) : x for (i,z) in enumerate(y)]
    return dydtreactor!(dy,v,t,domain,interfaces;p=p,sensitivity=false)
end

function jacobianytherm!(jac::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,ind::I,x::F,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:AbstractDomain,I<:Int64,F<:Float64}
    f(dy::X,x::Y) where {Y<:Real,X} = _dydttherm(dy,x,y,p,t,domain,interfaces,ind)
    jac[:,ind] = ForwardDiff.derivative(f,zeros(size(y)),x)
end

# function jacobianp!(d::W;cs::Q,V::Y,T::Y2,Us::Z3,Cvave::Z4,N::Z5,kfs::Z,krevs::X,wV::Q2,ratederiv::Q3) where {Q3,W<:Union{ConstantTPDomain,ConstantTVDomain},Z4<:Real,Z5<:Real,Z3<:AbstractArray,Q2<:AbstractArray,Q<:AbstractArray,Y2<:Real,Y<:Real,Z<:AbstractArray,X<:AbstractArray}
#     Nspcs = length(cs)
#     rxns = d.phase.reactions
#     Nrxns = length(rxns)
#     RTinv = 1.0/(R*T)
#     ratederiv .= 0.0
# 
#     for (j,rxn) in enumerate(rxns)
#         Nreact = length(rxn.reactantinds)
#         Nprod = length(rxn.productinds)
# 
#         if Nreact == 1
#             rind1 = rxn.reactantinds[1]
#             fderiv = cs[rind1]
#         elseif Nreact == 2
#             rind1,rind2 = rxn.reactantinds
#             fderiv = cs[rind1]*cs[rind2]
#         else
#             rind1,rind2,rind3 = rxn.reactantinds
#             fderiv = cs[rind1]*cs[rind2]*cs[rind3]
#         end
# 
#         if Nprod == 1
#             pind1 = rxn.productinds[1]
#             rderiv = krevs[j]/kfs[j]*cs[pind1]
#         elseif Nprod == 2
#             pind1,pind2 = rxn.productinds
#             rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]
#         else
#             pind1,pind2,pind3 = rxn.productinds
#             rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]*cs[pind3]
#         end
# 
#         flux = fderiv-rderiv
#         gderiv = rderiv*kfs[j]*RTinv
# 
#         deriv = zeros(Nspcs)
# 
#         deriv[rind1] += gderiv
#         if Nreact > 1
#             deriv[rind2] += gderiv
#             if Nreact > 2
#                 deriv[rind3] == gderiv
#             end
#         end
# 
#         deriv[pind1] -= gderiv
#         if Nprod > 1
#             deriv[pind2] -= gderiv
#             if Nprod > 2
#                 deriv[pind3] -= gderiv
#             end
#         end
# 
#         ratederiv[rind1,j] -= flux
#         ratederiv[rind1,Nrxns+1:Nrxns+Nspcs] .-= deriv
#         if Nreact > 1
#             ratederiv[rind2,j] -= flux
#             ratederiv[rind2,Nrxns+1:Nrxns+Nspcs] .-= deriv
#             if Nreact > 2
#                 ratederiv[rind3,j] -= flux
#                 ratederiv[rind3,Nrxns+1:Nrxns+Nspcs] .-= deriv
#             end
#         end
# 
#         ratederiv[pind1,j] += flux
#         ratederiv[pind1,Nrxns+1:Nrxns+Nspcs] .+= deriv
#         if Nprod > 1
#             ratederiv[pind2,j] += flux
#             ratederiv[pind2,Nrxns+1:Nrxns+Nspcs] .+= deriv
#             if Nprod > 2
#                 ratederiv[pind3,j] += flux
#                 ratederiv[pind3,Nrxns+1:Nrxns+Nspcs] .+= deriv
#             end
#         end
#     end
#     return V*ratederiv
# end
# 
# function jacobianp!(d::W; cs::Q,V::Y,T::Y2,Us::Z3,Cvave::Y3,N::Y2,kfs::Z,krevs::X,wV::Q2,ratederiv::Q3) where {W<:Union{ConstantVDomain,ParametrizedVDomain},Q3,Z3<:AbstractArray,Q<:AbstractArray,Q2<:AbstractArray,Y3<:Real,Y2<:Real,Y<:Real,Z<:AbstractArray,X<:AbstractArray}
#     Nspcs = length(cs)
#     rxns = d.phase.reactions
#     Nrxns = length(rxns)
#     RTinv = 1.0/(R*T)
#     ratederiv .= 0.0
# 
#     for (j,rxn) in enumerate(rxns)
#         Nreact = length(rxn.reactantinds)
#         Nprod = length(rxn.productinds)
# 
#         if Nreact == 1
#             rind1 = rxn.reactantinds[1]
#             fderiv = cs[rind1]
#         elseif Nreact == 2
#             rind1,rind2 = rxn.reactantinds
#             fderiv = cs[rind1]*cs[rind2]
#         else
#             rind1,rind2,rind3 = rxn.reactantinds
#             fderiv = cs[rind1]*cs[rind2]*cs[rind3]
#         end
# 
#         if Nprod == 1
#             pind1 = rxn.productinds[1]
#             rderiv = krevs[j]/kfs[j]*cs[pind1]
#         elseif Nprod == 2
#             pind1,pind2 = rxn.productinds
#             rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]
#         else
#             pind1,pind2,pind3 = rxn.productinds
#             rderiv = krevs[j]/kfs[j]*cs[pind1]*cs[pind2]*cs[pind3]
#         end
# 
#         flux = fderiv-rderiv
#         gderiv = rderiv*kfs[j]*RTinv
# 
#         deriv = zeros(Nspcs)
# 
#         deriv[rind1] += gderiv
#         if Nreact > 1
#             deriv[rind2] += gderiv
#             if Nreact > 2
#                 deriv[rind3] == gderiv
#             end
#         end
# 
#         deriv[pind1] -= gderiv
#         if Nprod > 1
#             deriv[pind2] -= gderiv
#             if Nprod > 2
#                 deriv[pind3] -= gderiv
#             end
#         end
# 
#         ratederiv[rind1,j] -= flux
#         ratederiv[rind1,Nrxns+1:Nrxns+Nspcs] .-= deriv
#         if Nreact > 1
#             ratederiv[rind2,j] -= flux
#             ratederiv[rind2,Nrxns+1:Nrxns+Nspcs] .-= deriv
#             if Nreact > 2
#                 ratederiv[rind3,j] -= flux
#                 ratederiv[rind3,Nrxns+1:Nrxns+Nspcs] .-= deriv
#             end
#         end
# 
#         ratederiv[pind1,j] += flux
#         ratederiv[pind1,Nrxns+1:Nrxns+Nspcs] .+= deriv
#         if Nprod > 1
#             ratederiv[pind2,j] += flux
#             ratederiv[pind2,Nrxns+1:Nrxns+Nspcs] .+= deriv
#             if Nprod > 2
#                 ratederiv[pind3,j] += flux
#                 ratederiv[pind3,Nrxns+1:Nrxns+Nspcs] .+= deriv
#             end
#         end
#     end
#     ratederiv *= V
#     #Temperature stuff
#     @views ratederiv[end,:] += (Us'*ratederiv[1:end-1,:])[1,:]
#     ratederiv[end,1:Nspcs] += wV
#     ratederiv[end,:] /= (N*Cvave)
#     return ratederiv
# end

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
# function jacobiany!(y::Array{T,1},t::T,domain::ConstantTPDomain,kfs::Array{T,1},krevs::Array{T,1},jac::P;zero::Bool=true) where {P<:AbstractArray,T<:Real,J<:Integer}
#     if zero
#         jac .= 0
#     end
#     N = sum(y)
#     V = N*R*domain.T/domain.P
#     cs = y./V
#     C = N/V
#     rxnarray = domain.rxnarray
#     Nrxns = size(rxnarray)[2]
#     Nspcs = length(y)
#     for i in 1:Nrxns
#         kf = kfs[i]
#         krev = krevs[i]
#         if rxnarray[2,i] == 0
#             jac[rxnarray[1,i],rxnarray[1,i]] -= kf
#             if rxnarray[5,i] == 0 
#                 jac[rxnarray[4,i],rxnarray[1,i]] += kf
#             elseif rxnarray[6,i] == 0
#                 jac[rxnarray[4,i],rxnarray[1,i]] += kf
#                 jac[rxnarray[5,i],rxnarray[1,i]] += kf
#             else
#                 jac[rxnarray[4,i],rxnarray[1,i]] += kf
#                 jac[rxnarray[5,i],rxnarray[1,i]] += kf
#                 jac[rxnarray[6,i],rxnarray[1,i]] += kf
#             end
#         elseif rxnarray[3,i] == 0
#             corr = -kf*cs[rxnarray[1,i]]*cs[rxnarray[2,i]]/C #correction for the partial of the volume term
#             if rxnarray[1,i] == rxnarray[2,i]
#                 deriv = 2*kf*cs[rxnarray[1,i]]
#                 jac[rxnarray[1,i],rxnarray[1,i]] -= 2.0*deriv
#                 for j in 1:Nspcs
#                     jac[rxnarray[1,i],j] -= 2.0*corr 
#                 end
#                 jac[rxnarray[4,i],rxnarray[1,i]] += deriv
#                 for j in 1:Nspcs 
#                     jac[rxnarray[4,i],j] += corr 
#                 end 
#                 if rxnarray[5,i] != 0
#                     jac[rxnarray[5,i],rxnarray[1,i]] += deriv 
#                     for j = 1:Nspcs 
#                         jac[rxnarray[5,i],j] += corr 
#                     end 
#                     if rxnarray[6,i] != 0 
#                         jac[rxnarray[6,i],rxnarray[1,i]] += deriv 
#                         for j = 1:Nspcs 
#                             jac[rxnarray[6,i],j] += corr 
#                         end 
#                     end 
#                 end 
#             else 
#                 #derivative with respect to reactant 1
#                 deriv = kf*cs[rxnarray[2,i]]
#                 jac[rxnarray[1,i],rxnarray[1,i]] -= deriv
#                 jac[rxnarray[2,i],rxnarray[1,i]] -= deriv
# 
#                 jac[rxnarray[4,i],rxnarray[1,i]] += deriv 
#                 if rxnarray[5,i] != 0 
#                     jac[rxnarray[5,i],rxnarray[1,i]] += deriv
#                     if rxnarray[6,i] != 0 
#                         jac[rxnarray[6,i],rxnarray[1,i]] += deriv 
#                     end 
#                 end 
# 
#                 #derivative with respect to reactant 2
#                 deriv = kf*cs[rxnarray[1,i]]
#                 jac[rxnarray[1,i],rxnarray[2,i]] -= deriv 
#                 jac[rxnarray[2,i],rxnarray[2,i]] -= deriv 
#                 for j = 1:Nspcs 
#                     jac[rxnarray[1,i],j] -= corr 
#                     jac[rxnarray[2,i],j] -= corr 
#                 end
#                 jac[rxnarray[4,i],rxnarray[2,j]] += deriv 
#                 if rxnarray[5,i] != 0 
#                     jac[rxnarray[5,i],rxnarray[2,i]] += deriv 
#                     for j = 1:Nspcs 
#                         jac[rxnarray[5,i],j] += corr 
#                     end 
#                     if rxnarray[6,i] != 0 
#                         jac[rxnarray[6,i],rxnarray[2,i]] += deriv 
#                         for j = 1:Nspcs 
#                             jac[rxnarray[6,i],j] += corr 
#                         end 
#                     end 
#                 end 
#             end
#         else
#             corr = -2.0*kf*cs[rxnarray[1,i]]*cs[rxnarray[2,i]]*cs[rxnarray[3,i]]/C
#             if (rxnarray[1,i] == rxnarray[2,i] && rxnarray[1,i] == rxnarray[3,i])
#                 deriv = 3.0*kf*cs[rxnarray[1,i]]*cs[rxnarray[1,i]]
#                 jac[rxnarray[1,i],rxnarray[1,i]] -= 3.0*deriv 
#                 for j = 1:Nspcs 
#                     jac[rxnarray[1,i],j] -= 3.0*corr 
#                 end 
#                 jac[rxnarray[4,i],rxnarray[1,i]] += deriv 
#                 for j = 1:Nspcs 
#                     jac[rxnarray[4,i],j] += corr 
#                 end 
#                 if rxnarray[5,i] != 0 
#                     jac[rxnarray[5,i],rxnarray[1,i]] += deriv 
#                     for j = 1:Nspcs 
#                         jac[rxnarray[5,i],j] += corr 
#                     end 
#                     if rxnarray[6,i] != 0 
#                         jac[rxnarray[6,i],rxnarray[1,i]] += deriv 
#                         for j = 1:Nspcs 
#                             jac[rxnarray[6,i],j] += corr 
#                         end 
#                     end 
#                 end 
#             elseif rxnarray[1,i] == rxnarray[2,i]
#                 #derivative with respect to reactant 1 
#                 deriv = 2.0*kf*cs[rxnarray[1,i]]*cs[rxnarray[3,i]]
#                 jac[rxnarray[1,i],rxnarray[1,i]] -= 2.0*deriv 
#                 jac[rxnarray[3,i],rxnarray[1,i]] -= deriv
# 
#                 jac[rxnarray[4,i],rxnarray[1,i]] += deriv 
#                 if rxnarray[5,i] != 0 
#                     jac[rxnarray[5,i],rxnarray[1,i]] += deriv 
#                     if rxnarray[6,i] != 0 
#                         jac[rxnarray[6,i],rxnarray[1,i]] += deriv 
#                     end 
#                 end 
# 
#                 #derivative with respect to reactant 3
#                 deriv = kf*cs[rxnarray]
# 
#             ind1,ind2,ind3 = rxn.reactantinds
#             corr = -2.0*kf*cs[ind1]*cs[ind2]*cs[ind3]/C
#             deriv = kf*cs[ind1]*cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind3,Nprod)
#             deriv = kf*cs[ind1]*cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*cs[ind3]*cs[ind2]
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
#         elseif Nprod == 2
#             ind1,ind2 = rxn.productinds
#             corr = -krev*cs[ind1]*cs[ind2]/C #correction for the partial of the volume term
#             deriv = krev*cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind1,Nreact)
#             for i in 1:length(domain.phase.species)
#                 jac[ind1,i] -= corr
#                 jac[ind2,i] -= corr
#                 spreadpartials!(jac,corr,rxn.reactantinds,i,Nreact)
#             end
#         elseif Nprod == 3
#             ind1,ind2,ind3 = rxn.productinds
#             corr = -2.0*krev*cs[ind1]*cs[ind2]*cs[ind3]/C
#             deriv = krev*cs[ind1]*cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind3,Nreact)
#             deriv = kf*cs[ind1]*cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*cs[ind3]*cs[ind2]
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
#     return jac
# end
# 
# function jacobiany!(y::Array{T,1},t::T,domain::ConstantTPDomain,kfs::Array{T,1},krevs::Array{T,1},jac::P;zero::Bool=true) where {P<:AbstractArray,T<:Real,J<:Integer}
#     if zero
#         jac .= 0
#     end
#     N = sum(y)
#     V = N*R*domain.T/domain.P
#     cs = y./V
#     C = N/V
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
#             corr = -kf*cs[ind1]*cs[ind2]/C #correction for the partial of the volume term
#             deriv = kf*cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*cs[ind2]
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
#             corr = -2.0*kf*cs[ind1]*cs[ind2]*cs[ind3]/C
#             deriv = kf*cs[ind1]*cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind3,Nprod)
#             deriv = kf*cs[ind1]*cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*cs[ind3]*cs[ind2]
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
#         elseif Nprod == 2
#             ind1,ind2 = rxn.productinds
#             corr = -krev*cs[ind1]*cs[ind2]/C #correction for the partial of the volume term
#             deriv = krev*cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind1,Nreact)
#             for i in 1:length(domain.phase.species)
#                 jac[ind1,i] -= corr
#                 jac[ind2,i] -= corr
#                 spreadpartials!(jac,corr,rxn.reactantinds,i,Nreact)
#             end
#         elseif Nprod == 3
#             ind1,ind2,ind3 = rxn.productinds
#             corr = -2.0*krev*cs[ind1]*cs[ind2]*cs[ind3]/C
#             deriv = krev*cs[ind1]*cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind3,Nreact)
#             deriv = kf*cs[ind1]*cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*cs[ind3]*cs[ind2]
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
#     return jac
# end
# function jacobiany!(y::Array{T,1},t::T,domain::ConstantTVDomain,kfs::Array{T,1},krevs::Array{T,1},jac::P;zero::Bool=true) where {P<:AbstractArray,T<:Real,J<:Integer}
#     if zero
#         jac .= 0
#     end
#     cs = y./domain.V
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
#             deriv = kf*cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind1,Nprod)
#         elseif Nreact == 3
#             ind1,ind2,ind3 = rxn.reactantinds
#             deriv = kf*state.cs[ind1]*cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind3,Nprod)
#             deriv = kf*cs[ind1]*cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind2,Nprod)
#             deriv = kf*cs[ind3]*cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             jac[ind3,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.productinds,ind1,Nprod)
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
#         elseif Nprod == 2
#             ind1,ind2 = rxn.productinds
#             deriv = krev*cs[ind1]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind1,Nreact)
#         elseif Nprod == 3
#             ind1,ind2,ind3 = rxn.productinds
#             deriv = krev*cs[ind1]*cs[ind2]
#             jac[ind1,ind3] -= deriv
#             jac[ind2,ind3] -= deriv
#             jac[ind3,ind3] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind3,Nreact)
#             deriv = kf*cs[ind1]*cs[ind3]
#             jac[ind1,ind2] -= deriv
#             jac[ind2,ind2] -= deriv
#             jac[ind3,ind2] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind2,Nreact)
#             deriv = kf*cs[ind3]*cs[ind2]
#             jac[ind1,ind1] -= deriv
#             jac[ind2,ind1] -= deriv
#             jac[ind3,ind1] -= deriv
#             spreadpartials!(jac,deriv,rxn.reactantinds,ind1,Nreact)
#         end
#     end
#     for ind in domain.constantspeciesinds
#         jac[ind,:] .= 0
#     end
#     return jac
# end
