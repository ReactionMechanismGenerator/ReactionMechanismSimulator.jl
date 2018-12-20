using DifferentialEquations
import DifferentialEquations.DiffEqBase: AbstractODESolution, HermiteInterpolation,AbstractDiffEqInterpolation

abstract type AbstractSimulation end
export AbstractSimulation

struct Simulation{Q<:AbstractODESolution,W<:AbstractDomain,L<:AbstractArray,G<:Function,G2<:AbstractArray} <: AbstractSimulation
    sol::Q
    domain::W
    names::L
    N::G
    Ns::G2
end

function Simulation(sol::Q,domain::W) where {Q<:AbstractODESolution,W<:AbstractDomain}
    names = getfield.(domain.phase.species,:name)
    Ns = sum(hcat(sol.interp.u...)[domain.indexes[1]:domain.indexes[2],:],dims=1)
    Nderivs = sum(hcat(sol.interp.du...)[domain.indexes[1]:domain.indexes[2],:],dims=1)
    N = HermiteInterpolation(sol.interp.t,Ns,Nderivs)
    F(t::T) where {T<:Real} = N(t,nothing,Val{0},sol.prob.p,:left)
    return Simulation(sol,domain,names,F,Ns)
end

export Simulation

length(p::T) where {T<:AbstractSimulation} = 1
export length

iterate(p::T) where {T<:AbstractSimulation} = p
export iterate

Broadcast.broadcastable(p::T) where {T<:AbstractSimulation} = Ref(p)
export broadcastable

spcindex(bsol::Z,name::Q) where {Z<:Simulation,Q<:AbstractString} = findfirst(isequal(name),getfield.(bsol.domain.phase.species,:name))

function molefractions(bsol::Q,name::W,t::E) where {Q<:AbstractSimulation, W<:String, E<:Real}
    @assert name in bsol.names
    ind = findfirst(isequal(name),bsol.names)
    return bsol.sol(t)[ind]/bsol.N(t)
end

function molefractions(bsol::Q, t::E) where {Q<:AbstractSimulation,E<:Real}
    return bsol.sol(t)[bsol.domain.indexes[1]:bsol.domain.indexes[2]]./bsol.N(t)
end

function molefractions(bsol::Q) where {Q<:AbstractSimulation}
    @views return hcat(bsol.sol.u...)[bsol.domain.indexes[1]:bsol.domain.indexes[2],:]./bsol.Ns
end

export molefractions

getT(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantTPDomain,ConstantTVDomain},K<:Real,Q,G,L} = bsol.domain.T
getT(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantVDomain,K<:Real,Q,G,L} = bsol.sol(t)[bsol.domain.indexes[3]]
export getT
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.N(t)*getT(bsol,t)*R /bsol.domain.P
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain},K<:Real,Q,G,L} = bsol.domain.V
export getV
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.domain.P
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantTVDomain,K<:Real,Q,G,L} = 1.0e6
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantVDomain,K<:Real,Q,G,L} = bsol.N(t)*R*getT(bsol,t)/getV(bsol,t)
export getP
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.domain.P/(R*bsol.domain.T)
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain},K<:Real,Q,G,L} = bsol.N(t)/bsol.domain.V
export getC
"""
calculates the rates of production/loss at a given time point
this outputs a sparse matrix of  num reactions xnum species containing the production/loss
rate of that species associated with that reaction
"""
function rops(bsol::Q,t::X) where {Q<:Simulation,X<:Real}
    ropmat = spzeros(length(bsol.domain.phase.reactions),length(bsol.domain.phase.species))
    xs = molefractions(bsol,t)
    T = getT(bsol,t)
    V = getV(bsol,t)
    P = getP(bsol,t)
    if :Gs in fieldnames(typeof(bsol.domain))
        Gs = bsol.domain.Gs
    else
        Gs = calcgibbs(bsol.domain.phase,T)
    end
    if :solvent in fieldnames(typeof(bsol.domain.phase)) && typeof(bsol.domain.phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if bsol.domain.phase.diffusionlimited
        diffs = getfield.(phase.species,:diffusion)(T=T,mu=mu,P=P)
    else
        diffs = Array{typeof(T),1}()
    end
    kfs,krevs = getkfkrevs(phase=bsol.domain.phase,V=V,T=T,P=P,C=1.0/V,N=1.0,ns=xs,Gs=Gs,diffs=diffs)
    cs = xs./V
    for (i,rxn) in enumerate(bsol.domain.phase.reactions)
        R = getrate(rxn,cs,kfs,krevs)
        ropmat[i,ind] = R*(count(isequal(ind),rxn.productinds)-count(isequal(ind),rxn.reactantinds))
    end
    return ropmat
end

"""
calculates the rates of production/loss at a given time point for a single species
this outputs a sparse vector of length num reactions containing the production/loss
rate associated with that reaction for the given species
"""
function rops(bsol::Y,name::X,t::Z) where {Y<:Simulation, X<:AbstractString, Z<:Real}
    rop = spzeros(length(bsol.domain.phase.reactions))
    xs = molefractions(bsol,t)
    T = getT(bsol,t)
    V = getV(bsol,t)
    P = getP(bsol,t)
    if :Gs in fieldnames(typeof(bsol.domain))
        Gs = bsol.domain.Gs
    else
        Gs = calcgibbs(bsol.domain.phase,T)
    end
    if :solvent in fieldnames(typeof(bsol.domain.phase)) && typeof(bsol.domain.phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if bsol.domain.phase.diffusionlimited
        diffs = getfield.(phase.species,:diffusion)(T=T,mu=mu,P=P)
    else
        diffs = Array{typeof(T),1}()
    end
    kfs,krevs = getkfkrevs(phase=bsol.domain.phase,V=V,T=T,P=P,C=1.0/V,N=1.0,ns=xs,Gs=Gs,diffs=diffs)
    cs = xs./V
    ind = findfirst(isequal(name),getfield.(bsol.domain.phase.species,:name))
    for (i,rxn) in enumerate(bsol.domain.phase.reactions)
        c = 0
        R = getrate(rxn,cs,kfs,krevs)
        c -= count(isequal(ind),rxn.reactantinds)
        c += count(isequal(ind),rxn.productinds)
        if c != 0
            rop[i] = c*R
        end
    end
    return rop
end
export rops

function getconcentrationsensitivity(bsol::Simulation{Q,W,L,G}, numerator::String, denominator::String, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain},K<:Real,Q,G,L}
    @assert numerator in bsol.names
    @assert denominator in bsol.names
    indnum = findfirst(isequal(numerator),bsol.names)
    inddeno = findfirst(isequal(denominator),bsol.names)
    Nvars = bsol.domain.indexes[end]-bsol.domain.indexes[1]+1
    Nrxns = length(bsol.domain.phase.reactions)
    arr = bsol.sol(t)
    s = arr[Nvars+Nrxns*Nvars+(inddeno-1)*Nvars+indnum]
    val = s/arr[indnum] #constant volume
    if t == 0
        return 0.0
    else
        return val
    end
end

function getconcentrationsensitivity(bsol::Simulation{Q,W,L,G}, numerator::String, denominator::String, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L}
    @assert numerator in bsol.names
    @assert denominator in bsol.names
    indnum = findfirst(isequal(numerator),bsol.names)
    inddeno = findfirst(isequal(denominator),bsol.names)
    Nvars = bsol.domain.indexes[end]-bsol.domain.indexes[1]+1
    Nrxns = length(bsol.domain.phase.reactions)
    arr = bsol.sol(t)
    s = arr[Nvars+Nrxns*Nvars+(inddeno-1)*Nvars+indnum]
    V = getV(bsol,t)
    c = arr[indnum]/V
    val =  (s-c*sum(arr[Nvars+Nrxns*Nvars+(inddeno-1)*Nvars+1:Nvars+Nrxns*Nvars+inddeno*Nvars])*R*bsol.domain.T/bsol.domain.P)/(c*V) #known T and P
    if t == 0
        return 0.0
    else
        return val
    end
end

function getconcentrationsensitivity(bsol::Simulation{Q,W,L,G}, numerator::String, denominator::Z, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain},K<:Real,Z<:Integer,Q,G,L}
    @assert numerator in bsol.names
    indnum = findfirst(isequal(numerator),bsol.names)
    inddeno = denominator
    Nvars = bsol.domain.indexes[end]-bsol.domain.indexes[1]+1
    Nrxns = length(bsol.domain.phase.reactions)
    arr = bsol.sol(t)
    s = arr[Nvars+(inddeno-1)*Nvars+indnum]
    T = getT(bsol,t)
    P = getP(bsol,t)
    C = getC(bsol,t)
    k = bsol.domain.phase.reactions[inddeno].kinetics(T=T,P=P,C=C)
    val = s*k/arr[indnum] #constant volume
    if t == 0
        return 0.0
    else
        return val
    end
end

function getconcentrationsensitivity(bsol::Simulation{Q,W,L,G}, numerator::String, denominator::Z, t::K) where {W<:ConstantTPDomain,K<:Real,Z<:Integer,Q,G,L}
    @assert numerator in bsol.names
    indnum = findfirst(isequal(numerator),bsol.names)
    inddeno = denominator
    Nvars = bsol.domain.indexes[end]-bsol.domain.indexes[1]+1
    Nrxns = length(bsol.domain.phase.reactions)
    arr = bsol.sol(t)
    s = arr[Nvars+(inddeno-1)*Nvars+indnum]
    V = getV(bsol,t)
    T = getT(bsol,t)
    P = getP(bsol,t)
    C = getC(bsol,t)
    c = arr[indnum]/V
    k = bsol.domain.phase.reactions[inddeno].kinetics(T=T,P=P,C=C)
    val = k*(s-c*sum(arr[Nvars+(inddeno-1)*Nvars+1:Nvars+inddeno*Nvars])*R*bsol.domain.T/bsol.domain.P)/(c*V) #known T and P
    if t == 0
        return 0.0
    else
        return val
    end
end

export getconcentrationsensitivity

"""
calculate the rates of all reactions at time t
"""
function rates(bsol::Q,t::X) where {Q<:Simulation,X<:Real}
    rates = zeros(length(bsol.domain.phase.reactions))
    xs = molefractions(bsol,t)
    T = getT(bsol,t)
    V = getV(bsol,t)
    P = getP(bsol,t)
    if :Gs in fieldnames(typeof(bsol.domain))
        Gs = bsol.domain.Gs
    else
        Gs = calcgibbs(bsol.domain.phase,T)
    end
    if :solvent in fieldnames(typeof(bsol.domain.phase)) && typeof(bsol.domain.phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if bsol.domain.phase.diffusionlimited
        diffs = getfield.(phase.species,:diffusion)(T=T,mu=mu,P=P)
    else
        diffs = Array{typeof(T),1}()
    end
    kfs,krevs = getkfkrevs(phase=bsol.domain.phase,V=V,T=T,P=P,C=1.0/V,N=1.0,ns=xs,Gs=Gs,diffs=diffs)
    cs = xs./V
    for i in 1:length(bsol.domain.phase.reactions)
        rates[i] =  getrate(bsol.domain.phase.reactions[i],cs,kfs,krevs)
    end
    return rates
end

"""
calculate the rates of all reactions at given times ts
defaults to using bsol.sol.t if ts is not supplied
"""
function rates(bsol::Q;ts::X=Array{Float64,1}()) where {Q<:Simulation,X<:AbstractArray}
    if length(ts) == 0
        ts = bsol.sol.t
    end
    return hcat([rates(bsol,t) for t in ts]...)
end

export rates
