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
export spcindex

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
getT(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ParametrizedVDomain},K<:Real,Q,G,L} = bsol.sol(t)[bsol.domain.indexes[3]]
getT(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ParametrizedTConstantVDomain,ParametrizedTPDomain},K<:Real,Q,G,L} = bsol.domain.T(t)
export getT
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.N(t)*getT(bsol,t)*R /bsol.domain.P
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain,ParametrizedTConstantVDomain},K<:Real,Q,G,L} = bsol.domain.V
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedVDomain,K<:Real,Q,G,L} = bsol.domain.V(t)
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedTPDomain,K<:Real,Q,G,L} = bsol.N(t)*getT(bsol,t)*R /bsol.domain.P(t)
export getV
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.domain.P
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantTVDomain,ParametrizedTConstantVDomain},K<:Real,Q,G,L} = 1.0e6
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ParametrizedVDomain},K<:Real,Q,G,L} = bsol.N(t)*R*getT(bsol,t)/getV(bsol,t)
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedTPDomain,K<:Real,Q,G,L} = bsol.domain.P(t)
export getP
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.domain.P/(R*bsol.domain.T)
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain,ParametrizedTConstantVDomain},K<:Real,Q,G,L} = bsol.N(t)/bsol.domain.V
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedVDomain,K<:Real,Q,G,L} = bsol.N(t)/bsol.domain.V(t)
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedTPDomain,K<:Real,Q,G,L} = bsol.domain.P(t)/(R*bsol.domain.T(t))
export getC
"""
calculates the rates of production/loss at a given time point
this outputs a sparse matrix of  num reactions xnum species containing the production/loss
rate of that species associated with that reaction
"""
function rops(bsol::Q,t::X) where {Q<:Simulation,X<:Real}
    ropmat = zeros(length(bsol.domain.phase.reactions),length(bsol.domain.phase.species))
    cs,kfs,krevs = calcthermo(bsol.domain,bsol.sol(t),t)[[2,9,10]]
    @simd for i in 1:length(bsol.domain.phase.reactions)
        rxn = bsol.domain.phase.reactions[i]
        R = getrate(rxn,cs,kfs,krevs)
        for ind in rxn.productinds
            ropmat[i,ind] += R
        end
        for ind in rxn.reactantinds
            ropmat[i,ind] -= R
        end
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
    cs,kfs,krevs = calcthermo(bsol.domain,bsol.sol(t),t)[[2,9,10]]
    ind = findfirst(isequal(name),getfield.(bsol.domain.phase.species,:name))
    @assert !isa(ind,Nothing) "species $name not in species array"
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
    cs,kfs,krevs = calcthermo(bsol.domain,bsol.sol(t),t)[[2,9,10]]
    return [getrate(rxn,cs,kfs,krevs) for rxn in bsol.domain.phase.reactions]
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
