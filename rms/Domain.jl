using Parameters
using LinearAlgebra
include("Constants.jl")
include("Phase.jl")
include("State.jl")
include("Interface.jl")

abstract type AbstractDomain end
export AbstractDomain

abstract type AbstractConstantKDomain <: AbstractDomain  end
export AbstractConstantKDomain

abstract type AbstractVariableKDomain <: AbstractDomain end
export AbstractVariableKDomain

@with_kw struct ConstantTPDomain{T<:AbstractState,N<:AbstractPhase,S<:Integer} <: AbstractConstantKDomain
    state::T
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Array{S,1}
    constantspeciesinds::Array{S,1}
end
function ConstantTPDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{EmptyInterface,1}(),constantspecies::Array{String,1}=Array{String,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface}
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    return ConstantTPDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index],constspcinds)
end
export ConstantTPDomain

@with_kw struct ConstantVDomain{T<:AbstractState,N<:AbstractPhase,S<:Integer} <: AbstractVariableKDomain
    state::T
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Array{S,1}
    constantspeciesinds::Array{S,1}
end
function ConstantVDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{AbstractInterface,1}(),constantspecies::Array{String,1}=Array{String,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface}
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    return ConstantVDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index,phase.species[end].index+1],constspcinds)
end
export ConstantVDomain

@with_kw struct ConstantTVDomain{T<:AbstractState,N<:AbstractPhase,S<:Integer} <: AbstractConstantKDomain
    state::T
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Array{S,1}
    constantspeciesinds::Array{S,1}
end
function ConstantTVDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{EmptyInterface,1}(),constantspecies::Array{String,1}=Array{String,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface}
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    return ConstantTVDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index],constspcinds)
end
export ConstantTVDomain


function calcthermo!(d::ConstantTPDomain{Z,W,Y},y::T,t::Q) where {Z<:MolarState,W<:IdealGas,Y<:Integer,T<:AbstractArray,Q<:AbstractFloat}
    d.state.ns = y[d.indexes[1]:d.indexes[2]]
    d.state.t = t
    d.state.N = sum(d.state.ns)
    @fastmath d.state.V = d.state.N*d.state.T*R/d.state.P
    @fastmath d.state.cs = d.state.ns./d.state.V
    @fastmath d.state.C = d.state.N/d.state.V
end

function calcthermo!(d::ConstantVDomain{Z,W,Y},y::T,t::Q) where {Z<:MolarState,W<:IdealGas,Y<:Integer,T<:AbstractArray,Q<:AbstractFloat}
    @inbounds d.state.ns = y[d.indexes[1]:d.indexes[2]]
    @inbounds d.state.T = y[d.indexes[3]]
    d.state.t = t
    d.state.N = sum(d.state.ns)
    @fastmath d.state.cs = d.state.ns./d.state.V
    @fastmath d.state.C = d.state.N/d.state.V
    @fastmath d.state.P = d.state.C*R*d.state.T
    recalcgibbsandinternal!(d.phase,d.state)
end

function calcthermo!(d::ConstantTVDomain{Z,W,Y},y::T,t::Q) where {Z<:MolarState,W<:IdealDiluteSolution,Y<:Integer,T<:AbstractArray,Q<:AbstractFloat}
    @inbounds d.state.ns = y[d.indexes[1]:d.indexes[2]]
    d.state.t = t
    d.state.N = sum(d.state.ns)
    @fastmath d.state.cs = d.state.ns./d.state.V
    @fastmath d.state.C = d.state.N/d.state.V
    d.state.mu = d.phase.solvent.mu(d.state.T)
end
export calcthermo!

function calcdomainderivatives!(d::T,dydt::Array{N,1}) where {T<:AbstractDomain,N<:AbstractFloat}
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
end

function calcdomainderivatives!(d::ConstantVDomain{Z,W,Y},dydt::Array{N,1}) where {Z<:MolarState,W<:IdealGas,N<:AbstractFloat,Y<:Integer}
    @fastmath @inbounds Cpave = mapreduce(x->getHeatCapacity(x.thermo,d.state.T)*d.state.ns[x.index],+,d.phase.species)/d.state.N
    @fastmath Cvave = Cpave-R
    @inbounds @fastmath dydt[d.indexes[3]] = -d.state.Us'*(dydt[d.indexes[1]:d.indexes[2]]/d.state.V)/(d.state.C*Cvave) #divide by V to cancel ωV to ω
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
end
export calcdomainderivatives!
