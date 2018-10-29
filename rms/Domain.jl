using Parameters
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
end
ConstantTPDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{EmptyInterface,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} = ConstantTPDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index])
export ConstantTPDomain

@with_kw struct ConstantVDomain{T<:AbstractState,N<:AbstractPhase,S<:Integer} <: AbstractVariableKDomain
    state::T
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Array{S,1}
end
ConstantVDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{AbstractInterface,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} = ConstantVDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index,0])
export ConstantVDomain

@with_kw struct ConstantTVDomain{T<:AbstractState,N<:AbstractPhase,S<:Integer} <: AbstractConstantKDomain
    state::T
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Array{S,1}
end
ConstantTVDomain(;state::T,phase::N,interfaces::Array{Q,1}=Array{EmptyInterface,1}()) where {T<:AbstractState,N<:AbstractPhase,Q<:AbstractInterface} = ConstantTVDomain(state,phase,interfaces,[phase.species[1].index,phase.species[end].index])
export ConstantTVDomain


function calcthermo!(d::ConstantTPDomain{Z,W,Y},y::T,t::Q) where {Z<:MolarState,W<:IdealGas,Y<:Integer,T<:AbstractArray,Q<:AbstractFloat}
    d.state.ns = y[d.indexes[1]:d.indexes[2]]
    d.state.t = t
    d.state.N = sum(d.state.ns)
    d.state.V = d.state.N*d.state.T*R/d.state.P
    d.state.cs = d.state.ns./d.state.V
    d.state.C = d.state.N/d.state.V
end

function calcthermo!(d::ConstantVDomain{Z,W,Y},y::T,t::Q) where {Z<:MolarState,W<:IdealGas,Y<:Integer,T<:AbstractArray,Q<:AbstractFloat}
    d.state.ns = y[d.indexes[1]:d.indexes[2]]
    d.state.T = y[d.indexes[3]]
    d.state.t = t
    d.state.N = sum(d.state.ns)
    d.state.cs = d.state.ns./d.state.V
    d.state.C = d.state.N/d.state.V
    d.state.P = d.state.C*R*d.state.T
    recalcgibbsandinternal!(d.phase,d.state)
end

function calcthermo!(d::ConstantTVDomain{Z,W,Y},y::T,t::Q) where {Z<:MolarState,W<:IdealDiluteSolution,Y<:Integer,T<:AbstractArray,Q<:AbstractFloat}
    d.state.ns = y[d.indexes[1]:d.indexes[2]]
    d.state.t = t
    d.state.N = sum(d.state.ns)
    d.state.cs = d.state.ns./d.state.V
    d.state.C = d.state.N/d.state.V
end
export calcthermo!


end
export ConstantTVDomain 
