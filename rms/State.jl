using Parameters
import Base: length
abstract type AbstractState end
export AbstractState

struct EmptyState <: AbstractState end
export EmptyState

@with_kw mutable struct MolarState{W<:AbstractFloat} <: AbstractState
    ns::Array{W,1}
    cs::Array{W,1} = Array{Float64,1}()
    T::W = 0.0
    P::W = 0.0
    V::W = 0.0
    C::W = 0.0
    t::W = 0.0
    N::W = 0.0
    mu::W = 0.0
    Gs::Array{W,1} = Array{Float64,1}()
    Hs::Array{W,1} = Array{Float64,1}()
    Us::Array{W,1} = Array{Float64,1}()
    Ss::Array{W,1} = Array{Float64,1}()
    Cps::Array{W,1} = Array{Float64,1}()
    diffusivity::Array{W,1} = Array{Float64,1}()
end
#Note that the main MolarState constructor is in PhaseState.jl
#to avoid importing AbstractPhase here
export MolarState

length(m::MolarState) = 1
export length

iterate(m::MolarState) = m
export iterate

Broadcast.broadcastable(st::MolarState) = Ref(st)
export broadcastable
