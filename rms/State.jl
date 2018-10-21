using Parameters

abstract type AbstractState end
export AbstractState

struct EmptyState <: AbstractState end
export EmptyState

@with_kw mutable struct ConcentrationState{T,Q,S,N,Z,J,Y<:AbstactFloat} <: AbstractState
    cs::Array{T,1}
    T::Q
    P::S
    V::N
    C::J
    t::Y
end
export ConcentrationState
