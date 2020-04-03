using LinearAlgebra

abstract type AbstractInterface end
export AbstractInterface

abstract type AbstractBoundaryInterface <: AbstractInterface end
export AbstractBoundaryInterface

abstract type AbstractInternalInterface <: AbstractInterface end
export AbstractInternalInterface

struct EmptyInterface <: AbstractInterface end
export EmptyInterface

@with_kw struct IdealGasCatalystInterface{T<:AbstractPhase,N<:AbstractPhase,Q<:AbstractReaction} <: AbstractInternalInterface
    gas::T
    catalyst::N
    reactions::Array{Q,1}
end
export IdealGasCatalystInterface

struct Inlet{Q<:Real,S,V<:AbstractArray,U<:Real,X<:Real} <: AbstractBoundaryInterface
    domain::S
    y::V
    F::Function
    T::U
    P::X
    H::Q
end

function Inlet(domain::V,conddict::Dict{String,X},F::Function) where {V,X<:Real,B<:Real}
    y = makespcsvector(domain.phase,conddict)
    T = conddict["T"]
    P = conddict["P"]
    yout = y./sum(y)
    H = dot(getEnthalpy.(getfield.(domain.phase.species,:thermo),T),yout)
    return Inlet(domain,yout,F,T,P,H)
end

export Inlet

struct Outlet{V} <: AbstractBoundaryInterface
    domain::V
    F::Function
end
export Outlet
