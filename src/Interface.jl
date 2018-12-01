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

struct Inlet{Q<:Real,V<:AbstractArray,U<:Real,X<:Real} <: AbstractBoundaryInterface
    F::V
    T::U
    P::X
    dUdt::Q
end

function Inlet(phase::V,conddict::Dict{String,X},F::B) where {V<:AbstractPhase,X<:Real,B<:Real}
    y = makespcsvector(phase,conddict)
    T = conddict["T"]
    P = conddict["P"]
    Fout = F.*y./sum(y)
    dUdt = dot(getEnthalpy.(getfield.(ph.species,:thermo),T),Fout)
    return Inlet(Fout,T,P,dUdt)
end

export Inlet

struct Outlet{V<:Real} <: AbstractBoundaryInterface
    F::V
end
export Outlet
