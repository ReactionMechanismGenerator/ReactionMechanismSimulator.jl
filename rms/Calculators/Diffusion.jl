using Parameters

include("../Constants.jl")

abstract type AbstractDiffusivity end

struct EmptyDiffusivity <: AbstractDiffusivity end

@with_kw struct ConstantDiffusivity{N<:Number} <: AbstractDiffusivity
    D::N = nothing
end
(c::ConstantDiffusivity)(;T::AbstractFloat=0.0,P::AbstractFloat=0.0,mu::AbstractFloat=0.0) = c.D

@with_kw struct StokesDiffusivity{N<:Number} <: AbstractDiffusivity
    r::N
end

StokesDiffusivity(;V::N=nothing) where {N<:Number} = StokesDiffusivity(r=(75.0*V/(pi*Na))^(1.0/3.0)/100.0)

(sd::StokesDiffusivity)(;T::N,mu::Q,P::R=0.0) where {N,R,Q<:Number} = kB*T/(6*pi*mu*sd.r)
