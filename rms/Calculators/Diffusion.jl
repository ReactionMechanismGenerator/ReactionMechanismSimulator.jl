using Parameters

abstract type AbstractDiffusivity end
export AbstractDiffusivity

struct EmptyDiffusivity <: AbstractDiffusivity end
export EmptyDiffusivity

@with_kw struct ConstantDiffusivity{N<:Number} <: AbstractDiffusivity
    D::N = nothing
end
(c::ConstantDiffusivity)(;T::AbstractFloat=0.0,P::AbstractFloat=0.0,mu::AbstractFloat=0.0) = c.D
export ConstantDiffusivity

@with_kw struct StokesDiffusivity{N<:Number} <: AbstractDiffusivity
    r::N
end
(sd::StokesDiffusivity)(;T::N,mu::Q,P::R=0.0) where {N,R,Q<:Number} = @fastmath kB*T/(6*Base.pi*mu*sd.r)
export StokesDiffusivity
