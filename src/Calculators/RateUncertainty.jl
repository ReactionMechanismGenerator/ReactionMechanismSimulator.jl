using Parameters

abstract type AbstractRateUncertainty end
export AbstractRateUncertainty

struct EmptyRateUncertainty <: AbstractRateUncertainty end
export EmptyRateUncertainty

@with_kw struct ConstantRateUncertainty <: AbstractRateUncertainty
    u::AbstractFloat = nothing
end
(c::ConstantRateUncertainty)(T::AbstractFloat=nothing) = c.u
export ConstantRateUncertainty

@with_kw struct EnergyRateUncertainty <: AbstractRateUncertainty
    dE::AbstractFloat=nothing
end
(er::EnergyRateUncertainty)(T::AbstractFloat=nothing) = exp(dE/(8.314*T))
export EnergyRateUncertainty
