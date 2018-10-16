using Parameters

abstract type AbstractRateUncertainty end

struct EmptyRateUncertainty <: AbstractRateUncertainty end

@with_kw struct ConstantRateUncertainty <: AbstractRateUncertainty
    u::AbstractFloat = nothing
end
(c::ConstantRateUncertainty)(T::AbstractFloat=nothing) = c.u

@with_kw struct EnergyRateUncertainty <: AbstractRateUncertainty
    dE::AbstractFloat=nothing
end
(er::EnergyRateUncertainty)(T::AbstractFloat=nothing) = exp(dE/(8.314*T))
