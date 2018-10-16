using Parameters

abstract type AbstractThermoUncertainty end

struct EmptyThermoUncertainty <: AbstractThermoUncertainty end

@with_kw struct ConstantThermoUncertainty <: AbstractThermoUncertainty
    u::AbstractFloat = nothing
end
(c::ConstantThermoUncertainty)(T::AbstractFloat=nothing) = c.u
