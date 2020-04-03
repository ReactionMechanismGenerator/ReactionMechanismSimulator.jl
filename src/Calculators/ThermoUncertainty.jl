using Parameters

abstract type AbstractThermoUncertainty end
export AbstractThermoUncertainty

struct EmptyThermoUncertainty <: AbstractThermoUncertainty end
export EmptyThermoUncertainty

@with_kw struct ConstantThermoUncertainty <: AbstractThermoUncertainty
    u::AbstractFloat = nothing
end
(c::ConstantThermoUncertainty)(T::AbstractFloat=nothing) = c.u
