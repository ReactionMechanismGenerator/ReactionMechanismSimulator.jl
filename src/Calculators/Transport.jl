abstract type AbstractTransportModel end
export AbstractTransportModel

struct EmptyTransportModel <: AbstractTransportModel end

@with_kw struct TransportModel{N1<:AbstractPotential} <: AbstractTransportModel
    m::Float64 # molecular mass
    potential::N1
    dipolemoment::Float64
    polarizability::Float64
end
export TransportModel
