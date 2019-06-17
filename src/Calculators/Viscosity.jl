using Parameters

abstract type AbstractViscosity end
export AbstractViscosity

struct EmptyViscosity <: AbstractViscosity end
export EmptyViscosity

@with_kw struct ConstantViscosity{T<:AbstractFloat} <: AbstractViscosity
    mu::T
end

(cv::ConstantViscosity)(T::Q) where {Q<:AbstractFloat} = cv.mu
export ConstantViscosity

@with_kw struct RiedelViscosity{N<:AbstractFloat} <: AbstractViscosity
    A::N
    B::N
    C::N
    D::N
    E::N
end

(rv::RiedelViscosity)(T::Q) where {Q<:Number} = @fastmath exp(rv.A+rv.B/T+rv.C*log(T)+rv.D*T^rv.E)
export RiedelViscosity
