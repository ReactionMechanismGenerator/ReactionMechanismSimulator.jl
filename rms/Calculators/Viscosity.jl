using Parameters

abstract type AbstractViscosity end
export AbstractViscosity

struct EmptyViscosity <: AbstractViscosity end
export EmptyViscosity

@with_kw struct ConstantViscosity{T<:Number}
    mu::T
end

(cv::ConstantViscosity)(T::Q) where {Q<:Number} = cv.mu
export ConstantViscosity

@with_kw struct RiedelViscosity{Q,V,X,Z,N<:Number}
    A::Q
    B::V
    C::X
    D::Z
    E::N
end

(rv::RiedelViscosity)(T::Q) where {Q<:Number} = exp(rv.A+rv.B/T+rv.C*log(T)+rv.D*T^rv.E)
export RiedelViscosity
