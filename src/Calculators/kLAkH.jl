using Parameters

abstract type AbstractHenryLawConstant end
export AbstractHenryLawConstant

struct EmptyHenryLawConstant <: AbstractHenryLawConstant end
(eh::EmptyHenryLawConstant)(;T::N) where {N<:Number} = Inf
export EmptyHenryLawConstant

@with_kw struct ConstantHenryLawConstant{N<:Number} <: AbstractHenryLawConstant
    kH::N
end
(ch::ConstantHenryLawConstant)(;T::N) where {N<:Number} = ch.kH
export ConstantHenryLawConstant

struct TemperatureDependentHenryLawConstant{N<:Function} <: AbstractHenryLawConstant
    kH::N
end

function TemperatureDependentHenryLawConstant(;Ts::N1,kHs::N2) where {N1<:AbstractArray,N2<:AbstractArray}
    return TemperatureDependentHenryLawConstant(getspline(Ts,kHs))
end

(th::TemperatureDependentHenryLawConstant)(;T::N) where {N<:Number} = th.kH(T)
export TemperatureDependentHenryLawConstant

abstract type AbstractLiquidVolumetricMassTransferCoefficient end
export AbstractLiquidVolumetricMassTransferCoefficient

struct EmptyLiquidVolumetricMassTransferCoefficient <: AbstractLiquidVolumetricMassTransferCoefficient end
(eh::EmptyLiquidVolumetricMassTransferCoefficient)(;T::N) where {N<:Number} = 0
export EmptyLiquidVolumetricMassTransferCoefficient

@with_kw struct ConstantLiquidVolumetricMassTransferCoefficient{N<:Number} <: AbstractLiquidVolumetricMassTransferCoefficient
    kLA::N
end
(cv::ConstantLiquidVolumetricMassTransferCoefficient)(;T::N) where {N<:Number} = cv.kLA
export ConstantLiquidVolumetricMassTransferCoefficient

struct TemperatureDependentLiquidVolumetricMassTransferCoefficient{N<:Function} <: AbstractLiquidVolumetricMassTransferCoefficient 
    kLA::N
end

function TemperatureDependentLiquidVolumetricMassTransferCoefficient(;Ts::N1,kLAs::N2) where {N1<:AbstractArray,N2<:AbstractArray}
    return TemperatureDependentLiquidVolumetricMassTransferCoefficient(getspline(Ts,kLAs)) 
end

(tv::TemperatureDependentLiquidVolumetricMassTransferCoefficient)(;T::N) where {N<:Number} = tv.kLA(T)
export TemperatureDependentLiquidVolumetricMassTransferCoefficient