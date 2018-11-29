using Parameters
import Base: length

using StaticArrays
abstract type AbstractReaction end
export AbstractReaction

@with_kw struct ElementaryReaction{T<:AbstractRate,Q<:Integer,V1<:StaticArray,V2<:StaticArray,V3<:StaticArray,V4<:StaticArray,V5<:StaticArray} <: AbstractReaction
    index::Q
    reactants::V1
    reactantinds::V2
    products::V3
    productinds::V4
    kinetics::T
    pairs::V5 = @SArray [@SArray [""]]
end
export ElementaryReaction

length(r::T) where {T<:AbstractReaction}= 1
export length
