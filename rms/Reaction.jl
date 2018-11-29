using Parameters
import Base: length,show,print,println
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

getrxnstr(rxn::T) where {T<:AbstractReaction} = join([join(getfield.(rxn.reactants,:name),"+"),join(getfield.(rxn.products,:name),"+")],"<=>")
show(io::IO,rxn::T) where {T<:AbstractReaction} = print(io,getrxnstr(rxn))
print(rxn::T) where {T<:AbstractReaction} = print(getrxnstr(rxn))
println(rxn::T) where {T<:AbstractReaction} = print(string(getrxnstr(rxn),"\n"))
length(r::T) where {T<:AbstractReaction}= 1
export length
