using Parameters
import Base: length,show,print,println
using StaticArrays
abstract type AbstractReaction end
export AbstractReaction

@with_kw struct ElementaryReaction{T<:AbstractRate,Q<:Integer,V1<:AbstractArray,V2<:AbstractArray,V3<:AbstractArray,V4<:AbstractArray,V5<:AbstractArray,V6<:AbstractArray,V7<:AbstractArray,V8<:AbstractArray,V9<:AbstractArray} <: AbstractReaction
    index::Q
    reactants::V1
    reactantinds::V2
    products::V3
    productinds::V4
    kinetics::T
    electronchange::Int64 = 0
    radicalchange::Int64 = -100
    reversible::Bool = true
    pairs::V5 = [[""]]
    fragmentbasedreactants::V6 = []
    fragmentbasedproducts::V7 = []
    fragmentbasedreactantinds::V8 = []
    fragmentbasedproductinds::V9 = []
    comment::String = ""
end
export ElementaryReaction

getrxnstr(rxn::T) where {T<:AbstractReaction} = join([join(getfield.(rxn.reactants,:name),"+"),join(getfield.(rxn.products,:name),"+")],"<=>")
export getrxnstr
show(io::IO,rxn::T) where {T<:AbstractReaction} = print(io,getrxnstr(rxn))
print(rxn::T) where {T<:AbstractReaction} = print(getrxnstr(rxn))
println(rxn::T) where {T<:AbstractReaction} = print(string(getrxnstr(rxn),"\n"))

"""
matches reactants to the products that most resemble them
returns a StaticArray of length 2 static arrays containing the names of each pair
"""
function getpairs(rxn::T) where {T<:AbstractReaction}
    if length(rxn.reactants) == 1
        if length(rxn.products) == 1
            return [ [rxn.reactants[1].name,rxn.products[1].name ] ]
        elseif length(rxn.products) == 2
            return [ [rxn.reactants[1].name,rxn.products[1].name ],  [rxn.reactants[1].name,rxn.products[2].name ]]
        elseif length(rxn.products) == 3
            return [ [rxn.reactants[1].name,rxn.products[1].name ], [rxn.reactants[1].name,rxn.products[2].name ], [rxn.reactants[1].name,rxn.products[3].name ]]
        end
    elseif length(rxn.reactants) == 2
        if length(rxn.products) == 1
            return [ [rxn.reactants[1].name,rxn.products[1].name ], [rxn.reactants[2].name,rxn.products[1].name ]]
        else
            return choosepairs(rxn)
        end
    elseif length(rxn.reactants) == 3
        if length(rxn.products) == 1
            return [ [rxn.reactants[1].name,rxn.products[1].name ], [rxn.reactants[2].name,rxn.products[1].name ], [rxn.reactants[3].name,rxn.products[1].name ]]
        else
            return choosepairs(rxn)
        end
    else
        return choosepairs(rxn)
    end
end
export getpairs

"""
returns a negative number representing compositional similarity between two species
this only uses the number of bonds and number of each atom type
"""
function getsimilarity(spc1::T,spc2::T2) where {T<:AbstractSpecies,T2<:AbstractSpecies}
    s = -abs(spc1.bondnum - spc2.bondnum)
    keys1 = keys(spc1.atomnums)
    keys2 = keys(spc2.atomnums)
    for x in keys1
        if x in keys2
            s -= abs(spc1.atomnums[x]-spc2.atomnums[x])
        else
            s -= spc1.atomnums[x]
        end
    end
    for y in keys2
        if y in keys1
            continue
        else
            s -= spc2.atomnums[y]
        end
    end
    return s
end

function choosepairs(rxn::T) where {T<:AbstractReaction}
    pairs = []
    for (i,r) in enumerate(rxn.reactants)
        sims = [getsimilarity(r,p) for p in rxn.products]
        if length(sims) == 2 && sims[1] == sims[2] #add both if the same
            push!(pairs,[r.name,rxn.products[1].name])
            push!(pairs,[r.name,rxn.products[2].name])
        else
            ind = argmax(sims)
            push!(pairs,[r.name,rxn.products[ind].name])
        end
    end
    return pairs
end

length(r::T) where {T<:AbstractReaction}= 1
export length
