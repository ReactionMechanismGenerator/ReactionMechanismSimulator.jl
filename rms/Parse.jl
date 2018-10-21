using Unitful
using YAML
include("Tools.jl")
include("Calculators.jl")
include("Species.jl")
include("Reaction.jl")
module Calc
    include("Calculators.jl")
end
module Spc
    include("Species.jl")
end
module Rxn
    include("Reaction.jl")
end

const unitsdict = Dict()
const allowedfcnlist = vcat(names(Calc),names(Spc),names(Rxn))

function parseentry(q,units=nothing)
    """
    parses single YAML lines into strings, numbers and arrays of Numbers in SI units
    """
    if isa(q,String)
        q = parsestring(q)
    end

    if isa(q,String)
        return q
    elseif isa(q,Quantity)
        return upreferred(q).val
    elseif isa(q,Number) || isa(q,AbstractArray)
        if units == nothing
            return q
        else
            return upreferred(Quantity(q,units)).val
        end
    else
        tq = typeof(q)
        throw(error("unable to parse $q of type $tq"))
    end

end

function parsestring(s)
    """
    Detects if a given string is actually a quantity input using metaprogramming
    the string is parsed by julia's parser and then analyzed to check if it is in
    the format of a quantity input 
    """
    ex = Meta.parse(s)
    if !isa(ex,Symbol) && !isa(ex,String) && length(ex.args) == 3 && ex.args[1] == :* &&
        (isa(ex.args[2],Number) || (ex.args[2].head == :vec && all([isa(x,Number) for x in ex.args[2].args]))) &&
        ex.args[3].args[1] == Symbol("@u_str") &&
        isa(ex.args[3].args[2],LineNumberNode) && ex.args[3].args[2].line == 1 &&
        ex.args[3].args[2].file == :none && isa(ex.args[3].args[3],String)
        return eval(ex)
    else
        return s
    end
end
