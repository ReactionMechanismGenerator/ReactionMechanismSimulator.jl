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
