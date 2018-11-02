using Parameters
using SpecialFunctions
include("Constants.jl")
include("Species.jl")
include("Reaction.jl")
include("Solvent.jl")
include("State.jl")

function MolarState(inputdict::Dict{Z,W},ph::Q) where {Z<:String,W<:AbstractFloat,Q<:AbstractPhase}
    n = length(ph.species)
    ns = zeros(n)
    D = Dict()
    for (key,val) in inputdict
        if !(key in ["T","P","t","V"])
            ns[ph.spcdict[key]] = val
        elseif isa(key,String)
            D[Symbol(key)] = val
        end
    end
    D[:ns] = ns
    ms = MolarState(;D...)
    ms.cs = zeros(n)
    ms.Gs = zeros(n)
    ms.Hs = zeros(n)
    ms.Us = zeros(n)
    if isa(ph,IdealDiluteSolution) #volume must be defined for constant V reactors
        @assert ms.V != 0.0 "Volume must be defined for IdealDiluteSolution Phase"
    end
    if :solvent in fieldnames(typeof(ph))
        ms.mu = ph.solvent.mu(ms.T)
    end
    if ph.diffusionlimited == true
        ms.diffusivity = map(x->x.diffusion(T=ms.T,mu=ms.mu),ph.species)
    else
        ms.diffusivity = zeros(n)
    end
    return ms
end
export MolarState

function recalcgibbs!(ph::T,st::MolarState) where {T<:AbstractPhase}
    map!(x->getGibbs(x.thermo,st.T),st.Gs,ph.species)
end
export recalcgibbs!

function recalcgibbsandinternal!(ph::T,st::MolarState) where {T<:IdealPhase}
    map!(x->getEnthalpy(x,st.T),st.Hs,ph.species)
    st.Us = st.Hs .- st.P*st.V
    st.Gs = st.Hs .- st.T.*map(x->getEntropy(x,st.T),ph.species)
end
export recalcgibbsandinternal!

function recalcgibbsandinternal!(ph::IdealGas,st::MolarState)
    map!(x->getEnthalpy(x.thermo,st.T),st.Hs,ph.species)
    st.Us = st.Hs .- R*st.T
    st.Gs = st.Hs .- st.T.*map(x->getEntropy(x.thermo,st.T),ph.species)
end
export recalcgibbsandinternal!

function getkf(rxn::ElementaryReaction,ph::T,st::MolarState) where {T<:AbstractPhase}
    return rxn.kinetics(T=st.T,P=st.P,C=st.C)
end
export getkf

"""
Calculates the diffusion limited rate coefficient
for 1 spc returns Inf
for 2 spc calculates using the Smolchowski equation
for >2 spc calculates using the Generalized Smolchowski equation
Equations from Flegg 2016
"""
function getDiffusiveRate(spcs::Q,st::MolarState) where {Q<:AbstractArray}
    if length(spcs) == 1
        return Inf
    elseif length(spcs) == 2
        kf = 4.0*Base.pi*(st.diffusivity[1]+st.diffusivity[2])*(spcs[1].radius+spcs[2].radius)*Na
    else
        N = length(spcs)
        a = (3.0*length(spcs)-5.0)/2.0
        Dinv = 1.0./st.diffusivity
        Dbar = 1.0./reverse(cumsum(Dinv))
        Dhat = st.diffusivity .+ Dbar
        deltaN = sum(Dinv)/sum(sum([[1.0/(st.diffusivity[i]*st.diffusivity[m]) for m in 1:N-1 if i>m] for i in 2:N]))
        kf = prod(Dhat[2:end].^1.5)*4*Base.pi^(a+1)/gamma(a)*(sum(getfield.(spcs,:radius))/sqrt(deltaN))^(2*a)*Na^(N-1)
    end
    return kf
end

function getKc(rxn::ElementaryReaction,ph::T,st::MolarState) where {T<:AbstractPhase}
    return exp(-(sum(st.Gs[rxn.productinds])-sum(st.Gs[rxn.reactantinds]))/(R*st.T))*(1.0e5/(R*st.T))^(length(rxn.productinds)-length(rxn.reactantinds))
end
export getKc

function getrate(rxn::ElementaryReaction,ph::IdealGas,st::MolarState)
    kf = getkf(rxn,ph,st)
    Kc = getKc(rxn,ph,st)
    krev = kf/Kc
    return kf*prod(st.cs[rxn.reactantinds]) - krev*prod(st.cs[rxn.productinds])
end
export getrate

function addreactionratecontribution!(y::Array{Q,1},rxn::ElementaryReaction,ph::IdealGas,st::MolarState) where {Q<:Number,T<:Integer}
    R = getrate(rxn,ph,st)
    for ind in rxn.reactantinds
        y[ind] -= R
    end
    for ind in rxn.productinds
        y[ind] += R
    end
end
export addreactionratecontribution!
