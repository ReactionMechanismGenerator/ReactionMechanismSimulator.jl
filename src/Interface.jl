using LinearAlgebra

abstract type AbstractInterface end
export AbstractInterface

abstract type AbstractBoundaryInterface <: AbstractInterface end
export AbstractBoundaryInterface

abstract type AbstractInternalInterface <: AbstractInterface end
export AbstractInternalInterface

struct EmptyInterface <: AbstractInterface end
export EmptyInterface

@with_kw struct IdealGasCatalystInterface{T<:AbstractPhase,N<:AbstractPhase,Q<:AbstractReaction} <: AbstractInternalInterface
    gas::T
    catalyst::N
    reactions::Array{Q,1}
"""
construct the stochiometric matrix for the reactions crossing both domains and the reaction molecule # change
"""
function getstoichmatrix(domain1,domain2,rxns)
    M = spzeros(length(rxns),length(domain1.phase.species)+length(domain2.phase.species))
    Nrp1 = zeros(length(rxns))
    Nrp2 = zeros(length(rxns))
    N1 = length(domain1.phase.species)
    spcs1 = domain1.phase.species
    spcs2 = domain2.phase.species
    for (i,rxn) in enumerate(rxns)
        Nrp1[i] = Float64(length([x for x in rxn.products if x in spcs1]) - length([x for x in rxn.reactants if x in spcs1]))
        Nrp2[i] = Float64(length([x for x in rxn.products if x in spcs2]) - length([x for x in rxn.reactants if x in spcs2]))
        for (j,r) in enumerate(rxn.reactants)
            isfirst = true
            ind = findfirst(isequal(r),domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r),domain2.phase.species)
            end
            M[i,isfirst ? ind : ind+N1] += 1
        end
        for (j,r) in enumerate(rxn.products)
            isfirst = true
            ind = findfirst(isequal(r),domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r),domain2.phase.species)
            end
            M[i,isfirst ? ind : ind+N1] -= 1
        end
    end
    return M,Nrp1,Nrp2
end

function getinterfacereactioninds(domain1,domain2,reactions)
    indices = zeros(Int64,(6,length(reactions)))
    N1 = length(domain1.phase.species)
    for (i,rxn) in enumerate(reactions)
        for (j,r) in enumerate(rxn.reactants)
            isfirst = true
            ind = findfirst(isequal(r),domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r),domain2.phase.species)
            end
            indices[j,i] = isfirst ? ind : ind+N1
        end
        for (j,r) in enumerate(rxn.products)
            isfirst = true
            ind = findfirst(isequal(r),domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r),domain2.phase.species)
            end
            indices[j+3,i] = isfirst ? ind : ind+N1
        end
    end
    return indices
end
end
export IdealGasCatalystInterface

struct Inlet{Q<:Real,S,V<:AbstractArray,U<:Real,X<:Real} <: AbstractBoundaryInterface
    domain::S
    y::V
    F::Function
    T::U
    P::X
    H::Q
end

function Inlet(domain::V,conddict::Dict{String,X},F::Function) where {V,X<:Real,B<:Real}
    y = makespcsvector(domain.phase,conddict)
    T = conddict["T"]
    P = conddict["P"]
    yout = y./sum(y)
    H = dot(getEnthalpy.(getfield.(domain.phase.species,:thermo),T),yout)
    return Inlet(domain,yout,F,T,P,H)
end

export Inlet

struct Outlet{V} <: AbstractBoundaryInterface
    domain::V
    F::Function
end
export Outlet
