using Parameters
using LinearAlgebra
using StaticArrays

abstract type AbstractDomain end
export AbstractDomain

abstract type AbstractConstantKDomain <: AbstractDomain  end
export AbstractConstantKDomain

abstract type AbstractVariableKDomain <: AbstractDomain end
export AbstractVariableKDomain

@with_kw struct ConstantTPDomain{N<:AbstractPhase,S<:Integer,W<:Real, W2<:Real, I<:Integer, Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Q #assumed to be in ascending order
    constantspeciesinds::Array{S,1}
    T::W
    P::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{UInt16,2}
    mu::W = 0.0
    diffusivity::Array{W,1} = Array{Float64,1}()
    jacobian::Array{W,2} = Array{Float64,2}(undef,(0,0))
    sensitivity::Bool = false
    jacuptodate::MArray{Tuple{1},Bool,1,1}=MVector(false)
    t::MArray{Tuple{1},W2,1,1}=MVector(0.0)
end
function ConstantTPDomain(;phase::E2,interfaces::Array{Q,1}=Array{EmptyInterface,1}(),initialconds::Dict{String,E},constantspecies::Array{String,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {E<:Real,E2<:AbstractPhase,Q<:AbstractInterface,W<:Real}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    if sensitivity
        y0 = zeros(length(phase.species)*(1+length(phase.species)+length(phase.reactions)))
    else
        y0 = zeros(length(phase.species))
    end
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end
    ns = y0
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase,T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if phase.diffusionlimited
        diffs = getfield.(phase.species,:diffusion)(T=T,mu=mu,P=P)
    else
        diffs = Array{typeof(T),1}()
    end
    C = P/(R*T)
    V = N*R*T/P
    kfs,krevs = getkfkrevs(phase=phase,T=T,P=P,C=C,N=N,ns=ns,Gs=Gs,diffs=diffs,V=V)
    if sparse
        jacobian=spzeros(typeof(T),length(phase.species),length(phase.species))
    else
        jacobian=zeros(typeof(T),length(phase.species),length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTPDomain(phase,interfaces,SVector(phase.species[1].index,phase.species[end].index),constspcinds,
        T,P,kfs,krevs,efficiencyinds,Gs,rxnarray,mu,diffs,jacobian,sensitivity,MVector(false),MVector(0.0)), y0
end
export ConstantTPDomain

@with_kw struct ConstantVDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Q #assumed to be in ascending order
    constantspeciesinds::Array{S,1}
    V::W
    rxnarray::Array{UInt16,2}
    jacobian::Array{W,2}
    sensitivity::Bool = false
    jacuptodate::MArray{Tuple{1},Bool,1,1}=MVector(false)
    t::MArray{Tuple{1},W2,1,1}=MVector(0.0)
end
function ConstantVDomain(;phase::Z,interfaces::Array{Q,1}=Array{EmptyInterface,1}(),initialconds::Dict{String,E},constantspecies::Array{String,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {E<:Real,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    N = sum(ns)
    if V == 0.0
        V = N*R*T/P
    elseif T == 0.0
        T = P*V/(R*N)
    elseif P == 0.0
        P = N*R*T/V
    else
        throw(error("ConstantVDomain overspecified with T,P and V"))
    end
    if sensitivity
        y0 = vcat(ns,T,zeros((length(ns)+1)*(length(ns)+length(phase.reactions))))
    else
        y0 = vcat(ns,T)
    end
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    if sparse
        jacobian=zeros(typeof(T),length(phase.species)+1,length(phase.species)+1)
    else
        jacobian=zeros(typeof(T),length(phase.species)+1,length(phase.species)+1)
    end
    rxnarray = getreactionindices(phase)
    return ConstantVDomain(phase,interfaces,SVector(phase.species[1].index,phase.species[end].index,phase.species[end].index+1),constspcinds,
    V,rxnarray,jacobian,sensitivity,MVector(false),MVector(0.0)), y0
end
export ConstantVDomain

@with_kw struct ConstantTVDomain{N<:AbstractPhase,S<:Integer,W<:Real, W2<:Real, I<:Integer, Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    interfaces::Array{AbstractInterface,1} = Array{AbstractInterface,1}()
    indexes::Q #assumed to be in ascending order
    constantspeciesinds::Array{S,1}
    T::W
    V::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{UInt16,2}
    mu::W = 0.0
    diffusivity::Array{W,1} = Array{Float64,1}()
    jacobian::Array{W,2} = Array{Float64,2}(undef,(0,0))
    sensitivity::Bool = false
    jacuptodate::MArray{Tuple{1},Bool,1,1}=MVector(false)
    t::MArray{Tuple{1},W2,1,1}=MVector(0.0)
end
function ConstantTVDomain(;phase::Z,interfaces::Array{Q,1}=Array{EmptyInterface,1}(),initialconds::Dict{String,E},constantspecies::Array{String,1}=Array{String,1}(),
    sparse=false,sensitivity=false) where {E<:Real, Z<:AbstractPhase,Q<:AbstractInterface,W<:Real}
    #set conditions and initialconditions
    T = 0.0
    V = 0.0
    P = 1.0e9
    if sensitivity
        y0 = zeros(length(phase.species)*(length(phase.species)+1+length(phase.reactions)))
    else
        y0 = zeros(length(phase.species))
    end
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end
    ns = y0
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase,T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if phase.diffusionlimited
        diffs = [x(T=T,mu=mu,P=P) for x in getfield.(phase.species,:diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    P = 1.0e9  #essentiallly assuming this is a liquid
    C = N/V
    kfs,krevs = getkfkrevs(phase=phase,T=T,P=P,C=C,N=N,ns=ns,Gs=Gs,diffs=diffs,V=V)
    if sparse
        jacobian=zeros(typeof(T),length(phase.species),length(phase.species))
    else
        jacobian=zeros(typeof(T),length(phase.species),length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTVDomain(phase,interfaces,SVector(phase.species[1].index,phase.species[end].index),constspcinds,
        T,V,kfs,krevs,efficiencyinds,Gs,rxnarray,mu,diffs,jacobian,sensitivity,MVector(false),MVector(0.0)), y0
end
export ConstantTVDomain


@inline function calcthermo(d::ConstantTPDomain{W,Y},y::J,t::Q) where {W<:IdealGas,Y<:Integer,J<:AbstractArray{Float64,1},Q<:Float64}
    if t != d.t[1]
        d.t[1] = t
        d.jacuptodate[1] = false
    end
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    V = N*d.T*R/d.P
    cs = ns./V
    C = N/V
    for ind in d.efficiencyinds #efficiency related rates may have changed
        d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V)
    end
    return ns,cs,d.T,d.P,V,C,N,d.mu,d.kfs,d.krevs,[],[],[],[],0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y},y::J,t::Q) where {W<:IdealGas,Y<:Integer,J<:AbstractArray,Q<:Real}
    if t != d.t[1]
        if isa(t,Float64)
            d.t[1] = t
        end
        d.jacuptodate[1] = false
    end
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    V = N*d.T*R/d.P
    cs = ns./V
    C = N/V
    kfs = convert(typeof(y),copy(d.kfs))
    krevs = convert(typeof(y),copy(d.krevs))
    for ind in d.efficiencyinds #efficiency related rates may have changed
        kfs[ind],krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V)
    end
    return ns,cs,d.T,d.P,V,C,N,d.mu,kfs,krevs,[],[],[],[],0.0
end

@inline function calcthermo(d::ConstantVDomain{W,Y},y::J,t::Q) where {W<:IdealGas,Y<:Integer,J<:AbstractArray,Q<:Real}
    if t != d.t[1]
        d.t[1] = t
        d.jacuptodate[1] = false
    end
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    N = sum(ns)
    cs = ns./d.V
    C = N/d.V
    P = C*R*T
    @views Us,Gs = calcenthalpyinternalgibbs(d.phase,T,P,d.V)[2:3]
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=mu,P=P)
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(phase=d.phase,T=T,P=P,C=C,N=N,ns=ns,Gs=Gs,diffs=diffs,V=d.V)
    @fastmath @inbounds f(spc::Species) = getHeatCapacity(spc.thermo,T)*ns[spc.index]
    @fastmath @inbounds Cvave = mapreduce(f,+,d.phase.species)/N - R
    return ns,cs,T,P,d.V,C,N,0.0,kfs,krevs,[],Us,Gs,diffs,Cvave
end

@inline function calcthermo(d::ConstantTVDomain{W,Y},y::J,t::Q) where {W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q<:Real}
    if t != d.t[1]
        d.t[1] = t
        d.jacuptodate[1] = false
    end
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.V
    C = N/d.V
    P = 1.0e9
    return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,[],[],[],[],0.0
end
export calcthermo

@inline function calcdomainderivatives!(d::Q,dydt::Array{Z7,1};T::Z4,Us::Array{Z,1},V::Z2,C::Z3,ns::Array{Z5,1},N::Z6,Cvave::Z8) where {Q<:AbstractDomain,Z8<:Real,Z7<:Real,W<:IdealGas,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5<:Real}
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
end

@inline function calcdomainderivatives!(d::ConstantVDomain{W,Y},dydt::Array{K,1};T::Z4,Us::Array{Z,1},V::Z2,C::Z3,ns::Array{Z5,1},N::Z6,Cvave::Z7) where {W<:IdealGas,Z7<:Real,K<:Real,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5<:Real}
    @views @fastmath @inbounds dydt[d.indexes[3]] = -dot(Us,dydt[d.indexes[1]:d.indexes[2]])/(N*Cvave) #divide by V to cancel ωV to ω
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
end
export calcdomainderivatives!

function getreactionindices(ig::Q) where {Q<:AbstractPhase}
    arr = zeros(UInt16,(6,length(ig.reactions)))
    for (i,rxn) in enumerate(ig.reactions)
        arr[1:length(rxn.reactantinds),i] = rxn.reactantinds
        arr[4:length(rxn.productinds)+3,i] = rxn.productinds
    end
    return arr
end
export getreactionindices
