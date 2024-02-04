using Parameters
using LinearAlgebra
using StaticArrays
using Calculus
using SciMLBase
using ForwardDiff
using Tracker
using ReverseDiff
using Logging

abstract type AbstractDomain end
export AbstractDomain

abstract type AbstractConstantKDomain <: AbstractDomain end
export AbstractConstantKDomain

abstract type AbstractVariableKDomain <: AbstractDomain end
export AbstractVariableKDomain

mutable struct ConstantTPDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::W
    P::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{Int64,2}
    mu::W
    diffusivity::Array{W,1}
    jacobian::Array{W,2}
    sensitivity::Bool
    alternativepformat::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end

function ConstantTPDomain(; phase::E2, initialconds::Dict{X,X2}, constantspecies::Array{X3,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false) where {E<:Real,E2<:AbstractPhase,Q<:AbstractInterface,W<:Real,X,X2,X3}
    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    y0 = zeros(length(phase.species) + 1)
    spnames = [x.name for x in phase.species]
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end


    @assert T != 0.0
    @assert P != 0.0
    ns = y0[1:end-1]
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase, T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if phase.diffusionlimited
        diffs = [x(T=T, mu=mu, P=P) for x in getfield.(phase.species, :diffusion)]
    else
        diffs = Array{typeof(T),1}()
    end
    C = P / (R * T)
    V = N * R * T / P
    y0[end] = V
    kfs, krevs = getkfkrevs(phase, T, P, C, N, ns, Gs, diffs, V, 0.0, 0.0)
    kfsp = deepcopy(kfs)
    for ind in efficiencyinds
        kfsp[ind] = 1.0
    end
    p = vcat(deepcopy(Gs), kfsp)
    if sparse
        jacobian = spzeros(typeof(T), length(phase.species), length(phase.species))
    else
        jacobian = zeros(typeof(T), length(phase.species), length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTPDomain(phase, [1, length(phase.species), length(phase.species) + 1], [1, length(phase.species) + length(phase.reactions)], constspcinds,
        T, P, kfs, krevs, efficiencyinds, Gs, rxnarray, mu, diffs, jacobian, sensitivity, false, MVector(false), MVector(0.0), p, Dict("V" => length(phase.species) + 1)), y0, p
end
export ConstantTPDomain

struct ConstantVDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    V::W
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ConstantVDomain(; phase::Z, initialconds::Dict{X,E}, constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false) where {E,X,X2,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert V != 0.0 || (T != 0.0 && P != 0.0)
    N = sum(ns)
    if V == 0.0
        V = N * R * T / P
    elseif T == 0.0
        T = P * V / (R * N)
    elseif P == 0.0
        P = N * R * T / V
    else
        throw(error("ConstantVDomain overspecified with T,P and V"))
    end
    y0 = vcat(ns, T, P)
    p = vcat(zeros(length(phase.species)), ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian = zeros(typeof(T), length(phase.species) + 2, length(phase.species) + 2)
    else
        jacobian = zeros(typeof(T), length(phase.species) + 2, length(phase.species) + 2)
    end
    rxnarray = getreactionindices(phase)
    return ConstantVDomain(phase, [1, length(phase.species), length(phase.species) + 1, length(phase.species) + 2], [1, length(phase.species) + length(phase.reactions)], constspcinds,
        V, efficiencyinds, rxnarray, jacobian, sensitivity, MVector(false), MVector(0.0), p, Dict("T" => length(phase.species) + 1, "P" => length(phase.species) + 2)), y0, p
end
export ConstantVDomain

struct ConstantPDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    P::W
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ConstantPDomain(; phase::Z, initialconds::Dict{X,E}, constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false) where {E,X,X2,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert P != 0.0 || (T != 0.0 && V != 0.0)
    N = sum(ns)
    if P == 0.0
        P = N * R * T / V
    elseif T == 0.0
        T = P * V / (R * N)
    elseif V == 0.0
        V = N * R * T / P
    else
        throw(error("ConstantPDomain overspecified with T,P and V"))
    end
    y0 = vcat(ns, T, V)
    p = vcat(zeros(length(phase.species)), ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian = zeros(typeof(T), length(phase.species) + 2, length(phase.species) + 2)
    else
        jacobian = zeros(typeof(T), length(phase.species) + 2, length(phase.species) + 2)
    end
    rxnarray = getreactionindices(phase)
    return ConstantPDomain(phase, [1, length(phase.species), length(phase.species) + 1, length(phase.species) + 2], [1, length(phase.species) + length(phase.reactions)], constspcinds,
        P, efficiencyinds, rxnarray, jacobian, sensitivity, MVector(false), MVector(0.0), p, Dict("T" => length(phase.species) + 1, "V" => length(phase.species) + 2)), y0, p
end
export ConstantPDomain

struct ParametrizedTPDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray,FT<:Function,FP<:Function} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::FT
    P::FP
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ParametrizedTPDomain(; phase::Z, initialconds::Dict{X,Any}, constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false) where {X,X2,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ts = 0.0
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        elseif key == "ts"
            ts = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert V != 0.0 || (T != 0.0 && P != 0.0)
    if isa(T, AbstractArray)
        Tfcn = getspline(ts, T)
    elseif isa(T, Function)
        Tfcn = T
    else
        throw(error("ParametrizedTPDomain must take \"T\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end
    if isa(P, AbstractArray)
        Pfcn = getspline(ts, P)
    elseif isa(P, Function)
        Pfcn = P
    else
        throw(error("ParametrizedTPDomain must take \"P\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end

    N = sum(ns)
    V = N * R * Tfcn(0.0) / Pfcn(0.0)
    y0 = zeros(length(phase.species) + 1)
    y0[1:length(phase.species)] = ns
    y0[length(phase.species)+1] = V
    p = vcat(zeros(length(phase.species)), ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian = zeros(typeof(V), length(phase.species) + 1, length(phase.species) + 1)
    else
        jacobian = zeros(typeof(V), length(phase.species) + 1, length(phase.species) + 1)
    end
    rxnarray = getreactionindices(phase)
    return ParametrizedTPDomain(phase, [1, length(phase.species), length(phase.species) + 1], [1, length(phase.species) + length(phase.reactions)], constspcinds,
        Tfcn, Pfcn, efficiencyinds, rxnarray, jacobian, sensitivity, MVector(false), MVector(0.0), p, Dict("V" => length(phase.species) + 1)), y0, p
end
export ParametrizedTPDomain

struct ParametrizedVDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray,FV<:Function} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    V::FV
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ParametrizedVDomain(; phase::Z, initialconds::Dict{X,Any}, constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false) where {X,X2,E<:Real,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ts = Array{Float64,1}()
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    @assert "V" in keys(initialconds)
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        elseif key == "ts"
            ts = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert isa(V, Function) || isa(V, AbstractArray)
    if isa(V, AbstractArray)
        Vfcn = getspline(ts, V)
    elseif isa(V, Function)
        Vfcn = V
    else
        throw(error("ParametrizedVDomain must take \"V\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end
    N = sum(ns)
    if T == 0.0
        T = P * Vfcn(0.0) / (R * N)
    elseif P == 0.0
        P = N * R * T / Vfcn(0.0)
    else
        ns *= (P * Vfcn(0.0) / (R * T)) / sum(ns) #automatically scale down moles if pressure specified
    end
    y0 = vcat(ns, T, P)
    p = vcat(zeros(length(phase.species)), ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian = zeros(typeof(T), length(phase.species) + 2, length(phase.species) + 2)
    else
        jacobian = zeros(typeof(T), length(phase.species) + 2, length(phase.species) + 2)
    end
    rxnarray = getreactionindices(phase)
    return ParametrizedVDomain(phase, [1, length(phase.species), length(phase.species) + 1, length(phase.species) + 2], [1, length(phase.species) + length(phase.reactions)], constspcinds,
        Vfcn, efficiencyinds, rxnarray, jacobian, sensitivity, MVector(false), MVector(0.0), p, Dict("T" => length(phase.species) + 1, "P" => length(phase.species) + 2)), y0, p
end
export ParametrizedVDomain

struct ParametrizedPDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray,FP<:Function} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    P::FP
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ParametrizedPDomain(; phase::Z, initialconds::Dict{X,Any}, constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false) where {X,X2,E<:Real,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ts = Array{Float64,1}()
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    @assert "P" in keys(initialconds)
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        elseif key == "ts"
            ts = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert isa(P, Function) || isa(P, AbstractArray)
    if isa(P, AbstractArray)
        Pfcn = getspline(ts, P)
    elseif isa(P, Function)
        Pfcn = P
    else
        throw(error("ParametrizedPDomain must take \"P\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end
    N = sum(ns)
    if T == 0.0
        T = Pfcn(0.0) * V / (R * N)
    elseif V == 0.0
        V = N * R * T / Pfcn(0.0)
    else
        ns *= (Pfcn(0.0) * V / (R * T)) / sum(ns) #automatically scale down moles if volume specified
    end
    y0 = vcat(ns, T, V)
    p = vcat(zeros(length(phase.species)), ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian = zeros(typeof(T), length(phase.species) + 2, length(phase.species) + 2)
    else
        jacobian = zeros(typeof(T), length(phase.species) + 2, length(phase.species) + 2)
    end
    rxnarray = getreactionindices(phase)
    return ParametrizedPDomain(phase, [1, length(phase.species), length(phase.species) + 1, length(phase.species) + 2], [1, length(phase.species) + length(phase.reactions)], constspcinds,
        Pfcn, efficiencyinds, rxnarray, jacobian, sensitivity, MVector(false), MVector(0.0), p, Dict("T" => length(phase.species) + 1, "V" => length(phase.species) + 2)), y0, p
end
export ParametrizedPDomain

mutable struct ConstantTVDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::W
    V::W
    phi::W
    d::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    kfsnondiff::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{Int64,2}
    mu::W
    diffusivity::Array{W,1}
    jacobian::Array{W,2}
    sensitivity::Bool
    alternativepformat::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ConstantTVDomain(; phase::Z, initialconds::Dict{X,E}, constantspecies::Array{X2,1}=Array{String,1}(),
    sparse=false, sensitivity=false) where {E,X,X2,Z<:AbstractPhase,Q<:AbstractInterface,W<:Real}
    #set conditions and initialconditions
    T = 0.0
    V = 0.0
    P = 1.0e8
    phi = 0.0
    d = 0.0
    y0 = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            throw(error("ConstantTVDomain cannot specify P"))
        elseif key == "V"
            V = val
        elseif key == "Phi"
            phi = val
        elseif key == "d"
            d = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end
    @assert T != 0.0
    @assert V != 0.0
    ns = y0
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase, T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if phase.diffusionlimited
        diffs = [x(T=T, mu=mu, P=P) for x in getfield.(phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    P = 1.0e8  #essentiallly assuming this is a liquid
    C = N/V
    dGrxns = -(phase.stoichmatrix*Gs).+phase.electronchange.*(phi*F)
    kfs,krevs = getkfkrevs(phase,T,P,C,N,ns,Gs,diffs,V,phi,d)
    kfsnondiff = getkfs(phase,T,P,C,ns,V,phi,dGrxns,d)
    p = vcat(Gs,kfsnondiff)
    if sparse
        jacobian = zeros(typeof(T), length(phase.species), length(phase.species))
    else
        jacobian = zeros(typeof(T), length(phase.species), length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTVDomain(phase,[1,length(phase.species)],[1,length(phase.species)+length(phase.reactions)],constspcinds,
        T,V,phi,d,kfs,krevs,kfsnondiff,efficiencyinds,Gs,rxnarray,mu,diffs,jacobian,sensitivity,false,MVector(false),MVector(0.0),p,Dict{String,Int64}()), y0, p
end
export ConstantTVDomain

struct ParametrizedTConstantVDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray,FT<:Function} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::FT
    V::W
    phi::W
    d::W
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ParametrizedTConstantVDomain(; phase::IdealDiluteSolution, initialconds::Dict{X,X3}, constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false) where {X,X2,X3,Q<:AbstractInterface}
    #set conditions and initialconditions
    T = 0.0
    P = 1.0e8 #essentiallly assuming this is a liquid
    V = 0.0
    phi = 0.0
    d = 0.0
    ts = Array{Float64,1}()
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    @assert "V" in keys(initialconds)
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            throw(error("ParametrizedTConstantVDomain cannot specify P"))
        elseif key == "V"
            V = val
        elseif key == "ts"
            ts = val
        elseif key == "Phi"
            phi = val
        elseif key == "d"
            d = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    if isa(T, AbstractArray)
        Tfcn = getspline(ts, T)
    elseif isa(T, Function)
        Tfcn = T
    else
        throw(error("ParametrizedTConstantVDomain must take \"T\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end
    N = sum(ns)
    y0 = zeros(length(phase.species))
    y0[1:length(phase.species)] = ns
    p = vcat(zeros(length(phase.species)), ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian = zeros(typeof(V), length(phase.species) + 1, length(phase.species) + 1)
    else
        jacobian = zeros(typeof(V), length(phase.species) + 1, length(phase.species) + 1)
    end
    rxnarray = getreactionindices(phase)
    return ParametrizedTConstantVDomain(phase,[1,length(phase.species)],[1,length(phase.species)+length(phase.reactions)],constspcinds,
    Tfcn,V,phi,d,efficiencyinds,rxnarray,jacobian,sensitivity,MVector(false),MVector(0.0),p,Dict{String,Int64}()), y0, p
end
export ParametrizedTConstantVDomain

mutable struct ConstantTAPhiDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::W
    A::W
    phi::W
    d::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{Int64,2}
    mu::W
    diffusivity::Array{W,1}
    jacobian::Array{W,2}
    sensitivity::Bool
    alternativepformat::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ConstantTAPhiDomain(; phase::E2, initialconds::Dict{X,X2}, constantspecies::Array{X3,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false, stationary::Bool=false) where {E<:Real,E2<:AbstractPhase,W<:Real,X,X2,X3}
    #set conditions and initialconditions
    T = 0.0
    A = 0.0
    phi = 0.0 #default 0.0
    d = 0.0
    y0 = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "A"
            A = val
        elseif key == "Phi"
            phi = val
        elseif key == "d"
            d = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end

    @assert A != 0.0
    @assert T != 0.0
    ns = y0
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase, T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    C = 0.0 #this currently shouldn't matter here, on a surface you shouldn't have pdep
    kfs,krevs = getkfkrevs(phase,T,0.0,C,N,ns,Gs,[],A,phi,d)
    p = vcat(Gs,kfs)
    if sparse
        jacobian = spzeros(typeof(T), length(phase.species), length(phase.species))
    else
        jacobian = zeros(typeof(T), length(phase.species), length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTAPhiDomain(phase,[1,length(phase.species)],[1,length(phase.species)+length(phase.reactions)],constspcinds,
        T,A,phi,d,kfs,krevs,efficiencyinds,Gs,rxnarray,mu,Array{Float64,1}(),jacobian,sensitivity,false,MVector(false),MVector(0.0),p,Dict{String,Int64}()), y0, p
end
export ConstantTAPhiDomain

"""
Fragment based constant T, rho domain is designed to simulate film growth with fragment based model.
It simulates the growth of mass of a swollen film caused by the chemical reaction between 
the functional groups (fragments) on the solid and 
the liquid phase reactants encapsulated in the film.
mass: the mass of the solid part in the swollen film
rho: the density of the solid part in the swollen film
"""
mutable struct FragmentBasedConstantTrhoDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::Float64
    rho::Float64
    A::Float64
    kfs::Array{W,1}
    krevs::Array{W,1}
    kfsnondiff::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{Int64,2}
    mu::W
    diffusivity::Array{W,1}
    jacobian::Array{W,2}
    sensitivity::Bool
    alternativepformat::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end

function FragmentBasedConstantTrhoDomain(; phase::Z, initialconds::Dict{X,E}, constantspecies::Array{X4,1}=Array{String,1}(),
    sparse::Bool=false, sensitivity::Bool=false) where {X,E,X1,E1,X3,X4,Z<:AbstractPhase}

    T = 0.0
    rho = 0.0
    mass = 0.0
    P = 1.0e8
    A = 0.0

    fragmentnames = getfield.(getphasespecies(phase), :name)

    y0 = zeros(length(fragmentnames) + 1) #track mass

    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "A"
            A = val
        elseif key == "rho"
            rho = val
        elseif key == "mass"
            mass = val
            y0[end] = val
        else
            ind = findfirst(isequal(key), fragmentnames)
            @assert typeof(ind) <: Integer "$key not found in fragment list: $fragmentnames"
            y0[ind] = val
        end
    end
    @assert T != 0.0
    @assert rho != 0.0
    @assert mass != 0.0
    @assert A != 0.0

    ns = y0[1:end-1]
    N = sum(ns)
    V = mass / rho

    if length(constantspecies) > 0
        for spc in constantspecies
            @assert spc in fragmentnames "$spc is not in the fragment list: $fragmentnames"
        end
        constspcinds = [findfirst(isequal(k), fragmentnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase, T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if phase.diffusionlimited
        diffs = [x(T=T, mu=mu, P=P) for x in getfield.(getphasespecies(phase), :diffusion)]
    else
        diffs = Array{Float64,1}()
    end

    C = N/V
    dGrxns = -(phase.stoichmatrix*Gs)
    kfs,krevs = getkfkrevs(phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    kfsnondiff = getkfs(phase,T,P,C,ns,V,0.0,dGrxns,0.0)
    p = vcat(Gs,kfsnondiff)
    if sparse
        jacobian = zeros(typeof(T), length(getphasespecies(phase)), length(getphasespecies(phase)))
    else
        jacobian = zeros(typeof(T), length(getphasespecies(phase)), length(getphasespecies(phase)))
    end
    rxnarray = getreactionindices(phase)
    return FragmentBasedConstantTrhoDomain(phase, [1, length(fragmentnames), length(fragmentnames) + 1], [1, length(phase.species) + length(phase.reactions)], constspcinds,
        T, rho, A, kfs, krevs, kfsnondiff, efficiencyinds, Gs, rxnarray, mu, diffs, jacobian, sensitivity, false, MVector(false), MVector(0.0), p, Dict("mass" => length(fragmentnames) + 1)), y0, p
end

export FragmentBasedConstantTrhoDomain

mutable struct ConstantTLiqFilmDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::W
    epsilon::W
    phi::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    kfsnondiff::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{Int64,2}
    mu::W
    diffusivity::Array{W,1}
    jacobian::Array{W,2}
    sensitivity::Bool
    alternativepformat::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end

function ConstantTLiqFilmDomain(; phase::Z, initialconds::Dict{X,E}, constantspecies::Array{X2,1}=Array{String,1}(),
    sparse=false, sensitivity=false) where {E,X,X2,Z<:AbstractPhase,Q<:AbstractInterface,W<:Real}
    #set conditions and initialconditions
    T = 0.0
    V = 0.0
    P = 1.0e8
    phi = 0.0
    epsilon = 0.0
    y0 = zeros(length(phase.species) + 1)
    spnames = [x.name for x in phase.species]
    for (key, val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            throw(error("ConstantTLiqFilmDomain cannot specify P"))
        elseif key == "V"
            V = val
            y0[end] = val
        elseif key == "epsilon"
            epsilon = val
        else
            ind = findfirst(isequal(key), spnames)
            @assert typeof(ind) <: Integer "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end
    @assert T != 0.0
    @assert V != 0.0
    @assert epsilon != 0.0

    ns = y0[1:end-1]
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species, :name)
        constspcinds = [findfirst(isequal(k), spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics) <: AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase, T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if phase.diffusionlimited
        diffs = [x(T=T, mu=mu, P=P) for x in getfield.(phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    P = 1.0e8  #essentiallly assuming this is a liquid
    C = N / V
    kfs, krevs = getkfkrevs(phase, T, P, C, N, ns, Gs, diffs, V, phi)
    kfsnondiff = getkfs(phase, T, P, C, ns, V, phi)
    p = vcat(Gs, kfsnondiff)
    if sparse
        jacobian = zeros(typeof(T), length(phase.species), length(phase.species))
    else
        jacobian = zeros(typeof(T), length(phase.species), length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTLiqFilmDomain(phase, [1, length(phase.species), length(phase.species) + 1], [1, length(phase.species) + length(phase.reactions)], constspcinds,
        T, epsilon, phi, kfs, krevs, kfsnondiff, efficiencyinds, Gs, rxnarray, mu, diffs, jacobian, sensitivity, false, MVector(false), MVector(0.0), p, Dict{String,Int64}(["V" => length(phase.species) + 1])), y0, p
end

export ConstantTLiqFilmDomain

@inline function calcthermo(d::ConstantTPDomain{W,Y}, y::J, t::Q, p::W3=SciMLBase.NullParameters()) where {W3<:SciMLBase.NullParameters,W<:IdealGas,Y<:Integer,J<:Array{Float64,1},Q} #no parameter input
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P * V / (R * d.T)
    cs = ns ./ V
    C = N / V
    for ind in d.efficiencyinds #efficiency related rates may have changed
        d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0)
    end
    return ns, cs, d.T, d.P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:Array{Float64,1},Q<:Float64} #uses parameter input
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P * V / (R * d.T)
    cs = ns ./ V
    C = N / V
    if !d.alternativepformat
        @views kfps = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        @views nothermochg = d.Gs == p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    else
        kfps = d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        nothermochg = d.Gs == d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    end
    @views nokfchg = count(d.kfs .!= kfps) <= length(d.efficiencyinds) && all(kfps[d.efficiencyinds] .== 1.0)
    if nothermochg && nokfchg
        for ind in d.efficiencyinds #efficiency related rates may have changed
            d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0;f=kfps[ind])
        end
        return ns, cs, d.T, d.P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
    elseif nothermochg
        d.kfs = kfps
        for ind in d.efficiencyinds #efficiency related rates may have changed
            d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0;f=kfps[ind])
        end
        return ns, cs, d.T, d.P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
    else #need to handle thermo changes
        d.kfs .= kfps
        if !d.alternativepformat
            d.Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        else
            d.Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        end
        krevs = getkfkrevs(d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0;kfs=d.kfs)[2]
        for ind in d.efficiencyinds #efficiency related rates may have changed
            d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0;f=kfps[ind])
        end
        return ns, cs, d.T, d.P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
    end
end

@inline function calcthermo(d::ConstantTPDomain{W,Y}, y::Array{W3,1}, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,W3<:ForwardDiff.Dual,Q} #Autodiff y
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P * V / (R * d.T)
    cs = ns ./ V
    C = N / V
    if !d.alternativepformat
        kfs = convert(typeof(y), p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    else
        kfs = convert(typeof(y), d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    end
    krevs = convert(typeof(y),getkfkrevs(d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;kfs=kfs)[2])
    for ind in d.efficiencyinds #efficiency related rates may have changed
        kfs[ind],krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;f=kfs[ind])
    end
    return ns, cs, d.T, d.P, V, C, N, d.mu, kfs, krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J,Q} #Autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P * V / (R * d.T)
    cs = ns ./ V
    C = N / V
    if !d.alternativepformat
        kfs = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    else
        kfs = d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    end
    krevs = getkfkrevs(d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;kfs=kfs)[2]
    for ind in d.efficiencyinds #efficiency related rates may have changed
        kfs[ind],krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;f=kfs[ind])
    end
    return ns, cs, d.T, d.P, V, C, N, d.mu, kfs, krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:Union{ReverseDiff.TrackedArray,Tracker.TrackedArray},Q} #Autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P * V / (R * d.T)
    cs = ns ./ V
    C = N / V
    kfs = similar(y, length(d.phase.reactions))
    krevs = similar(y, length(d.phase.reactions))
    if !d.alternativepformat
        kfs .= p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    else
        kfs .= d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    end
    krevs .= getkfkrevs(d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;kfs=kfs)[2]
    for ind in d.efficiencyinds #efficiency related rates may have changed
        kfs[ind],krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;f=kfs[ind])
    end
    return ns, cs, d.T, d.P, V, C, N, d.mu, kfs, krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Union{ReverseDiff.TrackedArray,Tracker.TrackedArray},W<:IdealGas,Y<:Integer,J,Q} #Tracker/reversediff
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P * V / (R * d.T)
    cs = ns ./ V
    C = N / V
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = [ind in d.efficiencyinds ? getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0)[1]*p[d.parameterindexes[1]-1+length(d.phase.species)+ind] : p[d.parameterindexes[1]-1+length(d.phase.species)+ind] for ind in 1:length(d.phase.reactions)]
    else
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = [ind in d.efficiencyinds ? getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0)[1]*d.p[length(d.phase.species)+ind]*p[d.parameterindexes[1]-1+length(d.phase.species)+ind] : p[d.parameterindexes[1]-1+length(d.phase.species)+ind] for ind in 1:length(d.phase.reactions)]
    end
    krevs = getkfkrevs(d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;kfs=kfs)[2]
    return ns,cs,d.T,d.P,V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Union{ReverseDiff.TrackedArray,Tracker.TrackedArray},W<:IdealGas,Y<:Integer,J<:Union{ReverseDiff.TrackedArray,Tracker.TrackedArray},Q} #Tracker/reversediff
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P * V / (R * d.T)
    cs = ns ./ V
    C = N / V
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = [ind in d.efficiencyinds ? getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0)[1]*p[d.parameterindexes[1]-1+length(d.phase.species)+ind] : p[d.parameterindexes[1]-1+length(d.phase.species)+ind] for ind in 1:length(d.phase.reactions)]
    else
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = [ind in d.efficiencyinds ? getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0)[1]*d.p[length(d.phase.species)+ind]*p[d.parameterindexes[1]-1+length(d.phase.species)+ind] : p[d.parameterindexes[1]-1+length(d.phase.species)+ind] for ind in 1:length(d.phase.reactions)]
    end
    krevs = getkfkrevs(d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;kfs=kfs)[2]
    return ns,cs,d.T,d.P,V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:SciMLBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P * d.V / (R * T)
    cs = ns ./ d.V
    C = N / d.V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    Cvave = 0.0
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Us = (hdivRT .- 1.0) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,d.V,0.0,0.0)
    return ns,cs,T,P,d.V,C,N,0.0,kfs,krevs,Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ConstantVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P * d.V / (R * T)
    cs = ns ./ d.V
    C = N / d.V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    Cvave = 0.0
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @views @fastmath hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Us = (hdivRT .- 1.0) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,d.V,0.0,0.0)
    return @views @fastmath ns,cs,T,P,d.V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ConstantVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P * d.V / (R * T)
    cs = ns ./ d.V
    C = N / d.V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    Cvave = 0.0
    cpdivR, hdivRT1, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Us = (hdivRT .- 1.0) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,d.V,0.0,0.0)
    return @views @fastmath ns,cs,T,P,d.V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ConstantPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:SciMLBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = d.P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Hs = hdivRT .* (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,d.P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return ns,cs,T,d.P,V,C,N,0.0,kfs,krevs,Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ConstantPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = d.P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Hs = hdivRT .* (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,d.P,C,N,ns,Gs,diffs,V,0.0,0.0)
    if p != SciMLBase.NullParameters()
        return @views @fastmath ns, cs, T, d.P, V, C, N, 0.0, kfs .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)], krevs .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)], Hs, Array{Float64,1}(), Gs, diffs, Cvave, cpdivR, 0.0
    else
        return ns, cs, T, d.P, V, C, N, 0.0, kfs, krevs, Hs, Array{Float64,1}(), Gs, diffs, Cvave, cpdivR, 0.0
    end
end

@inline function calcthermo(d::ConstantPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = d.P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR, hdivRT1, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Hs = hdivRT .* (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,d.P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return @views @fastmath ns,cs,T,d.P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:SciMLBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    V = d.V(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Us = (hdivRT .- 1.0) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return ns,cs,T,P,V,C,N,0.0,kfs,krevs,Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    V = d.V(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Us = (hdivRT .- 1.0) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    V = d.V(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR, hdivRT1, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Us = (hdivRT .- 1.0) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:SciMLBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Hs = hdivRT .* (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return ns,cs,T,P,V,C,N,0.0,kfs,krevs,Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Hs = hdivRT .* (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end
@inline function calcthermo(d::ParametrizedPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR, hdivRT1, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Hs = hdivRT .* (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedTConstantVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:SciMLBase.NullParameters,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q}
    V = d.V
    T = d.T(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ V
    C = N / V
    P = 1.0e8 #liquid phase
    Gs = zeros(length(d.phase.species))
    mu = d.phase.solvent.mu(T)
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=mu, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,d.phi,d.d)
    return ns,cs,T,P,V,C,N,mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ParametrizedTConstantVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q}
    V = d.V
    T = d.T(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ V
    C = N / V
    P = 1.0e8 #liquid phase
    Gs = zeros(length(d.phase.species))
    mu = d.phase.solvent.mu(T)
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=mu, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,d.phi,d.d)
    return @views @fastmath ns,cs,T,P,V,C,N,mu,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),d.phi
end

@inline function calcthermo(d::ParametrizedTConstantVDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q}
    V = d.V
    T = d.T(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ V
    C = N / V
    P = 1.0e8 #liquid phase
    Gs = zeros(length(d.phase.species))
    mu = d.phase.solvent.mu(T)
    cpdivR, hdivRT1, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=mu, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,d.phi,d.d)
    return @views @fastmath ns,cs,T,P,V,C,N,mu,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),d.phi
end

@inline function calcthermo(d::ParametrizedTPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:SciMLBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    T = d.T(t)
    @assert T < 10000.0
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return ns,cs,T,P,V,C,N,0.0,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ParametrizedTPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    T = d.T(t)
    @assert T < 10000.0
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR, hdivRT, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ParametrizedTPDomain{W,Y}, y::J, t::Q, p::W2=SciMLBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    T = d.T(t)
    @assert T < 10000.0
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = P * V / (R * T)
    cs = ns ./ V
    C = N / V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR, hdivRT1, sdivR = calcHSCpdless(d.phase.vecthermo, T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)] ./ (R * T)
    @fastmath Gs = (hdivRT .- sdivR) * (R * T)
    @fastmath Cvave = dot(cpdivR, ns)
    @fastmath Cvave *= R / N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = [x(T=T, mu=0.0, P=P) for x in getfield.(d.phase.species, :diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTVDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2<:SciMLBase.NullParameters,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ d.V
    C = N / d.V
    P = 1.0e8
    return ns, cs, d.T, P, d.V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
end

@inline function calcthermo(d::ConstantTVDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2<:Array{Float64,1},W<:IdealDiluteSolution,Y<:Integer,J<:Array{Float64,1},Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ d.V
    C = N / d.V
    P = 1.0e8
    if !d.alternativepformat
        @views nothermochg = d.Gs == p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        @views nokfchg = d.kfsnondiff == p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        if nothermochg && nokfchg
            return ns, cs, d.T, P, d.V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        elseif nothermochg
            d.kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,d.phi,d.d;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        else
            d.kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,d.phi,d.d;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        end
    else
        @views nothermochg = d.Gs == d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        @views nokfchg = d.kfsnondiff == d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        if nothermochg && nokfchg
            return ns, cs, d.T, P, d.V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        elseif nothermochg
            d.kfsnondiff .= d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,d.phi,d.d;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        else
            d.kfsnondiff .= d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs .= d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,d.phi,d.d;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        end
    end

end

@inline function calcthermo(d::ConstantTVDomain{W,Y}, y::Array{W2,1}, t::Q, p::Q2=SciMLBase.NullParameters()) where {W2<:ForwardDiff.Dual,Q2,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q} #autodiff y
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ d.V
    C = N / d.V
    P = 1.0e8
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = convert(typeof(y), p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    else
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = convert(typeof(y), d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    end
    kfs,krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,d.V,d.phi,d.d;kfs=kfsnondiff)
    return ns,cs,d.T,P,d.V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
end

@inline function calcthermo(d::ConstantTVDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q} #autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ d.V
    C = N / d.V
    P = 1.0e8
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    else
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    end
    kfs,krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,d.V,d.phi,d.d;kfs=kfsnondiff)
    return ns,cs,d.T,P,d.V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
end

@inline function calcthermo(d::ConstantTAPhiDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2<:SciMLBase.NullParameters,W<:IdealSurface,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ d.A
    C = N / d.A
    P = 0.0
    return ns, cs, d.T, P, d.A, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), d.Gs, Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
end

@inline function calcthermo(d::ConstantTAPhiDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2<:Array{Float64,1},W<:IdealSurface,Y<:Integer,J<:Array{Float64,1},Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ d.A
    C = N / d.A
    P = 0.0
    if !d.alternativepformat
        @views nothermochg = d.Gs == p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        if nothermochg
            return ns, cs, d.T, P, d.A, C, N, d.mu, p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)], d.krevs, Array{Float64,1}(), Array{Float64,1}(), d.Gs, Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        else
            d.kfs = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,d.phi,d.d;kfs=d.kfs)[2]
            return ns,cs,d.T,P,d.A,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),d.Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        end
    else
        @views nothermochg = d.Gs == d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        if nothermochg
            return ns, cs, d.T, P, d.A, C, N, d.mu, d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)], d.krevs, Array{Float64,1}(), Array{Float64,1}(), d.Gs, 0.0, Array{Float64,1}(), d.phi
        else
            d.kfs = d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,d.phi,d.d;kfs=d.kfs)[2]
            return ns,cs,d.T,P,d.A,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),d.Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        end
    end
end

@inline function calcthermo(d::ConstantTAPhiDomain{W,Y}, y::Array{W2,1}, t::Q, p::Q2=SciMLBase.NullParameters()) where {W2<:ForwardDiff.Dual,Q2,W<:IdealSurface,Y<:Integer,J<:AbstractArray,Q} #autodiff y
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ d.A
    C = N / d.A
    P = 0.0
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = convert(typeof(y), p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    else
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = convert(typeof(y), d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    end
    krevs = convert(typeof(y),getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,d.A,d.phi,d.d;kfs=kfs)[2])
    return ns,cs,d.T,P,d.A,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
end

@inline function calcthermo(d::ConstantTAPhiDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2,W<:IdealSurface,Y<:Integer,J<:AbstractArray,Q} #autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns ./ d.A
    C = N / d.A
    P = 0.0
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    else
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    end
    krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,d.A,d.phi,d.d;kfs=kfs)[2]
    return ns,cs,d.T,P,d.A,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
end

function calcthermo(d::FragmentBasedConstantTrhoDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2<:SciMLBase.NullParameters,W<:FragmentBasedIdealFilm,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    mass = y[d.indexes[3]]
    V = mass / d.rho
    cs = ns ./ V
    C = N / V
    P = 1.0e8
    return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
end

function calcthermo(d::FragmentBasedConstantTrhoDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2<:Array{Float64,1},W<:FragmentBasedIdealFilm,Y<:Integer,J<:Array{Float64,1},Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    mass = y[d.indexes[3]]
    V = mass / d.rho
    cs = ns ./ V
    C = N / V
    P = 1.0e8
    if !d.alternativepformat
        @views nothermochg = d.Gs == p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        @views nokfchg = d.kfsnondiff == p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        if nothermochg && nokfchg
            return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
        elseif nothermochg
            d.kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        else
            d.kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        end
    else
        @views nothermochg = d.Gs == d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        @views nokfchg = d.kfsnondiff == d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        if nothermochg && nokfchg
            return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
        elseif nothermochg
            d.kfsnondiff .= d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        else
            d.kfsnondiff .= d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs .= d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,V,0.0,0.0;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        end
    end

end

function calcthermo(d::FragmentBasedConstantTrhoDomain{W,Y}, y::Array{W2,1}, t::Q, p::Q2=SciMLBase.NullParameters()) where {W2<:ForwardDiff.Dual,Q2,W<:FragmentBasedIdealFilm,Y<:Integer,J<:AbstractArray,Q} #autodiff y
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    mass = y[d.indexes[3]]
    V = mass / d.rho
    cs = ns ./ V
    C = N / V
    P = 1.0e8
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = convert(typeof(y), p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    else
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = convert(typeof(y), d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    end
    kfs,krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;kfs=kfsnondiff)
    return ns,cs,d.T,P,V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

function calcthermo(d::FragmentBasedConstantTrhoDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2,W<:FragmentBasedIdealFilm,Y<:Integer,J<:AbstractArray,Q} #autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    mass = y[d.indexes[3]]
    V = mass / d.rho
    cs = ns ./ V
    C = N / V
    P = 1.0e8
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    else
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    end
    kfs, krevs = getkfkrevs(d.phase, d.T, P, C, N, ns, Gs, d.diffusivity, V, 0.0; kfs=kfsnondiff)
    return ns, cs, d.T, P, V, C, N, d.mu, kfs, krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), 0.0
end

@inline function calcthermo(d::ConstantTLiqFilmDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2<:Array{Float64,1},W<:IdealDiluteSolution,Y<:Integer,J<:Array{Float64,1},Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = sum(ns)
    cs = ns ./ V
    C = N / V
    P = 1.0e8
    if !d.alternativepformat
        @views nothermochg = d.Gs == p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        @views nokfchg = d.kfsnondiff == p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        if nothermochg && nokfchg
            return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        elseif nothermochg
            d.kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.kfs, d.krevs = getkfkrevs(d.phase, d.T, P, C, N, ns, d.Gs, d.diffusivity, V, d.phi; kfs=d.kfsnondiff)
            return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        else
            d.kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.kfs, d.krevs = getkfkrevs(d.phase, d.T, P, C, N, ns, d.Gs, d.diffusivity, V, d.phi; kfs=d.kfsnondiff)
            return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        end
    else
        @views nothermochg = d.Gs == d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        @views nokfchg = d.kfsnondiff == d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        if nothermochg && nokfchg
            return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        elseif nothermochg
            d.kfsnondiff .= d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.kfs, d.krevs = getkfkrevs(d.phase, d.T, P, C, N, ns, d.Gs, d.diffusivity, V, d.phi; kfs=d.kfsnondiff)
            return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        else
            d.kfsnondiff .= d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs .= d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.kfs, d.krevs = getkfkrevs(d.phase, d.T, P, C, N, ns, d.Gs, d.diffusivity, V, d.phi; kfs=d.kfsnondiff)
            return ns, cs, d.T, P, V, C, N, d.mu, d.kfs, d.krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
        end
    end

end

@inline function calcthermo(d::ConstantTLiqFilmDomain{W,Y}, y::Array{W2,1}, t::Q, p::Q2=SciMLBase.NullParameters()) where {W2<:ForwardDiff.Dual,Q2,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q} #autodiff y
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = sum(ns)
    cs = ns ./ V
    C = N / V
    P = 1.0e8
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = convert(typeof(y), p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    else
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = convert(typeof(y), d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    end
    kfs, krevs = getkfkrevs(d.phase, d.T, P, C, N, ns, Gs, d.diffusivity, V, d.phi; kfs=kfsnondiff)
    return ns, cs, d.T, P, V, C, N, d.mu, kfs, krevs, Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), 0.0, Array{Float64,1}(), d.phi
end

@inline function calcthermo(d::ConstantTLiqFilmDomain{W,Y}, y::J, t::Q, p::Q2=SciMLBase.NullParameters()) where {Q2,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q} #autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = sum(ns)
    cs = ns ./ V
    C = N / V
    P = 1.0e8
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    else
        Gs = d.p[1:length(d.phase.species)] .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = d.p[length(d.phase.species)+1:end] .* p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    end
    kfs,krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,V,0.0,0.0;kfs=kfsnondiff)
    return ns,cs,d.T,P,V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end
export calcthermo

@inline function calcdomainderivatives!(d::Q, dydt::Z7, interfaces::Z12; t::Z10, T::Z4, P::Z9, Us::Array{Z,1}, Hs::Array{Z11,1}, V::Z2, C::Z3, ns::Z5, N::Z6, Cvave::Z8) where {Q<:AbstractDomain,Z12,Z11,Z10,Z9,Z8<:Real,Z7,W<:IdealGas,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5}
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter, Inlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y .* inter.F(t)
        elseif isa(inter, Outlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.F(t) .* ns ./ N
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && d == inter.domain
            kLAs = map.(inter.kLAs, T)
            kHs = map.(inter.kHs, T)
            evap = kLAs .* ns
            cond = kLAs .* inter.molefractions .* inter.P ./ kHs * V

            dydt[d.indexes[1]:d.indexes[2]] .-= (evap .- cond)
        elseif isa(inter, VolumetricFlowRateInlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.Vin(t) * inter.cs
        elseif isa(inter, VolumetricFlowRateOutlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.Vout(t) * ns / V
        end
    end
end

@inline function calcdomainderivatives!(d::Q, dydt::Z7, interfaces::Z12; t::Z10, T::Z4, P::Z9, Us::Array{Z,1}, Hs::Array{Z11,1}, V::Z2, C::Z3, ns::Z5, N::Z6, Cvave::Z8) where {Q<:ConstantTPDomain,Z12,Z11,Z10,Z9,Z8<:Real,Z7,W<:IdealGas,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5}
    @views @fastmath @inbounds dydt[d.indexes[3]] = sum(dydt[d.indexes[1]:d.indexes[2]]) * R * T / P
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter, Inlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y .* inter.F(t)
            dydt[d.indexes[3]] += inter.F(t) * R * T / P
        elseif isa(inter, Outlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.F(t) .* ns ./ N
            dydt[d.indexes[3]] -= inter.F(t) * R * T / P
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && d == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs
            cond = kLAs .* inter.V .* cs * R * T ./ kHs
            net_evap = evap .- cond
            dydt[d.indexes[1]:d.indexes[2]] .+= net_evap
            dydt[d.indexes[3]] += sum(net_evap) * R * T / P
        elseif isa(inter, VolumetricFlowRateInlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.Vin(t) * inter.cs
            dydt[d.indexes[3]] += inter.Vin(t)
        elseif isa(inter, VolumetricFlowRateOutlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.Vout(t) * ns / V
            dydt[d.indexes[3]] -= inter.Vout(t)
        end
    end
    for inter in interfaces
        if isa(inter, VolumeMaintainingOutlet) && d == inter.domain #VolumeMaintainingOutlet has to be evaluated after dVdt has been modified by everything else
            @inbounds dVdt = dydt[d.indexes[3]]
            @inbounds flow = P * dVdt / (R * T)
            @views @inbounds dydt[d.indexes[1]:d.indexes[2]] .-= flow * ns / N
            @inbounds dydt[d.indexes[3]] -= dVdt
        end
    end
end

@inline function calcdomainderivatives!(d::ConstantVDomain{W,Y}, dydt::K, interfaces::Z12; t::Z10, T::Z4, P::Z9, Us::Z, Hs::Z11, V::Z2, C::Z3, ns::Z5, N::Z6, Cvave::Z7) where {Z12,Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5}
    @views @fastmath @inbounds dydt[d.indexes[3]] = -dot(Us, dydt[d.indexes[1]:d.indexes[2]]) / (N * Cvave) #divide by V to cancel ωV to ω
    @views @fastmath @inbounds dydt[d.indexes[4]] = sum(dydt[d.indexes[1]:d.indexes[2]]) * R * T / V + P / T * dydt[d.indexes[3]]
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter, Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y .* flow
            dTdt = flow * (inter.H - dot(Us, ns) / N) / (N * Cvave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow * R * T / V + P / T * dTdt
        elseif isa(inter, Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow .* ns ./ N
            dTdt = (P * V / N * flow) / (N * Cvave)
            dydt[d.indexes[3]] -= dTdt
            dydt[d.indexes[4]] -= flow * R * T / V + P / T * dTdt
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && d == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs
            cond = kLAs .* inter.V .* cs * R * T ./ kHs
            dydt[d.indexes[1]:d.indexes[2]] .+= (evap .- cond)

            flow = sum(evap)
            dTdt = flow * (inter.H - dot(Us, ns) / N) / (N * Cvave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow * R * T / V + P / T * dTdt

            flow = sum(cond)
            dTdt = (P * V / N * flow) / (N * Cvave)
            dydt[d.indexes[3]] -= dTdt
            dydt[d.indexes[4]] -= flow * R * T / V + P / T * dTdt
        elseif isa(inter, VolumetricFlowRateInlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.Vin(t) * inter.cs
            dTdt = inter.Vin(t) * (inter.Hpervolume - dot(Us, ns) / V) / (N * Cvave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += inter.Vin(t) * P / V + P / T * dTdt
        elseif isa(inter, VolumetricFlowRateOutlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.Vout(t) * ns / V
            dTdt = (P * inter.Vout(t)) / (N * Cvave)
            dydt[d.indexes[3]] -= dTdt
            dydt[d.indexes[4]] -= inter.Vout(t) * P / V + P / T * dTdt
        end
    end
end

@inline function calcdomainderivatives!(d::ConstantPDomain{W,Y}, dydt::K, interfaces::Z12; t::Z10, T::Z4, P::Z9, Us::Z, Hs::Z11, V::Z2, C::Z3, ns::Z5, N::Z6, Cvave::Z7) where {Z12,Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5}
    @fastmath Cpave = Cvave + R
    @views @fastmath @inbounds dydt[d.indexes[3]] = -dot(Hs, dydt[d.indexes[1]:d.indexes[2]]) / (N * Cpave) #divide by V to cancel ωV to ω
    @views @fastmath @inbounds dydt[d.indexes[4]] = sum(dydt[d.indexes[1]:d.indexes[2]]) * R * T / P + dydt[d.indexes[3]] * V / T
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter, Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y .* flow
            dTdt = flow * (inter.H - dot(Hs, ns) / N) / (N * Cpave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow * R * T / P + dTdt * V / T
        elseif isa(inter, Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow .* ns ./ N
            dydt[d.indexes[4]] -= flow * R * T / P
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && d == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs
            cond = kLAs .* inter.V .* cs * R * T ./ kHs
            dydt[d.indexes[1]:d.indexes[2]] .+= (evap .- cond)

            flow = sum(evap)
            dTdt = flow * (inter.H - dot(Hs, ns) / N) / (N * Cpave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow * R * T / P + dTdt * V / T

            flow = sum(cond)
            dydt[d.indexes[4]] -= flow * R * T / P
        elseif isa(inter, VolumetricFlowRateInlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.Vin(t) * inter.cs
            dTdt = inter.Vin(t) * (inter.Hpervolume - dot(Hs, ns) / V) / (N * Cpave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += inter.Vin(t)
        elseif isa(inter, VolumetricFlowRateOutlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.Vout(t) * ns / V
            dydt[d.indexes[4]] -= inter.Vout(t)
        end
    end
    for inter in interfaces
        if isa(inter, VolumeMaintainingOutlet) && d == inter.domain #VolumeMaintainingOutlet has to be evaluated after dVdt has been modified by everything else
            @inbounds dVdt = dydt[d.indexes[4]]
            @inbounds flow = P * dVdt / (R * T)
            @views @inbounds dydt[d.indexes[1]:d.indexes[2]] .-= flow * ns / N
            @inbounds dydt[d.indexes[4]] -= dVdt
        end
    end
end

@inline function calcdomainderivatives!(d::ParametrizedTPDomain{W,Y}, dydt::K, interfaces::Z12; t::Z10, T::Z4, P::Z9, Us::Z, Hs::Z11, V::Z2, C::Z3, ns::Z5, N::Z6, Cvave::Z7) where {Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5,Z12}
    @views @fastmath @inbounds dydt[d.indexes[3]] = sum(dydt[d.indexes[1]:d.indexes[2]]) * R * T / P + Calculus.derivative(d.T, t) * V / T - Calculus.derivative(d.P, t) * V / P
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter, Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y .* flow
            dydt[d.indexes[3]] += flow * R * T / P
        elseif isa(inter, Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow * ns ./ N
            dydt[d.indexes[3]] -= flow * R * T / P
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && d == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs
            cond = kLAs .* inter.V .* cs * R * T ./ kHs
            dydt[d.indexes[1]:d.indexes[2]] .+= (evap .- cond)

            flow = sum(evap)
            dydt[d.indexes[3]] += flow * R * T / P

            flow = sum(cond)
            dydt[d.indexes[3]] -= flow * R * T / P
        elseif isa(inter, VolumetricFlowRateInlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.Vin(t) * inter.cs
            dydt[d.indexes[3]] += inter.Vin(t)
        elseif isa(inter, VolumetricFlowRateOutlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.Vout(t) * ns / V
            dydt[d.indexes[3]] -= inter.Vout(t)
        end
    end
    for inter in interfaces
        if isa(inter, VolumeMaintainingOutlet) && d == inter.domain #VolumeMaintainingOutlet has to be evaluated after dVdt has been modified by everything else
            @inbounds dVdt = dydt[d.indexes[3]]
            @inbounds flow = P * dVdt / (R * T)
            @views @inbounds dydt[d.indexes[1]:d.indexes[2]] .-= flow * ns / N
            @inbounds dydt[d.indexes[3]] -= dVdt
        end
    end
end

@inline function calcdomainderivatives!(d::ParametrizedVDomain{W,Y}, dydt::K, interfaces::Z12; t::Z10, T::Z4, P::Z9, Us::Z, Hs::Z11, V::Z2, C::Z3, ns::Z5, N::Z6, Cvave::Z7) where {Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5,Z12}
    dVdt = Calculus.derivative(d.V, t)
    @views @fastmath @inbounds dydt[d.indexes[3]] = (-dot(Us, dydt[d.indexes[1]:d.indexes[2]]) - P * dVdt) / (N * Cvave) #divide by V to cancel ωV to ω
    @views @fastmath @inbounds dydt[d.indexes[4]] = sum(dydt[d.indexes[1]:d.indexes[2]]) * R * T / V + dydt[d.indexes[3]] * P / T - P / V * dVdt
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter, Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y .* flow
            dTdt = flow * (inter.H - dot(Us, ns) / N) / (N * Cvave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow * R * T / V + dTdt * P / T
        elseif isa(inter, Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow * ns ./ N
            dTdt = (P * V / N * flow) / (N * Cvave)
            dydt[d.indexes[3]] -= dTdt
            dydt[d.indexes[4]] -= flow * R * T / V + dTdt * P / T
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && d == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs
            cond = kLAs .* inter.V .* cs * R * T ./ kHs
            dydt[d.indexes[1]:d.indexes[2]] .+= (evap .- cond)

            flow = sum(evap)
            dTdt = flow * (inter.H - dot(Us, ns) / N) / (N * Cvave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow * R * T / V + dTdt * P / T

            flow = sum(cond)
            dTdt = (P * V / N * flow) / (N * Cvave)
            dydt[d.indexes[3]] -= dTdt
            dydt[d.indexes[4]] -= flow * R * T / V + dTdt * P / T
        elseif isa(inter, VolumetricFlowRateInlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.Vin(t) * inter.cs
            dTdt = inter.Vin(t) * (inter.Hpervolume - dot(Us, ns) / V) / (N * Cvave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += inter.Vin(t) * P / V + dTdt * P / T
        elseif isa(inter, VolumetricFlowRateOutlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.Vout(t) * ns / V
            dTdt = (P * inter.Vout(t)) / (N * Cvave)
            dydt[d.indexes[3]] -= dTdt
            dydt[d.indexes[4]] -= inter.Vout(t) * P / V + dTdt * P / T
        end
    end
end

@inline function calcdomainderivatives!(d::ParametrizedPDomain{W,Y}, dydt::K, interfaces::Z12; t::Z10, T::Z4, P::Z9, Us::Z, Hs::Z11, V::Z2, C::Z3, ns::Z5, N::Z6, Cvave::Z7) where {Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5,Z12}
    @fastmath Cpave = Cvave + R
    dPdt = Calculus.derivative(d.P, t)
    @views @fastmath @inbounds dydt[d.indexes[3]] = (-dot(Hs, dydt[d.indexes[1]:d.indexes[2]]) + V * dPdt) / (N * Cpave) #divide by V to cancel ωV to ω
    @views @fastmath @inbounds dydt[d.indexes[4]] = sum(dydt[d.indexes[1]:d.indexes[2]]) * R * T / P + dydt[d.indexes[3]] * V / T - dPdt * V / P
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter, Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y .* flow
            dTdt = flow * (inter.H - dot(Hs, ns) / N) / (N * Cpave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow * R * T / P + dTdt * V / T
        elseif isa(inter, Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow .* ns ./ N
            dydt[d.indexes[4]] -= flow * R * T / P
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && d == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs
            cond = kLAs .* inter.V .* cs * R * T ./ kHs
            dydt[d.indexes[1]:d.indexes[2]] .+= (evap .- cond)

            flow = sum(evap)
            dTdt = flow * (inter.H - dot(Hs, ns) / N) / (N * Cpave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow * R * T / P + dTdt * V / T

            flow = sum(cond)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow .* ns ./ N
            dydt[d.indexes[4]] -= flow * R * T / P
        elseif isa(inter, VolumetricFlowRateInlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.Vin(t) * inter.cs
            dTdt = inter.Vin(t) * (inter.Hpervolume - dot(Hs, ns) / V) / (N * Cpave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += inter.Vin(t)
        elseif isa(inter, VolumetricFlowRateOutlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.Vout(t) * ns / V
            dydt[d.indexes[4]] -= inter.Vout(t)
        end
    end
    for inter in interfaces
        if isa(inter, VolumeMaintainingOutlet) && d == inter.domain #VolumeMaintainingOutlet has to be evaluated after dVdt has been modified by everything else
            @inbounds dVdt = dydt[d.indexes[4]]
            @inbounds flow = P * dVdt / (R * T)
            @views @inbounds dydt[d.indexes[1]:d.indexes[2]] .-= flow * ns / N
            @inbounds dydt[d.indexes[4]] -= dVdt
        end
    end
end
export calcdomainderivatives!

@inline function jacobianyefficiencyderiv!(jac::S, domain::Union{ConstantTPDomain,ConstantPDomain,ParametrizedTPDomain,ParametrizedPDomain}, kinetics::AbstractFalloffRate, efficiencies::Dict{Int64,Float64}, rxnarray::Array{Int64,2}, rxnind::Int64, Vind::Int64, cs::Array{Float64,1}, T::Float64, V::Float64, C::Float64, Ceff::Float64, Kc::Float64) where {S<:AbstractArray}
    dkdCeff = _calcdkdCeff(kinetics, T, Ceff)
    for (spcind, efficiencyval) in efficiencies
        @fastmath dkdni = dkdCeff * (efficiencyval / V)
        _jacobianykderiv!(jac, spcind, dkdni, dkdni / Kc, rxnarray, rxnind, cs, V)
    end
    @fastmath dkdV = dkdCeff * (C - Ceff) / V
    _jacobianykderiv!(jac, Vind, dkdV, dkdV / Kc, rxnarray, rxnind, cs, V)
end

@inline function jacobianyefficiencyderiv!(jac::S, domain::Union{ConstantVDomain,ParametrizedVDomain,ConstantTVDomain,ParametrizedTConstantVDomain}, kinetics::AbstractFalloffRate, efficiencies::Dict{Int64,Float64}, rxnarray::AbstractArray, rxnind::Int64, Vind::Int64, cs::Array{Float64,1}, T::Float64, V::Float64, C::Float64, Ceff::Float64, Kc::Float64) where {S<:AbstractArray}
    dkdCeff = _calcdkdCeff(kinetics, T, Ceff)
    for (spcind, efficiencyval) in efficiencies
        @fastmath dkdni = dkdCeff * (efficiencyval / V)
        _jacobianykderiv!(jac, spcind, dkdni, dkdni / Kc, rxnarray, rxnind, cs, V)
    end
end

@inline function jacobianynsderiv!(jac::S, domain::Union{ConstantTPDomain,ParametrizedTPDomain}, rxnarray::Array{Int64,2}, efficiencyinds::Array{I,1}, cs::Array{Float64,1}, kfs::Array{Float64,1}, krevs::Array{Float64,1}, T::Float64, V::Float64, C::Float64) where {S<:AbstractArray,I<:Integer}
    @simd for rxnind = 1:size(rxnarray)[2]
        @inbounds _jacobianynswrtns!(jac, rxnarray, rxnind, cs, kfs[rxnind], krevs[rxnind])
        @inbounds _jacobianynswrtV!(jac, domain.indexes[3], rxnarray, rxnind, cs, kfs[rxnind], krevs[rxnind])
        if rxnind in efficiencyinds
            @inbounds @fastmath Kc = kfs[rxnind] / krevs[rxnind]
            @inbounds @fastmath Ceff = C + sum(cs[i] * val for (i, val) in domain.phase.reactions[rxnind].kinetics.efficiencies)
            @inbounds jacobianyefficiencyderiv!(jac, domain, domain.phase.reactions[rxnind].kinetics, domain.phase.reactions[rxnind].kinetics.efficiencies, rxnarray, rxnind, domain.indexes[3], cs, T, V, C, Ceff, Kc)
        end
    end
end

@inline function jacobianynsderiv!(jac::S, domain::Union{ConstantVDomain,ParametrizedVDomain,ConstantTVDomain,ConstantTAPhiDomain,ParametrizedTConstantVDomain}, rxnarray::Array{Int64,2}, efficiencyinds::Array{I,1}, cs::Array{Float64,1}, kfs::Array{Float64,1}, krevs::Array{Float64,1}, T::Float64, V::Float64, C::Float64) where {S<:AbstractArray,I<:Integer}
    @simd for rxnind = 1:size(rxnarray)[2]
        @inbounds _jacobianynswrtns!(jac, rxnarray, rxnind, cs, kfs[rxnind], krevs[rxnind])
        if rxnind in efficiencyinds
            @inbounds @fastmath Kc = kfs[rxnind] / krevs[rxnind]
            @inbounds @fastmath Ceff = C + sum(cs[i] * val for (i, val) in domain.phase.reactions[rxnind].kinetics.efficiencies)
            @inbounds jacobianyefficiencyderiv!(jac, domain, domain.phase.reactions[rxnind].kinetics, domain.phase.reactions[rxnind].kinetics.efficiencies, rxnarray, rxnind, 0, cs, T, V, C, Ceff, Kc)
        end
    end
end

@inline function jacobianynsderiv!(jac::S, domain::Union{ConstantPDomain,ParametrizedPDomain}, rxnarray::Array{Int64,2}, efficiencyinds::Array{I,1}, cs::Array{Float64,1}, kfs::Array{Float64,1}, krevs::Array{Float64,1}, T::Float64, V::Float64, C::Float64) where {S<:AbstractArray,I<:Integer}
    @simd for rxnind = 1:size(rxnarray)[2]
        @inbounds _jacobianynswrtns!(jac, rxnarray, rxnind, cs, kfs[rxnind], krevs[rxnind])
        @inbounds _jacobianynswrtV!(jac, domain.indexes[4], rxnarray, rxnind, cs, kfs[rxnind], krevs[rxnind])
        if rxnind in efficiencyinds
            @inbounds @fastmath Kc = kfs[rxnind] / krevs[rxnind]
            @inbounds @fastmath Ceff = C + sum(cs[i] * val for (i, val) in domain.phase.reactions[rxnind].kinetics.efficiencies)
            @inbounds jacobianyefficiencyderiv!(jac, domain, domain.phase.reactions[rxnind].kinetics, domain.phase.reactions[rxnind].kinetics.efficiencies, rxnarray, rxnind, domain.indexes[4], cs, T, V, C, Ceff, Kc)
        end
    end
end

@inline function jacobiany!(jac::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ConstantTPDomain}
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)
    jacobianynsderiv!(jac, domain, domain.rxnarray, domain.efficiencyinds, cs, kfs, krevs, T, V, C)

    # dV/dt = sum(dn/dt)*R*T/P
    # d/dni(dV/dt) = sum(d/dni(dn/dt))*R*T/P
    # d/dV(dV/dt) = sum(d/dV(dn/dt))*R*T/P
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @views @inbounds @fastmath jac[domain.indexes[3], i] = sum(jac[domain.indexes[1]:domain.indexes[2], i]) * R * T / P
    end
    @views @inbounds @fastmath jac[domain.indexes[3], domain.indexes[3]] = sum(jac[domain.indexes[1]:domain.indexes[2], domain.indexes[3]]) * R * T / P

    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind, :] .= 0.0
    end

    @simd for inter in interfaces

        if isa(inter, Inlet) && domain == inter.domain
            # dn/dt .+= inter.y.*inter.F(t)
            # dV/dt += inter.F(t)*R*T/P
            # d/dni(dn/dt) += 0
            # d/dV(dn/dt) += 0
            # d/dni(dV/dt) += 0
            # d/dV(dV/dt) += 0
            nothing

        elseif isa(inter, Outlet) && domain == inter.domain
            # dn/dt .-= inter.F(t).*ns./N = inter.F(t).*ns*(RT/PV)
            # dV/dt -= inter.F(t)*R*T/P
            # d/dni(dni/dt) -= inter.F(t)/N
            # d/dV(dn/dt) -= inter.F(t)*ns*(RT/PV)*(-1/V) = inter.F(t)*ns/N*(-1/V) = -inter.F(t)*ns/(N*V)
            # d/dni(dV/dt) -= 0
            # d/dV(dV/dt) -= 0

            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= flow / N
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[3]] .-= -flow .* ns / (N * V)
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && domain == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)

            # evaporation
            # dn/dt .+= kLAs.*inter.V.*inter.cs
            # dV/dt += sum(kLAs.*inter.V.*inter.cs)*R*T/P
            # d/dni(dn/dt) += 0
            # d/dV(dn/dt) += 0
            # d/dni(dV/dt) += 0
            # d/dV(dV/dt) += 0

            # condensation
            # dn/dt .-= kLAs.*inter.V.*cs*R*T./kHs
            # dV/dt -= sum(kLAs.*inter.V.*cs*R*T./kHs)*R*T/P
            # d/dni(dni/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]
            # d/dV(dn/dt) -= kLAs.*inter.V.*cs*(-1/V)*R*T./kHs
            # d/dni(dV/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]*R*T/P
            #              -= d/dni(dni/dt)*R*T/P
            # d/dV(dV/dt) -= sum(d/dV(dn/dt))*R*T/P
            @simd for i in domain.indexes[1]:domain.indexes[2]
                ddnidnidt = kLAs[i] * inter.V / V * R * T / kHs[i]
                @inbounds @fastmath jac[i, i] -= ddnidnidt
                @inbounds @fastmath jac[domain.indexes[3], i] -= ddnidnidt * R * T / P
            end
            ddVdndt = kLAs .* inter.V .* cs * (-1 / V) * R * T ./ kHs
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[3]] .-= ddVdndt
            @views @inbounds @fastmath jac[domain.indexes[3], domain.indexes[3]] -= sum(ddVdndt) * R * T / P
        elseif isa(inter, VolumetricFlowRateInlet) && domain == inter.domain
            # dn/dt += inter.Vin(t)*inter.cs
            # dV/dt += inter.Vin(t)
            # d/dn(dn/dt) += 0
            # d/dV(dn/dt) += 0
            # d/dn(dV/dt) += 0
            # d/dV(dV/dt) += 0
            nothing

        elseif isa(inter, VolumetricFlowRateOutlet) && domain == inter.domain
            # dn/dt .-= inter.Vout(t)*ns/V
            # dV/dt -= inter.Vout(t)
            # d/dni(dni/dt) -= inter.Vout(t)/V
            # d/dV(dn/dt) -= -inter.Vout(t)*ns/(V*V)
            # d/dni(dV/dt) -= 0
            # d/dV(dV/dt) -= 0

            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= inter.Vout(t) / V
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[3]] .-= -inter.Vout(t) * ns / (V * V)
        end
    end
end

@inline function jacobiany!(jac::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ConstantVDomain}
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)
    jacobianynsderiv!(jac, domain, domain.rxnarray, domain.efficiencyinds, cs, kfs, krevs, T, V, C)

    dydt = zeros(size(y))
    addreactionratecontributions!(dydt, domain.rxnarray, cs, kfs, krevs)
    dydt .*= V

    # dT/dt = -dot(Us,dn/dt)/(N*Cvave)
    # dP/dt = sum(dn/dt)*R*T/V + P/T*dT/dt
    # d/dni(dT/dt) = -dot(Us,d/dni(dn/dt))/(N*Cvave) + dT/dt/(N*Cvave) * d/dni(N*Cvave)
    #              = -dot(Us,d/dni(dn/dt))/(N*Cvave) + dT/dt/(Cvave) * d/dni(Cvave)
    #                Note: Cvave = dot(cpdivR,ns)*R/N-R
    #                Note: d/dni(Cvave) = cpdivR[i]*R/N
    # d/dni(dP/dt) = sum(d/dni(dn/dt))*R*T/V + P/T*d/dni(dT/dt)
    for i in domain.indexes[1]:domain.indexes[2]
        @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
        @views @inbounds @fastmath jac[domain.indexes[3], i] = -dot(Us, jac[domain.indexes[1]:domain.indexes[2], i]) / (N * Cvave) + dot(Us, dydt[domain.indexes[1]:domain.indexes[2]]) / (N * Cvave * Cvave) * dCvavedni
        @views @inbounds @fastmath jac[domain.indexes[4], i] = sum(jac[domain.indexes[1]:domain.indexes[2], i]) * R * T / V + P / T * jac[domain.indexes[3], i]
    end

    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind, :] .= 0.0
    end

    for inter in interfaces
        if isa(inter, Inlet) && domain == inter.domain
            # dn/dt .+= inter.y.*inter.F(t)
            # dTdt += inter.F(t)*(inter.H - dot(Us,ns)/N)/(N*Cvave)
            # dP/dt += inter.F(t)*R*T/V + P/T*dTdt
            # d/dni(dn/dt) += 0
            # d/dni(dTdt) += inter.F(t)*(-Us[i]/N)/(N*Cvave) - dTdt/(N*Cvave) * d/dni(N*Cvave)
            #             += inter.F(t)*(-Us[i]/N)/(N*Cvave) - dTdt/(Cvave) * d/dni(Cvave)
            #                Note: Cvave = dot(cpdivR,ns)*R/N-R
            #                Note: d/dni(Cvave) = cpdivR[i]*R/N
            # d/dni(dP/dt) += P/T*d/dni(dTdt)
            flow = inter.F(t)
            @fastmath dTdt = flow * (inter.H - dot(Us, ns) / N) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = flow * (-Us[i] / N) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += P / T * ddnidTdt
            end
        elseif isa(inter, Outlet) && domain == inter.domain
            # dn/dt .-= inter.F(t).*ns./N
            # dT/dt -= (P*V/N*inter.F(t))/(N*Cvave)
            # dP/dt -= inter.F(t)*R*T/V + P/T*dT/dt
            # d/dni(dni/dt) -= inter.F(t)/N
            # d/dni(dT/dt) -= dT/dt * (-1/(N*Cvave)) * d/dni(N*Cvave)
            #              -= dT/dt * (-1/(Cvave)) * d/dni(Cvave)
            #                 Note: Cvave = dot(cpdivR,ns)*R/N-R
            #                 Note: d/dni(Cvave) = cpdivR[i]*R/N
            #              -= -dT/dt * (dCvavedni/Cvave)
            # d/dni(dP/dt) -= P/T*d/dni(dT/dt)
            flow = inter.F(t)
            @fastmath dTdt = (P * V / N * flow) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= flow / N
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @fastmath ddnidTdt = -dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] -= ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] -= P / T * ddnidTdt
            end
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && domain == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs
            cond = kLAs .* inter.V .* cs * R * T ./ kHs

            # evaporation
            # dn/dt .+= kLAs.*inter.V.*inter.cs
            # dT/dt += sum(kLAs.*inter.V.*inter.cs)*(inter.H - dot(Us,ns)/N)/(N*Cvave)
            # dP/dt += sum(kLAs.*inter.V.*inter.cs)*R*T/V + P/T*dT/dt
            # d/dni(dni/dt) += 0
            # d/dni(dT/dt) += sum(kLAs.*inter.V.*inter.cs)*(-Us[i]/N)/(N*Cvave) - dT/dt * (-1/(N*Cvave)) * d/dni(N*Cvave)
            #              += sum(kLAs.*inter.V.*inter.cs)*(-Us[i]/N)/(N*Cvave) - dT/dt * (-1/(Cvave)) * d/dni(Cvave)
            #                 Note: Cvave = dot(cpdivR,ns)*R/N-R
            #                 Note: d/dni(Cvave) = cpdivR[i]*R/N
            #              += sum(kLAs.*inter.V.*inter.cs)*(-Us[i]/N)/(N*Cvave) - dT/dt * (dCvavedni/Cvave)
            # d/dni(dP/dt) += P/T*d/dni(dT/dt)
            flow = sum(evap)
            @fastmath dTdt = flow * (inter.H - dot(Us, ns) / N) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = flow * (-Us[i] / N) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += P / T * ddnidTdt
            end

            # condensation
            # dn/dt .-= kLAs.*inter.V.*cs*R*T./kHs
            # dT/dt -= (P*V/N*sum(kLAs.*inter.V.*cs*R*T./kHs))/(N*Cvave)
            # dP/dt -= sum(kLAs.*inter.V.*cs*R*T./kHs)*R*T/V + P/T*dT/dt
            # d/dni(dni/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]
            # d/dni(dT/dt) -= (P*V/N*kLAs[i]*inter.V/V*R*T/kHs[i])/(N*Cvave) + dT/dt * (-1/(N*Cvave)) * d/dni(N*Cvave)
            #              -= (P*V/N*kLAs[i]*inter.V/V*R*T/kHs[i])/(N*Cvave) + dT/dt * (-1/(Cvave)) * d/dni(Cvave)
            #                 Note: Cvave = dot(cpdivR,ns)*R/N-R
            #                 Note: d/dni(Cvave) = cpdivR[i]*R/N
            #              -= (P*V/N*kLAs[i]*inter.V/V*R*T/kHs[i])/(N*Cvave) - dT/dt * (dCvavedni/Cvave)
            #              -= (P*V/N*d/dni(dni/dt))/(N*Cvave) - dT/dt * (dCvavedni/Cvave)
            # d/dni(dP/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]*R*T/V + P/T*d/dni(dT/dt)
            #              -= d/dni(dni/dt)*R*T/V + P/T*d/dni(dT/dt)
            flow = sum(cond)
            @fastmath dTdt = (P * V / N * flow) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                ddnidnidt = kLAs[i] * inter.V / V * R * T / kHs[i]
                @inbounds @fastmath jac[i, i] -= ddnidnidt
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @fastmath ddnidTdt = (P * V / N * ddnidnidt) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] -= ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] -= ddnidnidt * R * T / V + P / T * ddnidTdt
            end
        elseif isa(inter, VolumetricFlowRateInlet) && domain == inter.domain
            # dn/dt += inter.Vin(t)*inter.cs
            # dT/dt += inter.Vin(t)*(inter.Hpervolume - dot(Us,ns)/V)/(N*Cvave)
            # dP/dt += inter.Vin(t)*P/V + P/T*dTdt
            # d/dni(dT/dt) += inter.Vin(t)*(-Us[i]/N)/(N*Cvave) - inter.Vin(t)*(inter.Hpervolume - dot(Us,ns)/V)/(N*Cvave)^2 * d/dni(N*Cvave)
            #              += inter.Vin(t)*(-Us[i]/N)/(N*Cvave) - dTdt * d/dni(N*Cvave)/(N*Cvave)
            #              += inter.Vin(t)*(-Us[i]/N)/(N*Cvave) - dTdt * d/dni(Cvave)/Cvave
            #                Note: Cvave = sum(ns.*cpdivR)*R/N - R
            #                Note: d/dni(Cvave) = cpdivR[i]*R/N
            # d/dni(dP/dt) += P/T * d/dni(dT/dt) 
            @fastmath dTdt = (inter.Vin(t) * (inter.Hpervolume - dot(Us, ns) / V)) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = inter.Vin(t) * (-Us[i] / N) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += P / T * ddnidTdt
            end

        elseif isa(inter, VolumetricFlowRateOutlet) && domain == inter.domain
            # dn/dt .-= inter.Vout(t)*ns/V
            # dT/dt -= (P*inter.Vout(t))/(N*Cvave)
            # dP/dt -= inter.Vout(t)*P/V + P/T*dTdt
            # d/dni(dni/dt) -= inter.Vout(t)/V
            # d/dni(dT/dt) -= dT/dt * (-1/(N*Cvave)) * d/dni(N*Cvave)
            #              -= dT/dt * (-1/(Cvave)) * d/dni(Cvave)
            #                 Note: Cvave = dot(cpdivR,ns)*R/N-R
            #                 Note: d/dni(Cvave) = cpdivR[i]*R/N
            #              -= -dT/dt * (dCvavedni/Cvave)
            # d/dni(dP/dt) -= P/T * d/dni(dT/dt)
            @fastmath dTdt = (P * inter.Vout(t)) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= inter.Vout(t) / V
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @fastmath ddnidTdt = -dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] -= ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] -= P / T * ddnidTdt
            end
        end
    end

    @inbounds jacobianytherm!(jac, y, p, t, domain, interfaces, domain.indexes[3], T, colorvec)
    @inbounds jacobianytherm!(jac, y, p, t, domain, interfaces, domain.indexes[4], P, colorvec)

    return jac
end

@inline function jacobiany!(jac::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ConstantPDomain}
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)
    jacobianynsderiv!(jac, domain, domain.rxnarray, domain.efficiencyinds, cs, kfs, krevs, T, V, C)

    dydt = zeros(size(y))
    addreactionratecontributions!(dydt, domain.rxnarray, cs, kfs, krevs)
    dydt .*= V

    # dT/dt = -dot(Hs,dn/dt)/(N*Cpave)
    # dV/dt = sum(dn/dt)*R*T/P + dT/dt*V/T
    # d/dni(dT/dt) = -dot(Hs,d/dni(dn/dt))/(N*Cpave) + dTdt * (-1/(N*Cpave)) * d/dni(N*Cpave)
    #              = -dot(Hs,d/dni(dn/dt))/(N*Cpave) + dTdt * (-1/(Cpave)) * d/dni(Cpave)
    #                Note: Cpave = dot(cpdivR,ns)*R/N-R
    #                Note: d/dni(Cpave) = cpdivR[i]*R/N
    #              = -dot(Hs,d/dni(dn/dt))/(N*Cpave) - dTdt * (dCpavedni/Cpave)
    # d/dV(dT/dt) = -dot(Hs,d/dV(dn/dt))/(N*Cpave) + dT/dt * (-1/(N*Cpave)) * d/dV(N*Cpave)
    #               Note: dN/dV = d(PV/RT)/dV = P/RT
    #               Note: Cpave = dot(cpdivR,ns)*R/N-R
    #               Note: d/dV(Cpave) = Cpave * (-1/N) * dN/dV = Cpave * (-1/N) * P/RT
    #               Note: d(N*Cpave)/dV = (dN/dV*Cpave + N*dCpave/dV) = P/RT*Cpave - Cpave * P/RT = 0
    #             = -dot(Hs,d/dV(dn/dt))/(N*Cpave)
    # d/dni(dV/dt) = sum(d/dni(dn/dt))*R*T/P + d/dni(dT/dt)*V/T
    # d/dV(dV/dt) = sum(d/dV(dn/dt))*R*T/P + d/dV(dT/dt)*V/T + dT/dt/T
    @fastmath Cpave = Cvave + R
    @views @inbounds @fastmath dTdt = -dot(Hs, dydt[domain.indexes[1]:domain.indexes[2]]) / (N * Cpave)
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @inbounds @fastmath dCpavedni = cpdivR[i] * R / N
        @views @inbounds @fastmath jac[domain.indexes[3], i] = -dot(Hs, jac[domain.indexes[1]:domain.indexes[2], i]) / (N * Cpave) - dTdt * (dCpavedni / Cpave)
        @views @inbounds @fastmath jac[domain.indexes[4], i] = sum(jac[domain.indexes[1]:domain.indexes[2], i]) * R * T / P + V / T * jac[domain.indexes[3], i]
    end
    @views @inbounds @fastmath jac[domain.indexes[3], domain.indexes[4]] = -dot(Hs, jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]]) / (N * Cpave)
    @views @inbounds @fastmath jac[domain.indexes[4], domain.indexes[4]] = sum(jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]]) * R * T / P + dTdt / T + V / T * jac[domain.indexes[3], domain.indexes[4]]

    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind, :] .= 0.0
    end

    @simd for inter in interfaces
        if isa(inter, Inlet) && domain == inter.domain
            # dn/dt .+= inter.y.*inter.F(t)
            # dT/dt += inter.F(t)*(inter.H - dot(Hs,ns)/N)/(N*Cpave)
            # dV/dt += inter.F(t)*R*T/P + dT/dt*V/T
            # d/dni(dn/dt) += 0
            # d/dV(dn/dt) += 0
            # d/dni(dT/dt) += inter.F(t)*(-Hs[i]/N)/(N*Cpave) - dTdt/(N*Cpave) * d/dni(N*Cpave)
            #              += inter.F(t)*(-Hs[i]/N)/(N*Cpave) - dTdt/(Cpave) * d/dni(Cpave)
            #                 Note: Cpave = dot(cpdivR,ns)*R/N-R
            #                 Note: d/dni(Cpave) = cpdivR[i]*R/N
            #              += inter.F(t)*(-Hs[i]/N)/(N*Cpave) - dTdt * (dCpavedni/Cpave)
            # d/dV(dT/dt) += inter.F(t)*(dot(Hs,ns)/N^2*dN/dV)/(N*Cpave) + dT/dt * (-1/(N*Cpave)) * d/dV(N*Cpave)
            #             += inter.F(t)*(dot(Hs,ns)/N^2*P/RT)/(N*Cpave) + 0
            #             += inter.F(t)*(dot(Hs,ns)/N^2*P/RT)/(N*Cpave)
            #             += inter.F(t)*(dot(Hs,ns)/N/V)/(N*Cpave)
            # d/dni(dV/dt) += d/dni(dT/dt)*V/T
            # d/dV(dV/dt) += d/dV(dT/dt)*V/T
            flow = inter.F(t)
            @fastmath H = dot(Hs, ns) / N
            @fastmath dTdt = flow * (inter.H - H) / (N * Cpave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCpavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = flow * (-Hs[i] / N) / (N * Cpave) - dTdt * (dCpavedni / Cpave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += V / T * ddnidTdt
            end
            @fastmath ddVdTdt = flow * H / V / (N * Cpave)
            @inbounds jac[domain.indexes[3], domain.indexes[4]] += ddVdTdt
            @inbounds @fastmath jac[domain.indexes[4], domain.indexes[4]] += V / T * ddVdTdt
        elseif isa(inter, Outlet) && domain == inter.domain
            # dn/dt .-= inter.F(t).*ns./N
            # dT/dt -= 0
            # dV/dt -= inter.F(t)*R*T/P
            # d/dni(dni/dt) -= inter.F(t)/N
            # d/dV(dni/dt) -= -inter.F(t).*ns./(N^2) * dN/dV
            #              -= -inter.F(t).*ns./(N^2) * P/RT
            #              -= -inter.F(t).*ns./N/V
            # d/dni(dT/dt) -= 0
            # d/dV(dT/dt) -= 0
            # d/dni(dV/dt) -= 0
            # d/dV(dV/dt) -= 0
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= flow / N
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]] .-= -flow * ns / N / V
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && domain == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs

            # evaporation
            # dn/dt .+= kLAs.*inter.V.*inter.cs
            # dT/dt += sum(kLAs.*inter.V.*inter.cs)*(inter.H - dot(Hs,ns)/N)/(N*Cpave)
            # dV/dt += sum(kLAs.*inter.V.*inter.cs)*R*T/P + dT/dt*V/T
            # d/dni(dn/dt) += 0
            # d/dV(dn/dt) += 0
            # d/dni(dT/dt) += sum(kLAs.*inter.V.*inter.cs)*(-Hs[i]/N)/(N*Cpave) - dTdt/(N*Cpave) * d/dni(N*Cpave)
            #              += sum(kLAs.*inter.V.*inter.cs)*(-Hs[i]/N)/(N*Cpave) - dTdt/(Cpave) * d/dni(Cpave)
            #                 Note: Cpave = dot(cpdivR,ns)*R/N-R
            #                 Note: d/dni(Cpave) = cpdivR[i]*R/N
            #              += sum(kLAs.*inter.V.*inter.cs)*(-Hs[i]/N)/(N*Cpave) - dTdt * (dCpavedni/Cpave)
            # d/dV(dT/dt) += sum(kLAs.*inter.V.*inter.cs)*(dot(Hs,ns)/N^2*dN/dV)/(N*Cpave) + dT/dt * (-1/(N*Cpave)) * d/dV(N*Cpave)
            #                Note: dN/dV = d(PV/RT)/dV = P/RT
            #                Note: Cpave = dot(cpdivR,ns)*R/N = dot(cpdivR,ns)*R*(RT/PV)
            #                Note: d/dV(Cpave) = dot(cpdivR,ns)*R*(RT/PV)*(-1/V) = -Cpave/V
            #                Note: d/dV(N*Cpave) = dN/dV*Cpave + N*dCpave/dV = P/RT*Cpave - (PV/RT)*Cpave/V = 0
            #             += sum(kLAs.*inter.V.*inter.cs)*(dot(Hs,ns)/N^2*P/RT)/(N*Cpave) + 0
            #             += sum(kLAs.*inter.V.*inter.cs)*(dot(Hs,ns)/N/V)/(N*Cpave)
            # d/dni(dV/dt) += d/dni(dT/dt)*V/T
            # d/dV(dV/dt) += d/dV(dT/dt)*V/T + dT/dt/T
            flow = sum(evap)
            @fastmath H = dot(Hs, ns) / N
            @fastmath dTdt = flow * (inter.H - H) / (N * Cpave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCpavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = flow * (-Hs[i] / N) / (N * Cpave) - dTdt * (dCpavedni / Cpave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += V / T * ddnidTdt
            end
            @fastmath ddVdTdt = flow * H / V / (N * Cpave)
            @inbounds jac[domain.indexes[3], domain.indexes[4]] += ddVdTdt
            @inbounds @fastmath jac[domain.indexes[4], domain.indexes[4]] += dTdt / T + V / T * ddVdTdt

            # condensation
            # dn/dt .-= kLAs.*inter.V.*cs*R*T./kHs
            # dT/dt -= 0
            # dV/dt -= sum(kLAs.*inter.V.*cs*R*T./kHs)*R*T/P
            # d/dni(dni/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]
            # d/dV(dn/dt) .-= kLAs.*inter.V.*cs*(-1/V)*R*T./kHs
            # d/dni(dT/dt) -= 0
            # d/dV(dT/dt) -= 0
            # d/dni(dV/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]*R*T/P
            #              -= ddnidnidt*R*T/P
            # d/dV(dV/dt) -= sum(d/dV(dn/dt))*R*T/P
            @simd for i in domain.indexes[1]:domain.indexes[2]
                ddnidnidt = kLAs[i] * inter.V / V * R * T / kHs[i]
                @inbounds @fastmath jac[i, i] -= ddnidnidt
                @inbounds @fastmath jac[domain.indexes[4], i] -= ddnidnidt * R * T / P
            end
            ddVdndt = kLAs .* inter.V .* cs * (-1 / V) * R * T ./ kHs
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]] .-= ddVdndt
            @views @inbounds @fastmath jac[domain.indexes[4], domain.indexes[4]] .-= sum(ddVdndt) * R * T / P
        elseif isa(inter, VolumetricFlowRateInlet) && domain == inter.domain
            # dn/dt .+= inter.Vin(t)*inter.cs
            # dT/dt += inter.Vin(t)*(inter.Hpervolume - dot(Hs,ns)/V)/(N*Cpave)
            # dV/dt += inter.Vin(t)
            # d/dn(dn/dt) = 0
            # d/dV(dn/dt) = 0
            # d/dni(dT/dt) += inter.Vin(t)*(-Hs[i]/N)/(N*Cpave) - dT/dt*(dCpavedni/Cpave)
            # d/dV(dT/dt) += inter.Vin(t)*(dot(Hs, ns)/V^2)/(N*Cpave) - (dT/dt)/(N*Cpave) * d/dV((PV/RT)*Cpave)
            #               Note: Cpave = dot(cpdivR,ns)*R/N = dot(cpdivR,ns)*R*(RT/PV)
            #               Note: d/dV(Cpave) = dot(cpdivR,ns)*R*(RT/PV)*(-1/V) = -Cpave/V
            #             = inter.Vin(t)*(dot(Hs, ns)/V^2)/(N*Cpave) - dT/dt/(N*Cpave) * ((P/RT)*Cpave - (PV/RT)*Cpave/V)
            #             = 0
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCpavedni = cpdivR[i] * R / N
                @inbounds jac[domain.indexes[3], i] += inter.Vin(t) * (-Hs[i] / N) / (N * Cpave) - dTdt * (dCpavedni / Cpave)
            end
        elseif isa(inter, VolumetricFlowRateOutlet) && domain == inter.domain
            # dn/dt .-= inter.Vout(t)*ns/V
            # dT/dt -= 0
            # dV/dt -= inter.Vout(t)
            # d/dni(dni/dt) -= inter.Vout(t)/V
            # d/dV(dn/dt) -= -inter.Vout(t)*ns/(V^2)
            # d/dni(dT/dt) -= 0
            # d/dV(dT/dt) -= 0
            # d/dni(dV/dt) -= 0
            # d/dV(dV/dt) -= 0
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= inter.Vout(t) / V
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]] .-= -inter.Vout(t) / (V * V) * ns
        end
    end

    @inbounds jacobianytherm!(jac, y, p, t, domain, interfaces, domain.indexes[3], T, colorvec)

    return jac
end

@inline function jacobiany!(jac::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ParametrizedTPDomain}
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)
    jacobianynsderiv!(jac, domain, domain.rxnarray, domain.efficiencyinds, cs, kfs, krevs, T, V, C)

    # dV/dt = sum(dn/dt)*R*T/P + Calculus.derivative(d.T,t)*V/T - Calculus.derivative(d.P,t)*V/P
    # d/dni(dV/dt) = sum(d/dni(dn/dt))*R*T/P
    # d/dV(dV/dt) = sum(d/dV(dn/dt))*R*T/P + Calculus.derivative(d.T,t)/T - Calculus.derivative(d.P,t)/P
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @views @inbounds @fastmath jac[domain.indexes[3], i] = sum(jac[domain.indexes[1]:domain.indexes[2], i]) * R * T / P
    end
    @views @inbounds @fastmath jac[domain.indexes[3], domain.indexes[3]] = sum(jac[domain.indexes[1]:domain.indexes[2], domain.indexes[3]]) * R * T / P + Calculus.derivative(domain.T, t) / T - Calculus.derivative(domain.P, t) / P

    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind, :] .= 0.0
    end

    @simd for inter in interfaces

        if isa(inter, Inlet) && domain == inter.domain
            # dn/dt .+= inter.y.*inter.F(t)
            # dV/dt += inter.F(t)*R*T/P
            # d/dni(dni/dt) += 0
            # d/dV(dni/dt) += 0
            # d/dni(dV/dt) += 0
            # d/dV(dV/dt) += 0
            nothing

        elseif isa(inter, Outlet) && domain == inter.domain
            # dn/dt .-= inter.F(t)*ns./N
            # dV/dt -= inter.F(t)*R*T/P
            # d/dni(dni/dt) -= inter.F(t)/N
            # d/dV(dn/dt) -= -inter.F(t)*ns/N^2 * dN/dV
            #                 Note: dN/dV = d(PV/RT)/dV = P/RT
            #              -= -inter.F(t)*ns/N^2 * P/RT
            #              -= -inter.F(t)*ns/N/V
            # d/dni(dV/dt) -= 0
            # d/dV(dV/dt) -= 0
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= flow / N
            end
            @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[3]] .-= -flow * ns / N / V
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && domain == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)

            # evaporaiton
            # dn/dt .+= kLAs.*inter.V.*inter.cs
            # dV/dt += sum(kLAs.*inter.V.*inter.cs)*R*T/P
            # d/dni(dni/dt) += 0
            # d/dV(dni/dt) += 0
            # d/dni(dV/dt) += 0
            # d/dV(dV/dt) += 0

            # condensation
            # dn/dt .-= kLAs.*inter.V.*cs*R*T./kHs
            # dV/dt -= sum(kLAs.*inter.V.*cs*R*T./kHs)*R*T/P
            # d/dni(dni/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]
            # d/dV(dn/dt) -= kLAs.*inter.V.*cs.*(-1/V)*R*T./kHs
            # d/dni(dV/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]*R*T/P
            #              -= d/dni(dni/dt)*R*T/P
            # d/dV(dV/dt) -= sum(d/dV(dn/dt))*R*T/P
            @simd for i in domain.indexes[1]:domain.indexes[2]
                ddnidnidt = kLAs[i] * inter.V / V * R * T / kHs[i]
                @inbounds @fastmath jac[i, i] -= ddnidnidt
                @inbounds @fastmath jac[domain.indexes[3], i] -= ddnidnidt * R * T / P
            end
            ddVdndt = kLAs .* inter.V .* cs .* (-1 / V) * R * T ./ kHs
            @views @inbounds @fastmath jac[domain.indexs[1]:domain.indexes[2], domain.indexes[3]] .-= ddVdndt
            @views @inbounds @fastmath jac[domain.indexes[3], domain.indexes[3]] -= sum(ddVdndt) * R * T / P
        elseif isa(inter, VolumetricFlowRateInlet) && domain == inter.domain
            # dn/dt .+= inter.Vin(t)*inter.cs
            # dV/dt += inter.Vin(t)
            nothing
        elseif isa(inter, VolumetricFlowRateOutlet) && domain == inter.domain
            # dn/dt .-= inter.Vout(t)*ns/V
            # dV/dt -= inter.Vout(t)
            # d/dni(dni/dt) -= inter.Vout(t)/V
            # d/dV(dni/dt) -= -inter.Vout(t)*ns/(V*V)
            # d/dni(dV/dt) -= 0
            # d/dV(dV/dt) -= 0
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= inter.Vout(t) / V
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[3]] .-= -inter.Vout(t) * ns / (V * V)
        end
    end

    return jac
end

@inline function jacobiany!(jac::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ParametrizedVDomain}
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)
    jacobianynsderiv!(jac, domain, domain.rxnarray, domain.efficiencyinds, cs, kfs, krevs, T, V, C)

    dydt = zeros(size(y))
    addreactionratecontributions!(dydt, domain.rxnarray, cs, kfs, krevs)
    dydt .*= V

    # dT/dt = (-dot(Us,dn/dt)-P*dV/dt(t))/(N*Cvave)
    # dP/dt = sum(dn/dt)*R*T/V + dT/dt*P/T - P/V*dVdt
    # d/dni(dT/dt) = -dot(Us, d/dni(dn/dt))/(N*Cvave) + dT/dt * (-1/N*Cvave) * d/dni(N*Cvave)
    #              = -dot(Us, d/dni(dn/dt))/(N*Cvave) + dT/dt * (-1/Cvave) * d/dni(Cvave)
    #                Note: Cvave = dot(cpdivR,ns)*R/N-R
    #                Note: d/dni(Cvave) = cpdivR[i]*R/N
    #              = -dot(Us, d/dni(dn/dt))/(N*Cvave) - dT/dt * (d/dni(Cvave)/Cvave)
    # d/dni(dP/dt) = sum(d/dni(dn/dt))*R*T/V + d/dni(dT/dt)*P/T
    dVdt = Calculus.derivative(domain.V, t)
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @fastmath @inbounds dCvavedni = cpdivR[i] * R / N
        @views @inbounds @fastmath jac[domain.indexes[3], i] = -dot(Us, jac[domain.indexes[1]:domain.indexes[2], i]) / (N * Cvave) - (-dot(Us, dydt[domain.indexes[1]:domain.indexes[2]]) - P * dVdt) / (N * Cvave * Cvave) * dCvavedni
        @views @inbounds @fastmath jac[domain.indexes[4], i] = sum(jac[domain.indexes[1]:domain.indexes[2], i]) * R * T / V + P / T * jac[domain.indexes[3], i]
    end

    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind, :] .= 0.0
    end

    @simd for inter in interfaces
        if isa(inter, Inlet) && domain == inter.domain
            # dn/dt .+= inter.y.*inter.F(t)
            # dT/dt += inter.F(t)*(inter.H - dot(Us,ns)/N)/(N*Cvave)
            # dP/dt += inter.F(t)*R*T/V + dT/dt*P/T
            # d/dni(dni/dt) += 0
            # d/dT(dn/dt) += 0
            # d/dP(dn/dt) += 0
            # d/dni(dT/dt) += inter.F(t)*(-Us[i]/N)/(N*Cvave) + dT/dt * (-1/N*Cvave) * d/dni(N*Cvave)
            #              += inter.F(t)*(-Us[i]/N)/(N*Cvave) + dT/dt * (-1/Cvave) * d/dni(Cvave)
            # d/dni(dP/dt) += d/dni(dT/dt)*P/T
            flow = inter.F(t)
            @fastmath dTdt = flow * (inter.H - dot(Us, ns) / N) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = flow * (-Us[i] / N) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += P / T * ddnidTdt
            end
        elseif isa(inter, Outlet) && domain == inter.domain
            # dn/dt .-= inter.F(t)*ns./N
            # dT/dt -= (P*V/N*inter.F(t))/(N*Cvave)
            # dP/dt -= inter.F(t)*R*T/V + dT/dt*P/T
            # d/dni(dni/dt) -= inter.F(t)/N
            # d/dni(dT/dt) -= dT/dt * (-1/N*Cvave) * d/dni(N*Cvave)
            #              -= dT/dt * (-1/Cvave) * d/dni(Cvave)
            # d/dni(dP/dt) -= d/dni(dT/dt)*P/T
            flow = inter.F(t)
            @fastmath dTdt = (P * V / N * flow) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= flow / N
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @fastmath ddnidTdt = -dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] -= ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] -= P / T * ddnidTdt
            end
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && domain == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs
            cond = kLAs .* inter.V .* cs * R * T ./ kHs

            # evaporation
            # dn/dt .+= kLAs.*inter.V.*inter.cs
            # dT/dt += sum(kLAs.*inter.V.*inter.cs)*(inter.H - dot(Us,ns)/N)/(N*Cvave)
            # dP/dt += sum(kLAs.*inter.V.*inter.cs)*R*T/V + dT/dt*P/T
            # d/dni(dni/dt) += 0
            # d/dni(dT/dt) += sum(kLAs.*inter.V.*inter.cs)*(-Us[i]/N)/(N*Cvave) + dT/dt * (-1/N*Cvave) * d/dni(N*Cvave)
            #              += sum(kLAs.*inter.V.*inter.cs)*(-Us[i]/N)/(N*Cvave) + dT/dt * (-1/Cvave) * d/dni(Cvave)
            # d/dni(dP/dt) += d/dni(dT/dt)*P/T
            flow = sum(evap)
            @fastmath dTdt = flow * (inter.H - dot(Us, ns) / N) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = flow * (-Us[i] / N) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += P / T * ddnidTdt
            end

            # condensation
            # dn/dt .-= kLAs.*inter.V.*cs*R*T./kHs
            # dT/dt -= (P*V/N*sum(kLAs.*inter.V.*cs*R*T./kHs))/(N*Cvave)
            # dP/dt -= sum(kLAs.*inter.V.*cs*R*T./kHs)*R*T/V + dT/dt*P/T
            # d/dni(dni/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]
            # d/dni(dT/dt) -= (P*V/N*ddnidnidt)/(N*Cvave) - dTdt*(dCvavedni/Cvave)
            # d/dni(dP/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]*R*T/V + P/T*ddnidTdt
            #              -= d/dni(dni/dt)*R*T/V + P/T*ddnidTdt
            flow = sum(cond)
            @fastmath dTdt = (P * V / N * flow) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                ddnidnidt = kLAs[i] * inter.V / V * R * T / kHs[i]
                @inbounds @fastmath jac[i, i] -= ddnidnidt
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @fastmath ddnidTdt = (P * V / N * ddnidnidt) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] -= ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] -= ddnidnidt * R * T / V + P / T * ddnidTdt
            end
        elseif isa(inter, VolumetricFlowRateInlet) && domain == inter.domain
            # dn/dt .+= inter.Vin(t)*inter.cs
            # dT/dt += inter.Vin(t)*(inter.Hpervolume - dot(Us,ns)/V)/(N*Cvave)
            # dP/dt += inter.Vin(t)*P/V + dT/dt*P/T
            # d/dn(dn/dt) += 0
            # d/dni(dT/dt) += inter.Vin(t)*(-Us[i]/V)/(N*Cvave) - inter.Vin(t)*(inter.Hpervolume - dot(Us,ns)/V)/(N*Cvave)^2 * d/dni(N*Cvave)
            #              = inter.Vin(t)*(-Us[i]/V)/(N*Cvave) - dT/dt * d/dni(N*Cvave) / (N*Cvave)
            #              = inter.Vin(t)*(-Us[i]/V)/(N*Cvave) - dT/dt * dCvave/dni / Cvave
            #                Note: Cvave = sum(ns.*cpdivR)*R/N - R
            #                Note: dCvave/dni = cpdivR[i]*R/N
            # d/dni(dP/dt) += d/dni(dT/dt)*P/T
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @fastmath ddnidTdt = (inter.Vout * P / N) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                ddnidTdt = inter.Vin(t) * (-Us[i] / V) / (N * Cvave) - dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += ddnidTdt * P / T
            end
        elseif isa(inter, VolumetricFlowRateOutlet) && domain == inter.domain
            # dn/dt .-= inter.Vout(t)*ns/V
            # dT/dt -= (P*inter.Vout(t))/(N*Cvave)
            # dP/dt -= inter.Vout(t)*P/V + dT/dt*P/T
            # d/dni(dni/dt) -= inter.Vout(t)/V
            # d/dni(dT/dt) -= dT/dt *(-1/(N*Cvave)) d/dni(N*Cvave)
            #              -= dT/dt *(-1/(Cvave)) d/dni(Cvave)
            #              -= - dT/dt * (dCvavedni/Cvave)
            # d/dni(dP/dt) -= P/T * d/dni(dT/dt)
            @fastmath dTdt = (P * inter.Vout(t)) / (N * Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= inter.Vout(t) / V
                @inbounds @fastmath dCvavedni = cpdivR[i] * R / N
                @fastmath ddnidTdt = -dTdt * (dCvavedni / Cvave)
                @inbounds jac[domain.indexes[3], i] -= ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] -= P / T * ddnidTdt
            end
        end
    end

    @inbounds jacobianytherm!(jac, y, p, t, domain, interfaces, domain.indexes[3], T, colorvec)
    @inbounds jacobianytherm!(jac, y, p, t, domain, interfaces, domain.indexes[4], P, colorvec)

    return jac
end

@inline function jacobiany!(jac::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ParametrizedPDomain}
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)
    jacobianynsderiv!(jac, domain, domain.rxnarray, domain.efficiencyinds, cs, kfs, krevs, T, V, C)

    dydt = zeros(size(y))
    addreactionratecontributions!(dydt, domain.rxnarray, cs, kfs, krevs)
    dydt .*= V

    # dT/dt = (-dot(Hs,dn/dt)+V*dPdt(t))/(N*Cpave)
    # dV/dt = sum(dn/dt)*R*T/P + dT/dt*V/T - dPdt(t)*V/P
    # d/dni(dT/dt) = -dot(Hs,d/dni(dn/dt))/(N*Cpave)-dTdt*(dCpavedni/Cpave)
    # d/dV(dT/dt) = (-dot(Hs,d/dV(dn/dt)) + dPdt(t))/(N*Cpave) + dT/dt * (-1/(N*Cpave)) * d/dV(N*Cpave)
    #               Note: d(N*Cpave)/dV = (dN/dV*Cpave + N*dCpave/dV) = P/RT*Cpave - Cpave * P/RT = 0
    #             = (-dot(Hs,d/dV(dn/dt)) + dPdt(t))/(N*Cpave) + dT/dt * (-1/(N*Cpave))
    # d/dni(dV/dt) = sum(d/dni(dn/dt))*R*T/P + V/T*d/dni(dT/dt)
    # d/dV(dV/dt) = sum(d/dV(dn/dt))*R*T/P + d/dV(dT/dt)*V/T + dT/dt/T - dPdt(t)/P
    @fastmath Cpave = Cvave + R
    dPdt = Calculus.derivative(domain.P, t)
    @views @inbounds @fastmath dTdt = (-dot(Hs, dydt[domain.indexes[1]:domain.indexes[2]]) + V * dPdt) / (N * Cpave)
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @inbounds @fastmath dCpavedni = cpdivR[i] * R / N
        @views @inbounds @fastmath jac[domain.indexes[3], i] = -dot(Hs, jac[domain.indexes[1]:domain.indexes[2], i]) / (N * Cpave) - dTdt * (dCpavedni / Cpave)
        @views @inbounds @fastmath jac[domain.indexes[4], i] = sum(jac[domain.indexes[1]:domain.indexes[2], i]) * R * T / P + V / T * jac[domain.indexes[3], i]
    end
    @views @inbounds @fastmath jac[domain.indexes[3], domain.indexes[4]] = (-dot(Hs, jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]]) + dPdt) / (N * Cpave)
    @views @inbounds @fastmath jac[domain.indexes[4], domain.indexes[4]] = sum(jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]]) * R * T / P + dTdt / T + V / T * jac[domain.indexes[3], domain.indexes[4]] - dPdt / P

    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind, :] .= 0.0
    end

    @simd for inter in interfaces
        if isa(inter, Inlet) && domain == inter.domain
            # dn/dt .+= inter.y.*inter.F(t)
            # dT/dt += inter.F(t)*(inter.H - dot(Hs,ns)/N)/(N*Cpave)
            # dV/dt += inter.F(t)*R*T/P + dT/dt*V/T
            # d/dni(dni/dt) += 0
            # d/dV(dni/dt) += 0
            # d/dni(dT/dt) += inter.F(t)*(-Hs[i]/N)/(N*Cpave) + dTdt * (-1/(N*Cpave)) * d/dni(N*Cpave)
            #              += inter.F(t)*(-Hs[i]/N)/(N*Cpave) + dTdt * (-1/(Cpave)) * d/dni(Cpave)
            #                 Note: Cpave = = dot(cpdivR,ns)*R/N
            #                 Note: dCpavedni = cpdivR[i]*R/N
            # d/dV(dT/dt) += inter.F(t)*(dot(Hs,ns)/N^2*dN/dV)/(N*Cpave) + dTdt * (-1/(N*Cpave)) * d/dV(N*Cpave)
            #                Note: dN/dV = d(PV/RT)/dV = P/RT
            #                Note: d(N*Cpave)/dV = (dN/dV*Cpave + N*dCpave/dV) = P/RT*Cpave - Cpave * P/RT = 0
            #             += inter.F(t)*(dot(Hs,ns)/N/V)/(N*Cpave)
            # d/dni(dV/dt) += V/T*d/dni(dT/dt)
            # d/dV(dV/dt) += d/dV(dT/dt)*V/T + dT/dt/T
            flow = inter.F(t)
            @fastmath H = dot(Hs, ns) / N
            @fastmath dTdt = flow * (inter.H - H) / (N * Cpave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCpavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = flow * (-Hs[i] / N) / (N * Cpave) - dTdt * (dCpavedni / Cpave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += V / T * ddnidTdt
            end
            @fastmath ddVdTdt = flow * H / V / (N * Cpave)
            @inbounds jac[domain.indexes[3], domain.indexes[4]] += ddVdTdt
            @inbounds @fastmath jac[domain.indexes[4], domain.indexes[4]] += dTdt / T + V / T * ddVdTdt
        elseif isa(inter, Outlet) && domain == inter.domain
            # dn/dt .-= inter.F(t).*ns./N
            # dT/dt -= 0
            # dV/dt -= inter.F(t)*R*T/P
            # d/dni(dni/dt) -= inter.F(t)/N
            # d/dV(dn/dt) -= -inter.F(t).*ns./N^2 * dN/dV
            #                Note: dN/dV = d(PV/RT)/dV = P/RT
            #             -= inter.F(t).*ns./N/V
            # d/dni(dT/dt) -= 0
            # d/dV(dT/dt) -= 0
            # d/dni(dV/dt) -= 0
            # d/dV(dV/dt) -= 0
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= flow / N
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]] .-= -flow * ns / N / V
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && domain == inter.domain
            kLAs = map.(inter.kLAs, inter.T)
            kHs = map.(inter.kHs, inter.T)
            evap = kLAs .* inter.V .* inter.cs

            # evaporation
            # dn/dt .+= kLAs.*inter.V.*inter.cs
            # dT/dt += sum(kLAs.*inter.V.*inter.cs)*(inter.H - dot(Hs,ns)/N)/(N*Cpave)
            # dV/dt += sum(kLAs.*inter.V.*inter.cs)*R*T/P + dT/dt*V/T
            # d/dni(dni/dt) += 0
            # d/dV(dni/dt) += 0
            # d/dni(dT/dt) += sum(kLAs.*inter.V.*inter.cs)*(-Hs[i]/N)/(N*Cpave) + dTdt * (-1/(N*Cpave)) * d/dni(N*Cpave)
            #              += sum(kLAs.*inter.V.*inter.cs)*(-Hs[i]/N)/(N*Cpave) + dTdt * (-1/(Cpave)) * d/dni(Cpave)
            #                 Note: Cpave = dot(cpdivR,ns)*R/N
            #                 Note: dCpavedni = cpdivR[i]*R/N
            # d/dV(dT/dt) += sum(kLAs.*inter.V.*inter.cs)*(dot(Hs,ns)/N^2*dN/dV)/(N*Cpave) + dTdt * (-1/(N*Cpave)) * d/dV(N*Cpave)
            #                Note: dN/dV = d(PV/RT)/dV = P/RT
            #                Note: Cpave = dot(cpdivR,ns)*R/N
            #                Note: d/dV(Cpave) = dot(cpdivR,ns)*R*(RT/PV)*(-1/V) = -Cpave/V
            #                Note: d(N*Cpave)/dV = (dN/dV*Cpave + N*dCpave/dV) = (P/RT*Cpave + PV/RT(-Cpave/V)) = 0
            #             += sum(kLAs.*inter.V.*inter.cs)*(dot(Hs,ns)/N/V)/(N*Cpave)
            # d/dni(dV/dt) += V/T*d/dni(dT/dt)
            # d/dV(dV/dt) += d/dV(dT/dt)*V/T + dT/dt/T
            flow = sum(evap)
            @fastmath H = dot(Hs, ns) / N
            @fastmath dTdt = flow * (inter.H - H) / (N * Cpave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCpavedni = cpdivR[i] * R / N
                @inbounds @fastmath ddnidTdt = flow * (-Hs[i] / N) / (N * Cpave) - dTdt * (dCpavedni / Cpave)
                @inbounds jac[domain.indexes[3], i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4], i] += V / T * ddnidTdt
            end
            @fastmath ddVdTdt = flow * H / V / (N * Cpave)
            @inbounds jac[domain.indexes[3], domain.indexes[4]] += ddVdTdt
            @inbounds @fastmath jac[domain.indexes[4], domain.indexes[4]] += dTdt / T + V / T * ddVdTdt

            # condensation
            # dn/dt .-= kLAs.*inter.V.*cs*R*T./kHs
            # dT/dt -= 0
            # dV/dt -= sum(kLAs.*inter.V.*cs*R*T./kHs)*R*T/P
            # d/dni(dni/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]
            # d/dV(dn/dt) -= kLAs.*inter.V.*cs.*(-1/V)*R*T./kHs
            # d/dni(dT/dt) -= 0
            # d/dV(dT/dt) -= 0
            # d/dni(dV/dt) -= kLAs[i]*inter.V/V*R*T/kHs[i]*R*T/P
            #              -= d/dni(dni/dt)*R*T/P
            # d/dV(dV/dt) -= sum(d/dV(dn/dt))*R*T/P
            @simd for i in domain.indexes[1]:domain.indexes[2]
                ddnidnidt = kLAs[i] * inter.V / V * R * T / kHs[i]
                @inbounds @fastmath jac[i, i] -= ddnidnidt
                @inbounds @fastmath jac[domain.indexes[4], i] -= ddnidnidt * R * T / P
            end
            ddVdndt = kLAs .* inter.V .* cs .* (-1 / V) * R * T ./ kHs
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]] .-= ddVdndt
            @inbounds @fastmath jac[domain.indexes[4], domain.indexes[4]] -= sum(ddVdndt) * R * T / P
        elseif isa(inter, VolumetricFlowRateInlet) && domain == inter.domain
            # dn/dt .+= inter.Vin(t)*inter.cs
            # dT/dt += inter.Vin(t)*(inter.Hpervolume - dot(Hs,ns)/V)/(N*Cpave)
            # dV/dt += inter.Vin(t)
            # d/dn(dn/dt) += 0
            # d/dV(dn/dt) += 0
            # d/dni(dT/dt) += inter.V(t)*(-Hs[i]/V)/(N*Cpave) - dT/dt/(N*Cpave) * d/dni(N*Cpave)
            #                Note: Cpave = dot(cpdivR,ns)*R/N
            #                Note: d/dni(Cpave) = cpdivR[i]*R/N
            #              = inter.V(t)*(-Hs[i]/V)/(N*Cpave) - dT/dt * (dCpavedni/Cpave)
            # d/dV(dT/dt) += inter.Vin(t)*(dot(Hs,ns)/V^2)/(N*Cpave) - dT/dt/(N*Cpave) * d/dV(PV/RT*Cpave)
            #             = inter.Vin(t)*(dot(Hs,ns)/V^2)/(N*Cpave) - dT/dt/(N*Cpave) * (P/RT*Cpave + PV/RT*dCpavedV)
            #               Note: Cpave = dot(cpdivR,ns)*R/N = dot(cpdivR,ns)*R*(RT/PV)
            #               Note: d/dV(Cpave) = dot(cpdivR,ns)*R*(RT/PV)*(-1/V) = -Cpave/V
            #             = inter.Vin(t)*(dot(Hs,ns)/V^2)/(N*Cpave) - dT/dt/(N*Cpave) * (P/RT*Cpave + PV/RT(-Cpave/V))
            #             = 0
            dTdt = inter.Vin(t) * (inter.Hpervolume - dot(Hs, ns) / V) / (N * Cpave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCpavedni = cpdivR[i] * R / N
                @inbounds @fastmath jac[domain.indexes[3], i] += inter.Vin(t) * (-Hs[i] / V) / (N * Cpave) - dTdt * (dCpavedni / Cpave)
            end
        elseif isa(inter, VolumetricFlowRateOutlet) && domain == inter.domain
            # dn/dt .-= inter.Vout(t)*ns/V
            # dT/dt -= 0
            # dV/dt -= inter.Vout(t)
            # d/dni(dni/dt) -= inter.Vout(t)/V
            # d/dV(dni/dt) -= -inter.Vout(t)*ns/V^2
            # d/dni(dT/dt) -= 0
            # d/dV(dT/dt) -= 0
            # d/dni(dV/dt) -= 0
            # d/dV(dV/dt) -= 0
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= inter.Vout(t) / V
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2], domain.indexes[4]] .-= -inter.Vout(t) / (V * V) * ns
        end
    end

    @inbounds jacobianytherm!(jac, y, p, t, domain, interfaces, domain.indexes[3], T, colorvec)

    return jac
end

@inline function jacobiany!(jac::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTVDomain,ConstantTAPhiDomain,ParametrizedTConstantVDomain}}
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, phi = calcthermo(domain, y, t, p)
    jacobianynsderiv!(jac, domain, domain.rxnarray, domain.efficiencyinds, cs, kfs, krevs, T, V, C)

    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind, :] .= 0.0
    end

    @simd for inter in interfaces
        if isa(inter, Outlet) && domain == inter.domain
            # dn/dt .-= inter.F(t).*ns./N
            # Note N = sum(ns) for liquid phase
            # d/dni(dni/dt) -= inter.F(t)/N - inter.F(t)*ns[i]/(N*N)
            # d/dnj(dni/dt) -= -inter.F(t)*ns[i]/(N*N)
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= flow / N
                @inbounds @fastmath jac[i, :] .+= flow * ns[i] / (N * N)
            end
        elseif isa(inter, kLAkHCondensationEvaporationWithReservoir) && domain == inter.domain
            kLAs = map.(inter.kLAs, T)

            # evaporation
            # dn/dt .-= kLAs.*ns
            # d/dni(dni/dt) -= kLAs[i]
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= kLAs[i]
            end

            # condensation
            # dn/dt .+= kLAs.*inter.molefractions*inter.P/kHs*V
            nothing
        elseif isa(inter, VolumetricFlowRateOutlet) && domain == inter.domain
            # dn/dt .-= inter.Vout(t)*ns/V
            # d/dni(dni/dt) -= inter.Vout(t)/V
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i, i] -= inter.Vout(t) / V
            end
        end
    end

    return jac
end

export jacobiany!

@inline function jacobianpnsderiv!(jacp::Q, y::U, p::W, t::Z, domain::D, rxnarray::Array{Int64,2}, cs::Array{Float64,1}, T::Float64, V::Float64, kfs::Array{Float64,1}, krevs::Array{Float64,1}, Nspcs::Int64, Nrxns::Int64) where {Q3<:AbstractArray,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTPDomain,ConstantTAPhiDomain}}
    @fastmath RTinv = 1.0 / (R * T)
    @simd for rxnind = 1:Nrxns
        if rxnind in domain.efficiencyinds
            if @inbounds rxnarray[2, rxnind] == 0
                @fastmath @inbounds fderiv = kfs[rxnind] * cs[rxnarray[1, rxnind]]
            elseif @inbounds rxnarray[3, rxnind] == 0
                @fastmath @inbounds fderiv = kfs[rxnind] * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]]
            elseif @inbounds rxnarray[4, rxnind] == 0
                @fastmath @inbounds fderiv = kfs[rxnind] * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]] * cs[rxnarray[3, rxnind]]
            else
                @fastmath @inbounds fderiv = kfs[rxnind] * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]] * cs[rxnarray[3, rxnind]] * cs[rxnarray[4, rxnind]]
            end

            if @inbounds rxnarray[6, rxnind] == 0
                @fastmath @inbounds rderiv = krevs[rxnind] * cs[rxnarray[5, rxnind]]
            elseif @inbounds rxnarray[7, rxnind] == 0
                @fastmath @inbounds rderiv = krevs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]]
            elseif @inbounds rxnarray[8, rxnind] == 0
                @fastmath @inbounds rderiv = krevs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]] * cs[rxnarray[7, rxnind]]
            else
                @fastmath @inbounds rderiv = krevs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]] * cs[rxnarray[7, rxnind]] * cs[rxnarray[8, rxnind]]
            end
        else
            if @inbounds rxnarray[2, rxnind] == 0
                @inbounds fderiv = cs[rxnarray[1, rxnind]]
            elseif @inbounds rxnarray[3, rxnind] == 0
                @fastmath @inbounds fderiv = cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]]
            elseif @inbounds rxnarray[4, rxnind] == 0
                @fastmath @inbounds fderiv = cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]] * cs[rxnarray[3, rxnind]]
            else
                @fastmath @inbounds fderiv = cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]] * cs[rxnarray[3, rxnind]] * cs[rxnarray[4, rxnind]]
            end

            if @inbounds rxnarray[6, rxnind] == 0
                @fastmath @inbounds rderiv = krevs[rxnind] / kfs[rxnind] * cs[rxnarray[5, rxnind]]
            elseif @inbounds rxnarray[7, rxnind] == 0
                @fastmath @inbounds rderiv = krevs[rxnind] / kfs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]]
            elseif @inbounds rxnarray[8, rxnind] == 0
                @fastmath @inbounds rderiv = krevs[rxnind] / kfs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]] * cs[rxnarray[7, rxnind]]
            else
                @fastmath @inbounds rderiv = krevs[rxnind] / kfs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]] * cs[rxnarray[7, rxnind]] * cs[rxnarray[8, rxnind]]
            end
        end

        flux = fderiv - rderiv
        _spreadreactantpartials!(jacp, flux, rxnarray, rxnind, Nspcs + rxnind)
        _spreadproductpartials!(jacp, -flux, rxnarray, rxnind, Nspcs + rxnind)

        @fastmath @inbounds gderiv = rderiv * kfs[rxnind] * RTinv

        @inbounds jacp[rxnarray[1, rxnind], rxnarray[1, rxnind]] -= gderiv
        @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[1, rxnind])
        if @inbounds rxnarray[2, rxnind] !== 0
            @inbounds jacp[rxnarray[2, rxnind], rxnarray[1, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[1, rxnind], rxnarray[2, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[2, rxnind], rxnarray[2, rxnind]] -= gderiv
            @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[2, rxnind])
            if @inbounds rxnarray[3, rxnind] !== 0
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[1, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[2, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[1, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[2, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[3, rxnind])
                if @inbounds rxnarray[4, rxnind] !== 0
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[1, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[2, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[3, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[1, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[2, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[3, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[4, rxnind])
                end
            end
        end

        @inbounds jacp[rxnarray[5, rxnind], rxnarray[5, rxnind]] -= gderiv
        @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[5, rxnind])
        if @inbounds rxnarray[6, rxnind] !== 0
            @inbounds jacp[rxnarray[6, rxnind], rxnarray[5, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[5, rxnind], rxnarray[6, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[6, rxnind], rxnarray[6, rxnind]] -= gderiv
            @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[6, rxnind])
            if @inbounds rxnarray[7, rxnind] !== 0
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[5, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[6, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[5, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[7, rxnind])
                if @inbounds rxnarray[8, rxnind] !== 0
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[5, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[6, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[7, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[5, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[6, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[7, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[8, rxnind])
                end
            end
        end

    end
    jacp .*= V
end

@inline function jacobianpnsderiv!(jacp::Q, y::U, p::W, t::Z, domain::D, rxnarray::Array{Int64,2}, cs::Array{Float64,1}, T::Float64, V::Float64, kfs::Array{Float64,1}, krevs::Array{Float64,1}, Nspcs::Int64, Nrxns::Int64) where {Q3<:AbstractArray,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantVDomain,ConstantPDomain,ParametrizedTPDomain,ParametrizedVDomain,ParametrizedPDomain,ParametrizedTConstantVDomain}}
    @fastmath RTinv = 1.0 / (R * T)
    @simd for rxnind = 1:Nrxns
        if @inbounds rxnarray[2, rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind] * cs[rxnarray[1, rxnind]]
        elseif @inbounds rxnarray[3, rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind] * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]]
        elseif @inbounds rxnarray[4, rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind] * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]] * cs[rxnarray[3, rxnind]]
        else
            @fastmath @inbounds fderiv = kfs[rxnind] * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]] * cs[rxnarray[3, rxnind]] * cs[rxnarray[4, rxnind]]
        end

        if @inbounds rxnarray[6, rxnind] == 0
            @fastmath @inbounds rderiv = krevs[rxnind] * cs[rxnarray[5, rxnind]]
        elseif @inbounds rxnarray[7, rxnind] == 0
            @fastmath @inbounds rderiv = krevs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]]
        elseif @inbounds rxnarray[8, rxnind] == 0
            @fastmath @inbounds rderiv = krevs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]] * cs[rxnarray[7, rxnind]]
        else
            @fastmath @inbounds rderiv = krevs[rxnind] * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]] * cs[rxnarray[7, rxnind]] * cs[rxnarray[8, rxnind]]
        end

        @fastmath flux = fderiv - rderiv
        _spreadreactantpartials!(jacp, flux, rxnarray, rxnind, Nspcs + rxnind)
        _spreadproductpartials!(jacp, -flux, rxnarray, rxnind, Nspcs + rxnind)

        @fastmath gderiv = rderiv * RTinv

        @inbounds jacp[rxnarray[1, rxnind], rxnarray[1, rxnind]] -= gderiv
        @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[1, rxnind])
        if @inbounds rxnarray[2, rxnind] !== 0
            @inbounds jacp[rxnarray[2, rxnind], rxnarray[1, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[1, rxnind], rxnarray[2, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[2, rxnind], rxnarray[2, rxnind]] -= gderiv
            @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[2, rxnind])
            if @inbounds rxnarray[3, rxnind] !== 0
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[1, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[2, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[1, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[2, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[3, rxnind])
                if @inbounds rxnarray[4, rxnind] !== 0
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[1, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[2, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[3, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[1, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[2, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[3, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[4, rxnind])
                end
            end
        end

        @inbounds jacp[rxnarray[5, rxnind], rxnarray[5, rxnind]] -= gderiv
        @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[5, rxnind])
        if @inbounds rxnarray[6, rxnind] !== 0
            @inbounds jacp[rxnarray[6, rxnind], rxnarray[5, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[5, rxnind], rxnarray[6, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[6, rxnind], rxnarray[6, rxnind]] -= gderiv
            @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[6, rxnind])
            if @inbounds rxnarray[7, rxnind] !== 0
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[5, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[6, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[5, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[7, rxnind])
                if @inbounds rxnarray[8, rxnind] !== 0
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[5, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[6, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[7, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[5, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[6, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[7, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[8, rxnind])
                end
            end
        end
    end
    jacp .*= V
end

@inline function jacobianpnsderiv!(jacp::Q, y::U, p::W, t::Z, domain::D, rxnarray::Array{Int64,2}, cs::Array{Float64,1}, T::Float64, V::Float64, kfs::Array{Float64,1}, krevs::Array{Float64,1}, Nspcs::Int64, Nrxns::Int64) where {Q3<:AbstractArray,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTVDomain}}
    @fastmath RTinv = 1.0 / (R * T)
    @simd for rxnind = 1:Nrxns
        if @inbounds rxnarray[2, rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind] * kfs[rxnind] / (p[Nspcs+rxnind] * p[Nspcs+rxnind]) * cs[rxnarray[1, rxnind]]
        elseif @inbounds rxnarray[3, rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind] * kfs[rxnind] / (p[Nspcs+rxnind] * p[Nspcs+rxnind]) * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]]
        elseif rxnarray[4, rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind] * kfs[rxnind] / (p[Nspcs+rxnind] * p[Nspcs+rxnind]) * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]] * cs[rxnarray[3, rxnind]]
        else
            @fastmath @inbounds fderiv = kfs[rxnind] * kfs[rxnind] / (p[Nspcs+rxnind] * p[Nspcs+rxnind]) * cs[rxnarray[1, rxnind]] * cs[rxnarray[2, rxnind]] * cs[rxnarray[3, rxnind]] * cs[rxnarray[4, rxnind]]
        end

        if @inbounds rxnarray[6, rxnind] == 0
            @fastmath @inbounds rderiv = kfs[rxnind] * krevs[rxnind] / (p[Nspcs+rxnind] * p[Nspcs+rxnind]) * cs[rxnarray[5, rxnind]]
        elseif @inbounds rxnarray[7, rxnind] == 0
            @fastmath @inbounds rderiv = kfs[rxnind] * krevs[rxnind] / (p[Nspcs+rxnind] * p[Nspcs+rxnind]) * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]]
        elseif @inbounds rxnarray[8, rxnind] == 0
            @fastmath @inbounds rderiv = kfs[rxnind] * krevs[rxnind] / (p[Nspcs+rxnind] * p[Nspcs+rxnind]) * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]] * cs[rxnarray[7, rxnind]]
        else
            @fastmath @inbounds rderiv = kfs[rxnind] * krevs[rxnind] / (p[Nspcs+rxnind] * p[Nspcs+rxnind]) * cs[rxnarray[5, rxnind]] * cs[rxnarray[6, rxnind]] * cs[rxnarray[7, rxnind]] * cs[rxnarray[8, rxnind]]
        end

        @fastmath flux = fderiv - rderiv
        _spreadreactantpartials!(jacp, flux, rxnarray, rxnind, Nspcs + rxnind)
        _spreadproductpartials!(jacp, -flux, rxnarray, rxnind, Nspcs + rxnind)

        @fastmath @inbounds gderiv = rderiv * (p[Nspcs+rxnind] * p[Nspcs+rxnind]) / kfs[rxnind] * RTinv

        @inbounds jacp[rxnarray[1, rxnind], rxnarray[1, rxnind]] -= gderiv
        @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[1, rxnind])
        if @inbounds rxnarray[2, rxnind] !== 0
            @inbounds jacp[rxnarray[2, rxnind], rxnarray[1, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[1, rxnind], rxnarray[2, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[2, rxnind], rxnarray[2, rxnind]] -= gderiv
            @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[2, rxnind])
            if @inbounds rxnarray[3, rxnind] !== 0
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[1, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[2, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[1, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[2, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3, rxnind], rxnarray[3, rxnind]] -= gderiv
                @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[3, rxnind])
                if @inbounds rxnarray[4, rxnind] !== 0
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[1, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[2, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[3, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[1, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[2, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[3, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[4, rxnind], rxnarray[4, rxnind]] -= gderiv
                    @inbounds _spreadreactantpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[4, rxnind])
                end
            end
        end

        @inbounds jacp[rxnarray[5, rxnind], rxnarray[5, rxnind]] -= gderiv
        @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[5, rxnind])
        if @inbounds rxnarray[6, rxnind] !== 0
            @inbounds jacp[rxnarray[6, rxnind], rxnarray[5, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[5, rxnind], rxnarray[6, rxnind]] -= gderiv
            @inbounds jacp[rxnarray[6, rxnind], rxnarray[6, rxnind]] -= gderiv
            @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[6, rxnind])
            if @inbounds rxnarray[7, rxnind] !== 0
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[5, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[6, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[5, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds jacp[rxnarray[7, rxnind], rxnarray[7, rxnind]] -= gderiv
                @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[7, rxnind])
                if @inbounds rxnarray[8, rxnind] !== 0
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[5, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[6, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[7, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[5, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[6, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[7, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds jacp[rxnarray[8, rxnind], rxnarray[8, rxnind]] -= gderiv
                    @inbounds _spreadproductpartials!(jacp, gderiv, rxnarray, rxnind, rxnarray[8, rxnind])
                end
            end
        end
    end
    jacp .*= V
end

function jacobianp!(jacp::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTPDomain,ParametrizedTPDomain}}
    jacp .= 0.0
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)

    Nspcs = length(cs)
    Nrxns = size(domain.rxnarray)[2]

    jacobianpnsderiv!(jacp, y, p, t, domain, domain.rxnarray, cs, T, V, kfs, krevs, Nspcs, Nrxns)

    @simd for i in 1:length(p)
        @views @inbounds @fastmath jacp[domain.indexes[3], i] = sum(jacp[domain.indexes[1]:domain.indexes[2], i]) * R * T / P
    end

    @simd for ind in domain.constantspeciesinds
        @inbounds jacp[ind, :] .= 0.0
    end
end

function jacobianp!(jacp::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantVDomain,ParametrizedVDomain}}
    jacp .= 0.0
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)

    Nspcs = length(cs)
    Nrxns = size(domain.rxnarray)[2]

    dydt = zeros(size(y))
    addreactionratecontributions!(dydt, domain.rxnarray, cs, kfs, krevs)
    dydt .*= V

    jacobianpnsderiv!(jacp, y, p, t, domain, domain.rxnarray, cs, T, V, kfs, krevs, Nspcs, Nrxns)

    @simd for i in 1:Nspcs
        @views @fastmath @inbounds jacp[domain.indexes[3], i] = -(dot(Us, jacp[domain.indexes[1]:domain.indexes[2], i]) + dydt[i]) / (N * Cvave)
        @views @fastmath @inbounds jacp[domain.indexes[4], i] = sum(jacp[domain.indexes[1]:domain.indexes[2], i]) * R * T / V + P / T * jacp[domain.indexes[3], i]
    end

    @simd for i in Nspcs+1:Nspcs+Nrxns
        @views @fastmath @inbounds jacp[domain.indexes[3], i] = -dot(Us, jacp[domain.indexes[1]:domain.indexes[2], i]) / (N * Cvave)
        @views @fastmath @inbounds jacp[domain.indexes[4], i] = sum(jacp[domain.indexes[1]:domain.indexes[2], i]) * R * T / V + P / T * jacp[domain.indexes[3], i]
    end

    @simd for ind in domain.constantspeciesinds
        @inbounds jacp[ind, :] .= 0.0
    end

    @simd for inter in interfaces
        if isa(inter, Inlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in 1:Nspcs
                ddGidTdt = flow * (-ns[i] / N) / (N * Cvave)
                jacp[domain.indexes[3], i] += ddGidTdt
                jacp[domain.indexes[4], i] += P / T * ddGidTdt
            end
        end
    end
end

function jacobianp!(jacp::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantPDomain,ParametrizedPDomain}}
    jacp .= 0.0
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)

    Nspcs = length(cs)
    Nrxns = size(domain.rxnarray)[2]

    dydt = zeros(size(y))
    addreactionratecontributions!(dydt, domain.rxnarray, cs, kfs, krevs)
    dydt .*= V

    jacobianpnsderiv!(jacp, y, p, t, domain, domain.rxnarray, cs, T, V, kfs, krevs, Nspcs, Nrxns)

    @fastmath Cpave = Cvave + R
    @simd for i in 1:Nspcs
        @views @fastmath @inbounds jacp[domain.indexes[3], i] = -(dot(Hs, jacp[domain.indexes[1]:domain.indexes[2], i]) + dydt[i]) / (N * Cpave) #divide by V to cancel ωV to ω
        @views @fastmath @inbounds jacp[domain.indexes[4], i] = sum(jacp[domain.indexes[1]:domain.indexes[2], i]) * R * T / P + jacp[domain.indexes[3], i] * V / T
    end

    @simd for i in Nspcs+1:Nspcs+Nrxns
        @views @fastmath @inbounds jacp[domain.indexes[3], i] = -(dot(Hs, jacp[domain.indexes[1]:domain.indexes[2], i])) / (N * Cpave) #divide by V to cancel ωV to ω
        @views @fastmath @inbounds jacp[domain.indexes[4], i] = sum(jacp[domain.indexes[1]:domain.indexes[2], i]) * R * T / P + jacp[domain.indexes[3], i] * V / T
    end

    @simd for ind in domain.constantspeciesinds
        @inbounds jacp[ind, :] .= 0.0
    end

    @simd for inter in interfaces
        if isa(inter, Inlet) && domain == inter.domain
            flow = inter.F(t)
            for i in 1:Nspcs
                ddGidTdt = flow * (-ns[i] / N) / (N * Cpave)
                jacp[domain.indexes[3], i] += ddGidTdt
                jacp[domain.indexes[4], i] += ddGidTdt * V / T
            end
        end
    end
end

function jacobianp!(jacp::Q, y::U, p::W, t::Z, domain::D, interfaces::Q3, colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTVDomain,ParametrizedTConstantVDomain,ConstantTAPhiDomain}}
    jacp .= 0.0
    ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR = calcthermo(domain, y, t, p)

    Nspcs = length(cs)
    Nrxns = size(domain.rxnarray)[2]

    jacobianpnsderiv!(jacp, y, p, t, domain, domain.rxnarray, cs, T, V, kfs, krevs, Nspcs, Nrxns)

    @simd for ind in domain.constantspeciesinds
        @inbounds jacp[ind, :] .= 0.0
    end
end

export jacobianp!

function getreactionindices(ig::Q) where {Q<:AbstractPhase}
    return deepcopy(ig.rxnarray)
end
export getreactionindices

@inline function getsensspcsrxns(domain::D, ind::Int64) where {D<:AbstractDomain}
    sensspcinds = Array{Int64,1}()
    sensrxninds = Array{Int64,1}()
    for rxnind in 1:size(domain.rxnarray)[2]
        if ind in @inbounds @view domain.rxnarray[:, rxnind]
            for spcind in @inbounds @view domain.rxnarray[:, rxnind]
                if !(spcind in sensspcinds) && (spcind !== 0)
                    push!(sensspcinds, spcind)
                end
            end
            push!(sensrxninds, rxnind)
        end
    end

    sensrxns = Array{ElementaryReaction,1}(undef, length(sensrxninds))
    sensspcs = Array{Species,1}(undef, length(sensspcinds))
    sensspcnames = Array{String,1}(undef, length(sensspcinds))
    senstooriginspcind = Array{Int64,1}(undef, length(sensspcinds))
    senstooriginrxnind = Array{Int64,1}(undef, length(sensrxninds))
    for (i, spcind) in enumerate(sensspcinds)
        spc = getphasespecies(domain.phase)[spcind]
        sensspcnames[i] = spc.name
        @inbounds sensspcs[i] = Species(
            name=spc.name,
            index=i,
            inchi=spc.inchi,
            smiles=spc.smiles,
            thermo=spc.thermo,
            atomnums=spc.atomnums,
            bondnum=spc.bondnum,
            diffusion=spc.diffusion,
            radius=spc.radius,
            radicalelectrons=spc.radicalelectrons,
            molecularweight=spc.molecularweight,
            henrylawconstant=spc.henrylawconstant,
            liquidvolumetricmasstransfercoefficient=spc.liquidvolumetricmasstransfercoefficient,
            comment=spc.comment,
            isfragment=spc.isfragment,
            isfragmentintermediate=spc.isfragmentintermediate
        )
        @inbounds senstooriginspcind[i] = spcind
    end

    for (i, rxnind) in enumerate(sensrxninds)
        rxn = domain.phase.reactions[rxnind]
        reactants = Array{Species,1}()
        reactantinds = Array{Int64,1}()
        @simd for reactant in rxn.reactants
            ind = findfirst(isequal(reactant.name), sensspcnames)
            @inbounds push!(reactants, sensspcs[ind])
            push!(reactantinds, ind)
        end
        products = Array{Species,1}()
        productinds = Array{Int64,1}()
        @simd for product in rxn.products
            ind = findfirst(isequal(product.name), sensspcnames)
            @inbounds push!(products, sensspcs[ind])
            push!(productinds, ind)
        end

        @inbounds sensrxns[i] = ElementaryReaction(
            index=i,
            reactants=SVector(reactants...),
            reactantinds=MVector(reactantinds...),
            products=SVector(products...),
            productinds=MVector(productinds...),
            kinetics=rxn.kinetics,
            radicalchange=rxn.radicalchange,
            pairs=rxn.pairs
        )
        @inbounds senstooriginrxnind[i] = rxnind
    end

    return sensspcs, sensrxns, sensspcnames, senstooriginspcind, senstooriginrxnind
end

@inline function getsensdomain(domain::D, ind::Int64) where {D<:AbstractDomain}

    sensspcs, sensrxns, sensspcnames, senstooriginspcind, senstooriginrxnind = getsensspcsrxns(domain, ind)

    initialconds = Dict{String,Float64}()

    for fn in fieldnames(typeof(domain))
        if fn in (:T, :P, :V)
            initialconds["$fn"] = getfield(domain, fn)
        end
    end

    d = split(repr(typeof(domain)), "{")[1]

    if occursin(".", d)
        d = Symbol(split(d, ".")[2])
    else
        d = Symbol(d)
    end

    if isa(domain.phase, IdealGas)
        return eval(d)(phase=IdealGas(sensspcs, sensrxns, name="phase"), initialconds=initialconds)[1], sensspcnames, senstooriginspcind, senstooriginrxnind
    else
        return eval(d)(phase=IdealDiluteSolution(sensspcs, sensrxns, domain.phase.solvent, name="phase"), initialconds=initialconds)[1], sensspcnames, senstooriginspcind, senstooriginrxnind
    end
end
