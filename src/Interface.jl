using LinearAlgebra
using SciMLBase
using LsqFit

abstract type AbstractInterface end
export AbstractInterface

abstract type AbstractBoundaryInterface <: AbstractInterface end
export AbstractBoundaryInterface

abstract type AbstractInternalInterface <: AbstractInterface end
export AbstractInternalInterface

abstract type AbstractReactiveInternalInterface <: AbstractInternalInterface end
export AbstractReactiveInternalInterface

struct EmptyInterface <: AbstractInterface end
export EmptyInterface

struct ReactiveInternalInterface{T,B,C,C2,N,Q<:AbstractReaction,X} <: AbstractReactiveInternalInterface
    domain1::T
    domain2::N
    reactions::Array{Q,1}
    veckinetics::X
    veckineticsinds::Array{Int64,1}
    rxnarray::B
    stoichmatrix::C
    Nrp1::C2
    Nrp2::C2
    A::Float64
    parameterindexes::Array{Int64,1}
    domaininds::Array{Int64,1}
    p::Array{Float64,1}
    reversibililty::Array{Bool,1}
    forwardability::Array{Bool,1}
end
function ReactiveInternalInterface(domain1,domain2,reactions,A)
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,products=rxn.products,
        productinds=rxn.productinds,kinetics=rxn.kinetics,electronchange=rxn.electronchange,radicalchange=rxn.radicalchange,reversible=rxn.reversible,forwardable=rxn.forwardable,pairs=rxn.pairs) for (i,rxn) in enumerate(rxns)]
    rxnarray = getinterfacereactioninds(domain1,domain2,rxns)
    M,Nrp1,Nrp2 = getstoichmatrix(domain1,domain2,reactions)
    reversibility = Array{Bool,1}(getfield.(rxns,:reversible))
    forwardability = Array{Bool,1}(getfield.(rxns,:forwardable))
    return ReactiveInternalInterface(domain1,domain2,
            rxns,vectuple,posinds,rxnarray,M,Nrp1,Nrp2,A,[1,length(reactions)],
            [0,1],ones(length(rxns)),reversibility,forwardability),ones(length(rxns))
end
export ReactiveInternalInterface

function getkfskrevs(ri::ReactiveInternalInterface,T1,T2,phi1,phi2,Gs1,Gs2,cstot::Array{Q,1}) where {Q}
    Gpart = ArrayPartition(Gs1,Gs2)
    dGrxn = -ri.stoichmatrix*Gpart
    kfs = getkfs(ri,T1,0.0,0.0,Array{Q,1}(),ri.A,phi1,dGrxns,0.0)
    Kc = getKc.(ri.reactions,ri.domain1.phase,ri.domain2.phase,Ref(Gs1),Ref(Gs2),dGrxns,T1,phi1)
    krevs = kfs./Kc
    return kfs,krevs
end

function evaluate(ri::ReactiveInternalInterface, dydt, domains, T1, T2, phi1, phi2, Gs1, Gs2, cstot, p::W) where {W<:SciMLBase.NullParameters}
    kfs, krevs = getkfskrevs(ri, T1, T2, phi1, phi2, Gs1, Gs2, cstot)
    addreactionratecontributions!(dydt, ri.rxnarray, cstot, kfs, krevs, ri.A)
end

function evaluate(ri::ReactiveInternalInterface, dydt, domains, T1, T2, phi1, phi2, Gs1, Gs2, cstot, p)
    kfs, krevs = getkfskrevs(ri, T1, T2, phi1, phi2, Gs1, Gs2, cstot)
    addreactionratecontributions!(dydt, ri.rxnarray, cstot, kfs .* p[ri.parameterindexes[1]:ri.parameterindexes[2]], krevs .* p[ri.parameterindexes[1]:ri.parameterindexes[2]], ri.A)
end
export evaluate


struct ReactiveInternalInterfaceConstantTPhi{J,N,B,B2,B3,C,C2,Q<:AbstractReaction} <: AbstractReactiveInternalInterface
    domain1::J
    domain2::N
    reactions::Array{Q,1}
    rxnarray::B
    stoichmatrix::C
    Nrp1::C2
    Nrp2::C2
    kfs::B2
    krevs::B3
    T::Float64
    A::Float64
    parameterindexes::Array{Int64,1}
    domaininds::Array{Int64,1}
    p::Array{Float64,1}
    reversibility::Array{Bool,1}
    forwardability::Array{Bool,1}
end
function ReactiveInternalInterfaceConstantTPhi(domain1,domain2,reactions,T,A,phi=0.0)
    @assert domain1.T == domain2.T 
    reactions = upgradekinetics(reactions,domain1,domain2)
    rxnarray = getinterfacereactioninds(domain1,domain2,reactions)
    M,Nrp1,Nrp2 = getstoichmatrix(domain1,domain2,reactions)
    Gpart = ArrayPartition(domain1.Gs,domain2.Gs)
    dGrxns = -M*Gpart
    kfs = getkf.(reactions,nothing,T,0.0,0.0,Ref([]),A,phi,dGrxns,0.0)
    Kc = getKcs(domain1.phase,domain2.phase,T,Nrp1,Nrp2,dGrxns)
    krevs = kfs./Kc
    reversibility = Array{Bool,1}(getfield.(reactions,:reversible))
    forwardability = Array{Bool,1}(getfield.(reactions,:forwardable))
    if isa(reactions,Vector{Any})
        reactions = convert(Vector{ElementaryReaction},reactions)
    end
    if isa(kfs, Vector{Any})
        kfs = convert(Vector{Float64}, kfs)
    end
    return ReactiveInternalInterfaceConstantTPhi(domain1, domain2, reactions,
        rxnarray, M, Nrp1, Nrp2, kfs, krevs, T, A, [1, length(reactions)],
        [0, 1], kfs[1:end], reversibility, forwardability), kfs[1:end]
end
export ReactiveInternalInterfaceConstantTPhi

function getkfskrevs(ri::ReactiveInternalInterfaceConstantTPhi, T1, T2, phi1, phi2, Gs1, Gs2, cstot)
    return ri.kfs, ri.krevs
end

function evaluate(ri::ReactiveInternalInterfaceConstantTPhi, dydt, domains, T1, T2, phi1, phi2, Gs1, Gs2, cstot, p::W) where {W<:SciMLBase.NullParameters}
    addreactionratecontributions!(dydt, ri.rxnarray, cstot, ri.kfs, ri.krevs, ri.A)
end

function evaluate(ri::ReactiveInternalInterfaceConstantTPhi, dydt, domains, T1, T2, phi1, phi2, Gs1, Gs2, cstot, p)
    if p[ri.parameterindexes[1]:ri.parameterindexes[2]] == ri.kfs
        kfs = ri.kfs
    else
        kfs = p[ri.parameterindexes[1]:ri.parameterindexes[2]]
    end
    if length(Gs1) == 0 || length(Gs2) == 0 || (all(Gs1 .== ri.domain1.Gs) && all(Gs2 .== ri.domain2.Gs))
        krevs = ri.krevs
    else
        Kc = getKcs(ri, T1, Gs1, Gs2)
        krevs = kfs ./ Kc
    end
    addreactionratecontributions!(dydt, ri.rxnarray, cstot, kfs, krevs, ri.A)
end
export evaluate

"""
construct the stochiometric matrix for the reactions crossing both domains and the reaction molecule # change
"""
function getstoichmatrix(domain1, domain2, rxns)
    M = spzeros(length(rxns), length(domain1.phase.species) + length(domain2.phase.species))
    Nrp1 = zeros(length(rxns))
    Nrp2 = zeros(length(rxns))
    N1 = length(domain1.phase.species)
    spcs1 = domain1.phase.species
    spcs2 = domain2.phase.species
    for (i, rxn) in enumerate(rxns)
        Nrp1[i] = Float64(length([x for x in rxn.products if x in spcs1]) - length([x for x in rxn.reactants if x in spcs1]))
        Nrp2[i] = Float64(length([x for x in rxn.products if x in spcs2]) - length([x for x in rxn.reactants if x in spcs2]))
        for (j, r) in enumerate(rxn.reactants)
            isfirst = true
            ind = findfirst(isequal(r), domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r), domain2.phase.species)
            end
            M[i, isfirst ? ind : ind + N1] += 1
        end
        for (j, r) in enumerate(rxn.products)
            isfirst = true
            ind = findfirst(isequal(r), domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r), domain2.phase.species)
            end
            M[i, isfirst ? ind : ind + N1] -= 1
        end
    end
    return M, Nrp1, Nrp2
end

function getinterfacereactioninds(domain1, domain2, reactions)
    indices = zeros(Int64, (8, length(reactions)))
    N1 = length(domain1.phase.species)
    for (i, rxn) in enumerate(reactions)
        for (j, r) in enumerate(rxn.reactants)
            isfirst = true
            ind = findfirst(isequal(r), domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r), domain2.phase.species)
            end
            indices[j, i] = isfirst ? ind : ind + N1
        end
        for (j, r) in enumerate(rxn.products)
            isfirst = true
            ind = findfirst(isequal(r), domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r), domain2.phase.species)
            end
            indices[j+4, i] = isfirst ? ind : ind + N1
        end
    end
    return indices
end

function upgradekinetics(rxns, domain1, domain2)
    domain1surf = hasproperty(domain1.phase, :sitedensity)
    domain2surf = hasproperty(domain2.phase, :sitedensity)
    @assert !(domain1surf && domain2surf)
    if domain1surf
        surfdomain = domain1
    elseif domain2surf
        surfdomain = domain2
    end
    newrxns = Array{ElementaryReaction,1}(undef, length(rxns))
    for (i, rxn) in enumerate(rxns)
        if isa(rxn.kinetics, StickingCoefficient)
            spc = [spc for spc in rxn.reactants if !(spc in surfdomain.phase.species)]
            @assert length(spc) == 1
            kin = stickingcoefficient2arrhenius(rxn.kinetics, surfdomain.phase.sitedensity, length(rxn.reactants) - 1, spc[1].molecularweight)
            newrxns[i] = ElementaryReaction(index=rxn.index, reactants=rxn.reactants, reactantinds=rxn.reactantinds, products=rxn.products,
                productinds=rxn.productinds, kinetics=kin, electronchange=rxn.electronchange, radicalchange=rxn.radicalchange, reversible=rxn.reversible, forwardable=rxn.forwardable, pairs=rxn.pairs)
        else
            newrxns[i] = rxn
        end
    end
    return [rxn for rxn in newrxns]
end

function stickingcoefficient2arrhenius(sc, sitedensity, N, mw; Tmin=300.0, Tmax=2000.0)
    mass = mw / Na
    ksc(T) = sc(T) / sitedensity^N * sqrt(kB * T / (2.0 * pi * mass))
    Ts = Array{Float64,1}(Tmin:10:Tmax)
    kscvals = ksc.(Ts)
    k(T, p) = log.(abs(p[1]) .* T .^ p[2] .* exp.(-p[3] ./ (R .* T)))
    p0 = [sc.A / sitedensity * sqrt(kB * 1000.0 / (2.0 * pi * mass)), 0.5, sc.Ea]
    fit = curve_fit(k, Ts, log.(kscvals), p0; x_tol=1e-18)
    @assert fit.converged
    p = fit.param
    p[1] = abs(p[1])
    return Arrhenius(; A=p[1], n=p[2], Ea=p[3])
end

struct Inlet{Q<:Real,S,V<:AbstractArray,U<:Real,X<:Real,FF<:Function} <: AbstractBoundaryInterface
    domain::S
    y::V
    F::FF
    T::U
    P::X
    H::Q
end

function Inlet(domain::V, conddict::Dict{X1,X}, F::FF) where {V,X1,X,B<:Real,FF<:Function}
    T = 0.0
    P = 0.0

    y = makespcsvector(domain.phase, conddict)
    yout = y ./ sum(y)

    if haskey(conddict, "T")
        T = conddict["T"]
    end
    if haskey(conddict, "P")
        P = conddict["P"]
    end

    if haskey(conddict, "Hin")
        H = conddict["Hin"]
    else
        @assert T != 0.0 && P != 0.0 "T and P need to be provided if Hin is not provided for Inlet"
        H = dot(getEnthalpy.(getfield.(domain.phase.species, :thermo), T), yout)
    end
    return Inlet(domain, yout, F, T, P, H)
end

export Inlet

struct Outlet{V,FF<:Function} <: AbstractBoundaryInterface
    domain::V
    F::FF
end
export Outlet

"""
kLAkHCondensationEvaporationWithReservoir adds evaporation and condensation to
(1) a liquid phase domain with a constant composition vapor resevoir, where number of moles, P, and T need to be specified, or
(2) a gas phase domain with a constant composition liquid resevoir, where number of moles, V, and T need to be specified.
kLA and kH are used to model cond/evap. 
kLA is liquid volumetric mass transfer coefficient with unit 1/s , and kH is Henry's law constant defined as gas phase partial pressure of solute over liquid phase concentration of solute.
"""

struct kLAkHCondensationEvaporationWithReservoir{S,V1<:AbstractArray,V2<:Real,V3<:Real,V4<:Real,V5<:AbstractArray,V6<:Real,V7<:AbstractArray,V8<:AbstractArray} <: AbstractBoundaryInterface
    domain::S
    molefractions::V1
    T::V2
    P::V3
    V::V4
    cs::V5
    H::V6
    kLAs::V7
    kHs::V8
end

function kLAkHCondensationEvaporationWithReservoir(domain::D, conddict::Dict{X1,X}) where {D,X1,X}
    y = makespcsvector(domain.phase, conddict)
    molefractions = y ./ sum(y)
    T = conddict["T"]
    H = dot(getEnthalpy.(getfield.(domain.phase.species, :thermo), T), molefractions)
    kLAs = [T -> kLA(T=T) for kLA in getfield.(domain.phase.species, :liquidvolumetricmasstransfercoefficient)]
    kHs = [T -> kH(T=T) for kH in getfield.(domain.phase.species, :henrylawconstant)]

    if isa(domain.phase, IdealDiluteSolution)
        if !haskey(conddict, "P")
            @error "P needs to be specified for the vapor resevoir over the liquid phase domain"
        end
        P = conddict["P"]
        V = sum(y) * R * T / P
        return kLAkHCondensationEvaporationWithReservoir(domain, molefractions, T, P, V, Array{Float64,1}(), H, kLAs, kHs)
    elseif isa(domain.phase, IdealGas)
        if !haskey(conddict, "V")
            @error "V needs to be specified for the liquid resevoir under the gas phase domain"
        end
        V = conddict["V"]
        cs = y ./ V
        return kLAkHCondensationEvaporationWithReservoir(domain, Array{Float64,1}(), T, 1e8, V, cs, H, kLAs, kHs)
    end
end

export kLAkHCondensationEvaporationWithReservoir

struct VolumetricFlowRateOutlet{V,F1<:Function} <: AbstractBoundaryInterface
    domain::V
    Vout::F1
end
export VolumetricFlowRateOutlet

struct VolumetricFlowRateInlet{Q<:Real,S,V<:AbstractArray,U<:Real,X<:Real,FF<:Function} <: AbstractBoundaryInterface
    domain::S
    cs::V
    Vin::FF
    T::U
    P::X
    Hpervolume::Q
end

"""
VolumetricFlowRateInlet adds a volumetric flow rate inlet to a domain.
conddict should provide the species composition in concentration, T, and P.
Hin is the inlet enthalpy per mol (J/mol) and is optional, if not provided, it will be calculated from T and P.
"""
function VolumetricFlowRateInlet(domain::V, conddict::Dict{X1,X}, Vin::FF) where {V,X1,X,FF<:Function}
    T = 0.0
    P = 0.0

    cs = makespcsvector(domain.phase, conddict)

    if haskey(conddict, "T")
        T = conddict["T"]
    end
    if haskey(conddict, "P")
        P = conddict["P"]
    end

    if haskey(conddict, "Hin")
        Hpervolume = conddict["Hin"] / sum(cs) # convert to J/m3
    else
        @assert T != 0.0 && P != 0.0 "T and P need to be provided if Hin is not provided for Inlet"
        Hpervolume = dot(getEnthalpy.(getfield.(domain.phase.species, :thermo), T), cs) # J/m3
    end
    return VolumetricFlowRateInlet(domain, cs, Vin, T, P, Hpervolume)
end

export VolumetricFlowRateInlet

"""
VolumeMaintainingOutlet is designed for gas phase domain such that the flow rate of this outlet will adjust to maintain the volume of the 
    domain to be constant. This is particularly useful to simulate any vapor-liquid phase system where the gas phase outlet
    is determined by the amount of evaporation.
"""
struct VolumeMaintainingOutlet{V} <: AbstractBoundaryInterface
    domain::V
end

export VolumeMaintainingOutlet

struct VaporLiquidMassTransferInternalInterfaceConstantT{D1,D2,B} <: AbstractInternalInterface
    domaingas::D1
    domainliq::D2
    ignoremasstransferspcnames::Array{String,1}
    ignoremasstransferspcinds::B
    kLAs::Array{Float64,1}
    kHs::Array{Float64,1}
    parameterindexes::Array{Int64,1}
    domaininds::Array{Int64,1}
    p::Array{Float64,1}
end

function VaporLiquidMassTransferInternalInterfaceConstantT(domaingas, domainliq, ignoremasstransferspcnames)
    @assert isa(domaingas.phase, IdealGas)
    @assert isa(domainliq.phase, IdealDiluteSolution)
    @assert getfield.(domaingas.phase.species, :name) == getfield.(domainliq.phase.species, :name)
    T = domainliq.T
    phase = domainliq.phase
    ignoremasstransferspcinds = getinterfaceignoremasstransferspcinds(domaingas, domainliq, ignoremasstransferspcnames)
    kLAs = [kLA(T=T) for kLA in getfield.(phase.species, :liquidvolumetricmasstransfercoefficient)]
    kHs = [kH(T=T) for kH in getfield.(phase.species, :henrylawconstant)]
    return VaporLiquidMassTransferInternalInterfaceConstantT(domaingas, domainliq, ignoremasstransferspcnames, ignoremasstransferspcinds, kLAs, kHs, [1, length(domainliq.phase.species)], [0, 0], ones(length(domainliq.phase.species))), ones(length(domainliq.phase.species))
end
export VaporLiquidMassTransferInternalInterfaceConstantT

function getkLAkHs(vl::VaporLiquidMassTransferInternalInterfaceConstantT, Tgas, Tliq)
    return vl.kLAs, vl.kHs
end

function evaluate(vl::VaporLiquidMassTransferInternalInterfaceConstantT, dydt, Vgas, Vliq, Tgas, Tliq, N1, N2, P1, P2, Cvave1, Cvave2, ns1, ns2, Us1, Us2, cstot, p)
    kLAs, kHs = getkLAkHs(vl, Tgas, Tliq)
    @views @inbounds @fastmath evap = kLAs * Vliq .* cstot[vl.domainliq.indexes[1]:vl.domainliq.indexes[2]] #evap_i = kLA_i * Vliq * cliq_i
    @views @inbounds @fastmath cond = kLAs * Vliq .* cstot[vl.domaingas.indexes[1]:vl.domaingas.indexes[2]] * R * Tgas ./ kHs #cond_i = kLA_i * Vliq * Pgas_i / kH_i, Pgas_i = cgas_i * R * Tgas
    netevap = (evap .- cond)
    @simd for ind in vl.ignoremasstransferspcinds
        @inbounds netevap[ind] = 0.0
    end
    @views @inbounds @fastmath dydt[vl.domaingas.indexes[1]:vl.domaingas.indexes[2]] .+= netevap
    @views @inbounds @fastmath dydt[vl.domainliq.indexes[1]:vl.domainliq.indexes[2]] .-= netevap
end
export evaluate

function getinterfaceignoremasstransferspcinds(domaingas, domainliq, ignoremasstransferspcnames)
    indices = zeros(Int64, length(ignoremasstransferspcnames))
    spcnamesliq = getfield.(domainliq.phase.species, :name)
    for (i, name) in enumerate(ignoremasstransferspcnames)
        indliq = findfirst(isequal(name), spcnamesliq)
        indices[i] = domainliq.indexes[1] - 1 + indliq
    end
    return indices
end

struct FragmentBasedReactiveFilmGrowthInterfaceConstantT{D1,D2,Q<:AbstractReaction,M1} <: AbstractReactiveInternalInterface
    domainfilm::D1
    domain2::D2
    reactions::Array{Q,1}
    rxnarray::Array{Int64,2}
    fragmentbasedrxnarray::Array{Int64,2}
    stoichmatrix::M1
    Nrp1::Array{Float64,1}
    Nrp2::Array{Float64,1}
    kfs::Array{Float64,1}
    krevs::Array{Float64,1}
    T::Float64
    parameterindexes::Array{Int64,1}
    domaininds::Array{Int64,1}
    p::Array{Float64,1}
    reversibililty::Array{Bool,1}
    forwardability::Array{Bool,1}
    Mws::Array{Float64,1}
end

function FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domain2, reactions)
    @assert isa(domainfilm.phase, FragmentBasedIdealFilm)
    @assert domainfilm.T == domain2.T
    T = domainfilm.T

    rxnarray, fragmentbasedrxnarray = getfragmentbasedinterfacereactioninds(domainfilm, domain2, reactions)

    kfs = getkf.(reactions, nothing, T, 0.0, 0.0, Ref([]), 0.0, 0.0)
    Kc = getKc.(reactions, domainfilm.phase, domain2.phase, Ref(domainfilm.Gs), Ref(domain2.Gs), T, 0.0)
    krevs = kfs ./ Kc

    M, Nrp1, Nrp2 = getstoichmatrix(domainfilm, domain2, reactions)
    reversibility = Array{Bool,1}(getfield.(reactions, :reversible))
    forwardability = Array{Bool,1}(getfield.(reactions, :forwardable))

    Mws = vcat(getfield.(domainfilm.phase.fragments, :molecularweight), getfield.(domain2.phase.species, :molecularweight))
    return FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm, domain2, reactions, rxnarray, fragmentbasedrxnarray, M, Nrp1, Nrp2, kfs, krevs, T, [1, length(reactions)],
        [0, 1], kfs[1:end], reversibility, forwardability, Mws), kfs[1:end]
end
export FragmentBasedReactiveFilmGrowthInterfaceConstantT

function getfragmentbasedinterfacereactioninds(domainfilm, domain2, reactions)
    # find maximum number of species in fragment-based reaction and names of fragment-based species
    maxnumfragmentbasedreacprod = 0
    for rxn in reactions
        numfragmentbasedreac = length(rxn.fragmentbasedreactants)
        numfragmentbasedprod = length(rxn.fragmentbasedproducts)
        maxnumfragmentbasedreacprod = max(maxnumfragmentbasedreacprod, numfragmentbasedreac, numfragmentbasedprod)
    end
    fragmentbasedrxnarray = zeros(Int64, (maxnumfragmentbasedreacprod * 2, length(reactions)))
    rxnarray = zeros(Int64, (8, length(reactions)))

    N1 = length(domainfilm.phase.fragments)
    for (i, rxn) in enumerate(reactions)
        for (j, r) in enumerate(rxn.reactants)
            if !r.isfragmentintermediate
                isfirst = true
                ind = findfirst(isequal(r), domainfilm.phase.fragments)
                if ind === nothing
                    isfirst = false
                    ind = findfirst(isequal(r), domain2.phase.species)
                end
                rxnarray[j, i] = isfirst ? ind : ind + N1
            end
        end
        for (j, r) in enumerate(rxn.fragmentbasedreactants)
            isfirst = true
            ind = findfirst(isequal(r), domainfilm.phase.fragments)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r), domain2.phase.species)
            end
            fragmentbasedrxnarray[j, i] = isfirst ? ind : ind + N1
        end
        for (j, r) in enumerate(rxn.products)
            if !r.isfragmentintermediate
                isfirst = true
                ind = findfirst(isequal(r), domainfilm.phase.fragments)
                if ind === nothing
                    isfirst = false
                    ind = findfirst(isequal(r), domain2.phase.species)
                end
                rxnarray[j+4, i] = isfirst ? ind : ind + N1
            end
        end
        for (j, r) in enumerate(rxn.fragmentbasedproducts)
            isfirst = true
            ind = findfirst(isequal(r), domainfilm.phase.fragments)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r), domain2.phase.species)
            end
            fragmentbasedrxnarray[j+maxnumfragmentbasedreacprod, i] = isfirst ? ind : ind + N1
        end
    end
    return rxnarray, fragmentbasedrxnarray
end

export getfragmentbasedinterfacereactioninds

function getkfskrevs(ri::FragmentBasedReactiveFilmGrowthInterfaceConstantT)
    return ri.kfs, ri.krevs
end

function evaluate(ri::FragmentBasedReactiveFilmGrowthInterfaceConstantT, dydt, Vfilm, cstot)
    kfs, krevs = getkfskrevs(ri)
    addreactionratecontributions!(dydt, ri.fragmentbasedrxnarray, ri.rxnarray, cstot, kfs, krevs, Vfilm, ri.domainfilm.indexes[3], ri.Mws, ri.domainfilm.indexes[1]:ri.domainfilm.indexes[2])
    if hasproperty(ri.domain2, :epsilon)
        epsilon = ri.domain2.epsilon
        dydt[ri.domain2.thermovariabledict["V"]] = dydt[ri.domainfilm.indexes[3]] / ri.domainfilm.rho / (1 - epsilon) * epsilon
    end
end

export evaluate