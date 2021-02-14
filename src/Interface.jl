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

function upgradekinetics(rxns,domain1,domain2)
    domain1surf = hasproperty(domain1.phase,:sitedensity)
    domain2surf = hasproperty(domain2.phase,:sitedensity)
    @assert !(domain1surf && domain2surf)
    if domain1surf
        surfdomain = domain1
    elseif domain2surf
        surfdomain = domain2
    end
    for (i,rxn) in enumerate(rxns)
        if isa(rxn.kinetics,StickingCoefficient)
            spc = [spc for spc in rxn.reactants if !(spc in surfdomain.phase.species)]
            @assert length(spc) == 1
            kin = stickingcoefficient2arrhenius(rxn.kinetics,surfdomain.phase.sitedensity,length(rxn.reactants)-1,spc[1].molecularweight)
            rxns[i] = ElementaryReaction(index=rxn.index,reactants=rxn.reactants,reactantinds=rxn.reactantinds,products=rxn.products,
                productinds=rxn.productinds,kinetics=kin,radicalchange=rxn.radicalchange,reversible=rxn.reversible,pairs=rxn.pairs)
        end
    end
    return rxns
end

function stickingcoefficient2arrhenius(sc,sitedensity,N,mw;Tmin=300.0,Tmax=2000.0)
    mass = mw/Na
    ksc(T) = sc(T)/sitedensity^N*sqrt(kB*T/(2.0*pi*mass))
    Ts = Array{Float64,1}(Tmin:10:Tmax);
    kscvals = ksc.(Ts)
    k(T,p) = log.(abs(p[1]).*T.^p[2].*exp.(-p[3]./(R.*T)))
    p0 = [sc.A/sitedensity*sqrt(kB*1000.0/(2.0*pi*mass)),0.5,sc.Ea]
    fit = curve_fit(k,Ts,log.(kscvals),p0;x_tol=1e-18)
    @assert fit.converged
    p = fit.param
    p[1] = abs(p[1])
    return Arrhenius(;A=p[1],n=p[2],Ea=p[3])
end

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
