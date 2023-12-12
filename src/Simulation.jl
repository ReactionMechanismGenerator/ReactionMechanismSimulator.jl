using SciMLBase
import SciMLBase: AbstractODESolution, HermiteInterpolation
using DiffEqSensitivity
using ForwardDiff
using OrdinaryDiffEq.PreallocationTools

abstract type AbstractSimulation end
export AbstractSimulation

struct Simulation{Q<:AbstractODESolution,W<:AbstractDomain,M,L<:AbstractArray,G<:Function,G2<:AbstractArray,G3,G4,G5} <: AbstractSimulation
    sol::Q
    domain::W
    interfaces::M
    names::L
    N::G
    Ns::G2
    species::G3
    reactions::G4
    p::G5
end

function Simulation(sol::Q, domain::W, interfaces=[], p=nothing) where {Q<:AbstractODESolution,W<:AbstractDomain}
    names = getfield.(getphasespecies(domain.phase), :name)
    Ns = sum(hcat(sol.interp.u...)[domain.indexes[1]:domain.indexes[2], :], dims=1)
    if hasproperty(sol.interp, :du)
        Nderivs = sum(hcat(sol.interp.du...)[domain.indexes[1]:domain.indexes[2], :], dims=1)
    else
        Nderivs = sum(hcat([sol(t, Val{1}) for t in sol.t]...)[domain.indexes[1]:domain.indexes[2], :], dims=1)
    end
    N = HermiteInterpolation(sol.interp.t, Ns, Nderivs)
    F(t::T) where {T<:Real} = N(t, nothing, Val{0}, sol.prob.p, :left)
    if p === nothing
        p = domain.p
    end
    return Simulation(sol, domain, interfaces, names, F, Ns, getphasespecies(domain.phase), domain.phase.reactions, p)
end

function Simulation(sol::Q, domain::W, reducedmodelmappings::ReducedModelMappings, interfaces=[], p=nothing) where {Q<:AbstractODESolution,W<:AbstractDomain}

    function unlumpsol(t::tt, sol::Q, domain::W, reducedmodelmappings::ReducedModelMappings, reducedmodelcache::ReducedModelCache, interfaces=[], p=nothing) where {tt<:Real,Q<:AbstractODESolution,W<:AbstractDomain}

        yunlumped = zeros(length(sol(0)) + length(reducedmodelmappings.qssindexes) - length(reducedmodelmappings.lumpedgroupmapping) + length(reducedmodelmappings.lumpedindexes))
        y = sol(t)

        qssc = get_tmp(reducedmodelcache.qssc, first(y) * t) .= 0.0

        @inbounds @views yunlumped[reducedmodelmappings.reducedindexes] .= y[1:end-length(domain.thermovariabledict)-length(reducedmodelmappings.lumpedgroupmapping)]
        for (i, group) in enumerate(reducedmodelmappings.lumpedgroupmapping)
            for (index, weight) in group
                @inbounds yunlumped[index] = weight * y[length(reducedmodelmappings.reducedindexes)+i]
            end
        end
        @inbounds @views yunlumped[end-length(domain.thermovariabledict)+1:end] .= y[end-length(domain.thermovariabledict)+1:end]

        ns, cs, T, P, V, C, N, mu, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR, phi = calcthermo(domain, yunlumped, t, p)

        reducedmodelmappings.qssc!(qssc, cs, kfs, krevs)
        @inbounds yunlumped[reducedmodelmappings.qssindexes] .= qssc .* V
        return yunlumped
    end

    qssc = dualcache(zeros(length(reducedmodelmappings.qssindexes)))
    reducedmodelcache = ReducedModelCache(nothing, nothing, qssc)

    unlumpsol(t::T) where {T<:Real} = unlumpsol(t, sol, domain, reducedmodelmappings, reducedmodelcache, interfaces, p)

    u = [unlumpsol(t) for t in sol.t]
    t = sol.t
    sol = SciMLBase.build_solution(sol.prob, sol.alg, t, u, retcode=:Success)

    species_range = 1:length(getphasespecies(domain.phase))
    names = getfield.(getphasespecies(domain.phase), :name)
    Ns = sum(hcat(sol.interp.u...)[species_range, :], dims=1)
    if hasproperty(sol.interp, :du)
        Nderivs = sum(hcat(sol.interp.du...)[species_range, :], dims=1)
    else
        Nderivs = sum(hcat([sol(t, Val{1}) for t in sol.t]...)[species_range, :], dims=1)
    end
    N = HermiteInterpolation(sol.interp.t, Ns, Nderivs)
    F(t::T) where {T<:Real} = N(t, nothing, Val{0}, sol.prob.p, :left)
    if p === nothing
        p = domain.p
    end
    return Simulation(sol, domain, interfaces, names, F, Ns, getphasespecies(domain.phase), domain.phase.reactions, p)
end

export Simulation

struct SystemSimulation{Q,B<:AbstractODESolution,X,Y,Z}
    sol::B
    sims::Q
    interfaces::X
    names::Array{String,1}
    species::Array{Z,1}
    reactions::Array{Y,1}
    p::Array{Float64,1}
end

function getinterfacesfordomain(domain, interfaces)
    return [inter for inter in interfaces if interofdomain(inter, domain)]
end

function interofdomain(inter, domain)
    if hasfield(typeof(inter), :domain)
        return domain == inter.domain
    elseif hasfield(typeof(inter), :domain1)
        return domain == inter.domain1 || domain == inter.domain2
    elseif hasfield(typeof(inter), :domaingas)
        return domain == inter.domaingas || domain == inter.domainliq
    elseif hasfield(typeof(inter), :domainfilm)
        return domain == inter.domainfilm || domain == inter.domain2
    end
end

function SystemSimulation(sol, domains, interfaces, p)
    sims = Tuple([Simulation(sol, domain, getinterfacesfordomain(domain, interfaces), p) for domain in domains])
    names = Array{String,1}()
    reactions = Array{ElementaryReaction,1}()
    species = Array{Species,1}()
    for sim in sims
        append!(names, sim.names)
        append!(species, getphasespecies(sim.domain.phase))
        append!(reactions, sim.domain.phase.reactions)
    end
    for inter in interfaces
        if hasproperty(inter, :reactions)
            append!(reactions, inter.reactions)
        end
    end
    return SystemSimulation(sol, sims, interfaces, names, species, reactions, p)
end
export SystemSimulation

length(p::T) where {T<:AbstractSimulation} = 1
export length

iterate(p::T) where {T<:AbstractSimulation} = p
export iterate

Broadcast.broadcastable(p::T) where {T<:AbstractSimulation} = Ref(p)
export broadcastable

spcindex(bsol::Z, name::Q) where {Z<:Simulation,Q<:AbstractString} = findfirst(isequal(name), getfield.(getphasespecies(bsol.domain.phase), :name))
export spcindex

function molefractions(bsol::Q, name::W, t::E) where {Q<:AbstractSimulation,W<:String,E<:Real}
    @assert name in bsol.names
    ind = findfirst(isequal(name), bsol.names) + bsol.domain.indexes[1] - 1
    return bsol.sol(t)[ind] / bsol.N(t)
end

function molefractions(bsol::Q, t::E) where {Q<:AbstractSimulation,E<:Real}
    return bsol.sol(t)[bsol.domain.indexes[1]:bsol.domain.indexes[2]] ./ bsol.N(t)
end

function molefractions(bsol::Q) where {Q<:AbstractSimulation}
    @views return hcat(bsol.sol.u...)[bsol.domain.indexes[1]:bsol.domain.indexes[2], :] ./ bsol.Ns
end

export molefractions

function concentrations(bsol::Q, name::W, t::E) where {Q<:AbstractSimulation,W<:String,E<:Real}
    @assert name in bsol.names
    ind = findfirst(isequal(name), bsol.names) + bsol.domain.indexes[1] - 1
    return bsol.sol(t)[ind] / getdomainsize(bsol, t)
end

function concentrations(bsol::Q, t::E) where {Q<:AbstractSimulation,E<:Real}
    return bsol.sol(t)[bsol.domain.indexes[1]:bsol.domain.indexes[2]] ./ getdomainsize(bsol, t)
end

function concentrations(bsol::Q) where {Q<:AbstractSimulation}
    @views return hcat(bsol.sol.u...)[bsol.domain.indexes[1]:bsol.domain.indexes[2], :] ./ getdomainsize.(bsol, bsol.sol.t)'
end

function concentrations(ssys::Q, name::W, t::E) where {Q<:SystemSimulation,W<:String,E<:Real}
    @assert name in ssys.names
    for sim in ssys.sims
        if name in sim.names
            return concentrations(sim, name, t)
        end
    end
end

function concentrations(ssys::Q, t::E) where {Q<:SystemSimulation,E<:Real}
    cstot = zeros(length(ssys.species))
    for sim in ssys.sims
        cstot[sim.domain.indexes[1]:sim.domain.indexes[2]] .= concentrations(sim, t)
    end
    return cstot
end

function concentrations(ssys::Q) where {Q<:SystemSimulation}
    cstots = zeros((length(ssys.species), length(ssys.sol.t)))
    for (i, t) in enumerate(ssys.sol.t)
        for sim in ssys.sims
            cstots[sim.domain.indexes[1]:sim.domain.indexes[2], i] .= concentrations(sim, t)
        end
    end
    return cstots
end

export concentrations

getT(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantTPDomain,ConstantTVDomain,FragmentBasedConstantTrhoDomain},K<:Real,Q,G,L} = bsol.domain.T
getT(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ParametrizedVDomain,ConstantPDomain,ParametrizedPDomain},K<:Real,Q,G,L} = bsol.sol(t)[bsol.domain.indexes[3]]
getT(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ParametrizedTConstantVDomain,ParametrizedTPDomain},K<:Real,Q,G,L} = bsol.domain.T(t)
export getT
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain,ParametrizedTConstantVDomain},K<:Real,Q,G,L} = bsol.domain.V
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedVDomain,K<:Real,Q,G,L} = bsol.domain.V(t)
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantTPDomain,ParametrizedTPDomain},K<:Real,Q,G,L} = bsol.sol(t)[bsol.domain.indexes[3]]
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ParametrizedPDomain,ConstantPDomain},K<:Real,Q,G,L} = bsol.sol(t)[bsol.domain.indexes[4]]
getV(bsol::Simulation{Q,W,L,G}, t::K) where {W<:FragmentBasedConstantTrhoDomain,K<:Real,Q,G,L} = bsol.sol(t)[bsol.domain.indexes[3]] / bsol.domain.rho
export getV
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantTPDomain,ConstantPDomain},K<:Real,Q,G,L} = bsol.domain.P
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantTVDomain,ParametrizedTConstantVDomain,FragmentBasedConstantTrhoDomain},K<:Real,Q,G,L} = 1.0e6
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ParametrizedTPDomain,ParametrizedPDomain},K<:Real,Q,G,L} = bsol.domain.P(t)
getP(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ParametrizedVDomain},K<:Real,Q,G,L} = bsol.sol(t)[bsol.domain.indexes[4]]
export getP
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantTPDomain,K<:Real,Q,G,L} = bsol.domain.P / (R * bsol.domain.T)
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain,ParametrizedTConstantVDomain},K<:Real,Q,G,L} = bsol.N(t) / bsol.domain.V
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:FragmentBasedConstantTrhoDomain,K<:Real,Q,G,L} = bsol.N(t) / getV(bsol, t)
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedVDomain,K<:Real,Q,G,L} = bsol.N(t) / bsol.domain.V(t)
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedTPDomain,K<:Real,Q,G,L} = bsol.domain.P(t) / (R * bsol.domain.T(t))
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ConstantPDomain,K<:Real,Q,G,L} = bsol.domain.P / (R * getT(bsol, t))
getC(bsol::Simulation{Q,W,L,G}, t::K) where {W<:ParametrizedPDomain,K<:Real,Q,G,L} = bsol.domain.P(t) / (R * getT(bsol, t))
export getC
getdomainsize(bsol, t) = getV(bsol, t)
getdomainsize(bsol::Simulation{Q,W,L,G}, t) where {W<:ConstantTAPhiDomain,Q,L,G} = bsol.domain.A

"""
calculates the rates of production/loss at a given time point
this outputs a sparse matrix of  num reactions xnum species containing the production/loss
rate of that species associated with that reaction
"""
function rops(bsol::Q, t::X) where {Q<:Simulation,X<:Real}
    ropmat = spzeros(length(bsol.domain.phase.reactions), length(getphasespecies(bsol.domain.phase)))
    cs, kfs, krevs = calcthermo(bsol.domain, bsol.sol(t), t)[[2, 9, 10]]
    V = getdomainsize(bsol, t)
    @simd for i in 1:length(bsol.domain.phase.reactions)
        rxn = bsol.domain.phase.reactions[i]
        R = getrate(rxn, cs, kfs, krevs) * V
        for ind in rxn.productinds
            ropmat[i, ind] += R
        end
        for ind in rxn.reactantinds
            ropmat[i, ind] -= R
        end
    end
    return ropmat
end

function rops(ssys::SystemSimulation, t)
    domains = getfield.(ssys.sims, :domain)
    Nrxns = sum([length(sim.domain.phase.reactions) for sim in ssys.sims]) + sum([length(inter.reactions) for inter in ssys.interfaces if hasproperty(inter, :reactions)])
    Nspcs = sum([length(getphasespecies(sim.domain.phase)) for sim in ssys.sims])
    cstot = zeros(Nspcs)
    vns = Array{Any,1}(undef, length(domains))
    vcs = Array{Any,1}(undef, length(domains))
    vT = Array{Any,1}(undef, length(domains))
    vP = Array{Any,1}(undef, length(domains))
    vV = Array{Any,1}(undef, length(domains))
    vC = Array{Any,1}(undef, length(domains))
    vN = Array{Any,1}(undef, length(domains))
    vmu = Array{Any,1}(undef, length(domains))
    vkfs = Array{Any,1}(undef, length(domains))
    vkrevs = Array{Any,1}(undef, length(domains))
    vHs = Array{Any,1}(undef, length(domains))
    vUs = Array{Any,1}(undef, length(domains))
    vGs = Array{Any,1}(undef, length(domains))
    vdiffs = Array{Any,1}(undef, length(domains))
    vCvave = Array{Any,1}(undef, length(domains))
    vphi = Array{Any,1}(undef, length(domains))
    ropmat = spzeros(Nrxns, Nspcs)
    start = 1
    for (k, sim) in enumerate(ssys.sims)
        vns[k], vcs[k], vT[k], vP[k], vV[k], vC[k], vN[k], vmu[k], vkfs[k], vkrevs[k], vHs[k], vUs[k], vGs[k], vdiffs[k], vCvave[k], vphi[k] = calcthermo(sim.domain, ssys.sol(t), t)
        cstot[sim.domain.indexes[1]:sim.domain.indexes[2]] = vcs[k]
        rops!(ropmat, sim.domain.rxnarray, cstot, vkfs[k], vkrevs[k], vV[k], start)
        start += length(vkfs[k])
    end
    for inter in ssys.interfaces
        if inter isa FragmentBasedReactiveFilmGrowthInterfaceConstantT
            kfs, krevs = getkfskrevs(inter)
            rops!(ropmat, inter.rxnarray, inter.fragmentbasedrxnarray, cstot, kfs, krevs, vV[inter.domaininds[1]], start)
            start += length(kfs)
        elseif hasproperty(inter, :reactions)
            kfs, krevs = getkfskrevs(inter, vT[inter.domaininds[1]], vT[inter.domaininds[2]], vphi[inter.domaininds[1]], vphi[inter.domaininds[2]], vGs[inter.domaininds[1]], vGs[inter.domaininds[2]], cstot)
            rops!(ropmat, inter.rxnarray, cstot, kfs, krevs, inter.A, start)
            start += length(kfs)
        end
    end
    return ropmat
end

"""
calculates the rates of production/loss at a given time point for a single species
this outputs a sparse vector of length num reactions containing the production/loss
rate associated with that reaction for the given species
"""
function rops(bsol::Y, name::X, t::Z) where {Y<:Simulation,X<:AbstractString,Z<:Real}
    rop = spzeros(length(bsol.domain.phase.reactions))
    cs, kfs, krevs = calcthermo(bsol.domain, bsol.sol(t), t)[[2, 9, 10]]
    V = getdomainsize(bsol, t)
    ind = findfirst(isequal(name), getfield.(getphasespecies(bsol.domain.phase), :name))
    @assert !isa(ind, Nothing) "species $name not in species array"
    for (i, rxn) in enumerate(bsol.domain.phase.reactions)
        c = 0
        R = getrate(rxn, cs, kfs, krevs) * V
        c -= count(isequal(ind), rxn.reactantinds)
        c += count(isequal(ind), rxn.productinds)
        if c != 0
            rop[i] = c * R
        end
    end
    return rop
end

function rops(ssys::SystemSimulation, name, t)
    domains = getfield.(ssys.sims, :domain)
    ind = findfirst(isequal(name), ssys.names)
    Nrxns = sum([length(sim.domain.phase.reactions) for sim in ssys.sims]) + sum(Vector{Int}([length(inter.reactions) for inter in ssys.interfaces if hasproperty(inter, :reactions)]))
    Nspcs = sum([length(getphasespecies(sim.domain.phase)) for sim in ssys.sims])
    cstot = zeros(Nspcs)
    vns = Array{Any,1}(undef, length(domains))
    vcs = Array{Any,1}(undef, length(domains))
    vT = Array{Any,1}(undef, length(domains))
    vP = Array{Any,1}(undef, length(domains))
    vV = Array{Any,1}(undef, length(domains))
    vC = Array{Any,1}(undef, length(domains))
    vN = Array{Any,1}(undef, length(domains))
    vmu = Array{Any,1}(undef, length(domains))
    vkfs = Array{Any,1}(undef, length(domains))
    vkrevs = Array{Any,1}(undef, length(domains))
    vHs = Array{Any,1}(undef, length(domains))
    vUs = Array{Any,1}(undef, length(domains))
    vGs = Array{Any,1}(undef, length(domains))
    vdiffs = Array{Any,1}(undef, length(domains))
    vCvave = Array{Any,1}(undef, length(domains))
    vphi = Array{Any,1}(undef, length(domains))
    ropvec = spzeros(Nrxns)
    start = 0
    for (k, sim) in enumerate(ssys.sims)
        vns[k], vcs[k], vT[k], vP[k], vV[k], vC[k], vN[k], vmu[k], vkfs[k], vkrevs[k], vHs[k], vUs[k], vGs[k], vdiffs[k], vCvave[k], vphi[k] = calcthermo(sim.domain, ssys.sol(t), t)
        cstot[sim.domain.indexes[1]:sim.domain.indexes[2]] = vcs[k]
        rops!(ropvec, sim.domain.rxnarray, cstot, vkfs[k], vkrevs[k], vV[k], start, ind)
        start += length(vkfs[k])
    end
    for inter in ssys.interfaces
        if hasproperty(inter, :reactions)
            kfs, krevs = getkfskrevs(inter, vT[inter.domaininds[1]], vT[inter.domaininds[2]], vphi[inter.domaininds[1]], vphi[inter.domaininds[2]], vGs[inter.domaininds[1]], vGs[inter.domaininds[2]], cstot)
            rops!(ropvec, inter.rxnarray, cstot, kfs, krevs, inter.A, start, ind)
            start += length(kfs)
        end
    end
    return ropvec
end

export rops

function getrate(rarray, cs, kfs, krevs, V, i)
    if @inbounds rarray[2, i] == 0
        @inbounds @fastmath fR = kfs[i] * cs[rarray[1, i]]
    elseif @inbounds rarray[3, i] == 0
        @inbounds @fastmath fR = kfs[i] * cs[rarray[1, i]] * cs[rarray[2, i]]
    elseif @inbounds rarray[4, i] == 0
        @inbounds @fastmath fR = kfs[i] * cs[rarray[1, i]] * cs[rarray[2, i]] * cs[rarray[3, i]]
    else
        @inbounds @fastmath fR = kfs[i] * cs[rarray[1, i]] * cs[rarray[2, i]] * cs[rarray[3, i]] * cs[rarray[4, i]]
    end
    if @inbounds rarray[6, i] == 0
        @inbounds @fastmath rR = krevs[i] * cs[rarray[5, i]]
    elseif @inbounds rarray[7, i] == 0
        @inbounds @fastmath rR = krevs[i] * cs[rarray[5, i]] * cs[rarray[6, i]]
    elseif rarray[8, i] == 0
        @inbounds @fastmath rR = krevs[i] * cs[rarray[5, i]] * cs[rarray[6, i]] * cs[rarray[7, i]]
    else
        @inbounds @fastmath rR = krevs[i] * cs[rarray[5, i]] * cs[rarray[6, i]] * cs[rarray[7, i]] * cs[rarray[8, i]]
    end
    @fastmath R = (fR - rR) * V
    return R
end

function rops!(ropmat, rarray::Array{Int64,2}, cs, kfs, krevs, V, start)
    for i = 1:length(kfs)
        R = getrate(rarray, cs, kfs, krevs, V, i)

        @inbounds @fastmath ropmat[i+start, rarray[1, i]] -= R
        if @inbounds rarray[2, i] != 0
            @inbounds @fastmath ropmat[i+start, rarray[2, i]] -= R
            if @inbounds rarray[3, i] != 0
                @inbounds @fastmath ropmat[i+start, rarray[3, i]] -= R
                if @inbounds rarray[4, i] != 0
                    @inbounds @fastmath ropmat[i+start, rarray[4, i]] -= R
                end
            end
        end
        @inbounds @fastmath ropmat[i+start, rarray[5, i]] += R
        if @inbounds rarray[6, i] != 0
            @inbounds @fastmath ropmat[i+start, rarray[6, i]] += R
            if @inbounds rarray[7, i] != 0
                @inbounds @fastmath ropmat[i+start, rarray[7, i]] += R
                if @inbounds rarray[8, i] != 0
                    @inbounds @fastmath ropmat[i+start, rarray[8, i]] += R
                end
            end
        end
    end
end

function rops!(ropvec, rarray::Array{Int64,2}, cs, kfs, krevs, V, start, ind)
    for i = 1:length(kfs)
        c = count(isequal(ind), rarray[5:8, i]) - count(isequal(ind), rarray[1:4, i])
        if c != 0.0
            R = getrate(rarray, cs, kfs, krevs, V, i)
            @fastmath @inbounds ropvec[i+start] = c * R
        end
    end
end


function rops!(ropmat, rarray::Array{Int64,2}, fragmentbasedrxnarray::Array{Int64,2}, cs, kfs, krevs, V, start)
    numfragmentbasedreacprod, numrxns = size(fragmentbasedrxnarray)
    half = Int(numfragmentbasedreacprod / 2)
    for i = 1:length(kfs)
        R = getrate(rarray, cs, kfs, krevs, V, i)

        for j = 1:half
            if fragmentbasedrxnarray[j, i] != 0
                @fastmath ropmat[fragmentbasedrxnarray[j, i]] -= R
            end
        end

        for j = half+1:numfragmentbasedreacprod
            if fragmentbasedrxnarray[j, i] != 0
                @fastmath ropmat[fragmentbasedrxnarray[j, i]] += R
            end
        end
    end
end

function rops!(ropvec, rarray::Array{Int64,2}, fragmentbasedrxnarray::Array{Int64,2}, cs, kfs, krevs, V, start, ind)
    numfragmentbasedreacprod, numrxns = size(fragmentbasedrxnarray)
    half = Int(numfragmentbasedreacprod / 2)
    for i = 1:length(kfs)
        @views c = count(isequal(ind), fragmentbasedrxnarray[half+1:end, i]) - count(isequal(ind), fragmentbasedrxnarray[1:half, i])
        if c != 0.0
            R = getrate(rarray, cs, kfs, krevs, V, i)
            @fastmath @inbounds ropvec[i+start] = c * R
        end
    end
end

"""
Calculates sensitivities with respect to `target` at the time point at the end of the simulation
The returned sensitivities are the normalized values

By default uses the InterpolatingAdjoint algorithm with vector Jacobian products calculated with ReverseDiffVJP(true)
this assumes no changes in code branching during simulation, if that were to become no longer true, the Tracker 
based alternative algorithm is slower, but avoids this concern. 
"""
function getadjointsensitivities(bsol::Q, target::String, solver::W; sensalg::W2=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)),
    abstol::Float64=1e-6, reltol::Float64=1e-3, normalize=true, kwargs...) where {Q,W,W2}
    @assert target in bsol.names || target in ["T", "V", "P", "mass"]

    pethane = 160

    if target in ["T", "V", "P", "mass"]
        if haskey(bsol.domain.thermovariabledict, target)
            ind = bsol.domain.thermovariabledict[target]
        else
            throw(error("$(bsol.domain) doesn't have $target in its thermovariables"))
        end
    else
        ind = findfirst(isequal(target), bsol.names)
        if isempty(bsol.interfaces)
            sensdomain, sensspcnames, senstooriginspcind, senstooriginrxnind = getsensdomain(bsol.domain, ind)
            if :thermovariabledict in fieldnames(typeof(bsol.domain))
                yinds = vcat(senstooriginspcind, collect(values(bsol.domain.thermovariabledict)))
            else
                yinds = vcat(senstooriginspcind)
            end
            pinds = vcat(senstooriginspcind, length(bsol.domain.phase.species) .+ senstooriginrxnind)
            ind = findfirst(isequal(target), sensspcnames)
        end
    end

    function sensg(y::X, p::Array{Y,1}, t::Z) where {Q,V,X,Y<:Float64,Z}
        sensy = y[yinds]
        sensp = p[pinds]
        dy = similar(sensy, length(sensy))
        return dydtreactor!(dy, sensy, t, sensdomain, [], p=sensp)[ind]
    end
    function sensg(y::Array{X,1}, p::Y, t::Z) where {Q,V,X<:Float64,Y,Z}
        sensy = y[yinds]
        sensp = p[pinds]
        dy = similar(sensp, length(sensy))
        return dydtreactor!(dy, sensy, t, sensdomain, [], p=sensp)[ind]
    end
    function sensg(y::Array{X,1}, p::Array{Y,1}, t::Z) where {Q,V,X<:Float64,Y<:Float64,Z}
        sensy = y[yinds]
        sensp = p[pinds]
        dy = similar(sensy, length(sensy))
        return dydtreactor!(dy, sensy, t, sensdomain, [], p=sensp)[ind]
    end
    function sensg(y::Array{X,1}, p::Array{Y,1}, t::Z) where {Q,V,X<:ForwardDiff.Dual,Y<:ForwardDiff.Dual,Z}
        sensy = y[yinds]
        sensp = p[pinds]
        dy = similar(sensy, length(sensy))
        return dydtreactor!(dy, sensy, t, sensdomain, [], p=sensp)[ind]
    end

    function g(y::X, p::Array{Y,1}, t::Z) where {Q,V,X,Y<:Float64,Z}
        dy = similar(y, length(y))
        return dydtreactor!(dy, y, t, bsol.domain, bsol.interfaces, p=p)[ind]
    end
    function g(y::Array{X,1}, p::Y, t::Z) where {Q,V,X<:Float64,Y,Z}
        dy = similar(p, length(y))
        return dydtreactor!(dy, y, t, bsol.domain, bsol.interfaces, p=p)[ind]
    end
    function g(y::Array{X,1}, p::Array{Y,1}, t::Z) where {Q,V,X<:Float64,Y<:Float64,Z}
        dy = zeros(length(y))
        return dydtreactor!(dy, y, t, bsol.domain, bsol.interfaces, p=p)[ind]
    end
    function g(y::Array{X,1}, p::Array{Y,1}, t::Z) where {Q,V,X<:ForwardDiff.Dual,Y<:ForwardDiff.Dual,Z}
        dy = similar(y, length(y))
        return dydtreactor!(dy, y, t, bsol.domain, bsol.interfaces, p=p)[ind]
    end

    dsensgdu(out, y, p, t) = ForwardDiff.gradient!(out, y -> sensg(y, p, t), y)
    dsensgdp(out, y, p, t) = ForwardDiff.gradient!(out, p -> sensg(y, p, t), p)
    dgdu(out, y, p, t) = ForwardDiff.gradient!(out, y -> g(y, p, t), y)
    dgdp(out, y, p, t) = ForwardDiff.gradient!(out, p -> g(y, p, t), p)
    dsensgdurevdiff(out, y, p, t) = ReverseDiff.gradient!(out, y -> sensg(y, p, t), y)
    dsensgdprevdiff(out, y, p, t) = ReverseDiff.gradient!(out, p -> sensg(y, p, t), p)
    dgdurevdiff(out, y, p, t) = ReverseDiff.gradient!(out, y -> g(y, p, t), y)
    dgdprevdiff(out, y, p, t) = ReverseDiff.gradient!(out, p -> g(y, p, t), p)

    if length(bsol.domain.p) <= pethane
        if target in ["T", "V", "P", "mass"] || !isempty(bsol.interfaces)
            du0, dpadj = adjoint_sensitivities(bsol.sol, solver, g, nothing, (dgdu, dgdp); sensealg=sensalg, abstol=abstol, reltol=reltol, kwargs...)
        else
            du0, dpadj = adjoint_sensitivities(bsol.sol, solver, sensg, nothing, (dsensgdu, dsensgdp); sensealg=sensalg, abstol=abstol, reltol=reltol, kwargs...)
        end
    else
        if target in ["T", "V", "P", "mass"] || !isempty(bsol.interfaces)
            du0, dpadj = adjoint_sensitivities(bsol.sol, solver, g, nothing, (dgdurevdiff, dgdprevdiff); sensealg=sensalg, abstol=abstol, reltol=reltol, kwargs...)
        else
            du0, dpadj = adjoint_sensitivities(bsol.sol, solver, sensg, nothing, (dsensgdurevdiff, dsensgdprevdiff); sensealg=sensalg, abstol=abstol, reltol=reltol, kwargs...)
        end
    end
    if normalize
        dpadj[length(bsol.domain.phase.species)+1:end] .*= bsol.domain.p[length(bsol.domain.phase.species)+1:end]
        if !(target in ["T", "V", "P", "mass"]) && isempty(bsol.interfaces)
            dpadj ./= bsol.sol(bsol.sol.t[end])[senstooriginspcind[ind]]
        elseif !(target in ["T", "V", "P", "mass"]) && !isempty(bsol.interfaces)
            dpadj ./= bsol.sol(bsol.sol.t[end])[ind]
        end
    end
    return dpadj
end

function getadjointsensitivities(syssim::Q, bsol::W3, target::String, solver::W; sensalg::W2=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)),
    abstol::Float64=1e-6, reltol::Float64=1e-3, normalize=true, kwargs...) where {Q,W,W2,W3}
    @assert target in bsol.names || target in ["T", "V", "P", "mass"]
    if target in ["T", "V", "P", "mass"]
        if haskey(bsol.domain.thermovariabledict, target)
            ind = bsol.domain.thermovariabledict[target]
        else
            throw(error("$(bsol.domain) doesn't have $target in its thermovariables"))
        end
    else
        ind = findfirst(isequal(target), bsol.names) + bsol.domain.indexes[1] - 1
    end
    domains = Tuple([x.domain for x in syssim.sims])
    function g(y::X, p::Array{Y,1}, t::Z) where {Q,V,X,Y<:Float64,Z}
        dy = similar(y, length(y))
        return dydtreactor!(dy, y, t, domains, syssim.interfaces, p=p)[ind]
    end
    function g(y::Array{X,1}, p::Y, t::Z) where {Q,V,X<:Float64,Y,Z}
        dy = similar(p, length(y))
        return dydtreactor!(dy, y, t, domains, syssim.interfaces, p=p)[ind]
    end
    function g(y::Array{Float64,1}, p::Array{Float64,1}, t::Z) where {Q,V,Z}
        dy = similar(p, length(y))
        return dydtreactor!(dy, y, t, domains, syssim.interfaces, p=p)[ind]
    end
    function g(y::Array{X,1}, p::Array{Y,1}, t::Z) where {Q,V,X<:ForwardDiff.Dual,Y<:ForwardDiff.Dual,Z}
        dy = similar(y, length(y))
        return dydtreactor!(dy, y, t, domains, syssim.interfaces, p=p)[ind]
    end
    dgdu(out, y, p, t) = ForwardDiff.gradient!(out, y -> g(y, p, t), y)
    dgdp(out, y, p, t) = ForwardDiff.gradient!(out, p -> g(y, p, t), p)
    du0, dpadj = adjoint_sensitivities(syssim.sol, solver, g, nothing, (dgdu, dgdp); sensealg=sensalg, abstol=abstol, reltol=reltol, kwargs...)
    if normalize
        for domain in domains
            dpadj[domain.parameterindexes[1]+length(domain.phase.species):domain.parameterindexes[2]] .*= syssim.p[domain.parameterindexes[1]+length(domain.phase.species):domain.parameterindexes[2]]
        end
        if !(target in ["T", "V", "P", "mass"])
            dpadj ./= bsol.sol(bsol.sol.t[end])[ind]
        end
    end
    return dpadj
end
export getadjointsensitivities

function getconcentrationsensitivity(bsol::Simulation{Q,W,L,G}, numerator::String, denominator::String, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain,ParametrizedTConstantVDomain},K<:Real,Q,G,L}
    @assert numerator in bsol.names
    @assert denominator in bsol.names
    indnum = findfirst(isequal(numerator), bsol.names) + bsol.domain.indexes[1] - 1
    inddeno = findfirst(isequal(denominator), bsol.names) + bsol.domain.parameterindexes[1] - 1
    Nvars = length(getphasespecies(bsol.domain.phase)) + length(bsol.domain.indexes) - 2
    Nrxns = length(bsol.domain.phase.reactions)
    x, dp = extract_local_sensitivities(bsol.sol, t)
    s = dp[inddeno][indnum]
    val = s / bsol.sol(t)[indnum] #constant volume
    if t == 0
        return 0.0
    else
        return val
    end
end

function getconcentrationsensitivity(bsol::Simulation{Q,W,L,G}, numerator::String, denominator::String, t::K) where {W<:Union{ConstantTPDomain,ParametrizedTPDomain,ConstantPDomain,ParametrizedPDomain,ParametrizedVDomain},K<:Real,Q,G,L}
    @assert numerator in bsol.names
    @assert denominator in bsol.names
    indnum = findfirst(isequal(numerator), bsol.names) + bsol.domain.indexes[1] - 1
    inddeno = findfirst(isequal(denominator), bsol.names) + bsol.domain.parameterindexes[1] - 1
    Nvars = length(getphasespecies(bsol.domain.phase)) + length(bsol.domain.indexes) - 2
    Nrxns = length(bsol.domain.phase.reactions)
    x, dp = extract_local_sensitivities(bsol.sol, t)
    svals = dp[inddeno][bsol.domain.indexes[1]:bsol.domain.indexes[2]]
    s = svals[indnum]
    V = getV(bsol, t)
    c = bsol.sol(t)[indnum] / V
    val = (s - c * sum(svals) * R * getT(bsol, t) / getP(bsol, t)) / (c * V) #known T and P
    if t == 0
        return 0.0
    else
        return val
    end
end

function getconcentrationsensitivity(bsol::Simulation{Q,W,L,G}, numerator::String, denominator::Z, t::K) where {W<:Union{ConstantVDomain,ConstantTVDomain,ParametrizedTConstantVDomain},K<:Real,Z<:Integer,Q,G,L}
    @assert numerator in bsol.names
    indnum = findfirst(isequal(numerator), bsol.names) + bsol.domain.indexes[1] - 1
    inddeno = denominator + bsol.domain.parameterindexes[1] - 1
    Nvars = length(getphasespecies(bsol.domain.phase)) + length(bsol.domain.indexes) - 2
    Nrxns = length(bsol.domain.phase.reactions)
    x, dp = extract_local_sensitivities(bsol.sol, t)
    s = dp[inddeno+length(bsol.domain.phase.species)][indnum]
    T = getT(bsol, t)
    P = getP(bsol, t)
    C = getC(bsol, t)
    k = bsol.domain.p[inddeno+length(bsol.domain.phase.species)]
    val = s * k / bsol.sol(t)[indnum] #constant volume
    if t == 0
        return 0.0
    else
        return val
    end
end

function getconcentrationsensitivity(bsol::Simulation{Q,W,L,G}, numerator::String, denominator::Z, t::K) where {W<:Union{ConstantTPDomain,ParametrizedTPDomain,ConstantPDomain,ParametrizedPDomain,ParametrizedVDomain},K<:Real,Z<:Integer,Q,G,L}
    @assert numerator in bsol.names
    indnum = findfirst(isequal(numerator), bsol.names) + bsol.domain.indexes[1] - 1
    inddeno = denominator + bsol.domain.parameterindexes[1] - 1
    Nvars = length(getphasespecies(bsol.domain.phase)) + length(bsol.domain.indexes) - 2
    Nrxns = length(bsol.domain.phase.reactions)
    x, dp = extract_local_sensitivities(bsol.sol, t)
    svals = dp[inddeno+length(bsol.domain.phase.species)][bsol.domain.indexes[1]:bsol.domain.indexes[2]]
    s = svals[indnum]
    V = getV(bsol, t)
    T = getT(bsol, t)
    P = getP(bsol, t)
    C = getC(bsol, t)
    c = bsol.sol(t)[indnum] / V
    k = bsol.domain.p[inddeno+length(bsol.domain.phase.species)]
    val = k * (s - c * sum(svals) * R * T / P) / (c * V) #known T and P
    if t == 0
        return 0.0
    else
        return val
    end
end

export getconcentrationsensitivity

"""
calculate the rates of all reactions at time t
"""
function rates(bsol::Q, t::X) where {Q<:Simulation,X<:Real}
    cs, kfs, krevs = calcthermo(bsol.domain, bsol.sol(t), t)[[2, 9, 10]]
    V = getdomainsize(bsol, t)
    return [getrate(rxn, cs, kfs, krevs) * V for rxn in bsol.domain.phase.reactions]
end

"""
calculate the rates of all reactions at given times ts
defaults to using bsol.sol.t if ts is not supplied
"""
function rates(bsol::Q; ts::X=Array{Float64,1}()) where {Q<:Simulation,X<:AbstractArray}
    if length(ts) == 0
        ts = bsol.sol.t
    end
    return hcat([rates(bsol, t) for t in ts]...)
end

"""
calculate the rates of all reactions at time t
"""
function rates(ssys::Q, t::X) where {Q<:SystemSimulation,X<:Real}
    rts = zeros(length(ssys.reactions))
    domains = getfield.(ssys.sims, :domain)
    Nrxns = sum([length(sim.domain.phase.reactions) for sim in ssys.sims]) + sum([hasproperty(inter, :reactions) ? length(inter.reactions) : 0 for inter in ssys.interfaces])
    Nspcs = sum([length(getphasespecies(sim.domain.phase)) for sim in ssys.sims])
    cstot = zeros(Nspcs)
    vns = Array{Any,1}(undef, length(domains))
    vcs = Array{Any,1}(undef, length(domains))
    vT = Array{Any,1}(undef, length(domains))
    vP = Array{Any,1}(undef, length(domains))
    vV = Array{Any,1}(undef, length(domains))
    vC = Array{Any,1}(undef, length(domains))
    vN = Array{Any,1}(undef, length(domains))
    vmu = Array{Any,1}(undef, length(domains))
    vkfs = Array{Any,1}(undef, length(domains))
    vkrevs = Array{Any,1}(undef, length(domains))
    vHs = Array{Any,1}(undef, length(domains))
    vUs = Array{Any,1}(undef, length(domains))
    vGs = Array{Any,1}(undef, length(domains))
    vdiffs = Array{Any,1}(undef, length(domains))
    vCvave = Array{Any,1}(undef, length(domains))
    vphi = Array{Any,1}(undef, length(domains))
    index = 1
    for (k, sim) in enumerate(ssys.sims)
        vns[k], vcs[k], vT[k], vP[k], vV[k], vC[k], vN[k], vmu[k], vkfs[k], vkrevs[k], vHs[k], vUs[k], vGs[k], vdiffs[k], vCvave[k], vphi[k] = calcthermo(sim.domain, ssys.sol(t), t)
        cstot[sim.domain.indexes[1]:sim.domain.indexes[2]] = vcs[k]
        rts[index:index+length(vkfs[k])-1] .= getrates(sim.domain.rxnarray, cstot, vkfs[k], vkrevs[k]) .* getdomainsize(sim, t)
        index += length(vkfs[k])
    end
    for inter in ssys.interfaces
        if inter isa FragmentBasedReactiveFilmGrowthInterfaceConstantT
            kfs, krevs = getkfskrevs(inter)
            @views rts[index:index+length(kfs)-1] = getrates(inter.rxnarray, cstot, kfs, krevs) .* vV[inter.domaininds[1]]
            index += length(kfs)
        elseif hasproperty(inter, :reactions)
            kfs, krevs = getkfskrevs(inter, vT[inter.domaininds[1]], vT[inter.domaininds[2]], vphi[inter.domaininds[1]], vphi[inter.domaininds[2]], vGs[inter.domaininds[1]], vGs[inter.domaininds[2]], cstot)
            rts[index:index+length(kfs)-1] = getrates(inter.rxnarray, cstot, kfs, krevs) .* inter.A
            index += length(kfs)
        end
    end
    return rts
end

"""
calculate the rates of all reactions at given times ts
defaults to using bsol.sol.t if ts is not supplied
"""
function rates(ssys::Q; ts::X=Array{Float64,1}()) where {Q<:SystemSimulation,X<:AbstractArray}
    if length(ts) == 0
        ts = ssys.sol.t
    end
    return hcat([rates(ssys, t) for t in ts]...)
end

export rates

"""
Save the simulation profile to a .csv file
"""
function save(sim::T, save_name::String) where {T<:Simulation}
    df = DataFrame(sim.sol)
    rename!(df, names(df)[sim.domain.indexes[1]:sim.domain.indexes[2]] .=> sim.names)
    for (thermovariable, index) in sim.domain.thermovariabledict
        rename!(df, names(df)[index] => thermovariable)
    end
    CSV.write(save_name, df)
end

function save(syss::T, save_name::String) where {T<:SystemSimulation}
    df = DataFrame(syss.sol)
    for sim in syss.sims
        rename!(df, names(df)[sim.domain.indexes[1]:sim.domain.indexes[2]] .=> sim.names .* "($(sim.domain.phase.name))")
        for (thermovariable, index) in sim.domain.thermovariabledict
            rename!(df, names(df)[index] => thermovariable)
        end
    end
    CSV.write(save_name, df)
end
export save

