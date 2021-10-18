#Time scale analysis
using DifferentialEquations
using FastGaussQuadrature

function getfractionbelow(histo,n,N)
    s = 0
    for w in histo.weights
        if w < n
            s += w
        else
            s += n
        end
    end
    return s/N
end

function getfractionleft(histo,tol,N)
    s = 0
    for i in 2:length(histo.edges[1])
        if histo.edges[1][i] < tol
            s += histo.weights[i-1]
        end
    end
    return s/N
end

function splitnum(histo,n)
    c = 0
    boo = histo.weights[1] < n
    for (i,x) in enumerate(histo.weights)
        if boo && x > n
            c += 1
            boo = false
        elseif !boo && x < n
            c += 1
            boo = true
        end
    end
    return c
end

function splitnuminds(histo,n)
    c = 0
    inds = Array{Int64,1}()
    boo = histo.weights[1] < n
    for (i,x) in enumerate(histo.weights)
        if boo && x > n
            c += 1
            boo = false
            push!(inds,i)
        elseif !boo && x < n
            c += 1
            boo = true
            push!(inds,i)
        end
    end
    return c,inds
end

"""
Identify the timescale that best divides the fastest relaxing species (usually radicals)
from the slower relaxing species (usually stable species)
This is done by analysis of the timescale histogram
If there are more than two timescale peaks this algorithm will pick
the timescale that separates the two fastest peaks

timescales are grained multiplicatively with taumax,taumin, with taures being the
multiplicative factor
usediag=true tells the algorithm to determine timescales from the diagonal of the jacobian
otherwise it will compute and use the eigenvalues, this hasn't been observed to have
a significant impact on dividing time scale prediction

Because this relies on constructing the timescale distribution of the species
It doesn't work for very small mechanisms with ~12 or less species and
may not be able to find a multimodal slice in these cases
"""
function getdividingtimescale(Jy;taumax=1e4,taumin=1e-18,taures=10.0^0.5,usediag=true)
    if usediag
        taus = 1.0./abs.(diag(Jy))
    else
        taus = 1.0./abs.(eigvals(Jy))
    end
    N = length(taus)
    histo = fit(StatsBase.Histogram, taus, 10.0.^(log10(taumin):log10(taures):log10(taumax)), closed=:left)
    Nmax = maximum(histo.weights)
    splits = zeros(Nmax)
    fracbelow = zeros(Nmax)
    for i = 1:Nmax
        splits[i] = splitnum(histo,i)
        fracbelow[i] = getfractionbelow(histo,i,N)
    end
    fracleft = cumsum(histo.weights)./N
    fourinds = Array{Int64,1}()
    fracleftvals = Array{Float64,1}()
    fracbelowvals = Array{Float64,1}()
    for i = 1:Nmax
        if splits[i] >= 4
            c,inds = splitnuminds(histo,i)
            ind = inds[2]
            fl = fracleft[ind]
            fb = fracbelow[i]
            push!(fourinds,ind)
            push!(fracleftvals,fl)
            push!(fracbelowvals,fb)
        end
    end
    if length(fourinds) == 0
        error("Could not find a multimodal slice")
    end
    ind = argmin((fracleftvals.- 0.5).^2 .+ (fracbelowvals .- 0.5).^2)
    tau = histo.edges[1][fourinds[ind]+1]
    return tau
end

function getdividingtimescale(sim,t;taumax=1e6,taumin=1e-18,taures=10.0^0.5,usediag=true)
    Jy = jacobiany(sim.sol(t),sim.domain.p,t,sim.domain,sim.interfaces);
    getdividingtimescale(Jy,taumax=taumax,taumin=taumin,taures=taures,usediag=usediag)
end
export getdividingtimescale

function jacobiany(sol,t,p)
    Jy = zeros(length(sol.prob.u0),length(sol.prob.u0))
    sol.prob.f.jac(Jy,sol(t),p,t)
    return Jy
end

function jacobianp(sol,t,p)
    Jp = zeros(length(sol.prob.u0),length(p))
    sol.prob.f.paramjac(Jp,sol(t),p,t)
    return Jp
end

function jacobianysparse(sol,t,p)
    Jy = spzeros(length(sol.prob.u0),length(sol.prob.u0))
    sol.prob.f.jac(Jy,sol(t),p,t)
    return Jy
end

function jacobianpsparse(sol,t,p)
    Jp = spzeros(length(sol.prob.u0),length(p))
    sol.prob.f.paramjac(Jp,sol(t),p,t)
    return Jp
end

function jacobianp(sol,t,p,ind)
   return ForwardDiff.jacobian(pv->sol.prob.f(sol(t),hcat(p[1:ind-1],[pv],p[ind+1:end]),t),p[ind])
end

function normalizefulltransitorysensitivities!(dSdt,sim::Simulation,t)
    y = sim.sol(t)
    ns = y[sim.domain.indexes[1]:sim.domain.indexes[2]]
    dSdt .*= sim.domain.p'
    dSdt[1:sim.domain.indexes[2],:] ./= ns
    return dSdt
end

function normalizefulltransitorysensitivities!(dSdt,ssys::SystemSimulation,t)
    y = ssys.sol(t)
    dSdt .*= ssys.p'
    for sim in ssys.sims
        ns = y[sim.domain.indexes[1]:sim.domain.indexes[2]]
        dSdt[sim.domain.indexes[1]:sim.domain.indexes[2],:] ./= ns
    end
    return dSdt
end

function normalizeparamtransitorysensitivities!(dSdt,sim::Simulation,t,ind)
    y = sim.sol(t)
    ns = y[sim.domain.indexes[1]:sim.domain.indexes[2]]
    dSdt .*= sim.domain.p[ind]
    dSdt[1:length(sim.domain.phase.species)] ./= ns
    return dSdt
end

function normalizeparamtransitorysensitivities!(dSdt,ssys::SystemSimulation,t,ind)
    y = ssys.sol(t)
    for sim in ssys.sims
        ns = y[sim.domain.indexes[1]:sim.domain.indexes[2]]
        dSdt .*= ssys.domain.p[ind]
        dSdt[sim.domain.indexes[1]:sim.domain.indexes[2]] ./= ns
    end
    return dSdt
end

function normalizeadjointtransitorysensitivities!(dSdt,sim::Simulation,t,ind)
    y = sim.sol(t)
    dSdt .*= sim.domain.p'
    if ind <= sim.domain.indexes[2]
        dSdt ./= y[ind]
    end
    return dSdt
end

function normalizeadjointtransitorysensitivities!(dSdt,ssys::SystemSimulation,t,ind)
    y = ssys.sol(t)
    for sim in ssys.sims
        dSdt[sim.domain.parameterindexes[1]:sim.domain.parameterindexes[2]] .*= ssys.domain.p[sim.domain.parameterindexes[1]:sim.domain.parameterindexes[2]]
        if ind >= sim.domain.indexes[1] && ind <= sim.domain.indexes[2]
            dSdt ./= y[ind]
        end
    end
    return dSdt
end

"""
Compute exact transitory sensitivities using Forward Sensitivity Analysis
"""
function transitorysensitivitiesfullexact(sim::Simulation,t;tau=NaN,
        normalized=true,solver=DifferentialEquations.CVODE_BDF(linear_solver=:GMRES),
        abstol=1e-16,reltol=1e-6)

    if isnan(tau)
        Jy = jacobiany(sim.sol,t,sim.p);
        tau = getdividingtimescale(Jy);
    end

    if tau == 0.0
        dSdt = jacobianp(sim.sol,t,sim.p);
    else
        tspan = (0.0,tau)
        react = Reactor(sim.domain,sim.sol(t),tspan,sim.interfaces;p=sim.sol.prob.p,forwardsensitivities=true);
        sol = solve(react.ode,solver,abstol=abstol,reltol=reltol)
        dSdt = reduce(hcat,DiffEqSensitivity.extract_local_sensitivities(sol,tau)[2])./tau
    end

    if normalized
        return normalizefulltransitorysensitivities!(dSdt,sim,t)
    else
        return dSdt
    end
end

"""
Compute exact transitory sensitivities using Forward Sensitivity Analysis
"""
function transitorysensitivitiesfullexact(ssys::SystemSimulation,t;tau=NaN,
        normalized=true,solver=DifferentialEquations.CVODE_BDF(linear_solver=:GMRES),
        abstol=1e-16,reltol=1e-6)

    if isnan(tau)
        Jy = jacobiany(ssys.sol,t,ssys.p);
        tau = getdividingtimescale(Jy);
    end

    if tau == 0.0
        dSdt = jacobianp(ssys.sol,t,ssys.p);
    else
        tspan = (0.0,tau)
        react = Reactor(ssys.domains,sim.sol(t),tspan,ssys.interfaces;p=ssys.p,forwardsensitivities=true);
        sol = solve(react.ode,solver,abstol=abstol,reltol=reltol)
        dSdt = reduce(hcat,DiffEqSensitivity.extract_local_sensitivities(sol,tau)[2])./tau
    end

    if normalized
        return normalizefulltransitorysensitivities!(dSdt,sim,t)
    else
        return dSdt
    end
end
export transitorysensitivitiesfullexact

"""
Compute approximate transitory sensitivities
Assumes constant jacobian y and p and then uses trapezoidal rule to
approximate the solution to the resulting forward sensitivity equations
requires one evaluation of both jacobiany and jacobianp
"""
function transitorysensitivitiesfulltrapezoidal(sim,t;tau=NaN,normalized=true)
    if tau == 0.0
        dSdt = jacobianp(sim.sol,t,sim.p);
    elseif isnan(tau)
        Jy = jacobiany(sim.sol,t,sim.p);
        tau = getdividingtimescale(Jy);
        Jp = jacobianp(sim.sol,t,sim.p);
        dSdt = (Jp .+ exp(tau*Jy)*Jp)./2.0
    else
        Jy = jacobiany(sim.sol,t,sim.p);
        Jp = jacobianp(sim.sol,t,sim.p);
        dSdt = (Jp .+ exp(tau*Jy)*Jp)./2.0
    end
    if normalized
        return normalizefulltransitorysensitivities!(dSdt,sim,t)
    else
        return dSdt
    end
end
export transitorysensitivitiesfulltrapezoidal

"""
Compute approximate transitory sensitivities with respect to one parameter
Assumes constant jacobian y and p and then uses trapezoidal rule to
approximate the solution to the resulting forward sensitivity equations
requires one evaluation of both jacobiany and jacobianp
"""
function transitorysensitivitiesparamtrapezoidal(sim,t,ind;tau=NaN,normalized=true)
    if tau == 0.0
        dSdt = jacobianp(sol,t,sim.p,ind);
    elseif isnan(tau)
        Jy = jacobiany(sim.sol,t,sim.p);
        tau = getdividingtimescale(Jy);
        Jp = jacobianp(sim.sol,t,sim.p,ind);
        dSdt = (Jp .+ exp(tau*Jy)*Jp)./2.0
    else
        Jy = jacobiany(sim.sol,t,sim.p);
        Jp = jacobianp(sim.sol,t,sim.p,ind);
        dSdt = (Jp .+ exp(tau*Jy)*Jp)./2.0
    end
    if normalized
        return normalizeparamtransitorysensitivities!(dSdt,sim,t)
    else
        return dSdt
    end
end
export transitorysensitivitiesparamtrapezoidal

"""
Compute exact transitory sensitivities with respect to one parameter
using forward sensitivity analysis
"""
function transitorysensitivitiesparamexact(sim::Simulation,t,ind;tau=NaN,
        normalized=true,solver=DifferentialEquations.CVODE_BDF(linear_solver=:GMRES),
        abstol=1e-16,reltol=1e-6)

    if isnan(tau)
        Jy = jacobiany(sim.sol,t,sim.p);
        tau = getdividingtimescale(Jy);
    end

    if tau == 0.0
        dSdt = jacobianp(sim.sol,t,sim.p,ind);
    else
        tspan = (0.0,tau)
        react = Reactor(sim.domain,sim.sol(t),tspan,sim.interfaces;p=sim.sol.prob.p,forwardsensitivities=true);
        fparam(u,p,t) = react.ode.f(u,hcat(sim.p[1:ind-1],[sim.p[ind]],sim.p[ind+1:end]),t)
        jacparam(u,p,tn) = jacobianp(sol,t+tn,p,ind)
        odeparam = remake(react.ode,f=fparam,paramjac=jacparam)
        sol = solve(odeparam,solver,abstol=abstol,reltol=reltol)
        dSdt = reduce(hcat,DiffEqSensitivity.extract_local_sensitivities(sol,tau)[2])./tau
    end

    if normalized
        return normalizeparamtransitorysensitivities!(dSdt,sim,t)
    else
        return dSdt
    end
end

"""
Compute exact transitory sensitivities with respect to one parameter
using forward sensitivity analysis
"""
function transitorysensitivitiesparamexact(ssys::SystemSimulation,t,ind;tau=NaN,
        normalized=true,solver=DifferentialEquations.CVODE_BDF(linear_solver=:GMRES),
        abstol=1e-16,reltol=1e-6)

    if isnan(tau)
        Jy = jacobiany(ssys.sol,t,ssys.p);
        tau = getdividingtimescale(Jy);
    end

    if tau == 0.0
        dSdt = jacobianp(ssys.sol,ssys.p,t,ind);
    else
        tspan = (0.0,tau)
        react = Reactor(ssys.domains,sim.sol(t),tspan,ssys.interfaces;p=ssys.p,forwardsensitivities=true);
        fparam(u,p,t) = react.ode.f(u,hcat(sim.sol.prob.p[1:ind-1],[sim.sol.prob.p[ind]],sim.sol.prob.p[ind+1:end]),t)
        jacparam(u,p,tn) = jacobianp(sol,t+tn,p,ind)
        odeparam = remake(react.ode,f=fparam,paramjac=jacparam)
        sol = solve(odeparam,solver,abstol=abstol,reltol=reltol)
        dSdt = reduce(hcat,DiffEqSensitivity.extract_local_sensitivities(sol,tau)[2])./tau
    end

    if normalized
        return normalizeparamtransitorysensitivities!(dSdt,sim,t)
    else
        return dSdt
    end
end
export transitorysensitivitiesparamexact

"""
compute exact transitory sensitivities using adjoint sensitivity analysis
abstol and reltol are for the adjoint simulation
appreciable speed up with minimal accuracy reduction can be achieved by
loosening these tolerances
does not run any jacobian computations if tau != 0, tau != NaN  is supplied and
the solver is jacobian free
"""
function transitorysensitivitiesadjointexact(sim::Simulation,t,name;tau=NaN,
        normalized=true,solver=DifferentialEquations.CVODE_BDF(),sensalg=InterpolatingAdjoint(),
        abstol=1e-16,reltol=1e-6)

    @assert name in sim.names || name in ["T","V","P"]
    if name in ["T","V","P"]
        if haskey(sim.domain.thermovariabledict, name)
            ind = sim.domain.thermovariabledict[name]
        else
            throw(error("$(sim.domain) doesn't have $name in its thermovariables"))
        end
    else
        ind = findfirst(isequal(name),sim.names)
    end

    if isnan(tau)
        Jy = jacobiany(sim.sol,t,sim.p);
        tau = getdividingtimescale(Jy);
    end

    if tau == 0.0
        dSdt = jacobianp(sim.sol,t,sim.p)[ind,:];
    else
        prob = remake(sim.sol.prob;u0=sim.sol(t),tspan=(0.0,tau));
        soladj = solve(prob,solver,reltol=reltol,abstol=abstol);
        simadj = Simulation(soladj,sim.domain);
        dSdt = getadjointsensitivities(simadj,name,solver;sensalg=sensalg,abstol=abstol,
            reltol=reltol,normalize=false)./tau;
    end

    if normalized
        return normalizeadjointtransitorysensitivities!(dSdt,sim,t,ind);
    else
        return dSdt
    end
end

"""
compute exact transitory sensitivities using adjoint sensitivity analysis
abstol and reltol are for the adjoint simulation
appreciable speed up with minimal accuracy reduction can be achieved by
loosening these tolerances
does not run any jacobian computations if tau != 0, tau != NaN  is supplied and
the solver is jacobian free
"""
function transitorysensitivitiesadjointexact(ssys::SystemSimulation,t,name;tau=NaN,
        normalized=true,solver=DifferentialEquations.CVODE_BDF(),sensalg=InterpolatingAdjoint(),
        abstol=1e-6,reltol=1e-3)

    @assert name in sim.names || name in ["T","V","P"]
    if name in ["T","V","P"]
        if haskey(sim.domain.thermovariabledict, name)
            ind = sim.domain.thermovariabledict[name]
        else
            throw(error("$(sim.domain) doesn't have $name in its thermovariables"))
        end
    else
        ind = findfirst(isequal(name),sim.names)
    end

    if isnan(tau)
        Jy = jacobianysparse(ssys.sol,t,ssys.p);
        tau = getdividingtimescale(Jy);
    end

    if tau == 0.0
        dSdt = jacobianp(ssys.sol,t,ssys.p)[ind,:];
    else
        prob = remake(ssys.sol.prob;u0=sim.sol(t),tspan=(0.0,tau));
        soladj = solve(prob,solver,reltol=reltol,abstol=abstol);
        ssysadj = SystemSimulation(soladj,[sim.domain for sim in ssys.sims],ssys.interfaces,ssys.p);
        dSdt = getadjointsensitivities(ssysadj,name,solver;sensalg=sensalg,abstol=abstol,reltol=reltol)./tau;
    end

    if normalized
        return normalizeadjointtransitorysensitivities!(dSdt,ssys,t,ind)
    else
        return dSdt
    end
end
export transitorysensitivitiesadjointexact
