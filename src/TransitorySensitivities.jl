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
