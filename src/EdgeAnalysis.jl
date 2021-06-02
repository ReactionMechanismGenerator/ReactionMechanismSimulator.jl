"""
Tools for model edge analysis for automatic mechanism generation
"""

using Logging
using Sundials
using SparseArrays
using DiffEqBase: build_solution

abstract type AbstractTerminationCriterion end

struct TerminationTime <: AbstractTerminationCriterion
    time::Float64
end
struct TerminationConversion <: AbstractTerminationCriterion
    species::Species
    conversion::Float64
end
struct TerminationRateRatio <: AbstractTerminationCriterion
    ratio::Float64
end

export TerminationTime
export TerminationConversion
export TerminationRateRatio

"""
Generate appropriate core Simulation/SystemSimulation object
"""
function tosim(sol,domains,inters,p)
    return Simulation(sol,domains,inters,p)
end

"""
Generate appropriate core Simulation/SystemSimulation object
"""
function tosim(sol,domains::Tuple,inters,p)
   return SystemSimulation(sol,domains,inters,p) 
end

"""
Generate appropriate edge Simulation/SystemSimulation object
"""
function getsim(inte,react,coreedgedomain,inters,p,coretoedgespcmap)
    ycoreedge = getycoreedge(inte.u,coretoedgespcmap,coreedgedomain.indexes[end])
    sol = build_solution(react.ode,inte.alg,[0.0,inte.t],[ycoreedge,ycoreedge])
    sim = Simulation(sol,coreedgedomain,inters,p)
    return sim
end

"""
Generate appropriate edge Simulation/SystemSimulation object
"""
function getsim(inte,react,coreedgedomains::Tuple,inters,p,coretoedgespcmap)
    ycoreedge = getycoreedge(inte.u,coretoedgespcmap,coreedgedomains[end].indexes[end])
    sol = build_solution(react.ode,inte.alg,[0.0,inte.t],[ycoreedge,ycoreedge])
    ssys = SystemSimulation(sol,coreedgedomains,inters,p)
    return ssys
end

"""
Generate a state vector appropriate for both the edge and core from the core state vector
"""
function getycoreedge(y,coretoedgespcmap,numedgespc)
    ycoreedge = zeros(numedgespc)
    for (coreind,edgeind) in coretoedgespcmap
        ycoreedge[edgeind] = y[coreind]
    end
    return ycoreedge
end