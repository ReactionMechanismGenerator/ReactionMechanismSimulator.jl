using DifferentialEquations
import DifferentialEquations.DiffEqBase: AbstractODESolution, HermiteInterpolation,AbstractDiffEqInterpolation

abstract type AbstractSolution end

struct BatchSolution{Q<:AbstractODESolution,W<:AbstractDomain,L<:AbstractArray,G<:Function} <: AbstractSolution
    sol::Q
    domain::W
    names::L
    N::G
end

function BatchSolution(sol::Q,domain::W) where {Q<:AbstractODESolution,W<:AbstractDomain}
    names = getfield.(domain.phase.species,:name)
    Ns = sum(hcat(sol.interp.u...)[domain.indexes[1]:domain.indexes[2],:],dims=1)
    Nderivs = sum(hcat(sol.interp.du...)[domain.indexes[1]:domain.indexes[2],:],dims=1)
    N = HermiteInterpolation(sol.interp.t,Ns,Nderivs)
    F(t::T) where {T<:Real} = N(t,nothing,Val{0},sol.prob.p,:left)
    return BatchSolution(sol,domain,names,F)
end
