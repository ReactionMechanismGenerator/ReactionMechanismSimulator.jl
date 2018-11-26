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

function molefractions(bsol::Q; name::W,t::E) where {Q<:AbstractSolution, W<:String, E<:Real}
    @assert name in bsol.names
    ind = findfirst(isequal(name),bsol.names)
    return bsol(t)[ind]/bsol.N(t)
end

function molefractions(bsol::Q; t::E) where {Q<:AbstractSolution,E<:Real}
    return bsol.sol(t)[bsol.domain.indexes[1]:bsol.domain.indexes[2]]./bsol.N(t)
end

function molefractions(bsol::Q) where {Q<:AbstractSolution}
    return bsol.sol.u./bsol.N.u
end

