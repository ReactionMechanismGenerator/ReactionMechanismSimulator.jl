using Parameters

abstract type AbstractDiffusivity end
export AbstractDiffusivity

struct EmptyDiffusivity <: AbstractDiffusivity end
(ed::EmptyDiffusivity)(;T::N,mu::Q,P::R=0.0) where {N,R,Q<:Number} = Inf
export EmptyDiffusivity

@with_kw struct ConstantDiffusivity{N<:Number} <: AbstractDiffusivity
    D::N = nothing
end
(c::ConstantDiffusivity)(;T::AbstractFloat=0.0,P::AbstractFloat=0.0,mu::AbstractFloat=0.0) = c.D
export ConstantDiffusivity

@with_kw struct StokesDiffusivity{N<:Number} <: AbstractDiffusivity
    r::N
end
(sd::StokesDiffusivity)(;T::N,mu::Q,P::R=0.0) where {N,R,Q<:Number} = @fastmath kB*T/(6*Base.pi*mu*sd.r)
export StokesDiffusivity

"""
ChapmanEnskogBinaryDiffusivity: Calculate binary diffusion coefficient for species i and species j using Chapman-Enskog theory
and Lorentz-Berthelot combining rules
See Eq. 11.4 in Kee et al. 2018
"""

function getChapmanEnskogBinaryDiffusivity(transporti::TransportModel, transportj::TransportModel; T::N1, P::N2) where {N1<:Number,N2<:Number}
    mi = transporti.m
    mj = transportj.m

    mij = (mi * mj) / (mi + mj)

    epsiloni = transporti.potential.epsilon
    epsilonj = transportj.potential.epsilon
    sigmai = transporti.potential.sigma
    sigmaj = transportj.potential.sigma

    reduced_polari = transporti.polarizability / sigmai^3 # Eq. 11.41 in Kee et al. 2018
    reduced_dipolej = transportj.dipolemoment / sqrt(epsilonj * sigmaj^3) # Eq. 11.42 in Kee et al. 2018
    xi = 1 + 0.25 * reduced_polari * reduced_dipolej * sqrt(epsiloni / epsilonj) # induction energy term Eq. 11.40 in Kee et al. 2018

    epsilonij = xi^2 * sqrt(epsiloni * epsilonj) # Eq. 11.38 in Kee et al. 2018
    sigmaij = 0.5 * (sigmai + sigmaj) * xi^(-1 / 6) # Eq. 11.39 in Kee et al. 2018

    omegaij = getcollisionintegral11(transporti.potential, transportj.potential, epsilonij; T=T)

    return 3 / 16 * sqrt(2 * pi * kB^3 * T^3 / mij) / (P * pi * sigmaij^2 * omegaij)
end

function getChapmanEnskogSelfDiffusivity(transport::TransportModel; T::N1, P::N2) where {N1<:Number,N2<:Number}
    m = transport.m
    epsilon = transport.potential.epsilon
    sigma = transport.potential.sigma

    omegaii = getcollisionintegral11(transport.potential, transport.potential, epsilon; T=T)

    return 3 / 8 * sqrt(pi * kB^3 * T^3 / m) / (P * pi * sigma^2 * omegaii)
end
