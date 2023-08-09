abstract type AbstractPotential end
export AbstractPotential

struct EmptyPotential <: AbstractPotential end
export EmptyPotential

"""
Lennard-Jones potential parameters
Used later in binary diffusion coefficient calculations
See Eq. 11.6 in Kee et al. 2018
"""
@with_kw struct LennardJonesPotential <: AbstractPotential
    epsilon::Float64 #L-J depth
    sigma::Float64 #L-J diameter
end
export LennardJonesPotential

function getcollisionintegral11(potentiali::LennardJonesPotential, potentialj::LennardJonesPotential, epsilonij::Float64; T::Number,)

    a1 = 1.0548
    a2 = 0.15504
    a3 = 0.55909
    a4 = 2.1705

    reduced_T = T * kB / epsilonij
    omegaij = a1 * reduced_T^(-a2) + (reduced_T + a3)^(-a4) # Eq. 11.6 in Kee et al. 2018

    return omegaij
end

function getcollisionintegral22(potentiali::LennardJonesPotential, potentialj::LennardJonesPotential, epsilonij::Float64; T::Number, )

    b1 = 1.0413
    b2 = 0.11930
    b3 = 0.43628
    b4 = 1.6041

    reduced_T = T * kB / epsilonij
    omegaij = b1 * reduced_T^(-b2) + (reduced_T + b3)^(-b4) # Eq. 11.6 in Kee et al. 2018

    return omegaij
end