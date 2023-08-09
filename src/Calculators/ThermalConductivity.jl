
function getChapmanEnskogThermalConductivity(transport::TransportModel, thermo::NASA; T::N1) where {N1<:Number}
    m = transport.m
    potential = transport.potential
    epsilon = potential.epsilon
    sigma = potential.sigma
    cp = getHeatCapacity(thermo, T)
    cv = cp - R

    omega = getcollisionintegral22(potential, potential, epsilon; T=T)
    return 25 / (32 * pi^0.5) * (kB * T / m)^0.5 * cv / (sigma^2 * Na * omega)
end
export ChapmanEnskogThermalConductivity
