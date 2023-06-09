abstract type AbstractDiffusiveFlux end

struct MixtureAveragedDiffusiveFlux <: AbstractDiffusiveFlux end

function mixtureaverageddiffusivities!(mixdiffs, phase::IdealGas, C, cs, T, P)
    mixdiffs .= 0.0
    
    sum_cs = sum(cs)
    for spcind in 1:length(phase.species)
        if sum_cs == cs[spcind] # pure species
            mixdiffs[spcind] = getChapmanEnskogSelfDiffusivity(phase.species[spcind].transport; T=T, P=P)
            return mixdiffs
        end
    end

    for spcind in 1:length(phase.species)
        cdivD = 0.0
        for otherspcind in 1:length(phase.species)
            if spcind != otherspcind
                bindarydiff = getChapmanEnskogBinaryDiffusivity(phase.species[spcind].transport, phase.species[otherspcind].transport; T=T, P=P)
                cdivD += cs[otherspcind] / bindarydiff
            end
        end
        mixdiffs[spcind] = (C - cs[spcind]) / cdivD
    end
    return mixdiffs
end

function mixtureaveragedthermalconductivity(phase::IdealGas, cs::N1,  T::N2, P::N3) where {N1<:AbstractArray, N2<:Number, N3<:Number}
    mixtureaveragedlambdamult = 0.0
    mixtureaveragedlambdadiv = 0.0
    C = P / (R * T)
    for spcind in 1:length(phase.species)
        spc = phase.species[spcind]
        lambda = getChapmanEnskogThermalConductivity(spc.transport, spc.thermo; T=T)
        molfrac = cs[spcind] / C
        mixtureaveragedlambdamult += molfrac * lambda
        mixtureaveragedlambdadiv += molfrac / lambda
    end
    return 0.5 * (mixtureaveragedlambdamult + 1.0 / mixtureaveragedlambdadiv)
end