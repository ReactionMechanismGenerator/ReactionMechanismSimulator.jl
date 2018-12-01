using Parameters
using MathProgBase
using Clp

"""
determine the moles of each species in phase at equilibrium
given an array of initial moles y and a temperature T
"""
function equilibrium(phase::Q,y::V,T::B) where {Q<:Union{IdealGas,IdealDiluteSolution},V<:AbstractArray,B<:Real}
    Gs = calcgibbs(phase,T)
    N = sum(y)
    y ./ N

    atoms = Array{String,1}()
    for spc in phase.species
        for x in keys(spc.atomnums)
            if !(x in atoms)
                push!(atoms,x)
            end
        end
    end

    constraintmat = zeros(length(atoms),length(phase.species)) #number of atoms in each species
    constrainteq = zeros(length(atoms)) #number of each atom type present

    for (i,spc) in enumerate(phase.species)
        for (j,atm) in enumerate(atoms)
            if atm in keys(spc.atomnums)
                constraintmat[j,i] = spc.atomnums[atm]
                constrainteq[j] += spc.atomnums[atm]*y[i]
            end
        end
    end

    out = linprog(Gs,constraintmat,'=',constrainteq,ClpSolver())
    return out.sol
end

"""
determine the moles of each species in phase at equilibrium
given an array of initial moles y and a temperature T
"""
function equilibrium(phase::Q,spcdict::Dict{String,V},T::B) where {Q<:Union{IdealGas,IdealDiluteSolution},V<:Real,B<:Real}
    Gs = calcgibbs(phase,T)
    y = makespcsvector(phase,spcdict)
    N = sum(y)
    y ./ N

    atoms = Array{String,1}()
    for spc in phase.species
        for x in keys(spc.atomnums)
            if !(x in atoms)
                push!(atoms,x)
            end
        end
    end

    constraintmat = zeros(length(atoms),length(phase.species)) #number of atoms in each species
    constrainteq = zeros(length(atoms)) #number of each atom type present

    for (i,spc) in enumerate(phase.species)
        for (j,atm) in enumerate(atoms)
            if atm in keys(spc.atomnums)
                constraintmat[j,i] = spc.atomnums[atm]
                constrainteq[j] += spc.atomnums[atm]*y[i]
            end
        end
    end

    out = linprog(Gs,constraintmat,'=',constrainteq,ClpSolver())
    return out.sol
end

export equilibrium
