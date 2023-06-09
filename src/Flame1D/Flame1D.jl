abstract type AbstractFlame end
abstract type AxisymmetricFlame <: AbstractFlame end
abstract type FlatFlame <: AxisymmetricFlame end

function getrho(cs, molecularweights::Array{Float64,1})
    return dot(cs, molecularweights)
end

"""
`Burner Flame`: struct for 1D pre-mixed burner stabilized flame
- `phase`: IdealGas phase
- `inlet`: Inlet1D struct contains the burner condition
    - `T`: burner temperature (K)
    - `cs`: Array{Float64,1} contains the inlet molar concentration of species (mol/m^3)
- `indexes`: Array{Int64,1} contains [start index of species, end index of species, index of temperature, index of axial velocity]
- `P`: pressure (Pa)
- `mdot`: net mass flux (kg/m^2/s)
- `u0`: initial axial velocity (m/s)
- `rxnarray`: Array{Int64, 2} contains the reaction array of shape (number of reactions, 8) with indexes of participating species
- `p`: Array{Float64,1} contains the parameters
- `variable_index_dict`: Dict{String,Int64} contains the index of the variable in the solution vector
"""

struct BurnerFlame{N1<:AbstractDiffusiveFlux} <: FlatFlame
    phase::IdealGas
    inlet::Inlet1D
    molecularweights::Array{Float64,1}
    indexes::Array{Int64,1}
    P::Float64
    mdot::Float64
    u0::Float64
    rxnarray::Array{Int64, 2}
    p::Array{Float64,1}
    variable_index_dict::Dict{String,Int64}
    diffusive_flux_model::N1
end

"""
`Burner Flame`: constructor for 1D pre-mixed burner stabilized flame
- `phase`: IdealGas phase
- `initialcond`: Dict{String,Float64} contains the initial conditions for the flame. 
    It takes in the following keys and values:
        - species name: premixed gas molar concentration (mol/m^3)
        - "P": pressure (Pa)
        - "T": the burner temperature (K)
        - "mdot": the net mass flux (kg/m^2/s)
"""

function BurnerFlame(;phase::IdealGas, initialcond::Dict, diffusive_flux_model::String="mix")
    @assert isa(phase, IdealGas) "phase must be an IdealGas phase"
    @assert "T" in keys(initialcond) "Initial condition must contain the burner temperature, T (K)"
    @assert "P" in keys(initialcond) "Initial condition must contain the pressure, P (Pa)"
    @assert "mdot" in keys(initialcond) "Initial condition must contain the net mass flux, mdot (kg/m^2/s)"
    @assert diffusive_flux_model in ["mix"] "diffusive flux model must be one of the following: mix" # currently only supporting mixture averaged diffusive flux
    
    # set inlet conditions
    T = 0.0
    P = 0.0
    mdot = 0.0

    cs = zeros(length(phase.species)) # inlet species molar concentration (mol/m^3)
    spcnames = getfield.(phase.species, :name)
    molecularweights = getfield.(phase.species, :molecularweight)

    for (key, value) in initialcond
        if key == "T"
            T = value
        elseif key == "P"
            P = value
        elseif key == "mdot"
            mdot = value
        else
            ind = findfirst(isequal(key), spcnames)
            @assert !isnothing(ind) "Species $key not found in the mechanism"
            cs[ind] = value
        end
    end

    rho = getrho(cs, molecularweights) # density (kg/m^3)
    u = mdot / rho # axial velocity (m/s)

    inlet = Inlet1D(T, cs)

    variable_index_dict = phase.spcdict
    variable_index_dict["T"] = length(spcnames) + 1
    variable_index_dict["u"] = length(spcnames) + 2

    p = vcat(zeros(length(phase.species)), ones(length(phase.reactions)))

    if diffusive_flux_model == "mix"
        diffusive_flux_model = MixtureAveragedDiffusiveFlux()
    end
    return BurnerFlame(phase, inlet, molecularweights, [1, length(phase.species), length(phase.species) + 1, length(phase.species) + 2], P, mdot, u, phase.rxnarray, p, variable_index_dict, diffusive_flux_model), p
end

export BurnerFlame

"""
`FreeFlame`: struct for 1D pre-mixed freely propagating adiabatic flame
- `phase`: IdealGas phase
- `inlet`: Inlet1D struct contains the inlet condition
    - `T`: inlet temperature (K)
    - `cs`: Array{Float64,1} contains the inlet molar concentration of species (mol/m^3)
- `indexes`: Array{Int64,1} contains [start index of species, end index of species, index of temperature, index of axial velocity]
- `P`: pressure (Pa)
- `rxnarray`: Array{Int64, 2} contains the reaction array of shape (number of reactions, 8) with indexes of participating species
- `p`: Array{Float64,1} contains the parameters
- `variable_index_dict`: Dict{String,Int64} contains the index of the variable in the solution vector
"""

struct FreeFlame{N1<:AbstractDiffusiveFlux} <: FlatFlame
    phase::IdealGas
    inlet::Inlet1D
    molecularweights::Array{Float64,1}
    indexes::Array{Int64,1}
    P::Float64
    rxnarray::Array{Int64, 2}
    p::Array{Float64,1}
    variable_index_dict::Dict{String,Int64}
    diffusive_flux_model::N1
end

"""
`FreeFlame`: constructor for 1D pre-mixed freely propagating adiabatic flame
- `phase`: IdealGas phase
- `initialcond`: Dict{String,Float64} contains the initial conditions for the flame. 
    It takes in the following keys and values:
        - species name: premixed gas molar concentration (mol/m^3)
        - "P": pressure (Pa)
        - "T": the inlet temperature (K)
"""

function FreeFlame(;phase::IdealGas, initialcond::Dict, diffusive_flux_model::String="mix")
    @assert isa(phase, IdealGas) "phase must be an IdealGas phase"
    @assert "T" in keys(initialcond) "Initial condition must contain the burner temperature, T (K)"
    @assert "P" in keys(initialcond) "Initial condition must contain the pressure, P (Pa)"
    @assert diffusive_flux_model in ["mix"] "diffusive flux model must be one of the following: mix" # currently only supporting mixture averaged diffusive flux
    
    # set inlet conditions
    T = 0.0
    P = 0.0

    cs = zeros(length(phase.species)) # inlet species molar concentration (mol/m^3)
    spcnames = getfield.(phase.species, :name)
    molecularweights = getfield.(phase.species, :molecularweight)

    for (key, value) in initialcond
        if key == "T"
            T = value
        elseif key == "P"
            P = value
        else
            ind = findfirst(isequal(key), spcnames)
            @assert !isnothing(ind) "Species $key not found in the mechanism"
            cs[ind] = value
        end
    end

    variable_index_dict = phase.spcdict
    variable_index_dict["T"] = length(spcnames) + 1
    variable_index_dict["u"] = length(spcnames) + 2

    p = vcat(zeros(length(phase.species)), ones(length(phase.reactions)))

    inlet = Inlet1D(T, cs)

    if diffusive_flux_model == "mix"
        diffusive_flux_model = MixtureAveragedDiffusiveFlux()
    end

    return FreeFlame(phase, inlet, molecularweights, [1, length(phase.species), length(phase.species) + 1, length(phase.species) + 2], P, phase.rxnarray, p, variable_index_dict, diffusive_flux_model), p
end

export FreeFlame

function calcthermo(flame::N1, y, t, p, cell_size) where {N1<:AxisymmetricFlame}
    @views cs = y[flame.indexes[1]:flame.indexes[2]]
    T = y[flame.variable_index_dict["T"]]
    P = flame.P
    V = cell_size
    C = P / (R * T)
    N = P * V / (R * T)
    ns = cs .* V

    cpdivR, hdivRT, sdivR = calcHSCpdless(flame.phase.vecthermo, T)
    Gs = (hdivRT .- sdivR) .* (R * T)
    Hs = hdivRT .* (R * T)
    Cvave = dot(cpdivR, ns) * R / N - R

    kfs, krevs = getkfkrevs(flame.phase, T, P, C, N, ns, Gs, Array{Float64,1}(), V, 0.0)
    return ns, cs, T, P, V, C, N, 0.0, kfs, krevs, Hs, Array{Float64,1}(), Gs, Array{Float64,1}(), Cvave, cpdivR, 0.0
end

function calcdomainderivatives!(flame::N1, rhs; t, T, P, Us, Hs, V, C, ns, N, Cvave) where {N1<:AxisymmetricFlame}
    Cpave = Cvave + R
    @views rhs[flame.variable_index_dict["T"]] = - dot(Hs, rhs[flame.indexes[1]:flame.indexes[2]])/(N * Cpave)
end