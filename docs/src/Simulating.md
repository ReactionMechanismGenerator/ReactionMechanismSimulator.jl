# Simulations

## Defining a Phase

A phase in RMS defines how thermodynamic and kinetic parameters for individual species and reactions
are calculated.  For the `IdealGas` case defining a phase is quite simple:  

```
ig = IdealGas(spcs,rxns,name="gas")
```

where `spcs` and `rxns` and lists of `AbstractSpecies` and `AbstractReaction` objects respectively.  
For an `IdealDiluteSolution` it is slightly more complicated because some of these properties (currently
 primarily diffusion limitations) can be dependent on the solvent:  

```
solv = Solvent("octane",RiedelViscosity(-98.805,3905.5,14.103,-2.5112e-5,2.0))
liq = IdealDiluteSolution(spcs,rxns,solv;name="phase",diffusionlimited=true)
```

## Defining Initial Conditions

Initial conditions for the simulation are defined within a dictionary in SI units.  The keys `"T"`,`"P"`
and `"V"` correspond to thermodynamic values (`"V"` being the extensive volume) while species names
correspond to numbers of moles (extensive).  

For example, since the initial `V` can be calculated implicitly a valid `IdealGas` initial condition might be:  

```
Dict(["T"=>1000.0,"P"=>1e5,"H2"=>0.67,"O2"=>0.33])
```

In the case of an `IdealDiluteSolution` where pressure is assumed to not affect the thermodynamic state
it might look more like:  
```
Dict(["T"=>450.0,"V"=>1.0e-6,"octane"=>6.154e-3,"oxygen"=>4.953e-6])
```

## Defining a Domain

A domain in RMS is a homogeneous volume that contains a single phase.  The `AbstractDomain` object defines how
the thermodynamics of the volume evolve with respect to time.  For example in a `ConstantTPDomain` the
temperature and pressure are defined in the domain object and held constant over the simulation and the volume
is integrated for while in a `ConstantVDomain` (a constant V adiabatic reactor) the volume is kept constant
and the temperature is integrated for.

When constructing a Domain object constant concentration species can be defined (list of species names).  
During integration the derivatives with respect to time of these species will be kept at zero.  

Sensitivity analysis can also be requested on the Domain object.  

`IdealDiluteSolution` example:  
```
domain,y0 = ConstantTVDomain(phase=liq,initialconds=initialconds,constantspecies=["oxygen"])
```

`IdealGas` example:  
```
domain,y0 = ConstantTPDomain(phase=ig,initialconds=initialconds;sensitivity=true)
```

## Defining a Reactor object

The Reactor object exists primarily to automatically sets up the ODEProblem object for you to solve.  
It takes the domain, the initial condition vector returned when constructing the domain and a time interval.  

Example:  
```
react = Reactor(domain,y0,(0.0,150.1))
```

## Solving a Reactor object

RMS purposefully exposes the solver interface provide users with all the options available from
Julia's DifferentialEquations package.  The ODEProblem object is a field of the Reactor
object `react.ode` and can be solved as the user desires.  

Example:  
```
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12)
```

In general CVODE_BDF tends to work well on these problems.  
