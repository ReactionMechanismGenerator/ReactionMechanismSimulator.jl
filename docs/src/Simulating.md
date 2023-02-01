# Simulations

## Defining a Phase

A phase in RMS defines how thermodynamic and kinetic parameters for individual species and reactions
are calculated.  For the `IdealGas` case defining a phase is quite simple:  

```
ig = IdealGas(spcs,rxns,name="gas")
```

where `spcs` and `rxns` are lists of `AbstractSpecies` and `AbstractReaction` objects respectively.  
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

## Defining an Interface object

In many cases you will have connections between multiple domains in your system or between one domain and something outside the system. In these cases you need to define interfaces.

Example Inlet:
```
Flow(t::Float64) = 10.0
inletdict= Dict(["T"=>800.0,"P"=>10.0e5,"O2"=>0.21, "N2"=>0.79])
inlet = Inlet(domain,inletdict,Flow
```

Example Reactive Interface:
```
inter,pinter = ReactiveInternalInterfaceConstantTPhi(domainliq,domaincat,interfacerxns,Tinter,areainter)
```

## Defining a Domain

A domain in RMS is a homogeneous volume that contains a single phase.  The `AbstractDomain` object defines how
the thermodynamics of the volume evolve with respect to time.  For example in a `ConstantTPDomain` the
temperature and pressure are defined in the domain object and held constant over the simulation and the volume
is integrated for while in a `ConstantVDomain` (a constant V adiabatic reactor) the volume is kept constant
and the temperature is integrated for.

When constructing a Domain object constant concentration species can be defined (list of species names).  
During integration the derivatives with respect to time of these species will be kept at zero.  

`IdealDiluteSolution` example:  
```
domain,y0,p = ConstantTVDomain(phase=liq,initialconds=initialconds,constantspecies=["oxygen"])
```

`IdealGas` example:  
```
domain,y0,p = ConstantTPDomain(phase=ig,initialconds=initialconds)
```


## Defining a Reactor object

The Reactor object exists primarily to automatically sets up the ODEProblem object for you to solve.  
For single domain reactors it takes the domain, the initial condition vector returned when constructing the domain, a time interval and an array of any interface objects and returns a Reactor object. If your system has multiple domains you need to pass the domains, initialconditions time interval, array of interfaces and parameters with the domains, associated initial conditions and associated parameters as tuples in consistent order from domains to interfaces.

Example single domain:  
```
react = Reactor(domain,y0,(0.0,150.1),interfaces,p=p)
```

Example multiple domains:
```
react,y0,p = Reactor((domainliq,domaincat), (y0liq,y0cat), (0.0, 1.0e5), [inter], (pliq,pcat,pinter))
```

## Solving a Reactor object

RMS purposefully exposes the solver interface provide users with all the options available from
Julia's DifferentialEquations package.  The ODEProblem object is a field of the Reactor
object `react.ode` and can be solved as the user desires. A recommended solver choice is stored in `react.recommendedsolver`. User can also specify their own choice of solver.

Forward sensitivity analysis can also be requested on the Reactor object by setting `forwardsensitivities=true`. Note that adjoint and threaded sensitivity analyses are usually much faster. Adjoint sensitivity analysis can be done as postprocessing analysis after the simulation is complete without a need to set `forwardsensitivities=true` during the simulation (this are discussed in the Analysis section). Threaded sensitivity analysis will be discussed in the next section.

Example:

```
sol = solve(react.ode,react.recommendedsolver,abstol=1e-20,reltol=1e-12)
```

```
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12;forwardsensitivities=true)
```

In general CVODE_BDF tends to work well on these problems.  

## Threaded Sensitivity Analysis
Instead of solving all of the sensitivity equations at once as is done in raw Forward sensitivity analysis we can
first solve the equations without sensitivity analysis alone. With an interpolatable solution to the original equations the sensitivities associated with each parameter are decoupled from the sensitivities of every other parameter and can be solved independently. Solving these groups of equations sequentially is often significantly faster than solving the equations together. However, by parallelizing the solution of these equations using multithreading it is possible to achieve dramatic speed ups. This approach is not always competitive with adjoint sensitivities as the adjoint approach requires solution of a much smaller system of equations, however, in practice this approach is often more robust especially for large systems.

In order to take advantage of multithreading you must set the number of threads before launching Julia see the
documentation <a href="https://docs.julialang.org/en/v1/manual/multi-threading/">here</a>. We provide two methods for this algorithm. The first method below calculates the sensitivities of all variables to all parameters and provides a solution object that can be used in the same way as the output of a solve on a Reactor with forwardsensitivities=true:

```
using Base.Threads
nthreads = Threads.nthreads()
println("Using $nthreads threads")

sol = threadedsensitivities(react; odesolver=react.recommendedsolver,senssolver=react.recommendedsolver,
        odekwargs=Dict([:abstol=>1e-20,:reltol=>1e-6]),senskwargs=Dict([:abstol=>1e-6,:reltol=>1e-3]))
```

We also provide a second method shown below. This variant allows the user to specify the indices of the parameters to calculate sensitivities for. However, because not all of the sensitivies are calculated the results cannot be formatted into a solution object that can be fed into a `Simulation` or `SystemSimulation` object. So instead a dictionary mapping the parameter index to an ODE solution of the associated sensitivity equations is provided. Note  that parameters are structured starting in species order with parameters corresponding to the Gibbs free energy of each species and continuing in reaction order with a parameter corresponding to each rate coefficient. The output solution objects are in species + thermodynamic variable order and are the raw sensitivities with no normalization.
```
using Base.Threads
nthreads = Threads.nthreads()
println("Using $nthreads threads")

soldict = threadedsensitivities(react,indices; odesolver=react.recommendedsolver,senssolver=react.recommendedsolver,
        odekwargs=Dict([:abstol=>1e-20,:reltol=>1e-6]),senskwargs=Dict([:abstol=>1e-6,:reltol=>1e-3]))
```
