# Analysis

## The Simulation Object

Since the solver interface is exposed the solution object generated only gives the raw moles
and raw sensitivity values, which typically are not what a user wants.  The Simulation object
combines the solution information and the domain information to calculate much more useful
properties of the solution.  

The Simulation object can be defined:  

```
bsol = Simulation(sol,domain,interfaces,p)
```

where `sol` is the ODESolution object output by the DifferentialEquations package, `domain` is the
domain `sol` corresponds to, interfaces is the array of interface objects and p is the parameter vector.

## The SystemSimulation Object

When a system involves multiple domains a single Simulation object is insufficient. For these systems you need to construct a SystemSimulation object. This works in much the same way except that the domains should be listed as a tuple and the domain and interface ordering should be the same as that used when constructing the Reactor object. In theory the SystemSimulation object should be able to be used in place of a Simulation object in most places (If this is not the case and should be the case please make an issue!).

```
ssys = SystemSimulation(sol,(domainliq,domaincat,),interfaces,p);
```

## Useful Properties

### Thermodynamic Properties

Any thermodynamic property (T,P,V,C) at any given time `t` can be calculated in the format
```
getT(bsol,t)
```
where T can be replaced by any of the other properties.  

### Mole Fractions

Mole fraction information can be retrieved using the `molefraction` function.  
`molefraction(bsol)` will give the mole fractions of all species at the times `bsol.sol.t`.  
`molefraction(bsol,t)` will the mole fractions of all species at time `t`.  
`molefraction(bsol,name,t)` will give the mole fraction of the species with name `name` at time `t`.  

### ROPs

Rates of production/consumption (associated with each reaction and each species at a given time)
can be retrived using the `rops` function.  
`rops(bsol,t)` will the matrix of rops for all species and all reactions at time `t`.  
`rops(bsol,name,t)` will give the array of rops for the species with name `name` at time `t`.  

### Concentration Sensitivities

Concentration sensitivities to rate coefficients and species gibbs free energies can be retrieved
using the `getconcentrationsensitivity(bsol,numerator,denominator,t)` function.  
Here `bsol` denotes the Simulation object, and `t` denotes the time.  The output is the sensitivity of
the `numerator` to the `denominator`.  For species concentration and species thermo this is the name of the
species for reactions this is the index of the reaction.  

### Rates

The function `rates` can be used to calculate the rates of all reactions at specific time points.  
`rates(bsol,t)` will give the array of reaction rates at time `t`
while `rates(bsol;ts=ts)` will give a matrix of reaction rates at all times in `ts`.  
Note that `ts` defaults to `bsol.sol.t`.  

### Adjoint Sensitivities

Sensitivity values to a target species or thermodynamic variable can be computed from a `Simulation` or `SystemSimulation` object using the `getadjointsensitivities(bsol::Q,target::String,solver::W;sensalg::W2=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)),abstol::Float64=1e-6,reltol::Float64=1e-3,normalize=true,kwargs...)` function. This computes the sensitivity with respect to the target at the final time point of bsol. It uses `solver`, `sensalg`, `abstol`, and `reltol` for the adjoint solve and by default will give the normalized sensitivity values (note these are molar sensitivities, concentration sensitivities can't be computed from a single adjoint pass). This is usually much faster than forward sensitivities.

### Transitory Sensitivities

Transitory sensitivity values can be computed using several different algorithms. `transitorysensitivitiesfullexact(sim::Simulation,t;tau=NaN,
        normalized=true,solver=Sundials.CVODE_BDF(linear_solver=:GMRES),
        abstol=1e-16,reltol=1e-6)` gives you the exact full matrix of transitory sensitivities using the forward sensitivity algorithm, while `transitorysensitivitiesfulltrapezoidal(sim,t;tau=NaN,normalized=true)` gives the approximate full matrix of transitory sensitivities using the trapezoidal method. `transitorysensitivitiesparamexact` and `transitorysensitivitiesparamtrapezoidal` are available for computing a single column of the matrix (with respect to a single parameter)). Lastly `transitorysensitivitiesadjointexact(sim::Simulation,t,name;tau=NaN,
        normalized=true,solver=Sundials.CVODE_BDF(),sensalg=InterpolatingAdjoint(),
        abstol=1e-16,reltol=1e-6)` is available for computing a single row of the matrix (sensitivity to a specified species with respect to all parameters). The adjoint algorithm is jacobian free if `tau` is specified and the solver is jacobian free. 

### Other Useful Properties

Please let us know on our Github issues page if we're missing any important
property calculators.  

## Crash Analysis

When thermodynamics and kinetics are poorly assigned mechanisms can become too stiff to simulate causing the solver to crash. In these
 cases it is very important to be able to efficiently identify potential offending thermochemistry and/or kinetics. Our crash analysis tool analyzes 
 the reaction and species fluxes to identify NaN quantities and quantities that are unusually large. A crash report can be generated from the associated `Simulation` or `SystemSimulation` object with `printcrashanalysis(analyzecrash(sim;tol=tol))`. In general, the correct `tol` depends on how much faster the offending chemistry is than the real chemistry. Rather than use the default tolerance one should adjust `tol` to achieve a reasonable amount of possible offending thermochemistry and kinetics, adjusting up to reduce the amount identified while decreasing `tol` will cause
the amount identified to increase. 

## Plotting

### Plotting Mole Fractions

Mole fractions can be plotting using the `plotmolefractions` function
`plotmolefractions(bsol, tf; t0=1e-15,N=1000,tol=0.01)` plots all species with molefraction greater than
`tol` at `N` logarithmically spaced time points between `t0` and `tf`.  
`plotmolefractions(bsol; tol=0.01)` plots all species with molefraction greater than
`tol` at the points in `bsol.sol.t`.  
`plotmolefractions(bsol,spcnames)` plots all the species with names in `spcnames` at the points in `bsol.sol.t`.  

### Plotting Forward Sensitivities

Sensitivities (normalized molar sensitivities) can be plotted using the `plotmaxthermoforwardsensitivity` and `plotmaxrateforwardsensitivity` functions.  
Both of these follow the same format:  
`plotmaxthermoforwardsensitivity(bsol, spcname; N=0, tol= 1e-2)`
`spcname` corresponds to the species sensitivities are being calculated for, `N` is the maximum number of
sensitive species/reactions plotted (0 corresponds to all of them), sensitive species/reactions with sensitivities
less than `tol` are not included in the plot.  Note that the thermo sensitivities are given in mol/kcal while the rate
sensitivities are fully non-dimensionalized (as displayed on the plots).  

### Plotting Adjoint Sensitivities

Adjoint sensitivities can be plotted using the `plotthermoadjointsensitivities(bsol::Y,name::X,dps::Z;N=0,tol=0.01)` and `plotrateadjointsensitivities(bsol::Y,name::X,dps::Z;N=0,tol=0.01)` functions where `dps` is the normalized adjoint sensitivity values.

### Plotting ROPs

ROPs can be plotted with the `plotrops` function.  
`plotrops(bsol,name,t;N=0,tol=0.01)` will plot the ROPs for species with name `name`
at time `t` including at most `N` reactions and not including reactions with absolute rates less than
tol * the largest absolute rate.  This is a bar plot comparing the reactions that contribute the most to
the production and loss of the species

The rops can be plotted with respect to time on a line plot using
`plotrops(bsol,name;rxnrates=Array{Float64,1}(),ts=Array{Float64,1}(),tol=0.05)`
in this case `rxnrates` corresponds to the matrix of reaction rates at each time point
(can be expensive to compute so if available can be reused), `ts` correpsonds to a set of
time points to plot at (otherwise defaults to `bsol.sol.t`), any reaction with flux smaller
than `tol` * the largest absolute rate at every time point is excluded from the plot.  

The analogous functions `plotradicalrops(bsol,t;N=0,tol=0.01)` and `plotradicalrops(bsol;rxnrates=Array{Float64,1}(),ts=Array{Float64,1}(),tol=0.05)` are available for plotting the rops for the sum of all radicals.

### Plotting Transitory Sensitivities

Transitory sensitivities and be combusted and plotted using `plotrxntransitorysensitivities(bsol,name,t;dSdt=nothing,tau=nothing,tol=1e-3,N=0,rxntol=1e-6)`
and `plotthermotransitorysensitivities(bsol,name,t;dSdt=nothing,tau=nothing,tol=1e-3,N=0)` where dSdt contains the transitory sensitivity values otherwise these values will be computed automatically using the trapezoidal method using the input `tau` or automatically selecting `tau` and the `rxntol` value. At most `N` reactions are included in the plot and the plot will not include reactions with absolute transitory sensitivities less than `tol` times the largest absolute value.

### Plotting Time Scales

The timescale distribution of a simulation at a point can be plotted using `plottimescales(sim,t;taumax=1e6,taumin=1e-18,taures=10.0^0.5,usediag=true)` or `plottimescales(Jy;taumax=1e6,taumin=1e-18,taures=10.0^0.5,usediag=true)` where `taumax`, `taumin` and `taures` control the bins. If `usediag=true` it will simply determine the timescales using the diagonal of the Jacobian otherwise it will compute and use the eigenvalues. In general we've observed no significant differences in distributions generated using the diagonal values vs the eigenvalues although there may be significant differences when the mechanism is small.

### Other Plots

While we are trying to build up a library of plotting functions that make mechanism analysis easier
we may not have the one you need.  However, the tools in the Useful Properties section should be enough
most of the time. We're happy to provide guidance on our Github issues page and add plotting functions you find useful.  

## Flux Diagrams

RMS generates flux diagrams with molecular images (or names when images aren't available) using the `getfluxdiagram` function:  
```
getfluxdiagram(bsol,t;centralspecieslist=Array{String,1}(),superimpose=false,
    maximumnodecount=50, maximumedgecount=50, concentrationtol=1e-6, speciesratetolerance=1e-6,
    maximumnodepenwidth=10.0,maximumedgepenwidth=10.0,radius=1,centralreactioncount=-1,outputdirectory="fluxdiagrams")
```
This generates a flux diagram at the time point `t`.  `centralspecieslist` denotes a list of species that must be included in
the diagram, `maximumnodecount` denotes the maximum number of species in the diagram, `maximumedgecount` denotes the maximum
number of connections between species, `concentrationtol` denotes the lowest fractional concentration to show in the diagram,
`speciesratetolerance` denotes the lowest fraction species rate to show in the diagram, `maximumNodePenWidth` denotes the
thickness of the border around a node a maximum concentration.  `maximumedgepenwidth` denotes the thickness of the edge at
maximum species rate, `radius` is the graph radius plotted around a central species.  

## Other Useful Functionality

`spcindex(bsol,name)` will give you the index of the species with name `name`.  
