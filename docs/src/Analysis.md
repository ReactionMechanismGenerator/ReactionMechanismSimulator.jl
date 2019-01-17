# Analysis

## The Simulation Object

Since the solver interface is exposed the solution object generated only gives the raw moles
and raw sensitivity values, which typically are not what a user wants.  The Simulation object
combines the solution information and the domain information to calculate much more useful
properties of the solution.  

The Simulation object can be defined:  

```
bsol = Simulation(sol,domain)
```

where `sol` is the ODESolution object output by the DifferentialEquations package and `domain` is the
domain `sol` corresponds to.  

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
while `rates(bsol;ts)` will give a matrix of reaction rates at all times in `ts`.  
Note that `ts` defaults to `bsol.sol.t`.  

### Other Useful Properties

Please let us know on our Github issues page if we're missing any important
property calculators.  

## Plotting

### Plotting Mole Fractions

Mole fractions can be plotting using the `plotmolefractions` function
`plotmolefractions(bsol, tf; t0=1e-15,N=1000,tol=0.01)` plots all species with molefraction greater than
`tol` at `N` logarithmically spaced time points between `t0` and `tf`.  
`plotmolefractions(bsol; tol=0.01)` plots all species with molefraction greater than
`tol` at the points in `bsol.sol.t`.  
`plotmolefractions(bsol,spcnames)` plots all the species with names in `spcnames` at the points in `bsol.sol.t`.  

### Plotting Concentration Sensitivity

Concentration sensitivities can be plotted using the `plotmaxthermosensitivity` and `plotmaxratesensitivity` functions.  
Both of these follow the same format:  
`plotmaxthermosensitivity(bsol, spcname; N=0, tol= 1e-2)`
`spcname` corresponds to the species sensitivities are being calculated for, `N` is the maximum number of
sensitive species/reactions plotted (0 corresponds to all of them), sensitive species/reactions with sensitivities
less than `tol` are not included in the plot.  Note that the thermo sensitivities are given in mol/kcal while the rate
sensitivities are fully non-dimensionalized (as displayed on the plots).  

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
