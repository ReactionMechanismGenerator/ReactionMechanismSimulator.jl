# Automatic Mechanism Analysis

## Running Automatic Mechanism Analysis

Automatic mechanism analysis can be run for a single species at a single time point using the function `analyzespc(sim,spcname,t;N=10,tol=1e-3,branchthreshold=0.9,
        pathbranchthreshold=0.2,branchtol=1e-2,steptol=1e-2,
        transitorysensitivitymethod=transitorysensitivitiesfulltrapezoidal,
        eliminate=true
        )`. This returns an array of ReactionAnalysis objects corresponding to each reaction found in the analysis (based on `N` and `tol` choices).

## ReactionAnalysis Object

The ReactionAnalysis object has seven attributes, `branchings` the array of potentially important Branching objects, `paths` the array of potentially important ReactionPath objects, `radprodlossfract` the fraction of production (+) or loss (-) of radicals that this reaction accounts for, `spcind` the index of the target species, `spcname` the name of the target species, `rxnind` the index of the reaction and `sens` the transitory sensitivity value. Can be dumpped to a string report by `getrxnanalysisstring(sim,ra;branchingcutoff=1e-2,radbranchfract=0.01)` or simply printed with `printrxnanalysis(sim,ra;branchingcutoff=1e-2,radbranchfract=0.01)` where the `branchingcutoff` is the fraction of the branching at which reactions will no longer show up as part of a branching and `radbranchfract` is the fraction of radical production/loss above which the value is displayed as important.

## Branching Object

The Branching object has three attributes, `spcind` the index of the target species, `rxninds` the array of the indices of the reactions in order of branching fraction, `branchingratios` the fraction of the branching accounted for by each reaction. Can be used to generate a flux diagram with `getfluxdiagram(bsol,t,b::Branching; branchtol=1.0e-2, kwargs...)` where `branchtol` is the fraction of the branching at which the product of a reaction will not be included in the flux diagram. Note this flux diagram includes all reaction between species part of the important branching reactions and not just the branching reactions themselves.

## ReactionPath Object  

The ReactionPath object has six important attributes, `forward` indices whether the path was generated following the flux forward from the target species or backwards from the target species, `spcsinds` is the array of the species indices followed in order along the flux path, `rxninds` is the array of reactions that connect these `spcsinds`, `spcind` is the index of the target species, `branchfracts` is the branching fraction of the reaction in `rxninds` for the species in `spcsinds`, `branchfract` is the fraction of the flux following from the start of the path to the end. Can be used to generate a flux diagram with `getfluxdiagram(bsol,t,rp::ReactionPath; radius=0, kwargs...)`. Note this flux diagram includes all reactions between species that are parts of the ReactionPath and not just the reaction path itself. 
