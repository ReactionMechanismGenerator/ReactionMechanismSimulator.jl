using Test
include("Parse.jl")
include("Reactor.jl")

#Constant T and P Ideal Gas
phaseDict = readinput("testing/superminimal.yml") #load mechanism dictionary
spcs = phaseDict["gas"]["Species"]; #mechanism dictionaries index:  phaseDict[phasename]["Species" or "Reactions"]
rxns = phaseDict["gas"]["Reactions"];

ig = IdealGas(spcs,rxns,name="gas") #Define the phase (how species thermodynamic and kinetic properties calculated)
initialconds = Dict(["T"=>1000.0,"P"=>1e5,"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp and Pressure
state = MolarState(initialconds,ig) #Define the initial state of the system
domain = ConstantTPDomain(state=state,phase=ig) #Define the domain (encodes how system thermodynamic properties calculated)

react = BatchSingleDomainReactor(domain,(0.0,150.11094)) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

spcnames = getfield.(ig.species,:name)
h2ind = findfirst(isequal("H2"),spcnames)
o2ind = findfirst(isequal("O2"),spcnames)
h2oind = findfirst(isequal("H2O"),spcnames)
y = sol(20.44002454)
N = sum(y)
@test y[h2ind]/N ≈ 0.412883111 rtol=1e-5 #from RMG simulator
@test y[o2ind]/N ≈ 0.200419093 rtol=1e-5
@test y[h2oind]/N ≈ 0.386618602 rtol=1e-5
