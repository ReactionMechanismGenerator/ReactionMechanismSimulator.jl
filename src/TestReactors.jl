using Test
using DiffEqBase
using Sundials

@testset "Test Reactors" begin

phaseDict = readinput("../src/testing/liquid_phase.rms")
spcs = phaseDict["phase"]["Species"]; #mechanism dictionaries index:  phaseDict[phasename]["Species" or "Reactions"]
rxns = phaseDict["phase"]["Reactions"];
solv = phaseDict["Solvents"][1];
liq = IdealDiluteSolution(spcs,rxns,solv;name="phase",diffusionlimited=true) #Define the phase (how species thermodynamic and kinetic properties calculated)

@testset "Test liquid phase Constant T Constant V reactor simulation" begin
#Constant T and V Ideal Dilute Liquid
initialconds = Dict(["T"=>450.0,"V"=>1.0e-6*1e6,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6]) #Set simulation Initial Temp and Pressure
domain,y0,p = ConstantTVDomain(phase=liq,initialconds=initialconds,constantspecies=["oxygen"]) #Define the domain (encodes how system thermodynamic properties calculated)
react = Reactor(domain,y0,(0.0,140000.01);p=p) #Create the reactor object

sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-8); #solve the ode associated with the reactor

spcnames = getfield.(liq.species,:name)
octaneind = findfirst(isequal("octane"),spcnames)
y = sol(32977.61568)
@test y[octaneind]/sum(y) ≈ 0.461599061 rtol=3e-2 #from RMG simulator I believe the slight difference is due to better calculation of diffusion limits in RMS

#analytic jacobian vs. ForwardDiff jacobian
t=32977.61568;
y=sol(t)
ja=jacobiany(y,p,t,domain,[],nothing);
j=jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

jpa=jacobianp(y,p,t,domain,[],nothing);
jp=jacobianpforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)
end;

@testset "Test liquid phase Parametrized T Constant V reactor jacobian" begin
#Parametrized T constant V Ideal Dilute Liquid
initialconds = Dict(["ts"=>[0.,600.,1200.],"T"=>[450.0,490.,500.],"V"=>1.0e-6*1e6,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6])
domain,y0,p = ParametrizedTConstantVDomain(phase=liq,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,140000.01),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-8); #solve the ode associated with the reactor

#analytic jacobian vs. ForwardDiff jacobian
t=32977.61568;
y=sol(t)
ja=jacobiany(y,p,t,domain,[],nothing);
j=jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

jpa=jacobianp(y,p,t,domain,[],nothing);
jp=jacobianpforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)
end;

@testset "Test ConstantTAPhi Reactor Simulation" begin
    phaseDict = readinput("../src/testing/ORR.rms")
    spcs = phaseDict["phase"]["Species"]; #mechanism dictionaries index:  phaseDict[phasename]["Species" or "Reactions"]
    rxns = phaseDict["phase"]["Reactions"];
    AdivV = 1.0
    is = IdealSurface(spcs,rxns,2.486e-5;name="surf");
    initialconds = Dict(["T"=>298.0,"A"=>1.0,"Phi"=>-0.5,"A*"=>0.5,"B*"=>0.5,"OHA"=>0.0,"O2"=>0.25*AdivV,"H+"=>1.81e-6*AdivV,"H2O2aq"=>0.0*AdivV,"H2O"=>5.55e4*AdivV])
    domain,y0,p =ConstantTAPhiDomain(phase=is,initialconds=initialconds;sensitivity=false,constantspecies=["H+","O2","H2O2aq","H2O"]);
    react = Reactor(domain,y0,(0.0,1e-5),p=p);
    sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-6);
    sim = Simulation(sol,domain);
    @test concentrations(sim,"OA",1e-6) ≈ 0.00013773350978007822 rtol=1e-5
    @test concentrations(sim,"A*",1e-6) ≈ 1.3559599927418817e-5 rtol=1e-5
    @test concentrations(sim,"OHA",1e-6) ≈ 0.4998487068452803 rtol=1e-5
end

#Use superminimal example to test
phaseDict = readinput("../src/testing/superminimal.rms") #load mechanism dictionary
spcs = phaseDict["phase"]["Species"]; #mechanism dictionaries index:  phaseDict[phasename]["Species" or "Reactions"]
rxns = phaseDict["phase"]["Reactions"];
ig = IdealGas(spcs,rxns,name="phase")

#Constant T and P Ideal Gas
@testset "Test constant T and P reactor simulation" begin
 #Define the phase (how species thermodynamic and kinetic properties calculated)
initialconds = Dict(["T"=>1000.0,"P"=>1e5,"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp and Pressure
domain,y0,p = ConstantTPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,150.11094);p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor
sim = Simulation(sol,domain);

spcnames = getfield.(ig.species,:name)
h2ind = findfirst(isequal("H2"),spcnames)
o2ind = findfirst(isequal("O2"),spcnames)
h2oind = findfirst(isequal("H2O"),spcnames)
y = sol(20.44002454)
N = sim.N(20.44002454)
@test y[h2ind]/N ≈ 0.412883111 rtol=1e-4 #from RMG simulator
@test y[o2ind]/N ≈ 0.200419093 rtol=1e-4
@test y[h2oind]/N ≈ 0.386618602 rtol=1e-4

#analytic jacobian vs. ForwardDiff jacobian
t=20.44002454;
y=sol(t)
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

jpa=jacobianp(y,p,t,domain,[],nothing);
jp=jacobianpforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)

#sensitivities
dps = getadjointsensitivities(sim,"H2",CVODE_BDF(linear_solver=:GMRES);sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)),abstol=1e-16,reltol=1e-6)
react2 = Reactor(domain,y0,(0.0,150.11094);p=p,forwardsensitivities=true)
sol2 = solve(react2.ode,CVODE_BDF(linear_solver=:GMRES),abstol=1e-21,reltol=1e-7); #solve the ode associated with the reactor
sim2 = Simulation(sol2,domain)

x,dp = extract_local_sensitivities(sol2,150.11094);
ind = findfirst(isequal("H2"),sim2.names)
dpvs = [v[ind] for v in dp]
dpvs[length(domain.phase.species)+1:end] .*= domain.p[length(domain.phase.species)+1:end]
dpvs ./= sol2(150.11094)[ind]
rerr = (dpvs .- dps')./dpvs
rerr = [isinf(x) ? 0.0 : x for x in rerr]
@test all((abs.(rerr) .> 1e-1).==false)
end;

#Constant T and P Ideal Gas
@testset "Test constant T and P reactor with interfaces simulation" begin
    #Define the phase (how species thermodynamic and kinetic properties calculated)
   initialconds = Dict(["T"=>1000.0,"P"=>1e5,"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp and Pressure
   domain,y0,p = ConstantTPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

   interfaces = [Inlet(domain,Dict{String,Float64}("H2"=>0.67,"O2"=>0.33,"T"=>1000.0,"P"=>1e5),x->0.001),
                Outlet(domain,x->0.001)]
   
   react = Reactor(domain,y0,(0.0,150.11094),interfaces;p=p) #Create the reactor object
   sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor
   sim = Simulation(sol,domain,interfaces);
   
   #analytic jacobian vs. ForwardDiff jacobian
   t = 20.44002454;
   y = sol(t)
   ja = jacobiany(y,p,t,domain,interfaces,nothing);
   j = jacobianyforwarddiff(y,p,t,domain,interfaces,nothing);
   @test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)
   
   jpa = jacobianp(y,p,t,domain,interfaces,nothing);
   jp = jacobianpforwarddiff(y,p,t,domain,interfaces,nothing);
   @test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)
   
   #sensitivities
   dps = getadjointsensitivities(sim,"H2",CVODE_BDF(linear_solver=:GMRES);sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)),abstol=1e-16,reltol=1e-6)
   react2 = Reactor(domain,y0,(0.0,150.11094),interfaces;p=p,forwardsensitivities=true)
   sol2 = solve(react2.ode,CVODE_BDF(linear_solver=:GMRES),abstol=1e-21,reltol=1e-7); #solve the ode associated with the reactor
   sim2 = Simulation(sol2,domain,interfaces)
   
   x,dp = extract_local_sensitivities(sol2,150.11094);
   ind = findfirst(isequal("H2"),sim2.names)
   dpvs = [v[ind] for v in dp]
   dpvs[length(domain.phase.species)+1:end] .*= domain.p[length(domain.phase.species)+1:end]
   dpvs ./= sol2(150.11094)[ind]
   rerr = (dpvs .- dps')./dpvs
   rerr = [isinf(x) ? 0.0 : x for x in rerr]
   @test all((abs.(rerr) .> 1e-1).==false)
   end;

#Constant V adiabatic Ideal Gas
#uses superminimal.yml mechanism
@testset "Constant volume adiabatic reactor simulation" begin
initialconds = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp and Pressure
domain,y0,p = ConstantVDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,0.101),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

ts = exp.(range(log(1e-15),length=10000,stop=log(0.1)))
IDT = ts[argmax(diff([sol(t)[end] for t in ts]))] #Ignition Delay Time based on argmax(dTdt(t))

@test IDT ≈ 0.038384723436228063 rtol=1e-5 #from Cantera simulation

#analytic jacobian vs. ForwardDiff jacobian
t=0.01;
y = sol(t);
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

jpa=jacobianp(y,p,t,domain,[],nothing);
jp=jacobianpforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)

#sensitivities
react = Reactor(domain,y0,(0.0,0.02),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor
sim = Simulation(sol,domain)
dps = getadjointsensitivities(sim,"H2",CVODE_BDF();sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)),abstol=1e-16,reltol=1e-6)
react2 = Reactor(domain,y0,(0.0,0.02);p=p,forwardsensitivities=true)
sol2 = solve(react2.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor
sim2 = Simulation(sol2,domain)

x,dp = extract_local_sensitivities(sol2,0.02);
ind = findfirst(isequal("H2"),sim2.names)
dpvs = [v[ind] for v in dp]
dpvs[length(domain.phase.species)+1:end] .*= domain.p[length(domain.phase.species)+1:end]
dpvs ./= sol2(0.02)[ind]
rerr = (dpvs .- dps')./dpvs
rerr = [isinf(x) ? 0.0 : x for x in rerr]
@test all((abs.(rerr) .> 2e-1).==false)
end;

#Constant P adiabatic Ideal Gas
#uses ethane.rms mechanism
@testset "Constant pressure adiabatic reactor simulation" begin

phaseDict = readinput("../src/testing/ethane.rms")
spcs = phaseDict["phase"]["Species"]
rxns = phaseDict["phase"]["Reactions"]
ig = IdealGas(spcs,rxns,name="phase")

initialconds = Dict(["T"=>1000.0,"P"=>2.0e5,"ethane"=>1.0,"Ar"=>1.0,"O2"=>3.5]) #Set simulation Initial Temp and Pressure
domain,y0,p = ConstantPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,1.0);p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor
sim = Simulation(sol,domain)

ts = exp.(range(log(1e-15),length=10000,stop=log(0.2)))
IDT = ts[argmax(diff([sol(t)[end] for t in ts]))] #Ignition Delay Time based on argmax(dTdt(t))

@test IDT ≈ 0.07324954954380769 rtol=1e-5

#analytic jacobian vs. ForwardDiff jacobian
t=0.01;
y=sol(t)
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

jpa=jacobianp(y,p,t,domain,[],nothing);
jp=jacobianpforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)

#transitory sensitivities
ethaneind = findfirst(isequal("ethane"),sim.names)
dSdt = transitorysensitivitiesfulltrapezoidal(sim,0.01)[ethaneind,:]
dSdttrapetrue = [0.0, 0.0, 0.0, -0.0, -0.0, 0.0, -0.0, 0.0, 0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, -0.0, 0.0, 0.0, -0.0, 0.0, 0.0, -0.0, -0.0, 0.0, -0.0, 0.0, -0.00044230108462476047, -1.0167340432041434e-8,
-2.5830363690594336e-13, 5.781731887234333e-15, 0.0, -4.33405803676164e-16, 0.0, 1.316845095798582e-15, 1.1998309904157391e-6, -5.949349200495832e-6, 3.144498190200055e-7, 1.8398290487975233e-5, 7.304780256982262e-8, 9.618451445647677e-8, 2.8366497171346418e-8, -1.3646961438660625e-6,
-8.81014064332219e-11, -3.9627574438125974e-12, 5.02754059580806e-9, -0.000239399114357028, -0.032852326043993726, -0.03798141508230999, -0.0033779935066780968, 7.105745741532718e-9, 1.41200350561073e-10, 1.0225614778464501e-8, 7.271674880751787e-10, 1.060317157357307e-5, 6.554625893343905e-10, -0.0009646230043544521, 3.66979777421103e-10, -2.601600506168768e-8, 8.244482600996483e-10, 9.87139208820034e-8, 7.087374979649977e-8, 3.7494699026570754e-6, -0.0016360860266715445, 7.467968950652317e-10,
-2.1377430063571472e-8, -1.5313395782786743e-6, -3.448141709131119e-8, -2.0205662027871653e-13, -6.167915597532084e-5, -1.0232310996457154e-8, 8.773160423995249e-13, -8.272836567202069e-12, -1.3849538738933464e-7, -2.826741537733756e-7, -1.956160665061754e-9, 1.4249212654326557e-10, 1.402586389960076e-13, 2.6192412351324203e-15, 6.209421075830585e-13, 7.563798137184527e-6, 1.0737707383431534e-17, 5.5720952223744194e-8, 6.707604770381715e-12, 4.083583009282063e-16, -4.9694659616291374e-11,
1.2353451155062566e-10, -2.9313710880975173e-9, 2.5209077489733237e-15, -3.6059531795945054e-8, -9.2629682425454e-7, -5.2677674006964645e-9, -2.053858494046743e-12, 2.4582858957008768e-15, -7.0550681708350215e-9, -1.4020276284380951e-6, -0.0019885261294226055, 5.222914011007864e-6, -0.0023894000351949375, 5.892619594144724e-9, -6.593316167799491e-12, 6.798280145259288e-14, -1.2498073552250886e-10, -9.721907208900016e-7, 1.8620274111607517e-10, 2.705888506814654e-7, 2.4326666876319783e-44,
9.630742927972707e-42, 2.435319911780137e-16, 1.8705658290774946e-14, 1.5145954056536096e-19, -6.585108903726748e-6, -8.587734090666396e-26, 2.890488130028411e-9, 9.637715763594455e-14, -6.489237417276642e-12, 3.037266296314503e-10, 0.0016140972380424994, 5.6414646097567624e-5, 1.9729816422841313e-11, 2.0730292250075028e-8, 1.1906724981252052e-7, 2.4661025783105124e-14, -2.805358029028032e-9, 1.3198518652025385e-6, 2.568038801682934e-15, 6.79630190955828e-14, -9.311987491283682e-10,
1.3580367145556702e-14, -4.143118184892303e-11, -1.605904738012428e-11, 3.722078919943519e-5, 4.849224291533013e-6, -0.0016554597624476694, -5.9196933275912326e-8, 1.7327013959131697e-9, 8.521863840979788e-9, -2.456204066534639e-6, -4.934842425166022e-6, -2.576882034943144e-8, 3.4288679720870346e-7, -1.1635869304276698e-9, 3.6501026379166063e-13, 2.568655079625506e-14, -8.203639185708493e-15, 1.3851394699501025e-14, 3.905955571296776e-12, 1.0508113259625332e-5, -3.7536060928405004e-7,
0.00023524067977726568, -0.032351483992083895, 3.8435972250994037e-13, -5.452697238181546e-10, 9.281350122952164e-10, 9.181092007735427e-6, 1.1203368555596499e-5, -2.0769630397597524e-5, 3.4957691201187006e-7, 1.4873642363068299e-7, 7.540628688149077e-13, 2.6673444167019607e-12]
@test all((abs.(dSdt.-dSdttrapetrue) .> 1e-4.*abs.(dSdttrapetrue).+1e-16).==false)
end;

#Parametrized T and P Ideal Gas
#uses superminimal.yml mechanism
@testset "Parametrized T and P reactor jacobian" begin
initialconds = Dict(["T"=>[1000.0,1400.0,1500.0],"ts"=>[0.,100.,200.],"P"=>[1.0e5,1.8e5,2.0e5],"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp/Pressure and Volume (function/array)
domain,y0,p = ParametrizedTPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,0.101),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

#analytic jacobian vs. ForwardDiff jacobian
t=0.1;
y = sol(t);
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

jpa=jacobianp(y,p,t,domain,[],nothing);
jp=jacobianpforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)
end;

#Parametrized V adiabatic Ideal Gas
#uses superminimal.yml mechanism
@testset "Parametrized volume adiabatic reactor jacobian" begin
initialconds = Dict(["T"=>1000.0,"ts"=>[0.,100.,200.],"V"=>[0.01,0.05,0.1],"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp/Pressure and Volume (function/array)
domain,y0,p = ParametrizedVDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,0.101),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

#analytic jacobian vs. ForwardDiff jacobian
t=0.1;
y = sol(t);
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

jpa=jacobianp(y,p,t,domain,[],nothing);
jp=jacobianpforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)
end;

#Parametrized P adiabatic Ideal Gas
#uses ethane.yml mechanism
@testset "Parametrized pressure adiabatic reactor jacobian" begin
initialconds = Dict(["T"=>1000.0,"ts"=>[0.,100.,200.],"P"=>[2.0e5,3.0e5,5.0e5],"H2"=>0.67,"O2"=>0.33])
domain,y0,p = ParametrizedPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,0.101),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

#analytic jacobian vs. ForwardDiff jacobian
t=0.1;
y = sol(t);
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

jpa=jacobianp(y,p,t,domain,[],nothing);
jp=jacobianpforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(jpa.-jp) .> 1e-4.*abs.(jp).+1e-16).==false)
end;



@testset "Multi-domain ConstantV sensitivity analysis" begin
    phaseDict = readinput("../src/testing/superminimal.rms")
    spcs = phaseDict["phase"]["Species"]
    rxns = phaseDict["phase"]["Reactions"]
    ig = IdealGas(spcs,rxns,name="phase")
    
    initialcondsV = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV,y0V,pV = ConstantVDomain(phase=ig,initialconds=initialcondsV) #Define the domain (encodes how system thermodynamic properties calculated)
    
    reactV = Reactor(domainV,y0V,(0.0,0.037);p=pV) #Create the reactor object
    solV = solve(reactV.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor
    simV = Simulation(solV,domainV)
    
    initialcondsV1 = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV1,y0V1,pV1 = ConstantVDomain(phase=ig,initialconds=initialcondsV1) #Define the domain (encodes how system thermodynamic properties calculated)
    initialcondsV2 = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV2,y0V2,pV2 = ConstantVDomain(phase=ig,initialconds=initialcondsV2) #Define the domain (encodes how system thermodynamic properties calculated)
    
    react,y0,p = Reactor((domainV1,domainV2),(y0V1,y0V2),(0.0,0.037),[],(pV1,pV2));
    sol = solve(react.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6);
    sysim = SystemSimulation(sol,(domainV1,domainV2),[],p);
    
    t = 0.03
    @test sol(t)[1:length(spcs)] ≈ solV(t)[1:end-2] rtol=1e-5
    @test sol(t)[length(spcs)+1:end-4] ≈ solV(t)[1:end-2] rtol=1e-5
    
    dpsV = getadjointsensitivities(simV,"H2",CVODE_BDF();sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)),abstol=1e-16,reltol=1e-6)
    dps = getadjointsensitivities(sysim,sysim.sims[1],"H2",CVODE_BDF();sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false)),abstol=1e-16,reltol=1e-6)
    @test dpsV ≈ dps rtol=1e-4
end;

@testset "Multi-domain ConstantV and ConstantTP simulation" begin
    phaseDict = readinput("../src/testing/superminimal.rms")
    spcs = phaseDict["phase"]["Species"]
    rxns = phaseDict["phase"]["Reactions"]
    ig = IdealGas(spcs,rxns,name="phase")
    
    initialcondsTP = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainTP,y0TP,pTP = ConstantTPDomain(phase=ig,initialconds=initialcondsTP) #Define the domain (encodes how system thermodynamic properties calculated)
    
    reactTP = Reactor(domainTP,y0TP,(0.0,0.04);p=pTP) #Create the reactor object
    solTP = solve(reactTP.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor
    
    initialcondsV = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV,y0V,pV = ConstantVDomain(phase=ig,initialconds=initialcondsV) #Define the domain (encodes how system thermodynamic properties calculated)
    
    reactV = Reactor(domainV,y0V,(0.0,0.04);p=pV) #Create the reactor object
    solV = solve(reactV.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor
    
    initialcondsTP = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainTP,y0TP,pTP = ConstantTPDomain(phase=ig,initialconds=initialcondsTP) #Define the domain (encodes how system thermodynamic properties calculated)
    initialcondsV = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV,y0V,pV = ConstantVDomain(phase=ig,initialconds=initialcondsV) #Define the domain (encodes how system thermodynamic properties calculated)
    
    react,y0,p = Reactor((domainTP,domainV),(y0TP,y0V),(0.0,0.04),[],(pTP,pV));
    sol = solve(react.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6);
    
    t = 0.03
    @test sol(t)[1:length(spcs)] ≈ solTP(t)[1:end-1] rtol=1e-5
    print(length(sol(t)[length(spcs)+1:end-3]))
    print(length(solV(t)[1:end-2]))
    @test sol(t)[length(spcs)+1:end-3] ≈ solV(t)[1:end-2] rtol=1e-5
end;

@testset "Multi-domain Gas-Surface ConstantTP and ConstantTAPhi Simulation" begin
    phaseDict = readinput("../src/testing/ch4o2cat.rms")
    gasspcs = phaseDict["gas"]["Species"]; 
    gasrxns = phaseDict["gas"]["Reactions"];
    surfacespcs = phaseDict["surface"]["Species"]
    surfacerxns = phaseDict["surface"]["Reactions"]
    interfacerxns = phaseDict[Set(["gas","surface"])]["Reactions"];
    
    ig = IdealGas(gasspcs,gasrxns;name="gas"); 
    cat = IdealSurface(surfacespcs,surfacerxns,2.486e-5;name="surface");
    
    initialconds = Dict(["T"=>800.0,"P"=>1.0e5,"O2"=>0.2,"N2"=>0.7,"CH4"=>0.1]); 
    domaingas,y0gas,pgas = ConstantTPDomain(phase=ig,initialconds=initialconds,); 
    
    V = 8.314*800.0/1.0e5
    A = 1.0e5*V
    initialconds = Dict(["T"=>800.0,"A"=>A,"vacantX"=>cat.sitedensity*A]); 
    domaincat,y0cat,pcat = ConstantTAPhiDomain(phase=cat,initialconds=initialconds,); 
    
    inter,pinter = ReactiveInternalInterfaceConstantTPhi(domaingas,domaincat,interfacerxns,800.0,A);
    
    react,y0,p = Reactor((domaingas,domaincat),(y0gas,y0cat),(0.0,0.1),(inter,),(pgas,pcat,pinter));
    
    sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-6);
    
    ssys = SystemSimulation(sol,(domaingas,domaincat,),(inter,),p);
    
    @test concentrations(ssys,"OX",0.5e-5) ≈ 8.033191655902819e-6 rtol=1e-5
    @test molefractions(ssys.sims[1],"H2O",0.5e-5) ≈ 0.10899527627867926 rtol=1e-5
end;

@testset "Multi-domain Gas-Surface ConstantV and ConstantTAPhi Simulation" begin
    phaseDict = readinput("../src/testing/ch4o2cat.rms")
    gasspcs = phaseDict["gas"]["Species"]; 
    gasrxns = phaseDict["gas"]["Reactions"];
    surfacespcs = phaseDict["surface"]["Species"]
    surfacerxns = phaseDict["surface"]["Reactions"]
    interfacerxns = phaseDict[Set(["gas","surface"])]["Reactions"];
    
    ig = IdealGas(gasspcs,gasrxns;name="gas"); 
    cat = IdealSurface(surfacespcs,surfacerxns,2.486e-5;name="surface");
    
    initialconds = Dict(["T"=>800.0,"P"=>1.0e5,"O2"=>0.2,"N2"=>0.7,"CH4"=>0.1]); 
    domaingas,y0gas,pgas = ConstantVDomain(phase=ig,initialconds=initialconds,); 
    
    V = 8.314*800.0/1.0e5
    A = 1.0e5*V
    initialconds = Dict(["T"=>800.0,"A"=>A,"vacantX"=>cat.sitedensity*A]); 
    domaincat,y0cat,pcat = ConstantTAPhiDomain(phase=cat,initialconds=initialconds,); 
    
    inter,pinter = ReactiveInternalInterface(domaingas,domaincat,interfacerxns,A);
    
    react,y0,p = Reactor((domaingas,domaincat),(y0gas,y0cat),(0.0,0.1),(inter,),(pgas,pcat,pinter));
    
    sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-6);
    
    ssys = SystemSimulation(sol,(domaingas,domaincat,),(inter,),p);
    
    @test concentrations(ssys,"OX",0.5e-5) ≈ 1.9165723392283484e-5 rtol=1e-5
    @test molefractions(ssys.sims[1],"H2O",1e-3) ≈ 0.12732345278036702 rtol=1e-5
end;
end;