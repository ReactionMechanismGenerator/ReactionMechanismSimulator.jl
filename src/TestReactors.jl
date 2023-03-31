using Test
using SciMLBase
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

@testset "Test liquid phase Constant T Constant V reactor with volumetric flow rate inlet and outlet simulation" begin
    #Constant T and V Ideal Dilute Liquid
    initialconds = Dict(["T"=>450.0,"V"=>1.0e-6*1e6,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6]) #Set simulation Initial Temp and Pressure
    domain,y0,p = ConstantTVDomain(phase=liq,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

    inlet = VolumetricFlowRateInlet(domain,Dict(["T"=>450.0,"P"=>1.e5,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6]),x->1.0) #set the inlet flow rate and conditions
    outlet = VolumetricFlowRateOutlet(domain,x->1.0) #set the outlet flow rate
    interfaces = [inlet,outlet]
    react = Reactor(domain,y0,(0.0,140000.01),interfaces;p=p) #Create the reactor object
    
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

    @testset "Test liquid phase Constant T Constant V reactor simulation with kLAkHCondensationEvaporationWithReservoir" begin
        phaseDict = readinput("../src/testing/constantkLAkH.rms")

        spcs = phaseDict["phase"]["Species"]
        rxns = phaseDict["phase"]["Reactions"]
        solvent = phaseDict["Solvents"][1]

        liq = IdealDiluteSolution(spcs,rxns,solvent;name="phase",diffusionlimited=true)

        initialconds = Dict(["T"=>450.0,"V"=>1.0e-6*1e6,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6])
        domain,y0,p = ConstantTVDomain(phase=liq,initialconds=initialconds)
        conds = Dict(["T"=>450.0,"P"=>1.e5,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6])
        interfaces = [kLAkHCondensationEvaporationWithReservoir(domain,conds)]
        react = Reactor(domain,y0,(0.0,140000.01),interfaces;p=p)

        sol1 = solve(react.ode,react.recommendedsolver,abstol=1e-18,reltol=1e-6);

        phaseDict = readinput("../src/testing/TdependentkLAkH.rms")
        spcs = phaseDict["phase"]["Species"]
        rxns = phaseDict["phase"]["Reactions"]
        solvent = phaseDict["Solvents"][1]
        liq = IdealDiluteSolution(spcs,rxns,solvent;name="phase",diffusionlimited=true)

        initialconds = Dict(["T"=>450.0,"V"=>1.0e-6*1e6,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6])
        domain,y0,p = ConstantTVDomain(phase=liq,initialconds=initialconds)
        conds = Dict(["T"=>450.0,"P"=>1.e5,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6])
        interfaces = [kLAkHCondensationEvaporationWithReservoir(domain,conds)]
        react = Reactor(domain,y0,(0.0,140000.01),interfaces;p=p) #Create the reactor object

        sol2 = solve(react.ode,react.recommendedsolver,abstol=1e-18,reltol=1e-6);
        
        spcnames = getfield.(liq.species,:name)
        octaneind = findfirst(isequal("octane"),spcnames)
        @test sol1(32977.61568) ≈ sol2(32977.61568)

        end;

    @testset "Test vapor-liquid phase multi-domain reactor simulation with VaporLiquidMassTransferInternalInterfaceConstantT and VolumeMaintainingOutlet interface" begin
        input_file = "../src/testing/TdependentkLAkH.rms"
        phaseDict = readinput(input_file)
        liqspcs = phaseDict["phase"]["Species"]
        liqrxns = phaseDict["phase"]["Reactions"]
        solvent = phaseDict["Solvents"][1]
        gasspcs = liqspcs
        gasrxns = []
        liqspcnames = getfield.(liqspcs,:name)
        gasspcnames = getfield.(gasspcs,:name)

        gas = IdealGas(gasspcs,gasrxns;name="gas"); 
        liq = IdealDiluteSolution(liqspcs,liqrxns,solvent;name="liq",diffusionlimited=true)

        Vliq = 1.0
        Vgas = 1.0
        T = 25 + 273.15
        octaneconc = 6478
        tf = 3600*24

        initialconds = Dict("octane"=>octaneconc*Vliq,"T"=>T,"V"=>Vliq)
        domainliq,y0liq,pliq = ConstantTVDomain(phase=liq,initialconds=initialconds);

        P = 101300.0
        N2N = P*Vgas/R/T*0.8
        oxygenN = P*Vgas/R/T*0.2
        initialconds = Dict("N2"=>N2N,"oxygen"=>oxygenN,"T"=>T,"P"=>P)
        domaingas,y0gas,pgas = ConstantTPDomain(phase=gas,initialconds=initialconds);
        initialconds = Dict("N2"=>N2N,"oxygen"=>oxygenN,"T"=>T,"P"=>P)
        inletgas = Inlet(domaingas,initialconds,x->42)
        outletgas = VolumeMaintainingOutlet(domaingas)

        ignoremasstransferspcnames = ["octane"] #only allowing octane to be consumed via reactions, to avoid the liquid phase to dry out completely
        vl,pinter = VaporLiquidMassTransferInternalInterfaceConstantT(domaingas,domainliq,ignoremasstransferspcnames);

        domains = (domainliq,domaingas)
        interfaces = [vl,inletgas,outletgas]
        react,y0,p = Reactor(domains,(y0liq,y0gas),(0.0,tf),interfaces,(pliq,pgas,pinter));
        sol = solve(react.ode,react.recommendedsolver,abstol=1e-18,reltol=1e-6);

        name = "oxygen"
        ind = findfirst(x->x==name,liqspcnames)
        @test sol(sol.t[end])[ind] ≈ 0.11758959354431776 rtol=1e-5 #test there are oxygen dissolved into the liquid 

    end

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
solthreaded = threadedsensitivities(Reactor(domain,y0,(0.0,150.11094);p=p,forwardsensitivities=false);
        odekwargs=Dict([:abstol=>1e-21,:reltol=>1e-7]),senskwargs=Dict([:abstol=>1e-21,:reltol=>1e-7]))
simthreaded = Simulation(solthreaded,domain)

x,dp = extract_local_sensitivities(sol2,150.11094);
xth,dpth = extract_local_sensitivities(solthreaded,150.11094);

ind = findfirst(isequal("H2"),sim2.names)
dpvs = [v[ind] for v in dp]
dpvsth = [v[ind] for v in dpth]
dpvs[length(domain.phase.species)+1:end] .*= domain.p[length(domain.phase.species)+1:end]
dpvsth[length(domain.phase.species)+1:end] .*= domain.p[length(domain.phase.species)+1:end]
dpvs ./= sol2(150.11094)[ind]
dpvsth ./= solthreaded(150.11094)[ind]

rerr = (dpvs .- dps')./dpvs
rerrth = (dpvsth .- dps')./dpvsth
rerr = [isinf(x) ? 0.0 : x for x in rerr]
rerrth = [isinf(x) ? 0.0 : x for x in rerrth]
@test all((abs.(rerr) .> 1e-1).==false)
@test all((abs.(rerrth) .> 1e-1).==false)
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
(dSdt,tau) = transitorysensitivitiesfulltrapezoidal(sim,0.01)
dSdttrape = dSdt[ethaneind,:]
dSdttrapetrue = [2.9575870145096627e-7, 4.515275851956882e-11, 4.988521948240434e-11, -2.5799274590483593e-11, -3.4555131774031248e-12, 0.0, -2.1492670519414238e-7, 7.46732761666698e-9, 0.0, 0.0, -2.40138848721267e-7, 3.9219318929084755e-7, -1.4739312254182412e-8, 5.92445424297349e-7, -7.125209636476018e-8, -1.1352622702113613e-7, 7.946413524885843e-8, 6.705472353867515e-12, -4.557442966614505e-10, 7.418526742019254e-9, 8.438347173720314e-11, -1.2366741309434584e-10, -4.225417298096224e-13, 3.365005913661073e-11, -1.0695252156578608e-12, 2.0373377649853602e-16, -0.00044230108462490786, -1.0167340432033213e-8, -2.583036369059463e-13, 5.781731887238174e-15, 0.0, -4.334058036761624e-16, 0.0, 1.3168450957987063e-15, 1.1998309904180759e-6, -5.949349200507399e-6, 3.1444981902056855e-7, 1.8398290487996734e-5, 7.304780256981038e-8, 9.618451445619252e-8, 2.8366497171337888e-8, -1.364696143868719e-6, -8.810140643319245e-11, -3.962757443811278e-12, 5.027540595806555e-9,
-0.00023939911435707088, -0.032852326043979294, -0.03798141508237119, -0.0033779935066780868, 7.105745741537494e-9, 1.412003505611414e-10, 1.0225614778472479e-8, 7.271674880756403e-10, 1.0603171573593312e-5, 6.554625893344657e-10, -0.0009646230043543857, 3.6697977742128043e-10, -2.6016005061672334e-8, 8.244482601003305e-10, 9.871392088202154e-8, 7.087374979644391e-8, 3.7494699026628695e-6, -0.0016360860266714356, 7.46796895065408e-10, -2.1377430063578516e-8, -1.5313395782779889e-6,
-3.44814170913112e-8, -2.0205662027875063e-13, -6.167915597541991e-5, -1.0232310996433263e-8, 8.773160423979728e-13, -8.272836567208742e-12, -1.3849538738942805e-7, -2.8267415377402184e-7, -1.956160665062943e-9, 1.4249212654328193e-10, 1.4025863899596462e-13, 2.6192412351312244e-15, 6.209421075823885e-13, 7.563798137179707e-6, 1.0737707383411033e-17, 5.5720952223799324e-8, 6.707604770379481e-12, 4.083583009271183e-16, -4.9694659616325786e-11, 1.2353451155074604e-10, -2.9313710880962973e-9, 2.5209077489688276e-15, -3.605953179593259e-8, -9.262968242527441e-7, -5.267767400715836e-9, -2.0538584940509884e-12, 2.4582858956997566e-15, -7.055068170841879e-9, -1.4020276284380862e-6, -0.001988526129422606, 5.222914010991184e-6, -0.002389400035195657, 5.892619594146276e-9, -6.593316167800227e-12, 6.798280145273607e-14, -1.2498073552259318e-10, -9.721907208893787e-7, 1.8620274111689625e-10, 2.705888506809334e-7, 2.4326666876274853e-44, 9.630742927954911e-42, 2.4353199117753496e-16, 1.8705658712224538e-14, 1.5145954056505608e-19, -6.5851089037225425e-6, -8.587734090686863e-26, 2.8904881300280135e-9, 9.637715763640123e-14, -6.489237417275701e-12, 3.037266296308308e-10, 0.0016140972380429778, 5.6414646097585933e-5, 1.972981642283126e-11, 2.073029225006578e-8, 1.190672498122864e-7, 2.4661025783165292e-14, -2.8053580290286356e-9, 1.3198518652028483e-6, 2.5680388016794715e-15, 6.7963019095535e-14, -9.311987491323173e-10, 1.3580367145534529e-14, -4.143118184898953e-11,
-1.605904738012102e-11, 3.72207891995559e-5, 4.849224291528469e-6, -0.0016554597624482278, -5.9196933275908064e-8, 1.73270139591348e-9, 8.521863840987097e-9, -2.4562040665363023e-6, -4.934842425177314e-6, -2.5768820349447168e-8, 3.4288679720876376e-7, -1.1635869304274178e-9, 3.650102637910092e-13, 2.568655079627211e-14, -8.203639185708453e-15, 1.3851394699508101e-14, 3.905955571296345e-12, 1.0508113259645755e-5, -3.7536060928465414e-7, 0.00023524067977726278, -0.03235148399206836, 3.843597225108039e-13, -5.452697163403704e-10, 9.281350122766843e-10, 9.181092007716816e-6, 1.1203368555573775e-5, -2.0769630397555456e-5, 3.495769120111611e-7, 1.4873642363083278e-7, 7.540628688138892e-13, 2.667344416697606e-12]
@test all((abs.(dSdttrape.-dSdttrapetrue) .> 1e-4.*abs.(dSdttrapetrue).+1e-16).==false)

out = analyzespc(sim,"ethane",0.01)
s = getrxnanalysisstring(sim,out[1])
@test s == "Analyzing ethane sensitivity to HO2+ethane<=>H2O2+C2H5 at a value of -0.0379814 \n\nKey branching for HO2 \nHO2+HO2<=>O2+H2O2 had a branching ratio of 0.632358 \nHO2+ethane<=>H2O2+C2H5 had a branching ratio of 0.358444 \n\nKey branching for ethane \nOH+ethane<=>H2O+C2H5 had a branching ratio of 0.576145 \nHO2+ethane<=>H2O2+C2H5 had a branching ratio of 0.328068 \nH+ethane<=>H2+C2H5 had a branching ratio of 0.0313989 \nO2+ethane<=>HO2+C2H5 had a branching ratio of 0.0291454 \nCH3+CH3<=>ethane had a branching ratio of 0.0165147 \nCH3+ethane<=>CH4+C2H5 had a branching ratio of 0.0145334 \n\nAssociated key reaction path in ethane loss direction \nHO2+ethane<=>H2O2+C2H5 at path branching of 0.328068 \nHO2+C2H4<=>O2+C2H5 at path branching of 0.919341 \n\nAssociated key reaction path in ethane loss direction \nHO2+ethane<=>H2O2+C2H5 at path branching of 0.328068 \nHO2+C2H4<=>O2+C2H5 at path branching of 0.919341 \n\nAssociated key reaction path in ethane loss direction \nHO2+ethane<=>H2O2+C2H5 at path branching of 0.328068 \nHO2+C2H4<=>O2+C2H5 at path branching of 0.919341 \n\n\n"

phaseDict = readinput("../src/testing/ethane.rms")
spcs = phaseDict["phase"]["Species"]
rxns = phaseDict["phase"]["Reactions"]
ig1 = IdealGas(spcs,[],name="phase1")
ig2 = IdealGas(spcs,rxns,name="phase2")

initialcondsV1 = Dict(["T"=>1000.0,"P"=>2.0e5,"ethane"=>1.0,"Ar"=>1.0,"O2"=>3.5])
domainV1,y0V1,pV1 = ConstantPDomain(phase=ig1,initialconds=initialcondsV1) #Define the domain (encodes how system thermodynamic properties calculated)
initialcondsV2 = Dict(["T"=>1000.0,"P"=>2.0e5,"ethane"=>1.0,"Ar"=>1.0,"O2"=>3.5])
domainV2,y0V2,pV2 = ConstantPDomain(phase=ig2,initialconds=initialcondsV2) #Define the domain (encodes how system thermodynamic properties calculated)

react,y0,p = Reactor((domainV1,domainV2),(y0V1,y0V2),(0.0,1.0),[],(pV1,pV2));
sol = solve(react.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6);
sysim = SystemSimulation(sol,(domainV1,domainV2),[],p);

out = analyzespc(sysim.sims[2],"ethane",0.01)
s = getrxnanalysisstring(sysim.sims[2],out[1])
@test s == "Analyzing ethane sensitivity to HO2+ethane<=>H2O2+C2H5 at a value of -0.0379814 \n\nKey branching for HO2 \nHO2+HO2<=>O2+H2O2 had a branching ratio of 0.632358 \nHO2+ethane<=>H2O2+C2H5 had a branching ratio of 0.358444 \n\nKey branching for ethane \nOH+ethane<=>H2O+C2H5 had a branching ratio of 0.576145 \nHO2+ethane<=>H2O2+C2H5 had a branching ratio of 0.328068 \nH+ethane<=>H2+C2H5 had a branching ratio of 0.0313989 \nO2+ethane<=>HO2+C2H5 had a branching ratio of 0.0291454 \nCH3+CH3<=>ethane had a branching ratio of 0.0165147 \nCH3+ethane<=>CH4+C2H5 had a branching ratio of 0.0145334 \n\nAssociated key reaction path in ethane loss direction \nHO2+ethane<=>H2O2+C2H5 at path branching of 0.328068 \nHO2+C2H4<=>O2+C2H5 at path branching of 0.919341 \n\nAssociated key reaction path in ethane loss direction \nHO2+ethane<=>H2O2+C2H5 at path branching of 0.328068 \nHO2+C2H4<=>O2+C2H5 at path branching of 0.919341 \n\nAssociated key reaction path in ethane loss direction \nHO2+ethane<=>H2O2+C2H5 at path branching of 0.328068 \nHO2+C2H4<=>O2+C2H5 at path branching of 0.919341 \n\n\n"
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

@testset "Multi-domain Liquid-Film ConstantTV and FragmentBasedConstantTrhoDomain Simulation" begin
    path = "../src/testing/minimal_fragment_based_film_growth_model.rms"
    phaseDict = readinput(path);
    liqspcs = phaseDict["liquid"]["Species"];
    liqrxns = phaseDict["liquid"]["Reactions"];
    filmspcs = phaseDict["film"]["Species"];
    filmrxns = phaseDict["film"]["Reactions"];
    interfacerxns = phaseDict[Set(["liquid","film"])]["Reactions"];
    solvent = phaseDict["Solvents"][1];

    liq = IdealDiluteSolution(liqspcs,liqrxns,solvent;name="liquid",diffusionlimited=true);
    film = FragmentBasedIdealFilm(filmspcs,filmrxns;name="film",diffusionlimited=false);

    initialconds = Dict(["T"=>298.0,"V"=>1.0,"1,3-BUTADIENE(L)"=>10000.0,"2-BUTENE(L)"=>1000.0,
                        "CYCLOPENTADIENE(L)"=>100.0]);
    domainliq,y0liq,pliq = ConstantTVDomain(phase=liq,initialconds=initialconds,constantspecies=["1,3-BUTADIENE(L)","2-BUTENE(L)","CYCLOPENTADIENE(L)"]);

    initialconds = Dict(["T"=>298.0, "A"=>1.0, "rho"=>1.0, "mass"=>1.0,
                         "AR"=>0.01, "KR"=>0.01, "AH"=>40000.0, "CDB"=>10000.0,
                         ]);
    domainfilm,y0film,pfilm = FragmentBasedConstantTrhoDomain(phase=film,initialconds=initialconds);

    inter,pinter = FragmentBasedReactiveFilmGrowthInterfaceConstantT(domainfilm,domainliq,interfacerxns);

    react,y0,p = Reactor((domainfilm,domainliq),(y0film,y0liq),(0.0,0.1),(inter,),(pfilm,pliq,pinter));
    
    sol = solve(react.ode,react.recommendedsolver,abstol=1e-20,reltol=1e-6);
    
    ssys = SystemSimulation(sol,(domainfilm,domainliq,),(inter,),p);
    
    @test concentrations(ssys, "CDB", 0.05) ≈ 10000.322038151966 rtol=1e-5
    @test concentrations(ssys, "KR", 0.05) ≈ 1.008673416492296e-9 rtol=1e-5
end

@testset "Test Crash Analysis" begin
 #Define the phase (how species thermodynamic and kinetic properties calculated)
    phaseDict = readinput("../src/testing/superminimal_broken.rms")
    spcs = phaseDict["phase"]["Species"];
    rxns = phaseDict["phase"]["Reactions"];
    ig = IdealGas(spcs,rxns;name="gas");
    initialconds = Dict(["T"=>1000.0,"P"=>1e5,"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp and Pressure
    domain,y0,p = ConstantTPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

    react = Reactor(domain,y0,(0.0,150.11094);p=p) #Create the reactor object
    sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor
    sim = Simulation(sol,domain,[],p);

    dmech = analyzecrash(sim)
    @test length(dmech.rxns) == 1
    @test length(dmech.spcs) == 4
    @test dmech.rxns[1].index == 11
end;
end;