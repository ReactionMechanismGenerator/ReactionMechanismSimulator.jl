using Test

@testset "Model reduction and simulation of reduced model" begin
	phaseDict = readinput("../src/testing/minimal.rms") 
	spcs = phaseDict["phase"]["Species"];
	rxns = phaseDict["phase"]["Reactions"];
	ig = IdealGas(spcs,rxns,name="phase")

	qssnames = ["[CH2][CH]C=C"]
	isomergroups = [Dict("[CH2]CC=C"=>0.5,"[CH]=CCC"=>0.5)]

	reducedmodelmappings = generateqsscmapping(ig,qssnames,isomergroups)
	
	initialconds = Dict(["T"=>1350.0,"P"=>1.0e5,"ethane"=>1.0]);
	domain,y0,p = ConstantTPDomain(phase=ig,initialconds=initialconds);
	react = Reactor(domain,y0,(0.0,1e-6),reducedmodelmappings,[];p=p);
	sol = solve(react.ode,react.recommendedsolver,abstol=1e-18,reltol=1e-6);
	unlumpedbsol = Simulation(sol,domain,reducedmodelmappings,[],p);

	react = Reactor(domain,y0,(0.0,1e-6);p=p);
	sol = solve(react.ode,react.recommendedsolver,abstol=1e-18,reltol=1e-6);
	bsol = Simulation(sol,domain,[],p);

	spcnames = getfield.(spcs,:name)
	ethaneind = findfirst(isequal("ethane"),spcnames)
	@test unlumpedbsol.sol(1e-6)[ethaneind] ≈ bsol.sol(1e-6)[ethaneind]
	
	isomer1_ind = findfirst(isequal("[CH2]CC=C"),spcnames)
	isomer2_ind = findfirst(isequal("[CH]=CCC"),spcnames)
	
	@test unlumpedbsol.sol(1e-6)[isomer1_ind] ≈ unlumpedbsol.sol(1e-6)[isomer2_ind]
end