using Test
using SciMLBase
using Sundials

@testset "Species Selection" begin
    outcore = readinput("../src/testing/minimal.rms")
    outedge = readinput("../src/testing/minimal_edge.rms")
    corespcs = outcore["phase"]["Species"]
    edgeishspcs = outedge["phase"]["Species"]
    corerxns = outcore["phase"]["Reactions"]
    edgeishrxns = outedge["phase"]["Reactions"]
    
    coreedgespcs = copy(corespcs)
    index = length(coreedgespcs)
    for spc in edgeishspcs
        if nothing === findfirst(x->x.name==spc.name,coreedgespcs)
            newspc = Species(;name=spc.name,index=index+1,inchi=spc.inchi,smiles=spc.smiles,
                adjlist=spc.adjlist,thermo=spc.thermo,atomnums=spc.atomnums,diffusion=spc.diffusion,
                radius=spc.radius,radicalelectrons=spc.radicalelectrons,molecularweight=spc.molecularweight)
            index += 1
            push!(coreedgespcs,newspc)
        end
    end
    coreedgespcsnames = getfield.(coreedgespcs,:name)

    coreedgerxns = copy(corerxns)
    index = length(coreedgerxns)
    for rxn in edgeishrxns
        out = findfirst(x->getrxnstr(x)==getrxnstr(rxn),coreedgerxns)
        if nothing === out
            inds = findall(x->getrxnstr(x)==getrxnstr(rxn),edgeishrxns)
            for ind in inds
                rxnout = edgeishrxns[ind]
                reactants = [coreedgespcs[findfirst(x->spc.name==x,coreedgespcsnames)] for spc in rxnout.reactants]
                products = [coreedgespcs[findfirst(x->spc.name==x,coreedgespcsnames)] for spc in rxnout.products]
                reactantinds = [findfirst(x->spc.name==x,coreedgespcsnames) for spc in rxnout.reactants]
                productinds = [findfirst(x->spc.name==x,coreedgespcsnames) for spc in rxnout.products]
                newrxn = ElementaryReaction(;index=index+1,reactants=reactants,reactantinds=reactantinds,products=products,
                    productinds=productinds,kinetics=rxnout.kinetics,electronchange=rxnout.electronchange,
                    radicalchange=rxnout.radicalchange,reversible=rxnout.reversible,pairs=rxnout.pairs)
                push!(coreedgerxns,newrxn)
            end
        end
    end
    
    coregas = IdealGas(corespcs,corerxns);
    coreedgegas = IdealGas(coreedgespcs,coreedgerxns);
    
    initialconds = Dict(["T"=>1350.0,"P"=>1.0e5,"ethane"=>1.0]);
    spc = coregas.species[5] #ethane
    termination = [TerminationConversion(spc,0.9),TerminationTime(1e6)];
    coredomain,y0,corep = ConstantTPDomain(phase=coregas,initialconds=initialconds);
    react = Reactor(coredomain,y0,(0.0,1e6);p=corep);
    coreedgedomain,coreedgey0,coreedgep = ConstantTPDomain(phase=coreedgegas,initialconds=initialconds);
    reactedge = Reactor(coreedgedomain,coreedgey0,(0.0,1e6);p=coreedgep);
    (terminated,resurrected,invalidobjects,unimolecularthreshold,bimolecularthreshold,
    trimolecularthreshold,maxedgespeciesrateratios) = selectobjects(react,reactedge,coreedgedomain,[],coredomain,
        [],corep,coreedgep,0.03,0.03,false,true,5,0.005,1.0,1.0,true,termination,1.0e8,Dict(),20)
    @test terminated == false
    @test invalidobjects[1].name == "[CH2]CCC"
    @test unimolecularthreshold[5] == true
    @test maxedgespeciesrateratios[5] ≈ 0.008333897672002308 rtol=1e-3
    
end;