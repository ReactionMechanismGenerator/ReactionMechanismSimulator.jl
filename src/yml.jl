function convertchemkin2yml(chemkinpath;spcdictpath="",output="chem.rms")
    if spcdictpath != ""
        spcs,rxns = chemkin.loadChemkinFile(chemkinpath,dictionaryPath=spcdictpath)
    else
        spcs,rxns = chemkin.loadChemkinFile(chemkinpath)
    end
    writeyml(spcs,rxns;path=output)
end

function writeyml(spcs,rxns;path="chem.rms")
    D = getmechdict(spcs,rxns)
    yaml.dump(D,stream=pybuiltin("file")(path,"w"))
end

function getmechdict(spcs,rxns)
    D = Dict([])
    D["Units"] = Dict([])
    D["Reactions"] = []
    D["Phases"] = [Dict([])]
    D["Phases"][1]["name"] = "phase"
    D["Phases"][1]["Species"] = [obj2dict(x,spcs) for x in spcs]
    D["Reactions"] = [obj2dict(x,spcs) for x in rxns]
    return D
end

function getradicals(obj::T) where {T}
    sm = obj.molecule[1].toSMILES()
    if sm == "[O][O]"
        return 0
    else
        return obj.molecule[1].multiplicity-1
    end
end

function obj2dict(obj,spcs;label="solvent")
    D = Dict([])
    if pybuiltin("isinstance")(obj,species.Species)
        D["name"] = obj.label
        D["type"] = "Species"
        if length(obj.molecule) == 0
            println(obj)
            println(obj.label)
        end
        D["smiles"] = obj.molecule[1].toSMILES()
        D["thermo"] = obj2dict(obj.thermo,spcs)
        if D["smiles"] != "[O][O]"
            D["radicalelectrons"] = obj.molecule[1].multiplicity-1
        else
            D["radicalelectrons"] = 0
        end
    elseif pybuiltin("isinstance")(obj,nasa.NASA)
        D["polys"] = [obj2dict(k,spcs) for k in obj.polynomials]
        D["type"] = "NASA"
    elseif pybuiltin("isinstance")(obj,nasa.NASAPolynomial)
        D["type"] = "NASApolynomial"
        D["coefs"] = PyVector(obj.coeffs)
        D["Tmax"] = obj.Tmax.value_si
        D["Tmin"] = obj.Tmin.value_si
    elseif pybuiltin("isinstance")(obj,reaction.Reaction)
        D["reactants"] = [x.label for x in obj.reactants]
        D["products"] = [x.label for x in obj.products]
        D["kinetics"] = obj2dict(obj.kinetics,spcs)
        D["type"] = "ElementaryReaction"
        D["radicalchange"] = sum([getradicals(x) for x in obj.products])-sum([getradicals(x) for x in obj.reactants])
    elseif pybuiltin("isinstance")(obj,arrhenius.Arrhenius)
        D["type"] = "Arrhenius"
        D["A"] = obj.A.value_si
        D["Ea"] = obj.Ea.value_si
        D["n"] = obj.n.value_si
    elseif pybuiltin("isinstance")(obj,arrhenius.PDepArrhenius)
        D["type"] = "PdepArrhenius"
        D["Ps"] = PyVector(obj.pressures.value_si)
        D["arrs"] = [obj2dict(x,spcs) for x in obj.arrhenius]
    elseif pybuiltin("isinstance")(obj,arrhenius.MultiArrhenius)
        D["type"] = "MultiArrhenius"
        D["arrs"] = [obj2dict(x,spcs) for x in obj.arrhenius]
    elseif pybuiltin("isinstance")(obj,arrhenius.MultiPDepArrhenius)
        D["type"] = "MultiPdepArrhenius"
        D["parrs"] = [obj2dict(x,spcs) for x in obj.arrhenius]
    elseif pybuiltin("isinstance")(obj,falloff.ThirdBody)
        D["type"] = "ThirdBody"
        D["arr"] = obj2dict(obj.arrheniusLow,spcs)
        D["efficiencies"] = Dict([spcs[i].label=>float(val) for (i,val) in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1])
    elseif pybuiltin("isinstance")(obj,falloff.Lindemann)
        D["type"] = "Lindemann"
        D["arrhigh"] = obj2dict(obj.arrheniusHigh,spcs)
        D["arrlow"] = obj2dict(obj.arrheniusLow,spcs)
        D["efficiencies"] = Dict([spcs[i].label=>float(val) for (i,val) in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1])
    elseif pybuiltin("isinstance")(obj,falloff.Troe)
        D["type"] = "Troe"
        D["arrhigh"] = obj2dict(obj.arrheniusHigh,spcs)
        D["arrlow"] = obj2dict(obj.arrheniusLow,spcs)
        D["efficiencies"] = Dict([spcs[i].label=>float(val) for (i,val) in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1])
        D["a"] = obj.alpha
        D["T1"] = obj.T1.value_si
        if !isa(obj.T2,Nothing)
            D["T2"] = obj.T2.value_si
        else
            D["T2"] = 0.0
        end
        D["T3"] = obj.T3.value_si
    elseif pybuiltin("isinstance")(obj,chebyshev.Chebyshev)
        D["type"] = "Chebyshev"
        n,m = size(obj.coeffs.value_si)
        D["coefs"] = PyVector([PyVector(obj.coeffs.value_si[i,:]) for i in 1:n])
        D["Tmin"] = obj.Tmin.value_si
        D["Tmax"] = obj.Tmax.value_si
        D["Pmin"] = obj.Pmin.value_si
        D["Pmax"] = obj.Pmax.value_si
    elseif pybuiltin("isinstance")(obj,wilhoit.Wilhoit)
        D["type"] = "Wilhoit"
        D["coefs"] = [obj.a0,obj.a1,obj.a2,obj.a3]
        D["Cp0"] = obj.Cp0.value_si
        D["Cpinf"] = obj.CpInf.value_si
        D["H0"] = obj.H0.value_si
        D["S0"] = obj.S0.value_si
        D["B"] = obj.B.value_si
    elseif pybuiltin("isinstance")(obj,solvation.SolventData)
        D["type"] = "Solvent"
        D["name"] = label
        dsub = dict()
        dsub["type"] = "RiedelViscosity"
        dsub["A"] = obj.A.value_si
        dsub["B"] = obj.B.value_si
        dsub["C"] = obj.C.value_si
        dsub["D"] = obj.D.value_si
        dsub["E"] = obj.E.value_si
        D["mu"] = dsub
    elseif pybuiltin("is")(obj,pybuiltin("None"))
        return pybuiltin("None")
    else
        throw(error("object $obj not understood"))
    end
    return D
end
