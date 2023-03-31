using YAML
using Unitful

function convertchemkin2yml(chemkinpath;spcdictpath="",output="chem.rms")
    if spcdictpath != ""
        spcs,rxns = chemkin.load_chemkin_file(chemkinpath,dictionary_path=spcdictpath)
    else
        spcs,rxns = chemkin.load_chemkin_file(chemkinpath)
    end
    D = getmechdictfromchemkin(spcs,rxns)
    writeyml(D;path=output)
end

function writeyml(D;path="chem.rms")
    YAML.write_file(path,D)
end

function getmechdictfromchemkin(spcs,rxns)
    names = [x.label for x in spcs]
    for (i,name) in enumerate(names)
        if count(names.==name) > 1
            names[i] = string(name,"-",count(names.==name))
        end
    end
    D = Dict([])
    D["Units"] = Dict([])
    D["Reactions"] = []
    D["Phases"] = [Dict([])]
    D["Phases"][1]["name"] = "phase"
    D["Phases"][1]["Species"] = [obj2dict(x,spcs,names) for x in spcs]
    D["Reactions"] = [obj2dict(x,spcs,names) for x in rxns]
    return D
end

function convertcanterayml2rmsyml(path;output="chem.rms")
    canterayml = YAML.load_file(path)
    D = getmechdictfromcantera(canterayml)
    writeyml(D;path=output)
end

function getmechdictfromcantera(canterayml)
    unsupportedphases = [
        "binary-solution-tabulated",
        "compound-lattice",
        "constant-density",
        "Debye-Huckel",
        "edge",
        "electron-cloud",
        "fixed-stoichiometry",
        "HMW-electrolyte",
        # "ideal-gas",
        "ideal-molal-solution",
        "ideal-condensed",
        "ideal-solution-VPSS",
        # "ideal-surface",
        "ions-from-neutral-molecule",
        "lattice",
        "liquid-water-IAPWS95",
        "Margules",
        "Maskell-solid-solution",
        "Peng-Robinson",
        "plasma",
        "pure-fluid",
        "Redlich-Kister",
        "Redlich-Kwong",
    ]
    
    spcs = canterayml["species"]
    rxns = canterayml["reactions"]
    phases = canterayml["phases"]
    units = canterayml["units"] #user specified units
    units = Dict(key=>uparse(value) for (key,value) in units) #convert to Unitful units
    
    #prevent duplicate names
    spc_names = [spc["name"] for spc in spcs]
    for (i,name) in enumerate(spc_names)
        if count(spc_names.==name) > 1
            k = 2
            new_name = name
            while new_name in spc_names
                new_name = string(name,"-",k)
            end
            spc_names[i] = new_name
        end
    end
    
    D = Dict([])
    D["Units"] = Dict([]) #will be converting all values to SI
    D["Phases"] = []
    for phase in phases
        if phase["thermo"] in unsupportedphases
            @error "Currently not supporting $(phase["thermo"])"
        end
        phasedict = Dict()
        phasedict["name"] = phase["name"]
        phasedict["Species"] = [canteradict2rmsdict(spc,spcs,spc_names,units,:species) for spc in spcs if spc["name"] in phase["species"]]
        push!(D["Phases"],phasedict)
    end
    D["Reactions"] = [canteradict2rmsdict(rxn,spcs,spc_names,units,:reaction) for rxn in rxns]
    return D
end

function getradicals(obj::T) where {T}
    sm = obj.molecule[1].to_smiles()
    if sm == "[O][O]"
        return 0
    else
        return obj.molecule[1].multiplicity-1
    end
end

function obj2dict(obj,spcs,names;label="solvent")
    D = Dict([])
    if pybuiltin("isinstance")(obj,species.Species)
        D["name"] = names[findall(x->x==obj,spcs)[1]] 
        D["type"] = "Species"
        if length(obj.molecule) == 0
            println(obj)
            println(obj.label)
        end
        D["smiles"] = obj.molecule[1].to_smiles()
        D["thermo"] = obj2dict(obj.thermo,spcs,names)
        if D["smiles"] != "[O][O]"
            D["radicalelectrons"] = obj.molecule[1].multiplicity-1
        else
            D["radicalelectrons"] = 0
        end
    elseif pybuiltin("isinstance")(obj,nasa.NASA)
        D["polys"] = [obj2dict(k,spcs,names) for k in obj.polynomials]
        D["type"] = "NASA"
    elseif pybuiltin("isinstance")(obj,nasa.NASAPolynomial)
        D["type"] = "NASApolynomial"
        D["coefs"] = PyVector(obj.coeffs)
        D["Tmax"] = obj.Tmax.value_si
        D["Tmin"] = obj.Tmin.value_si
    elseif pybuiltin("isinstance")(obj,reaction.Reaction)
        D["reactants"] = [names[findall(y->y==x,spcs)[1]]  for x in obj.reactants]
        D["products"] = [names[findall(y->y==x,spcs)[1]] for x in obj.products]
        D["kinetics"] = obj2dict(obj.kinetics,spcs,names)
        D["type"] = "ElementaryReaction"
        D["radicalchange"] = sum([getradicals(x) for x in obj.products])-sum([getradicals(x) for x in obj.reactants])
        D["reversible"] = obj.reversible
        D["forwardable"] = obj.forwardable
    elseif pybuiltin("isinstance")(obj,arrhenius.Arrhenius)
        D["type"] = "Arrhenius"
        D["A"] = obj.A.value_si
        D["Ea"] = obj.Ea.value_si
        D["n"] = obj.n.value_si
    elseif pybuiltin("isinstance")(obj,arrhenius.PDepArrhenius)
        D["type"] = "PdepArrhenius"
        D["Ps"] = PyVector(obj.pressures.value_si)
        D["arrs"] = [obj2dict(x,spcs,names) for x in obj.arrhenius]
    elseif pybuiltin("isinstance")(obj,arrhenius.MultiArrhenius)
        D["type"] = "MultiArrhenius"
        D["arrs"] = [obj2dict(x,spcs,names) for x in obj.arrhenius]
    elseif pybuiltin("isinstance")(obj,arrhenius.MultiPDepArrhenius)
        D["type"] = "MultiPdepArrhenius"
        D["parrs"] = [obj2dict(x,spcs,names) for x in obj.arrhenius]
    elseif pybuiltin("isinstance")(obj,falloff.ThirdBody)
        D["type"] = "ThirdBody"
        D["arr"] = obj2dict(obj.arrheniusLow,spcs,names)
        D["efficiencies"] = Dict([spcs[i].label=>float(val) for (i,val) in enumerate(obj.get_effective_collider_efficiencies(spcs)) if val != 1])
    elseif pybuiltin("isinstance")(obj,falloff.Lindemann)
        D["type"] = "Lindemann"
        D["arrhigh"] = obj2dict(obj.arrheniusHigh,spcs,names)
        D["arrlow"] = obj2dict(obj.arrheniusLow,spcs,names)
        D["efficiencies"] = Dict([spcs[i].label=>float(val) for (i,val) in enumerate(obj.get_effective_collider_efficiencies(spcs)) if val != 1])
    elseif pybuiltin("isinstance")(obj,falloff.Troe)
        D["type"] = "Troe"
        D["arrhigh"] = obj2dict(obj.arrheniusHigh,spcs,names)
        D["arrlow"] = obj2dict(obj.arrheniusLow,spcs,names)
        D["efficiencies"] = Dict([spcs[i].label=>float(val) for (i,val) in enumerate(obj.get_effective_collider_efficiencies(spcs)) if val != 1])
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

function canteradict2rmsdict(canteradict,spcs,names,units,dict_type;numreactants=nothing,polyindex=nothing)
    D = Dict([])
    if dict_type == :species #species
        D["name"] = names[findall(x->x==canteradict,spcs)[1]] 
        D["type"] = "Species"
        D["smiles"] = ""
        D["thermo"] = canteradict2rmsdict(canteradict["thermo"],spcs,names,units,:thermo)
        if haskey(canteradict,"equation-of-state")
            model = canteradict["equation-of-state"]["model"]
            @error "Currently not supporting $(model) thermo model"
        end
    elseif dict_type == :thermo
        if canteradict["model"] == "NASA7"
            D["type"] = "NASA"
            D["polys"] = [canteradict2rmsdict(canteradict,spcs,names,units,:NASApolynomial;polyindex=i) for i in 1:(length(canteradict["temperature-ranges"])-1)]
        else
             @error "Currently only support NASA7 thermo model from Cantera"
        end
    elseif dict_type == :NASApolynomial
        D["type"] = "NASApolynomial"
        D["Tmax"] = canteradict["temperature-ranges"][polyindex+1]
        D["Tmin"] = canteradict["temperature-ranges"][polyindex]
        D["coefs"] = canteradict["data"][polyindex,:][1]
    elseif dict_type == :reaction #reaction
        equation = canteradict["equation"]
        reversible = true
        forwardable = true
        if occursin(" <=> ",equation) #reversible
            reactants, products = split(canteradict["equation"]," <=> ")
        elseif occursin(" => ",equation) #irreversible
            reactants, products = split(canteradict["equation"]," => ")
            reversible = false
        elseif occursin(" <= ",equation) #not forwardable
            reactants, products = split(canteradict["equation"]," <= ")
            forwardable = false
        end
        reactants = _interpretstoichstring(reactants,names)
        products = _interpretstoichstring(products,names)
        kinetics = canteradict2rmsdict(canteradict,spcs,names,units,:kinetics,numreactants=length(reactants))
        D["reactants"] = reactants
        D["products"] = products
        D["kinetics"] = kinetics
        D["reversible"] = reversible
        D["forwardable"] = forwardable
        D["type"] = "ElementaryReaction"
    elseif dict_type == :kinetics
        if haskey(canteradict,"type")
            kinetics_type = canteradict["type"]
            if kinetics_type == "three-body"
                D["type"] = "ThirdBody"
                D["arr"] = canteradict2rmsdict(canteradict["rate-constant"],spcs,names,units,:arrhenius,numreactants=numreactants+1)
                D["efficiencies"] = get(canteradict,"efficiencies",Dict([]))
            elseif kinetics_type == "falloff"
                if haskey(canteradict,"Troe")
                    D["type"] = "Troe"
                    D["arrhigh"] = canteradict2rmsdict(canteradict["high-P-rate-constant"],spcs,names,units,:arrhenius,numreactants=numreactants)
                    D["arrlow"] = canteradict2rmsdict(canteradict["low-P-rate-constant"],spcs,names,units,:arrhenius,numreactants=numreactants+1)
                    D["efficiencies"] = get(canteradict,"efficiencies",Dict([]))
                    D["a"] = canteradict["Troe"]["A"]
                    D["T1"] = canteradict["Troe"]["T1"]
                    if haskey(canteradict["Troe"],"T2")
                        D["T2"] = canteradict["Troe"]["T2"]
                    else
                        D["T2"] = 0.0
                    end
                    D["T3"] = canteradict["Troe"]["T3"]
                else #Lindemann
                    D["type"] = "Lindemann"
                    D["arrhigh"] = canteradict2rmsdict(canteradict["high-P-rate-constant"],spcs,names,units,:arrhenius,numreactants=numreactants)
                    D["arrlow"] = canteradict2rmsdict(canteradict["low-P-rate-constant"],spcs,names,units,:arrhenius,numreactants=numreactants+1)
                    D["efficiencies"] = get(canteradict,"efficiencies",Dict([]))
                end
            elseif kinetics_type == "pressure-dependent-Arrhenius"
                D["type"] = "PdepArrhenius"
                D["Ps"] = [kinetics["P"] for kinetics in canteradict["rate-constants"]]
                D["arrs"] = [canteradict2rmsdict(kinetics,spcs,names,units,:arrhenius,numreactants=numreactants) for kinetics in canteradict["rate-constants"]]
            elseif kinetics_type == "Chebyshev"
                D["type"] = "Chebyshev"
                D["coefs"] = [canteradict["data"][:,i] for i in 1:size(canteradict["data"],2)]
                D["Tmin"] = tosivalue(canteradict["temperature-range"][1],units=units,value_type=:temperature)
                D["Tmax"] = tosivalue(canteradict["temperature-range"][end],units=units,value_type=:temperature)
                D["Pmin"] = tosivalue(canteradict["pressure-range"][1],units=units,value_type=:pressure)
                D["Pmax"] = tosivalue(canteradict["pressure-range"][end],units=units,value_type=:pressure)
            else
                @error "Currently not supporting $(kinetics_type)"
            end   
        elseif haskey(canteradict,"sticking-coefficient")
            D = canteradict2rmsdict(canteradict["sticking-coefficient"],spcs,names,units,:arrhenius,numreactants=0)
            D["type"] = "StickingCoefficient"
        elseif haskey(canteradict,"rate-constant") #arrhenius
            D = canteradict2rmsdict(canteradict["rate-constant"],spcs,names,units,:arrhenius,numreactants=numreactants)
        else
            @error "Currently not supporting $(canteradict) type kinetics"
        end
        if haskey(canteradict,"coverage-dependencies")
            @error "Currently not supporting coverage dependencies"
        end
    elseif dict_type == :arrhenius
        D["type"] = "Arrhenius"
        D["A"] = tosivalue(canteradict["A"],units=units,value_type=:Afactor,numreactants=numreactants)
        D["Ea"] = tosivalue(canteradict["Ea"],units=units,value_type=:activationenergy)
        D["n"] = canteradict["b"]
    end
    return D
end

function _interpretstoichstring(spcs,names)
    spc_names = Vector{String}()
    
    spcs = split(spcs," ")
    for (ind,item) in enumerate(spcs)
        stoich = tryparse(Int64,item)
        if stoich != nothing
            spc = spcs[ind+1]
            append!(spc_names,[spc for i in 1:(stoich-1)])
        else
            if item in names
                spc = item
                push!(spc_names,spc)
            end
        end
    end
    return spc_names
end

function tosivalue(value;units=nothing,value_type=nothing,numreactants=nothing)
    
    if value isa String #unit specified for this individual value
        value,unit = split(value," ")
        value = parse(Float64,value)
        unit = uparse(unit)
    else 
        if value_type == :temperature
            unit = get(units,"temperature",Unitful.K)
        elseif value_type == :pressure
            unit = get(units,"pressure",Unitful.Pa)
        elseif value_type == :activationenergy
            unit = get(units,"activation-energy",nothing)
            if unit == nothing
                quantity_unit = get(units,"quantity",Unitful.mol)
                energy_unit = get(units,"energy",Unitful.J)
                unit = (energy_unit)/(quantity_unit)
            end
        elseif value_type == :Afactor
            time_unit = get(units,"time",Unitful.s)
            length_unit = get(units,"length",Unitful.m)
            quantity_unit = get(units,"quantity",Unitful.mol)
            if numreactants == 1
                unit = time_unit^(-1)
            elseif numreactants == 2
                unit = (length_unit)^3/(quantity_unit*time_unit)
            elseif numreactants == 3
                unit = (length_unit)^6/((quantity_unit)^2*(time_unit))
            elseif numreactants == 0
                return value
            end
        end
    end
    return upreferred(value*unit).val
end

