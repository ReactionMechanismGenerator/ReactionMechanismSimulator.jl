import yaml
import os
from rmgpy.chemkin import loadChemkinFile
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.thermo.wilhoit import *
from rmgpy.kinetics.arrhenius import *
from rmgpy.kinetics.falloff import *
from rmgpy.kinetics.chebyshev import *
from rmgpy.data.solvation import *

def convertchemkin2yml(chemkinpath,spcdictpath=None,output="chem.rms"):
    if spcdictpath:
        spcs,rxns = loadChemkinFile(chemkinpath,dictionaryPath=spcdictpath)
    else:
        spcs,rxns = loadChemkinFile(chemkinpath)
    writeyml(spcs,rxns,path=output)

def writeyml(spcs,rxns,path="chem.yml"):
    D = getmechdict(spcs,rxns)
    yaml.dump(D,stream=file(path,'w'))

def getmechdict(spcs,rxns):
    D = dict()
    D["Units"] = dict()
    D["Reactions"] = []
    D["Phases"] = [dict()]
    D["Phases"][0]["name"] = "phase"
    D["Phases"][0]["Species"] = [obj2dict(x,spcs) for x in spcs]
    D["Reactions"] = [obj2dict(x,spcs) for x in rxns]
    return D

def getradicals(spc):
    if spc.molecule[0].toSMILES() == "[O][O]":
        return 0
    else:
        return spc.molecule[0].multiplicity-1

def obj2dict(obj,spcs,label="solvent"):
    D = dict()
    if isinstance(obj,Species):
        D["name"] = obj.label
        D["type"] = "Species"
        D["smiles"] = obj.molecule[0].toSMILES()
        D["thermo"] = obj2dict(obj.thermo,spcs)
        if D["smiles"]  == "[O][O]":
            D["radicalelectrons"] = 0
        else:
            D["radicalelectrons"] = obj.molecule[0].multiplicity-1
    elif isinstance(obj,NASA):
        D["polys"] = [obj2dict(k,spcs) for k in obj.polynomials]
        D["type"] = "NASA"
    elif isinstance(obj,NASAPolynomial):
        D["type"] = "NASApolynomial"
        D["coefs"] = obj.coeffs.tolist()
        D["Tmax"] = obj.Tmax.value_si
        D["Tmin"] = obj.Tmin.value_si
    elif isinstance(obj,Reaction):
        D["reactants"] = [x.label for x in obj.reactants]
        D["products"] = [x.label for x in obj.products]
        D["kinetics"] = obj2dict(obj.kinetics,spcs)
        D["type"] = "ElementaryReaction"
        D["radicalchange"] = sum([getradicals(x) for x in obj.products]) -  sum([getradicals(x) for x in obj.reactants])
    elif isinstance(obj,Arrhenius):
        D["type"] = "Arrhenius"
        D["A"] = obj.A.value_si
        D["Ea"] = obj.Ea.value_si
        D["n"] = obj.n.value_si
    elif isinstance(obj,PDepArrhenius):
        D["type"] = "PdepArrhenius"
        D["Ps"] = obj.pressures.value_si.tolist()
        D["arrs"] = [obj2dict(x,spcs) for x in obj.arrhenius]
    elif isinstance(obj,MultiArrhenius):
        D["type"] = "MultiArrhenius"
        D["arrs"] = [obj2dict(x,spcs) for x in obj.arrhenius]
    elif isinstance(obj,MultiPDepArrhenius):
        D["type"] = "MultiPdepArrhenius"
        D["parrs"] = [obj2dict(x,spcs) for x in obj.arrhenius]
    elif isinstance(obj,ThirdBody):
        D["type"] = "ThirdBody"
        D["arr"] = obj2dict(obj.arrheniusLow,spcs)
        D["efficiencies"] = {spcs[i].label:float(val) for i,val in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1}
    elif isinstance(obj,Lindemann):
        D["type"] = "Lindemann"
        D["arrhigh"] = obj2dict(obj.arrheniusHigh,spcs)
        D["arrlow"] = obj2dict(obj.arrheniusLow,spcs)
        D["efficiencies"] = {spcs[i].label:float(val) for i,val in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1}
    elif isinstance(obj,Troe):
        D["type"] = "Troe"
        D["arrhigh"] = obj2dict(obj.arrheniusHigh,spcs)
        D["arrlow"] = obj2dict(obj.arrheniusLow,spcs)
        D["efficiencies"] = {spcs[i].label:float(val) for i,val in enumerate(obj.getEffectiveColliderEfficiencies(spcs)) if val != 1}
        D["a"] = obj.alpha
        D["T1"] = obj.T1.value_si
        if obj.T2:
            D["T2"] = obj.T2.value_si
        else:
            D["T2"] = 0.0
        D["T3"] = obj.T3.value_si
    elif isinstance(obj,Chebyshev):
        D["type"] = "Chebyshev"
        D["coefs"] = obj.coeffs.value_si.tolist()
        D["Tmin"] = obj.Tmin.value_si
        D["Tmax"] = obj.Tmax.value_si
        D["Pmin"] = obj.Pmin.value_si
        D["Pmax"] = obj.Pmax.value_si
    elif isinstance(obj,Wilhoit):
        D["type"] = "Wilhoit"
        D["coefs"] = [obj.a0,obj.a1,obj.a2,obj.a3]
        D["Cp0"] = obj.Cp0.value_si
        D["Cpinf"] = obj.CpInf.value_si
        D["H0"] = obj.H0.value_si
        D["S0"] = obj.S0.value_si
        D["B"] = obj.B.value_si
    elif isinstance(obj,SolventData):
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
    elif obj is None:
        return None
    else:
        raise ValueError
    return D
