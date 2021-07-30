# <img align="top" src="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl/blob/master/logos/rms-logo-small.png"> RMS - Reaction Mechanism Simulator

[![Build status](https://img.shields.io/travis/ReactionMechanismGenerator/ReactionMechanismSimulator.jl/master.svg)](https://travis-ci.org/ReactionMechanismGenerator/ReactionMechanismSimulator.jl)
[![codecov](https://codecov.io/gh/ReactionMechanismGenerator/ReactionMechanismSimulator.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ReactionMechanismGenerator/ReactionMechanismSimulator.jl)

## Description
RMS is a Julia package designed for simulating and analyzing large chemical reaction mechanisms. 

RMS has been used in many applications: 
* Combustion:
  * Ignition quality tester
  * Rapid compression machine
  * Shock tube
  * Flow tube
* Pharmaceutical degradation
* Polymer film growth
* Gas phase catalysis
* Electrocatalytic reduction of Nitrogen to ammonia
* Solid electrolyte interfaces in batteries
* Liquid oxidation
* Pyrolysis of heavy oils

## Features
Ideal gas, dilute liquid and ideal surface phases. 
Wide selection of domains including but not limited to constant temperature and pressure, constant volume, parameterized volume, constant temperature and volume and constant temperature, potential and area. All of these have analytic jacobians! Easy to add more!
Domains can be coupled to fixed interfaces such as inlets and oulets and also to dynamic interfaces such as a surface-gas reactive interface between surface and gas phase domains.
Diffusion limited rates.
Forward and adjoint sensitivity analysis for all reactors.  
Flux diagrams with molecular images (if molecular information is provided).  
Handy plotting and other solution analysis tools.  
Easy to add new phases, domains, interfaces and other new features. 

## Installation

RMS can be installed with:
```
using Pkg
Pkg.add("ReactionMechanismSimulator")
Pkg.build("ReactionMechanismSimulator")

using ReactionMechanismSimulator
```

Detailed instructions and documentation are currently available in the <a href="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl/wiki">wiki</a>.
