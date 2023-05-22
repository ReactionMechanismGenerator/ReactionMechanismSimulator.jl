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
* Ideal gas, dilute liquid and ideal surface phases. 
* Wide selection of domains including but not limited to constant temperature and pressure, constant volume, parameterized volume, constant temperature and volume and constant temperature, potential and area. All of these have analytic jacobians! Easy to add more!
* Domains can be coupled to fixed interfaces such as inlets and outlets and also to dynamic interfaces such as surface-gas reactive interfaces between surface and gas phase domains.
* Diffusion limited rates.
* Forward and adjoint sensitivity analysis for all reactors.  
* Flux diagrams with molecular images (if molecular information is provided).  
* Handy plotting and other solution analysis tools.  
* Easy to add new phases, domains, interfaces and other new features. 

## How to cite
Please include the following citations for ReactionMechanismSimulator.jl in general and for transitory sensitivities and the automatic mechanism analysis toolkit respectively. 

- <div class="csl-entry">Johnson, M. S., Pang, H.-W., Payne, A. M., &#38; Green, W. H. (2023). <i>ReactionMechanismSimulator.jl: A Modern Approach to Chemical Kinetic Mechanism Simulation and Analysis</i>. https://doi.org/10.26434/CHEMRXIV-2023-TJ34T</div>
- <div class="csl-entry">Johnson, M. S., McGill, C. J., &#38; Green, W. H. (2022). <i>Transitory Sensitivity in Automatic Chemical Kinetic Mechanism Analysis</i>. https://doi.org/10.26434/CHEMRXIV-2022-ZSFJC</div>

## Installation

RMS can be installed with:
```
using Pkg
Pkg.add("ReactionMechanismSimulator")
Pkg.build("ReactionMechanismSimulator")

using ReactionMechanismSimulator
```

Detailed instructions and documentation are currently available in the <a href="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl/wiki">wiki</a>.
