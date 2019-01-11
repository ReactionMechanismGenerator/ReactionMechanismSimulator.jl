# RMS - Reaction Mechanism Simulator

## Description
RMS is a Julia package designed for simulating and analyzing large chemical reaction mechanisms.  

## Features
Ideal gas and dilute liquid phases.  
Constant T and P and constant V adiabatic ideal gas reactors.  
Constant T and V dilute liquid reactors.  
Diffusion limited rates. 
Sensitivity analysis for all reactors.  
Flux diagrams with molecular images (if molecular information is provided).  
Handy plotting and other solution analysis tools.  
Easy to add new features.  

## Installation

RMS currently has dependencies that are only available in python 2.  So currently in order to use RMS it is necessary for your PyCall to reference a python 2 installation.  

This can be done with PyCall uninstalled:  
```
using Pkg
ENV["PYTHON"] = "absoulte path to python 2 executable"
Pkg.add("PyCall")
Pkg.build("PyCall")
```
