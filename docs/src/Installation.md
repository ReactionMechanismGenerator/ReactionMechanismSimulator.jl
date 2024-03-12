# Install Julia
RMS is written in Juila language. So before stepping any further, Julia needs to be installed. The download links can be found at [download Julia](https://julialang.org/downloads/). More instructions can be found from the [instruction page](https://julialang.org/downloads/platform/).

## Standard Installation
With julia RMS can be installed with:

```
using Pkg
Pkg.add("ReactionMechanismSimulator")
Pkg.build("ReactionMechanismSimulator")
```

## Developer Installation
Clone RMS to your machine in an appropriate location we will refer to as `RMS_PATH``:  
```
git clone https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl.git
```
then you can install with 
```
import Pkg
Pkg.develop(Pkg.PackageSpec(path=RMS_PATH)) 
Pkg.build("ReactionMechanismSimulator")
```

## Testing RMS
Unit and functional tests for RMS can be run with: 
```
import Pkg
Pkg.test("ReactionMechanismSimulator")
```

## pyrms
We also provide a python wrapper for RMS, [pyrms](https://github.com/ReactionMechanismGenerator/pyrms). Installation instructions are available on its github page. 

## Julia-Python Linking
The above instructions will automatically handle Julia-Python linking. However, in some cases it can be useful to use python from a specific conda environment. For these cases we provide instructions to relink Julia to a different Anaconda Python environment where `PATH_TO_YOUR_ENV` is the path to the Anaconda environment and `PATH_TO_PYTHON` is the path to the associated Python executable: 

```
import Pkg
ENV["CONDA_JL_HOME"] = PATH_TO_YOUR_ENV
ENV["PYTHON"] = PATH_TO_PYTHON
Pkg.add("PyCall")
Pkg.build("PyCall")                            
```
