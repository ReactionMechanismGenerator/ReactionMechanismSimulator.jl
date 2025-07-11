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

## Julia-Python Linking
The above instructions will automatically handle Julia-Python linking. However, in some cases it can be useful to use python from a specific conda environment. For these cases we provide instructions to relink Julia to a different Anaconda Python environment.

First, makesure `pyjuliacall` is installed in your Conda environment:
```
conda install -y conda-forge:pyjuliacall
```
Then, activate your Conda environment:
```
conda activate your_conda_environment
```

```
conda env config vars set JULIA_CONDAPKG_BACKEND=Null
conda env config vars set JULIA_CONDAPKG_EXE=$(which conda)
conda env config vars set JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python
conda env config vars set PYTHON_JULIAPKG_EXE=$(which julia)
conda env config vars set PYTHON_JULIAPKG_PROJECT=$CONDA_PREFIX/julia_env

export JULIA_CONDAPKG_BACKEND=Null
export JULIA_CONDAPKG_EXE="$(which conda)"
export JULIA_PYTHONCALL_EXE="$CONDA_PREFIX/bin/python"
export PYTHON_JULIAPKG_EXE="$(which julia)"
export PYTHON_JULIAPKG_PROJECT="$CONDA_PREFIX/julia_env"
```
Then you install RMS into the Julia project location within your Conda environment:
```
julia
using Pkg
Pkg.activate(ENV["PYTHON_JULIAPKG_PROJECT"])
Pkg.add("ReactionMechanismSimulator") # Standard Installation
Pkg.add(Pkg.PackageSpec(path=RMS_PATH)) # Developer Installation
Pkg.instantiate()
Pkg.build("ReactionMechanismSimulator")
```
To test whether Julia-Python linking is set up correctly:
```
python
from juliacall import Main
import sys
Main.seval('Base.identify_package("ReactionMechanismSimulator")')
```
You should see an output in the following format:
```
Julia: ReactionMechanismSimulator [32-digit uuid]
```
If the output is empty, the linking is not set up properly.


## Testing RMS
Unit and functional tests for RMS can be run with:
```
import Pkg
Pkg.test("ReactionMechanismSimulator")
```

## pyrms
We also provide a python wrapper for RMS, [pyrms](https://github.com/ReactionMechanismGenerator/pyrms). Installation instructions are available on its github page.
