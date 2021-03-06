[![DOI](https://zenodo.org/badge/206770757.svg)](https://zenodo.org/badge/latestdoi/206770757)


# Earth system data cubes unravel global multivariate dynamics" 

Example code used to produce the results and figures in Mahecha, Gans, 
et al. (2020) "Earth system data cubes unravel global multivariate 
dynamics" Earth System Dynamics.

## Installation

To install run this code you need to install the unregistered dependencies first:

````julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/esa-esdl/ESDL.jl",rev="ESDL_Paper"))
````

After that you can run

````julia
Pkg.develop(PackageSpec(url="https://github.com/esa-esdl/PaperCode.jl"))
````

which will create a clone of the respoitory structure in your julia development folder.
To re-run all studies and reproduce the plots, you can run

````julia
using ESDLPaperCode
run_all()
````

or

````julia
run_case_study(i)
````

to run a single case study. Note that this will first download the necessary variables from the data cube (approx 20GB) which might take some time.
You should also have some free space on your hard drive since the script will store a large amount of temporary data, so an addition 100GB of free space should be available.

If you prefer to run the code interactively as jupyter notebooks, make sure
to activate the ESDLPaperCode environment to have all necessary packages available:

````julia
using Pkg
Pkg.activate("ESDLPaperCode")
````
