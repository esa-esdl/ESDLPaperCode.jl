# Case studies from the ESDL Paper (Mahecha et. al)

This repository contains the scripts to reproduce all results and figures from the case studies
in the ESDL overview paper (citation to follow).

## Installation

To install run this code you need to install the unregistered dependencies first:

````julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/esa-esdl/ESDL.jl",rev="v0.7.5"))
Pkg.add(PackageSpec(url="https://github.com/esa-esdl/ESDLPlots.jl",rev="v0.2.3"))
````

After that you can run

````julia
Pkg.add(PackageSpec(url="https://github.com/esa-esdl/PaperCode.jl")
````
