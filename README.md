# QMCIntegrationInR.jl


add package by typing
```julia
pkg> add https://github.com/Philippe1123/QMCIntegrationInR.jl
```

First cd to
```command
cd .julia/packages/QMCIntegrationInR
```

To run the simulation
```julia
include("IntRQMC_D_MultiLevel_v2.jl")
```



IntRQMC_D_MultiLevel_v4.jl   ->   Mapping of function to the unit cube  with Sobol' sequences and HKKN lattice sequences

IntRQMC_D_MultiLevel_v5.jl   ->   Integration of the function in R with Sobol' sequences and HKKN lattice sequences

IntRQMC_D_MultiLevel_v51.jl   ->   Integration of the function in R with Sobol' sequences and HKKN lattice sequences fixed number of samples, increasing domain

IntRQMC_D_MultiLevel_v52.jl   ->   Integration of the function in R with Sobol' sequences and HKKN lattice sequences , fixed domain, increasing number of samples

IntRQMC_D_MultiLevel_v53.jl   ->   Integration of the function in R with Sobol' sequences and HKKN lattice sequences (read of points from file)

IntRQMC_D_MultiLevel_v6.jl   ->   Improved mapping of function to the unit cube (see paper) with Sobol' sequences and HKKN lattice sequences (with increasing size box)

IntRQMC_D_MultiLevel_v61.jl   ->   Improved mapping of function to the unit cube (see paper) with Sobol' sequences and HKKN lattice sequences (box is full domein R in practice [-20,20]^s ~ R^s)


IntRQMC_D_MultiLevel_v7.jl   ->   Integration of the function in R with Niedereiter Xing sequences and HKKN lattice sequences

** IntRQMC_D_MultiLevel_v8.jl   ->   Integration of the output of a 1D PDE (log diffusion) in R with Sobol' sequences and HKKN lattice sequences **

IntRQMC_D_MultiLevel_v81.jl   ->   Integration of the output of a 1D PDE (log diffusion) in R with Sobol' sequences and HKKN lattice sequences, fixed number of samples, increasing domain

IntRQMC_D_MultiLevel_v82.jl   ->   Integration of the output of a 1D PDE (log diffusion) in R with Sobol' sequences and HKKN lattice sequences, fixed domain, increasing number of samples