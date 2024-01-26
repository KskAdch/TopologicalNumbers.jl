# TopologicalNumbers.jl: A Julia package for topological number computation

![TopologicalNumbers logo](docs/src/assets/logo_title.svg#gh-light-mode-only)
![TopologicalNumbers logo](docs/src/assets/logo_title_dark.svg#gh-dark-mode-only)

---


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/dev/)
[![Build Status](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl)
[![DOI](https://zenodo.org/badge/666049458.svg)](https://zenodo.org/badge/latestdoi/666049458)

## Overview

TopologicalNumbers.jl is a Julia package designed to compute topological numbers, 
such as the first and second Chern numbers and $\mathbb{Z}_2$ numbers, 
using a numerical approach based on the Fukui-Hatsugai-Suzuki method or the Shiozaki method, 
or method of calculating the Weyl nodes.  
This package includes the following functions:

- Computation of the dispersion relation.
- Provides numerical calculation methods for various types of topological numbers.
- Computation of the phase diagram.
- Compute Pfaffian and tridiagonarize skew-symmetric matrix (migration to Julia from [PFAPACK](https://pypi.org/project/pfapack/)).
- Utility functions for plotting.
- Support parallel computing using `MPI`.


The correspondence between the spatial dimension of the system and the supported topological numbers is as follows.

|Dimension|Function                                                                                            |
|---------|----------------------------------------------------------------------------------------------------|
|1D       |Calculation of Berry Phases ($\mathbb{Z}$)<br>                                                                   |
|2D       |Calculation of local Berry Fluxes ($\mathbb{Z}$)<br>Calculation of first Chern numbers ($\mathbb{Z}$)<br>Calculation of $\mathbb{Z}_2$ numbers ($\mathbb{Z}_2$)        |
|3D       |Calculation of Weyl nodes ($\mathbb{Z}$)<br> Calculation of first Chern numbers in sliced Surface ($\mathbb{Z}$)<br>Finding Weyl points ($\mathbb{Z}$)<br>|
|4D       |Calculation of second Chern numbers ($\mathbb{Z}$)<br>                                                                  |



This software is released under the MIT License, please see the LICENSE file for more details.  
It is confirmed to work on Julia 1.6 (LTS) and 1.10.


## Installation

To install TopologicalNumbers.jl, run the following command:

```julia
pkg> add TopologicalNumbers
```

Alternatively, you can use:

```julia
julia> using Pkg
julia> Pkg.add("TopologicalNumbers")
```



## Examples

### The Su-Schriffer-Heeger (SSH) model

Here's a simple example of the SSH Hamiltonian:

```julia
julia> using TopologicalNumbers
julia> function H₀(k, p)
            [
                0 p[1]+p[2]*exp(-im * k)
                p[1]+p[2]*exp(im * k) 0
            ]
        end
```

The band structure is computed as follows:

```julia
julia> H(k) = H₀(k, (0.9, 1.0))
julia> showBand(H; value=false, disp=true)
```

![Band structure of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/1a499692-681c-4166-b99c-0fd5cdcd506f)

Next, we can calculate the winding numbers using `BPProblem`:

```julia
julia> prob = BPProblem(H);
julia> sol = solve(prob)
```

The output is:

```julia
BPSolution{Vector{Int64}, Int64}([1, 1], 0)
```

The first argument `TopologicalNumber` in the named tuple is a vector that stores the winding number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the winding numbers for each band (mod 2).
`Total` is a quantity that should always return zero.

You can access these values as follows:

```julia
julia> sol.TopologicalNumber
2-element Vector{Int64}:
 1
 1

julia> sol.Total
0
```


A one-dimensional phase diagram is given by:

```julia
julia> param = range(-2.0, 2.0, length=1001)

julia> prob = BPProblem(H₀);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = -2.0:0.004:2.0, nums = [1 1; 1 1; … ; 1 1; 1 1])
```

![One-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/7a6bec77-9140-4257-ba66-8280eef4fe1d)


### Haldane model

Hamiltonian of Haldane model is given by:

```julia
julia> function H₀(k, p) # landau
           k1, k2 = k
           J = 1.0
           K = 1.0
           ϕ, M = p

           h0 = 2K * cos(ϕ) * (cos(k1) + cos(k2) + cos(k1 + k2))
           hx = J * (1 + cos(k1) + cos(k2))
           hy = J * (-sin(k1) + sin(k2))
           hz = M - 2K * sin(ϕ) * (sin(k1) + sin(k2) - sin(k1 + k2))

           s0 = [1 0; 0 1]
           sx = [0 1; 1 0]
           sy = [0 -im; im 0]
           sz = [1 0; 0 -1]

           h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
       end
```

The band structure is computed as follows:

```julia
julia> H(k) = H₀(k, (π/3, 0.5))
julia> showBand(H; value=false, disp=true)
```



![Band structure of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/5f1da3f6-c0d8-4bd7-84cf-a591bf59a137)


Then we can compute the Chern numbers using `FCProblem`:

```julia
julia> prob = FCProblem(H);
julia> sol = solve(prob)
```

The output is:

```julia
FCSolution{Vector{Int64}, Int64}([1, -1], 0)
```

The first argument `TopologicalNumber` in the named tuple is a vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.



A one-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (1, p, 2.5));
julia> param = range(-π, π, length=1000);

julia> prob = FCProblem(H);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = -3.141592653589793:0.006289474781961547:3.141592653589793, nums = [0 0; 0 0; … ; 0 0; 0 0])
```

![One-dimensional phase diagram of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/bb9682bc-55a2-41c8-ba1f-b90763e9233f)


Also, two-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (1, p[1], p[2]));
julia> param1 = range(-π, π, length=100);
julia> param2 = range(-6.0, 6.0, length=100);

julia> prob = FCProblem(H);
julia> sol = calcPhaseDiagram(prob, param1, param2; plot=true)
(param1 = -3.141592653589793:0.06346651825433926:3.141592653589793, param2 = -6.0:0.12121212121212122:6.0, nums = [0 0 … 0 0; 0 0 … 0 0;;; 0 0 … 0 0; 0 0 … 0 0;;; 0 0 … 0 0; 0 0 … 0 0;;; … ;;; 0 0 … 0 0; 0 0 … 0 0;;; 0 0 … 0 0; 0 0 … 0 0;;; 0 0 … 0 0; 0 0 … 0 0])
```

![Two-dimensional phase diagram of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/a61ca3b0-28b9-44a8-b50f-67547a453ebe)


### The Bernevig-Hughes-Zhang (BHZ) model

As an example of a two-dimensional topological insulator, the BHZ model is presented here:

```julia
julia> function H₀(k, p) # BHZ
    k1, k2 = k
    tₛₚ = 1
    t₁ = ϵ₁ = 2
    ϵ₂, t₂ = p

    R0 = -t₁*(cos(k1) + cos(k2)) + ϵ₁/2
    R3 = 2tₛₚ*sin(k2)
    R4 = 2tₛₚ*sin(k1)
    R5 = -t₂*(cos(k1) + cos(k2)) + ϵ₂/2

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    a0 = kron(s0, s0)
    a1 = kron(sx, sx)
    a2 = kron(sx, sy)
    a3 = kron(sx, sz)
    a4 = kron(sy, s0)
    a5 = kron(sz, s0)

    R0 .* a0 .+ R3 .* a3 .+ R4 .* a4 .+ R5 .* a5
end
```

To calculate the dispersion, execute:

```julia
julia> H(k) = H₀(k, (2, 2))
julia> showBand(H; value=false, disp=true)
```

![Dispersion of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/de14907c-777f-4667-810b-54c10888dfa1)



Next, we can compute the $\mathbb{Z}_2$ numbers using `Z2Problem`:

```julia
julia> prob = Z2Problem(H);
julia> sol = solve(prob)
```


The output is:

```julia
Z2Solution{Vector{Int64}, Nothing, Int64}([1, 1], nothing, 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the $\mathbb{Z}_2$ number for Energy bands below and above some filling condition that you selected in the options (the default is the half-filling). 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the $\mathbb{Z}_2$ numbers for each pair of two energy bands.
`Total` is a quantity that should always return zero.



A one-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (p, 0.25));
julia> param = range(-2, 2, length=1000);

julia> prob = Z2Problem(H);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = -2.0:0.004004004004004004:2.0, nums = [0 0; 0 0; … ; 0 0; 0 0])
```

![One-dimensional phase diagram of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/5d6d5364-68d0-4423-8ecf-49bf0538af63)


Also, two-dimensional phase diagram is given by:

```julia
julia> param1 = range(-2, 2, length=100);
julia> param2 = range(-0.5, 0.5, length=100);

julia> prob = Z2Problem(H₀);
julia> calcPhaseDiagram(prob, param1, param2; plot=true)
```


![Two-dimensional phase diagram of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/802eedbe-c893-44b4-8267-d80e1745415a)



### The four-dimensional lattice Dirac model

As an example of a four-dimensional topological insulator, the lattice Dirac model is presented here:

```julia
julia> function H₀(k, p) # lattice Dirac
            k1, k2, k3, k4 = k
            m = p

            # Define Pauli matrices and Gamma matrices
            σ₀ = [1 0; 0 1]
            σ₁ = [0 1; 1 0]
            σ₂ = [0 -im; im 0]
            σ₃ = [1 0; 0 -1]
            g1 = kron(σ₁, σ₀)
            g2 = kron(σ₂, σ₀)
            g3 = kron(σ₃, σ₁)
            g4 = kron(σ₃, σ₂)
            g5 = kron(σ₃, σ₃)

            h1 = m + cos(k1) + cos(k2) + cos(k3) + cos(k4)
            h2 = sin(k1)
            h3 = sin(k2)
            h4 = sin(k3)
            h5 = sin(k4)

            # Return the Hamiltonian matrix
            h1 .* g1 .+ h2 .* g2 .+ h3 .* g3 .+ h4 .* g4 .+ h5 .* g5
        end
```
You can also use our preset Hamiltonian function `LatticeDirac` to define the same Hamiltonian matrix as follows:

```julia
julia> H₀(k, p) = LatticeDirac(k, p)
```

Then we can compute the second Chern numbers using `SCProblem`:

```julia
julia> H(k) = H₀(k, -3.0)

julia> prob = SCProblem(H);
julia> sol = solve(prob)
```

The output is:

```julia
SCSolution{Float64}(0.9793607631927376)
```

The argument `TopologicalNumber` in the named tuple stores the second Chern number with some filling condition that you selected in the options (the default is the half-filling). 


A phase diagram is given by:

```julia
julia> param = range(-4.9, 4.9, length=50);
julia> prob = SCProblem(H₀);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
```

![Dense phase diagram of lattice Dirac model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/9449651c-46f4-4141-ac87-161e9e5fbf28)


If you want to use a parallel environment, you can utilize `MPI.jl`.
Let's create a file named `test.jl` with the following content:
```julia
using TopologicalNumbers
using MPI

H₀(k, p) = LatticeDirac(k, p)
H(k) = H₀(k, -3.0)

param = range(-4.9, 4.9, length=10)

prob = SCProblem(H₀)
result = calcPhaseDiagram(prob, param; parallel=UseMPI(MPI), progress=true)

plot1D(result; labels=true, disp=false, pdf=true)
```
You can perform calculations using `mpirun` (for example, with `4` cores) as follows:
```bash
mpirun -np 4 julia --project test.jl
```





Please see [Documentation](https://kskadch.github.io/TopologicalNumbers.jl/dev/) for more details.

