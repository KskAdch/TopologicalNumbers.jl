# TopologicalNumbers.jl: A Julia package for calculating topological numbers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/dev/)
[![Build Status](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl)

## Overview

TopologicalNumbers.jl is a Julia package designed to calculate topological numbers, such as the Chern numbers and $\mathbb{Z}_2$ numbers, 
using a numerical approach based on the Fukui-Hatsugai-Suzuki method or the Shiozaki method.  
This package mainly includes the following functions:

- `showBand` to calculate the dispersion relation,
- `calcBerryPhase` to calculate the winding numbers in the one-dimensional case,
- `calcChern` to calculate the first Chern numbers in the two-dimensional case,
- `calcZ2` to calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case,
- `calcPhaseDiagram` to calculate the phase diagram using the several methods,
- `calcBerryFlux` to calculate the Berry flux in the two-dimensional case.


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

![Band structure of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/a586aa22-6c79-454e-a82f-6f5056d98f6c)

Next, we can calculate the winding numbers using `calcBerryPhase`:

```julia
julia> calcBerryPhase(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the winding number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the winding numbers for each band (mod 2).
`Total` is a quantity that should always return zero.



One-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (p, 1.0))

julia> param = range(-2.0, 2.0, length=1001)
julia> calcPhaseDiagram(H, param, "BerryPhase"; plot=true)
```

![One-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/2b53e455-83ee-42d5-9824-84120c2be093)

Also, two-dimensional phase diagram is given by:

```julia
julia> param = range(-2.0, 2.0, length=101)
julia> calcPhaseDiagram(H₀, param, param, "BerryPhase"; plot=true)
```

![Two-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/0ceef1a3-01fd-4e8b-9f01-4a4932039d26)



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


![Band structure of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/b05dd8a2-8047-43d6-a835-d4c18348ec97)


Then we can compute the Chern numbers using `calcChern`:

```julia
julia> calcChern(H)
```

The output is:

```julia
(TopologicalNumber = [1, -1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.



One-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (p, 2.5))

julia> param = range(-π, π, length=1000)
julia> calcPhaseDiagram(H, param, "Chern"; plot=true)
```

![One-dimensional phase diagram of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/8af528e5-20df-4b9e-ad1b-cc11de0b902a)

Also, two-dimensional phase diagram is given by:

```julia
julia> param1 = range(-π, π, length=100)
julia> param2 = range(-6.0, 6.0, length=100)
julia> calcPhaseDiagram(H₀, param1, param2, "Chern"; plot=true)
```

![Two-dimensional phase diagram of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/bafd50e9-2752-4f81-9b9d-e271a8d8061f)


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

![Dispersion of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/a9cf9768-6920-45e6-89bd-ed7ec434152c)


Next, we can compute the $\mathbb{Z}_2$ numbers using `calcZ2`:

```julia
julia> calcZ2(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the $\mathbb{Z}_2$ number for each pair of two energy bands. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the $\mathbb{Z}_2$ numbers for each pair of two energy bands.
`Total` is a quantity that should always return zero.


One-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (p, 0.25))

julia> param = range(-2, 2, length=1000)
julia> calcPhaseDiagram(H, param, "Z2"; plot=true)
```

![One-dimensional phase diagram of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/8e27a9d9-f52a-4f24-9d9e-c1254edabdcc)


Also, two-dimensional phase diagram is given by:

```julia
julia> param1 = range(-2, 2, length=100)
julia> param2 = range(-0.5, 0.5, length=100)
julia> calcPhaseDiagram(H₀, param1, param2, "Z2"; plot=true)
```


![Two-dimensional phase diagram of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/f8a36504-372b-4e23-b7e9-02ada709bdc4)



Please see [Documentation](https://kskadch.github.io/TopologicalNumbers.jl/dev/) for more details.

