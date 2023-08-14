# TopologicalNumbers.jl: A Julia package for calculating topological numbers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/dev/)
[![Build Status](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl)

## Overview

TopologicalNumbers.jl is a Julia package designed to calculate topological numbers, such as the Chern numbers and $\mathbb{Z}_2$ numbers, 
using a numerical approach based on the Fukui-Hatsugai-Suzuki method or the Shiozaki method.  
This package includes the following functions:

- `showBand` to calculate the dispersion relation,
- `calcBerryPhase` to calculate the winding numbers in the one-dimensional case,
- `calcChern` to calculate the first Chern numbers in the two-dimensional case,
- `calcZ2` to calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case.


This software is released under the MIT License, please see the LICENSE file for more details.  
It is confirmed to work on Julia 1.6 (LTS) and 1.9.


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
julia> function H(k) # set SSH Hamiltonian function of wavenumber k
    g = 0.9
    
    [
        0 g+exp(-im*k)
        g+exp(im*k) 0
    ]
end
```

The band structure is computed as follows:

```julia
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
julia> function H0(k, p)
            [
                0 p[1]+p[2]*exp(-im * k)
                p[1]+p[2]*exp(im * k) 0
            ]
        end
julia> H(k, p) = H0(k, (p, 1.0))

julia> param = range(-2.0, 2.0, length=1001)
julia> calcPhaseDiagram(H, param, "BerryPhase"; plot=true)
```

![One-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/2b53e455-83ee-42d5-9824-84120c2be093)

Also, two-dimensional phase diagram is given by:

```julia
julia> param = range(-2.0, 2.0, length=101)
julia> calcPhaseDiagram(H0, param, param, "BerryPhase"; plot=true)
```

![Two-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/0ceef1a3-01fd-4e8b-9f01-4a4932039d26)


### Two-dimensional square lattice model with flux

A two-dimensional example is presented here:

```julia
julia> function H(k) # landau
    k1, k2 = k
    t = 1

    Hsize = 6
    Hmat = zeros(ComplexF64, Hsize, Hsize)

    for i in 1:Hsize
        Hmat[i, i] = -2*cos(k2-2pi*i/Hsize)
    end

    for i in 1:Hsize-1
        Hmat[i, i+1] = -t
        Hmat[i+1, i] = -t
    end

    Hmat[1, Hsize] = -t*exp(-im*k1)
    Hmat[Hsize, 1] = -t*exp(im*k1)
    
    Hmat
end
```

To calculate the dispersion, run:

```julia
julia> showBand(H; value=false, disp=true)
```

![Dispersion of 2D square lattice with flux model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/8470121f-bdad-4960-9848-7ade1ae805d3)


Then we can compute the Chern numbers using `calcChern`:

```julia
julia> calcChern(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1, -2, -2, 1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.


One-dimensional phase diagram is given by:

```julia
julia> function H(k, p)
    k1, k2 = k
    t = p

    Hsize = 6
    Hmat = zeros(ComplexF64, Hsize, Hsize)

    for i in 1:Hsize
        Hmat[i, i] = -2 * cos(k2 - 2pi * i / Hsize)
    end

    for i in 1:Hsize-1
        Hmat[i, i+1] = -t
        Hmat[i+1, i] = -t
    end

    Hmat[1, Hsize] = -t * exp(-im * k1)
    Hmat[Hsize, 1] = -t * exp(im * k1)

    Hmat
end

julia> param = range(-2.0, 2.0, length=500)
julia> calcPhaseDiagram(H, param, "Chern"; plot=true)
```

![One-dimensional phase diagram](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/f9c179c3-1275-4640-ac21-0d10737fcaf7)




### The Bernevig-Hughes-Zhang (BHZ) model

As an example of a two-dimensional topological insulator, the BHZ model is presented here:

```julia
julia> function H(k) # BHZ
    k1, k2 = k
    t = 1

    R0 = -2(cos(k1) + cos(k2)) + 1
    R3 = 2sin(k2)
    R4 = 2sin(k1)
    R5 = -2t*(cos(k1) + cos(k2)) + 1

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

    R0*a0+R3*a3+R4*a4+R5*a5
end
```

To calculate the dispersion, execute:

```julia
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

The first argument `TopologicalNumber` in the named tuple is an vector that stores the $\mathbb{Z}_2$ number for each each pair of two energy bands. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the $\mathbb{Z}_2$ numbers for each pair of two energy bands.
`Total` is a quantity that should always return zero.


Please see [Documentation](https://kskadch.github.io/TopologicalNumbers.jl/dev/) for more details.
