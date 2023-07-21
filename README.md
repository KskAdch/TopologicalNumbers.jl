# TopologicalNumbers.jl: A Julia package for calculating topological numbers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/dev/)
[![Build Status](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl)

## Overview

TopologicalNumbers.jl is a Julia package designed to calculate topological numbers, such as the Chern numbers and $\mathbb{Z}_2$ numbers, 
using a numerical approach based on the Fukui-Hatsugai-Suzuki method or the Shiozaki method.  
This package includes the following functions:

- `Dispersion` to calculate the dispersion relation,
- `QuantizedBerryPhase` to calculate the winding numbers in the one-dimensional case,
- `FirstChern` to calculate the first Chern numbers in the two-dimensional case,
- `Z2Invariants2D` to calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case.


This software is released under the MIT License, please see the LICENSE file for more details.  
It is confirmed to work on Julia 1.6 (LTS) and 1.9.


## Installation

To install TopologicalNumbers.jl, run the following command:

```julia
]add https://github.com/KskAdch/TopologicalNumbers.jl
```

Alternatively, you can use:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/KskAdch/TopologicalNumbers.jl")
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
julia> Dispersion(H, 1)
```

![Band structure of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/ec08e558-3b0c-4ab0-9c4f-99b977b20142)

In this case, `1` signifies the dimension of the wavenumber space.

Next, we can calculate the winding numbers using `QuantizedBerryPhase`:

```julia
julia> QuantizedBerryPhase(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the winding number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the winding numbers for each band (mod 2).
`Total` is a quantity that should always return zero.


### Two-dimensional square lattice with flux model

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

To calculate the Dispersion, run:

```julia
julia> Dispersion(H, 2)
```

![Dispersion of 2D square lattice with flux model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/abf47c01-94f3-4c66-8b54-bb243ce48b5f)


Then we can compute the Chern numbers using `FirstChern`:

```julia
julia> FirstChern(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1, -2, -2, 1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.




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

To calculate the Dispersion, execute:

```julia
julia> Dispersion(H, 2)
```

![Dispersion of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/7a7b67a3-7efc-44e6-8e28-bb607e17f7f5)


Next, we can compute the $\mathbb{Z}_2$ numbers using `Z2Invariants2D`:

```julia
julia> Z2Invariants2D(H)
```

The output is:

```julia
(TopologicalNumber = [1 1 1 1; 1 1 1 1], Total = 0)
```

This implies... (edit required)

`Total` is a value that should consistently return zero.


Please see [Documentation](https://kskadch.github.io/TopologicalNumbers.jl/dev/) for more details.
