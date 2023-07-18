```@meta
CurrentModule = TopologicalNumbers
```

# TopologicalNumbers.jl

Documentation for [TopologicalNumbers](https://github.com/KskAdch/TopologicalNumbers.jl).

```@index
```

```@autodocs
Modules = [TopologicalNumbers]
```



# TopologicalNumbers.jl: A Julia package for calculating topological numbers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/dev/)
[![Build Status](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl)

## Abstract

TopologicalNumbers.jl is a Julia package to calculate topological numbers, including the Chern numbers and $\mathbb{Z}_2$ numbers using a numerical approach based on the Fukui-Hatsugai-Suzuki method or the Shiozaki method.
This package serves the following functions:

- `Dispersion` to calculate the dispersion relation,
- `QuantizedBerryPhase` to calculate the winding numbers in the one-dimensional case,
- `FirstChern` to calculate the first Chern numbers in the two-dimensional case,
- `Z2invariants2D` to calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case.


This software is released under the MIT License, see LICENSE.
We checked that it works on Julia 1.6 (LTS) and 1.9.


## Install

You can install TopologicalNumbers.jl with:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/KskAdch/TopologicalNumbers.jl")
```



## Examples

### The Su-Schriffer-Heeger (SSH) model

Here is a simple example of the SSH Hamiltonian:

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

A band structure is here:
```
julia> Dispersion(H, 1)
```

![Band structure of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/ec08e558-3b0c-4ab0-9c4f-99b977b20142)

Here, `1` is the dimension of the wavenumber space.

Then we can calculate the winding numbers using `QuantizedBerryPhase`:

```julia
julia> QuantizedBerryPhase(H)
```

Output is 

```julia
(TopologicalNumber = [1, 1], Total = 0)
```

This means...(edit required)


### Two-dimensional square lattice with flux model

Two-dimensional example is here:

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

To calcurate Dispersion, run:

```julia
julia> Dispersion(H, 2)
```

![Dispersion of 2D square lattice with flux model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/abf47c01-94f3-4c66-8b54-bb243ce48b5f)


Then we can calculate the Chern numbers using `FirstChern`:

```julia
julia> FirstChern(H)
```

Output is 

```julia
(TopologicalNumber = [1, 1, -2, -2, 1, 1], Total = 0)
```

This means...(edit required)
