```@meta
CurrentModule = TopologicalNumbers
```

# TopologicalNumbers.jl

Documentation for [TopologicalNumbers](https://github.com/KskAdch/TopologicalNumbers.jl).



[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/dev/)
[![Build Status](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl)

## Overview

TopologicalNumbers.jl is a Julia package designed to calculate topological numbers, such as the Chern numbers and $\mathbb{Z}_2$ numbers, 
using a numerical approach based on the Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite), the Shiozaki method [Shiozaki2023discrete](@cite) or method of calculating the Weyl node [Hirayama2017;Yang2011](@cite).

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

!!! note "Note"
    If you are using a headless server, you may receive an error related to the `GLMakie` package. To resolve this issue, please refer to the [`Makie` documentation](https://docs.makie.org/stable/documentation/headless/) or the [`GLMakie` troubleshooting guide](https://github.com/MakieOrg/Makie.jl/tree/master/GLMakie#troubleshooting-opengl).



## Package Features

This package includes the following functions:

- `showBand` to calculate the dispersion relation,
- `calcBerryPhase` to calculate the winding numbers in the one-dimensional case,
- `calcChern` to calculate the first Chern numbers in the two-dimensional case,
- `calcZ2` to calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case,
- `calcPhaseDiagram` to calculate the phase diagram using the several methods,
- `calcBerryFlux` to calculate the Berry flux in the two-dimensional case,
- `calcWeylNode` to calculate the Weyl node in the three-dimensional case,
- `calcChernSurface` to calculate the Chern numbers in the three-dimensional case,
- `findWeylPoint` to find the Weyl points in the three-dimensional case.