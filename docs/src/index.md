```@meta
CurrentModule = TopologicalNumbers
```

# TopologicalNumbers.jl


```@raw html

<picture>
  <!-- ダークモード用 -->
  <source
    srcset="assets/logo_title_dark.svg"
    media="(prefers-color-scheme: dark)"
  />
  <!-- ライトモード用 -->
  <source srcset="assets/logo_title.svg" />
  <img
    src="assets/logo_title.svg"
    alt="TopologicalNumbers logo"
  />
</picture>

```


---


Documentation for [TopologicalNumbers](https://github.com/KskAdch/TopologicalNumbers.jl).



[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KskAdch.github.io/TopologicalNumbers.jl/dev/)
[![Build Status](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KskAdch/TopologicalNumbers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KskAdch/TopologicalNumbers.jl)

## Overview

TopologicalNumbers.jl is a Julia package designed to compute topological numbers, such as the Chern numbers and $\mathbb{Z}_2$ numbers, 
using a numerical approach based on the Fukui-Hatsugai-Suzuki method [Fukui2005Chern,Mochol-Grzelak2018Efficient](@cite), 
the Shiozaki method [Fukui2007Quantum,Shiozaki2023discrete](@cite), 
or method of calculating the Weyl nodes [Hirayama2018Topological,Yang2011Quantum,Hirayama2015Weyl,Du2017Emergence](@cite).

This software is released under the MIT License, please see the LICENSE file for more details.  
It is confirmed to work on Julia 1.10 (LTS) and 1.11.


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



## Package Features

This package includes the following functions:

- Calculation of the dispersion relation.
- Provides numerical calculation methods for various types of topological numbers.
- Calculation of the phase diagram.
- Compute Pfaffian and tridiagonarize skew-symmetric matrix (migration to Julia from [PFAPACK](https://pypi.org/project/pfapack/) [Wimmer2012Algorithm](@cite)).
- Utility functions for plotting.
- Support parallel computing using `MPI`.


The correspondence between the spatial dimension of the system and the supported topological numbers is as follows.

```@raw html
<table>
    <tr>
        <th>Dimension</th>
        <th>Function</th>
    </tr>
    <tr>
        <td>1D</td>
        <td>Calculation of Berry Phases ($\mathbb{Z}$)<br></td>
    </tr>
    <tr>
        <td>2D</td>
        <td>Calculation of local Berry Fluxes ($\mathbb{Z}$)<br>Calculation of first Chern numbers ($\mathbb{Z}$)<br>Calculation of $\mathbb{Z}_{2}$ numbers ($\mathbb{Z}_{2}$)<br></td>
    </tr>
    <tr>
        <td>3D</td>
        <td>Calculation of Weyl nodes ($\mathbb{Z}$)<br> Calculation of first Chern numbers in sliced Surface ($\mathbb{Z}$)<br>Finding Weyl points ($\mathbb{Z}$)<br></td>
    </tr>
    <tr>
        <td>4D</td>
        <td>Calculation of second Chern numbers ($\mathbb{Z}$)<br></td>
    </tr>
</table>
```

## Contributing

User feedback is appreciated. Please create a [GitHub Issue](https://github.com/KskAdch/TopologicalNumbers.jl/issues) to report bugs or suggest new features, or a [Pull Request](https://github.com/KskAdch/TopologicalNumbers.jl/pulls) to contribute to the package. 
