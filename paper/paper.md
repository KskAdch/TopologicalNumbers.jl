---
title: 'TopologicalNumbers.jl: A Julia package for calculating topological numbers'
tags:
  - Julia
  - condensed matter physics
  - solid state physics
  - topological number
  - topological insulator
  - Berry phase
  - Chern number
  - Z2 number
  - phase diagram
  - Weyl node
  - Weyl point
  - pfaffian
  - skew-symmetric matrix
authors:
  - name: Keisuke Adachi
    orcid: 0009-0004-0195-7952
    equal-contrib: true
    corresponding: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Minoru Kanega
    orcid: 0009-0008-4623-8010
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    corresponding: true
    affiliation: 2
affiliations:
 - name: Department of Physics, Ibaraki University, Mito, Ibaraki, Japan
   index: 1
 - name: Department of Physics, Chiba University, Chiba, Japan
   index: 2
bibliography: paper.bib

---

# Summary
Topological insulators have been of considerable interest in the last decades [@Hasan2010Colloquium;@Qi2011Topological]. 
These materials are new states of matter that are insulating in the bulk but have conducting states on their surfaces.
The conducting states on the surface are protected by the topology of the bulk band structure, 
and topological numbers, such as first Chern number, second one, Z2 number, etc., 
are used to characterize them.
As a typical example, 
the quantum Hall effect has the quantized Hall conductivity, 
which can be calculated by the first Chern number. 
Other topological numbers similarly become important physical quantities that characterize the system, 
depending on the system dimension and symmetry classes [@Ryu2010Topological].

To obtain the topological numbers,
we often need numerical calculations,
and it may require an enormous amount of computation before convergence is achieved. 
Therefore, creating tools to easily calculate these numbers will lead to advances in research for topological phase of matters. 
So far, several methods [@Fukui2005Chern;@Fukui2007Quantum;@Mochol-Grzelak2018Efficient;@Shiozaki2023discrete] have been reported that suggest that some topological numbers can be computed efficiently. 
However, since each method is specialized for a specific dimension or symmetry class,
it is necessary to implement the algorithm for each problem, respectively.
Our project of `TopologicalNumbers.jl` aims to provide a package that can easily and efficiently compute topological numbers in various dimensions and symmetry classes, comprehensively.



# Statement of need
`TopologicalNumbers.jl` is an open-source Julia package for computing various topological numbers. 
This package currently includes various methods for computing the numbers.
The first one is the Fukui-Hatsugai-Suzuki (FHS) method [@Fukui2005Chern] for computing first Chern numbers in two-dimensional solid state systems.
First Chern numbers are obtained by integrating the Berry curvature, 
derived from the eigenstates of the Hamiltonian, in the Brillouin zone.
FHS method enables us to calculate the numbers efficiently by discretizing Berry curvature in the Brillouin zone.
Based on FHS method, several calculation methods have been proposed to compute various topological numbers. 
One is the method of second Chern number calculation in four-dimensional systems [@Mochol-Grzelak2018Efficient].
Z2 numbers can also be calculated in two-dimensional systems with time-reversal symmetry [@Fukui2007Quantum;@Shiozaki2023discrete].
FHS method is also applied to find Weyl points and Weyl node in three-dimensional systems [@Hirayama2018Topological;@Yang2011Quantum;@Hirayama2015Weyl;@Du2017Emergence].



Currently, there is no comprehensive Julia package that implements all these calculation methods. 
Users can easily calculate topological numbers using these methods. 
The simplest requirement for users is to provide a function of the Hamiltonian with wave numbers as arguments. 
By creating a corresponding `Problem` and calling the `solve` function (`solve(Problem)`), calculations can be executed. 
At present, each `Problem` has one implemented `Algorithm`, limiting users' choices in computational methods. 
However, `solve` is designed to accept an `Algorithm` as an argument (`solve(Problem, Algorithm)`), allowing for future extensions. 
For example, methods based on numerical integration are planned for implementation. 
The package also offers a `calcPhaseDiagram` function, enabling the computation of topological numbers in one-dimensional or two-dimensional parameter spaces by providing a `Problem` (`solve(Problem)`).



For the calculation of Z2 invariants, which require the computation of pfaffian, 
`PFAPACK` has been ported to Julia and is utilized in our package. 
`PFAPACK` is a Fortran/C++/Python library for calculating the pfaffian of skew-symmetric matrices, as noted in [@Wimmer2012Algorithm], and our package includes a pure-Julia implementation of all the functions originally provided. 
While `SkewLinearAlgebra.jl` exists as an official Julia package for computing pfaffians of real skew-symmetric matrices, 
`TopologicalNumbers.jl` is the first to offer a pure-Julia implementation for handling complex skew-symmetric matrices. 
Additionally, several utility functions are available, such as `showBand`/`plot1D`/`plot2D` for visualizing energy band structures and phase diagrams. 
We also provide various model Hamiltonians (e.g., `SSH`, `Haldane`) to enable users to quickly check the functionality and learn how to use these features. 
Moreover, the package supports parallel computing using `MPI.jl`.




<!-- There is no Julia package yet that comprehensively implements these calculation methods.  -->
<!-- The basic topological numbers in this package can be calculated if only the Hamiltonian is given.  -->
<!-- The calculations can be performed with a minimum number of arguments, 
making them easy to use even for Julia beginners and beginning students of condensed matter physics. 
It is also easy for researchers to use because it is designed with many optional arguments so that it can be used for general-purpose calculations. 
It is designed to be more accessible and with clear documentation. -->


# Acknowledgements
The authors are grateful to Takahiro Fukui for fruitful discussions.
M.K. was supported by JST, the establishment of university fellowships towards the creation of science technology innovation, Grant No. JPMJFS2107.


# References


<!-- 

This package includes the following functions:

- Calculation of the dispersion relation.
- Provides numerical calculation methods for various types of topological numbers.
- Calculation of the phase diagram.
- Compute Pfaffian and tridiagonarize skew-symmetric matrix (migration to Julia from [PFAPACK](https://pypi.org/project/pfapack/) [Wimmer2012Algorithm](@cite)).
- Utility functions for plotting.
- Support parallel computing using `MPI`.


The correspondence between the spatial dimension of the system and the supported topological numbers is as follows.


+-------------------+-----------------------------------------------------------------------+
| Dimension         | Function                                                              |
|                   |                                                                       |
+:=================:+:=====================================================================:+
| 0D                | - Calculation of Weyl nodes ($\mathbb{Z}$)                            |
+-------------------+-----------------------------------------------------------------------+
| 1D                | - Calculation of Berry Phases ($\mathbb{Z}$)                          |
+-------------------+-----------------------------------------------------------------------+
| 2D                | - Calculation of local Berry Fluxes ($\mathbb{Z}$)                    |
|                   | - Calculation of first Chern numbers ($\mathbb{Z}$)                   |
|                   | - Calculation of $\mathbb{Z}_2$ numbers ($\mathbb{Z}_2$)              |
+-------------------+-----------------------------------------------------------------------+
| 3D                | - Calculation of Weyl nodes ($\mathbb{Z}$)                            |
|                   | - Calculation of first Chern numbers in sliced Surface ($\mathbb{Z}$) |
|                   | - Finding Weyl points ($\mathbb{Z}$)                                  |
+-------------------+-----------------------------------------------------------------------+
| 4D                | - Calculation of second Chern numbers ($\mathbb{Z}$)                  |
+-------------------+-----------------------------------------------------------------------+



+-------------------+-----------------------------------------------------------------------+
| Dimension         | Function                                                              |
|                   |                                                                       |
+:=================:+:=====================================================================:+
| 0D                | Calculation of Weyl nodes ($\mathbb{Z}$)                              |
+-------------------+-----------------------------------------------------------------------+
| 1D                | Calculation of Berry Phases ($\mathbb{Z}$)                            |
+-------------------+-----------------------------------------------------------------------+
| 2D                | Calculation of local Berry Fluxes ($\mathbb{Z}$)                      |
|                   | Calculation of first Chern numbers ($\mathbb{Z}$)                     |
|                   | Calculation of $\mathbb{Z}_2$ numbers ($\mathbb{Z}_2$)                |
+-------------------+-----------------------------------------------------------------------+
| 3D                | Calculation of Weyl nodes ($\mathbb{Z}$)                              |
|                   | Calculation of first Chern numbers in sliced Surface ($\mathbb{Z}$)   |
|                   | Finding Weyl points ($\mathbb{Z}$)                                    |
+-------------------+-----------------------------------------------------------------------+
| 4D                | Calculation of second Chern numbers ($\mathbb{Z}$)                    |
+-------------------+-----------------------------------------------------------------------+ -->