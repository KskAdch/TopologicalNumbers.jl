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
These materials show new states of matter that are insulating in the bulk but have conducting states on their surfaces. 
The conducting states on the surface are protected by the topology of the bulk band structure, and topological numbers, 
such as first Chern number, second Chern number, Z2 number, etc., are used to characterize them. 
As a typical example, the quantum Hall effect has the quantized Hall conductivity, which can be calculated by the first Chern number. 
Other topological numbers similarly become important physical quantities that characterize the system, 
depending on the system dimension and symmetry classes [@Ryu2010Topological].

To obtain the topological numbers, we often need numerical calculations, 
and it may require an enormous amount of computation before convergence is achieved. 
Therefore, creating tools to easily calculate these numbers will lead to advances in research for the topological phase of matters. 
So far, several methods have been reported that suggest that some topological numbers can be computed efficiently [@Fukui2005Chern;@Fukui2007Quantum;@Mochol-Grzelak2018Efficient;@Shiozaki2023discrete]. 
However, since each method is specialized for a specific dimension or symmetry class, 
it is necessary to implement the algorithm for each problem, respectively. 
Our project, `TopologicalNumbers.jl`, aims to provide a package that can easily and efficiently compute topological numbers comprehensively in various dimensions and symmetry classes.



# Statement of need
`TopologicalNumbers.jl` is an open-source Julia package for computing various topological numbers. 
This package currently includes various methods for computing topological numbers. 
The first is the Fukui-Hatsugai-Suzuki (FHS) method for computing first Chern numbers in two-dimensional solid-state systems [@Fukui2005Chern]. 
First Chern numbers are obtained by integrating the Berry curvature, 
derived from the eigenstates of the Hamiltonian, 
in the Brillouin zone. 
The FHS method enables us to calculate the numbers efficiently by discretizing Berry curvature in the Brillouin zone. 
Based on the FHS method, several calculation methods have been proposed to compute various topological numbers. 
One is the method of second Chern number calculation in four-dimensional systems [@Mochol-Grzelak2018Efficient]. 
Z2 numbers can also be calculated in two-dimensional systems with time-reversal symmetry [@Fukui2007Quantum;@Shiozaki2023discrete]. 
The FHS method is also applied to find Weyl points and Weyl nodes in three-dimensional systems [@Hirayama2018Topological;@Yang2011Quantum;@Hirayama2015Weyl;@Du2017Emergence].



Currently, there is no comprehensive Julia package that implements all these calculation methods. 
Users can easily calculate topological numbers using these methods included in our package. 
In the simplest case, users only need to provide a function of the Hamiltonian with wave numbers as arguments. 
Calculations can be executed by creating a corresponding `Problem` and calling the `solve` function (`solve(Problem)`). 
The package also offers a `calcPhaseDiagram` function, 
enabling the computation of topological numbers in one-dimensional or two-dimensional parameter spaces by providing a `Problem` (`solve(Problem)`).



For the calculation of Z2 invariants, which require the computation of pfaffian, we have ported `PFAPACK` to Julia. 
`PFAPACK` is a Fortran/C++/Python library for calculating the pfaffian of skew-symmetric matrices [@Wimmer2012Algorithm], 
and our package includes a pure-Julia implementation of all the functions originally provided. 
While `SkewLinearAlgebra.jl` exists as an official Julia package for computing pfaffians of real skew-symmetric matrices, 
`TopologicalNumbers.jl` is the first official package to offer a pure-Julia implementation for handling complex skew-symmetric matrices. 
Additionally, several utility functions are available, such as `showBand`/`plot1D`/`plot2D` for visualizing energy band structures and phase diagrams. 
We also provide various model Hamiltonians (e.g., `SSH`, `Haldane`) to enable users to quickly check the functionality and learn how to use these features. 
Moreover, the package supports parallel computing using `MPI.jl`. 
Consequently, `TopologicalNumbers.jl` is the first comprehensive Julia package for computing topological numbers in solid-state systems, 
and we believe that it will be useful for researchers in the field of solid-state physics.



# Acknowledgements
The authors are grateful to Takahiro Fukui for fruitful discussions. 
M.K. was supported by JST, 
the establishment of university fellowships towards the creation of science technology innovation, 
Grant No. JPMJFS2107.