---
title: 'TopologicalNumbers.jl: A Julia package for topological number computation'
tags:
  - Julia
  - condensed matter physics
  - solid state physics
  - topological number
  - topological insulator
  - Berry phase
  - Chern number
  - $\mathbb{Z}_2$ number
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

`TopologicalNumbers.jl` is an open-source Julia package designed to calculate topological invariants — mathematical quantities that characterize the properties of materials in condensed matter physics. 
These invariants, such as the Chern number and the $\mathbb{Z}_2$ invariant, are crucial for understanding exotic materials like topological insulators and superconductors, which have potential applications in advanced electronics, spintronics, and quantum computing [@Nayak2008NonAbelian;@Hasan2010Colloquium;@Qi2011Topological]. 
This package provides researchers and educators with an easy-to-use and efficient toolset to compute these invariants across various dimensions and symmetry classes, facilitating the exploration and discovery of new topological phases of matter.


# Statement of need

Understanding the properties of materials is essential in solid-state physics. 
For example, electrical conductivity is a key physical quantity indicating how well a material conducts electric current. 
Typically, when a weak electric field is applied to a material, if quantum eigenstates exist in the bulk into which electrons can transition, the material exhibits finite electrical conductivity and behaves as a metal.
Conversely, if such states do not exist, the electrical conductivity is low, and the material behaves as an insulator.
Since the 1980s, it has been revealed that certain materials exhibit states in which the bulk is insulating while the material’s surface has conducting electronic states [@Hasan2010Colloquium;@Qi2011Topological]. 
These materials are known as topological electronic systems, including quantum Hall insulators and topological insulators. 
Due to these novel properties, extensive research has been conducted to identify candidate materials and evaluate their characteristics.


The features of surface conducting states are determined by the topology of quantum eigenstates in momentum space. 
Topological numbers, such as the first Chern number, the second Chern number, and the $\mathbb{Z}_2$ invariant, are used to characterize these properties [@Thouless1982Quantized;@Kane2005Z_2]. 
A typical example is the quantum Hall effect, where applying a weak electric field to a two-dimensional material results in a quantized finite electrical conductivity (Hall conductivity) perpendicular to the applied field [@Thouless1982Quantized]. 
The Hall conductivity $\sigma_{xy}$ is characterized by the first Chern number $\nu \in \mathbb{Z}$ and is given by $\sigma_{xy} = \frac{e^{2}}{h} \nu$, where $e$ is the elementary charge and $h$ is Planck’s constant. 
Other topological numbers similarly serve as important physical quantities that characterize systems, depending on their dimensions and symmetry classes [@Ryu2010Topological].


Obtaining topological numbers often requires extensive numerical calculations, which may demand considerable computational effort before achieving convergence. 
Therefore, creating tools that simplify the computation of these quantities will advance research on topological phases of matter. 
Several methods have been developed to enable efficient computation of certain topological numbers [@Fukui2005Chern;@Fukui2007Quantum;@Mochol-Grzelak2018Efficient;@Shiozaki2023discrete]. 
However, since each method is typically specialized for specific dimensions or symmetry classes, one must often implement algorithms separately for each problem. 
Our project, `TopologicalNumbers.jl`, aims to provide a package that can efficiently and easily compute topological numbers across various dimensions and symmetry classes.


This package currently includes several methods for calculating topological numbers. 
The first is the Fukui–Hatsugai–Suzuki (FHS) method for computing the first Chern number in two-dimensional solid-state systems [@Fukui2005Chern]. 
The first Chern number is obtained by integrating the Berry curvature, derived from the Hamiltonian’s eigenstates, over the Brillouin zone. 
The FHS method enables efficient computation by discretizing the Berry curvature in the Brillouin zone. 
Several methods based on the FHS approach have been proposed to compute various topological numbers.
One such method calculates the second Chern number in four-dimensional systems [@Mochol-Grzelak2018Efficient]. 
The $\mathbb{Z}_2$ invariant can be computed in two-dimensional systems with time-reversal symmetry [@Fukui2007Quantum;@Shiozaki2023discrete]. 
The FHS method also applies to identifying Weyl points and Weyl nodes in three-dimensional systems [@Yang2011Quantum;@Hirayama2015Weyl;@Du2017Emergence;@Hirayama2018Topological].


Currently, there is no comprehensive Julia package that implements all these calculation methods. 
On other platforms, software packages using different approaches — such as those based on Wannier charge centers [@Soluyanov2011Computing] or Wilson loops [@Yu2011Equivalent] — are available. 
For example, `Z2Pack` [@Gresch2017Z2Pack] is a Python-based tool widely used for calculating the $\mathbb{Z}_2$ invariant and the first Chern number. 
`WannierTools` [@Wu2018WannierTools] offers powerful features for analyzing topological materials but is implemented in Fortran, which may pose a steep learning curve for some users.


`TopologicalNumbers.jl` distinguishes itself by providing an efficient, pure Julia implementation. Julia is known for its high performance and user-friendly syntax. 
This package supports various topological invariants across multiple dimensions and symmetry classes, including the first and second Chern numbers and the $\mathbb{Z}_2$ invariant. 
It also offers parallel computing capabilities through `MPI.jl`, enhancing computational efficiency for large-scale problems. 
By leveraging Julia’s multiple dispatch feature, we adopt a consistent interface using the `Problem`, `Algorithm`, and `solve` style — similar to `DifferentialEquations.jl` [@Rackauckas2017DifferentialEquationsjl] — to improve extensibility. 
With these features, `TopologicalNumbers.jl` achieves a unique balance of performance, usability, maintainability, and extensibility, providing an alternative perspective rather than directly competing with other libraries.


Additionally, to compute the $\mathbb{Z}_2$ invariant, which requires calculating the Pfaffian, we have ported `PFAPACK` to Julia. 
`PFAPACK` is a Fortran/C++/Python library for computing the Pfaffian of skew-symmetric matrices [@Wimmer2012Algorithm].
Our package includes pure Julia implementations of all originally provided functions. 
While `SkewLinearAlgebra.jl` exists as an official Julia package for computing the Pfaffian of real skew-symmetric matrices, to our knowledge, `TopologicalNumbers.jl` is the first official package to offer a pure Julia implementation that handles complex skew-symmetric matrices. 


# Usage

Users can easily compute topological numbers using the various methods included in this package.
In the simplest case, they need only provide a function that returns the Hamiltonian matrix as a function of the wave numbers. 
Computations are performed by creating the corresponding `Problem` instance and calling the `solve` function (`solve(Problem)`). 
The package also provides the `calcPhaseDiagram` function, which enables the computation of topological numbers in one-dimensional or two-dimensional parameter spaces by specifying the `Problem` and parameter ranges (`calcPhaseDiagram(Problem, range...)`).


Moreover, utility functions such as `showBand`, `plot1D`, and `plot2D` are available for visualizing energy band structures and phase diagrams. 
We also offer various model Hamiltonians, such as the Su–Schrieffer–Heeger (SSH) model [@Su1979Solitons] and the Haldane model [@Haldane1988Model], allowing users to quickly test functionalities and learn how to use these features.


# Acknowledgements
The authors are grateful to Takahiro Fukui for fruitful discussions. 
M. K. was supported by JST, the establishment of university fellowships towards the creation of science technology innovation (Grant No. JPMJFS2107), and by JST SPRING (Grant No. JPMJSP2109).
