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
and the topological numbers, such as the first Chern number, second one, Z2 number, etc., 
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
So far, several methods [@Fukui2005;@Fukui2007Quantum;@Mochol-Grzelak2018Efficient;@Shiozaki2023discrete] have been reported that suggest that some topological numbers can be computed efficiently. 
However, since each method is specialized for a specific dimension or symmetry class,
it is necessary to implement the algorithm for each problem.
Our project of `TopologicalNumbers.jl` aims to provide a package that can easily and efficiently compute topological numbers in various dimensions and symmetry classes, comprehensively.



# Statement of need
`TopologicalNumbers.jl` is an open-source Julia package for computing various topological numbers. 
These numbers are very important physical quantities in the study of condensed matter physics. 
However, their calculation can be very computationally intensive. 
In addition, the computation of topological numbers is often not the essence of the study. 
Most researchers in the field of condensed matter physics calculate topological numbers as an assistance to their research. 
Therefore, this package is very useful and practical for them. 
The package is also useful for beginning students of condensed matter physics.

There are currently various methods for computing topological numbers, 
the first of efficient calculation methods was the Fukui-Hatsugai-Suzuki [@Fukui2005] in 2005 for method of computing two-dimensional Chern numbers. 
In general, topological numbers are obtained by integrating the eigenstates of the Hamiltonian in the Brillouin zone, 
and this method can efficiently calculate them by discretizing the Brillouin zone, 
which is the integral range. 
This method can be useful in a practical computation for more complicated systems with a topological order for which a number of data points of the wave functions cannot easily be increased. 
Also, methods have been proposed to compute various topological invariants using this method. 
One is the method of [@Shiozaki2023discrete], which computes the Z2 numbers in two-dimensional with time-reversal symmetry in 2023. 
This method does not require any gauge fixing conditions and is quantized for any discrete approximation of the Brillouin zone. 
It is also used for methods to find Weyl points and Weyl node [@Hirayama2017;@Yang2011;@Hirayama2015;@Du2017] in three-dimensional.

There is no Julia package yet that comprehensively implements these methods. 
The basic topological numbers in this package can be calculated if only the Hamiltonian is given. 
The calculations can be performed with a minimum number of arguments, 
making them easy to use even for Julia beginners and beginning students of condensed matter physics. 
It is also easy for researchers to use because it is designed with many optional arguments so that it can be used for general-purpose calculations. 
It is designed to be more accessible and with clear documentation.


# Acknowledgements
The authors are grateful to Takahiro Fukui for fruitful discussions.
M.K. was supported by JST, the establishment of university fellowships towards the creation of science technology innovation, Grant No. JPMJFS2107.


# References




Non-negative matrix factorization (NMF) is a widely used algorithm to decompose non-negative matrix data into factor matrices. Due to the interpretability of its non-negativity and the convenience of using decomposition results as clustering, there are many applications of NMF in image processing, audio processing, and bioinformatics (Cichocki et al., 2009).
NMF has been applied to matrix data but there is a growing demand to apply NMF to more heterogeneous non-negative data such as multiple matrices and tensors (high-dimensional arrays), which are higher-order data structures than matrices (Cichocki et al., 2009). To meet these requirements, we originally developed nnTensor, which is an R/CRAN package to perform some non-negative matrix/tensor decomposition algorithms (https://cran.r-project.org/web/packages/nnTensor/index.html).