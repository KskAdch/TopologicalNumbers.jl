---
title: 'TopologicalNumbers.jl: A Julia package for calculating topological numbers'
tags:
  - Julia
  - 
  - 
  - 
  - 
authors:
  - name: Keisuke Adachi
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    corresponding: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Minoru Kanega
    orcid: 0009-0008-4623-8010
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    corresponding: true
    affiliation: 2
#   - given-names: Ludwig
#     dropping-particle: van
#     surname: Beethoven
#     affiliation: 3
affiliations:
 - name: Department of Physics, Ibaraki University, Mito, Ibaraki, Japan
   index: 1
 - name: Department of Physics, Chiba University, Chiba, Japan
   index: 2
#  - name: Independent Researcher, Country
#    index: 3
# date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
The topological number is a physical quantity that plays an important role in characterizing phase transitions in condensed matter physics. Recently, topological phase transitions have been of considerable
interest in recent condensed matter physics. A typical example is the integer quantum Hall transition, where the quantized Hall conductivity can be calculated by the Chern number, one of the topological numbers. The quantization of such physical quantities is closely related to the concept of topology. However, the calculation of this topological number often requires numerical calculations. In addition, it may require an enormous amount of computation before convergence is achieved. Creating tools to easily calculate such topological numbers will lead to advances in condensed matter physics.


# Statement of need
`TopologicalNumbers.jl` is a package for computing various topological invariants. These quantities are very important physical quantities in the study of condensed matter physics. However, their computation can be very computationally intensive. In addition, the calculation of topological invariants is often not the essence of the research. Therefore, this package is very useful for researchers in the field of condensed matter physics to easily and efficiently calculate topological invariants.

There are currently various methods for computing topological invariants, the first of which was the Fukui-Hatsugai-Suzuki method [@Fukui2005] in 2005 for computing two-dimensional Chern numbers. In general, topological invariants are obtained by integrating the eigenstates of the Hamiltonian in the Brillouin zone, and this method can efficiently calculate them by discretizing the range of integration. This method can be useful in a practical computation for more complicated systems with a topological order for which a number of data points of the wave functions cannot easily be increased. Also, methods have been proposed to compute various topological invariants using this method. One is the method of [@Shiozaki2023discrete], which computes the Z2 number in two dimensions in 2023. This method does not require any gauge fixing conditions and is quantized for any discrete approximation of the Brillouin zone. It is also used for methods to find three-dimensional Weil points [@Hirayama2017,@Yang2011,@Hirayama2015,@Du2017].

There is no Julia package yet that comprehensively implements these methods. This package is easy for beginners to use because calculations can be done with a minimum of arguments. It is also easy for researchers to use because it is designed with many optional arguments so that it can be used for general-purpose calculations. It is designed to be more accessible and with clear documentation.


# Mathematics



# Citations



# Figures



# Acknowledgements
The authors are grateful to Takahiro Fukui for fruitful discussions.
M.K. was supported by JST, the establishment of university fellowships towards the creation of science technology innovation, Grant No. JPMJFS2107.


# References
