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



There is no Julia package yet that comprehensively implements these calculation methods. 
ユーザーはこれらのメソッドを用いて、簡単にトポロジカル数を計算することができます。
最も簡単には、ユーザーが与えるべき情報は波数を引数に持つハミルトニアンの関数のみです。
対応する`Problem`を作成し、`solve`関数を呼び出す (`solve(Problem)`) ことで、計算が実行されます。
現在、それぞれの`Problem`に対して実装されている`Algorithm`が1種類ずつのため、ユーザーは計算手法を選ぶことができませんが、
将来的な拡張性のために`solve`は`Algorithm`を引数に取ることができます (`solve(Problem, Algorithm)`)。
例えば、数値積分による計算手法を実装予定です。
また、相図を計算する`calcPhaseDiagram`関数も提供しており、`Problem`を与えることで1次元/2次元パラメータ空間におけるトポロジカル数の値を計算することができます(`solve(Problem)`)。



また、Z2不変量の計算に必要なpfaffianの計算については、`PFAPACK`をJuliaに移植し、利用しています。
PFAPACKは、skew-symmetric行列のpfaffianを計算するためのFortran/C++/Pythonライブラリであり [@Wimmer2012Algorithm]、
元々提供されているすべての関数のpure-Julia実装が我々のパッケージに含まれています。
現在までにreal skew-symmetric matricesに対するpfaffianを計算するJulia公式パッケージとして`SkewLinearAlgebra.jl`がありますが、
complex skew-symmetric matricesを取り扱えるpure-Juliaの公式ライブラリは`TopologicalNumbers.jl`が初めてです。
エネルギーバンド構造や相図を可視化する`showBand`/`plot1D`/`plot2D`、
ユーザーがこれらの機能の動作をすぐに確認し、使い方を学ぶことができるようないくつかのモデルハミルトニアン(`SSH`/`Haldane`など)を提供するなど、いくつかのユーティリティ関数も利用可能です。
また、`MPI.jl`を用いた並列計算にも対応しています。




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