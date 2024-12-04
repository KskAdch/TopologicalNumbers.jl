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

`TopologicalNumbers.jl` is an open-source Julia package designed to calculate topological invariants, which are mathematical quantities characterizing the properties of materialsin condensed matter physics. 
These invariants, such as the Chern number and the $\mathbb{Z}_2$ invariant, are crucial for understanding exotic materials like topological insulators and superconductors, which have potential applications in advanced electronics, spintronics, and quantum computing [Ref.]. 
This package provides researchers and educators with an easy-to-use and efficient toolset to compute these invariants across various dimensions and symmetry classes, facilitating the exploration and discovery of new topological phases of matter.

`TopologicalNumbers.jl`は、凝縮系物理学における物質の性質を特徴付ける数学的量であるトポロジカル不変量を計算するための、オープンソースのJuliaパッケージです。
これらの不変量、例えばチャーン数や$\mathbb{Z}_2$不変量は、先進的なエレクトロニクスやスピントロニクス、量子コンピューティングへの応用が期待されているトポロジカル絶縁体や超伝導体のようなエキゾチックな物質を理解する上で重要です[Ref.]。
本パッケージは、研究者や教育者がこれらの不変量を様々な次元や対称性クラスで計算できる、使いやすく効率的なツールセットを提供し、新たなトポロジカル相の探索と発見を促進します。


# Statement of need

Understanding the properties of materials is essential in solid-state physics. 
For example, electrical conductivity is a key physical quantity that indicates how well a material conducts electricity. 
Typically, when a weak electrostatic field is applied to a material, if there exist quantum states in the bulk where electrons can transition, the material exhibits finite electrical conductivity and behaves as a metal. 
Conversely, if such states do not exist, the electrical conductivity is low, and the material behaves as an insulator.
Since the 1980s, topological electronic systems—including quantum Hall insulators and topological insulators—have attracted significant attention because they exhibit new states where the bulk is insulating, but the material’s surface possesses conducting electronic states [@Hasan2010Colloquium;@Qi2011Topological]. 
Due to these novel properties, extensive research has been conducted to explore candidate materials and evaluate their characteristics.

固体物理学において、物質がどのような特性を持つか理解することは重要です。
例えば、電気伝導度は物質が電気をどの程度通すかを示す重要な物理量です。
通常、物質に弱い静電場をかけたとき、物質中の大域的空間（バルク）において電子が遷移可能な量子状態が存在すれば、電気伝導度は有限の値を持ち、電気が流れる金属となります。
逆に、存在しなければ電気伝導度は低い値となり、電気を流さない絶縁体となります。
しかし、1980年以降、量子ホール絶縁体やトポロジカル絶縁体などを含むトポロジカル電子系は、バルクでは絶縁体でありながら、物質表面には電子の伝導状態を持つ新しい物質状態を示し、今日に至るまで大きな関心を集めています[@Hasan2010Colloquium;@Qi2011Topological]。
この新奇な特性から、トポロジカル電子系の候補物質の探索や特性評価が盛んに行われてきました。


The features of surface conducting states are determined by the topology of eigen quantum states in momentum space, dictated by the material’s Hamiltonian. 
Topological numbers, such as the first Chern number, second Chern number, and the $\mathbb{Z}2$ invariant, are used to characterize these properties [@Thouless1982Quantized;@Kane2005Z_2]. 
A typical example is the quantum Hall effect, where applying a weak electric field to a two-dimensional material results in a quantized finite electrical conductivity (Hall conductivity) in the direction perpendicular to the applied field [@Thouless1982Quantized]. 
The Hall conductivity $\sigma{xy}$ is characterized by the integer-valued first Chern number $\nu$ and is given by $\sigma_{xy} = \frac{e^{2}}{h} \nu$, where $e$ is the elementary charge and $h$ is Planck’s constant. 
Other topological numbers similarly serve as important physical quantities characterizing systems, depending on their dimensions and symmetry classes [@Ryu2010Topological].

表面の伝導状態の特徴は、物質のハミルトニアンから決定される固有量子状態が波数空間上でどのようなトポロジーを持つかによって決まり、第一チャーン数、第二チャーン数、$\mathbb{Z}_2$不変量などのトポロジカル数がそれらを特徴付けるために用いられます[@Thouless1982Quantized;@Kane2005Z_2]。
典型的な例として、2次元物質に弱い電場をかけた際、かけた電場と直交する方向に量子化された有限の電気伝導度（ホール伝導度）が現れる現象を量子ホール効果と言います[@Thouless1982Quantized]。
ホール伝導度$\sigma_{xy}$は、整数値を取る第一チャーン数$\nu$によって特徴づけられ、$\sigma_{xy}=\frac{e^{2}}{h}\nu$と与えられることが知られています。
ここで、$e$は電気素量、$h$はプランク定数です。
他のトポロジカル数も同様に、系の次元や持つ対称性クラスに応じて、系を特徴付ける重要な物理量となります [@Ryu2010Topological]。


To obtain the topological numbers, we often need numerical calculations, and it may require an enormous amount of computation before convergence is achieved. 
Therefore, creating tools to easily compute these numbers will lead to advances in research on the topological phase of matter. 
So far, several methods have been reported that suggest that some topological numbers can be computed efficiently [@Fukui2005Chern;@Fukui2007Quantum;@Mochol-Grzelak2018Efficient;@Shiozaki2023discrete]. 
However, since each method is specialized for specific dimensions or symmetry classes, it is necessary to implement the algorithm for each problem, respectively. 
Our project, `TopologicalNumbers.jl`, aims to provide a package that can easily and efficiently compute topological numbers comprehensively in various dimensions and symmetry classes.


This package currently includes several methods for calculating topological numbers. 
The first is the Fukui-Hatsugai-Suzuki (FHS) method for computing the first Chern number in two-dimensional solid-state systems [@Fukui2005Chern]. 
The first Chern number is obtained by integrating the Berry curvature, derived from the eigenstates of the Hamiltonian, in the Brillouin zone. 
The FHS method enables us to compute the numbers efficiently by discretizing the Berry curvature in the Brillouin zone. 
Based on the FHS method, several calculation methods have been proposed to compute various topological numbers. 
One such method is for calculating the second Chern number calculation in four-dimensional systems [@Mochol-Grzelak2018Efficient]. 
The $\mathbb{Z}_2$ invariant can also be calculated in two-dimensional systems with time-reversal symmetry [@Fukui2007Quantum;@Shiozaki2023discrete]. 
The FHS method is also applied to find Weyl points and Weyl nodes in three-dimensional systems [@Hirayama2018Topological;@Yang2011Quantum;@Hirayama2015Weyl;@Du2017Emergence].


Currently, there is no comprehensive Julia package that implements all these calculation methods. 
On other platforms, software packages utilizing different approaches—such as methods based on Wannier charge centers [@Soluyanov2011Computing] or Wilson loops [@Yu2011Equivalent]—are available. 
For example, `Z2Pack` [@Gresch2017Z2Pack] is a Python-based tool widely used for calculating the $\mathbb{Z}_2$ invariant and the first Chern number. 
`WannierTools` [@Wu2018WannierTools] offers powerful features for analyzing topological materials using Wannier functions but is implemented in Fortran, which may present a steep learning curve for some users.

現在、これらすべての計算法を実装した包括的なJuliaパッケージは存在しません。
一方、他のプラットフォームでは、本パッケージと異なるワニエ電荷中心[@Soluyanov2011Computing]またはWilson ループ[@Yu2011Equivalent]と呼ばれる量を用いた手法で実装されたソフトウェアパッケージが利用可能です。
例えば、`Z2Pack` [@Gresch2017Z2Pack]は、$\mathbb{Z}_2$不変量や第一チャーン数の計算に広く用いられるPythonベースのツールです。
また、`WannierTools` [@Wu2018WannierTools]は、ワニエ関数を用いたトポロジカル物質の解析に強力な機能を提供しますが、Fortranで実装されており、一部のユーザーにとっては習得のハードルが高い可能性があります。


`TopologicalNumbers.jl` distinguishes itself by providing an efficient, pure Julia implementation framework within the Julia programming language, known for its high performance and user-friendly syntax. 
Our package supports various topological invariants—including the first and second Chern numbers and the $\mathbb{Z}_2$ invariant—across multiple dimensions and symmetry classes. 
It also offers parallel computing capabilities through `MPI.jl`, enhancing computational efficiency for large-scale problems. 
Furthermore, by leveraging Julia’s multiple dispatch feature, we adopt a consistent interface using the `Problem`, `Algorithm`, and `solve` style, similar to `DifferentialEquations.jl` [@Rackauckas2017DifferentialEquationsjl], enhancing extensibility. 
By combining these features, `TopologicalNumbers.jl` offers a unique balance of performance, usability, maintainability, and extensibility, providing an alternative perspective rather than competing directly with other libraries.

`TopologicalNumbers.jl`は、高性能でユーザーフレンドリーな構文で知られるJuliaプログラミング言語内で、効率的な純粋Julia実装のフレームワークを提供することで際立っています。
我々のパッケージは、第一チャーン数および第二チャーン数、$\mathbb{Z}_2$不変量など、様々なトポロジカル不変量を複数の次元や対称性クラスにわたってサポートしています。
また、`MPI.jl`を通じた並列計算機能を提供し、大規模な問題に対する計算効率を向上させます。
さらに、Juliaのマルチディスパッチ機能を活用して、`DifferentialEquations.jl`[@Rackauckas2017DifferentialEquationsjl]のような`Problem`、`Algorithm`、`solve`スタイルを採用することで一貫したインターフェースを提供し、拡張性を高めています。
これらの機能を組み合わせることで、`TopologicalNumbers.jl`は性能、使いやすさ、メンテナンス性、拡張性のユニークなバランスを提供し、他のライブラリと競合するのではなく、代替のツールとして異なる視点を提供します。


Additionally, to compute the $\mathbb{Z}_2$ invariant, which requires calculating the Pfaffian, we have ported `PFAPACK` to Julia. 
`PFAPACK` is a Fortran/C++/Python library for computing the Pfaffian of skew-symmetric matrices [@Wimmer2012Algorithm], and our package includes pure Julia implementations of all originally provided functions. 
While `SkewLinearAlgebra.jl` exists as an official Julia package for computing the Pfaffian of real skew-symmetric matrices, `TopologicalNumbers.jl` is the first official package to offer a pure Julia implementation handling complex skew-symmetric matrices. 


# Usage

Users can easily compute topological numbers using the verious methods included in this package.
In the simplest case, users only need to provide a function of the Hamiltonian matrix with wave numbers as arguments. 
Computations are performed by creating the corresponding `Problem` and calling the `solve` function (`solve(Problem)`). 
The package also provides the `calcPhaseDiagram` function, which allows the computation of topological numbers in one-dimensional or two-dimensional parameter spaces by specifying the `Problem` and parameter ranges (`calcPhaseDiagram(Problem, range...)`).


Furthermore, utility functions such as `showBand`, `plot1D`, and `plot2D` are available for visualizing energy band structures and phase diagrams. 
We also offer various model Hamiltonians (e.g., `SSH`, `Haldane`) enabling users to quickly test functionalities and learn how to use these features.


# Acknowledgements
The authors are grateful to Takahiro Fukui for fruitful discussions. 
M. K. was supported by JST, the establishment of university fellowships towards the creation of science technology innovation (Grant No. JPMJFS2107), and by JST SPRING (Grant No. JPMJSP2109).