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

`TopologicalNumbers.jl` is an open-source Julia package designed to calculate topological invariants—mathematical quantities that characterize the properties of materials—in condensed matter physics. 
These invariants, such as the Chern number and $\mathbb{Z}_2$ invariants, are crucial for understanding exotic materials like topological insulators and superconductors which have potential applications in advanced electronics, spintronics, and quantum computing. 
Our package provides an accessible and efficient toolset for researchers and educators to compute these invariants across various dimensions and symmetry classes, facilitating the exploration and discovery of new topological phases of matter.

TopologicalNumbers.jlは、凝縮系物理学における物質の性質を特徴付ける数学的量であるトポロジカル不変量を計算するための、オープンソースのJuliaパッケージです。
これらの不変量、例えばチャーン数や$\mathbb{Z}_2$不変量は、トポロジカル絶縁体や超伝導体のようなエキゾチックな物質を理解する上で重要であり、これらは先進的なエレクトロニクスやスピントロニクス、量子コンピューティングへの応用が期待されています[Ref.]。
本パッケージは、研究者や教育者がこれらの不変量を様々な次元や対称性クラスで計算できる、使いやすく効率的なツールセットを提供し、新たなトポロジカル相の探索と発見を促進します。



# Statement of need

固体物理学において、物質がどのような特性を持つか理解することは重要です。
例えば、電気伝導度は物質が電気をどの程度通すかを示す重要な物理量であり、通常は、物質に弱い静電場をかけたとき物質中の大域的空間（バルク）において電子が遷移可能であるような量子状態が存在すれば有限の値を持ち、電気が流れる金属となります。
逆に存在しなければ低い値を持ち、電気を流さない絶縁体となります。
しかし、1980年以降、今日に至るまで大きな関心を集めている量子ホール絶縁体やトポロジカル絶縁体などを含むトポロジカル電子系は、バルクでは絶縁体でありながら、物質表面には電子の伝導状態を持つ新しい物質状態を示します[Ref.]。
このような新奇特性から、トポロジカル電子系の候補物質の探索や特性評価がこれまで盛んに行われて来ています[@Hasan2010Colloquium;@Qi2011Topological]。


表面の伝導状態の特徴は、物質のハミルトニアンから決定される固有量子状態が波数空間上でどのようなトポロジーを持つかによって決定されており、第一チャーン数、第二チャーン数、$\mathbb{Z}_2$不変量などのトポロジカル数がそれらを特徴付けるために用いられます[Ref.]。
典型的な例として、2次元物質に対して弱い電場をかけた際、かけた電場と直交方向に量子化された有限の電気伝導度（ホール伝導度）が現れる現象を量子ホール効果と言います[Ref.]。
ホール伝導度$\sigma_{xy}$は整数値を取る第一チャーン数$\nu$によって特徴づけられ、$\sigma_{xy}=\frac{e^{2}}{h}\nu$と与えられることが知られています。
ここで、$e$は電気素量、$h$はプランク定数です。
他のトポロジカル数も同様に、系の次元や、系の持つ対称性のクラスに応じて、系を特徴付ける重要な物理量となります [@Ryu2010Topological]。


To obtain the topological numbers, we often need numerical calculations, and it may require an enormous amount of computation before convergence is achieved. 
Therefore, creating tools to easily compute these numbers will lead to advances in research on the topological phase of matter. 
So far, several methods have been reported that suggest that some topological numbers can be computed efficiently [@Fukui2005Chern;@Fukui2007Quantum;@Mochol-Grzelak2018Efficient;@Shiozaki2023discrete]. 
However, since each method is specialized for a specific dimension or symmetry class, it is necessary to implement the algorithm for each problem, respectively. 
Our project, `TopologicalNumbers.jl`, aims to provide a package that can easily and efficiently compute topological numbers comprehensively in various dimensions and symmetry classes.

トポロジカル数を得るためには、しばしば数値計算が必要であり、収束するまでに膨大な計算量を要することがあります。
したがって、これらの数を容易に計算できるツールを作成することは、物質のトポロジカル相の研究の進展につながります。
これまで、いくつかのトポロジカル数を効率的に計算できることを示唆する手法が報告されています[@Fukui2005Chern;@Fukui2007Quantum;@Mochol-Grzelak2018Efficient;@Shiozaki2023discrete]。
しかし、各手法は特定の次元や対称性クラスに特化しているため、それぞれの問題に対してアルゴリズムを実装する必要があります。
我々のプロジェクトであるTopologicalNumbers.jlは、様々な次元や対称性クラスにおいて、トポロジカル数を容易かつ効率的に包括的に計算できるパッケージを提供することを目指しています。


This package currently includes various methods for calculating topological numbers. 
The first is the Fukui-Hatsugai-Suzuki (FHS) method for computing first Chern numbers in two-dimensional solid-state systems [@Fukui2005Chern]. 
First Chern numbers are obtained by integrating the Berry curvature, derived from the eigenstates of the Hamiltonian, in the Brillouin zone. 
The FHS method enables us to compute the numbers efficiently by discretizing the Berry curvature in the Brillouin zone. 
Based on the FHS method, several calculation methods have been proposed to compute various topological numbers. 
One is the method of second Chern number calculation in four-dimensional systems [@Mochol-Grzelak2018Efficient]. 
$\mathbb{Z}_2$ invariants can also be calculated in two-dimensional systems with time-reversal symmetry [@Fukui2007Quantum;@Shiozaki2023discrete]. 
The FHS method is also applied to find Weyl points and Weyl nodes in three-dimensional systems [@Hirayama2018Topological;@Yang2011Quantum;@Hirayama2015Weyl;@Du2017Emergence].

このパッケージは現在、トポロジカル数を計算するための様々な手法を含んでいます。
最初に挙げられるのは、二次元固体物理系における第一チャーン数を計算するための福井・初貝・鈴木（FHS）法です [@Fukui2005Chern]。
第一チャーン数は、ハミルトニアンの固有状態から導かれるベリー曲率をブリルアンゾーン内で積分することで得られます。
FHS法は、ブリルアンゾーン内でベリー曲率を離散化することで、数を効率的に計算することを可能にします。
FHS法に基づいて、様々なトポロジカル数を計算するためのいくつかの計算法が提案されています。その一つは、四次元系における第二チャーン数の計算法です [@Mochol-Grzelak2018Efficient]。
また、時間反転対称性を持つ二次元系では、$\mathbb{Z}_2$不変量も計算できます [@Fukui2007Quantum;@Shiozaki2023discrete]。
FHS法は、三次元系におけるワイル点やワイルノードを見つけるためにも適用されています [@Hirayama2018Topological;@Yang2011Quantum;@Hirayama2015Weyl;@Du2017Emergence]。


Currently, there is no comprehensive Julia package that implements all these calculation methods. 
Various software packages are available for computing topological invariants in condensed matter systems. 
For instance, `Z2Pack` [@Gresch2017Z2Pack] is a Python-based tool widely used for calculating $\mathbb{Z}_2$ invariants and Chern numbers, primarily focusing on two-dimensional systems and requiring extensive setup. 
`WannierTools` [@Wu2018WannierTools] offers powerful features for analyzing topological materials using Wannier functions but is implemented in Fortran, which may present a steeper learning curve for some users.

現在、これらすべての計算法を実装した包括的なJuliaパッケージは存在しません。
一方、他のプラットフォームでは、本パッケージと異なるワニエ電荷中心[Ref.]またはWilson loop[Ref.]と呼ばれる量を用いた手法で実装されたいくつかのソフトウェアパッケージが利用可能です。
例えば、Z2Pack [@Gresch2017Z2Pack]は、$\mathbb{Z}_2$不変量や第一チャーン数の計算に広く用いられるPythonベースのツールです。
また、WannierTools [@Wu2018WannierTools]は、ワニエ関数を用いたトポロジカル物質の解析に強力な機能を提供しますが、Fortranで実装されており、一部のユーザーにとっては習得のハードルが高い可能性があります。


`TopologicalNumbers.jl` sets itself apart by providing a comprehensive and efficient framework within the Julia programming language, known for its high performance and user-friendly syntax. 
Our package supports a wide range of topological invariants—including first and second Chern numbers and $\mathbb{Z}_2$ invariants—across multiple dimensions and symmetry classes. 
It also offers parallel computing capabilities through `MPI.jl`, enhancing computational efficiency for large-scale problems. 
By combining these features, `TopologicalNumbers.jl` fills a gap in the current ecosystem of computational tools for topological materials, offering a unique balance of performance, usability, and extensibility.

TopologicalNumbers.jlは、高性能でユーザーフレンドリーな構文で知られるJuliaプログラミング言語内で、効率的なpure-Julia implementationのフレームワークを提供することで際立っています。
我々のパッケージは、第一チャーン数および第二チャーン数や$\mathbb{Z}_2$不変量など、様々なトポロジカル不変量を複数の次元や対称性クラスにわたってサポートしています。
また、MPI.jlを通じた並列計算機能を提供し、大規模な問題に対する計算効率を向上させます。
さらに、Juliaのマルチディスパッチ機能を活用し、`DifferentialEquations.jl`[Ref.]のような`Problem`、`Algorithm`、`solve`のような一貫したインターフェースを提供することで、拡張性を高めています。
これらの機能を組み合わせることで、TopologicalNumbers.jlは性能、使いやすさ、メンテナンス性、拡張性のユニークなバランスを提供し、他のライブラリと必ずしも競合するのではなく、代替のツールと異なる視点を提供します。


Additionally, for the calculation of $\mathbb{Z}_2$ invariants, which require the computation of Pfaffians, we have ported `PFAPACK` to Julia. 
`PFAPACK` is a Fortran/C++/Python library for computing the Pfaffian of skew-symmetric matrices [@Wimmer2012Algorithm], and our package includes a pure-Julia implementation of all the functions originally provided. 
While `SkewLinearAlgebra.jl` exists as an official Julia package for computing Pfaffians of real skew-symmetric matrices, TopologicalNumbers.jl is the first official package to offer a pure-Julia implementation for handling complex skew-symmetric matrices. 

また、Pfaffianの計算を必要とする$\mathbb{Z}_2$不変量の計算のために、我々はPFAPACKをJuliaに移植しました。
PFAPACKは、反対称行列のPfaffianを計算するためのFortran/C++/Pythonライブラリであり [@Wimmer2012Algorithm]、本パッケージには元々提供されていたすべての関数の純粋なJulia実装が含まれています。
SkewLinearAlgebra.jlは実対称行列のPfaffianを計算するための公式のJuliaパッケージとして存在しますが、TopologicalNumbers.jlは複素数の反対称行列を扱うための純粋なJulia実装を提供する最初の公式パッケージです。



# Usage

Users can easily calculate topological numbers using these methods included in our package.
In the simplest case, users only need to provide a function of the Hamiltonian matrix with wave numbers as arguments. 
Computations can be performed by creating a corresponding `Problem` and calling the `solve` function (`solve(Problem)`). 
The package also offers a `calcPhaseDiagram` function, enabling the computation of topological numbers in one-dimensional or two-dimensional parameter spaces by providing a `Problem` and parameter ranges (`calcPhaseDiagram(Problem, range...)`).

ユーザーは、本パッケージに含まれる数々の手法を用いて、容易にトポロジカル数を計算できます。
最も単純な場合、ユーザーは波数を引数とするハミルトニアン行列の関数を提供するだけで済みます。
計算は、対応するProblemを作成し、solve関数（solve(Problem)）を呼び出すことで実行できます。
また、本パッケージはcalcPhaseDiagram関数を提供しており、Problemとパラメータ範囲を指定することで、一次元または二次元のパラメータ空間におけるトポロジカル数の計算を可能にします（calcPhaseDiagram(Problem, range...)）。


Furthermore, several utility functions are available, such as `showBand`, `plot1D`, and `plot2D` for visualizing energy band structures and phase diagrams. 
We also provide various model Hamiltonians (e.g., `SSH`, `Haldane`) to enable users to quickly check the functionality and learn how to use these features.

さらに、エネルギーバンド構造や相図を可視化するためのshowBand、plot1D、plot2Dなどのユーティリティ関数も利用可能です。
また、ユーザーが機能を迅速にチェックし、これらの機能の使い方を学ぶことができるように、SSH、Haldaneなどの様々なモデルハミルトニアンも提供しています。



# Acknowledgements
The authors are grateful to Takahiro Fukui for fruitful discussions. 
M. K. was supported by JST, the establishment of university fellowships towards the creation of science technology innovation (Grant No. JPMJFS2107), and by JST SPRING (Grant No. JPMJSP2109).