module TopologicalNumbers

include("./pfaffian.jl")
export pfaffian

# algorithms
include("./algorithms.jl")
export BP, FHS, FHS2, Shio
export FHSlocal2, FHSsurface, FHSlocal3, Evar

# problems
include("./problems.jl")
export BPProblem, FCProblem, SCProblem, Z2Problem
export LBFProblem, WCSProblem, WNProblem, WPProblem

# solutions
export BPSolution, FCSolution, SCSolution, Z2Solution
export LBFSolution, WCSSolution, WNSolution, WPSolution

# parallel environment
include("./parallelEnv.jl")
export UseSingleThread, UseMPI

# models
include("./models.jl")
export SSH, KitaevChain
export Flux2d, Haldane, KitaevHoneycomb
export ThoulessPump, KaneMele, BHZ
export LatticeDirac


include("./packages.jl")
include("./params.jl")
include("./showBand.jl")
include("./calcZ2.jl")
include("./calcChern.jl")
include("./calcBerryPhase.jl")
include("./phaseDiagram.jl")
include("./plot.jl")
include("./calcBerryFlux.jl")
include("./calcWeylNode.jl")
include("./calcChernSurface.jl")
include("./findWeylPoint.jl")
include("./SecondChern.jl")

# main functions
export solve
export calcPhaseDiagram
export showBand
export plot1D, plot2D

# Old methods
export calcZ2, calcChern, calcSecondChern, calcBerryPhase
export calcBerryFlux
export calcWeylNode, calcChernSurface, findWeylPoint

end
