abstract type TopologicalNumbersAlgorithms end

# Algorithms for calculating the Berry phase
abstract type BerryPhaseAlgorithms <: TopologicalNumbersAlgorithms end
# struct Int1DBP <: Z2Algorithms end
struct BP <: BerryPhaseAlgorithms end

# Algorithms for calculating the first Chern number
abstract type FirstChernAlgorithms <: TopologicalNumbersAlgorithms end
# struct IntFChern <: FirstChernAlgorithms end
struct FHS <: FirstChernAlgorithms end

# Algorithms for calculating the second Chern number
abstract type SecondChernAlgorithms <: TopologicalNumbersAlgorithms end
# struct IntSChern <: SecondChernAlgorithms end
struct FHS2 <: SecondChernAlgorithms end

# Algorithms for calculating the Z2 invariant
abstract type Z2Algorithms <: TopologicalNumbersAlgorithms end
# struct Int2DZ2 <: Z2Algorithms end
struct Shio <: Z2Algorithms end

# Algorithms for calculating the local Berry flux
abstract type BerryFluxAlgorithms <: TopologicalNumbersAlgorithms end
struct FHSlocal2 <: BerryFluxAlgorithms end

# Algorithms for finding and calculating the Weyl points
abstract type WeylPointsAlgorithms <: TopologicalNumbersAlgorithms end
struct FHSsurface <: WeylPointsAlgorithms end
struct FHSlocal3 <: WeylPointsAlgorithms end
struct Evar <: WeylPointsAlgorithms end