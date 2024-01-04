@doc raw"""
"""
abstract type TopologicalNumbersAlgorithms end

# Algorithms for calculating the Berry phase
@doc raw"""
"""
abstract type BerryPhaseAlgorithms <: TopologicalNumbersAlgorithms end
# struct Int1DBP <: Z2Algorithms end

@doc raw"""
"""
struct BP <: BerryPhaseAlgorithms end

# Algorithms for calculating the first Chern number
@doc raw"""
"""
abstract type FirstChernAlgorithms <: TopologicalNumbersAlgorithms end
# struct IntFChern <: FirstChernAlgorithms end

@doc raw"""
"""
struct FHS <: FirstChernAlgorithms end


# Algorithms for calculating the second Chern number
@doc raw"""
"""
abstract type SecondChernAlgorithms <: TopologicalNumbersAlgorithms end
# struct IntSChern <: SecondChernAlgorithms end

@doc raw"""
"""
struct FHS2 <: SecondChernAlgorithms end

# Algorithms for calculating the Z2 invariant
@doc raw"""
"""
abstract type Z2Algorithms <: TopologicalNumbersAlgorithms end
# struct Int2DZ2 <: Z2Algorithms end

@doc raw"""
"""
struct Shio <: Z2Algorithms end


# Algorithms for calculating the local Berry flux
@doc raw"""
"""
abstract type BerryFluxAlgorithms <: TopologicalNumbersAlgorithms end

@doc raw"""
"""
struct FHSlocal2 <: BerryFluxAlgorithms end


# Algorithms for finding and calculating the Weyl points
@doc raw"""
"""
abstract type WeylPointsAlgorithms <: TopologicalNumbersAlgorithms end

@doc raw"""
"""
struct FHSsurface <: WeylPointsAlgorithms end

@doc raw"""
"""
struct FHSlocal3 <: WeylPointsAlgorithms end

@doc raw"""
"""
struct Evar <: WeylPointsAlgorithms end