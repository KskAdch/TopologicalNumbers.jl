abstract type TopologicalNumbersProblems end

# Problem for calculating the first Chern number
@kwdef struct BPProblem <: TopologicalNumbersProblems
    H
    N::T2 = 51
    gapless = 0.0
    rounds = true
end

# Problem for calculating the first Chern number
@kwdef struct FCProblem <: TopologicalNumbersProblems
    H
    N::T2 = 51
    gapless = 0.0
    rounds = true
end

# Problem for calculating the second Chern number
struct SCProblem <: TopologicalNumbersProblems
    H
    N::T2 = 30
    Nfill = nothing
    RV = true
end

# Problem for calculating the Z2 invariant
struct Z2Problem <: TopologicalNumbersProblems
    H
    N::T2 = 50
    rounds = true
    TR = false
end



abstract type TopologicalNumbersSolutions end

# Solution for calculating the first Chern number
@kwdef struct BPSolution <: TopologicalNumbersSolutions
    TopologicalNumber = nothing
    Total = nothing
end

# Solution for calculating the first Chern number
@kwdef struct FCSolution <: TopologicalNumbersSolutions
    TopologicalNumber = nothing
    Total = nothing
end

# Solution for calculating the second Chern number
struct SCSolution <: TopologicalNumbersSolutions
    TopologicalNumber = nothing
end

# Solution for calculating the Z2 invariant
struct Z2Solution <: TopologicalNumbersSolutions
    TopologicalNumber = nothing
    TRTopologicalNumber = nothing
    Total = nothing
end
