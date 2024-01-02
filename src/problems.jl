abstract type TopologicalNumbersProblems end

# Problem for calculating the first Chern number
@kwdef struct BPProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Real,T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 51
    gapless::T3 = 0.0
    rounds::T4 = true
end

# Problem for calculating the first Chern number
@kwdef struct FCProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Real,T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 51
    gapless::T3 = 0.0
    rounds::T4 = true
end

# Problem for calculating the second Chern number
@kwdef struct SCProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Union{Real,Nothing},T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 30
    Nfill::T3 = nothing
    RV::T4 = true
end

# Problem for calculating the Z2 invariant
@kwdef struct Z2Problem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 50
    rounds::T3 = true
    TR::T3 = false
end



abstract type TopologicalNumbersSolutions end

# Solution for calculating the first Chern number
@kwdef struct BPSolution{T} <: TopologicalNumbersSolutions
    TopologicalNumber::T = nothing
    Total::T = nothing
end

# Solution for calculating the first Chern number
@kwdef struct FCSolution{T} <: TopologicalNumbersSolutions
    TopologicalNumber::T = nothing
    Total::T = nothing
end

# Solution for calculating the second Chern number
@kwdef struct SCSolution{T} <: TopologicalNumbersSolutions
    TopologicalNumber::T = nothing
end

# Solution for calculating the Z2 invariant
@kwdef struct Z2Solution{T} <: TopologicalNumbersSolutions
    TopologicalNumber::T = nothing
    TRTopologicalNumber::T = nothing
    Total::T = nothing
end
