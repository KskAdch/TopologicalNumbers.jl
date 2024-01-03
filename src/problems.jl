abstract type TopologicalNumbersProblems end

# Problem for calculating the Berry phase
@kwdef struct BPProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Real,T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 51
    gapless::T3 = 0.0
    rounds::T4 = true
end
# default
BPProblem(H) = BPProblem(; H=H)
BPProblem(H, N) = BPProblem(; H=H, N=N)

# Problem for calculating the first Chern number
@kwdef struct FCProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Real,T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 51
    gapless::T3 = 0.0
    rounds::T4 = true
end
# default
FCProblem(H) = FCProblem(; H=H)
FCProblem(H, N) = FCProblem(; H=H, N=N)

# Problem for calculating the second Chern number
@kwdef struct SCProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Union{Real,Nothing},T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 30
    Nfill::T3 = nothing
    RV::T4 = true
end
# default
SCProblem(H) = SCProblem(; H=H)
SCProblem(H, N) = SCProblem(; H=H, N=N)

# Problem for calculating the Z2 invariant
@kwdef struct Z2Problem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 50
    rounds::T3 = true
    TR::T3 = false
end
# default
Z2Problem(H) = Z2Problem(; H=H)
Z2Problem(H, N) = Z2Problem(; H=H, N=N)



abstract type TopologicalNumbersSolutions end

# Solution for calculating the first Chern number
@kwdef struct BPSolution{T1,T2} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    Total::T2 = nothing
end

# Solution for calculating the first Chern number
@kwdef struct FCSolution{T1,T2} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    Total::T2 = nothing
end

# Solution for calculating the second Chern number
@kwdef struct SCSolution{T} <: TopologicalNumbersSolutions
    TopologicalNumber::T = nothing
end

# Solution for calculating the Z2 invariant
@kwdef struct Z2Solution{T1,T2,T3} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    TRTopologicalNumber::T2 = nothing
    Total::T3 = nothing
end
