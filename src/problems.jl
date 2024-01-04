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

# Problem for calculating the local Berry flux
@kwdef struct LBFProblem{T1<:Function,T2<:AbstractVector,T3<:Union{Tuple,AbstractVector,Int},T4<:Real,T5<:Bool} <: TopologicalNumbersProblems
    H::T1
    n::T2
    N::T3 = 51
    gapless::T4 = 0.0
    rounds::T5 = true
end
# default
LBFProblem(H, n) = LBFProblem(; H=H, n=n)
LBFProblem(H, n, N) = LBFProblem(; H=H, n=n, N=N)

# Problem for finding and calculating the Weyl points
@kwdef struct WCSProblem{T1<:Function,T2<:String,T3<:Int,T4<:Union{Tuple,AbstractVector,Int},T5<:Real,T6<:Bool} <: TopologicalNumbersProblems
    H::T1
    kn::T2
    kn_mesh::T3 = 51
    N::T4 = 51
    gapless::T5 = 0.0
    rounds::T6 = true
end
# default
WCSProblem(H, kn) = WCSProblem(; H=H, kn=kn)
WCSProblem(H, kn, N::T) where {T<:Int} = WCSProblem(; H=H, kn=kn, kn_mesh=N, N=N)
WCSProblem(H, kn, N1, N2) = WCSProblem(; H=H, kn=kn, kn_mesh=N1, N=N2)

# Problem for finding and calculating the Weyl points
@kwdef struct WNProblem{T1<:Function,T2<:AbstractVector,T3<:Int,T4<:Real,T5<:Bool} <: TopologicalNumbersProblems
    H::T1
    n::T2
    N::T3 = 51
    gapless::T4 = 0.0
    rounds::T5 = true
end
# default
WNProblem(H, n) = WNProblem(; H=H, n=n)
WNProblem(H, n, N) = WNProblem(; H=H, n=n, N=N)

# Problem for finding and calculating the Weyl points
@kwdef struct WPProblem{T1<:Function,T2<:Int,T3<:AbstractVector,T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 10
    gapless::T3 = [1e-1, 1e-2, 1e-3, 1e-4]
    rounds::T4 = true
end
# default
WPProblem(H) = WPProblem(; H=H)
WPProblem(H, N) = WPProblem(; H=H, N=N)



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

# Solution for calculating the local Berry flux
@kwdef struct LBFSolution{T1,T2} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    n::T2 = nothing
end

# Solution for finding and calculating the Weyl points
@kwdef struct WCSSolution{T1,T2,T3} <: TopologicalNumbersSolutions
    kn::T1 = nothing
    param::T2 = nothing
    nums::T3 = nothing
end

# Solution for finding and calculating the Weyl points
@kwdef struct WNSolution{T1,T2,T3} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    n::T2 = nothing
    N::T3 = nothing
end

# Solution for finding and calculating the Weyl points
@kwdef struct WPSolution{T1,T2,T3} <: TopologicalNumbersSolutions
    WeylPoint::T1 = nothing
    N::T2 = nothing
    Nodes::T3 = nothing
end
