function psi!(i, v::TemporalBerryPhase, p::Params) # wave function
    @unpack Ham, N = p

    k = 2pi * (i - 1) / N .+ 2pi * 1e-5
    eigens = eigen!(Ham(k))
    v.psi1[:, :] .= eigens.vectors
    v.Evec1[:] .= eigens.values
end

@views function U!(i, v::TemporalBerryPhase, p::Params) # link variable
    @unpack N, gapless, Hs = p

    if i == 1
        psi!(i, v, p)
        v.psiN1 .= v.psi1
    end

    if i != N
        v.psi0 .= v.psi1
        v.Evec0 .= v.Evec1
        psi!(i + 1, v, p)
    end

    l = 1
    while l <= Hs
        l0 = Hs - count(v.Evec0 .> (gapless + v.Evec0[l]))

        if i == N
            if l == l0
                v.Link[l:l0] .= dot(v.psi1[:, l:l0], v.psiN1[:, l:l0])
            else
                v.Link[l:l0] .= det(v.psi1[:, l:l0]' * v.psiN1[:, l:l0])
            end
        else
            if l == l0
                v.Link[l:l0] .= dot(v.psi0[:, l:l0], v.psi1[:, l:l0])
            else
                v.Link[l:l0] .= det(v.psi0[:, l:l0]' * v.psi1[:, l:l0])
            end
        end

        l = 1 + l0
    end
end

# @views function BerryPhase_round!(TopologicalNumber, p::Params) # berry phase
#     @unpack N, Hs = p
#     Link = zeros(ComplexF64, Hs)

#     Evec0 = zeros(Hs)
#     Evec1 = zeros(Hs)

#     psi0 = zeros(ComplexF64, Hs, Hs)
#     psi1 = zeros(ComplexF64, Hs, Hs)
#     psiN1 = zeros(ComplexF64, Hs, Hs)

#     phi = zeros(Hs)

#     TN = zeros(Hs)

#     for i in 1:N
#         U!(Link, Evec0, Evec1, psi0, psi1, psiN1, i, p)
#         L!(phi, Link, p)
#         TN[:] .+= phi[:]
#     end

#     TopologicalNumber .= [abs(rem(round(Int, TN[i] / pi), 2)) for i in 1:Hs]
# end

@views function L!(phi, v::TemporalBerryPhase, p::Params) # lattice field strength
    @unpack N, Hs = p

    phi .= angle.(v.Link)
end

@doc raw"""
"""
@views function BerryPhase!(TopologicalNumber, p::Params) # berry phase
    @unpack N, Hs = p
    Link = zeros(ComplexF64, Hs)

    Evec0 = zeros(Hs)
    Evec1 = zeros(Hs)

    psi0 = zeros(ComplexF64, Hs, Hs)
    psi1 = zeros(ComplexF64, Hs, Hs)
    psiN1 = zeros(ComplexF64, Hs, Hs)

    v = TemporalBerryPhase(Link, Evec0, Evec1, psi0, psi1, psiN1)

    phi = zeros(Hs)

    TN = zeros(Hs)

    for i in 1:N
        U!(i, v, p)
        L!(phi, v, p)
        TN[:] .+= phi[:]
    end

    for i in 1:Hs
        TopologicalNumber[i] = 1 - abs(1 - rem(abs(TN[i]) / pi, 2))
    end
end

@doc raw"""

 Calculate the winding numbers in the one-dimensional case.

    calcBerryPhase(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)


 Arguments
 - `Hamiltonian::Function`: the Hamiltonian matrix function with one-dimensional wavenumber `k` as an argument.
 - `N::Int`: the number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - `gapless::Real`: the threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - `rounds::Bool`: an option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.


# Definition

 The Berry phase of the $n$th band $\nu_{n}$ is defined by
```math
\nu_{n}=\frac{1}{\pi}\sum_{k\in\mathrm{BZ}}U_{n}(k)
```
 The range $\mathrm{BZ}$(Brillouin Zone) is $k\in[0,2\pi]$. $U_{n,i}(k)$ is the link variable at wavenumber $k$. $e_{1}$ is the unit vector.
```math
U_{n}(k)=\braket{\Psi_{n}(k)|\Psi_{n}(k+e_{1})}
```
 $\ket{\Psi_{n}(k)}$ is the wave function of the $n$th band.
"""
function calcBerryPhase(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

    Hs = size(Hamiltonian(0.0), 1)
    p = Params(; Ham=Hamiltonian, N, gapless, rounds, Hs, dim=1)

    TopologicalNumber = zeros(Hs)
    BerryPhase!(TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
        Total = rem(sum(TopologicalNumber), 2)
    else
        Total = abs(sum(TopologicalNumber))
        Total = rem(1 - abs(1 - Total), 2)
    end

    (; TopologicalNumber, Total)
end


@doc raw"""

 Calculate the winding numbers in the one-dimensional case.

    solve(prob::BPProblem, alg::T1=BP(); parallel::T2=UseSingleThread()) where {T1<:BerryPhaseAlgorithms,T2<:TopologicalNumbersParallel}


 Arguments
 - `Hamiltonian::Function`: the Hamiltonian matrix function with one-dimensional wavenumber `k` as an argument.
 - `N::Int`: the number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - `gapless::Real`: the threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - `rounds::Bool`: an option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.


# Definition

 The Berry phase of the $n$th band $\nu_{n}$ is defined by
```math
\nu_{n}=\frac{1}{\pi}\sum_{k\in\mathrm{BZ}}U_{n}(k)
```
 The range $\mathrm{BZ}$(Brillouin Zone) is $k\in[0,2\pi]$. $U_{n,i}(k)$ is the link variable at wavenumber $k$. $e_{1}$ is the unit vector.
```math
U_{n}(k)=\braket{\Psi_{n}(k)|\Psi_{n}(k+e_{1})}
```
 $\ket{\Psi_{n}(k)}$ is the wave function of the $n$th band.
"""
function solve(
    prob::BPProblem,
    alg::T1=BP();
    parallel::T2=UseSingleThread()
) where {T1<:BerryPhaseAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack H, N, gapless, rounds = prob

    Hs = size(H(0.0), 1)
    p = Params(; Ham=H, N, gapless, rounds, Hs, dim=1)

    TopologicalNumber = zeros(Hs)
    BerryPhase!(TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
        Total = rem(sum(TopologicalNumber), 2)
    else
        Total = abs(sum(TopologicalNumber))
        Total = rem(1 - abs(1 - Total), 2)
    end

    BPSolution(; TopologicalNumber, Total)
end
