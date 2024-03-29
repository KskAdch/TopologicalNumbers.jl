function psimat_square!(n, psimat, Evec, p::Params) # wave function □
    @unpack Ham, N = p

    n10 = n .+ [1, 0]
    n11 = n .+ [1, 1]
    n01 = n .+ [0, 1]

    n10 .= [mod(n10[i], N) for i in 1:2]
    n11 .= [mod(n11[i], N) for i in 1:2]
    n01 .= [mod(n01[i], N) for i in 1:2]

    k1 = 2pi * n / N
    k2 = 2pi * n10 / N
    k3 = 2pi * n11 / N
    k4 = 2pi * n01 / N

    eigens = eigen!(Ham(k1))
    psimat[1, :, :] .= eigens.vectors
    Evec[:] .= eigens.values

    psimat[2, :, :] .= eigen!(Ham(k2)).vectors
    psimat[3, :, :] .= eigen!(Ham(k3)).vectors
    return psimat[4, :, :] .= eigen!(Ham(k4)).vectors
end

@views function Linkmat_square!(psimat, Evec, Linkmat, p::Params)
    @unpack gapless, Hs = p

    l = 1
    while l <= Hs
        l0 = Hs - count(Evec .> (gapless + Evec[l]))

        if l == l0
            Linkmat[1, l:l0] .= dot(psimat[1, :, l:l0], psimat[2, :, l:l0])
            Linkmat[2, l:l0] .= dot(psimat[2, :, l:l0], psimat[3, :, l:l0])
            Linkmat[3, l:l0] .= dot(psimat[4, :, l:l0], psimat[3, :, l:l0])
            Linkmat[4, l:l0] .= dot(psimat[1, :, l:l0], psimat[4, :, l:l0])
        else
            Linkmat[1, l:l0] .= det(psimat[1, :, l:l0]' * psimat[2, :, l:l0])
            Linkmat[2, l:l0] .= det(psimat[2, :, l:l0]' * psimat[3, :, l:l0])
            Linkmat[3, l:l0] .= det(psimat[4, :, l:l0]' * psimat[3, :, l:l0])
            Linkmat[4, l:l0] .= det(psimat[1, :, l:l0]' * psimat[4, :, l:l0])
        end

        l = 1 + l0
    end
end

@views function F!(Linkmat, TopologicalNumber, p::Params) # lattice field strength
    @unpack rounds, Hs = p

    dphi = zeros(Hs)

    TopologicalNumber[:] = [
        angle(Linkmat[1, l] * Linkmat[2, l] * conj(Linkmat[3, l]) * conj(Linkmat[4, l])) for
        l in 1:Hs
    ]
    dphi[:] = [
        angle(Linkmat[1, l]) + angle(Linkmat[2, l]) - angle(Linkmat[3, l]) -
        angle(Linkmat[4, l]) for l in 1:Hs
    ]

    return TopologicalNumber .= (TopologicalNumber - dphi) ./ 2pi
end

@doc raw"""

 Calculate the Berry flux in the two-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    calcBerryFlux(Hamiltonian::Function, n::Vector{Int64}; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

 Arguments
 - Hamiltionian::Function: the Hamiltonian matrix with one-dimensional wavenumber `k` as an argument.
 - n::Vector{Int64}: The wavenumber($2\pi n/N$) when calculating Berry flux.
 - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.


# Definition
 The Berry flux at the wavenumber $\bm{k}$ of the $n$th band $F_{n}(\bm{k})$ is defined by
```math
F_{n}(\bm{k})=f_{n}(\bm{k})-df_{n}(\bm{k})
```
```math
f_{n}(\bm{k})=\frac{1}{2\pi}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,1}^{*}(\bm{k})\right]\right]
```
```math
df_{n}(\bm{k})=\frac{1}{2\pi}\mathrm{Im}\left[\mathrm{Log}[U_{n,1}(\bm{k})]+\mathrm{Log}[U_{n,2}(\bm{k}+\bm{e}_{1})]-\mathrm{Log}[U_{n,1}(\bm{k}+\bm{e}_{2})]-\mathrm{Log}[U_{n,1}(\bm{k})]\right]
```
 $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $n$th band.
"""
function calcBerryFlux(
    Hamiltonian::Function, n::Vector{Int64}; N::Int=51, gapless::Real=0.0, rounds::Bool=true
)
    Hs = size(Hamiltonian(n), 1)
    p = Params(; Ham=Hamiltonian, N, gapless, rounds, Hs, dim=2)

    psimat = zeros(ComplexF64, 4, Hs, Hs)
    Evec = zeros(Hs)
    Linkmat = zeros(ComplexF64, 4, Hs)
    TopologicalNumber = zeros(Hs)

    n .= [mod(n[i], N) for i in 1:2]

    psimat_square!(n, psimat, Evec, p)
    Linkmat_square!(psimat, Evec, Linkmat, p)
    F!(Linkmat, TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
    end

    return (; TopologicalNumber, n)
end

@doc raw"""
Calculate the Berry flux in the two-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    solve(prob::LBFProblem, alg::T1=FHSlocal2(); parallel::T2=UseSingleThread()) where {T1<:BerryFluxAlgorithms,T2<:TopologicalNumbersParallel}

# Arguments
- `prob::LBFProblem`: The LBFProblem struct that contains the Hamiltonian matrix function in the wave number space and other parameters.
- `alg::T1=FHSlocal2()`: The algorithm to use for calculating the second Chern numbers. Default is `FHSlocal2` algorithm.
- `parallel::T2=UseSingleThread()`: The parallelization strategy to use. Default is to use a single thread.

# Returns
- `LBFSolution`: A struct that contains the calculated Berry flux.

# Examples

```julia
julia> H(k) = Flux2d(k, (6, 1))
julia> n = [0, 0]
julia> prob = LBFProblem(H, n)
julia> result = solve(prob)
LBFSolution{Vector{Int64}, Vector{Int64}}([0, 0, 0, 0, -1, 0], [0, 0])
julia> result.TopologicalNumber
6-element Vector{Int64}:
  0
  0
  0
  0
 -1
  0
```

"""
function solve(
    prob::LBFProblem, alg::T1=FHSlocal2(); parallel::T2=UseSingleThread()
) where {T1<:BerryFluxAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack H, n, N, gapless, rounds = prob

    Hs = size(H(n), 1)
    p = Params(; Ham=H, N, gapless, rounds, Hs, dim=2)

    psimat = zeros(ComplexF64, 4, Hs, Hs)
    Evec = zeros(Hs)
    Linkmat = zeros(ComplexF64, 4, Hs)
    TopologicalNumber = zeros(Hs)

    n .= [mod(n[i], N) for i in 1:2]

    psimat_square!(n, psimat, Evec, p)
    Linkmat_square!(psimat, Evec, Linkmat, p)
    F!(Linkmat, TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
    end

    return LBFSolution(; TopologicalNumber, n)
end
