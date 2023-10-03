function psimat!(n, psimat, Evec, p) # wave function □
    @unpack Hamiltonian, N = p

    if n[1] == N-1 && n[2] == N-1
        n10 = [0, n[2]]
        n11 = [0, 0]
        n01 = [n[1], 0]
    elseif n[1] == N-1
        n10 = [0, n[2]]
        n11 = [0, n[2]+1]
        n01 = [n[1], n[2]+1]
    elseif n[2] == N-1
        n10 = [n[1]+1, n[2]]
        n11 = [n[1]+1, 0]
        n01 = [n[1], 0]
    else
        n10 = n .+ [1, 0]
        n11 = n .+ [1, 1]
        n01 = n .+ [0, 1]
    end

    k1 = 2pi * n / N
    k2 = 2pi * n10 / N
    k3 = 2pi * n11 / N
    k4 = 2pi * n01 / N

    eigens = eigen!(Hamiltonian(k1))
    psimat[1, :, :] .= eigens.vectors
    Evec[:] .= eigens.values

    psimat[2, :, :] .= eigen!(Hamiltonian(k2)).vectors
    psimat[3, :, :] .= eigen!(Hamiltonian(k3)).vectors
    psimat[4, :, :] .= eigen!(Hamiltonian(k4)).vectors
end

@views function Linkmat!(psimat, Evec, Linkmat, p)
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
            Linkmat[1, l:l0] .= det(psimat[1, :, l:l0]' * psimat[1, :, l:l0])
            Linkmat[2, l:l0] .= det(psimat[2, :, l:l0]' * psimat[2, :, l:l0])
            Linkmat[3, l:l0] .= det(psimat[4, :, l:l0]' * psimat[4, :, l:l0])
            Linkmat[4, l:l0] .= det(psimat[1, :, l:l0]' * psimat[1, :, l:l0])
        end

        l = 1 + l0
    end
end

@views function F(Linkmat, p) # lattice field strength
    @unpack rounds, Hs = p

    phi = zeros(Hs)
    dphi = zeros(Hs)

    phi[:] = [imag(log(Linkmat[1, l] * Linkmat[2, l] * conj(Linkmat[3, l]) * conj(Linkmat[4, l]))) for l in 1:Hs]
    dphi[:] = [imag(log(Linkmat[1, l]) + log(Linkmat[2, l]) - log(Linkmat[3, l]) - log(Linkmat[4, l])) for l in 1:Hs]

    if rounds == true
        phi[:] = [round(Int, (phi[i] - dphi[i]) / 2pi) for i in 1:Hs]
    else
        phi .= (phi - dphi) / 2pi
    end
end

@doc raw"""

 Calculate the Berry flux in the two-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    calcBerryFlux(Hamiltonian::Function, n::Vector{Int64}; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

 Arguments
 - Hamiltionian::Function: the Hamiltonian matrix with one-dimensional wavenumber `k` as an argument.
 - n::Vector{Int64}: The wavenumber($2\pi*n/N$) when calculating Berry flux.
 - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.


# Definition
 The Berry flux at the wavenumber $\bm{k}$ of the $n$th band $F_{n}(\bm{k})$ is defined by
```math
F_{n}(\bm{k})=f_{n}(\bm{k})-df_{n}(\bm{k})
```
```math
f_{n}(\bm{k})=\frac{1}{2\pi}\Im\log\left(U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,1}^{*}(\bm{k})\right)
```
```math
df_{n}(\bm{k})=\frac{1}{2\pi}\Im\left(\logU_{n,1}(\bm{k})+\logU_{n,2}(\bm{k}+\bm{e}_{1})-\logU_{n,1}(\bm{k}+\bm{e}_{2})-\logU_{n,1}(\bm{k})\right)
```
 $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is primitive translation vector.
```math
U_{n,i}(\bm{k})=\bra{\Psi_{n}(\bm{k})}\ket{\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $n$th band.
"""
function calcBerryFlux(Hamiltonian::Function, n::Vector{Int64}; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

    Hs = size(Hamiltonian(n))[1]
    p = Params(; Hamiltonian, N, gapless, rounds, Hs, dim=2)

    psimat = zeros(ComplexF64, 4, Hs, Hs)
    Evec = zeros(Hs)
    Linkmat = zeros(ComplexF64,4, Hs)

    n[1] = mod(n[1], N)
    n[2] = mod(n[2], N)

    psimat!(n, psimat, Evec, p)
    Linkmat!(psimat, Evec, Linkmat, p)
    F(Linkmat, p)
end