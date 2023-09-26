function psi!(n, psimat, Evec, p) # wave function
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

    eigens1 = eigen!(Hamiltonian(k1))
    psimat[1, :, :] .= eigens1.vectors
    Evec[:] .= eigens1.values

    psimat[2, :, :] .= eigen!(Hamiltonian(k2)).vectors
    psimat[3, :, :] .= eigen!(Hamiltonian(k3)).vectors
    psimat[4, :, :] .= eigen!(Hamiltonian(k4)).vectors
end

@views function Link!(psimat, Evec, Linkmat, p)
    @unpack gapless, Hs = p

    l = 1
    while l <= Hs
        l0 = Hs - count(Evec .> (gapless + Evec[l]))

        if l == l0
            for i in 1:3
                Linkmat[i, l:l0] .= dot(psimat[i, :, l:l0], psimat[i+1, :, l:l0])
            end
            Linkmat[4, l:l0] .= dot(psimat[4, :, l:l0], psimat[1, :, l:l0])
        else
            for i in 1:3
                Linkmat[i, l:l0] .= det(psimat[i, :, l:l0]' * psimat[i+1, :, l:l0])
            end
            Linkmat[4, l:l0] .= det(psimat[4, :, l:l0]' * psimat[1, :, l:l0])
        end

        l = 1 + l0
    end
end

@views function F(Linkmat, p) # lattice field strength
    @unpack rounds, Hs = p

    phi = zeros(Hs)
    dphi = zeros(Hs)

    phi[:] = [imag(log(Linkmat[1, l] * Linkmat[2, l] * Linkmat[3, l] * Linkmat[4, l])) for l in 1:Hs]
    dphi[:] = [imag(log(Linkmat[1, l]) + log(Linkmat[2, l]) + log(Linkmat[3, l]) + log(Linkmat[4, l])) for l in 1:Hs]

    if rounds == true
        phi[:] = [round(Int, (phi[i] - dphi[i]) / 2pi) for i in 1:Hs]
    else
        phi .= (phi - dphi) / 2pi
    end
end

@doc raw"""

 Calculate the BerryFlux in the two-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    calcBerryFlux(Hamiltonian::Function, n::Vector{Int64}; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

 Arguments
 - Hamiltionian::Function: the Hamiltonian matrix with one-dimensional wavenumber `k` as an argument.
 - n::Vector{Int64}: The wavenumber when calculating Berry flux.
 - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.


# Definition
 The Berry flux at the wavenumber $\bm{k}$ of the $n$th band $F_{n}(\bm{k})$ is defined by
```math
F_{n}(\bm{k})=\frac{1}{2\pi i}\left(\partial_{k_{1}}A_{n,2}(\bm{k})-\partial_{k_{2}}A_{n,1}(\bm{k})\right)
```
 $A_{n,i}(\bm{k})$ is the Berry connection at wavenumber $\bm{k}$.
```math
A_{n,i}(\bm{k})=\bra{\Psi_{n}(\bm{k})}\partial_{k_{i}}\ket{\Psi_{n}(\bm{k})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $n$th band.
"""
function calcBerryFlux(Hamiltonian::Function, n::Vector{Int64}; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

    Hs = size(Hamiltonian(n))[1]
    p = (; Hamiltonian, N, gapless, rounds, Hs)

    psimat = zeros(ComplexF64, 4, Hs, Hs)
    Evec = zeros(Hs)
    Linkmat = zeros(ComplexF64,4, Hs)

    if n[1] >= N
        n[1] = mod(n[1], N)
    end
    
    if n[2] >= N
        n[2] = mod(n[2], N)
    end

    psi!(n, psimat, Evec, p)
    Link!(psimat, Evec, Linkmat, p)
    F(Linkmat, p)
end