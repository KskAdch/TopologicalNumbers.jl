function psi!(i, psi1, Evec1, p) # wave function
    @unpack Hamiltonian, N = p

    k = 2pi * (i - 1) / N
    eigens = eigen!(Hamiltonian(k))
    psi1[:, :] .= eigens.vectors
    Evec1[:] .= eigens.values
end

@views function U!(Link, Evec0, Evec1, psi0, psi1, psiN1, i, p) # link variable
    @unpack N, gapless, Hs = p

    if i == 1
        psi!(i, psi1, Evec1, p)
        psiN1 .= psi1
    end

    if i != N
        psi0 .= psi1
        Evec0 .= Evec1
        psi!(i + 1, psi1, Evec1, p)
    end

    l = 1
    while l <= Hs
        l0 = Hs - count(Evec0 .> (gapless + Evec0[l]))

        if i == N
            if l == l0
                Link[l:l0] .= dot(psi1[:, l:l0], psiN1[:, l:l0])
            else
                Link[l:l0] .= det(psi1[:, l:l0]' * psiN1[:, l:l0])
            end
        else
            if l == l0
                Link[l:l0] .= dot(psi0[:, l:l0], psi1[:, l:l0])
            else
                Link[l:l0] .= det(psi0[:, l:l0]' * psi1[:, l:l0])
            end
        end

        l = 1 + l0
    end
end

@views function BerryPhase_round!(TopologicalNumber, p) # berry phase
    @unpack N, Hs = p
    Link = zeros(ComplexF64, Hs)

    Evec0 = zeros(Hs)
    Evec1 = zeros(Hs)

    psi0 = zeros(ComplexF64, Hs, Hs)
    psi1 = zeros(ComplexF64, Hs, Hs)
    psiN1 = zeros(ComplexF64, Hs, Hs)

    phi = zeros(Hs)

    TN = zeros(Hs)

    for i in 1:N
        U!(Link, Evec0, Evec1, psi0, psi1, psiN1, i, p)
        F!(phi, Link, p)
        TN[:] .+= phi[:]
    end

    TopologicalNumber .= [abs(rem(round(Int, TN[i] / pi), 2)) for i in 1:Hs]
end

@views function BerryPhase!(TopologicalNumber, p) # berry phase
    @unpack N, Hs = p
    Link = zeros(ComplexF64, Hs)

    Evec0 = zeros(Hs)
    Evec1 = zeros(Hs)

    psi0 = zeros(ComplexF64, Hs, Hs)
    psi1 = zeros(ComplexF64, Hs, Hs)
    psiN1 = zeros(ComplexF64, Hs, Hs)

    phi = zeros(Hs)

    TN = zeros(Hs)

    for i in 1:N
        U!(Link, Evec0, Evec1, psi0, psi1, psiN1, i, p)
        F!(phi, Link, p)
        TN[:] .+= phi[:]
    end

    for i in 1:Hs
        while abs(TN[i]) > 1.5pi
            TN[i] = abs(TN[i]) - 2pi
        end
    end
    TopologicalNumber .= [abs(rem(TN[i] / pi, 2)) for i in 1:Hs]
end

@views function F!(phi, Link, p) # lattice field strength
    @unpack N, Hs = p

    phi[:] .= [imag(log(Link[l])) for l in 1:Hs]
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
U_{n,i}(k)=\braket{\Psi_{n}(k)|\Psi_{n}(k+e_{1})}
```
 $\ket{\Psi_{n}(k)}$ is the wave function of the $n$th band.
"""
function calcBerryPhase(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

    Hs = size(Hamiltonian(0.0))[1]
    p = (; Hamiltonian, N, gapless, rounds, Hs, dim=1)

    if rounds == true
        TopologicalNumber = zeros(Int64, Hs)

        BerryPhase_round!(TopologicalNumber, p)

        Total = rem(sum(TopologicalNumber), 2)
    else
        TopologicalNumber = zeros(Hs)

        BerryPhase!(TopologicalNumber, p)

        Total = abs(sum(TopologicalNumber))
        while Total > 1.5
            Total -= 2
        end
        Total = abs(rem(Total, 2))
    end

    (; TopologicalNumber, Total)
end