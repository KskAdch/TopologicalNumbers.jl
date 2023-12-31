function psi_j!(j, psi_1, Evec1, p::Params) # wave function
    @unpack Hamiltonian, N = p
    for i in 1:N
        k = [i-1, j-1] * 2pi / N .+ 2pi * [1e-5, 1e-5]
        eigens = eigen!(Hamiltonian(k))
        psi_1[i, :, :] .= eigens.vectors
        Evec1[i, :] .= eigens.values
    end
end

@views function Link!(psi00, psi01, psi10, Enevec, link10, link01, gapless, Hs)
    l = 1
    while l <= Hs
        l0 = Hs - count(Enevec .> (gapless + Enevec[l]))

        if l == l0
            link10[l:l0] .= dot(psi00[:, l:l0], psi10[:, l:l0])
            link01[l:l0] .= dot(psi00[:, l:l0], psi01[:, l:l0])
        else
            link10[l:l0] .= det(psi00[:, l:l0]' * psi10[:, l:l0])
            link01[l:l0] .= det(psi00[:, l:l0]' * psi01[:, l:l0])
        end

        l = 1 + l0
    end
end

@views function U!(Link0, Link1, LinkN, link10, link01, psi_0, psi_1, psi_N, Evec0, Evec1, psi00, psi10, psi01, Enevec, j, p::Params) # link variable
    @unpack N, gapless, Hs = p

    if j != 1
        Link0 .= Link1
    end

    if j != N
        if j == 1

            psi_j!(1, psi_1, Evec1, p)
            psi_N .= psi_1
            psi_0 .= psi_1
            Evec0 .= Evec1
            psi_j!(2, psi_1, Evec1, p)

            for i in 1:N

                psi00[:, :] = psi_0[i, :, :]
                if i == N
                    psi10[:, :] = psi_0[1, :, :]
                    psi01[:, :] = psi_1[N, :, :]
                else
                    psi10[:, :] = psi_0[i+1, :, :]
                    psi01[:, :] = psi_1[i, :, :]
                end

                Enevec[:] = Evec1[i, :]
                Link!(psi00, psi01, psi10, Enevec, link10, link01, gapless, Hs)

                Link0[:, :, i] .= [link10 link01]
                LinkN .= Link0
            end
        end

        psi_0 .= psi_1
        Evec0 .= Evec1

        if j != N - 1
            psi_j!(j + 2, psi_1, Evec1, p)
        end

        for i in 1:N

            psi00[:, :] = psi_0[i, :, :]
            if i == N && j == N - 1
                psi10[:, :] = psi_0[1, :, :]
                psi01[:, :] = psi_N[N, :, :]
            elseif i == N
                psi10[:, :] = psi_0[1, :, :]
                psi01[:, :] = psi_1[N, :, :]
            elseif j == N - 1
                psi10[:, :] = psi_0[i+1, :, :]
                psi01[:, :] = psi_N[i, :, :]
            else
                psi10[:, :] = psi_0[i+1, :, :]
                psi01[:, :] = psi_1[i, :, :]
            end

            Enevec[:] = Evec0[i, :]
            Link!(psi00, psi01, psi10, Enevec, link10, link01, gapless, Hs)

            Link1[:, :, i] .= [link10 link01]
        end
    end
end

@views function F!(phi, dphi, i, j, Link0, Link1, LinkN, p::Params) # lattice field strength
    @unpack N, rounds, Hs = p

    if i == N && j == N
        phi[:] = [angle(Link0[l, 1, N] * Link0[l, 2, 1] * conj(LinkN[l, 1, N]) * conj(Link0[l, 2, N])) for l in 1:Hs]
        dphi[:] = [angle(Link0[l, 1, N]) + angle(Link0[l, 2, 1]) - angle(LinkN[l, 1, N]) - angle(Link0[l, 2, N]) for l in 1:Hs]
    elseif i == N
        phi[:] = [angle(Link0[l, 1, N] * Link0[l, 2, 1] * conj(Link1[l, 1, N]) * conj(Link0[l, 2, N])) for l in 1:Hs]
        dphi[:] = [angle(Link0[l, 1, N]) + angle(Link0[l, 2, 1]) - angle(Link1[l, 1, N]) - angle(Link0[l, 2, N]) for l in 1:Hs]
    elseif j == N
        phi[:] = [angle(Link0[l, 1, i] * Link0[l, 2, i+1] * conj(LinkN[l, 1, i]) * conj(Link0[l, 2, i])) for l in 1:Hs]
        dphi[:] = [angle(Link0[l, 1, i]) + angle(Link0[l, 2, i+1]) - angle(LinkN[l, 1, i]) - angle(Link0[l, 2, i]) for l in 1:Hs]
    else
        phi[:] = [angle(Link0[l, 1, i] * Link0[l, 2, i+1] * conj(Link1[l, 1, i]) * conj(Link0[l, 2, i])) for l in 1:Hs]
        dphi[:] = [angle(Link0[l, 1, i]) + angle(Link0[l, 2, i+1]) - angle(Link1[l, 1, i]) - angle(Link0[l, 2, i]) for l in 1:Hs]
    end

    phi .= (phi - dphi) ./ 2pi
end

@views function ChernPhase!(TopologicalNumber, p::Params) # chern number # Bug
    @unpack N, Hs = p
    TopologicalNumber[:] .= zero(Float64)
    Link0 = zeros(ComplexF64, Hs, 2, N)
    Link1 = zeros(ComplexF64, Hs, 2, N)
    LinkN = zeros(ComplexF64, Hs, 2, N)
    link1 = zeros(ComplexF64, Hs)
    link2 = zeros(ComplexF64, Hs)

    psi_0 = zeros(ComplexF64, N, Hs, Hs)
    psi_1 = zeros(ComplexF64, N, Hs, Hs)
    psi_N = zeros(ComplexF64, N, Hs, Hs)
    Evec0 = zeros(N, Hs)
    Evec1 = zeros(N, Hs)

    psi00 = zeros(ComplexF64, Hs, Hs)
    psi10 = zeros(ComplexF64, Hs, Hs)
    psi01 = zeros(ComplexF64, Hs, Hs)
    Enevec = zeros(Hs)

    phi = zeros(Hs)
    dphi = zeros(Hs)

    for j in 1:N
        U!(Link0, Link1, LinkN, link1, link2, psi_0, psi_1, psi_N, Evec0, Evec1, psi00, psi10, psi01, Enevec, j, p)
        for i in 1:N
            F!(phi, dphi, i, j, Link0, Link1, LinkN, p)
            TopologicalNumber[:] .+= phi[:]
        end
    end
end

@doc raw"""

 Calculate the first Chern numbers in the two-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    calcChern(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

 Arguments
 - Hamiltionian::Function: The Hamiltonian matrix with two-dimensional wavenumber `k` as an argument.
 - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.


# Definition
 The first Chern number of the $n$th band $\nu_{n}$ is defined by
```math
\nu_{n}=\frac{1}{2\pi}\sum_{\bm{k}\in\mathrm{BZ}}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,2}^{*}(\bm{k})\right]\right]
```
 The range $\mathrm{BZ}$(Brillouin Zone) is $\bm{k}\in[0,2\pi]^{2}$. $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $n$th band.
"""
function calcChern(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

    Hs = size(Hamiltonian(zeros(2)), 1)
    p = Params(; Hamiltonian, N, gapless, rounds, Hs, dim=2)

    TopologicalNumber = zeros(Hs)

    ChernPhase!(TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
    end

    Total = sum(TopologicalNumber)

    (; TopologicalNumber, Total)
end
