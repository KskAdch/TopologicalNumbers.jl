function psi_j!(j, v::TemporalFirstChern, p::Params) # wave function
    @unpack Ham, N = p
    for i in 1:N
        v.k[1] = (i - 1) * 2pi / N + 2pi * 1e-5
        v.k[2] = (j - 1) * 2pi / N + 2pi * 1e-5
        eigens = eigen!(Ham(v.k))
        v.psi_1[i, :, :] .= eigens.vectors
        v.Evec1[i, :] .= eigens.values
    end
end

@views function Link!(v::TemporalFirstChern, gapless, Hs)
    l = 1
    while l <= Hs
        l0 = Hs - count(v.Enevec .> (gapless + v.Enevec[l]))

        if l == l0
            v.link10[l:l0] .= dot(v.psi00[:, l:l0], v.psi10[:, l:l0])
            v.link01[l:l0] .= dot(v.psi00[:, l:l0], v.psi01[:, l:l0])
        else
            v.link10[l:l0] .= det(v.psi00[:, l:l0]' * v.psi10[:, l:l0])
            v.link01[l:l0] .= det(v.psi00[:, l:l0]' * v.psi01[:, l:l0])
        end

        l = 1 + l0
    end
end

@views function U!(j, v::TemporalFirstChern, p::Params) # link variable
    @unpack N, gapless, Hs = p

    if j != 1
        v.Link0 .= v.Link1
    end

    if j != N
        if j == 1

            psi_j!(1, v, p)
            v.psi_N .= v.psi_1
            v.psi_0 .= v.psi_1
            v.Evec0 .= v.Evec1
            psi_j!(2, v, p)

            for i in 1:N

                v.psi00[:, :] = v.psi_0[i, :, :]
                if i == N
                    v.psi10[:, :] = v.psi_0[1, :, :]
                    v.psi01[:, :] = v.psi_1[N, :, :]
                else
                    v.psi10[:, :] = v.psi_0[i+1, :, :]
                    v.psi01[:, :] = v.psi_1[i, :, :]
                end

                v.Enevec[:] = v.Evec1[i, :]
                Link!(v, gapless, Hs)

                v.Link0[:, :, i] .= [v.link10 v.link01]
                v.LinkN .= v.Link0
            end
        end

        v.psi_0 .= v.psi_1
        v.Evec0 .= v.Evec1

        if j != N - 1
            psi_j!(j + 2, v, p)
        end

        for i in 1:N

            v.psi00[:, :] = v.psi_0[i, :, :]
            if i == N && j == N - 1
                v.psi10[:, :] = v.psi_0[1, :, :]
                v.psi01[:, :] = v.psi_N[N, :, :]
            elseif i == N
                v.psi10[:, :] = v.psi_0[1, :, :]
                v.psi01[:, :] = v.psi_1[N, :, :]
            elseif j == N - 1
                v.psi10[:, :] = v.psi_0[i+1, :, :]
                v.psi01[:, :] = v.psi_N[i, :, :]
            else
                v.psi10[:, :] = v.psi_0[i+1, :, :]
                v.psi01[:, :] = v.psi_1[i, :, :]
            end

            v.Enevec[:] = v.Evec0[i, :]
            Link!(v, gapless, Hs)

            v.Link1[:, 1, i] .= v.link10
            v.Link1[:, 2, i] .= v.link01
        end
    end
end

@views function F!(phi, dphi, i, j, v::TemporalFirstChern, p::Params) # lattice field strength
    @unpack N, rounds, Hs = p

    if i != N && j != N
        for l in 1:Hs
            phi[l] = angle(v.Link0[l, 1, i] * v.Link0[l, 2, i+1] * conj(v.Link1[l, 1, i]) * conj(v.Link0[l, 2, i]))
            dphi[l] = angle(v.Link0[l, 1, i]) + angle(v.Link0[l, 2, i+1]) - angle(v.Link1[l, 1, i]) - angle(v.Link0[l, 2, i])
        end
    elseif i == N && j != N
        for l in 1:Hs
            phi[l] = angle(v.Link0[l, 1, N] * v.Link0[l, 2, 1] * conj(v.Link1[l, 1, N]) * conj(v.Link0[l, 2, N]))
            dphi[l] = angle(v.Link0[l, 1, N]) + angle(v.Link0[l, 2, 1]) - angle(v.Link1[l, 1, N]) - angle(v.Link0[l, 2, N])
        end
    elseif i != N && j == N
        for l in 1:Hs
            phi[l] = angle(v.Link0[l, 1, i] * v.Link0[l, 2, i+1] * conj(v.LinkN[l, 1, i]) * conj(v.Link0[l, 2, i]))
            dphi[l] = angle(v.Link0[l, 1, i]) + angle(v.Link0[l, 2, i+1]) - angle(v.LinkN[l, 1, i]) - angle(v.Link0[l, 2, i])
        end
    elseif i == N && j == N
        for l in 1:Hs
            phi[l] = angle(v.Link0[l, 1, N] * v.Link0[l, 2, 1] * conj(v.LinkN[l, 1, N]) * conj(v.Link0[l, 2, N]))
            dphi[l] = angle(v.Link0[l, 1, N]) + angle(v.Link0[l, 2, 1]) - angle(v.LinkN[l, 1, N]) - angle(v.Link0[l, 2, N])
        end
    end

    phi .= (phi .- dphi) ./ 2pi
end

@views function ChernPhase!(TopologicalNumber, p::Params) # chern number # Bug
    @unpack N, Hs = p
    TopologicalNumber[:] .= zero(Float64)

    k = zeros(2)
    Link0 = zeros(ComplexF64, Hs, 2, N)
    Link1 = zeros(ComplexF64, Hs, 2, N)
    LinkN = zeros(ComplexF64, Hs, 2, N)
    link10 = zeros(ComplexF64, Hs)
    link01 = zeros(ComplexF64, Hs)

    psi_0 = zeros(ComplexF64, N, Hs, Hs)
    psi_1 = zeros(ComplexF64, N, Hs, Hs)
    psi_N = zeros(ComplexF64, N, Hs, Hs)
    Evec0 = zeros(N, Hs)
    Evec1 = zeros(N, Hs)

    psi00 = zeros(ComplexF64, Hs, Hs)
    psi10 = zeros(ComplexF64, Hs, Hs)
    psi01 = zeros(ComplexF64, Hs, Hs)
    Enevec = zeros(Hs)

    v = TemporalFirstChern(k, Link0, Link1, LinkN, link10, link01, psi_0, psi_1, psi_N, Evec0, Evec1, psi00, psi10, psi01, Enevec)

    phi = zeros(Hs)
    dphi = zeros(Hs)

    for j in 1:N
        U!(j, v, p)
        for i in 1:N
            F!(phi, dphi, i, j, v, p)
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
    p = Params(; Ham=Hamiltonian, N, gapless, rounds, Hs, dim=2)

    TopologicalNumber = zeros(Hs)

    ChernPhase!(TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
    end

    Total = sum(TopologicalNumber)

    (; TopologicalNumber, Total)
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
function solve(prob::FCProblem, alg::T1=FHS(); parallel::T2=UseSingleThread()) where {T1<:FirstChernAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack H, N, gapless, rounds = prob

    Hs = size(H(zeros(2)), 1)
    p = Params(; Ham=H, N, gapless, rounds, Hs, dim=2)

    TopologicalNumber = zeros(Hs)

    ChernPhase!(TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
    end

    Total = sum(TopologicalNumber)

    FCSolution(; TopologicalNumber, Total)
end