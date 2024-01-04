function psi_j!(j, v::TemporalZ2, p::Params) # wave function
    @unpack Ham, N = p
    for i in 1:N
        v.k[1] = (i - 1) * 2pi / N
        v.k[2] = (j - 1) * 2pi / N
        v.psi_1[i, :, :] .= eigen!(Ham(v.k)).vectors
    end
end

@views function Link!(v::TemporalZ2)
    for l in 1:v.Hshalf
        v.link10[l] = det(v.psi00[:, 2l-1:2l]' * v.psi10[:, 2l-1:2l])
        v.link01[l] = det(v.psi00[:, 2l-1:2l]' * v.psi01[:, 2l-1:2l])
    end
end

@views function U!(j, v::TemporalZ2, p::Params) # link variable
    @unpack N = p

    if j != 1
        v.Link1 .= v.Link2
    end

    if j != N
        if j == 1

            psi_j!(1, v, p)
            v.psi_0 .= v.psi_1
            v.psi_N .= v.psi_1
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

                Link!(v)

                v.Link1[:, 1, i] .= v.link10
                v.Link1[:, 2, i] .= v.link01
                v.LinkN1 .= v.Link1
            end

            v.w00 .= v.psi_0[1, :, :]' * v.T * conj(v.psi_0[1, :, :])
            v.wp0 .= v.psi_0[v.Nhalf, :, :]' * v.T * conj(v.psi_0[v.Nhalf, :, :])
        end

        v.psi_0 .= v.psi_1

        if j != N - 1
            psi_j!(j + 2, v, p)
        end

        if j == v.Nhalf - 1
            v.w0p .= v.psi_0[1, :, :]' * v.T * conj(v.psi_0[1, :, :])
            v.wpp .= v.psi_0[v.Nhalf, :, :]' * v.T * conj(v.psi_0[v.Nhalf, :, :])
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

            Link!(v)

            v.Link2[:, 1, i] .= v.link10
            v.Link2[:, 2, i] .= v.link01
        end
    end
end

@views function F!(phi, Px0, Pxp, i, j, v::TemporalZ2, p::Params) # lattice field strength
    @unpack N = p

    if i != N && j != N
        for l in 1:v.Hshalf
            phi[l] = angle(v.Link1[l, 1, i] * v.Link1[l, 2, i+1] * conj(v.Link2[l, 1, i]) * conj(v.Link1[l, 2, i]))
        end
    elseif i == N && j != N
        for l in 1:v.Hshalf
            phi[l] = angle(v.Link1[l, 1, N] * v.Link1[l, 2, 1] * conj(v.Link2[l, 1, N]) * conj(v.Link1[l, 2, N]))
        end
    elseif i != N && j == N
        for l in 1:v.Hshalf
            phi[l] = angle(v.Link1[l, 1, i] * v.Link1[l, 2, i+1] * conj(v.LinkN1[l, 1, i]) * conj(v.Link1[l, 2, i]))
        end
    elseif i == N && j == N
        for l in 1:v.Hshalf
            phi[l] = angle(v.Link1[l, 1, N] * v.Link1[l, 2, 1] * conj(v.LinkN1[l, 1, N]) * conj(v.Link1[l, 2, N]))
        end
    end

    if j == 1 && i < v.Nhalf
        for l in 1:v.Hshalf
            Px0[l] -= angle(v.Link1[l, 1, i])
        end
    elseif j == v.Nhalf && i < v.Nhalf
        for l in 1:v.Hshalf
            Pxp[l] -= angle(v.Link1[l, 1, i])
        end
    end
end

# @views function Z2Phase_round!(TopologicalNumber, p::Params) # chern number
#     @unpack N, Hs, rounds = p
#     Nhalf = N ÷ 2 + 1
#     Hshalf = Hs ÷ 2

#     Link1 = zeros(ComplexF64, Hshalf, 2, N)
#     Link2 = zeros(ComplexF64, Hshalf, 2, N)
#     LinkN1 = zeros(ComplexF64, Hshalf, 2, N)
#     link10 = zeros(ComplexF64, Hshalf)
#     link01 = zeros(ComplexF64, Hshalf)

#     psi_0 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
#     psi_1 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
#     psi_N = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)

#     psi00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     psi10 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     psi01 = zeros(ComplexF64, 2Hshalf, 2Hshalf)

#     phi = zeros(Hshalf)
#     Px0 = zeros(Hshalf)
#     Pxp = zeros(Hshalf)

#     s0 = Matrix{ComplexF64}(I, Hshalf, Hshalf)
#     sy = [0 -1; 1 0]
#     T = kron(s0, sy)

#     w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

#     TN = zeros(Hshalf)

#     for j in 1:Nhalf
#         U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, j, Hshalf, Nhalf, p)
#         for i in 1:N
#             F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, Hshalf, Nhalf, p)
#             if j < Nhalf
#                 TN[:] .+= phi[:]
#             end
#         end
#     end

#     for l in 1:Hshalf
#         Px0[l] += angle((w00[2l-1, 2l]) / (wp0[2l-1, 2l]))
#         Pxp[l] += angle((w0p[2l-1, 2l]) / (wpp[2l-1, 2l]))
#     end

#     # for l in 1:Hshalf
#     #     Px0[l] += imag(log((w00[2l-1, 2l]) / (wp0[2l-1, 2l])))
#     #     Pxp[l] += imag(log((w0p[2l-1, 2l]) / (wpp[2l-1, 2l])))
#     # end

#     for l in 1:Hshalf
#         if TN[l] - 2Px0[l] + 2Pxp[l] !== NaN
#             TopologicalNumber[l] = abs(rem(round(Int, (TN[l] - 2Px0[l] + 2Pxp[l]) / 2pi), 2))
#         end
#     end
# end

@doc raw"""
"""
@views function Z2Phase!(TopologicalNumber, p::Params) # chern number
    @unpack N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1
    Hshalf = Hs ÷ 2

    k = zeros(2)

    Link1 = zeros(ComplexF64, Hshalf, 2, N)
    Link2 = zeros(ComplexF64, Hshalf, 2, N)
    LinkN1 = zeros(ComplexF64, Hshalf, 2, N)
    link10 = zeros(ComplexF64, Hshalf)
    link01 = zeros(ComplexF64, Hshalf)

    psi_0 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    psi_1 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    psi_N = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)

    psi00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    psi10 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    psi01 = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    phi = zeros(Hshalf)
    Px0 = zeros(Hshalf)
    Pxp = zeros(Hshalf)

    s0 = Matrix{ComplexF64}(I, Hshalf, Hshalf)
    sy = [0 -1; 1 0] # imaginary???
    T = kron(s0, sy)

    w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    v = TemporalZ2(k, T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, Nhalf, Hshalf)

    TN = zeros(Hshalf)

    for j in 1:Nhalf
        U!(j, v, p)
        for i in 1:N
            F!(phi, Px0, Pxp, i, j, v, p)
            if j < Nhalf
                TN[:] .+= phi[:]
            end
        end
    end

    for l in 1:Hshalf
        Px0[l] += angle((v.w00[2l-1, 2l]) / (v.wp0[2l-1, 2l]))
        Pxp[l] += angle((v.w0p[2l-1, 2l]) / (v.wpp[2l-1, 2l]))

        if TN[l] - 2Px0[l] + 2Pxp[l] !== NaN
            TopologicalNumber[l] = 1 - abs(1 - rem(abs(TN[l] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
        end
    end
end

# @views function Z2Phase_round!(TopologicalNumber, TRTopologicalNumber, p::Params) # chern number
#     @unpack N, Hs, rounds = p
#     Nhalf = N ÷ 2 + 1
#     Hshalf = Hs ÷ 2

#     Link1 = zeros(ComplexF64, Hshalf, 2, N)
#     Link2 = zeros(ComplexF64, Hshalf, 2, N)
#     LinkN1 = zeros(ComplexF64, Hshalf, 2, N)
#     link10 = zeros(ComplexF64, Hshalf)
#     link01 = zeros(ComplexF64, Hshalf)

#     psi_0 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
#     psi_1 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
#     psi_N = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)

#     psi00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     psi10 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     psi01 = zeros(ComplexF64, 2Hshalf, 2Hshalf)

#     phi = zeros(Hshalf)
#     Px0 = zeros(Hshalf)
#     Pxp = zeros(Hshalf)

#     s0 = Matrix{ComplexF64}(I, Hshalf, Hshalf)
#     sy = [0 -1; 1 0]
#     T = kron(s0, sy)

#     w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

#     TN = zeros(Hshalf, 2)

#     for j in 1:N
#         U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, j, Hshalf, Nhalf, p)
#         for i in 1:N
#             F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, Hshalf, Nhalf, p)
#             if j < Nhalf
#                 TN[:, 1] .+= phi[:]
#             else
#                 TN[:, 2] .+= phi[:]
#             end
#         end
#     end

#     for l in 1:Hshalf
#         Px0[l] += angle((w00[2l-1, 2l]) / (wp0[2l-1, 2l]))
#         Pxp[l] += angle((w0p[2l-1, 2l]) / (wpp[2l-1, 2l]))
#     end

#     # for l in 1:Hshalf
#     #     Px0[l] += imag(log((w00[2l-1, 2l]) / (wp0[2l-1, 2l])))
#     #     Pxp[l] += imag(log((w0p[2l-1, 2l]) / (wpp[2l-1, 2l])))
#     # end

#     for l in 1:Hshalf
#         if TN[l, 1] - 2Px0[l] + 2Pxp[l] !== NaN
#             TopologicalNumber[l] = abs(rem(round(Int, (TN[l, 1] - 2Px0[l] + 2Pxp[l]) / 2pi), 2))
#             TRTopologicalNumber[l] = abs(rem(round(Int, (TN[l, 2] - 2Px0[l] + 2Pxp[l]) / 2pi), 2))
#         end
#     end
# end

@doc raw"""
"""
@views function Z2Phase!(TopologicalNumber, TRTopologicalNumber, p::Params) # chern number
    @unpack N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1
    Hshalf = Hs ÷ 2

    k = zeros(2)

    Link1 = zeros(ComplexF64, Hshalf, 2, N)
    Link2 = zeros(ComplexF64, Hshalf, 2, N)
    LinkN1 = zeros(ComplexF64, Hshalf, 2, N)
    link10 = zeros(ComplexF64, Hshalf)
    link01 = zeros(ComplexF64, Hshalf)

    psi_0 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    psi_1 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    psi_N = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)

    psi00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    psi10 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    psi01 = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    phi = zeros(Hshalf)
    Px0 = zeros(Hshalf)
    Pxp = zeros(Hshalf)

    s0 = Matrix{ComplexF64}(I, Hshalf, Hshalf)
    sy = [0 -1; 1 0] # imaginary???
    T = kron(s0, sy)

    w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    v = TemporalZ2(k, T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, Nhalf, Hshalf)

    TN = zeros(Hshalf, 2)

    for j in 1:N
        U!(j, v, p)
        for i in 1:N
            F!(phi, Px0, Pxp, i, j, v, p)
            if j < Nhalf
                TN[:, 1] .+= phi[:]
            else
                TN[:, 2] .+= phi[:]
            end
        end
    end

    for l in 1:Hshalf
        Px0[l] += angle((v.w00[2l-1, 2l]) / (v.wp0[2l-1, 2l]))
        Pxp[l] += angle((v.w0p[2l-1, 2l]) / (v.wpp[2l-1, 2l]))

        if TN[l] - 2Px0[l] + 2Pxp[l] !== NaN
            TopologicalNumber[l] = 1 - abs(1 - rem(abs(TN[l, 1] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
            TRTopologicalNumber[l] = 1 - abs(1 - rem(abs(TN[l, 2] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
        end
    end
end

@doc raw"""

 Calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case with reference to Shiozaki method [Shiozaki2023discrete](@cite).

    calcZ2(Hamiltonian::Function; N::Int=50, rounds::Bool=true, TR::Bool=false)

 Arguments
 - `Hamiltonian::Function` is a matrix with one-dimensional wavenumber `k` as an argument.
 - `N::Int` is the number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - `rounds::Bool` is an option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.

# Definition
 The $\mathbb{Z}_{2}$ number of the $2n$th (and $2n-1$th) band $\nu_{n}$ is defined by
```math
\nu_{n}=F_{n}-\left(P_{n}(0)-P_{n}(\pi)\right)
```
 $F_{n}$ is the Berry flux of the $n$th band in the $\mathrm{BZ}'$. The range $\mathrm{BZ}'$ is $\bm{k}\in[0,2\pi]\times[0,\pi]$ half of BZ(Brillouin Zone).
```math
F_{n}=\frac{1}{2\pi}\sum_{\bm{k}\in\mathrm{BZ}'}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,1}^{*}(\bm{k})\right]\right]
```
 $P_{n}(k_{2})$ is the time-reversal polarization at wavenumber $k_{2}$.
```math
P_{n}(k_{2})=\frac{1}{2\pi}\frac{\mathrm{Pf}[\omega(0,k_{2})]}{\mathrm{Pf}[\omega(\pi,k_{2})]}\sum_{k_{1}=0}^{\pi-e_{1}}U_{n,1}(\bm{k})
```
 $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $2n$th (and $2n-1$th) band. $\omega(\bm{k})$ is the unitary matrix given by
```math
\omega(\bm{k})=\bra{\Psi(-\bm{k})}T\ket{\Psi(\bm{k})}
```
 $T$ is the time-reversal operator.
"""
function calcZ2(Hamiltonian::Function; N::Int=50, rounds::Bool=true, TR::Bool=false)

    Hs = size(Hamiltonian(zeros(2)), 1)
    Hshalf = Hs ÷ 2
    p = Params(; Ham=Hamiltonian, N, Hs, gapless=0.0, rounds, dim=2)

    TopologicalNumber = zeros(Hshalf)
    if TR == false
        Z2Phase!(TopologicalNumber, p)

        if rounds == true
            TopologicalNumber = round.(Int, TopologicalNumber)
            Total = rem(sum(TopologicalNumber), 2)
        else
            Total = abs(sum(TopologicalNumber))
            Total = abs(rem(Total, 2))
        end

        (; TopologicalNumber, Total)
    else
        TRTopologicalNumber = zeros(Hshalf)
        Z2Phase!(TopologicalNumber, TRTopologicalNumber, p)

        if rounds == true
            TopologicalNumber = round.(Int, TopologicalNumber)
            TRTopologicalNumber = round.(Int, TRTopologicalNumber)
            Total = rem(sum(TopologicalNumber), 2)
        else
            Total = abs(sum(TopologicalNumber))
            Total = rem(1 - abs(1 - Total), 2)
        end

        (; TopologicalNumber, TRTopologicalNumber, Total)
    end
end


@doc raw"""

 Calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case with reference to Shiozaki method [Shiozaki2023discrete](@cite).

    solve(prob::Z2Problem, alg::T1=Shio(); parallel::T2=UseSingleThread()) where {T1<:Z2Algorithms,T2<:TopologicalNumbersParallel}

 Arguments
 - `Hamiltonian::Function` is a matrix with one-dimensional wavenumber `k` as an argument.
 - `N::Int` is the number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - `rounds::Bool` is an option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.

# Definition
 The $\mathbb{Z}_{2}$ number of the $2n$th (and $2n-1$th) band $\nu_{n}$ is defined by
```math
\nu_{n}=F_{n}-\left(P_{n}(0)-P_{n}(\pi)\right)
```
 $F_{n}$ is the Berry flux of the $n$th band in the $\mathrm{BZ}'$. The range $\mathrm{BZ}'$ is $\bm{k}\in[0,2\pi]\times[0,\pi]$ half of BZ(Brillouin Zone).
```math
F_{n}=\frac{1}{2\pi}\sum_{\bm{k}\in\mathrm{BZ}'}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,1}^{*}(\bm{k})\right]\right]
```
 $P_{n}(k_{2})$ is the time-reversal polarization at wavenumber $k_{2}$.
```math
P_{n}(k_{2})=\frac{1}{2\pi}\frac{\mathrm{Pf}[\omega(0,k_{2})]}{\mathrm{Pf}[\omega(\pi,k_{2})]}\sum_{k_{1}=0}^{\pi-e_{1}}U_{n,1}(\bm{k})
```
 $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $2n$th (and $2n-1$th) band. $\omega(\bm{k})$ is the unitary matrix given by
```math
\omega(\bm{k})=\bra{\Psi(-\bm{k})}T\ket{\Psi(\bm{k})}
```
 $T$ is the time-reversal operator.
"""
function solve(
    prob::Z2Problem,
    alg::T1=Shio();
    parallel::T2=UseSingleThread()
) where {T1<:Z2Algorithms,T2<:TopologicalNumbersParallel}
    @unpack H, N, rounds, TR = prob

    Hs = size(H(zeros(2)), 1)
    Hshalf = Hs ÷ 2
    p = Params(; Ham=H, N, Hs, rounds, dim=2)

    TopologicalNumber = zeros(Hshalf)
    if TR == false
        Z2Phase!(TopologicalNumber, p)

        if rounds == true
            TopologicalNumber = round.(Int, TopologicalNumber)
            Total = rem(sum(TopologicalNumber), 2)
        else
            Total = abs(sum(TopologicalNumber))
            Total = abs(rem(Total, 2))
        end

        Z2Solution(; TopologicalNumber, Total)
    else
        TRTopologicalNumber = zeros(Hshalf)
        Z2Phase!(TopologicalNumber, TRTopologicalNumber, p)

        if rounds == true
            TopologicalNumber = round.(Int, TopologicalNumber)
            TRTopologicalNumber = round.(Int, TRTopologicalNumber)
            Total = rem(sum(TopologicalNumber), 2)
        else
            Total = abs(sum(TopologicalNumber))
            Total = rem(1 - abs(1 - Total), 2)
        end

        Z2Solution(; TopologicalNumber, TRTopologicalNumber, Total)
    end
end