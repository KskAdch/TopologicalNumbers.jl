function psi_j!(j, psi_1, p::Params) # wave function
    @unpack Hamiltonian, N = p
    for i in 1:N
        k = [i-1, j-1] * 2pi / N .+ 2pi * [1e-5, 1e-5]
        psi_1[i, :, :] .= eigen!(Hamiltonian(k)).vectors
    end
end

@views function Link!(psi00, psi10, psi01, link10, link01, Hshalf)
    for l in 1:Hshalf
        link10[l] = det(psi00[:, 2l-1:2l]' * psi10[:, 2l-1:2l])
        link01[l] = det(psi00[:, 2l-1:2l]' * psi01[:, 2l-1:2l])
    end
end

@views function U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, j, Hshalf, Nhalf, p::Params) # link variable
    @unpack N = p

    if j != 1
        Link1 .= Link2
    end

    if j != N
        if j == 1

            psi_j!(1, psi_1, p)
            psi_0 .= psi_1
            psi_N .= psi_1
            psi_j!(2, psi_1, p)

            for i in 1:N

                psi00[:, :] = psi_0[i, :, :]
                if i == N
                    psi10[:, :] = psi_0[1, :, :]
                    psi01[:, :] = psi_1[N, :, :]
                else
                    psi10[:, :] = psi_0[i+1, :, :]
                    psi01[:, :] = psi_1[i, :, :]
                end

                Link!(psi00, psi10, psi01, link10, link01, Hshalf)

                Link1[:, :, i] .= [link10 link01]
                LinkN1 .= Link1
            end

            w00 .= psi_0[1, :, :]' * T * conj(psi_0[1, :, :])
            wp0 .= psi_0[Nhalf, :, :]' * T * conj(psi_0[Nhalf, :, :])
        end

        psi_0 .= psi_1

        if j != N - 1
            psi_j!(j + 2, psi_1, p)
        end

        if j == Nhalf - 1
            w0p .= psi_0[1, :, :]' * T * conj(psi_0[1, :, :])
            wpp .= psi_0[Nhalf, :, :]' * T * conj(psi_0[Nhalf, :, :])
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

            Link!(psi00, psi10, psi01, link10, link01, Hshalf)

            Link2[:, :, i] .= [link10 link01]
        end
    end
end

@views function F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, Hshalf, Nhalf, p::Params) # lattice field strength
    @unpack N = p

    if i == N && j == N
        phi[:] = [angle(Link1[l, 1, N] * Link1[l, 2, 1] * conj(LinkN1[l, 1, N]) * conj(Link1[l, 2, N])) for l in 1:Hshalf]
    elseif i == N
        phi[:] = [angle(Link1[l, 1, N] * Link1[l, 2, 1] * conj(Link2[l, 1, N]) * conj(Link1[l, 2, N])) for l in 1:Hshalf]
    elseif j == N
        phi[:] = [angle(Link1[l, 1, i] * Link1[l, 2, i+1] * conj(LinkN1[l, 1, i]) * conj(Link1[l, 2, i])) for l in 1:Hshalf]
    else
        phi[:] = [angle(Link1[l, 1, i] * Link1[l, 2, i+1] * conj(Link2[l, 1, i]) * conj(Link1[l, 2, i])) for l in 1:Hshalf]
    end

    # if i == N && j == N
    #     phi[:] = [imag(log(Link1[l, 1, N] * Link1[l, 2, 1] * conj(LinkN1[l, 1, N]) * conj(Link1[l, 2, N]))) for l in 1:Hshalf]
    # elseif i == N
    #     phi[:] = [imag(log(Link1[l, 1, N] * Link1[l, 2, 1] * conj(Link2[l, 1, N]) * conj(Link1[l, 2, N]))) for l in 1:Hshalf]
    # elseif j == N
    #     phi[:] = [imag(log(Link1[l, 1, i] * Link1[l, 2, i+1] * conj(LinkN1[l, 1, i]) * conj(Link1[l, 2, i]))) for l in 1:Hshalf]
    # else
    #     phi[:] = [imag(log(Link1[l, 1, i] * Link1[l, 2, i+1] * conj(Link2[l, 1, i]) * conj(Link1[l, 2, i]))) for l in 1:Hshalf]
    # end

    if j == 1 && i < Nhalf
        Px0 .-= [angle(Link1[l, 1, i]) for l in 1:Hshalf]
    elseif j == Nhalf && i < Nhalf
        Pxp .-= [angle(Link1[l, 1, i]) for l in 1:Hshalf]
    end

    # if j == 1 && i < Nhalf
    #     Px0 .-= [imag(log(Link1[l, 1, i])) for l in 1:Hshalf]
    # elseif j == Nhalf && i < Nhalf
    #     Pxp .-= [imag(log(Link1[l, 1, i])) for l in 1:Hshalf]
    # end
end

@views function Z2Phase_round!(TopologicalNumber, p::Params) # chern number
    @unpack N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1
    Hshalf = Hs ÷ 2

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
    sy = [0 -1; 1 0]
    T = kron(s0, sy)

    w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    TN = zeros(Hshalf)

    for j in 1:Nhalf
        U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, j, Hshalf, Nhalf, p)
        for i in 1:N
            F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, Hshalf, Nhalf, p)
            if j < Nhalf
                TN[:] .+= phi[:]
            end
        end
    end

    for l in 1:Hshalf
        Px0[l] += angle((w00[2l-1, 2l]) / (wp0[2l-1, 2l]))
        Pxp[l] += angle((w0p[2l-1, 2l]) / (wpp[2l-1, 2l]))
    end

    # for l in 1:Hshalf
    #     Px0[l] += imag(log((w00[2l-1, 2l]) / (wp0[2l-1, 2l])))
    #     Pxp[l] += imag(log((w0p[2l-1, 2l]) / (wpp[2l-1, 2l])))
    # end

    for l in 1:Hshalf
        if TN[l] - 2Px0[l] + 2Pxp[l] !== NaN
            TopologicalNumber[l] = abs(rem(round(Int, (TN[l] - 2Px0[l] + 2Pxp[l]) / 2pi), 2))
        end
    end
end

@views function Z2Phase!(TopologicalNumber, p::Params) # chern number
    @unpack N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1
    Hshalf = Hs ÷ 2

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
    sy = [0 -1; 1 0]
    T = kron(s0, sy)

    w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    TN = zeros(Hshalf)

    for j in 1:Nhalf
        U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, j, Hshalf, Nhalf, p)
        for i in 1:N
            F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, Hshalf, Nhalf, p)
            if j < Nhalf
                TN[:] .+= phi[:]
            end
        end
    end

    for l in 1:Hshalf
        Px0[l] += angle((w00[2l-1, 2l]) / (wp0[2l-1, 2l]))
        Pxp[l] += angle((w0p[2l-1, 2l]) / (wpp[2l-1, 2l]))
    end

    # for l in 1:Hshalf
    #     Px0[l] += imag(log((w00[2l-1, 2l]) / (wp0[2l-1, 2l])))
    #     Pxp[l] += imag(log((w0p[2l-1, 2l]) / (wpp[2l-1, 2l])))
    # end

    for l in 1:Hshalf
        if TN[l] - 2Px0[l] + 2Pxp[l] !== NaN
            TopologicalNumber[l] = abs(rem((TN[l] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
        end
    end
end

@views function Z2Phase_round!(TopologicalNumber, TRTopologicalNumber, p::Params) # chern number
    @unpack N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1
    Hshalf = Hs ÷ 2

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
    sy = [0 -1; 1 0]
    T = kron(s0, sy)

    w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    TN = zeros(Hshalf, 2)

    for j in 1:N
        U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, j, Hshalf, Nhalf, p)
        for i in 1:N
            F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, Hshalf, Nhalf, p)
            if j < Nhalf
                TN[:, 1] .+= phi[:]
            else
                TN[:, 2] .+= phi[:]
            end
        end
    end

    for l in 1:Hshalf
        Px0[l] += angle((w00[2l-1, 2l]) / (wp0[2l-1, 2l]))
        Pxp[l] += angle((w0p[2l-1, 2l]) / (wpp[2l-1, 2l]))
    end

    # for l in 1:Hshalf
    #     Px0[l] += imag(log((w00[2l-1, 2l]) / (wp0[2l-1, 2l])))
    #     Pxp[l] += imag(log((w0p[2l-1, 2l]) / (wpp[2l-1, 2l])))
    # end

    for l in 1:Hshalf
        if TN[l, 1] - 2Px0[l] + 2Pxp[l] !== NaN
            TopologicalNumber[l] = abs(rem(round(Int, (TN[l, 1] - 2Px0[l] + 2Pxp[l]) / 2pi), 2))
            TRTopologicalNumber[l] = abs(rem(round(Int, (TN[l, 2] - 2Px0[l] + 2Pxp[l]) / 2pi), 2))
        end
    end
end

@views function Z2Phase!(TopologicalNumber, TRTopologicalNumber, p::Params) # chern number
    @unpack N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1
    Hshalf = Hs ÷ 2

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
    sy = [0 -1; 1 0]
    T = kron(s0, sy)

    w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    TN = zeros(Hshalf, 2)

    for j in 1:N
        U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, j, Hshalf, Nhalf, p)
        for i in 1:N
            F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, Hshalf, Nhalf, p)
            if j < Nhalf
                TN[:, 1] .+= phi[:]
            else
                TN[:, 2] .+= phi[:]
            end
        end
    end

    for l in 1:Hshalf
        Px0[l] += angle((w00[2l-1, 2l]) / (wp0[2l-1, 2l]))
        Pxp[l] += angle((w0p[2l-1, 2l]) / (wpp[2l-1, 2l]))
    end
    
    # for l in 1:Hshalf
    #     Px0[l] += imag(log((w00[2l-1, 2l]) / (wp0[2l-1, 2l])))
    #     Pxp[l] += imag(log((w0p[2l-1, 2l]) / (wpp[2l-1, 2l])))
    # end

    for l in 1:Hshalf
        if TN[l, 1] - 2Px0[l] + 2Pxp[l] !== NaN
            TopologicalNumber[l] = abs(rem((TN[l, 1] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
            TRTopologicalNumber[l] = abs(rem((TN[l, 2] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
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

    n0 = zeros(2)
    Hs = size(Hamiltonian(n0))[1]
    Hshalf = Hs ÷ 2
    p = Params(; Hamiltonian, N, Hs, gapless=0.0, rounds, dim=2)

    if TR == false
        if rounds == true
            TopologicalNumber = zeros(Int, Hshalf)

            Z2Phase_round!(TopologicalNumber, p)

            Total = rem(sum(TopologicalNumber), 2)
        else
            TopologicalNumber = zeros(Hshalf)

            Z2Phase!(TopologicalNumber, p)

            Total = abs(sum(TopologicalNumber))
            while Total > 1.5
                Total -= 2
            end
            Total = abs(rem(Total, 2))
        end

        (; TopologicalNumber, Total)
    else
        if rounds == true
            TopologicalNumber = zeros(Int, Hshalf)
            TRTopologicalNumber = zeros(Int, Hshalf)

            Z2Phase_round!(TopologicalNumber, TRTopologicalNumber, p)

            Total = rem(sum(TopologicalNumber), 2)
        else
            TopologicalNumber = zeros(Hshalf)
            TRTopologicalNumber = zeros(Hshalf)

            Z2Phase!(TopologicalNumber, TRTopologicalNumber, p)

            Total = abs(sum(TopologicalNumber))
            while Total > 1.5
                Total -= 2
            end
            Total = abs(rem(Total, 2))
        end

        (; TopologicalNumber, TRTopologicalNumber, Total)
    end
end
