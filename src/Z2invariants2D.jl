function Z2invariants2D(Hamiltonian::Function; N::Int=50, rounds::Bool=true)

    # function psi_j!(j, psi_2, Evec2, p) # wave function
    function psi_j!(j, psi_2, p) # wave function
        @unpack Hamiltonian, N, Hs = p
        for i in 1:N
            k = [i - 1, j - 1] * 2pi / N
            eigens = eigen!(Hamiltonian(k))
            psi_2[i, :, :] .= eigens.vectors
            # Evec2[i, :] .= eigens.values
        end
    end

    # @views function U!(w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link1, link2, psi_1, psi_2, psi_N1, Evec1, Evec2, psi0, psi1, psi2, Enevec, j, p) # link variable
    @views function U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link1, link2, psi_1, psi_2, psi_N1, psi0, psi1, psi2, j, p) # link variable
        @unpack N, Nhalf, Hshalf = p

        if j != 1
            Link1 .= Link2
        end

        if j != N
            if j == 1

                psi_j!(1, psi_2, p)
                psi_1 .= psi_2
                psi_N1 .= psi_2
                # Evec1 .= Evec2
                psi_j!(2, psi_2, p)

                for i in 1:N

                    psi0[:, :] = psi_1[i, :, :]
                    if i == N
                        psi1[:, :] = psi_1[1, :, :]
                        psi2[:, :] = psi_2[N, :, :]
                    else
                        psi1[:, :] = psi_1[i+1, :, :]
                        psi2[:, :] = psi_2[i, :, :]
                    end

                    for l in 1:Hshalf
                        link1[l] = det(psi0[:, 2l-1:2l]' * psi1[:, 2l-1:2l])
                        link2[l] = det(psi0[:, 2l-1:2l]' * psi2[:, 2l-1:2l])
                    end

                    Link1[:, :, i] .= [link1 link2]
                    LinkN1 .= Link1
                end

                w00 .= psi_1[1, :, :]' * T * conj(psi_1[1, :, :])
                wp0 .= psi_1[Nhalf, :, :]' * T * conj(psi_1[Nhalf, :, :])
            end

            psi_1 .= psi_2
            # Evec1 .= Evec2

            if j != N - 1
                psi_j!(j + 2, psi_2, p)
            end

            if j == Nhalf - 1
                w0p .= psi_1[1, :, :]' * T * conj(psi_1[1, :, :])
                wpp .= psi_1[Nhalf, :, :]' * T * conj(psi_1[Nhalf, :, :])
            end

            for i in 1:N

                psi0[:, :] = psi_1[i, :, :]
                if i == N && j == N - 1
                    psi1[:, :] = psi_1[1, :, :]
                    psi2[:, :] = psi_N1[N, :, :]
                elseif i == N
                    psi1[:, :] = psi_1[1, :, :]
                    psi2[:, :] = psi_2[N, :, :]
                elseif j == N - 1
                    psi1[:, :] = psi_1[i+1, :, :]
                    psi2[:, :] = psi_N1[i, :, :]
                else
                    psi1[:, :] = psi_1[i+1, :, :]
                    psi2[:, :] = psi_2[i, :, :]
                end

                for l in 1:Hshalf
                    link1[l] = det(psi0[:, 2l-1:2l]' * psi1[:, 2l-1:2l])
                    link2[l] = det(psi0[:, 2l-1:2l]' * psi2[:, 2l-1:2l])
                end

                Link2[:, :, i] .= [link1 link2]
            end
        end
    end

    @views function F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, p) # lattice field strength
        @unpack N, Nhalf, Hshalf = p

        if i == N && j == N
            phi[:] = [imag(log(Link1[l, 1, N] * Link1[l, 2, 1] * conj(LinkN1[l, 1, N]) * conj(Link1[l, 2, N]))) for l in 1:Hshalf]
        elseif i == N
            phi[:] = [imag(log(Link1[l, 1, N] * Link1[l, 2, 1] * conj(Link2[l, 1, N]) * conj(Link1[l, 2, N]))) for l in 1:Hshalf]
        elseif j == N
            phi[:] = [imag(log(Link1[l, 1, i] * Link1[l, 2, i+1] * conj(LinkN1[l, 1, i]) * conj(Link1[l, 2, i]))) for l in 1:Hshalf]
        else
            phi[:] = [imag(log(Link1[l, 1, i] * Link1[l, 2, i+1] * conj(Link2[l, 1, i]) * conj(Link1[l, 2, i]))) for l in 1:Hshalf]
        end

        if j == 1 && i < Nhalf
            Px0 .-= [imag(log(Link1[l, 1, i])) for l in 1:Hshalf]
        elseif j == Nhalf && i < Nhalf
            Pxp .-= [imag(log(Link1[l, 1, i])) for l in 1:Hshalf]
        end
    end

    @views function Phase!(TopologicalNumber, p) # chern number
        @unpack N, Nhalf, Hs, Hshalf = p
        Link1 = zeros(ComplexF64, Hshalf, 2, N)
        Link2 = zeros(ComplexF64, Hshalf, 2, N)
        LinkN1 = zeros(ComplexF64, Hshalf, 2, N)
        link1 = zeros(ComplexF64, Hshalf)
        link2 = zeros(ComplexF64, Hshalf)

        psi_1 = zeros(ComplexF64, N, Hs, Hs)
        psi_2 = zeros(ComplexF64, N, Hs, Hs)
        psi_N1 = zeros(ComplexF64, N, Hs, Hs)
        # Evec1 = zeros(N, Hs)
        # Evec2 = zeros(N, Hs)

        psi0 = zeros(ComplexF64, Hs, Hs)
        psi1 = zeros(ComplexF64, Hs, Hs)
        psi2 = zeros(ComplexF64, Hs, Hs)
        # Enevec = zeros(Hs)

        phi = zeros(Hshalf)
        Px0 = zeros(Hshalf)
        Pxp = zeros(Hshalf)

        s0 = Matrix{ComplexF64}(I, Int(Hs / 2), Int(Hs / 2))
        sy = [0 -1; 1 0]
        T = kron(s0, sy)

        w00 = zeros(ComplexF64, Hs, Hs)
        w0p = zeros(ComplexF64, Hs, Hs)
        wp0 = zeros(ComplexF64, Hs, Hs)
        wpp = zeros(ComplexF64, Hs, Hs)

        TN = zeros(Hshalf, 2)

        for j in 1:N
            U!(T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link1, link2, psi_1, psi_2, psi_N1, psi0, psi1, psi2, j, p)
            for i in 1:N
                F!(phi, Px0, Pxp, i, j, Link1, Link2, LinkN1, p)
                if j < Nhalf
                    TN[:, 1] .+= phi[:]
                else
                    TN[:, 2] .+= phi[:]
                end
            end
        end

        for l in 1:Hshalf
            Px0[l] += imag(log((w00[2l-1, 2l]) / (wp0[2l-1, 2l])))
            Pxp[l] += imag(log((w0p[2l-1, 2l]) / (wpp[2l-1, 2l])))
        end

        for l in 1:Hshalf
            for TRS in 1:2
                if TN[l, TRS] - 2Px0[l] + 2Pxp[l] !== NaN
                    if rounds == true
                        TopologicalNumber[TRS, 2l-1:2l] .= abs(rem(round(Int, (TN[l, TRS] - 2Px0[l] + 2Pxp[l]) / 2pi), 2))
                    else
                        TopologicalNumber[TRS, 2l-1:2l] .= abs(rem((TN[l, TRS] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
                    end
                end
            end
        end
    end

    function main(Hamiltonian, N, rounds)

        Nhalf = N ÷ 2 + 1
        n0 = zeros(2)
        Hs = size(Hamiltonian(n0))[1]
        Hshalf = Hs ÷ 2
        p = (; Hamiltonian, N, Nhalf, Hs, Hshalf)

        if rounds == true
            TopologicalNumber = zeros(Int, 2, Hs)
        else
            TopologicalNumber = zeros(2, Hs)
        end

        Phase!(TopologicalNumber, p)

        if rounds == true
            Total = rem(sum(TopologicalNumber[1, 1:2:Hs]), 2)
        else
            Total = sum(TopologicalNumber[1, 1:2:Hs])
            if Total > 1.5
                Total -= 2
            end
            Total = rem(Total, 2)
        end

        (; TopologicalNumber, Total)
    end

    main(Hamiltonian, N, rounds)
end