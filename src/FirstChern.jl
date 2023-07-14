function FirstChern(Hamiltonian::Function; N::Int=51, gapless::Real=0.0)

    function psi_j!(j, psi_2, Evec2, p) # wave function
        @unpack Hamiltonian, N, Hs = p
        for i in 1:N
            k = [i - 1, j - 1] * 2pi / N
            eigens = eigen(Hamiltonian(k))
            psi_2[i, :, :] .= eigens.vectors
            Evec2[i, :] .= eigens.values
        end
    end

    @views function U!(Link1, Link2, LinkN1, link1, link2, psi_1, psi_2, psi_N1, Evec1, Evec2, psi0, psi1, psi2, Enevec, j, p) # link variable
        @unpack N, gapless, Hs = p

        if j != 1
            Link1 .= Link2
        end

        if j != N
            if j == 1

                psi_j!(1, psi_2, Evec2, p)
                psi_1 .= psi_2
                psi_N1 .= psi_2
                Evec1 .= Evec2
                psi_j!(2, psi_2, Evec2, p)

                for i in 1:N

                    psi0[:, :] = psi_1[i, :, :]
                    if i == N
                        psi1[:, :] = psi_1[1, :, :]
                        psi2[:, :] = psi_2[N, :, :]
                    else
                        psi1[:, :] = psi_1[i+1, :, :]
                        psi2[:, :] = psi_2[i, :, :]
                    end

                    Enevec[:] = Evec1[i, :]
                    l = 1
                    while l <= Hs
                        Enevector = Enevec[Enevec.>Enevec[l]] .- Enevec[l]
                        l0 = l + count(Enevector .<= gapless)

                        link1[l:l0] .= det(psi0[:, l:l0]' * psi1[:, l:l0])
                        link2[l:l0] .= det(psi0[:, l:l0]' * psi2[:, l:l0])

                        l = 1 + l0
                    end

                    Link1[:, :, i] .= [link1 link2]
                    LinkN1 .= Link1
                end
            end

            psi_1 .= psi_2
            Evec1 .= Evec2

            if j != N - 1
                psi_j!(j + 2, psi_2, Evec2, p)
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

                Enevec[:] = Evec1[i, :]
                l = 1
                while l <= Hs
                    Enevector = Enevec[Enevec.>Enevec[l]] .- Enevec[l]
                    l0 = l + count(Enevector .<= gapless)

                    link1[l:l0] .= det(psi0[:, l:l0]' * psi1[:, l:l0])
                    link2[l:l0] .= det(psi0[:, l:l0]' * psi2[:, l:l0])

                    l = 1 + l0
                end

                Link2[:, :, i] .= [link1 link2]
            end
        end
    end

    @views function F!(phi, dphi, i, j, Link1, Link2, LinkN1, p) # lattice field strength
        @unpack N, Hs = p

        if i == N && j == N
            phi[:] = [imag(log(Link1[l, 1, N] * Link1[l, 2, 1] * conj(LinkN1[l, 1, N]) * conj(Link1[l, 2, N]))) for l in 1:Hs]
            dphi[:] = [imag(log(Link1[l, 1, N])) + imag(log(Link1[l, 2, 1])) - imag(log(LinkN1[l, 1, N])) - imag(log(Link1[l, 2, N])) for l in 1:Hs]
        elseif i == N
            phi[:] = [imag(log(Link1[l, 1, N] * Link1[l, 2, 1] * conj(Link2[l, 1, N]) * conj(Link1[l, 2, N]))) for l in 1:Hs]
            dphi[:] = [imag(log(Link1[l, 1, N])) + imag(log(Link1[l, 2, 1])) - imag(log(Link2[l, 1, N])) - imag(log(Link1[l, 2, N])) for l in 1:Hs]
        elseif j == N
            phi[:] = [imag(log(Link1[l, 1, i] * Link1[l, 2, i+1] * conj(LinkN1[l, 1, i]) * conj(Link1[l, 2, i]))) for l in 1:Hs]
            dphi[:] = [imag(log(Link1[l, 1, i])) + imag(log(Link1[l, 2, i+1])) - imag(log(LinkN1[l, 1, i])) - imag(log(Link1[l, 2, i])) for l in 1:Hs]
        else
            phi[:] = [imag(log(Link1[l, 1, i] * Link1[l, 2, i+1] * conj(Link2[l, 1, i]) * conj(Link1[l, 2, i]))) for l in 1:Hs]
            dphi[:] = [imag(log(Link1[l, 1, i])) + imag(log(Link1[l, 2, i+1])) - imag(log(Link2[l, 1, i])) - imag(log(Link1[l, 2, i])) for l in 1:Hs]
        end

        phi[:] = [round(Int, (phi[i] - dphi[i]) / 2pi) for i in 1:Hs]
    end

    @views function Phase!(TopologicalNumber, p) # chern number
        @unpack N, Hs = p
        Link1 = zeros(ComplexF64, Hs, 2, N)
        Link2 = zeros(ComplexF64, Hs, 2, N)
        LinkN1 = zeros(ComplexF64, Hs, 2, N)
        link1 = zeros(ComplexF64, Hs)
        link2 = zeros(ComplexF64, Hs)

        psi_1 = zeros(ComplexF64, N, Hs, Hs)
        psi_2 = zeros(ComplexF64, N, Hs, Hs)
        psi_N1 = zeros(ComplexF64, N, Hs, Hs)
        Evec1 = zeros(N, Hs)
        Evec2 = zeros(N, Hs)

        psi0 = zeros(ComplexF64, Hs, Hs)
        psi1 = zeros(ComplexF64, Hs, Hs)
        psi2 = zeros(ComplexF64, Hs, Hs)
        Enevec = zeros(Hs)

        phi = zeros(Hs)
        dphi = zeros(Hs)

        for j in 1:N
            U!(Link1, Link2, LinkN1, link1, link2, psi_1, psi_2, psi_N1, Evec1, Evec2, psi0, psi1, psi2, Enevec, j, p)
            for i in 1:N
                F!(phi, dphi, i, j, Link1, Link2, LinkN1, p)
                TopologicalNumber[:] .+= phi[:]
            end
        end
    end

    function main(Hamiltonian, N, gapless)

        n0 = zeros(2)
        Hs = size(Hamiltonian(n0))[1]
        p = (; Hamiltonian, N, gapless, Hs)

        TopologicalNumber = zeros(Int, Hs)

        Phase!(TopologicalNumber, p)

        Total = sum(TopologicalNumber)

        (; TopologicalNumber, Total)
    end

    main(Hamiltonian, N, gapless)
end