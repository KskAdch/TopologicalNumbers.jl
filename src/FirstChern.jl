function FirstChern(Hamiltonian::Function; N::Int=51, gapless::Real=0.0)

    function psi_j!(j, psi_1, Evec1, p) # wave function
        @unpack Hamiltonian, N, Hs = p
        for i in 1:N
            k = [i - 1, j - 1] * 2pi / N
            eigens = eigen!(Hamiltonian(k))
            psi_1[i, :, :] .= eigens.vectors
            Evec1[i, :] .= eigens.values
        end
    end

    @views function U!(Link0, Link1, LinkN, link10, link01, psi_0, psi_1, psi_N, Evec0, Evec1, psi00, psi10, psi01, Enevec, j, p) # link variable
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
                    l = 1
                    while l <= Hs
                        l0 = Hs - count(Enevec .> (gapless + Enevec[l]))

                        if l == l0
                            link10[l:l0] .= psi00[:, l:l0] ⋅ psi10[:, l:l0]
                            link01[l:l0] .= psi00[:, l:l0] ⋅ psi01[:, l:l0]
                        else
                            link10[l:l0] .= det(psi00[:, l:l0]' * psi10[:, l:l0])
                            link01[l:l0] .= det(psi00[:, l:l0]' * psi01[:, l:l0])
                        end

                        l = 1 + l0
                    end

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
                if i == N && j == N-1
                    psi10[:, :] = psi_0[1, :, :]
                    psi01[:, :] = psi_N[N, :, :]
                elseif i == N
                    psi10[:, :] = psi_0[1, :, :]
                    psi01[:, :] = psi_1[N, :, :]
                elseif j == N-1
                    psi10[:, :] = psi_0[i+1, :, :]
                    psi01[:, :] = psi_N[i, :, :]
                else
                    psi10[:, :] = psi_0[i+1, :, :]
                    psi01[:, :] = psi_1[i, :, :]
                end

                Enevec[:] = Evec1[i, :]
                l = 1
                while l <= Hs
                    l0 = Hs - count(Enevec .> (gapless + Enevec[l]))

                    if l == l0
                        link10[l:l0] .= psi00[:, l:l0] ⋅ psi10[:, l:l0]
                        link01[l:l0] .= psi00[:, l:l0] ⋅ psi01[:, l:l0]
                    else
                        link10[l:l0] .= det(psi00[:, l:l0]' * psi10[:, l:l0])
                        link01[l:l0] .= det(psi00[:, l:l0]' * psi01[:, l:l0])
                    end

                    l = 1 + l0
                end

                Link1[:, :, i] .= [link10 link01]
            end
        end
    end

    @views function F!(phi, dphi, i, j, Link0, Link1, LinkN, p) # lattice field strength
        @unpack N, Hs = p

        if i == N && j == N
            phi[:] = [imag(log(Link0[l, 1, N]*Link0[l, 2, 1]*conj(LinkN[l, 1, N])*conj(Link0[l, 2, N]))) for l in 1:Hs]
            dphi[:] = [imag(log(Link0[l, 1, N]))+imag(log(Link0[l, 2, 1]))-imag(log(LinkN[l, 1, N]))-imag(log(Link0[l, 2, N])) for l in 1:Hs]
        elseif i == N
            phi[:] = [imag(log(Link0[l, 1, N]*Link0[l, 2, 1]*conj(Link1[l, 1, N])*conj(Link0[l, 2, N]))) for l in 1:Hs]
            dphi[:] = [imag(log(Link0[l, 1, N]))+imag(log(Link0[l, 2, 1]))-imag(log(Link1[l, 1, N]))-imag(log(Link0[l, 2, N])) for l in 1:Hs]
        elseif j == N
            phi[:] = [imag(log(Link0[l, 1, i]*Link0[l, 2, i+1]*conj(LinkN[l, 1, i])*conj(Link0[l, 2, i]))) for l in 1:Hs]
            dphi[:] = [imag(log(Link0[l, 1, i]))+imag(log(Link0[l, 2, i+1]))-imag(log(LinkN[l, 1, i]))-imag(log(Link0[l, 2, i])) for l in 1:Hs]
        else
            phi[:] = [imag(log(Link0[l, 1, i]*Link0[l, 2, i+1]*conj(Link1[l, 1, i])*conj(Link0[l, 2, i]))) for l in 1:Hs]
            dphi[:] = [imag(log(Link0[l, 1, i]))+imag(log(Link0[l, 2, i+1]))-imag(log(Link1[l, 1, i]))-imag(log(Link0[l, 2, i])) for l in 1:Hs]
        end

        if rounds == true
            phi[:] = [round(Int, (phi[i] - dphi[i]) / 2pi) for i in 1:Hs]
        else
            phi .= (phi - dphi) / 2pi
        end
    end

    @views function Phase!(TopologicalNumber, p) # chern number
        @unpack N, Hs = p
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

    function main(Hamiltonian, N, gapless, rounds)

        n0 = zeros(2)
        Hs = size(Hamiltonian(n0))[1]
        p = (; Hamiltonian, N, gapless, Hs)

        if rounds == true
            TopologicalNumber = zeros(Int, Hs)
        else
            TopologicalNumber = zeros(Hs)
        end

        Phase!(TopologicalNumber, p)

        Total = sum(TopologicalNumber)

        (; TopologicalNumber, Total)
    end

    main(Hamiltonian, N, gapless, rounds)
end