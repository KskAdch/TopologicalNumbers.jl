function QuantizedBerryPhase(Hamiltonian::Function; N::Int=51, gapless::Real=0.0)

    function psi!(i, psi1, Evec1, p) # wave function
        @unpack Hamiltonian, N, Hs = p

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
                    Link[l:l0] .= dot(psi1[:, l:l0],  psiN1[:, l:l0])
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

    @views function F!(phi, Link, p) # lattice field strength
        @unpack N, Hs = p

        phi[:] = [imag(log(Link[l])) for l in 1:Hs]
    end

    @views function Phase!(TopologicalNumber, p) # berry phase
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

        TopologicalNumber .= [round(Int, abs(rem(TN[i] / pi, 2))) for i in 1:Hs]
    end

    function main(Hamiltonian, N, gapless)

        Hs = size(Hamiltonian(0.0))[1]
        p = (; Hamiltonian, N, gapless, Hs)

        TopologicalNumber = zeros(Int, Hs)

        Phase!(TopologicalNumber, p)

        Total = rem(sum(TopologicalNumber), 2)

        (; TopologicalNumber, Total)
    end

    main(Hamiltonian, N, gapless)
end