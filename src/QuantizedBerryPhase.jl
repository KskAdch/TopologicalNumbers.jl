function QuantizedBerryPhase(Hamiltonian::Function; N::Int=51, gapless::Float64=0.0)

    function psi!(i, psi1, Evec1, p) # wave function
        (; Hamiltonian, N, Hs) = p

        nrang = range(0, N - 1, length=N)

        n = nrang[i]
        k = 2pi * n / N
        eigens = eigen(Hamiltonian(k))
        psi1[:, :] .= eigens.vectors
        Evec1[:] .= eigens.values
    end

    @views function U!(Link, Evec0, Evec1, psi0, psi1, psiN1, i, p) # link variable
        (; N, gapless, Hs) = p

        if i == 1
            psi!(i, psi1, Evec1, p)
            psiN1 .= psi1
        end

        if i != N
            psi0 .= psi1
            Evec0 .= Evec1
            psi!(i + 1, psi1, Evec1, p)
        end

        Evec0[:]
        l = 1
        while l <= Hs
            Enevec = Evec0[Evec0.>Evec0[l]] .- Evec0[l]
            l0 = l + count(Enevec .<= gapless)

            if i == N
                Link[l:l0] .= det(psi1[:, l:l0]' * psiN1[:, l:l0])
            else
                Link[l:l0] .= det(psi0[:, l:l0]' * psi1[:, l:l0])
            end

            l = 1 + l0
        end
    end

    @views function F!(phi, i, Link, p) # lattice field strength
        (; N, Hs) = p

        phi[:] = [imag(log(Link[l])) for l in 1:Hs]
    end

    @views function Phase!(TopologicalNumber, p) # berry phase
        (; N, Hs) = p
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
            F!(phi, i, Link, p)
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