function Dispersion(Hamiltonian::Function, dim::Int; N::Int=51, labels::Bool=true)
    # GLMakie.activate!(inline=false)

    function Ene1D(nrang, p) # 1D Energy
        @unpack Hamiltonian, N, Hs = p

        Ene = zeros(N, Hs)

        for i in 1:N
            n = nrang[i]
            k = 2pi * n / (N)
            Ene[i, :] = eigvals(Hamiltonian(k))
        end
        Ene
    end

    function Ene2D(nrang, p) # 2D Energy
        @unpack Hamiltonian, N, Hs = p

        Ene = zeros(N, N, Hs)

        for j in 1:N
            n2 = nrang[j]
            for i in 1:N
                n1 = nrang[i]
                n = [n1, n2]
                k = 2pi * n / (N)
                Ene[i, j, :] = eigvals(Hamiltonian(k))
            end
        end
        Ene
    end

    function diagram(p)
        @unpack dim, Hs = p

        nrang = range(0, N, length=N)

        fig = Figure()

        if dim == 1

            if labels == true
                Axis(fig[1, 1], xlabel="k", ylabel="Energy")
            else
                Axis(fig[1, 1])
            end

            Ene = Ene1D(nrang, p)

            for i in 1:Hs
                lines!(nrang, Ene[:, i])
            end
        elseif dim == 2

            if labels == true
                Axis3(fig[1, 1], xlabel="k_1", ylabel="k_2", zlabel="Energy")
            else
                Axis3(fig[1, 1])
            end

            Ene = Ene2D(nrang, p)

            for i in 1:Hs
                surface!(nrang, nrang, Ene[:, :, i])
            end
        end

        # save("Dispersion.png", fig)
        # save("Dispersion.svg", fig)
        # save("Dispersion.pdf", fig)
        fig
    end

    if dim == 1
        n0 = 0.0
        Hs = size(Hamiltonian(n0))[1]
    elseif dim == 2
        n0 = zeros(2)
        Hs = size(Hamiltonian(n0))[1]
    end

    p = (; Hamiltonian, dim, N, Hs)
    diagram(p)
end