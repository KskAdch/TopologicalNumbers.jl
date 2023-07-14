function Dispersion(Hamiltonian::Function, dim::Int; N::Int=51, labels::Bool=true)
    # GLMakie.activate!(inline=false)

    function Ene1D(p) # 1D Energy
        @unpack Hamiltonian, N, Hs = p

        Ene = zeros(N, Hs)

        for i in 1:N
            k = 2pi*(i - 1) / (N - 1)
            Ene[i, :] = eigvals!(Hamiltonian(k))
        end
        Ene
    end

    function Ene2D(p) # 2D Energy
    @unpack Hamiltonian, N, Hs = p

        Ene = zeros(N, N, Hs)

        for j in 1:N
            for i in 1:N
                k = 2pi*[i - 1, j - 1] / (N - 1)
                Ene[i, j, :] = eigvals!(Hamiltonian(k))
            end
        end
        Ene
    end

    function diagram(p)
        @unpack dim, Hs = p

        nrang = range(0, 2pi, length=N)

        fig = Figure()

        if dim == 1

            if labels == true
                Axis(fig[1, 1], xlabel="k", ylabel="Eₖ")
            else
                Axis(fig[1, 1], xlabel="", ylabel="")
            end

            Ene = Ene1D(p)

            for i in 1:Hs
                lines!(nrang, Ene[:, i])
            end
        elseif dim == 2

            if labels == true
                Axis3(fig[1, 1], xlabel="k₁", ylabel="k₂", zlabel="Eₖ")
            else
                Axis3(fig[1, 1], xlabel = "", ylabel = "", zlabel = "")
            end

            Ene = Ene2D(p)

            for i in 1:Hs
                surface!(nrang, nrang, Ene[:, :, i])
            end
        end

        save("Dispersion.png", fig)
        save("Dispersion.svg", fig)
        save("Dispersion.pdf", fig)
        fig
    end

    if dim == 1
        Hs = size(Hamiltonian(0.0))[1]
    elseif dim == 2
        Hs = size(Hamiltonian(zeros(2)))[1]
    end

    p = (; Hamiltonian, dim, N, Hs)
    diagram(p)
end