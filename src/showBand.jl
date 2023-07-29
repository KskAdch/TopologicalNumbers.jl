function Ene1D(p) # 1D Energy
    @unpack Hamiltonian, N, Hs = p

    Ene = zeros(N, Hs)

    for i in 1:N
        k = 2pi * (i - 1) / (N - 1)
        Ene[i, :] = eigvals!(Hamiltonian(k))
    end
    Ene
end

function Ene2D(p) # 2D Energy
    @unpack Hamiltonian, N, Hs = p

    Ene = zeros(N, N, Hs)

    for j in 1:N
        for i in 1:N
            k = 2pi * [i - 1, j - 1] / (N - 1)
            Ene[i, j, :] = eigvals!(Hamiltonian(k))
        end
    end
    Ene
end

function diagram(p)
    @unpack dim, N, labels, Hs = p

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
            Axis3(fig[1, 1], xlabel="", ylabel="", zlabel="")
        end

        Ene = Ene2D(p)

        for i in 1:Hs
            surface!(nrang, nrang, Ene[:, :, i])
        end
    end

    # save("Band.png", fig)
    # save("Band.svg", fig)
    # save("Band.pdf", fig)
    fig
end

@doc raw"""

 Drawing the band structure of the Hamiltonian.

    showBand(Hamiltonian::Function; N::Int=51, labels::Bool=true)

 Arguments
 - `Hamiltonian::Function`: the Hamiltonian matrix function of wave number $\bm k$.
 - `N::Int`: the number of divisions in the wave number space.
 - `labels::Bool`: whether to display the labels of the figure.

```math
```
"""
function showBand(Hamiltonian::Function; N::Int=51, labels::Bool=true)
    # GLMakie.activate!(inline=false)

    dim = Hs = 0

    try
        Hs = size(Hamiltonian(0.0))[1]
        dim = 1
    catch
        Hs = size(Hamiltonian(zeros(2)))[1]
        dim = 2
    end

    p = (; Hamiltonian, dim, N, labels, Hs)
    diagram(p)
end