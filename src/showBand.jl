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
            Axis(fig[1, 1], xticks = ([0, pi, 2pi], ["0", "π", "2π"]), xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(2), xlabel="k", ylabel="Eₖ")
        else
            Axis(fig[1, 1], xlabelvisible = false, ylabelvisible = false)
        end

        Ene = Ene1D(p)

        for i in 1:Hs
            lines!(nrang, Ene[:, i])
        end
    elseif dim == 2

        if labels == true
            Axis3(fig[1, 1], xticks = ([0, pi, 2pi], ["0", "π", "2π"]), yticks = ([0, pi, 2pi], ["0", "π", "2π"]), xlabel="k₁", ylabel="k₂", zlabel="Eₖ")
        else
            Axis3(fig[1, 1], xlabelvisible = false, ylabelvisible = false, zlabelvisible = false)
        end

        Ene = Ene2D(p)

        for i in 1:Hs
            surface!(nrang, nrang, Ene[:, :, i])
        end
    end
    
    Ene, fig
end

function output(Ene, fig, p)
    @unpack value, png, pdf, svg, fig3d = p
    if png == true
        save("Band.png", fig)
    end
    if pdf == true
        save("Band.pdf", fig)
    end
    if svg == true
        save("Band.svg", fig)
    end
    if fig3d == true
        GLMakie.activate!()
        display(fig)
    end
    if value == true
        return Ene
    end
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
function showBand(Hamiltonian::Function; N::Int=51, labels::Bool=true, value::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, fig3d::Bool=false)

    dim = Hs = 0

    try
        Hs = size(Hamiltonian(0))[1]
        dim = 1
    catch
        Hs = size(Hamiltonian(zeros(2)))[1]
        dim = 2
    end

    p = (; Hamiltonian, dim, N, labels, Hs, value, png, pdf, svg, fig3d)
    Ene, fig = diagram(p)
    output(Ene, fig, p)
end