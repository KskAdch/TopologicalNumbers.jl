function Ene1D(p) # 1D Energy
    @unpack Hamiltonian, N, Hs = p

    k = range(-π, π, length=N)
    Ene = zeros(N, Hs)

    for i in 1:N
        Ene[i, :] = eigvals!(Hamiltonian(k[i]))
    end
    k, Ene
end

function Ene2D(p) # 2D Energy
    @unpack Hamiltonian, N, Hs = p

    k = range(-π, π, length=N)
    k0 = zeros(2)
    Ene = zeros(N, N, Hs)

    for j in 1:N
        k0[2] = k[j]
        for i in 1:N
            k0[1] = k[i]
            Ene[i, j, :] .= eigvals!(Hamiltonian(k0))
        end
    end
    k = hcat(k, k)
    k, Ene
end

function diagram(p)
    @unpack dim, N, labels, Hs = p

    nrang = range(-π, π, length=N)

    fig = Figure()

    if dim == 1

        if labels == true
            Axis(fig[1, 1], xticks=([-π, 0, π], ["-π", "0", "π"]), xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(2), xlabel="k", ylabel="Eₖ")
        else
            Axis(fig[1, 1], xlabelvisible=false, ylabelvisible=false)
        end

        k, Ene = Ene1D(p)

        for i in 1:Hs
            lines!(nrang, Ene[:, i])
        end
    elseif dim == 2

        if labels == true
            Axis3(fig[1, 1], xticks=([-π, 0, π], ["-π", "0", "π"]), yticks=([-π, 0, π], ["-π", "0", "π"]), xlabel="k₁", ylabel="k₂", zlabel="Eₖ")
        else
            Axis3(fig[1, 1], xlabelvisible=false, ylabelvisible=false, zlabelvisible=false)
        end

        k, Ene = Ene2D(p)

        for i in 1:Hs
            surface!(nrang, nrang, Ene[:, :, i])
        end
    end

    k, Ene, fig
end

function output(k, Ene, fig, p)
    @unpack value, disp, png, pdf, svg = p
    if png == true
        CairoMakie.activate!()
        save("Band.png", fig)
    end
    if pdf == true
        CairoMakie.activate!()
        save("Band.pdf", fig)
    end
    if svg == true
        CairoMakie.activate!()
        save("Band.svg", fig)
    end
    if disp == true
        GLMakie.activate!()
        display(fig)
    end
    if value == true
        return (; k, Ene)
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
function showBand(Hamiltonian::Function; N::Int=51, labels::Bool=true, value::Bool=true, disp::Bool=false, png::Bool=false, pdf::Bool=false, svg::Bool=false)

    dim = Hs = 0

    try
        Hs = size(Hamiltonian(0))[1]
        dim = 1
    catch
        Hs = size(Hamiltonian(zeros(2)))[1]
        dim = 2
    end

    p = (; Hamiltonian, dim, N, labels, Hs, value, disp, png, pdf, svg)
    k, Ene, fig = diagram(p)
    output(k, Ene, fig, p)
end