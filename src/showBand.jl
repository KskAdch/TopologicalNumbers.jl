function Ene1D(p::Params) # 1D Energy
    @unpack Hamiltonian, N, Hs = p

    k = range(-π, π, length=N)
    Ene = zeros(N, Hs)

    for i in 1:N
        Ene[i, :] = eigvals!(Hamiltonian(k[i]))
    end
    k, Ene
end

function Ene2D(p::Params) # 2D Energy
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

function Ene4D(p::Params) # 4D Energy
    @unpack Hamiltonian, N, Hs = p

    k = range(-π, π, length=N)
    k0 = zeros(4)
    Ene = zeros(N, N, N, N, Hs)

    for m in 1:N
        k0[4] = k[m]
        for l in 1:N
            k0[3] = k[l]
            for j in 1:N
                k0[2] = k[j]
                for i in 1:N
                    k0[1] = k[i]
                    Ene[i, j, l, m, :] .= eigvals!(Hamiltonian(k0))
                end
            end
        end
    end
    k = hcat(k, k)
    k = hcat(k, k)
    k = hcat(k, k)
    k, Ene
end

function diagram(p::Params, p_out)
    @unpack dim, N, Hs = p
    @unpack labels, disp, png, pdf, svg = p_out

    nrang = range(-π, π, length=N)

    fig = figure()

    if dim == 1

        ax = fig.add_subplot(111)

        if labels == true
            ax.set_xlabel(L"k")
            ax.set_ylabel(L"E_k")
            ax.set_xticks([-π, 0, π])
            ax.set_xticklabels([L"-\pi", L"0", L"\pi"])
        end
        ax.grid()

        k, Ene = Ene1D(p)

        if disp == true || png == true || pdf == true || svg == true
            for i in 1:Hs
                ax.plot(nrang, Ene[:, i])
            end
        else
            plotclose()
        end
    elseif dim == 2

        ax = fig.add_subplot(111, projection="3d")

        if labels == true
            ax.set_xlabel(L"k_1")
            ax.set_ylabel(L"k_2")
            ax.set_zlabel(L"E_k")
            ax.set_xticks([-π, 0, π])
            ax.set_xticklabels([L"-\pi", L"0", L"\pi"])
            ax.set_yticks([-π, 0, π])
            ax.set_yticklabels([L"-\pi", L"0", L"\pi"])
        end
        ax.grid()

        k, Ene = Ene2D(p)

        if disp == true || png == true || pdf == true || svg == true
            X, Y = complex.(nrang', nrang) |> z -> (real.(z), imag.(z))

            for i in 1:Hs
                ax.plot_surface(X, Y, Ene[:, :, i], shade=true, antialiased=false)
            end
        else
            plotclose()
        end
    elseif dim == 4

        k, Ene = Ene4D(p)
        plotclose()

    end

    return k, Ene, fig
end

function output(k, Ene, fig, p)
    @unpack value, disp, png, pdf, svg, filename = p
    if png == true
        savefig(filename * ".png")
    end
    if pdf == true
        savefig(filename * ".pdf")
    end
    if svg == true
        savefig(filename * ".svg")
    end
    if disp == true
        # plotshow()
        display(fig)
    end
    if value == true
        return (; k, Ene)
    end
end

@doc raw"""

 Drawing the band structure of the Hamiltonian.

    showBand(Hamiltonian::Function; N::Int=51, labels::Bool=true, value::Bool=true, disp::Bool=false, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="Band")

 Arguments
 - `Hamiltonian::Function`: the Hamiltonian matrix function of wave number $\bm k$.
 - `N::Int`: the number of divisions in the wave number space.
 - `labels::Bool`: whether to display the labels of the figure.
 - `value::Bool`: whether to output the values of the wave number and the energy in the matrix form.
 - `disp::Bool`: whether to display the figure.
 - `png::Bool`: whether to save the figure as a PNG file.
 - `pdf::Bool`: whether to save the figure as a PDF file.
 - `svg::Bool`: whether to save the figure as a SVG file.
 - `filename::String`: the name of the output file.

```math
```
"""
function showBand(Hamiltonian::Function; N::Int=51, labels::Bool=true, value::Bool=true, disp::Bool=false, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="Band")

    dim = Hs = 0

    try
        Hs = size(Hamiltonian(0.0), 1)
        dim = 1
    catch
        try
            Hs = size(Hamiltonian(zeros(2)), 1)
            dim = 2
        catch
            try
                Hs = size(Hamiltonian(zeros(3)), 1)
                dim = 3
            catch
                Hs = size(Hamiltonian(zeros(4)), 1)
                dim = 4
            end
        end
    end

    p = Params(; Hamiltonian, dim, N, Hs, gapless=0.0, rounds=false)
    p_out = (; labels, value, disp, png, pdf, svg, filename)

    k, Ene, fig = diagram(p, p_out)

    output(k, Ene, fig, p_out)
end