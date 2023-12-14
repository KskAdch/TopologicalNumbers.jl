using LinearAlgebra
# using PythonPlot
using LaTeXStrings
using StaticArrays
using Accessors

⊗(x, y) = kron(x, y)

function H(k, p)
    k1, k2, k3, k4 = k

    h1 = p.m + cos(k1) + cos(k2) + cos(k3) + cos(k4)
    h2 = sin(k1)
    h3 = sin(k2)
    h4 = sin(k3)
    h5 = sin(k4)

    return h1 .* p.g1 .+ h2 .* p.g2 .+ h3 .* p.g3 .+ h4 .* p.g4 .+ h5 .* p.g5
end

function linkUx(i, j, l, m, evec)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m] * P
    return adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m]
end

function linkUx_inv(i, j, l, m, evec)
    P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i+1, j, l, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i+1, j, l, m]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m])
    return P * inv(adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m]) * P
end

function linkUy(i, j, l, m, evec)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m] * P
    return adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m]
end

function linkUy_inv(i, j, l, m, evec)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j+1, l, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j+1, l, m]) * evec[i, j, l, m]
    return inv(adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m])
end

function linkUz(i, j, l, m, evec)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m] * P
    return adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m]
end

function linkUz_inv(i, j, l, m, evec)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l+1, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j, l+1, m]) * evec[i, j, l, m]
    return inv(adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m])
end

function linkUw(i, j, l, m, evec)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1] * P
    return adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1]
end

function linkUw_inv(i, j, l, m, evec)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m+1]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j, l, m+1]) * evec[i, j, l, m]
    return inv(adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1])
end

function Fxy(i, j, l, m, evec)
    return log(linkUx(i, j, l, m, evec) * linkUy(i + 1, j, l, m, evec) * linkUx_inv(i, j + 1, l, m, evec) * linkUy_inv(i, j, l, m, evec))
end
function Fzw(i, j, l, m, evec)
    return log(linkUz(i, j, l, m, evec) * linkUw(i, j, l + 1, m, evec) * linkUz_inv(i, j, l, m + 1, evec) * linkUw_inv(i, j, l, m, evec))
end
function Fwx(i, j, l, m, evec)
    return log(linkUw(i, j, l, m, evec) * linkUx(i, j, l, m + 1, evec) * linkUw_inv(i + 1, j, l, m, evec) * linkUx_inv(i, j, l, m, evec))
end
function Fzy(i, j, l, m, evec)
    return log(linkUz(i, j, l, m, evec) * linkUy(i, j, l + 1, m, evec) * linkUz_inv(i, j + 1, l, m, evec) * linkUy_inv(i, j, l, m, evec))
end
function Fzx(i, j, l, m, evec)
    return log(linkUz(i, j, l, m, evec) * linkUx(i, j, l + 1, m, evec) * linkUz_inv(i + 1, j, l, m, evec) * linkUx_inv(i, j, l, m, evec))
end
function Fyw(i, j, l, m, evec)
    return log(linkUy(i, j, l, m, evec) * linkUw(i, j + 1, l, m, evec) * linkUy_inv(i, j, l, m + 1, evec) * linkUw_inv(i, j, l, m, evec))
end

function chernF(i, j, l, m, evec)
    # display(linkUz(i, j, l, m, evec) * linkUw(i, j, l + 1, m, evec) * linkUz_inv(i, j, l, m + 1, evec) * linkUw_inv(i, j, l, m, evec))
    # println(i, " ", j, " ", l, " ", m)
    tr(Fxy(i, j, l, m, evec) * Fzw(i, j, l, m, evec) + Fwx(i, j, l, m, evec) * Fzy(i, j, l, m, evec) + Fzx(i, j, l, m, evec) * Fyw(i, j, l, m, evec))
end

function makeFigure(mList, ChernList)

    fig = figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.plot(mList, ChernList)
    ax.set_xlabel(L"m")
    ax.set_ylabel("Second Chern number")
    ax.grid(true)
    savefig("secondChern.png")
    plotshow()
end

function setParams()
    σ₀ = [1 0; 0 1]
    σ₁ = [0 1; 1 0]
    σ₂ = [0 -im; im 0]
    σ₃ = [1 0; 0 -1]
    # σ₀ = @MMatrix [1 0; 0 1]
    # σ₁ = @MMatrix [0 1; 1 0]
    # σ₂ = @MMatrix [0 -im; im 0]
    # σ₃ = @MMatrix [1 0; 0 -1]
    g1 = σ₁ ⊗ σ₀
    g2 = σ₂ ⊗ σ₀
    g3 = σ₂ ⊗ σ₁
    g4 = σ₂ ⊗ σ₂
    g5 = σ₃ ⊗ σ₃
    (; g1, g2, g3, g4, g5, m=0.0)
end

function main()
    Nx = 50
    Ny = 50
    Nz = 50
    Nw = 50

    Kxrange = range(-pi, pi, length=Nx)
    Kyrange = range(-pi, pi, length=Ny)
    Kzrange = range(-pi, pi, length=Nz)
    Kwrange = range(-pi, pi, length=Nw)
    r = (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange)

    mN = 10
    mList = range(-5.0, 5.0, length=mN)
    # mList = [-3.0]
    ChernList = zeros(ComplexF64, mN)

    p = setParams()
    # k = @MVector [0.0, 0.0, 0.0, 0.0]
    k = [0.0, 0.0, 0.0, 0.0]

    # evec = [@SMatrix zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

    for i in eachindex(mList)

        # evec = [@MMatrix zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
        evec = [zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
        # evec = [zeros(ComplexF64, 4, 4) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

        Chern = 0.0 + 0.0im

        @reset p.m = mList[i]

        for m in eachindex(Kwrange), l in eachindex(Kzrange), j in eachindex(Kyrange), i in eachindex(Kxrange)
            # k = @SVector [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
            k = [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
            evec[i, j, l, m] = @view(eigvecs(H(k, p))[:, 1:2])
            # evec[i, j, l, m] = eigvecs(H(k, p))
        end

        for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange)
            evec[l, j, i, Nz+1] = evec[l, j, i, 1]
        end
        for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange)
            evec[l, j, Nz+1, i] = evec[l, j, 1, i]
        end
        for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange)
            evec[l, Ny+1, j, i] = evec[l, 1, j, i]
        end
        for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange)
            evec[Nz+1, l, j, i] = evec[1, l, j, i]
        end

        for i in eachindex(Kxrange), j in eachindex(Kyrange)
            evec[j, i, Nz+1, Nw+1] = evec[j, i, 1, 1]
        end
        for i in eachindex(Kxrange), j in eachindex(Kyrange)
            evec[j, Nz+1, i, Nw+1] = evec[j, 1, i, 1]
        end
        for i in eachindex(Kxrange), j in eachindex(Kyrange)
            evec[Ny+1, j, i, Nw+1] = evec[1, j, i, 1]
        end
        for i in eachindex(Kxrange), j in eachindex(Kyrange)
            evec[j, Nz+1, Nw+1, i] = evec[j, 1, 1, i]
        end
        for i in eachindex(Kxrange), j in eachindex(Kyrange)
            evec[Nz+1, j, Nw+1, i] = evec[1, j, 1, i]
        end
        for i in eachindex(Kxrange), j in eachindex(Kyrange)
            evec[Ny+1, Nz+1, j, i] = evec[1, 1, j, i]
        end

        for i in eachindex(Kxrange)
            evec[i, Ny+1, Nz+1, Nw+1] = evec[i, 1, 1, 1]
        end
        for i in eachindex(Kxrange)
            evec[Nz+1, i, Ny+1, Nw+1] = evec[1, i, 1, 1]
        end
        for i in eachindex(Kxrange)
            evec[Ny+1, Nz+1, i, Nw+1] = evec[1, 1, i, 1]
        end
        for i in eachindex(Kxrange)
            evec[Nz+1, Ny+1, Nw+1, i] = evec[1, 1, 1, i]
        end

        evec[Nx+1, Ny+1, Nz+1, Nw+1] = evec[1, 1, 1, 1]

        # for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1
        #     if evec[i, j, l, m] == zeros(ComplexF64, 4, 2)
        #         println(i, " ", j, " ", l, " ", m)
        #     end
        # end
        # k = [Kxrange[1], Kyrange[1], Kzrange[1], Kwrange[1]]
        # display(evec[1, 1, 1, 1]' * H(k, p) * evec[1, 1, 1, 1])# == diagm(eigvals(H(k, p))))
        # display(diagm(eigvals(H(k, p))))
        # display(evec[1, 1, 1, 1]' * H(k, p) * evec[1, 1, 1, 1] ≈ diagm(eigvals(H(k, p))))
        # P = evec[1, 1, 1, 1][:, 1] * adjoint(evec[1, 1, 1, 1][:, 1]) + evec[1, 1, 1, 1][:, 2] * adjoint(evec[1, 1, 1, 1][:, 2])
        # display(P)
        # display(P^2)
        # display(P ≈ P^2)
        # display(H(k, p))
        # display(P * H(k, p) * P)
        # display(eigvals(H(k, p)))
        # display(eigvals(P * H(k, p) * P))
        # display(adjoint(evec[1, 1, 1, 1][:, 1]) * evec[1, 1, 1, 1][:, 1])
        # display(adjoint(evec[1, 1, 1, 1][:, 1]) * H(k, p) * evec[1, 1, 1, 1][:, 1])

        for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange), m in eachindex(Kwrange)
            Chern += chernF(i, j, l, m, evec)
        end

        Chern /= 4(pi^2)
        display(p.m)
        display(Chern)
        ChernList[i] = Chern
    end

    # makeFigure(mList, ChernList)

end

@time main()