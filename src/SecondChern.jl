using LinearAlgebra
# using PythonPlot
using LaTeXStrings
using StaticArrays
using Accessors
using Parameters

@with_kw mutable struct Temporal{T}
    P::T
end

⊗(x, y) = kron(x, y)

function H(p)
    k1, k2, k3, k4 = p.sys.k
    b = p.basis

    h1 = p.m + cos(k1) + cos(k2) + cos(k3) + cos(k4)
    h2 = sin(k1)
    h3 = sin(k2)
    h4 = sin(k3)
    h5 = sin(k4)

    return h1 .* b.g1 .+ h2 .* b.g2 .+ h3 .* b.g3 .+ h4 .* b.g4 .+ h5 .* b.g5
end

function linkUx(i, j, l, m, s)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m] * P
    # return adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m]
    return s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i+1, j, l, m] * s.P
end

function linkUx_inv(i, j, l, m, s)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i+1, j, l, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i+1, j, l, m]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m])
    return s.P * inv(adjoint(s.evec[i, j, l, m]) * s.evec[i+1, j, l, m]) * s.P
end

function linkUy(i, j, l, m, s)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m] * P
    # return adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m]
    return s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i, j+1, l, m] * s.P
end

function linkUy_inv(i, j, l, m, s)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j+1, l, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j+1, l, m]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m])
    return s.P * inv(adjoint(s.evec[i, j, l, m]) * s.evec[i, j+1, l, m]) * s.P
end

function linkUz(i, j, l, m, s)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m] * P
    # return adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m]
    return s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l+1, m] * s.P
end

function linkUz_inv(i, j, l, m, s)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l+1, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j, l+1, m]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m])
    return s.P * inv(adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l+1, m]) * s.P
end

function linkUw(i, j, l, m, s)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1] * P
    # return adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1]
    return s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l, m+1] * s.P
end

function linkUw_inv(i, j, l, m, s)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m+1]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j, l, m+1]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1])
    return s.P * inv(adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l, m+1]) * s.P
end

function Fxy(i, j, l, m, s)
    return log(linkUx(i, j, l, m, s) * linkUy(i + 1, j, l, m, s) * linkUx_inv(i, j + 1, l, m, s) * linkUy_inv(i, j, l, m, s))
end
function Fzw(i, j, l, m, s)
    return log(linkUz(i, j, l, m, s) * linkUw(i, j, l + 1, m, s) * linkUz_inv(i, j, l, m + 1, s) * linkUw_inv(i, j, l, m, s))
end
function Fwx(i, j, l, m, s)
    return log(linkUw(i, j, l, m, s) * linkUx(i, j, l, m + 1, s) * linkUw_inv(i + 1, j, l, m, s) * linkUx_inv(i, j, l, m, s))
end
function Fzy(i, j, l, m, s)
    return log(linkUz(i, j, l, m, s) * linkUy(i, j, l + 1, m, s) * linkUz_inv(i, j + 1, l, m, s) * linkUy_inv(i, j, l, m, s))
end
function Fzx(i, j, l, m, s)
    return log(linkUz(i, j, l, m, s) * linkUx(i, j, l + 1, m, s) * linkUz_inv(i + 1, j, l, m, s) * linkUx_inv(i, j, l, m, s))
end
function Fyw(i, j, l, m, s)
    return log(linkUy(i, j, l, m, s) * linkUw(i, j + 1, l, m, s) * linkUy_inv(i, j, l, m + 1, s) * linkUw_inv(i, j, l, m, s))
end

function chernF!(i, j, l, m, p)
    s = p.sys
    evec = s.evec
    # display(linkUz(i, j, l, m, evec) * linkUw(i, j, l + 1, m, evec) * linkUz_inv(i, j, l, m + 1, evec) * linkUw_inv(i, j, l, m, evec))
    # println(i, " ", j, " ", l, " ", m)
    s.P = s.evec[i, j, l, m][:, 1] * adjoint(s.evec[i, j, l, m][:, 1]) + s.evec[i, j, l, m][:, 2] * adjoint(s.evec[i, j, l, m][:, 2])
    s.chern += tr(Fxy(i, j, l, m, s) * Fzw(i, j, l, m, s) + Fwx(i, j, l, m, s) * Fzy(i, j, l, m, s) + Fzx(i, j, l, m, s) * Fyw(i, j, l, m, s))
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

    basis = (; g1, g2, g3, g4, g5)

    Nx = 50
    Ny = 50
    Nz = 50
    Nw = 50

    Kxrange = range(-pi, pi, length=Nx)
    Kyrange = range(-pi, pi, length=Ny)
    Kzrange = range(-pi, pi, length=Nz)
    Kwrange = range(-pi, pi, length=Nw)

    mN = 10
    mList = range(-5.0, 5.0, length=mN)
    # mList = [-3.0]
    ChernList = zeros(ComplexF64, mN)

    r = (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange)


    # k = @MVector [0.0, 0.0, 0.0, 0.0]
    k = [0.0, 0.0, 0.0, 0.0]

    # evec = [@MMatrix zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
    # evec = [zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
    evec = [zeros(ComplexF64, 4, 4) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

    P = zeros(ComplexF64, 4, 4)
    sys = Temporal(; chern=zero(ComplexF64), k, evec, P)

    (; m=zero(Float64), mN, mList, ChernList, basis, r, sys)
end

function setBasis!(p)
    (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange) = p.r
    s = p.sys

    for m in eachindex(Kwrange), l in eachindex(Kzrange), j in eachindex(Kyrange), i in eachindex(Kxrange)
        # k = @SVector [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
        s.k .= [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
        # s.evec[i, j, l, m] = @view(eigvecs(H(k, p))[:, 1:2])
        s.evec[i, j, l, m] .= eigvecs(H(p))
    end

    for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange)
        s.evec[l, j, i, Nz+1] = s.evec[l, j, i, 1]
    end
    for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange)
        s.evec[l, j, Nz+1, i] = s.evec[l, j, 1, i]
    end
    for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange)
        s.evec[l, Ny+1, j, i] = s.evec[l, 1, j, i]
    end
    for i in eachindex(Kxrange), j in eachindex(Kyrange), l in eachindex(Kzrange)
        s.evec[Nz+1, l, j, i] = s.evec[1, l, j, i]
    end

    for i in eachindex(Kxrange), j in eachindex(Kyrange)
        s.evec[j, i, Nz+1, Nw+1] = s.evec[j, i, 1, 1]
    end
    for i in eachindex(Kxrange), j in eachindex(Kyrange)
        s.evec[j, Nz+1, i, Nw+1] = s.evec[j, 1, i, 1]
    end
    for i in eachindex(Kxrange), j in eachindex(Kyrange)
        s.evec[Ny+1, j, i, Nw+1] = s.evec[1, j, i, 1]
    end
    for i in eachindex(Kxrange), j in eachindex(Kyrange)
        s.evec[j, Nz+1, Nw+1, i] = s.evec[j, 1, 1, i]
    end
    for i in eachindex(Kxrange), j in eachindex(Kyrange)
        s.evec[Nz+1, j, Nw+1, i] = s.evec[1, j, 1, i]
    end
    for i in eachindex(Kxrange), j in eachindex(Kyrange)
        s.evec[Ny+1, Nz+1, j, i] = s.evec[1, 1, j, i]
    end

    for i in eachindex(Kxrange)
        s.evec[i, Ny+1, Nz+1, Nw+1] = s.evec[i, 1, 1, 1]
    end
    for i in eachindex(Kxrange)
        s.evec[Nz+1, i, Ny+1, Nw+1] = s.evec[1, i, 1, 1]
    end
    for i in eachindex(Kxrange)
        s.evec[Ny+1, Nz+1, i, Nw+1] = s.evec[1, 1, i, 1]
    end
    for i in eachindex(Kxrange)
        s.evec[Nz+1, Ny+1, Nw+1, i] = s.evec[1, 1, 1, i]
    end

    s.evec[Nx+1, Ny+1, Nz+1, Nw+1] = s.evec[1, 1, 1, 1]
end

function updateChern!(p)
    r = p.r
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange), l in eachindex(r.Kzrange), m in eachindex(r.Kwrange)
        chernF!(i, j, l, m, p)
    end
end

function main()

    p = setParams()

    s = p.sys

    for i in eachindex(p.mList)

        s.chern = 0.0 + 0.0im

        @reset p.m = p.mList[i]

        setBasis!(p)

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

        updateChern!(p)

        s.chern /= 4(pi^2)
        display(p.m)
        display(s.chern)
        p.ChernList[i] = s.chern
    end

    makeFigure(mList, ChernList)

end

@time main()