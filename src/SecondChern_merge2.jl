
using LinearAlgebra
using PythonPlot
using LaTeXStrings
using StaticArrays
using Accessors
using Parameters

@with_kw mutable struct Temporal{T1,T2,T3,T4,T5}
    chern::T1
    k::T2
    evec::T3
    allP::T3
    P::T4
    H::T4
    Fxy::T5
    Fzw::T5
    Fwx::T5
    Fzy::T5
    Fzx::T5
    Fyw::T5
end

⊗(x, y) = kron(x, y)

function H!(p)
    k1, k2, k3, k4 = p.sys.k
    b = p.basis

    h1 = p.m + cos(k1) + cos(k2) + cos(k3) + cos(k4)
    h2 = sin(k1)
    h3 = sin(k2)
    h4 = sin(k3)
    h5 = sin(k4)

    p.sys.H = h1 .* b.g1 .+ h2 .* b.g2 .+ h3 .* b.g3 .+ h4 .* b.g4 .+ h5 .* b.g5
end


function Fxy!(i, j, l, m, s, p)
    y = s.evec[i, j, l, m]' * s.allP[i, j+1, l, m] * s.allP[i+1, j+1, l, m] * s.allP[i+1, j, l, m] * s.evec[i, j, l, m]
    z = y * p.r.Nfill / tr(y)
    s.Fxy = log(z)
end
function Fzw!(i, j, l, m, s, p)
    y = s.evec[i, j, l, m]' * s.allP[i, j, l, m+1] * s.allP[i, j, l+1, m+1] * s.allP[i, j, l+1, m] * s.evec[i, j, l, m]
    z = y * p.r.Nfill / tr(y)
    s.Fzw = log(z)
end
function Fwx!(i, j, l, m, s, p)
    y = s.evec[i, j, l, m]' * s.allP[i+1, j, l, m] * s.allP[i+1, j, l, m+1] * s.allP[i, j, l, m+1] * s.evec[i, j, l, m]
    z = y * p.r.Nfill / tr(y)
    s.Fwx = log(z)
end
function Fzy!(i, j, l, m, s, p)
    y = s.evec[i, j, l, m]' * s.allP[i, j+1, l, m] * s.allP[i, j+1, l+1, m] * s.allP[i, j, l+1, m] * s.evec[i, j, l, m]
    z = y * p.r.Nfill / tr(y)
    s.Fzy = log(z)
end
function Fzx!(i, j, l, m, s, p)
    y = s.evec[i, j, l, m]' * s.allP[i+1, j, l, m] * s.allP[i+1, j, l+1, m] * s.allP[i, j, l+1, m] * s.evec[i, j, l, m]
    z = y * p.r.Nfill / tr(y)
    s.Fzx = log(z)
end
function Fyw!(i, j, l, m, s, p)
    y = s.evec[i, j, l, m]' * s.allP[i, j, l, m+1] * s.allP[i, j+1, l, m+1] * s.allP[i, j+1, l, m] * s.evec[i, j, l, m]
    z = y * p.r.Nfill / tr(y)
    s.Fyw = log(z)
end

function chernF!(i, j, l, m, p)
    s = p.sys

    # println("i = ", i, ", j = ", j, ", l = ", l, ", m = ", m)

    Fxy!(i, j, l, m, s, p)
    Fzw!(i, j, l, m, s, p)
    Fwx!(i, j, l, m, s, p)
    Fzy!(i, j, l, m, s, p)
    Fzx!(i, j, l, m, s, p)
    Fyw!(i, j, l, m, s, p)
    s.Fxy = s.Fxy * s.Fzw + s.Fwx * s.Fzy + s.Fzx * s.Fyw
    s.chern += tr(s.Fxy)
end

function makeFigure(mList, ChernList)

    fig = figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.plot(mList, ChernList, "o--")
    ax.set_xlabel(L"m")
    ax.set_ylabel("Second Chern number")
    ax.grid(true)
    tight_layout()
    savefig("secondChern.png")
    plotshow()
end

function setParams()
    σ₀ = [1 0; 0 1]
    σ₁ = [0 1; 1 0]
    σ₂ = [0 -im; im 0]
    σ₃ = [1 0; 0 -1]
    g1 = σ₁ ⊗ σ₀
    g2 = σ₂ ⊗ σ₀
    g3 = σ₃ ⊗ σ₁
    g4 = σ₃ ⊗ σ₂
    g5 = σ₃ ⊗ σ₃

    basis = (; g1, g2, g3, g4, g5)

    Nx = 30
    Ny = 30
    Nz = 30
    Nw = 30

    Kxrange = range(-pi, pi, length=Nx + 1)
    Kyrange = range(-pi, pi, length=Ny + 1)
    Kzrange = range(-pi, pi, length=Nz + 1)
    Kwrange = range(-pi, pi, length=Nw + 1)

    threshold = 1e-15

    NH = 4
    Nfill = 2

    mN = 10
    mList = range(-5.0, 5.0, length=mN)
    # mList = [-3.0]
    ChernList = zeros(ComplexF64, mN)

    r = (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange, threshold, NH, Nfill)


    k = [0.0, 0.0, 0.0, 0.0]

    evec = [zeros(ComplexF64, NH, Nfill) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

    allP = [zeros(ComplexF64, NH, NH) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

    P = zeros(ComplexF64, NH, NH)
    H = zeros(ComplexF64, NH, NH)
    Fxy = zeros(ComplexF64, Nfill, Nfill)
    Fzw = zeros(ComplexF64, Nfill, Nfill)
    Fwx = zeros(ComplexF64, Nfill, Nfill)
    Fzy = zeros(ComplexF64, Nfill, Nfill)
    Fzx = zeros(ComplexF64, Nfill, Nfill)
    Fyw = zeros(ComplexF64, Nfill, Nfill)

    sys = Temporal(; chern=zero(ComplexF64), k, evec, allP, P, H, Fxy, Fzw, Fwx, Fzy, Fzx, Fyw)

    (; m=zero(Float64), mN, mList, ChernList, basis, r, sys)
end

function setBasis!(p)
    (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange, Nfill) = p.r
    s = p.sys

    for m in eachindex(Kwrange), l in eachindex(Kzrange), j in eachindex(Kyrange), i in eachindex(Kxrange)
        s.k = [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
        H!(p)
        s.H = eigvecs(s.H)
        s.evec[i, j, l, m] = @view s.H[:, 1:Nfill]
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

# Projector
function projector!(i, j, l, m, p)
    s = p.sys
    s.P = s.evec[i, j, l, m] * s.evec[i, j, l, m]'
end

# Compute all projectors
function setProjectors!(p)
    s = p.sys
    r = p.r
    for m in eachindex(r.Kwrange), l in eachindex(r.Kzrange), j in eachindex(r.Kyrange), i in eachindex(r.Kxrange)
        projector!(i, j, l, m, p)
        s.allP[i, j, l, m] = s.P
    end
end

function updateChern!(p)
    r = p.r
    for i in eachindex(r.Kxrange)[1:end-1], j in eachindex(r.Kyrange)[1:end-1], l in eachindex(r.Kzrange)[1:end-1], m in eachindex(r.Kwrange)[1:end-1]
        chernF!(i, j, l, m, p)
    end
end

function threshold_zero!(A::T, p) where {T<:Union{AbstractMatrix,AbstractVector}}
    for i in eachindex(A)
        if abs(A[i]) < p.r.threshold
            A[i] = 0
        elseif abs(real(A[i])) < p.r.threshold
            A[i] = imag(A[i]) * im
        elseif abs(imag(A[i])) < p.r.threshold
            A[i] = real(A[i])
        end
    end
end
function threshold_zero!(A::T, p) where {T<:Number}
    if abs(A) < p.r.threshold
        A = 0.0
    elseif abs(real(A)) < p.r.threshold
        A = imag(A[i]) * im
    elseif abs(imag(A)) < p.r.threshold
        A = real(A)
    end
    return A
end



function main()

    p = setParams()

    s = p.sys

    for i in eachindex(p.mList)

        s.chern = 0.0 + 0.0im

        @reset p.m = p.mList[i]

        setBasis!(p)

        setProjectors!(p)

        updateChern!(p)

        s.chern /= 4pi^2
        println("m: ", p.m, ", Second Chern number: ", s.chern)

        p.ChernList[i] = s.chern
    end

    # makeFigure(p.mList, p.ChernList)

end

@time main()