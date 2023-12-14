using LinearAlgebra
using PythonPlot
using LaTeXStrings
using StaticArrays
using Accessors
using Parameters

@with_kw mutable struct Temporal{T1,T2,T3,T4}
    chern::T1
    k::T2
    evec::T3
    P::T4
    H::T4
    Ux::T4
    Uy::T4
    Uz::T4
    Uw::T4
    Ux_inv::T4
    Uy_inv::T4
    Uz_inv::T4
    Uw_inv::T4
    Fxy::T4
    Fzw::T4
    Fwx::T4
    Fzy::T4
    Fzx::T4
    Fyw::T4
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

    p.sys.H = h1 .* b.g1 .+ h2 .* b.g2 .+ h3 .* b.g3 .+ h4 .* b.g4 .+ h5 .* b.g5
    threshold_zero!(p.sys.H, p)
end

function linkUx!(i, j, l, m, s, p)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m] * P
    # return adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m]
    s.Ux = s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i+1, j, l, m] * s.P
    threshold_zero!(s.Ux, p)
end

function linkUx_inv!(i, j, l, m, s, p)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i+1, j, l, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i+1, j, l, m]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i+1, j, l, m])
    s.Ux_inv = s.P * inv(adjoint(s.evec[i, j, l, m]) * s.evec[i+1, j, l, m]) * s.P
    threshold_zero!(s.Ux_inv, p)
end

function linkUy!(i, j, l, m, s, p)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m] * P
    # return adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m]
    s.Uy = s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i, j+1, l, m] * s.P
    threshold_zero!(s.Uy, p)
end

function linkUy_inv!(i, j, l, m, s, p)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j+1, l, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j+1, l, m]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i, j+1, l, m])
    s.Uy_inv = s.P * inv(adjoint(s.evec[i, j, l, m]) * s.evec[i, j+1, l, m]) * s.P
    threshold_zero!(s.Uy_inv, p)
end

function linkUz!(i, j, l, m, s, p)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m] * P
    # return adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m]
    s.Uz = s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l+1, m] * s.P
    threshold_zero!(s.Uz, p)
end

function linkUz_inv!(i, j, l, m, s, p)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l+1, m]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j, l+1, m]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i, j, l+1, m])
    s.Uz_inv = s.P * inv(adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l+1, m]) * s.P
    threshold_zero!(s.Uz_inv, p)
end

function linkUw!(i, j, l, m, s, p)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1] * P
    # return adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1]
    s.Uw = s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l, m+1] * s.P
    threshold_zero!(s.Uw, p)
end

function linkUw_inv!(i, j, l, m, s, p)
    # P = evec[i, j, l, m][:, 1] * adjoint(evec[i, j, l, m][:, 1]) + evec[i, j, l, m][:, 2] * adjoint(evec[i, j, l, m][:, 2])
    # return P * adjoint(evec[i, j, l, m+1]) * evec[i, j, l, m] * P
    # return adjoint(evec[i, j, l, m+1]) * evec[i, j, l, m]
    # return inv(adjoint(evec[i, j, l, m]) * evec[i, j, l, m+1])
    s.Uw_inv = s.P * inv(adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l, m+1]) * s.P
    threshold_zero!(s.Uw_inv, p)
end

function Fxy!(i, j, l, m, s, p)
    linkUx!(i, j, l, m, s, p)
    linkUy!(i + 1, j, l, m, s, p)
    linkUx_inv!(i, j + 1, l, m, s, p)
    linkUy_inv!(i, j, l, m, s, p)
    s.Fxy = log(s.Ux * s.Uy * s.Ux_inv * s.Uy_inv)
    threshold_zero!(s.Fxy, p)
end
function Fzw!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s, p)
    linkUw!(i, j, l + 1, m, s, p)
    linkUz_inv!(i, j, l, m + 1, s, p)
    linkUw_inv!(i, j, l, m, s, p)
    s.Fzw = log(s.Uz * s.Uw * s.Uz_inv * s.Uw_inv)
    threshold_zero!(s.Fzw, p)
end
function Fwx!(i, j, l, m, s, p)
    linkUw!(i, j, l, m, s, p)
    linkUx!(i, j, l, m + 1, s, p)
    linkUw_inv!(i + 1, j, l, m, s, p)
    linkUx_inv!(i, j, l, m, s, p)
    s.Fwx = log(s.Uw * s.Ux * s.Uw_inv * s.Ux_inv)
    threshold_zero!(s.Fwx, p)
end
function Fzy!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s, p)
    linkUy!(i, j, l + 1, m, s, p)
    linkUz_inv!(i, j + 1, l, m, s, p)
    linkUy_inv!(i, j, l, m, s, p)
    s.Fzy = log(s.Uz * s.Uy * s.Uz_inv * s.Uy_inv)
    threshold_zero!(s.Fzy, p)
end
function Fzx!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s, p)
    linkUx!(i, j, l + 1, m, s, p)
    linkUz_inv!(i + 1, j, l, m, s, p)
    linkUx_inv!(i, j, l, m, s, p)
    s.Fzx = log(s.Uz * s.Ux * s.Uz_inv * s.Ux_inv)
    threshold_zero!(s.Fzx, p)
end
function Fyw!(i, j, l, m, s, p)
    linkUy!(i, j, l, m, s, p)
    linkUw!(i, j + 1, l, m, s, p)
    linkUy_inv!(i, j, l, m + 1, s, p)
    linkUw_inv!(i, j, l, m, s, p)
    s.Fyw = log(s.Uy * s.Uw * s.Uy_inv * s.Uw_inv)
    threshold_zero!(s.Fyw, p)
end

function chernF!(i, j, l, m, p)
    s = p.sys
    # display(linkUz(i, j, l, m, evec) * linkUw(i, j, l + 1, m, evec) * linkUz_inv(i, j, l, m + 1, evec) * linkUw_inv(i, j, l, m, evec))
    # println(i, " ", j, " ", l, " ", m)
    s.P = s.evec[i, j, l, m][:, 1] * adjoint(s.evec[i, j, l, m][:, 1]) + s.evec[i, j, l, m][:, 2] * adjoint(s.evec[i, j, l, m][:, 2])
    threshold_zero!(s.P, p)
    Fxy!(i, j, l, m, s, p)
    Fzw!(i, j, l, m, s, p)
    Fwx!(i, j, l, m, s, p)
    Fzy!(i, j, l, m, s, p)
    Fzx!(i, j, l, m, s, p)
    Fyw!(i, j, l, m, s, p)
    s.Fxy = s.Fxy * s.Fzw + s.Fwx * s.Fzy + s.Fzx * s.Fyw
    threshold_zero!(s.Fxy, p)
    s.chern += tr(s.Fxy)
end

function makeFigure(mList, ChernList)

    fig = figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.plot(mList, ChernList)
    ax.set_xlabel(L"m")
    ax.set_ylabel("Second Chern number")
    ax.grid(true)
    tight_layout()
    savefig("secondChern.png")
    plotshow()
end

function setParams()
    # σ₀ = [1 0; 0 1]
    # σ₁ = [0 1; 1 0]
    # σ₂ = [0 -im; im 0]
    # σ₃ = [1 0; 0 -1]
    σ₀ = @MMatrix [1 0; 0 1]
    σ₁ = @MMatrix [0 1; 1 0]
    σ₂ = @MMatrix [0 -im; im 0]
    σ₃ = @MMatrix [1 0; 0 -1]
    g1 = σ₁ ⊗ σ₀
    g2 = σ₂ ⊗ σ₀
    g3 = σ₂ ⊗ σ₁
    g4 = σ₂ ⊗ σ₂
    g5 = σ₃ ⊗ σ₃

    basis = (; g1, g2, g3, g4, g5)

    Nx = 10
    Ny = 10
    Nz = 10
    Nw = 10

    Kxrange = range(-pi, pi, length=Nx)
    Kyrange = range(-pi, pi, length=Ny)
    Kzrange = range(-pi, pi, length=Nz)
    Kwrange = range(-pi, pi, length=Nw)

    threshold = 1e-12

    mN = 1 # 10
    # mList = range(-5.0, 5.0, length=mN)
    mList = [-3.0]
    ChernList = zeros(ComplexF64, mN)

    r = (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange, threshold)


    k = @MVector [0.0, 0.0, 0.0, 0.0]
    # k = [0.0, 0.0, 0.0, 0.0]

    # evec = [@MMatrix zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
    # evec = [zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
    # evec = [zeros(ComplexF64, 4, 4) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
    evec = [@MMatrix zeros(ComplexF64, 4, 4) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

    # P = zeros(ComplexF64, 4, 4)
    P = @MMatrix zeros(ComplexF64, 4, 4)
    H = @MMatrix zeros(ComplexF64, 4, 4)
    Ux = @MMatrix zeros(ComplexF64, 4, 4)
    Uy = @MMatrix zeros(ComplexF64, 4, 4)
    Uz = @MMatrix zeros(ComplexF64, 4, 4)
    Uw = @MMatrix zeros(ComplexF64, 4, 4)
    Ux_inv = @MMatrix zeros(ComplexF64, 4, 4)
    Uy_inv = @MMatrix zeros(ComplexF64, 4, 4)
    Uz_inv = @MMatrix zeros(ComplexF64, 4, 4)
    Uw_inv = @MMatrix zeros(ComplexF64, 4, 4)
    Fxy = @MMatrix zeros(ComplexF64, 4, 4)
    Fzw = @MMatrix zeros(ComplexF64, 4, 4)
    Fwx = @MMatrix zeros(ComplexF64, 4, 4)
    Fzy = @MMatrix zeros(ComplexF64, 4, 4)
    Fzx = @MMatrix zeros(ComplexF64, 4, 4)
    Fyw = @MMatrix zeros(ComplexF64, 4, 4)

    sys = Temporal(; chern=zero(ComplexF64), k, evec, P, H, Ux, Uy, Uz, Uw, Ux_inv, Uy_inv, Uz_inv, Uw_inv, Fxy, Fzw, Fwx, Fzy, Fzx, Fyw)

    (; m=zero(Float64), mN, mList, ChernList, basis, r, sys)
end

function setBasis!(p)
    (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange) = p.r
    s = p.sys

    for m in eachindex(Kwrange), l in eachindex(Kzrange), j in eachindex(Kyrange), i in eachindex(Kxrange)
        s.k = @MVector [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
        # s.k = [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
        # s.evec[i, j, l, m] = @view(eigvecs(H(k, p))[:, 1:2])
        s.evec[i, j, l, m] = eigvecs(H(p))
        threshold_zero!(s.evec[i, j, l, m], p)
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

function test(p)
    s = p.sys

    for i in 1:p.r.Nx+1, j in 1:p.r.Ny+1, l in 1:p.r.Nz+1, m in 1:p.r.Nw+1
        if s.evec[i, j, l, m] == zeros(ComplexF64, 4, 2)
            println(i, " ", j, " ", l, " ", m)
        end
    end
    s.k = @MVector [-pi, -pi, -pi, -pi]

    println("[1] Check for eigen values")
    A = s.evec[1, 1, 1, 1]' * H(p) * s.evec[1, 1, 1, 1]
    threshold_zero!(A, p)
    display(A)# == diagm(eigvals(H(k, p))))
    # display(typeof(H(p)))
    # display(typeof(eigvals(H(p))))
    B = MMatrix(diagm(eigvals(H(p))))
    # display(typeof(B))
    threshold_zero!(B, p)
    display(B)
    C = adjoint(s.evec[1, 1, 1, 1][:, 1]) * H(p) * s.evec[1, 1, 1, 1][:, 1]
    C = threshold_zero!(C, p)
    display(C)
    display(A ≈ B)

    println("[2] Check for projection operator")
    # s.P = (
    #     s.evec[1, 1, 1, 1][:, 1] * adjoint(s.evec[1, 1, 1, 1][:, 1]) + s.evec[1, 1, 1, 1][:, 2] * adjoint(s.evec[1, 1, 1, 1][:, 2])
    #     + s.evec[2, 1, 1, 1][:, 1] * adjoint(s.evec[2, 1, 1, 1][:, 1]) + s.evec[2, 1, 1, 1][:, 2] * adjoint(s.evec[2, 1, 1, 1][:, 2])
    # )
    s.P = s.evec[1, 1, 1, 1][:, 1] * adjoint(s.evec[1, 1, 1, 1][:, 1]) + s.evec[1, 1, 1, 1][:, 2] * adjoint(s.evec[1, 1, 1, 1][:, 2])
    threshold_zero!(s.P, p)
    display(s.P)
    display(s.P^2)
    display(s.P ≈ s.P^2)

    println("[3] Check for projected Hamiltonian")
    A = H(p)
    threshold_zero!(A, p)
    display(A)
    B = s.P * H(p) * s.P
    threshold_zero!(B, p)
    display(B)

    println("[4] Check for projected eigen values")
    A = MVector(eigvals(A))
    threshold_zero!(A, p)
    display(A)
    B = MVector(eigvals(B)) # projected
    threshold_zero!(B, p)
    display(B)

    println("[5] Check for orthonormality of eigen vectors")
    A = adjoint(s.evec[1, 1, 1, 1][:, 1]) * s.evec[1, 1, 1, 1][:, 1]
    A = threshold_zero!(A, p)
    display(A)
    B = adjoint(s.evec[1, 1, 1, 1][:, 1]) * s.evec[1, 1, 1, 1][:, 2]
    B = threshold_zero!(B, p)
    display(B)
end

function main()

    p = setParams()

    s = p.sys

    for i in eachindex(p.mList)

        s.chern = 0.0 + 0.0im

        @reset p.m = p.mList[i]

        setBasis!(p)

        # test(p) # Bug check

        updateChern!(p)

        s.chern /= 4(pi^2)
        display(p.m)
        display(s.chern)
        p.ChernList[i] = s.chern
    end

    makeFigure(p.mList, p.ChernList)

end

@time main()