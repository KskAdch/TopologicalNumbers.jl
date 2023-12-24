using LinearAlgebra
using PythonPlot
using LaTeXStrings
using StaticArrays
using Accessors
using Parameters

@with_kw mutable struct Temporal{T1,T2,T3,T4,T5,T6,T7}
    chern::T1
    k::T2
    evec::T3
    allP::T4
    P::T5
    H::T5
    Ux::T6
    Uy::T6
    Uz::T6
    Uw::T6
    Ux_inv::T6
    Uy_inv::T6
    Uz_inv::T6
    Uw_inv::T6
    Fxy::T6
    Fzw::T6
    Fwx::T6
    Fzy::T6
    Fzx::T6
    Fyw::T6
    Ftemp::T7
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
    # threshold_zero!(p.sys.H, p)
    p.sys.H
end

function H!(p)
    k1, k2, k3, k4 = p.sys.k
    b = p.basis

    h1 = p.m + cos(k1) + cos(k2) + cos(k3) + cos(k4)
    h2 = sin(k1)
    h3 = sin(k2)
    h4 = sin(k3)
    h5 = sin(k4)

    p.sys.H = h1 .* b.g1 .+ h2 .* b.g2 .+ h3 .* b.g3 .+ h4 .* b.g4 .+ h5 .* b.g5
    # threshold_zero!(p.sys.H, p)
end

function linkUx!(i, j, l, m, s, p)
    s.Ux = adjoint(s.evec[i, j, l, m]) * s.evec[i+1, j, l, m]
    # threshold_zero!(s.Ux, p)
end

function linkUx_inv!(i, j, l, m, s, p)
    # s.Ux_inv = adjoint(s.evec[i, j, l, m]) * s.evec[i+1, j, l, m]
    s.Ux_inv = adjoint(s.evec[i+1, j, l, m]) * s.evec[i, j, l, m]
    # threshold_zero!(s.Ux_inv, p)
    # s.Ux_inv = inv(s.Ux_inv)
    # threshold_zero!(s.Ux_inv, p)
end

function linkUy!(i, j, l, m, s, p)
    s.Uy = adjoint(s.evec[i, j, l, m]) * s.evec[i, j+1, l, m]
    # threshold_zero!(s.Uy, p)
end

function linkUy_inv!(i, j, l, m, s, p)
    # s.Uy_inv = adjoint(s.evec[i, j, l, m]) * s.evec[i, j+1, l, m]
    s.Uy_inv = adjoint(s.evec[i, j+1, l, m]) * s.evec[i, j, l, m]
    # threshold_zero!(s.Uy_inv, p)
    # s.Uy_inv = inv(s.Uy_inv)
    # threshold_zero!(s.Uy_inv, p)
end

function linkUz!(i, j, l, m, s, p)
    s.Uz = adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l+1, m]
    # threshold_zero!(s.Uz, p)
end

function linkUz_inv!(i, j, l, m, s, p)
    # s.Uz_inv = adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l+1, m]
    s.Uz_inv = adjoint(s.evec[i, j, l+1, m]) * s.evec[i, j, l, m]
    # threshold_zero!(s.Uz_inv, p)
    # s.Uz_inv = inv(s.Uz_inv)
    # threshold_zero!(s.Uz_inv, p)
end

function linkUw!(i, j, l, m, s, p)
    # s.Uw = s.P * adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l, m+1] * s.P
    s.Uw = adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l, m+1]
    # threshold_zero!(s.Uw, p)
end

function linkUw_inv!(i, j, l, m, s, p)
    # s.Uw_inv = adjoint(s.evec[i, j, l, m]) * s.evec[i, j, l, m+1]
    s.Uw_inv = adjoint(s.evec[i, j, l, m+1]) * s.evec[i, j, l, m]
    # threshold_zero!(s.Uw_inv, p)
    # s.Uw_inv = inv(s.Uw_inv)
    # threshold_zero!(s.Uw_inv, p)
end

function Fxy!(i, j, l, m, s, p)
    linkUx!(i, j, l, m, s, p)
    linkUy!(i + 1, j, l, m, s, p)
    linkUx_inv!(i, j + 1, l, m, s, p)
    linkUy_inv!(i, j, l, m, s, p)
    s.Fxy = s.Ux * s.Uy * s.Ux_inv * s.Uy_inv
    # threshold_zero!(s.Ux, p)
    s.Ftemp .= s.Fxy * p.r.Nfill / tr(s.Fxy)
    # display(s.Fxy)
    s.Fxy = log(s.Ftemp)
    # threshold_zero!(s.Fxy, p)
    # display(s.Ux * s.Uy * s.Ux_inv * s.Uy_inv)
    # display(s.Fxy)
end
function Fzw!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s, p)
    linkUw!(i, j, l + 1, m, s, p)
    linkUz_inv!(i, j, l, m + 1, s, p)
    linkUw_inv!(i, j, l, m, s, p)
    s.Fzw = s.Uz * s.Uw * s.Uz_inv * s.Uw_inv
    # threshold_zero!(s.Uz, p)
    s.Ftemp .= s.Fzw * p.r.Nfill / tr(s.Fzw)
    # display(s.Fzw)
    s.Fzw = log(s.Ftemp)
    # threshold_zero!(s.Fzw, p)
end
function Fwx!(i, j, l, m, s, p)
    linkUw!(i, j, l, m, s, p)
    linkUx!(i, j, l, m + 1, s, p)
    linkUw_inv!(i + 1, j, l, m, s, p)
    linkUx_inv!(i, j, l, m, s, p)
    s.Fwx = s.Uw * s.Ux * s.Uw_inv * s.Ux_inv
    # threshold_zero!(s.Uw, p)
    s.Ftemp .= s.Fwx * p.r.Nfill / tr(s.Fwx)
    # display(s.Fwx)
    s.Fwx = log(s.Ftemp)
    # threshold_zero!(s.Fwx, p)
end
function Fzy!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s, p)
    linkUy!(i, j, l + 1, m, s, p)
    linkUz_inv!(i, j + 1, l, m, s, p)
    linkUy_inv!(i, j, l, m, s, p)
    s.Fzy = s.Uz * s.Uy * s.Uz_inv * s.Uy_inv
    # threshold_zero!(s.Uz, p)
    s.Ftemp .= s.Fzy * p.r.Nfill / tr(s.Fzy)
    # display(s.Fzy)
    s.Fzy = log(s.Ftemp)
    # threshold_zero!(s.Fzy, p)
end
function Fzx!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s, p)
    linkUx!(i, j, l + 1, m, s, p)
    linkUz_inv!(i + 1, j, l, m, s, p)
    linkUx_inv!(i, j, l, m, s, p)
    s.Fzx = s.Uz * s.Ux * s.Uz_inv * s.Ux_inv
    # threshold_zero!(s.Uz, p)
    s.Ftemp .= s.Fzx * p.r.Nfill / tr(s.Fzx)
    # display(s.Fzx)
    s.Fzx = log(s.Ftemp)
    # threshold_zero!(s.Fzx, p)
end
function Fyw!(i, j, l, m, s, p)
    linkUy!(i, j, l, m, s, p)
    linkUw!(i, j + 1, l, m, s, p)
    linkUy_inv!(i, j, l, m + 1, s, p)
    linkUw_inv!(i, j, l, m, s, p)
    s.Fyw = s.Uy * s.Uw * s.Uy_inv * s.Uw_inv
    # threshold_zero!(s.Uy, p)
    s.Ftemp .= s.Fyw * p.r.Nfill / tr(s.Fyw)
    # display(s.Fyw)
    s.Fyw = log(s.Ftemp)
    # threshold_zero!(s.Fyw, p)
end

function chernF!(i, j, l, m, p)
    s = p.sys

    # println("i: ", i, ", j: ", j, ", l: ", l, ", m: ", m)

    # s.P = s.evec[i, j, l, m][:, 1] * adjoint(s.evec[i, j, l, m][:, 1]) + s.evec[i, j, l, m][:, 2] * adjoint(s.evec[i, j, l, m][:, 2])
    # threshold_zero!(s.P, p)
    Fxy!(i, j, l, m, s, p)
    Fzw!(i, j, l, m, s, p)
    Fwx!(i, j, l, m, s, p)
    Fzy!(i, j, l, m, s, p)
    Fzx!(i, j, l, m, s, p)
    Fyw!(i, j, l, m, s, p)
    s.Fxy = s.Fxy * s.Fzw + s.Fwx * s.Fzy + s.Fzx * s.Fyw
    # threshold_zero!(s.Fxy, p)
    s.chern += tr(s.Fxy)
    # s.chern += tr(s.Fxy * s.Fzw + s.Fwx * s.Fzy + s.Fzx * s.Fyw)
    # display(s.Fzw)
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
    # σ₀ = [1 0; 0 1]
    # σ₁ = [0 1; 1 0]
    # σ₂ = [0 -im; im 0]
    # σ₃ = [1 0; 0 -1]
    σ₀ = @SMatrix [1 0; 0 1]
    σ₁ = @SMatrix [0 1; 1 0]
    σ₂ = @SMatrix [0 -im; im 0]
    σ₃ = @SMatrix [1 0; 0 -1]
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

    mN = 1
    # mN = 10
    # mList = range(-5.0, 5.0, length=mN)
    mList = [-3.0]
    ChernList = zeros(ComplexF64, mN)

    r = (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange, threshold, NH, Nfill)


    k = @SVector [0.0, 0.0, 0.0, 0.0]
    # k = [0.0, 0.0, 0.0, 0.0]

    # evec = [@MMatrix zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
    # evec = [zeros(ComplexF64, 4, 2) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
    evec = [@SMatrix zeros(ComplexF64, NH, Nfill) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]
    # evec = [@MMatrix zeros(ComplexF64, 4, 4) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

    allP = [@SMatrix zeros(ComplexF64, NH, NH) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

    P = @SMatrix zeros(ComplexF64, NH, NH)
    H = @SMatrix zeros(ComplexF64, NH, NH)
    Ux = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Uy = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Uz = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Uw = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Ux_inv = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Uy_inv = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Uz_inv = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Uw_inv = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Fxy = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Fzw = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Fwx = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Fzy = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Fzx = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Fyw = @SMatrix zeros(ComplexF64, Nfill, Nfill)
    Ftemp = zeros(ComplexF64, Nfill, Nfill)

    sys = Temporal(; chern=zero(ComplexF64), k, evec, allP, P, H, Ux, Uy, Uz, Uw, Ux_inv, Uy_inv, Uz_inv, Uw_inv, Fxy, Fzw, Fwx, Fzy, Fzx, Fyw, Ftemp)

    (; m=zero(Float64), mN, mList, ChernList, basis, r, sys)
end

function warn_zeroEigen!(p)
    if any(abs.(eigvals(H(p))) .< 1e-10)
        println("Zero eigen value detected!")
    end
end

function setBasis!(p)
    (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange, Nfill) = p.r
    s = p.sys

    for m in eachindex(Kwrange), l in eachindex(Kzrange), j in eachindex(Kyrange), i in eachindex(Kxrange)
        s.k = @SVector [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
        # s.k = [Kxrange[i], Kyrange[j], Kzrange[l], Kwrange[m]]
        # s.evec[i, j, l, m] = @view(eigvecs(H(k, p))[:, 1:2])
        H!(p)
        # warn_zeroEigen!(p)
        s.H = eigvecs(s.H)
        # threshold_zero!(s.H, p)
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
    # y = vecs[i, j, l, m]
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
    for m in eachindex(r.Kwrange)[1:end-1], l in eachindex(r.Kzrange)[1:end-1], j in eachindex(r.Kyrange)[1:end-1], i in eachindex(r.Kxrange)[1:end-1]
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
        if s.evec[i, j, l, m] == zeros(ComplexF64, 4, 4)
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

function plotband(Kxrange, A)
    fig = figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    for i in 1:size(A, 1)
        ax.plot(Kxrange, A[i, :])
    end
    ax.set_xlabel(L"k_x")
    ax.set_ylabel(L"E")
    ax.grid(true)
    tight_layout()
    savefig("band.png")
    plotshow()
end

function band(p)
    s = p.sys
    A = zeros(p.r.NH, p.r.Nx)
    for i in eachindex(p.r.Kxrange)
        s.k = [p.r.Kxrange[i], 0.0, 0.0, 0.0]
        H!(p)
        A[:, i] = eigvals(s.H)
    end
    plotband(p.r.Kxrange, A)
end

function main()

    p = setParams()

    s = p.sys

    for i in eachindex(p.mList)

        s.chern = 0.0 + 0.0im

        @reset p.m = p.mList[i]

        setBasis!(p)

        # setProjectors!(p)

        # band(p)

        # test(p) # Bug check

        updateChern!(p)

        s.chern /= 4pi^2
        println("m: ", p.m, ", Second Chern number: ", s.chern)

        p.ChernList[i] = s.chern
    end

    # makeFigure(p.mList, p.ChernList)

end

@time main()
@time main()



# inv使わないver.
# 840.959875 seconds (2.39 G allocations: 373.897 GiB, 7.07% gc time)

# inv使うver.
# 907.860957 seconds (2.68 G allocations: 579.615 GiB, 8.86% gc time)


# Benchmarks (inv使うver.)

# 1回目
# m: -3.0, Second Chern number: 0.9792956713277512 - 1.0237402388473857e-16im
# 112.240905 seconds (268.87 M allocations: 58.056 GiB, 11.05% gc time)

# m: -3.0, Second Chern number: 0.9792956713277512 - 1.0237402388473857e-16im
# 105.275859 seconds (268.85 M allocations: 58.055 GiB, 9.87% gc time)

# 2回目
# m: -3.0, Second Chern number: 0.9792956713277515 - 1.0237402388473857e-16im
#  97.521954 seconds (268.85 M allocations: 58.055 GiB, 10.00% gc time)

# m: -3.0, Second Chern number: 0.9792956713277515 - 1.0237402388473857e-16im
# 105.103001 seconds (268.85 M allocations: 58.055 GiB, 9.37% gc time)

# 3回目
# m: -3.0, Second Chern number: 0.9793607631927379 - 3.90140613426141e-16im
# 110.723529 seconds (251.33 M allocations: 31.571 GiB, 4.39% gc time)
# m: -3.0, Second Chern number: 0.9793607631927379 - 3.90140613426141e-16im
#  96.559988 seconds (225.91 M allocations: 29.996 GiB, 4.40% gc time)