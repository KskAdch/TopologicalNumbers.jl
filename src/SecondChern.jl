using LinearAlgebra
using PythonPlot
using LaTeXStrings
using StaticArrays
using Accessors
using Parameters

# Define a mutable struct named Temporal with multiple type parameters
@with_kw mutable struct Temporal{T1,T2,T3,T4,T5,T6}
    chern::T1      # Chern number
    k::T2          # Momentum vector
    evec::T3       # Eigenvectors
    P::T4          # Projection operator
    H::T4          # Hamiltonian
    Ux::T5         # Unitary matrices in the x direction
    Uy::T5         # Unitary matrices in the y direction
    Uz::T5         # Unitary matrices in the z direction
    Uw::T5         # Unitary matrices in the w direction
    Ux_inv::T5     # Inverse unitary matrices in the x direction
    Uy_inv::T5     # Inverse unitary matrices in the y direction
    Uz_inv::T5     # Inverse unitary matrices in the z direction
    Uw_inv::T5     # Inverse unitary matrices in the w direction
    Fxy::T5        # Field strength tensor components
    Fzw::T5
    Fwx::T5
    Fzy::T5
    Fzx::T5
    Fyw::T5
    Ftemp::T6      # Temporary variable for LinearAlgebra.log calculations
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
end

function linkUx!(i, j, l, m, s)
    s.Ux = s.evec[i, j, l, m]' * s.evec[i+1, j, l, m]
end

function linkUx_inv!(i, j, l, m, s)
    s.Ux_inv = s.evec[i+1, j, l, m]' * s.evec[i, j, l, m]
end

function linkUy!(i, j, l, m, s)
    s.Uy = s.evec[i, j, l, m]' * s.evec[i, j+1, l, m]
end

function linkUy_inv!(i, j, l, m, s)
    s.Uy_inv = s.evec[i, j+1, l, m]' * s.evec[i, j, l, m]
end

function linkUz!(i, j, l, m, s)
    s.Uz = s.evec[i, j, l, m]' * s.evec[i, j, l+1, m]
end

function linkUz_inv!(i, j, l, m, s)
    s.Uz_inv = s.evec[i, j, l+1, m]' * s.evec[i, j, l, m]
end

function linkUw!(i, j, l, m, s)
    s.Uw = s.evec[i, j, l, m]' * s.evec[i, j, l, m+1]
end

function linkUw_inv!(i, j, l, m, s)
    s.Uw_inv = s.evec[i, j, l, m+1]' * s.evec[i, j, l, m]

end

function Fxy!(i, j, l, m, s, p)
    linkUx!(i, j, l, m, s)
    linkUy!(i + 1, j, l, m, s)
    linkUx_inv!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m, s)
    s.Fxy = s.Ux * s.Uy * s.Ux_inv * s.Uy_inv
    s.Ftemp .= s.Fxy .* p.r.Nfill ./ tr(s.Fxy)
    s.Fxy = log(s.Ftemp)

end
function Fzw!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s)
    linkUw!(i, j, l + 1, m, s)
    linkUz_inv!(i, j, l, m + 1, s)
    linkUw_inv!(i, j, l, m, s)
    s.Fzw = s.Uz * s.Uw * s.Uz_inv * s.Uw_inv
    s.Ftemp .= s.Fzw * p.r.Nfill / tr(s.Fzw)
    s.Fzw = log(s.Ftemp)
end
function Fwx!(i, j, l, m, s, p)
    linkUw!(i, j, l, m, s)
    linkUx!(i, j, l, m + 1, s)
    linkUw_inv!(i + 1, j, l, m, s)
    linkUx_inv!(i, j, l, m, s)
    s.Fwx = s.Uw * s.Ux * s.Uw_inv * s.Ux_inv
    s.Ftemp .= s.Fwx * p.r.Nfill / tr(s.Fwx)
    s.Fwx = log(s.Ftemp)
end
function Fzy!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s)
    linkUy!(i, j, l + 1, m, s)
    linkUz_inv!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m, s)
    s.Fzy = s.Uz * s.Uy * s.Uz_inv * s.Uy_inv
    s.Ftemp .= s.Fzy * p.r.Nfill / tr(s.Fzy)
    s.Fzy = log(s.Ftemp)
end
function Fzx!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s)
    linkUx!(i, j, l + 1, m, s)
    linkUz_inv!(i + 1, j, l, m, s)
    linkUx_inv!(i, j, l, m, s)
    s.Fzx = s.Uz * s.Ux * s.Uz_inv * s.Ux_inv
    s.Ftemp .= s.Fzx * p.r.Nfill / tr(s.Fzx)
    s.Fzx = log(s.Ftemp)
end
function Fyw!(i, j, l, m, s, p)
    linkUy!(i, j, l, m, s)
    linkUw!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m + 1, s)
    linkUw_inv!(i, j, l, m, s)
    s.Fyw = s.Uy * s.Uw * s.Uy_inv * s.Uw_inv
    s.Ftemp .= s.Fyw * p.r.Nfill / tr(s.Fyw)
    s.Fyw = log(s.Ftemp)
end

function chernF!(i, j, l, m, p)
    s = p.sys

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

    NH = 4
    Nfill = 2

    mN = 1
    # mN = 10
    # mList = range(-5.0, 5.0, length=mN)
    mList = [-3.0]
    ChernList = zeros(ComplexF64, mN)

    r = (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange, NH, Nfill)


    k = @SVector zeros(4)

    evec = [@SMatrix zeros(ComplexF64, NH, Nfill) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

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
    Ftemp = zeros(ComplexF64, Nfill, Nfill) # Temporal variable for LinearAlgebra.log calculation

    chern = zero(ComplexF64)
    sys = Temporal(; chern, k, evec, P, H, Ux, Uy, Uz, Uw, Ux_inv, Uy_inv, Uz_inv, Uw_inv, Fxy, Fzw, Fwx, Fzy, Fzx, Fyw, Ftemp)

    m = zero(Float64)
    (; m, mN, mList, ChernList, basis, r, sys)
end

# Boundary gauge fixing
function boundaryGauge(r, s)
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange), l in eachindex(r.Kzrange)
        s.evec[l, j, i, r.Nz+1] = s.evec[l, j, i, 1]
    end
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange), l in eachindex(r.Kzrange)
        s.evec[l, j, r.Nz+1, i] = s.evec[l, j, 1, i]
    end
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange), l in eachindex(r.Kzrange)
        s.evec[l, r.Ny+1, j, i] = s.evec[l, 1, j, i]
    end
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange), l in eachindex(r.Kzrange)
        s.evec[r.Nz+1, l, j, i] = s.evec[1, l, j, i]
    end

    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange)
        s.evec[j, i, r.Nz+1, r.Nw+1] = s.evec[j, i, 1, 1]
    end
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange)
        s.evec[j, r.Nz+1, i, r.Nw+1] = s.evec[j, 1, i, 1]
    end
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange)
        s.evec[r.Ny+1, j, i, r.Nw+1] = s.evec[1, j, i, 1]
    end
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange)
        s.evec[j, r.Nz+1, r.Nw+1, i] = s.evec[j, 1, 1, i]
    end
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange)
        s.evec[r.Nz+1, j, r.Nw+1, i] = s.evec[1, j, 1, i]
    end
    for i in eachindex(r.Kxrange), j in eachindex(r.Kyrange)
        s.evec[r.Ny+1, r.Nz+1, j, i] = s.evec[1, 1, j, i]
    end

    for i in eachindex(r.Kxrange)
        s.evec[i, r.Ny+1, r.Nz+1, r.Nw+1] = s.evec[i, 1, 1, 1]
    end
    for i in eachindex(r.Kxrange)
        s.evec[r.Nz+1, i, r.Ny+1, r.Nw+1] = s.evec[1, i, 1, 1]
    end
    for i in eachindex(r.Kxrange)
        s.evec[r.Ny+1, r.Nz+1, i, r.Nw+1] = s.evec[1, 1, i, 1]
    end
    for i in eachindex(r.Kxrange)
        s.evec[r.Nz+1, r.Ny+1, r.Nw+1, i] = s.evec[1, 1, 1, i]
    end

    s.evec[r.Nx+1, r.Ny+1, r.Nz+1, r.Nw+1] = s.evec[1, 1, 1, 1]
end

function setBasis!(p)
    r = p.r
    s = p.sys

    for m in eachindex(r.Kwrange), l in eachindex(r.Kzrange), j in eachindex(r.Kyrange), i in eachindex(r.Kxrange)
        s.k = @SVector [r.Kxrange[i], r.Kyrange[j], r.Kzrange[l], r.Kwrange[m]]
        H!(p)
        s.H = eigvecs(s.H)
        s.evec[i, j, l, m] = @view s.H[:, 1:r.Nfill]
    end

    boundaryGauge(p.r, s)
end

function updateChern!(p)
    r = p.r
    for m in eachindex(r.Kwrange)[1:end-1], l in eachindex(r.Kzrange)[1:end-1], j in eachindex(r.Kyrange)[1:end-1], i in eachindex(r.Kxrange)[1:end-1]
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

        updateChern!(p)

        s.chern /= 4pi^2
        println("m: ", p.m, ", Second Chern number: ", s.chern)

        p.ChernList[i] = s.chern
    end

    # makeFigure(p.mList, p.ChernList)

end

@time main()
@time main()

# using Profile
# main()
# Profile.clear_malloc_data()
# main()



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

# 4回目
# m: -3.0, Second Chern number: 0.9793607631927379 - 3.90140613426141e-16im
#  84.573190 seconds (226.39 M allocations: 30.263 GiB, 5.72% gc time)

# m: -3.0, Second Chern number: 0.9793607631927379 - 3.90140613426141e-16im
#  69.130391 seconds (200.98 M allocations: 28.689 GiB, 5.96% gc time)