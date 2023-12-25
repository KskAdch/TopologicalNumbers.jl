using LinearAlgebra
using PythonPlot
using LaTeXStrings
using StaticArrays
using Accessors
using Parameters

# Define a mutable struct named Temporal with multiple type parameters
@with_kw mutable struct Temporal{T1,T2,T3,T4,T5,T6}
    chern::T1      # Second Chern number
    k::T2          # Momentum vector
    evec::T3       # Eigenvectors
    H::T4          # Hamiltonian
    Ux::T5         # Link variables
    Uy::T5
    Uz::T5
    Uw::T5
    Ux_inv::T5     # Inverse Link variables
    Uy_inv::T5
    Uz_inv::T5
    Uw_inv::T5
    Fxy::T5        # Field strength tensor components
    Fzw::T5
    Fwx::T5
    Fzy::T5
    Fzx::T5
    Fyw::T5
    Ftemp::T6      # Temporary variable for LinearAlgebra.log calculations
end

# Define a tensor product function
⊗(x, y) = kron(x, y)

# Functions to calculate and update Link variables and their inverses
# for the links in the x, y, z, and w directions
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

# Functions to calculate field strength tensors Fxy, Fzw, Fwx, Fzy, Fzx, Fyw
function Fxy!(i, j, l, m, s, p)
    linkUx!(i, j, l, m, s)
    linkUy!(i + 1, j, l, m, s)
    linkUx_inv!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fxy = s.Ux * s.Uy * s.Ux_inv * s.Uy_inv
    s.Ftemp .= s.Fxy .* p.r.Nfill ./ tr(s.Fxy)
    s.Fxy = log(s.Ftemp)

end
function Fzw!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s)
    linkUw!(i, j, l + 1, m, s)
    linkUz_inv!(i, j, l, m + 1, s)
    linkUw_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fzw = s.Uz * s.Uw * s.Uz_inv * s.Uw_inv
    s.Ftemp .= s.Fzw * p.r.Nfill / tr(s.Fzw)
    s.Fzw = log(s.Ftemp)
end
function Fwx!(i, j, l, m, s, p)
    linkUw!(i, j, l, m, s)
    linkUx!(i, j, l, m + 1, s)
    linkUw_inv!(i + 1, j, l, m, s)
    linkUx_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fwx = s.Uw * s.Ux * s.Uw_inv * s.Ux_inv
    s.Ftemp .= s.Fwx * p.r.Nfill / tr(s.Fwx)
    s.Fwx = log(s.Ftemp)
end
function Fzy!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s)
    linkUy!(i, j, l + 1, m, s)
    linkUz_inv!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fzy = s.Uz * s.Uy * s.Uz_inv * s.Uy_inv
    s.Ftemp .= s.Fzy * p.r.Nfill / tr(s.Fzy)
    s.Fzy = log(s.Ftemp)
end
function Fzx!(i, j, l, m, s, p)
    linkUz!(i, j, l, m, s)
    linkUx!(i, j, l + 1, m, s)
    linkUz_inv!(i + 1, j, l, m, s)
    linkUx_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fzx = s.Uz * s.Ux * s.Uz_inv * s.Ux_inv
    s.Ftemp .= s.Fzx * p.r.Nfill / tr(s.Fzx)
    s.Fzx = log(s.Ftemp)
end
function Fyw!(i, j, l, m, s, p)
    linkUy!(i, j, l, m, s)
    linkUw!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m + 1, s)
    linkUw_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fyw = s.Uy * s.Uw * s.Uy_inv * s.Uw_inv
    s.Ftemp .= s.Fyw * p.r.Nfill / tr(s.Fyw)
    s.Fyw = log(s.Ftemp)
end

# Function to calculate the second Chern number at local k-point using field strength tensors
function chernF!(i, j, l, m, p)
    s = p.sys

    Fxy!(i, j, l, m, s, p)
    Fzw!(i, j, l, m, s, p)
    Fwx!(i, j, l, m, s, p)
    Fzy!(i, j, l, m, s, p)
    Fzx!(i, j, l, m, s, p)
    Fyw!(i, j, l, m, s, p)
    # Combine components to calculate the Chern number
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

# Function to set temporal parameters for the simulation
function setParams()
    # Define Pauli matrices and Gamma matrices
    # σ₀ = @SMatrix [1 0; 0 1]
    # σ₁ = @SMatrix [0 1; 1 0]
    # σ₂ = @SMatrix [0 -im; im 0]
    # σ₃ = @SMatrix [1 0; 0 -1]
    σ₀ = [1 0; 0 1]
    σ₁ = [0 1; 1 0]
    σ₂ = [0 -im; im 0]
    σ₃ = [1 0; 0 -1]
    g1 = σ₁ ⊗ σ₀
    g2 = σ₂ ⊗ σ₀
    g3 = σ₃ ⊗ σ₁
    g4 = σ₃ ⊗ σ₂
    g5 = σ₃ ⊗ σ₃

    basis = (; g1, g2, g3, g4, g5) # Define a basis set

    # Define the ranges and sizes for the grid in the x, y, z, and w directions
    Nx = 30
    Ny = 30
    Nz = 30
    Nw = 30

    Kxrange = range(-pi, pi, length=Nx + 1)
    Kyrange = range(-pi, pi, length=Ny + 1)
    Kzrange = range(-pi, pi, length=Nz + 1)
    Kwrange = range(-pi, pi, length=Nw + 1)

    NH = 4 # Hamiltonian size
    Nfill = 2 # Filling number

    # mN = 1 # Number of m values
    mN = 50 # Number of m values
    mList = range(-4.9, 4.9, length=mN) # List of m values
    # mList = [-3.0] # List of m values
    ChernList = zeros(ComplexF64, mN) # Initialize Chern number list

    r = (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange, NH, Nfill) # Parameters for the ranges


    # k = @SVector zeros(4) # Initialize momentum vector
    k = zeros(4) # Initialize momentum vector

    # Initialize matrices for eigenvectors, Hamiltonian, unitary matrices, and field strength tensors
    evec = [@SMatrix zeros(ComplexF64, NH, Nfill) for i in 1:Nx+1, j in 1:Ny+1, l in 1:Nz+1, m in 1:Nw+1]

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

    chern = zero(ComplexF64) # Initialize Chern number
    sys = Temporal(; chern, k, evec, H, Ux, Uy, Uz, Uw, Ux_inv, Uy_inv, Uz_inv, Uw_inv, Fxy, Fzw, Fwx, Fzy, Fzx, Fyw, Ftemp) # Create Temporal struct

    m = zero(Float64) # Initialize m value
    (; m, mN, mList, ChernList, basis, r, sys) # Return parameters
end

# Function to fix the gauge at the boundaries
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

# Function to set the basis for calculations
function setBasis!(p)
    r = p.r
    s = p.sys

    # Loop over the grid to calculate the eigenvectors
    for m in eachindex(r.Kwrange), l in eachindex(r.Kzrange), j in eachindex(r.Kyrange), i in eachindex(r.Kxrange)
        # s.k = @SVector [r.Kxrange[i], r.Kyrange[j], r.Kzrange[l], r.Kwrange[m]]
        s.k .= [r.Kxrange[i], r.Kyrange[j], r.Kzrange[l], r.Kwrange[m]]
        # H!(p)
        # s.H = eigvecs(s.H)
        s.H = eigvecs(H(p))
        s.evec[i, j, l, m] = @view s.H[:, 1:r.Nfill]
    end

    boundaryGauge(p.r, s) # Apply boundary conditions
end

# Function to update the second Chern number
function updateChern!(p)
    r = p.r
    # Loop over the grid and calculate the Chern number
    for m in eachindex(r.Kwrange)[1:end-1], l in eachindex(r.Kzrange)[1:end-1], j in eachindex(r.Kyrange)[1:end-1], i in eachindex(r.Kxrange)[1:end-1]
        chernF!(i, j, l, m, p)
    end
end

# Main function to execute the simulation
function main()

    p = setParams() # Set parameters

    s = p.sys

    # Loop over m values to calculate the Chern number for each
    for i in eachindex(p.mList)

        s.chern = 0.0 + 0.0im # Reset second Chern number

        @reset p.m = p.mList[i] # Set current m value

        setBasis!(p) # Set the basis

        updateChern!(p) # Update the second Chern number

        s.chern /= 4pi^2 # Normalize the second Chern number
        println("m: ", p.m, ", Second Chern number: ", s.chern)

        p.ChernList[i] = s.chern # Store the calculated Chern number
    end

    makeFigure(p.mList, p.ChernList)

end

function calcSecondChern()
    main()
end

@time main()
# @time main()


# @views function ChernPhase!(TopologicalNumber, p::Params) # chern number # Bug
#     @unpack N, Hs = p
#     TopologicalNumber[:] .= zero(Float64)
#     Link0 = zeros(ComplexF64, Hs, 2, N)
#     Link1 = zeros(ComplexF64, Hs, 2, N)
#     LinkN = zeros(ComplexF64, Hs, 2, N)
#     link1 = zeros(ComplexF64, Hs)
#     link2 = zeros(ComplexF64, Hs)

#     psi_0 = zeros(ComplexF64, N, Hs, Hs)
#     psi_1 = zeros(ComplexF64, N, Hs, Hs)
#     psi_N = zeros(ComplexF64, N, Hs, Hs)
#     Evec0 = zeros(N, Hs)
#     Evec1 = zeros(N, Hs)

#     psi00 = zeros(ComplexF64, Hs, Hs)
#     psi10 = zeros(ComplexF64, Hs, Hs)
#     psi01 = zeros(ComplexF64, Hs, Hs)
#     Enevec = zeros(Hs)

#     phi = zeros(Hs)
#     dphi = zeros(Hs)

#     for j in 1:N
#         U!(Link0, Link1, LinkN, link1, link2, psi_0, psi_1, psi_N, Evec0, Evec1, psi00, psi10, psi01, Enevec, j, p)
#         for i in 1:N
#             F!(phi, dphi, i, j, Link0, Link1, LinkN, p)
#             TopologicalNumber[:] .+= phi[:]
#         end
#     end
# end

# @doc raw"""

#  Calculate the first Chern numbers in the two-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

#     calcChern(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

#  Arguments
#  - Hamiltionian::Function: The Hamiltonian matrix with two-dimensional wavenumber `k` as an argument.
#  - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
#  - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
#  - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.


# # Definition
#  The first Chern number of the $n$th band $\nu_{n}$ is defined by
# ```math
# \nu_{n}=\frac{1}{2\pi}\sum_{\bm{k}\in\mathrm{BZ}}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,2}^{*}(\bm{k})\right]\right]
# ```
#  The range $\mathrm{BZ}$(Brillouin Zone) is $\bm{k}\in[0,2\pi]^{2}$. $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
# ```math
# U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
# ```
#  $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $n$th band.
# """
# function calcChern(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

#     Hs = size(Hamiltonian(zeros(2)), 1)
#     p = Params(; Hamiltonian, N, gapless, rounds, Hs, dim=2)

#     TopologicalNumber = zeros(Hs)

#     ChernPhase!(TopologicalNumber, p)

#     if rounds == true
#         TopologicalNumber = round.(Int, TopologicalNumber)
#     end

#     Total = sum(TopologicalNumber)

#     (; TopologicalNumber, Total)
# end