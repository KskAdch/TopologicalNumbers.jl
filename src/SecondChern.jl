# Functions to calculate and update Link variables and their inverses
# for the links in the x, y, z, and w directions
function linkUx!(i, j, l, m, s::TemporalSecondChern)
    s.Ux = s.evec[i, j, l, m]' * s.evec[i+1, j, l, m]
end
function linkUx_inv!(i, j, l, m, s::TemporalSecondChern)
    s.Ux_inv = s.evec[i+1, j, l, m]' * s.evec[i, j, l, m]
end

function linkUy!(i, j, l, m, s::TemporalSecondChern)
    s.Uy = s.evec[i, j, l, m]' * s.evec[i, j+1, l, m]
end
function linkUy_inv!(i, j, l, m, s::TemporalSecondChern)
    s.Uy_inv = s.evec[i, j+1, l, m]' * s.evec[i, j, l, m]
end

function linkUz!(i, j, l, m, s::TemporalSecondChern)
    s.Uz = s.evec[i, j, l, m]' * s.evec[i, j, l+1, m]
end
function linkUz_inv!(i, j, l, m, s::TemporalSecondChern)
    s.Uz_inv = s.evec[i, j, l+1, m]' * s.evec[i, j, l, m]
end

function linkUw!(i, j, l, m, s::TemporalSecondChern)
    s.Uw = s.evec[i, j, l, m]' * s.evec[i, j, l, m+1]
end
function linkUw_inv!(i, j, l, m, s::TemporalSecondChern)
    s.Uw_inv = s.evec[i, j, l, m+1]' * s.evec[i, j, l, m]

end

# Functions to calculate field strength tensors Fxy, Fzw, Fwx, Fzy, Fzx, Fyw
function Fxy!(i, j, l, m, s::TemporalSecondChern, Nfill)
    linkUx!(i, j, l, m, s)
    linkUy!(i + 1, j, l, m, s)
    linkUx_inv!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fxy = s.Ux * s.Uy * s.Ux_inv * s.Uy_inv
    s.Ftemp .= s.Fxy .* Nfill ./ tr(s.Fxy)
    s.Fxy = log(s.Ftemp)

end
function Fzw!(i, j, l, m, s::TemporalSecondChern, Nfill)
    linkUz!(i, j, l, m, s)
    linkUw!(i, j, l + 1, m, s)
    linkUz_inv!(i, j, l, m + 1, s)
    linkUw_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fzw = s.Uz * s.Uw * s.Uz_inv * s.Uw_inv
    s.Ftemp .= s.Fzw * Nfill / tr(s.Fzw)
    s.Fzw = log(s.Ftemp)
end
function Fwx!(i, j, l, m, s::TemporalSecondChern, Nfill)
    linkUw!(i, j, l, m, s)
    linkUx!(i, j, l, m + 1, s)
    linkUw_inv!(i + 1, j, l, m, s)
    linkUx_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fwx = s.Uw * s.Ux * s.Uw_inv * s.Ux_inv
    s.Ftemp .= s.Fwx * Nfill / tr(s.Fwx)
    s.Fwx = log(s.Ftemp)
end
function Fzy!(i, j, l, m, s::TemporalSecondChern, Nfill)
    linkUz!(i, j, l, m, s)
    linkUy!(i, j, l + 1, m, s)
    linkUz_inv!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fzy = s.Uz * s.Uy * s.Uz_inv * s.Uy_inv
    s.Ftemp .= s.Fzy * Nfill / tr(s.Fzy)
    s.Fzy = log(s.Ftemp)
end
function Fzx!(i, j, l, m, s::TemporalSecondChern, Nfill)
    linkUz!(i, j, l, m, s)
    linkUx!(i, j, l + 1, m, s)
    linkUz_inv!(i + 1, j, l, m, s)
    linkUx_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fzx = s.Uz * s.Ux * s.Uz_inv * s.Ux_inv
    s.Ftemp .= s.Fzx * Nfill / tr(s.Fzx)
    s.Fzx = log(s.Ftemp)
end
function Fyw!(i, j, l, m, s::TemporalSecondChern, Nfill)
    linkUy!(i, j, l, m, s)
    linkUw!(i, j + 1, l, m, s)
    linkUy_inv!(i, j, l, m + 1, s)
    linkUw_inv!(i, j, l, m, s)
    # Calculate Fxy and adjust for normalization
    s.Fyw = s.Uy * s.Uw * s.Uy_inv * s.Uw_inv
    s.Ftemp .= s.Fyw * Nfill / tr(s.Fyw)
    s.Fyw = log(s.Ftemp)
end

# Function to calculate the second Chern number at local k-point using field strength tensors
function chernF!(i, j, l, m, v)
    s = v.sys
    @unpack Nfill = v.r

    Fxy!(i, j, l, m, s, Nfill)
    Fzw!(i, j, l, m, s, Nfill)
    Fwx!(i, j, l, m, s, Nfill)
    Fzy!(i, j, l, m, s, Nfill)
    Fzx!(i, j, l, m, s, Nfill)
    Fyw!(i, j, l, m, s, Nfill)
    # Combine components to calculate the Chern number
    s.Fxy = s.Fxy * s.Fzw + s.Fwx * s.Fzy + s.Fzx * s.Fyw
    s.chern += tr(s.Fxy)
end

# Function to set temporal parameters for the simulation
function setParams(p)

    # Define the ranges and sizes for the grid in the x, y, z, and w directions
    Nx = p.N[1]
    Ny = p.N[2]
    Nz = p.N[3]
    Nw = p.N[4]

    Kxrange = range(-pi, pi, length=Nx + 1)
    Kyrange = range(-pi, pi, length=Ny + 1)
    Kzrange = range(-pi, pi, length=Nz + 1)
    Kwrange = range(-pi, pi, length=Nw + 1)

    NH = p.Hs # Hamiltonian size
    Nfill = p.Nfill # Filling number

    r = (; Nx, Ny, Nz, Nw, Kxrange, Kyrange, Kzrange, Kwrange, NH, Nfill) # Parameters for the ranges

    k = @SVector zeros(4) # Initialize momentum vector

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
    sys = TemporalSecondChern(; chern, k, evec, H, Ux, Uy, Uz, Uw, Ux_inv, Uy_inv, Uz_inv, Uw_inv, Fxy, Fzw, Fwx, Fzy, Fzx, Fyw, Ftemp) # Create Temporal struct

    (; r, sys) # Return parameters
end

# Function to fix the gauge at the boundaries
function boundaryGauge(r, s::TemporalSecondChern)
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
function setBasis!(v, p)
    r = v.r
    s = v.sys

    # Loop over the grid to calculate the eigenvectors
    for m in eachindex(r.Kwrange), l in eachindex(r.Kzrange), j in eachindex(r.Kyrange), i in eachindex(r.Kxrange)
        s.k = @SVector [r.Kxrange[i], r.Kyrange[j], r.Kzrange[l], r.Kwrange[m]]
        # H!(p)
        # s.H = eigvecs(s.H)
        s.H = eigvecs(p.Ham(s.k))
        s.evec[i, j, l, m] = @view s.H[:, 1:r.Nfill]
    end

    boundaryGauge(r, s) # Apply boundary conditions
end

# Function to update the second Chern number
function updateChern!(v, ::UseSingleThread)
    r = v.r
    # Loop over the grid and calculate the Chern number
    for m in eachindex(r.Kwrange)[1:end-1], l in eachindex(r.Kzrange)[1:end-1], j in eachindex(r.Kyrange)[1:end-1], i in eachindex(r.Kxrange)[1:end-1]
        chernF!(i, j, l, m, v)
    end
end

# Function to update the second Chern number
function updateChern!(v, mod::UseMPI)
    r = v.r

    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(r.Kwrange) - 1, nprocs)[myrank+1]

    # Loop over the grid and calculate the Chern number
    for m in idxs
        for l in eachindex(r.Kzrange)[1:end-1], j in eachindex(r.Kyrange)[1:end-1], i in eachindex(r.Kxrange)[1:end-1]
            chernF!(i, j, l, m, v)
        end
    end
    v.sys.chern = mod.MPI.Allreduce(v.sys.chern, mod.MPI.SUM, comm)

    mod.MPI.Barrier(comm)
end

function warn_finiteImaginary(x)
    if abs(imag(x)) > 1e-10
        warn("Imaginary part of the second Chern number has a finite value.")
    end
end

# Main function to execute the simulation
@doc raw"""
"""
function SecondChernPhase(p; parallel::T=UseSingleThread()) where {T<:TopologicalNumbersParallel}
    # @unpack N, Hs = p

    v = setParams(p) # Set parameters

    setBasis!(v, p) # Set the basis

    SecondChernPhase!(v; parallel) # Update the second Chern number

    v.sys.chern # Return the second Chern number

end

# Main function to execute the simulation
@doc raw"""
"""
function SecondChernPhase!(v; parallel::T=UseSingleThread()) where {T<:TopologicalNumbersParallel}
    s = v.sys

    updateChern!(v, parallel) # Update the second Chern number

    s.chern /= 4pi^2 # Normalize the second Chern number

end


# Old method
@doc raw"""

 Calculate the second Chern numbers in the four-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    calcSecondChern(Hamiltonian::Function; Nfill::T1=nothing, N::T2=(30, 30, 30, 30), returnRealValue::Bool=true, parallel::T3=UseSingleThread()) where {T1<:Union{Int,Nothing},T2<:Union{AbstractVector,Tuple},T3<:TopologicalNumbersParallel}

 Arguments
 - `Hamiltionian`: The Hamiltonian matrix with two-dimensional wavenumber `k` as an argument.
 - `Nfill::T1`: The filling number. The default value is `Hs ÷ 2`, where `Hs` is the size of the Hamiltonian matrix.
 - `N::T2`: The numbers of meshes when discretizing the Brillouin Zone. Each element of `N` is the number of meshes in the x, y, z, and w directions, respectively.
 - `returnRealValue::Bool`: An option to return the value of the topological number by an real value. The topological number returns a value of type `Float64` when `true`, and a value of type `ComplexF64` when `false`.


# Definition

# Examples

"""
function calcSecondChern(
    Hamiltonian::Function;
    Nfill::T1=nothing,
    N::T2=(30, 30, 30, 30),
    returnRealValue::Bool=true,
    parallel::T3=UseSingleThread()
) where {T1<:Union{Int,Nothing},T2<:Union{AbstractVector,Tuple},T3<:TopologicalNumbersParallel}

    Hs = size(Hamiltonian(zeros(4)), 1)
    if isnothing(Nfill)
        Nfill = Hs ÷ 2 # Half filling
    end
    p = Params(; Ham=Hamiltonian, Nfill, N, gapless=0.0, rounds=returnRealValue, Hs, dim=4)

    TopologicalNumber = SecondChernPhase(p; parallel)
    warn_finiteImaginary(TopologicalNumber)

    if returnRealValue == true
        TopologicalNumber = real(TopologicalNumber)
    end

    (; TopologicalNumber)
end

@doc raw"""

 Calculate the second Chern numbers in the four-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    solve(prob::SCProblem, alg::T1=FHS2(); parallel::T2=UseSingleThread()) where {T1<:SecondChernAlgorithms,T2<:TopologicalNumbersParallel}

 Arguments
 - `Hamiltionian`: The Hamiltonian matrix with two-dimensional wavenumber `k` as an argument.
 - `Nfill::T1`: The filling number. The default value is `Hs ÷ 2`, where `Hs` is the size of the Hamiltonian matrix.
 - `N::T2`: The numbers of meshes when discretizing the Brillouin Zone. Each element of `N` is the number of meshes in the x, y, z, and w directions, respectively.
 - `returnRealValue::Bool`: An option to return the value of the topological number by an real value. The topological number returns a value of type `Float64` when `true`, and a value of type `ComplexF64` when `false`.


# Definition

# Examples

"""
function solve(
    prob::SCProblem,
    alg::T1=FHS2();
    parallel::T2=UseSingleThread()
) where {T1<:SecondChernAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack H, N, Nfill, RV = prob

    Hs = size(H(zeros(4)), 1)
    if N isa Int
        N = (N, N, N, N)
    end
    if isnothing(Nfill)
        Nfill = Hs ÷ 2 # Half filling
    end
    p = Params(; Ham=H, Nfill, N, rounds=false, returnRealValue=RV, Hs, dim=4)

    TopologicalNumber = SecondChernPhase(p; parallel)
    warn_finiteImaginary(TopologicalNumber)

    if p.returnRealValue == true
        TopologicalNumber = real(TopologicalNumber)
    end

    SCSolution(; TopologicalNumber)
end