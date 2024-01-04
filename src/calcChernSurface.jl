function nums_k1(H, kn_range, N, gapless, rounds, Hs; parallel::T=UseSingleThread()) where {T<:TopologicalNumbersParallel}
    H0(k, p) = H([p, k[1], k[2]])
    Hamiltonian(k) = H0(k, 0.0)

    param = Params(; Ham=Hamiltonian, dim=2, N, gapless, rounds, Hs)

    calc_data1D(H0, kn_range, "Chern", parallel, param)
end

function nums_k2(H, kn_range, N, gapless, rounds, Hs; parallel::T=UseSingleThread()) where {T<:TopologicalNumbersParallel}
    H0(k, p) = H([k[1], p, k[2]])
    Hamiltonian(k) = H0(k, 0.0)

    param = Params(; Ham=Hamiltonian, dim=2, N, gapless, rounds, Hs)

    calc_data1D(H0, kn_range, "Chern", parallel, param)
end

function nums_k3(H, kn_range, N, gapless, rounds, Hs; parallel::T=UseSingleThread()) where {T<:TopologicalNumbersParallel}
    H0(k, p) = H([k[1], k[2], p])
    Hamiltonian(k) = H0(k, 0.0)

    param = Params(; Ham=Hamiltonian, dim=2, N, gapless, rounds, Hs)

    calc_data1D(H0, kn_range, "Chern", parallel, param)
end

@doc raw"""
    calcChernSurface(H::Function, kn::String; kn_mesh::Int=51, N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false)

 Arguments
 - Hamiltionian::Function: The Hamiltonian matrix with three-dimensional wavenumber `k` as an argument.
 - kn::String: Compute the Chern number of the plane perpendicular to the `"kn"` direction in Brillouin zone (`"k1"`, `"k2"`, `"k3"`).
 - kn_mesh::T: Number of mesh in `"kn"` direction.
 - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.
"""
function calcChernSurface(
    H::Function,
    kn::String;
    kn_mesh::Int=51,
    N::Int=51,
    gapless::Real=0.0,
    rounds::Bool=true,
    parallel::T=UseSingleThread(),
    plot::Bool=false
) where {T<:TopologicalNumbersParallel}

    Hs = size(H(zeros(3)), 1)
    kn_range = 2pi * (0:(kn_mesh-1)) / kn_mesh .+ 2pi * 1e-5

    if kn == "k1"
        nums = nums_k1(H, kn_range, N, gapless, rounds, Hs; parallel)
    elseif kn == "k2"
        nums = nums_k2(H, kn_range, N, gapless, rounds, Hs; parallel)
    elseif kn == "k3"
        nums = nums_k3(H, kn_range, N, gapless, rounds, Hs; parallel)
    else
        throw(ArgumentError("Unknown keyword $kn"))
    end

    if plot == true
        plot1D(nums, kn_range)
    end

    (; kn, param=kn_range, nums)
end


@doc raw"""
    solve(prob::WCSProblem, alg::T1=FHSsurface(); parallel::T2=UseSingleThread(), plot::Bool=false) where {T1<:WeylPointsAlgorithms,T2<:TopologicalNumbersParallel}

 Arguments
 - Hamiltionian::Function: The Hamiltonian matrix with three-dimensional wavenumber `k` as an argument.
 - kn::String: Compute the Chern number of the plane perpendicular to the `"kn"` direction in Brillouin zone (`"k1"`, `"k2"`, `"k3"`).
 - kn_mesh::T: Number of mesh in `"kn"` direction.
 - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.
"""
function solve(
    prob::WCSProblem,
    alg::T1=FHSsurface();
    parallel::T2=UseSingleThread(),
    plot::Bool=false
) where {T1<:WeylPointsAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack H, kn, kn_mesh, N, gapless, rounds = prob

    Hs = size(H(zeros(3)), 1)
    kn_range = 2pi * (0:(kn_mesh-1)) / kn_mesh .+ 2pi * 1e-5

    if kn == "k1"
        nums = nums_k1(H, kn_range, N, gapless, rounds, Hs; parallel)
    elseif kn == "k2"
        nums = nums_k2(H, kn_range, N, gapless, rounds, Hs; parallel)
    elseif kn == "k3"
        nums = nums_k3(H, kn_range, N, gapless, rounds, Hs; parallel)
    else
        throw(ArgumentError("Unknown keyword $kn"))
    end

    if plot == true
        plot1D(nums, kn_range)
    end

    WCSSolution(; kn, param=kn_range, nums)
end

