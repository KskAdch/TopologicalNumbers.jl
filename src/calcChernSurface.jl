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
Calculate the sliced first Chern numbers in the three-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    solve(prob::WCSProblem, alg::T1=FHSsurface(); parallel::T2=UseSingleThread()) where {T1<:WeylPointsAlgorithms,T2<:TopologicalNumbersParallel}

# Arguments
- `prob::WCSProblem`: The WCSProblem struct that contains the Hamiltonian matrix function in the wave number space and other parameters.
- `alg::T1=FHSsurface()`: The algorithm to use for calculating the sliced first Chern numbers. Default is `FHSsurface` algorithm.
- `parallel::T2=UseSingleThread()`: The parallelization strategy to use. Default is to use a single thread.

# Returns
- `WCSSolution`: A struct that contains the calculated sliced first Chern numbers.

# Examples

```julia
julia> function H₀(k, p) # Weyl
    k1, k2, k3 = k
    t1, t2, t3, m, k0 = p

    h0 = 0
    hx = 2t1*(cos(k1) - cos(k0)) + m*(2 - cos(k2) - cos(k3))
    hy = 2t2*sin(k2)
    hz = 2t3*sin(k3)

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
end
julia> p0 = (1, 1, 1, 2, 2pi*2/5);
julia> H(k) = H₀(k, p0);
julia> prob = WCSProblem(H, "k1");
julia> sol = solve(prob)
WCSSolution{String, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}, Matrix{Int64}}("k1", 6.283185307179587e-5:0.12319971190548208:6.160048427127176, [0 0; 0 0; … ; 0 0; 0 0])
julia> sol.nums
51×2 Matrix{Int64}:
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 1  -1
 1  -1
 1  -1
 1  -1
 1  -1
 1  -1
 1  -1
 1  -1
 1  -1
 1  -1
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
 0   0
```

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

