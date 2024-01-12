# Weyl semimetal

A three-dimensional example is presented here:

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
```

The location of the Weyl point can be determined by fixing one component of the wavenumber vector and calculating the Chern number as follows:

```julia
julia> p0 = (1, 1, 1, 2, 2pi*2/5);
julia> H(k) = H₀(k, p0);
julia> prob = WCSProblem(H, "k1");
julia> sol = solve(prob; plot = true)
WCSSolution{String, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}, Matrix{Int64}}("k1", 6.283185307179587e-5:0.12319971190548208:6.160048427127176, [0 0; 0 0; … ; 0 0; 0 0])
```

The first argument `kn` in the named tuple is a fixed component of the wavenumber vector. 
The second argument `param` stores a range of fixed wavenumber vectors.
Take $2\pi$ from $0$ by default.
The third argument `nums` is a matrix that stores the Chern number of each band in each `kn`.

You can access these values as follows:

```julia
julia> sol.kn
"k1"

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

julia> sol.param
6.283185307179587e-5:0.12319971190548208:6.160048427127176
```

![Chern surface of k1](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/02ba7790-4997-495a-9b74-48964271410d)

```julia
julia> prob = WCSProblem(H, "k2");
julia> sol = solve(prob; plot = true)
```

![Chern surface of k2](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/35aeb1af-41b0-4b91-9603-81a2d58772c1)

```julia
julia> prob = WCSProblem(H, "k3");
julia> sol = solve(prob; plot = true)
```

![Chern surface of k3](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/2827affe-0c57-42f9-828c-2740818ea425)

Or, Weyl points can be found as follows:

```julia
julia> prob = WPProblem(H);
julia> result = solve(prob)
WPSolution{Vector{Vector{Vector{Int64}}}, Int64, Vector{Vector{Int64}}}([[[4000, 0, 0], [6000, 0, 0]], [[4000, 0, 0], [6000, 0, 0]]], 10000, [[1, -1], [-1, 1]])
```

The second argument `N` in the named tuple is the number of Brillouin zone divisions.
The first argument `WeylPoint` represents the the wavenumber vector of position with the Weyl points.
The third argument `Nodes` stores the node at the Weyl points.


`WeylPoint` can be converted to a wavenumber vector as follows:

```julia
julia> 2pi*result.WeylPoint[1] / result.N .- pi*[ones(3), ones(3)]
2-element Vector{Vector{Float64}}:
 [-0.6283185307179586, -3.141592653589793, -3.141592653589793]
 [0.6283185307179586, -3.141592653589793, -3.141592653589793]
```


If you already know the wavenumber vector of the Weyl point, you can calculate the node as follows:

```julia
julia> H(k) = H₀(k .- 2pi*1e-5, p0)
julia> prob = WNProblem(H, result.WeylPoint[1][1], result.N);
julia> sol = solve(prob)
WNSolution{Vector{Int64}, Vector{Int64}, Int64}([1, -1], [4000, 0, 0], 10000)
```
