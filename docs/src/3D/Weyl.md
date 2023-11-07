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
julia> p0 = (1, 1, 1, 2, 2pi*2/5)
julia> H(k) = H₀(k, p0)
julia> calcChernSurface(H, "k1"; plot = true)
```

```julia
julia> (kn = "k1", param = 6.283185307179587e-5:0.12319971190548208:6.160048427127176, nums = [0 0; 0 0; … ; 0 0; 0 0])
```

The first argument `kn` in the named tuple is a fixed component of the wavenumber vector. 
The second argument `param` stores a range of fixed wavenumber vectors.
Take $2\pi$ from $0$ by default.
The third argument `nums` is a matrix that stores the Chern number of each band in each `kn`.

![Chern surface of k1](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/47e4a36a-4c67-4158-ada1-31a957d19b87)

```julia
julia> calcChernSurface(H, "k2"; plot = true)
```

![Chern surface of k2](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/f44c1499-93dc-488e-b398-884d41ca1a6a)

```julia
julia> calcChernSurface(H, "k3"; plot = true)
```

![Chern surface of k3](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/d4a98efe-8fa1-4960-9ccc-c6b71be62661)

Or, Weyl points can be found as follows:

```julia
julia> result = findWeylPoint(H)
julia> 2pi*result.WeylPoint[1] / result.N
```

```julia
julia> [[2.5132741228718345, 0.0, 0.0], [3.7699111843077517, 0.0, 0.0]]
```

The values returned by findWeylPoint are as follows:

```julia
julia> result
```

```julia
julia> (WeylPoint = [[[4000, 0, 0], [6000, 0, 0]], [[4000, 0, 0], [6000, 0, 0]]], N = 10000, Nodes = [[1, -1], [-1, 1]])
```

The second argument `N` in the named tuple is the number of Brillouin zone divisions.
The first argument `WeylPoint` represents the the wavenumber vector of position with the Weyl point.
The third argument `Nodes` stores the node at the Wel point.