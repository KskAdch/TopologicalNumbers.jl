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

![Chern surface of k1]()

```julia
julia> calcChernSurface(H, "k2"; plot = true)
```

![Chern surface of k2]()

```julia
julia> calcChernSurface(H, "k3"; plot = true)
```

![Chern surface of k3]()

Or, Weyl points can be found as follows:

```julia
julia> result = findWeylPoint(H)
julia> 2pi*result.WeylPoint / result.N
```