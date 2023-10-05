# The Kane-Mele model

As an example of a two-dimensional topological insulator, the Kane-Mele model is presented here:

```julia
julia> function H₀(k, p) # Kane-Mele
    k1, k2 = k
    t, λₛₒ = p

    R1 = 0
    R2 = 0
    R3 = 2λₛₒ*(sin(k1) - sin(k2) - sin(k1-k2))
    R4 = -t*(sin(k1) + sin(k2))
    R0 = -t*(cos(k1) + cos(k2) + 1)

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    a1 = kron(sz, sx)
    a2 = kron(sz, sy)
    a3 = kron(sz, sz)
    a4 = kron(sy, s0)
    a0 = kron(sx, s0)

    R1 .* a1 .+ R2 .* a2 .+ R3 .* a3 .+ R4 .* a4 .+ R0 .* a0
end
```

To calculate the dispersion, execute:

```julia
julia> H(k) = H₀(k, (1, 0.5))
julia> showBand(H; value=false, disp=true)
```

![Dispersion of Kane-Mele model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/2e82e405-aaa3-4851-bd71-5ba12b1c6b3c)


Next, we can compute the $\mathbb{Z}_2$ numbers using `calcZ2`:

```julia
julia> calcZ2(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the $\mathbb{Z}_2$ number for each pair of two energy bands. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the $\mathbb{Z}_2$ numbers for each pair of two energy bands.
`Total` is a quantity that should always return zero.


One-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (p, 0.5))

julia> param = range(-1.0, 1.0, length=1001)
julia> calcPhaseDiagram(H, param, "Z2"; plot=true)
```

![One-dimensional phase diagram of Kane-Mele model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/a6b58620-9deb-4594-9c32-4a728a79c7d0)


Also, two-dimensional phase diagram is given by:

```julia
julia> param = range(-1.0, 1.0, length=101)
julia> calcPhaseDiagram(H₀, param, param, "Z2"; plot=true)
```


![Two-dimensional phase diagram of Kane-Mele model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/a5a3e842-ec96-4da5-8e0c-f03a2b551ed3)