# The Bernevig-Hughes-Zhang (BHZ) model

As an example of a two-dimensional topological insulator, the BHZ model is presented here:

```julia
julia> using LinearAlgebra
julia> function H₀(k, p) # BHZ
    k1, k2 = k
    tₛₚ = 1
    t₁ = ϵ₁ = 2
    ϵ₂, t₂ = p

    ϵ = -t₁*(cos(k1) + cos(k2)) + ϵ₁/2
    R1 = 0
    R2 = 0
    R3 = 2tₛₚ*sin(k2)
    R4 = 2tₛₚ*sin(k1)
    R0 = -t₂*(cos(k1) + cos(k2)) + ϵ₂/2

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    I0 = Matrix{Int64}(I, 4, 4)
    a1 = kron(sz, sx)
    a2 = kron(sz, sy)
    a3 = kron(sz, sz)
    a4 = kron(sy, s0)
    a0 = kron(sx, s0)

    ϵ .* I0 .+ R1 .* a1 .+ R2 .* a2 .+ R3 .* a3 .+ R4 .* a4 .+ R0 .* a0
end
```

To calculate the dispersion, execute:

```julia
julia> H(k) = H₀(k, (2, 2))
julia> showBand(H; value=false, disp=true)
```

![Dispersion of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/de14907c-777f-4667-810b-54c10888dfa1)


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
julia> H(k, p) = H₀(k, (p, 0.25))

julia> param = range(-2, 2, length=1000)
julia> calcPhaseDiagram(H, param, "Z2"; plot=true)
```

![One-dimensional phase diagram of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/5d6d5364-68d0-4423-8ecf-49bf0538af63)


Also, two-dimensional phase diagram is given by:

```julia
julia> param1 = range(-2, 2, length=100)
julia> param2 = range(-0.5, 0.5, length=100)
julia> calcPhaseDiagram(H₀, param1, param2, "Z2"; plot=true)
```


![Two-dimensional phase diagram of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/802eedbe-c893-44b4-8267-d80e1745415a)