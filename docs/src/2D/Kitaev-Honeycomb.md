# Kitaev honeycomb model

Hamiltonian of Kitaev honeycomb model is given by:

```julia
julia> function H₀(k, p) # Kitaev
            k1, k2 = k
            K₁ = 1
            K₂ = p

            hx = -K₂ * (sin(k2) - sin(k1) + sin(k1 - k2))
            hy = -K₁ * (sin(k1) + sin(k2))
            hz = K₁ * (cos(k1) + cos(k2) + 1)

            sx = [0 1; 1 0]
            sy = [0 -im; im 0]
            sz = [1 0; 0 -1]

            hx .* sx .+ hy .* sy .+ hz .* sz
       end
```

The band structure is computed as follows:

```julia
julia> H(k) = H₀(k, 0.5)
julia> showBand(H; value=false, disp=true)
```


![Band structure of Kitaev honeycomb model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/68e8a814-7f5d-4a9c-9abb-2dad90e808e9)


Then we can compute the Chern numbers using `calcChern`:

```julia
julia> calcChern(H)
```

The output is:

```julia
(TopologicalNumber = [1, -1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.



One-dimensional phase diagram is given by:

```julia
julia> param = range(-1, 1, length=1000)
julia> calcPhaseDiagram(H₀, param, "Chern"; plot=true)
```

![One-dimensional phase diagram of Kitaev honeycomb model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/1af15ce6-6274-4816-b0c6-50a8762c18a6)

<!-- Also, two-dimensional phase diagram is given by:

```julia
julia> param1 = range(-1.0, 1.0, length=100)
julia> param2 = range(-1, 1, length=2)
julia> calcPhaseDiagram(H₀, param1, param2, "Chern"; plot=true)
```

![Two-dimensional phase diagram of Kitaev honeycomb model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/7801d67e-faf6-435b-aa2a-20fb721274b1) -->
