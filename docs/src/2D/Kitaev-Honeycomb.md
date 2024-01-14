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


Then we can compute the Chern numbers using `FCProblem`:

```julia
julia> prob = FCProblem(H);
julia> sol = solve(prob)
```

The output is:

```julia
FCSolution{Vector{Int64}, Int64}([-1, 1], 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.

You can access these values as follows:

```julia
julia> sol.TopologicalNumber
2-element Vector{Int64}:
 -1
  1

julia> sol.Total
0
```


One-dimensional phase diagram is given by:

```julia
julia> param = range(-1, 1, length=1000);
julia> calcPhaseDiagram(H₀, param, "Chern"; plot=true)
(param = -1.0:0.002002002002002002:1.0, nums = [1 -1; 1 -1; … ; -1 1; -1 1])
```

![One-dimensional phase diagram of Kitaev honeycomb model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/1af15ce6-6274-4816-b0c6-50a8762c18a6)

