# Haldane model

Hamiltonian of Haldane model is given by:

```julia
julia> function H₀(k, p) # landau
    k1, k2 = k
    J = 1.0
    K = 1.0
    ϕ, M = p

    h0 = 2K * cos(ϕ) * (cos(k1) + cos(k2) + cos(k1 + k2))
    hx = -J * (1 + cos(k1) + cos(k2))
    hy = -J * (-sin(k1) + sin(k2))
    hz = M + 2K * sin(ϕ) * (sin(k1) + sin(k2) - sin(k1 + k2))

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
end
```

The band structure is computed as follows:

```julia
julia> H(k) = H₀(k, (π / 2, 1.0))
julia> showBand(H; value=false, disp=true)
```


![Band structure of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/c9c0aa00-c412-4a4c-be94-d080a62bd14a)


Then we can compute the Chern numbers using `calcChern`:

```julia
julia> calcChern(H)
```

The output is:

```julia
(TopologicalNumber = [-1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.



One-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (p, 2.5))

julia> param = range(-π, π, length=1000)
julia> calcPhaseDiagram(H, param, "Chern"; plot=true)
```

![One-dimensional phase diagram of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/bbbbd989-b2e2-4073-812f-a7c81bb253ba)

Also, two-dimensional phase diagram is given by:

```julia
julia> param1 = range(-π, π, length=100)
julia> param2 = range(-6.0, 6.0, length=100)
julia> calcPhaseDiagram(H₀, param1, param2, "Chern"; plot=true)
```

![Two-dimensional phase diagram of Haldane model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/81e4648d-4b6a-49e9-bb33-081557495a20)