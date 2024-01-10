# The Su-Schriffer-Heeger (SSH) model

Here's a simple example of the SSH Hamiltonian:

```julia
julia> using TopologicalNumbers
julia> function H₀(k, p)
           t₁ = 1
           t₂ = p

           [
               0 t₁+t₂*exp(-im * k)
               t₁+t₂*exp(im * k) 0
           ]
       end
```

The band structure is computed as follows:

```julia
julia> H(k) = H₀(k, 1.1)
julia> showBand(H; value=false, disp=true)
```

![Band structure of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/1a499692-681c-4166-b99c-0fd5cdcd506f)

Next, we can calculate the winding numbers using `calcBerryPhase`:

```julia
julia> calcBerryPhase(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the winding number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the winding numbers for each band (mod 2).
`Total` is a quantity that should always return zero.



One-dimensional phase diagram is given by:

```julia
julia> param = range(-2.0, 2.0, length=1001)
julia> calcPhaseDiagram(H₀, param, "BerryPhase"; plot=true)
```

![One-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/7a6bec77-9140-4257-ba66-8280eef4fe1d)

<!-- Also, two-dimensional phase diagram is given by:

```julia
julia> param = range(-2.0, 2.0, length=101)
julia> calcPhaseDiagram(H₀, param, param, "BerryPhase"; plot=true)
```

![Two-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/adac3cbf-64ce-4324-964f-f42a66948fd4) -->