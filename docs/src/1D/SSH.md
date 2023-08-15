# The Su-Schriffer-Heeger (SSH) model

Here's a simple example of the SSH Hamiltonian:

```julia
julia> using TopologicalNumbers
julia> function H₀(k, p)
            [
                0 p[1]+p[2]*exp(-im * k)
                p[1]+p[2]*exp(im * k) 0
            ]
        end
```

The band structure is computed as follows:

```julia
julia> H(k) = H₀(k, (0.9, 1.0))
julia> showBand(H; value=false, disp=true)
```

![Band structure of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/a586aa22-6c79-454e-a82f-6f5056d98f6c)

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
julia> H(k, p) = H₀(k, (p, 1.0))

julia> param = range(-2.0, 2.0, length=1001)
julia> calcPhaseDiagram(H, param, "BerryPhase"; plot=true)
```

![One-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/2b53e455-83ee-42d5-9824-84120c2be093)

Also, two-dimensional phase diagram is given by:

```julia
julia> param = range(-2.0, 2.0, length=101)
julia> calcPhaseDiagram(H₀, param, param, "BerryPhase"; plot=true)
```

![Two-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/0ceef1a3-01fd-4e8b-9f01-4a4932039d26)