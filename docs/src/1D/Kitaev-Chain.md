# The Kitaev chain model

Here's a simple example of the Kitaev chain Hamiltonian:

```julia
julia> using TopologicalNumbers
julia> function H₀(k, p)
            μ, Δ = p
            t = 1

            [
                -μ-2t*cos(k) 2im*Δ*sin(k)
                -2im*Δ*sin(k) μ+2t*cos(k)
            ]
        end
```

The band structure is computed as follows:

```julia
julia> H(k) = H₀(k, (1.0, 0.5))
julia> showBand(H; value=false, disp=true)
```

![Band structure of Kitaev chain model]()

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

![One-dimensional phase diagram of Kitaev chain model]()

Also, two-dimensional phase diagram is given by:

```julia
julia> param = range(-2.0, 2.0, length=101)
julia> calcPhaseDiagram(H₀, param, param, "BerryPhase"; plot=true)
```

![Two-dimensional phase diagram of Kitaev chain model]()