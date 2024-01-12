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

Next, we can calculate the winding numbers using `BPProblem`:

```julia
julia> prob = BPProblem(H);
julia> sol = solve(prob)
```

The output is:

```julia
BPSolution{Vector{Int64}, Int64}([1, 1], 0)
```

The first argument `TopologicalNumber` in the named tuple is a vector that stores the winding number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the winding numbers for each band (mod 2).
`Total` is a quantity that should always return zero.

You can access these values as follows:

```julia
julia> sol.TopologicalNumber
2-element Vector{Int64}:
 1
 1

julia> sol.Total
0
```



A one-dimensional phase diagram is given by:

```julia
julia> param = range(-2.0, 2.0, length=1001)

julia> prob = BPProblem(H₀);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = -2.0:0.004:2.0, nums = [1 1; 1 1; … ; 1 1; 1 1])
```

![One-dimensional phase diagram of SSH model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/7a6bec77-9140-4257-ba66-8280eef4fe1d)
