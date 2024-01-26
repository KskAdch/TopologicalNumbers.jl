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
You can also use our preset Hamiltonian function `KitaevChain` to define the same Hamiltonian matrix as follows:

```julia
julia> H₀(k, p) = KitaevChain(k, p)
```

The band structure is computed as follows:

```julia
julia> H(k) = H₀(k, (-1.0, 0.5))
julia> showBand(H; value=false, disp=true)
```

![Band structure of Kitaev chain model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/7a11c4b0-1949-457d-8cee-9b57f43af2f1)

Next, we can calculate the winding numbers using `BPProblem`:

```julia
julia> prob = BPProblem(H);
julia> sol = solve(prob)
```

The output is:

```julia
BPSolution{Vector{Int64}, Int64}([1, 1], 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the winding number for each band. 
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



One-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (p, 1.0));
julia> param = range(-3.0, 0, length=601);

julia> prob = BPProblem(H);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = -3.0:0.005:0.0, nums = [0 0; 0 0; … ; 1 1; 1 1])
```

![One-dimensional phase diagram of Kitaev chain model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/6b0798c3-ffa0-4bfb-b724-92aee1eda384)

Also, two-dimensional phase diagram is given by:

```julia
julia> param1 = range(-3.0, 3.0, length=101);
julia> param2 = range(-1.0, 1.0, length=101);

julia> prob = BPProblem(H₀);
julia> calcPhaseDiagram(prob, param1, param2; plot=true)
(param1 = -3.0:0.06:3.0, param2 = -1.0:0.02:1.0, nums = [0 0 … 0 0; 0 0 … 0 0;;; 0 0 … 0 0; 0 0 … 0 0;;; 0 0 … 0 0; 0 0 … 0 0;;; … ;;; 0 0 … 0 0; 0 0 … 0 0;;; 0 0 … 0 0; 0 0 … 0 0;;; 0 0 … 0 0; 0 0 … 0 0])
```

![Two-dimensional phase diagram of Kitaev chain model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/937f9c32-f63e-414e-bc12-92c0008d11d2)