# The Kane-Mele model

As an example of a two-dimensional topological insulator, the Kane-Mele model is presented here:

```julia
julia> function H₀(k, p) # Kane-Mele
           k1, k2 = k
           t = 1
           λₛₒ = p

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
julia> H(k) = H₀(k, 0.5)
julia> showBand(H; value=false, disp=true)
```

![Dispersion of Kane-Mele model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/2a1f8488-0e5b-4d79-be68-88bb2d744910)


Next, we can compute the $\mathbb{Z}_2$ numbers using `Z2Problem`:

```julia
julia> prob = Z2Problem(H);
julia> sol = solve(prob)
```

The output is:

```julia
Z2Solution{Vector{Int64}, Nothing, Int64}([1, 1], nothing, 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the $\mathbb{Z}_2$ number for Energy bands below and above some filling condition that you selected in the options (the default is the half-filling).
The vector is arranged in order of bands, starting from the lower energy.
The second argument `Total` stores the total of the $\mathbb{Z}_2$ numbers for Energy bands below and above some filling condition (mod2).
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
julia> param = range(-1.0, 1.0, length=1001)

julia> prob = Z2Problem(H₀);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = -1.0:0.002:1.0, nums = [1 1; 1 1; … ; 1 1; 1 1])
```

![One-dimensional phase diagram of Kane-Mele model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/779bbbb4-78c8-4599-9aba-acfe22553036)
