# The Thouless pump model

As an example of a two-dimensional topological insulator, the Thouless pump model is presented here:

```julia
julia> function H₀(k, p) # Thouless pump
           k1, t = k
           m₀, Δ = p

           R1 = -Δ*sin(k1)
           R2 = -Δ*(1-cos(k1))
           R3 = m₀*sin(t)
           R4 = -(1+cos(t))*sin(k1)
           R0 = -(1-cos(t))-(1+cos(t))*cos(k1)

           s0 = [1 0; 0 1]
           sx = [0 1; 1 0]
           sy = [0 -im; im 0]
           sz = [1 0; 0 -1]

           a1 = kron(sy, sx)
           a2 = kron(sy, sy)
           a3 = kron(sz, sz)
           a4 = kron(s0, sy)
           a0 = kron(s0, sx)

           R1 .* a1 .+ R2 .* a2 .+ R3 .* a3 .+ R4 .* a4 .+ R0 .* a0
       end
```

To calculate the dispersion, execute:

```julia
julia> H(k) = H₀(k, (-1.0, 0.5))
julia> showBand(H; value=false, disp=true)
```

![Dispersion of Thouless pump model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/f159045f-4437-455f-8e78-ba989c95a28e)


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
julia> H(k, p) = H₀(k, (p, 0.25));
julia> param = range(-1.0, 0.0, length=1001)

julia> prob = Z2Problem(H);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = -1.0:0.001:0.0, nums = [1 1; 1 1; … ; 0 0; 0 0])
```

![One-dimensional phase diagram of Thouless pump model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/f4524c76-f252-45d8-bc41-655cca5112a3)


Also, two-dimensional phase diagram is given by:

```julia
julia> param1 = range(-1.0, -0.1, length=91);
julia> param2 = range(-0.5, 0.5, length=101);

julia> prob = Z2Problem(H₀);
julia> calcPhaseDiagram(prob, param1, param2; plot=true)
```


![Two-dimensional phase diagram of Thouless pump model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/3c9327b0-3761-46b5-b4b0-62390a0aa40c)
