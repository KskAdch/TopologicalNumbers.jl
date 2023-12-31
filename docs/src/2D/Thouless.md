# The Thouless pump model

<!-- As an example of a two-dimensional topological insulator, the Kane-Mele model is presented here: -->

```julia
julia> function H₀(k, p) # Thouless pump
           k, t = k
           m₀, Δ = p

           R1 = -Δ*sin(k)
           R2 = -Δ*(1-cos(k))
           R3 = m₀*sin(t)
           R4 = -(1+cos(t))*sin(k)
           R0 = -(1-cos(t))-(1+cos(t))*cos(k)

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

<!-- To calculate the dispersion, execute:

```julia
julia> H(k) = H₀(k, 0.5)
julia> showBand(H; value=false, disp=true)
```

![Dispersion of Kane-Mele model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/2a1f8488-0e5b-4d79-be68-88bb2d744910)


Next, we can compute the $\mathbb{Z}_2$ numbers using `calcZ2`:

```julia
julia> calcZ2(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the $\mathbb{Z}_2$ number for each pair of two energy bands. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the $\mathbb{Z}_2$ numbers for each pair of two energy bands.
`Total` is a quantity that should always return zero.


One-dimensional phase diagram is given by:

```julia
julia> param = range(-1.0, 1.0, length=1001)
julia> calcPhaseDiagram(H₀, param, "Z2"; plot=true)
```

![One-dimensional phase diagram of Kane-Mele model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/779bbbb4-78c8-4599-9aba-acfe22553036) -->


<!-- Also, two-dimensional phase diagram is given by:

```julia
julia> param = range(-1.0, 1.0, length=101)
julia> calcPhaseDiagram(H₀, param, param, "Z2"; plot=true)
```


![Two-dimensional phase diagram of Kane-Mele model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/4fe6d699-e3f9-45cd-b176-b4f1216d39d8) -->