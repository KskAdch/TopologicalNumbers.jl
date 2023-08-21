# Two-dimensional square lattice model with flux

A two-dimensional example is presented here:

```julia
julia> function H₀(k, p) # landau
    k1, k2 = k
    t = 1

    Hsize = 6
    Hmat = zeros(ComplexF64, Hsize, Hsize)

    ϕ = 2π * p / Hsize

    for i in 1:Hsize
        Hmat[i, i] = -2t * cos(k2 - i * ϕ)
    end

    for i in 1:Hsize-1
        Hmat[i, i+1] = -t
        Hmat[i+1, i] = -t
    end

    Hmat[1, Hsize] = -t * exp(-im * k1)
    Hmat[Hsize, 1] = -t * exp(im * k1)

    Hmat
end
```

To calculate the dispersion, run:

```julia
julia> H(k) = H₀(k, 1)
julia> showBand(H; value=false, disp=true)
```

![Dispersion of 2D square lattice with flux model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/8470121f-bdad-4960-9848-7ade1ae805d3)


Then we can compute the Chern numbers using `calcChern`:

```julia
julia> calcChern(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1, -2, -2, 1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.


One-dimensional phase diagram is given by:

```julia
julia> param = 1:6
julia> calcPhaseDiagram(H₀, param, "Chern"; plot=true)
```
%%%dummy%%%
![One-dimensional phase diagram](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/f9c179c3-1275-4640-ac21-0d10737fcaf7)
%%%dummy%%%
