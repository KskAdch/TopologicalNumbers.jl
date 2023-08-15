# The Bernevig-Hughes-Zhang (BHZ) model

As an example of a two-dimensional topological insulator, the BHZ model is presented here:

```julia
julia> function H(k) # BHZ
    k1, k2 = k
    t = 1

    R0 = -2(cos(k1) + cos(k2)) + 1
    R3 = 2sin(k2)
    R4 = 2sin(k1)
    R5 = -2t*(cos(k1) + cos(k2)) + 1

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    a0 = kron(s0, s0)
    a1 = kron(sx, sx)
    a2 = kron(sx, sy)
    a3 = kron(sx, sz)
    a4 = kron(sy, s0)
    a5 = kron(sz, s0)

    R0*a0+R3*a3+R4*a4+R5*a5
end
```

To calculate the dispersion, execute:

```julia
julia> showBand(H; value=false, disp=true)
```

![Dispersion of BHZ model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/a9cf9768-6920-45e6-89bd-ed7ec434152c)


Next, we can compute the $\mathbb{Z}_2$ numbers using `calcZ2`:

```julia
julia> calcZ2(H)
```

The output is:

```julia
(TopologicalNumber = [1, 1], Total = 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the $\mathbb{Z}_2$ number for each each pair of two energy bands. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the $\mathbb{Z}_2$ numbers for each pair of two energy bands.
`Total` is a quantity that should always return zero.

`Total` is a value that should consistently return zero.