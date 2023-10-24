# Weyl semimetal

A three-dimensional example is presented here:

```julia
julia> function H₀(k) # Weyl
    k1, k2, k3 = k
    t1 = 1
    t2 = 1
    t3 = 1
    m = 2
    k0 = 2pi*2/5

    h0 = 0
    hx = 2t1*(cos(k1) - cos(k0)) + m*(2 - cos(k2) - cos(k3))
    hy = 2t2*sin(k2)
    hz = 2t3*sin(k3)

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
end
```
