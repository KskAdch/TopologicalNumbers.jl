# Two-dimensional square lattice model with flux

A two-dimensional example is presented here:

```julia
julia> function H₀(k, p)
           k1, k2 = k
           Hsize, ν = p
           t = 1

           Hmat = zeros(ComplexF64, Hsize, Hsize)

           ϕ = 2π * ν / Hsize

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
You can also use our preset Hamiltonian function `Flux2d` to define the same Hamiltonian matrix as follows:

```julia
julia> H₀(k, p) = Flux2d(k, p)
```

To calculate the dispersion, run:

```julia
julia> H(k) = H₀(k, (6, 1))
julia> showBand(H; value=false, disp=true)
```

![Dispersion of 2D square lattice with flux model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/630d8fc6-e1ee-4f0c-855e-80ee1b8c115f)


Then we can compute the Chern numbers using `FCProblem`:

```julia
julia> prob = FCProblem(H);
julia> sol = solve(prob)
```

The output is:

```julia
FCSolution{Vector{Int64}, Int64}([1, 1, -2, -2, 1, 1], 0)
```

The first argument `TopologicalNumber` in the named tuple is an vector that stores the first Chern number for each band. 
The vector is arranged in order of bands, starting from the one with the lowest energy.
The second argument `Total` stores the total of the first Chern numbers for each band.
`Total` is a quantity that should always return zero.

You can access these values as follows:

```julia
julia> sol.TopologicalNumber
6-element Vector{Int64}:
  1
  1
 -2
 -2
  1
  1

julia> sol.Total
0
```


One-dimensional phase diagram is given by:

```julia
julia> H(k, p) = H₀(k, (6, p));
julia> param = 0:6;

julia> prob = FCProblem(H);
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = 0:6, nums = [0 0 … 0 0; 1 1 … 1 1; … ; -1 -1 … -1 -1; 0 0 … 0 0])
```

![One-dimensional phase diagram](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139373570/42f0532e-03b5-4d4f-a8e1-4777a9777d13)
