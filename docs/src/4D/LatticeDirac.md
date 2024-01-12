# Lattice Dirac model

Hamiltonian of lattice Dirac model is given by:

```julia
julia> function H₀(k, p) # lattice Dirac
            k1, k2, k3, k4 = k
            m = p

            # Define Pauli matrices and Gamma matrices
            σ₀ = [1 0; 0 1]
            σ₁ = [0 1; 1 0]
            σ₂ = [0 -im; im 0]
            σ₃ = [1 0; 0 -1]
            g1 = kron(σ₁, σ₀)
            g2 = kron(σ₂, σ₀)
            g3 = kron(σ₃, σ₁)
            g4 = kron(σ₃, σ₂)
            g5 = kron(σ₃, σ₃)

            h1 = m + cos(k1) + cos(k2) + cos(k3) + cos(k4)
            h2 = sin(k1)
            h3 = sin(k2)
            h4 = sin(k3)
            h5 = sin(k4)

            # Return the Hamiltonian matrix
            h1 .* g1 .+ h2 .* g2 .+ h3 .* g3 .+ h4 .* g4 .+ h5 .* g5
        end
```
You can also use our preset Hamiltonian function `LatticeDirac` to define the same Hamiltonian matrix as follows:

```julia
julia> H₀(k, p) = LatticeDirac(k, p)
```

Then we can compute the second Chern numbers using `SCProblem`:

```julia
julia> H(k) = H₀(k, -3.0)

julia> prob = SCProblem(H);
julia> sol = solve(prob)
```

The output is:

```julia
SCSolution{Float64}(0.9793607631927376)
```

The argument `TopologicalNumber` in the named tuple stores the second Chern number with some filling condition that you selected in the options (the default is the half-filling). 

You can access this value as follows:

```julia
julia> sol.TopologicalNumber
0.9793607631927376
```


A phase diagram is given by:

```julia
julia> param = range(-4.9, 4.9, length=10);

julia> N = 30; # number of k-points in each direction
julia> prob = SCProblem(H₀, N); # N is optional argument
julia> sol = calcPhaseDiagram(prob, param; plot=true)
(param = -4.9:1.0888888888888888:4.9, nums = [0.00024580085568788514, 0.8920579621583358, 0.9779212824560908, -2.8575405041314244, -2.915655604396968, 2.9165403378212695, 2.8604644187734776, -0.9777674289766198, -0.886041497183358, -0.0002442110189681556])
```

![Phase diagram of lattice Dirac model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/4e967ed2-9011-4e88-85ab-3b64deaaf09a)

A more dense diagram can also be obtained:

```julia
julia> param = range(-4.9, 4.9, length=50)
julia> sol = calcPhaseDiagram(prob, param; plot=true)
```

![Dense phase diagram of lattice Dirac model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/9449651c-46f4-4141-ac87-161e9e5fbf28)

Since the system dimension is high,
a computational cost is high comparing with other low dimensional cases.



If you want to use a parallel environment, you can utilize `MPI.jl`.
Let's create a file named `test.jl` with the following content:
```julia
using TopologicalNumbers
using MPI

H₀(k, p) = LatticeDirac(k, p)
H(k) = H₀(k, -3.0)

param = range(-4.9, 4.9, length=10)

prob = SCProblem(H₀)
result = calcPhaseDiagram(prob, param; parallel=UseMPI(MPI), progress=true)

plot1D(result; labels=true, disp=false, pdf=true)
```
You can perform calculations using `mpirun` (for example, with `4` cores) as follows:
```bash
mpirun -np 4 julia --project test.jl
```
For more details, refer to the `MPI.jl` document at [https://juliaparallel.org/MPI.jl/stable/](https://juliaparallel.org/MPI.jl/stable/).