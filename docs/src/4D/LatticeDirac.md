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

            # Update the Hamiltonian matrix in place
            h1 .* g1 .+ h2 .* g2 .+ h3 .* g3 .+ h4 .* g4 .+ h5 .* g5
        end
```

Then we can compute the second Chern numbers using `calcSecondChern`:

```julia
julia> H(k) = H₀(k, -3.0)
julia> calcSecondChern(H)
```

The output is:

```julia
(TopologicalNumber = 0.9793607631927379,)
```

The argument `TopologicalNumber` in the named tuple stores the second Chern number with some filling condition that you selected in options (the default is the half-filling). 


Phase diagram is given by:

```julia
julia> param = range(-4.9, 4.9, length=10)
julia> calcPhaseDiagram(H₀, param, SecondChern_FHS(); N=30, plot=true)
```

![Phase diagram of lattice Dirac model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/4e967ed2-9011-4e88-85ab-3b64deaaf09a)

A more dense diagram can also be obtained:

```julia
julia> param = range(-4.9, 4.9, length=50)
julia> calcPhaseDiagram(H₀, param, SecondChern_FHS(); N=30, plot=true)
```

![Dense phase diagram of lattice Dirac model](https://github.com/KskAdch/TopologicalNumbers.jl/assets/139110206/9449651c-46f4-4141-ac87-161e9e5fbf28)

Since the system dimension is high,
a computational cost is high comparing with other low dimensional cases.



If you want to use a parallel environment, you can utilize MPI.jl.
Let's create a file named test.jl with the following content:
```julia
using TopologicalNumbers
using MPI

H₀(k, p) = LatticeDirac(k, p)
H(k) = H₀(k, -3.0)

param = range(-4.9, 4.9, length=10)
result = calcPhaseDiagram(H₀, param, SecondChern_FHS(); N=30, parallel=UseMPI(MPI), progress=true)
plot1D(result; labels=true, disp=false, pdf=true)
```
You can perform calculations using `mpirun` (for example, with `4` cores) as follows:
```bash
mpirun -np 4 julia --project test.jl
```
For more details, refer to the `MPI.jl` document at [https://juliaparallel.org/MPI.jl/stable/](https://juliaparallel.org/MPI.jl/stable/).