@doc raw"""

Hamiltonian of the Su-Schrieffer-Heeger model.

     SSH(k::T1, p::T2) where {T1<:Real,T2<:Real}

# Arguments
 - `k::T1`: one-dimensional wavenumber `k`.
 - `p::T2`: parameter defined as below.


# Definition
 Hamiltonian of the Su-Schrieffer-Heeger model is defined as
```math
H(k)=\begin{pmatrix}
    0                 & t_{1}+t_{2}e^{-ik} \\
    t_{1}+t_{2}e^{ik} & 0
\end{pmatrix}
```
where $t_{1}$ and $t_{2}$ are the amplitudes of the nearest neighbor hopping in the tight binding model.
Nondimensionalize with $t_{1}=1$ and $t_{2}=$`p`.

# Examples

```julia
julia> SSH(π/3, 0.5)
2×2 Matrix{ComplexF64}:
  0.0+0.0im       1.25-0.433013im
 1.25+0.433013im   0.0+0.0im
```
"""
function SSH(k::T1, p::T2) where {T1<:Real,T2<:Real}
    t₁ = 1
    t₂ = p

    [
        0 t₁+t₂*exp(-im * k)
        t₁+t₂*exp(im * k) 0
    ]
end

@doc raw"""

Hamiltonian of the Kitaev chain model.

     KitaevChain(k::T1, p::T2) where {T1<:Real, T2<:Union{AbstractVector,Tuple}}

# Arguments
 - `k::T1`: one-dimensional wavenumber `k`.
 - `p::T2`: parameters defined as below.


# Definition
 Hamiltonian of the Kitaev chain model is defined as
```math
H(k)=\begin{pmatrix}
    -\mu-2t\cos(k)   & 2i\Delta\sin(k) \\
    -2i\Delta\sin(k) & \mu+2t\cos(k)
\end{pmatrix}
```
where $t$ is the amplitude of the nearest neighbor hopping in the tight binding model.
$\Delta$ is the pairing with the nearest neighbor site, $\mu$ is the chemical potential.
Nondimensionalize with $t=1$ and $\mu,\Delta=$`p`.

# Examples

```julia
julia> SSH(π/3, 0.5)
2×2 Matrix{ComplexF64}:
 -0.5+0.0im       0.0+0.866025im
  0.0-0.866025im  0.5+0.0im
```
"""
function KitaevChain(k::T1, p::T2) where {T1<:Real,T2<:Union{AbstractVector,Tuple}}
    μ, Δ = p
    t = 1

    [
        -μ-2t*cos(k) 2im*Δ*sin(k)
        -2im*Δ*sin(k) μ+2t*cos(k)
    ]
end

@doc raw"""

Hamiltonian of the two-dimensional square lattice with flux model.

     Flux2d(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple}, T2<:Union{AbstractVector,Tuple}}

# Arguments
 - `k::T1`: one-dimensional wavenumber `k`.
 - `p::T2`: parameters defined as below.


# Definition
 Hamiltonian of the two-dimensional square lattice with flux model is defined as
```math
H(k)=\begin{pmatrix}
    -2t\cos(k_{2}-2\pi m\frac{1}{n}) & -t                               & 0      & 0  & \cdots                               & e^{-ik_{1}}           \\
    -t                               & -2t\cos(k_{2}-2\pi m\frac{2}{n}) & -t     & 0  & \cdots                               & 0                     \\
    0                                & -t                               & \ddots & -t & \cdots                               & \vdots                \\
    \vdots                           & \vdots                           &        & -t & -2t\cos(k_{2}-2\pi m\frac{(n-1)}{n}) & -t                    \\
    e^{ik_{1}}                       & 0                                & \cdots & 0  & -t                                   & -2t\cos(k_{2}-2\pi m)
\end{pmatrix}
```
where $t$ is the amplitude of the nearest neighbor hopping in the tight binding model.
$n$ and $m$ are the size of the unit cell of the phase produced by the application of the static magnetic field.
Nondimensionalize with $t=1$ and $n,m=$`p`.

# Examples

```julia
julia> Flux2d([0, π/3], [4, 1])
4×4 Matrix{ComplexF64}:
 -1.73205+0.0im  -1.0+0.0im      0.0+0.0im  -1.0+0.0im
     -1.0+0.0im   1.0+0.0im     -1.0+0.0im   0.0+0.0im
      0.0+0.0im  -1.0+0.0im  1.73205+0.0im  -1.0+0.0im
     -1.0-0.0im   0.0+0.0im     -1.0+0.0im  -1.0+0.0im
```
"""
function Flux2d(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple},T2<:Union{AbstractVector,Tuple}}
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

@doc raw"""

Hamiltonian of the Haldane model.

     Haldane(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple}, T2<:Union{AbstractVector,Tuple}}

# Arguments
 - `k::T1`: two-dimensional wavenumber `k`.
 - `p::T2`: parameters defined as below.


# Definition
 Hamiltonian of the Haldane model is defined as
```math
H(\bm{k})=\begin{pmatrix}
    h_{0}(\bm{k})+h_{3}(\bm{k})  & h_{1}(\bm{k})-ih_{2}(\bm{k}) \\
    h_{1}(\bm{k})+ih_{2}(\bm{k}) & h_{0}(\bm{k})-h_{3}(\bm{k})
\end{pmatrix}
```

```math
\begin{align*}
     h_{0}(\bm{k})=&2t_{2}\cos(\phi)(\sin(k_{1})+\sin(k_{2})+\sin(k_{1}+k_{2})) \\
     h_{1}(\bm{k})=&-t_{1}(1+\cos(k_{1})+\cos(k_{2})) \\
     h_{2}(\bm{k})=&-t_{1}(-\sin(k_{1})+\sin(k_{2})) \\
     h_{3}(\bm{k})=&m+2t_{2}\sin(\phi)(\sin(k_{1})+\sin(k_{2})-\sin(k_{1}+k_{2}))
\end{align*}
```
where $t_{1}$ and $t_{2}$ are the amplitudes of the nearest and the next nearest neighbor hopping in the tight binding model.
$\phi$ is the phase produced by the application of the static magnetic field.
Nondimensionalize with $t_{1}=1$ and $t_{2},\phi,m=$`p`.

# Examples

```julia
julia> Haldane([0, π/3], [0.5, π/2, 0.3])
2×2 Matrix{ComplexF64}:
  0.3-0.0im       -2.5+0.866025im
 -2.5-0.866025im  -0.3-0.0im
```
"""
function Haldane(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple}, T2<:Union{AbstractVector,Tuple}}
    k1, k2 = k
    t₁ = 1
    t₂, ϕ, m = p

    h0 = 2t₂ * cos(ϕ) * (cos(k1) + cos(k2) + cos(k1 + k2))
    hx = -t₁ * (1 + cos(k1) + cos(k2))
    hy = -t₁ * (-sin(k1) + sin(k2))
    hz = m + 2t₂ * sin(ϕ) * (sin(k1) + sin(k2) - sin(k1 + k2))

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
end

@doc raw"""

Hamiltonian of the Kitaev honeycomb model.

     KitaevHoneycomb(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple}, T2<:Real}

# Arguments
 - `k::T1`: two-dimensional wavenumber `k`.
 - `p::T2`: parameter defined as below.


# Definition
 Hamiltonian of the Kitaev honeycomb model is defined as
```math
H(\bm{k})=\begin{pmatrix}
    h_{3}(\bm{k})                & h_{1}(\bm{k})-ih_{2}(\bm{k}) \\
    h_{1}(\bm{k})+ih_{2}(\bm{k}) & -h_{3}(\bm{k})
\end{pmatrix}
```

```math
\begin{align*}
     h_{1}(\bm{k})=&-K_{2}(\sin(k_{2})-\sin(k_{1})+\sin(k_{1}-k_{2})) \\
     h_{2}(\bm{k})=&-K_{1}(\sin(k_{1})+\sin(k_{2})) \\
     h_{3}(\bm{k})=&K_{1}(\cos(k_{1})+\cos(k_{2})+1)
\end{align*}
```
where $K_{1}$ is the magnitude of Kitaev interaction.
$K_{2}$ is the magnitude of spin triples term due to magnetic field.
Nondimensionalize with $K_{1}=1$ and $K_{2}=$`p`.

# Examples

```julia
julia> KitaevHoneycomb([0, π/3], 0.5)
2×2 Matrix{ComplexF64}:
 2.5-0.0im        0.0+0.866025im
 0.0-0.866025im  -2.5-0.0im
```
"""
function KitaevHoneycomb(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple},T2<:Real}
    k1, k2 = k
    K₁ = 1
    K₂ = p

    hx = -K₂ * (sin(k2) - sin(k1) + sin(k1 - k2))
    hy = -K₁ * (sin(k1) + sin(k2))
    hz = K₁ * (cos(k1) + cos(k2) + 1)

    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    hx .* sx .+ hy .* sy .+ hz .* sz
end

@doc raw"""

Hamiltonian of the Thouless pumping model.

     ThoulessPump(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple}, T2<:Union{AbstractVector,Tuple}}

# Arguments
 - `k::T1`: two-dimensional wavenumber `k`.
 - `p::T2`: parameters defined as below.


# Definition
 Hamiltonian of the Thouless pumping model is defined as
```math
H(k)=
```
where,,,
"""
function ThoulessPump(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple},T2<:Union{AbstractVector,Tuple}}
end

@doc raw"""

Hamiltonian of the Kane--Mele model.

     KaneMele(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple}, T2<:Real}

# Arguments
 - `k::T1`: two-dimensional wavenumber `k`.
 - `p::T2`: parameter defined as below.


# Definition
 Hamiltonian of the Kane--Mele model is defined as
```math
H(k)=\begin{pmatrix}
    h_{3}(\bm{k})                & 0                            & h_{5}(\bm{k})-ih_{4}(\bm{k}) & 0                            \\
    0                            & -h_{3}(\bm{k})               & 0                            & h_{5}(\bm{k})-ih_{4}(\bm{k}) \\
    h_{5}(\bm{k})+ih_{4}(\bm{k}) & 0                            & -h_{3}(\bm{k})               & 0                            \\
    0                            & h_{5}(\bm{k})+ih_{4}(\bm{k}) & 0                            & h_{3}(\bm{k})
\end{pmatrix}
```

```math
\begin{align*}
     h_{3}(\bm{k})=&2\lambda_{\mathrm{SO}}(\sin(k_{1})-\sin(k_{2})-\sin(k_{1}-k_{2})) \\
     h_{4}(\bm{k})=&-t(\sin(k_{1})+\sin(k_{2})) \\
     h_{5}(\bm{k})=&-t(\cos(k_{1})+\cos(k_{2})+1)
\end{align*}
```
where $t$ is the amplitude of the nearest neighbor hopping in the tight binding model.
$\lambda_{\mathrm{SO}}$ is the magnitude of spin-orbit interaction.
Nondimensionalize with $t=1$ and $\lambda_{\mathrm{SO}}=$`p`.

# Examples

```julia
julia> KaneMele([0, π/3], 0.5)
4×4 Matrix{ComplexF64}:
  0.0+0.0im        0.0+0.0im       -2.5+0.866025im   0.0+0.0im
  0.0+0.0im        0.0+0.0im        0.0+0.0im       -2.5+0.866025im
 -2.5-0.866025im   0.0+0.0im        0.0+0.0im        0.0+0.0im
  0.0+0.0im       -2.5-0.866025im   0.0+0.0im        0.0+0.0im
```
"""
function KaneMele(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple},T2<:Real}
    k1, k2 = k
    t = 1
    λₛₒ = p

    R1 = 0
    R2 = 0
    R3 = 2λₛₒ * (sin(k1) - sin(k2) - sin(k1 - k2))
    R4 = -t * (sin(k1) + sin(k2))
    R0 = -t * (cos(k1) + cos(k2) + 1)

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


@doc raw"""

Hamiltonian of the Bernevig--Hughes--Zhang (BHZ) model.

     BHZ(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple}, T2<:Union{AbstractVector,Tuple}}

# Arguments
 - `k::T1`: two-dimensional wavenumber `k`.
 - `p::T2`: parameters defined as below.


# Definition
 Hamiltonian of the BHZ model is defined as
```math
H(k)=\begin{pmatrix}
    h_{0}(\bm{k})+h_{3}(\bm{k})  & 0                            & h_{5}(\bm{k})-ih_{4}(\bm{k}) & 0                            \\
    0                            & h_{0}(\bm{k})-h_{3}(\bm{k})  & 0                            & h_{5}(\bm{k})-ih_{4}(\bm{k}) \\
    h_{5}(\bm{k})+ih_{4}(\bm{k}) & 0                            & h_{0}(\bm{k})-h_{3}(\bm{k})  & 0                            \\
    0                            & h_{5}(\bm{k})+ih_{4}(\bm{k}) & 0                            & h_{0}(\bm{k})+h_{3}(\bm{k})
\end{pmatrix}
```

```math
\begin{align*}
     h_{0}(\bm{k})=&-(t_{ss}-t_{pp})(\cos(k_{1})+\cos(k_{2}))+\frac{(\epsilon_{s}+\epsilon_{p})}{2} \\
     h_{3}(\bm{k})=&2t_{sp}\sin(k_{2}) \\
     h_{4}(\bm{k})=&2t_{sp}\sin(k_{1}) \\
     h_{5}(\bm{k})=&-(t_{ss}+t_{pp})(\cos(k_{1})+\cos(k_{2}))+\frac{(\epsilon_{s}-\epsilon_{p})}{2}
\end{align*}
```
where $t_{sp}$ is the amplitude of the mixed hopping between s- and p-orbitals in the tight binding model.
$t_{ss}$ is the hopping between s-orbitals, $t_{pp}$ is the hopping between p-orbitals, and $t_{1}=t_{s}-t_{p}$, $t_{2}=t_{s}+t_{p}$.
$\epsilon_{s}$ and $\epsilon_{p}$ are the site potential terms for the s- and p-orbitals, $\epsilon_{1}\epsilon_{s}s+\epsilon_{p}$, $\epsilon_{2}=\epsilon_{s}-\epsilon_{p}.
Nondimensionalize with $t_{sp}=t_{1}=\epsilon_{1}=1$ and $\epsilon_{2},t_{2}=$`p`.

# Examples

```julia
julia> BHZ([0, π/3], [0.5, 0.5])
4×4 Matrix{ComplexF64}:
 0.732051+0.0im       0.0+0.0im      -0.5+0.0im       0.0+0.0im
      0.0+0.0im  -2.73205+0.0im       0.0+0.0im      -0.5+0.0im
     -0.5+0.0im       0.0+0.0im  -2.73205+0.0im       0.0+0.0im
      0.0+0.0im      -0.5+0.0im       0.0+0.0im  0.732051+0.0im
"""
function BHZ(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple},T2<:Union{AbstractVector,Tuple}}
    k1, k2 = k
    tₛₚ = 1
    t₁ = ϵ₁ = 1
    ϵ₂, t₂ = p

    ϵ = -t₁ * (cos(k1) + cos(k2)) + ϵ₁ / 2
    R1 = 0
    R2 = 0
    R3 = 2tₛₚ * sin(k2)
    R4 = 2tₛₚ * sin(k1)
    R0 = -t₂ * (cos(k1) + cos(k2)) + ϵ₂ / 2

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    I0 = Matrix{Int64}(I, 4, 4)
    a1 = kron(sz, sx)
    a2 = kron(sz, sy)
    a3 = kron(sz, sz)
    a4 = kron(sy, s0)
    a0 = kron(sx, s0)

    ϵ .* I0 .+ R1 .* a1 .+ R2 .* a2 .+ R3 .* a3 .+ R4 .* a4 .+ R0 .* a0
end


@doc raw"""

Hamiltonian of the lattice Dirac model.

     LatticeDirac(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple},T2<:Real}

# Arguments
 - `k::T1`: four-dimensional wavenumber `k`.
 - `p::T2`: parameter defined as below.


# Definition
 Hamiltonian of the lattice Dirac model is defined as
```math
H(k)=\begin{pmatrix}
    h_{4}(\bm{k})                & h_{2}(\bm{k})-ih_{3}(\bm{k}) & h_{0}(\bm{k})-ih_{1}(\bm{k})  & 0                             \\
    h_{2}(\bm{k})+ih_{3}(\bm{k}) & -h_{4}(\bm{k})               & 0                             & h_{0}(\bm{k})-ih_{1}(\bm{k})  \\
    h_{0}(\bm{k})+ih_{1}(\bm{k}) & 0                            & -h_{4}(\bm{k})                & -h_{2}(\bm{k})+ih_{3}(\bm{k}) \\
    0                            & h_{0}(\bm{k})+ih_{1}(\bm{k}) & -h_{2}(\bm{k})-ih_{3}(\bm{k}) & h_{4}(\bm{k})
\end{pmatrix}
```

```math
\begin{align*}
     h_{0}(\bm{k})=&m+c(cos(k_{1})+cos(k_{2})+cos(k_{3})+cos(k_{4})) \\
     h_{1}(\bm{k})=&sin(k_{1}) \\
     h_{2}(\bm{k})=&sin(k_{2}) \\
     h_{3}(\bm{k})=&sin(k_{3}) \\
     h_{4}(\bm{k})=&sin(k_{4})
\end{align*}
```
where,,,

# Examples

```julia
julia> LatticeDirac([0, 0, 0, π/2], 0.5)
4×4 Matrix{ComplexF64}:
 1.0+0.0im   0.0+0.0im   3.5+0.0im  0.0+0.0im
 0.0+0.0im  -1.0+0.0im   0.0+0.0im  3.5+0.0im
 3.5+0.0im   0.0+0.0im  -1.0+0.0im  0.0+0.0im
 0.0+0.0im   3.5+0.0im   0.0+0.0im  1.0+0.0im
"""
function LatticeDirac(k::T1, p::T2) where {T1<:Union{AbstractVector,Tuple},T2<:Real}
    k1, k2, k3, k4 = k
    c = 1
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

    h1 = m + c*(cos(k1) + cos(k2) + cos(k3) + cos(k4))
    h2 = sin(k1)
    h3 = sin(k2)
    h4 = sin(k3)
    h5 = sin(k4)

    # Update the Hamiltonian matrix in place
    h1 .* g1 .+ h2 .* g2 .+ h3 .* g3 .+ h4 .* g4 .+ h5 .* g5
end