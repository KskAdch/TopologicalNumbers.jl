@doc raw"""

Hamiltonian of the Su-Schrieffer-Heeger model.

     SSH(k, p)

 Arguments
 - k::Float64: one-dimensional wavenumber `k`.
 - p::Tuple{Float64, Float64}: parameters defined as below.


# Definition
 Hamiltonian of the Su-Schrieffer-Heeger model is defined as
```math
H(k)=
```
where,,,
"""
function SSH(k, p)
    [
        0 p[1]+p[2]*exp(-im * k)
        p[1]+p[2]*exp(im * k) 0
    ]
end

@doc raw"""

Hamiltonian of the Kitaev chain model.

     KitaevChain(k, p)

 Arguments
 - k::Float64: one-dimensional wavenumber `k`.
 - p::Tuple{Float64, Float64}: parameters defined as below.


# Definition
 Hamiltonian of the Kitaev chain model is defined as
```math
H(k)=
```
where,,,
"""
function KitaevChain(k, p)
    μ, Δ = p
    t = 1

    [
        -μ-2t*cos(k) 2im*Δ*sin(k)
        -2im*Δ*sin(k) μ+2t*cos(k)
    ]
end

@doc raw"""

Hamiltonian of the two-dimensional square lattice with flux model.

     Flux2d(k, p)

 Arguments
 - k::Float64: one-dimensional wavenumber `k`.
 - p::Tuple{Float64, Float64}: parameters defined as below.


# Definition
 Hamiltonian of the two-dimensional square lattice with flux model is defined as
```math
H(k)=
```
where,,,
"""
function Flux2d(k, p)
    k1, k2 = k
    Hsize, ν = p
    t = 1

    # Hsize = 6
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

     Haldane(k, p)

 Arguments
 - k::Float64: two-dimensional wavenumber `k`.
 - p::Tuple{Float64, Float64}: parameters defined as below.


# Definition
 Hamiltonian of the Haldane model is defined as
```math
H(k)=
```
where,,,
"""
function Haldane(k, p)
    k1, k2 = k
    J = 1.0
    K = 1.0
    ϕ, M = p

    h0 = 2K * cos(ϕ) * (cos(k1) + cos(k2) + cos(k1 + k2))
    hx = -J * (1 + cos(k1) + cos(k2))
    hy = -J * (-sin(k1) + sin(k2))
    hz = M + 2K * sin(ϕ) * (sin(k1) + sin(k2) - sin(k1 + k2))

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
end

@doc raw"""

Hamiltonian of the Kitaev honeycomb model.

     KitaevHoneycomb(k, p)

 Arguments
 - k::Float64: two-dimensional wavenumber `k`.
 - p::Tuple{Float64, Float64}: parameters defined as below.


# Definition
 Hamiltonian of the Kitaev honeycomb model is defined as
```math
H(k)=
```
where,,,
"""
function KitaevHoneycomb(k, p)
    k1, k2 = k
    K, κ = p

    hx = -κ * (sin(k2) - sin(k1) + sin(k1 - k2))
    hy = K * (sin(k1) + sin(k2))
    hz = K * (cos(k1) + cos(k2) + 1)

    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    hx .* sx .+ hy .* sy .+ hz .* sz
end

@doc raw"""

Hamiltonian of the Thouless pumping model.

     ThoulessPump(k, p)

 Arguments
 - k::Float64: two-dimensional wavenumber `k`.
 - p::Tuple{Float64, Float64}: parameters defined as below.


# Definition
 Hamiltonian of the Thouless pumping model is defined as
```math
H(k)=
```
where,,,
"""
function ThoulessPump(k, p)
end

@doc raw"""

Hamiltonian of the Kane--Mele model.

     KaneMele(k, p)

 Arguments
 - k::Float64: two-dimensional wavenumber `k`.
 - p::Tuple{Float64, Float64}: parameters defined as below.


# Definition
 Hamiltonian of the Kane--Mele model is defined as
```math
H(k)=
```
where,,,
"""
function KaneMele(k, p)
    k1, k2 = k
    t, λₛₒ = p

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


@doc raw"""

Hamiltonian of the Bernevig--Hughes--Zhang (BHZ) model.

     BHZ(k, p)

 Arguments
 - k::Float64: two-dimensional wavenumber `k`.
 - p::Tuple{Float64, Float64}: parameters defined as below.


# Definition
 Hamiltonian of the BHZ model is defined as
```math
H(k)=
```
where,,,
"""
function BHZ(k, p)
    k1, k2 = k
    tₛₚ = 1
    t₁ = ϵ₁ = 2
    ϵ₂, t₂ = p

    ϵ = -t₁*(cos(k1) + cos(k2)) + ϵ₁/2
    R1 = 0
    R2 = 0
    R3 = 2tₛₚ*sin(k2)
    R4 = 2tₛₚ*sin(k1)
    R0 = -t₂*(cos(k1) + cos(k2)) + ϵ₂/2

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
