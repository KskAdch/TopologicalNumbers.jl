abstract type TopologicalNumbersProblems end

# Problem for calculating the Berry phase
@doc raw"""
    struct BPProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Real,T4<:Bool} <: TopologicalNumbersProblems

The `BPProblem` struct represents a problem for calculating Berry phase.

# Fields
- `H::T1`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 1.
- `N::T2`: The number of points for one direction in the Brillouin zone. Default is 51.
- `gapless::T3`: The threshold for considering a band as gapless. Default is 0.0.
- `rounds::T4`: A boolean indicating whether to round a returned variable. Default is true.

# Example
```julia
julia> 
```
"""
@kwdef struct BPProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Real,T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 51
    gapless::T3 = 0.0
    rounds::T4 = true
end

# default
@doc raw"""
    BPProblem(H)

Constructs a Berry phase problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 1.

# Returns
A `BPProblem` object.

# Example
```julia
julia> 
```

"""
BPProblem(H) = BPProblem(; H=H)

@doc raw"""
    BPProblem(H, N)

Constructs a Berry phase problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 1.
- `N`: The number of points for one direction in the Brillouin zone.

# Returns
A `BPProblem` object.

# Example
```julia
julia> 
```

"""
BPProblem(H, N) = BPProblem(; H=H, N=N)


# Problem for calculating the first Chern number
@doc raw"""
    FCProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Real,T4<:Bool} <: TopologicalNumbersProblems

A struct representing a problem for calculating the first Chern number.

# Fields
- `H::T1`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 2.
- `N::T2`: The number of points for one direction in the Brillouin zone. Default is 51.
- `gapless::T3`: The threshold for considering a band as gapless. Default is 0.0.
- `rounds::T4`: A boolean indicating whether to round a returned variable. Default is true.

# Example
```julia
julia> 
```
"""
@kwdef struct FCProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Real,T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 51
    gapless::T3 = 0.0
    rounds::T4 = true
end

# default
@doc raw"""
    FCProblem(H)

Constructs a first Chern number problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 2.

# Returns
A `FCProblem` object.

# Example
```julia
julia> 
```
"""
FCProblem(H) = FCProblem(; H=H)

@doc raw"""
    FCProblem(H, N)

Constructs a first Chern number problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 2.
- `N`: The number of points for one direction in the Brillouin zone.

# Returns
A `FCProblem` object.

# Example
```julia
julia> 
```
"""
FCProblem(H, N) = FCProblem(; H=H, N=N)


# Problem for calculating the second Chern number
@doc raw"""
    SCProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Union{Real,Nothing},T4<:Bool} <: TopologicalNumbersProblems

A struct representing a problem for calculating the second Chern number.

# Fields
- `H::T1`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 4.
- `N::T2`: The number of points in the Brillouin zone. Type of `Int` means the uniform mesh. You can also specify the mesh by giving a tuple or a vector. Default is 30.
- `Nfill::T3`: The number of filled bands. `nothing` means the half-filling. Default is `nothing`.
- `RV::T4`: A boolean indicating whether to return a `real` value. Default is true.

# Example
```julia
julia> 
```
"""
@kwdef struct SCProblem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Union{Real,Nothing},T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 30
    Nfill::T3 = nothing
    RV::T4 = true
end

# default

@doc raw"""
    SCProblem(H)

Constructs a second Chern number problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 4.

# Returns
A `SCProblem` object.

# Example
```julia
julia> 
```
"""
SCProblem(H) = SCProblem(; H=H)

@doc raw"""
    SCProblem(H, N)

Constructs a second Chern number problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 4.
- `N`: The number of points in the Brillouin zone. Type of `Int` means the uniform mesh. You can also specify the mesh by giving a tuple or a vector.

# Returns
A `SCProblem` object.

# Example
```julia
julia> 
```
"""
SCProblem(H, N) = SCProblem(; H=H, N=N)


# Problem for calculating the Z2 invariant
@doc raw"""
    Z2Problem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Bool} <: TopologicalNumbersProblems

A struct representing a problem for calculating the Z2 number.

# Fields
- `H::T1`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 2.
- `N::T2`: The number of points for one direction in the Brillouin zone. Default is 50.
- `rounds::T3`: A boolean indicating whether to round a returned variable. Default is true.
- `TR::T3`: A boolean indicating whether to calculate the remaining part of the Brillouin zone. If `true`, `solve` returns an additional result `TRTopologicalNumber`. If the calculation is done nomally, `TRTopologicalNumber` is equal to `TopologicalNumber`. Default is false.

# Example
```julia
julia> 
```
"""
@kwdef struct Z2Problem{T1<:Function,T2<:Union{Tuple,AbstractVector,Int},T3<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 50
    rounds::T3 = true
    TR::T3 = false
end
# default

@doc raw"""
Z2Problem(H)

Constructs a Z2 number problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 2.

# Returns
A `Z2Problem` object.

# Example
```julia
julia> 
```
"""
Z2Problem(H) = Z2Problem(; H=H)

@doc raw"""
    Z2Problem(H, N)

Constructs a Z2 number problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k, p)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector and `p` contains parameter. Dimension of `k` must be 2.
- `N`: The number of points for one direction in the Brillouin zone.

# Returns
A `Z2Problem` object.

# Example
```julia
julia> 
```
"""
Z2Problem(H, N) = Z2Problem(; H=H, N=N)


# Problem for calculating the local Berry flux
@doc raw"""
    LBFProblem{T1<:Function,T2<:AbstractVector,T3<:Union{Tuple,AbstractVector,Int},T4<:Real,T5<:Bool} <: TopologicalNumbersProblems

A struct representing a problem for calculating the $k$-local value of Berry flux.

# Fields
- `H::T1`: The Hamiltonian function `H=H(k)` that defines the system. `k` is an `AbstractVector` (or a `Tuple`) of the wavenumber vector. Dimension of `k` must be 2.
- `n::T2`: An `AbstractVector` (or a `Tuple`) including two elements of `Int`, which represents wavenumber ($2\pi n/N$) when calculating Berry flux. Dimension of `n` must be 2.
- `N::T3`: The number of points for one direction in the Brillouin zone. Default is 51.
- `gapless::T4`: The threshold for considering a band as gapless. Default is 0.0.
- `rounds::T5`: A boolean indicating whether to round a returned variable. Default is true.

# Example
```julia
julia> 
```
"""
@kwdef struct LBFProblem{T1<:Function,T2<:AbstractVector,T3<:Union{Tuple,AbstractVector,Int},T4<:Real,T5<:Bool} <: TopologicalNumbersProblems
    H::T1
    n::T2
    N::T3 = 51
    gapless::T4 = 0.0
    rounds::T5 = true
end
# default

@doc raw"""
    LBFProblem(H, n)

Constructs a local Berry flux problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 2.
- `n`: An `AbstractVector` (or a `Tuple`) including two elements of `Int`, which represents wavenumber ($2\pi n/N$) when calculating Berry flux. Dimension of `n` must be 2.

# Returns
A `LBFProblem` object.

# Example
```julia
julia> 
```
"""
LBFProblem(H, n) = LBFProblem(; H=H, n=n)

@doc raw"""
    LBFProblem(H, n, N)

Constructs a local Berry flux problem with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 2.
- `n`: An `AbstractVector` (or a `Tuple`) including two elements of `Int`, which represents wavenumber ($2\pi n/N$) when calculating Berry flux. Dimension of `n` must be 2.
- `N`: The number of points for one direction in the Brillouin zone.

# Returns
A `LBFProblem` object.

# Example
```julia
julia> 
```
"""
LBFProblem(H, n, N) = LBFProblem(; H=H, n=n, N=N)


# Problem for finding and calculating the Weyl points
@doc raw"""
    WCSProblem{T1<:Function,T2<:String,T3<:Int,T4<:Union{Tuple,AbstractVector,Int},T5<:Real,T6<:Bool} <: TopologicalNumbersProblems

A struct representing a problem for finding and calculating the Weyl points.

# Fields
- `H::T1`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `kn::T2`: Compute the Chern number of the plane perpendicular to the `"kn"` direction in Brillouin zone (`"k1"`, `"k2"`, `"k3"`).
- `kn_mesh::T3`: Number of mesh in `"kn"` direction. Default is 51.
- `N::T4`: The number of points for one direction in the Brillouin zone. Default is 51.
- `gapless::T5`: The threshold for considering a band as gapless. Default is 0.0.
- `rounds::T6`: A boolean indicating whether to round a returned variable. Default is true.

# Example
```julia
julia> 
```
"""
@kwdef struct WCSProblem{T1<:Function,T2<:String,T3<:Int,T4<:Union{Tuple,AbstractVector,Int},T5<:Real,T6<:Bool} <: TopologicalNumbersProblems
    H::T1
    kn::T2
    kn_mesh::T3 = 51
    N::T4 = 51
    gapless::T5 = 0.0
    rounds::T6 = true
end
# default

@doc raw"""
    WCSProblem(H, kn)

Constructs a problem for finding and calculating the Weyl points with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `kn`: Compute the Chern number of the plane perpendicular to the `"kn"` direction in Brillouin zone (`"k1"`, `"k2"`, `"k3"`).

# Returns
A `WCSProblem` object.

# Example
```julia
julia> 
```
"""
WCSProblem(H, kn) = WCSProblem(; H=H, kn=kn)

@doc raw"""
    WCSProblem(H, kn, N::T) where {T<:Int}

Constructs a problem for finding and calculating the Weyl points with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `kn`: Compute the Chern number of the plane perpendicular to the `"kn"` direction in Brillouin zone (`"k1"`, `"k2"`, `"k3"`).
- `N`: The number of points for one direction in the Brillouin zone.

# Returns
A `WCSProblem` object.

# Example
```julia
julia> 
```
"""
WCSProblem(H, kn, N::T) where {T<:Int} = WCSProblem(; H=H, kn=kn, kn_mesh=N, N=N)

@doc raw"""
    WCSProblem(H, kn, N1, N2)

Constructs a problem for finding and calculating the Weyl points with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `kn`: Compute the Chern number of the plane perpendicular to the `"kn"` direction in Brillouin zone (`"k1"`, `"k2"`, `"k3"`).
- `N1`: Number of mesh in `"kn"` direction.
- `N2`: The number of points for one direction in the Brillouin zone.

# Returns
A `WCSProblem` object.

# Example
```julia
julia> 
```
"""
WCSProblem(H, kn, N1, N2) = WCSProblem(; H=H, kn=kn, kn_mesh=N1, N=N2)


# Problem for finding and calculating the Weyl points
@doc raw"""
    WNProblem{T1<:Function,T2<:AbstractVector,T3<:Int,T4<:Real,T5<:Bool} <: TopologicalNumbersProblems

A struct representing a problem for calculating the Weyl nodes.

# Fields
- `H::T1`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `n::T2`: An `AbstractVector` (or a `Tuple`) including two elements of `Int`, which represents wavenumber ($2\pi n/N$). Dimension of `n` must be 3.
- `N::T3`: The number of points for one direction in the Brillouin zone. Default is 51.
- `gapless::T4`: The threshold for considering a band as gapless. Default is 0.0.
- `rounds::T5`: A boolean indicating whether to round a returned variable. Default is true.

# Example
```julia
julia> 
```
"""
@kwdef struct WNProblem{T1<:Function,T2<:AbstractVector,T3<:Int,T4<:Real,T5<:Bool} <: TopologicalNumbersProblems
    H::T1
    n::T2
    N::T3 = 51
    gapless::T4 = 0.0
    rounds::T5 = true
end
# default

@doc raw"""
    WNProblem(H, n)

Constructs a problem for calculating the Weyl nodes with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `n`: An `AbstractVector` (or a `Tuple`) including two elements of `Int`, which represents wavenumber ($2\pi n/N$). Dimension of `n` must be 3.

# Returns
A `WNProblem` object.

# Example
```julia
julia> 
```
"""
WNProblem(H, n) = WNProblem(; H=H, n=n)

@doc raw"""
    WNProblem(H, n, N)

Constructs a problem for calculating the Weyl nodes with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `n`: An `AbstractVector` (or a `Tuple`) including two elements of `Int`, which represents wavenumber ($2\pi n/N$). Dimension of `n` must be 3.
- `N`: The number of points for one direction in the Brillouin zone.

# Returns
A `WNProblem` object.

# Example
```julia
julia> 
```
"""
WNProblem(H, n, N) = WNProblem(; H=H, n=n, N=N)


# Problem for finding and calculating the Weyl points
@doc raw"""
    WPProblem{T1<:Function,T2<:Int,T3<:AbstractVector,T4<:Bool} <: TopologicalNumbersProblems

A struct representing a problem for calculating the Weyl points.

# Fields
- `H::T1`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `N::T2`: The number of meshes when discretizing the Brillouin Zone. The $n$th iteration divides the Brillouin zone into $1/N^n$. Default is 10.
- `gapless::T3` The threshold that determines the state to be degenerate. The $n$th iteration adopts the threshold value of the $n$th value of the vector. The number of iterations can be varied by the length of the vector. Default is `[1e-1, 1e-2, 1e-3, 1e-4]`.
- `rounds::T4`: A boolean indicating whether to round a returned variable. Default is true.

# Example
```julia
julia> 
```
"""
@kwdef struct WPProblem{T1<:Function,T2<:Int,T3<:AbstractVector,T4<:Bool} <: TopologicalNumbersProblems
    H::T1
    N::T2 = 10
    gapless::T3 = [1e-1, 1e-2, 1e-3, 1e-4]
    rounds::T4 = true
end
# default

@doc raw"""
    WPProblem(H)

Constructs a problem for calculating the Weyl points with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.

# Returns
A `WPProblem` object.

# Example
```julia
julia> 
```
"""
WPProblem(H) = WPProblem(; H=H)

@doc raw"""
    WPProblem(H, N)

Constructs a problem for calculating the Weyl points with the default parameters.

# Arguments
- `H`: The Hamiltonian function `H=H(k)` that defines the system. `k` is a abstract vector (or a tuple) of the wavenumber vector. Dimension of `k` must be 3.
- `N`: The number of meshes when discretizing the Brillouin Zone. The $n$th iteration divides the Brillouin zone into $1/N^n$.

# Returns
A `WPProblem` object.

# Example
```julia
julia> 
```
"""
WPProblem(H, N) = WPProblem(; H=H, N=N)



abstract type TopologicalNumbersSolutions end

# Solution for calculating the first Chern number
@doc raw"""
    BPSolution{T1,T2} <: TopologicalNumbersSolutions

The `BPSolution` struct represents a solution for calculating Berry phase.

# Fields
- `TopologicalNumber::T1`: The Berry phase for each energy bands.
- `Total::T2`: The total Berry phase.
"""
@kwdef struct BPSolution{T1,T2} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    Total::T2 = nothing
end

# Solution for calculating the first Chern number
@doc raw"""
    FCSolution{T1,T2} <: TopologicalNumbersSolutions

The `FCSolution` struct represents a solution for calculating the first Chern number.

# Fields
- `TopologicalNumber::T1`: The first Chern number for each energy bands.
- `Total::T2`: The total of the first Chern number.
"""
@kwdef struct FCSolution{T1,T2} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    Total::T2 = nothing
end

# Solution for calculating the second Chern number
@doc raw"""
    SCSolution{T} <: TopologicalNumbersSolutions

The `SCSolution` struct represents a solution for calculating the second Chern number.

# Fields
- `TopologicalNumber::T`: The second Chern number.
"""
@kwdef struct SCSolution{T} <: TopologicalNumbersSolutions
    TopologicalNumber::T = nothing
end

# Solution for calculating the Z2 invariant
@doc raw"""
    Z2Solution{T1,T2,T3} <: TopologicalNumbersSolutions

The `Z2Solution` struct represents a solution for calculating Z2 number.

# Fields
- `TopologicalNumber::T1`: The Z2 number for each pair of energy bands.
- `TRTopologicalNumber::T2`: The Z2 number for the remaining part of the Brillouin zone. This field is only returned when `TR` is `true` in `Z2Problem`.
- `Total::T3`: The total Z2 number.
"""
@kwdef struct Z2Solution{T1,T2,T3} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    TRTopologicalNumber::T2 = nothing
    Total::T3 = nothing
end

# Solution for calculating the local Berry flux
@doc raw"""
    LBFSolution{T1,T2} <: TopologicalNumbersSolutions

The `LBFSolution` struct represents a solution for calculating the $k$-local value of Berry flux.

# Fields
- `TopologicalNumber::T1`: The local Berry flux for each energy bands.
- `n::T2`: An `AbstractVector` (or a `Tuple`) including two elements of `Int`, which represents wavenumber ($2\pi n/N$) when calculating Berry flux.
"""
@kwdef struct LBFSolution{T1,T2} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    n::T2 = nothing
end

# Solution for finding and calculating the Weyl points
@doc raw"""
    WCSSolution{T1,T2,T3} <: TopologicalNumbersSolutions

The `WCSSolution` struct represents a solution for finding and calculating the Weyl points.

# Fields
- `kn::T1`: Compute the Chern number of the plane perpendicular to the `"kn"` direction in Brillouin zone (`"k1"`, `"k2"`, `"k3"`).
- `param::T2`: Wavenumber parameters in `"kn"` direction.
- `nums::T3`: The Chern number for each energy bands and each wavenumber parameters.
"""
@kwdef struct WCSSolution{T1,T2,T3} <: TopologicalNumbersSolutions
    kn::T1 = nothing
    param::T2 = nothing
    nums::T3 = nothing
end

# Solution for finding and calculating the Weyl points
@doc raw"""
    WNSolution{T1,T2,T3} <: TopologicalNumbersSolutions

The `WNSolution` struct represents a solution for calculating the Weyl nodes.

# Fields
- `TopologicalNumber::T1`: 
- `n::T2`: 
- `N::T3`: 
"""
@kwdef struct WNSolution{T1,T2,T3} <: TopologicalNumbersSolutions
    TopologicalNumber::T1 = nothing
    n::T2 = nothing
    N::T3 = nothing
end

# Solution for finding and calculating the Weyl points
@doc raw"""
    WPSolution{T1,T2,T3} <: TopologicalNumbersSolutions

The `WPSolution` struct represents a solution for calculating the Weyl points.

# Fields
- `WeylPoint::T1`: 
- `N::T2`: 
- `Nodes::T3`:
"""
@kwdef struct WPSolution{T1,T2,T3} <: TopologicalNumbersSolutions
    WeylPoint::T1 = nothing
    N::T2 = nothing
    Nodes::T3 = nothing
end
