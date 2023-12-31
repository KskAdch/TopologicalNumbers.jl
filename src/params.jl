
@with_kw mutable struct Params{T1<:Int,T2<:Union{Int,Tuple,AbstractVector},T3<:AbstractFloat,T4<:Bool}
    Hamiltonian::Function
    dim::T1
    Nfill::T1 = 0
    Hs::T1
    N::T2
    gapless::T3
    rounds::T4
end

# Define a mutable struct named TemporalSecondChern with multiple type parameters
@with_kw mutable struct TemporalSecondChern{T1,T2,T3,T4,T5,T6}
    chern::T1      # Second Chern number
    k::T2          # Momentum vector
    evec::T3       # Eigenvectors
    H::T4          # Hamiltonian or eigen vectors
    Ux::T5         # Link variables
    Uy::T5
    Uz::T5
    Uw::T5
    Ux_inv::T5     # Inverse Link variables
    Uy_inv::T5
    Uz_inv::T5
    Uw_inv::T5
    Fxy::T5        # Field strength tensor components
    Fzw::T5
    Fwx::T5
    Fzy::T5
    Fzx::T5
    Fyw::T5
    Ftemp::T6      # Temporary variable for LinearAlgebra.log calculations
end