@with_kw mutable struct Params{T1<:Int,T2<:AbstractFloat,T3<:Bool}
    Hamiltonian::Function
    dim::T1
    N::T1
    Hs::T1
    gapless::T2
    rounds::T3
end

# Define a mutable struct named TemporalSecondChern with multiple type parameters
@with_kw mutable struct TemporalSecondChern{T1,T2,T3,T4,T5,T6}
    chern::T1      # Second Chern number
    k::T2          # Momentum vector
    evec::T3       # Eigenvectors
    H::T4          # Hamiltonian
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