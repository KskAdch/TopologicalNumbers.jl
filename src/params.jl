
@kwdef struct Params{T1<:Int,T2<:Union{Int,Tuple,AbstractVector},T3<:AbstractFloat,T4<:Bool}
    Ham::Function
    dim::T1
    Nfill::T1 = 1
    Hs::T1
    N::T2
    gapless::T3 = 0.0
    rounds::T4 = true
    returnRealValue::T4 = true
end

# Define a mutable struct named TemporalSecondChern with multiple type parameters
@kwdef mutable struct TemporalSecondChern{T1,T2,T3,T4,T5,T6}
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


@kwdef struct TemporalBerryPhase{L,EV,PS}
    Link::L
    Evec0::EV
    Evec1::EV
    psi0::PS
    psi1::PS
    psiN1::PS
end


@kwdef struct TemporalFirstChern{K,L1,L2,EV1,EV2,PS1,PS2}
    k::K
    Link0::L1
    Link1::L1
    LinkN::L1
    link10::L2
    link01::L2
    psi_0::PS1
    psi_1::PS1
    psi_N::PS1
    Evec0::EV1
    Evec1::EV1
    psi00::PS2
    psi10::PS2
    psi01::PS2
    Enevec::EV2
end


@kwdef struct TemporalZ2{K,B,W,L1,L2,PS1,PS2,N}
    k::K
    T::B
    w00::W
    w0p::W
    wp0::W
    wpp::W
    Link1::L1
    Link2::L1
    LinkN1::L1
    link10::L2
    link01::L2
    psi_0::PS1
    psi_1::PS1
    psi_N::PS1
    psi00::PS2
    psi10::PS2
    psi01::PS2
    Nhalf::N
    Hshalf::N
end