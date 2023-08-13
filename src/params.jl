@with_kw mutable struct Params{T1<:Int,T2<:AbstractFloat,T3<:Bool}
    Hamiltonian::Function
    dim::T1
    N::T1
    Hs::T1
    gapless::T2
    rounds::T3
end