function calcChernSurface(H::Function, kn::String, kn_range::T; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}
    
    Hs = size(H(zeros(3)))[1]

    if kn == "k1"
        k0(k, p) = [p, k[1], k[2]]
    elseif kn == "k2"
        k0(k, p) = [k[1], p, k[2]]
    elseif kn == "k3"
        k0(k, p) = [k[1], k[2], p]
    else
        throw(ArgumentError("Unknown keyword $kn"))
    end
    
    H0(k, p) = H(k0(k, p))
    Hamiltonian(k) = H0(k, 0.0)

    param = Params(; Hamiltonian, dim=2, N, gapless, rounds, Hs)

    nums = calc_data1D(H0, param_range, alg, param)
    
    if plot == true
        plot1D(nums, param_range)
    end

    (; kn=kn_range, nums)
end