function calcChernSurface(Hamiltonian::Function, kn::String, kn_range::T; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}
    
    if kn == "k1"
        k0(k, p) = [p, k[1], k[2]]
    elseif kn == "k2"
        k0(k, p) = [k[1], p, k[2]]
    elseif kn == "k3"
        k0(k, p) = [k[1], k[2], p]
    else
        throw(ArgumentError("Unknown keyword $kn"))
    end
    
    H0(k, p) = Hamiltonian(k0(k, p))

    Hs = size(Hamiltonian(zeros(3)))[1]
    param = Params(; Hamiltonian, dim=2, N, gapless, rounds, Hs)
    # capcPhaseDigram(Hamiltonian, kn_range, "Chern")
    calc_data1D!(nums, param)
    
    if plot == true
        plot1D(nums, param_range)
    end

    (; param=param_range, nums)
end