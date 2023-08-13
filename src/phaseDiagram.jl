function update1D!(nums, num0, H, alg!, range1::T, p::Params) where {T<:AbstractVector}
    for i in eachindex(range1)
        Ham0(k) = H(k, range1[i])
        p.Hamiltonian = Ham0
        alg!(num0, p)
        nums[:, i] .= num0
    end
end

function update2D!(nums, num0, H, alg!, range1::T, range2::T, p::Params) where {T<:AbstractVector}
    for j in eachindex(range1), i in eachindex(range2)
        param = (range1[i], range2[j])
        Ham0(k) = H(k, param)
        p.Hamiltonian = Ham0
        alg!(num0, p)
        nums[:, i, j] .= num0
    end
end

@doc raw"""
    calcPhaseDiagram(H::Function, param_range::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}
"""
function calcPhaseDiagram(H::Function, param_range::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}

    dim = Hs = 0
    Hamiltonian(k) = H(k, 0.0)

    try
        Hs = size(Hamiltonian(0.0))[1]
        dim = 1
    catch
        Hs = size(Hamiltonian(zeros(2)))[1]
        dim = 2
    end

    p = Params(; Hamiltonian, dim, N, gapless, rounds, Hs)

    nums = zeros(Float64, Hs, size(param_range, 1))

    if rounds == true
        num0 = zeros(Int64, Hs)

        if alg == "BerryPhase"
            algorithm! = BerryPhase_round!
        elseif alg == "Z2"
            algorithm! = Z2Phase_round!
            nums = zeros(Float64, Hs ÷ 2, size(param_range, 1))
            num0 = zeros(Int64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        end

        update1D!(nums, num0, H, algorithm!, param_range, p)
        nums = Int.(transpose(nums))
    elseif rounds == false
        num0 = zeros(Float64, Hs)

        if alg == "BerryPhase"
            algorithm! = BerryPhase!
        elseif alg == "Z2"
            algorithm! = Z2Phase!
            nums = zeros(Float64, Hs ÷ 2, size(param_range, 1))
            num0 = zeros(Float64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        end

        update1D!(nums, num0, H, algorithm!, param_range, p)
        nums = transpose(nums)
    end

    if plot == true
        plot1D(nums, param_range)
    end

    nums
end



@doc raw"""
    calcPhaseDiagram(H::Function, param_range1::T, param_range2::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}
"""
function calcPhaseDiagram(H::Function, param_range1::T, param_range2::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}

    dim = Hs = 0
    Hamiltonian(k) = H(k, (0.0, 0.0))

    try
        Hs = size(Hamiltonian(0.0))[1]
        dim = 1
    catch
        Hs = size(Hamiltonian(zeros(2)))[1]
        dim = 2
    end

    p = Params(; Hamiltonian, dim, N, gapless, rounds, Hs)

    nums = zeros(Float64, Hs, size(param_range1, 1), size(param_range2, 1))

    if rounds == true
        num0 = zeros(Int64, Hs)

        if alg == "BerryPhase"
            algorithm! = BerryPhase_round!
        elseif alg == "Z2"
            algorithm! = Z2Phase_round!
            nums = zeros(Float64, Hs ÷ 2, size(param_range1, 1), size(param_range2, 1))
            num0 = zeros(Int64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        end

        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)
        nums = Int.(nums)
    elseif rounds == false
        num0 = zeros(Float64, Hs)

        if alg == "BerryPhase"
            algorithm! = BerryPhase!
        elseif alg == "Z2"
            algorithm! = Z2Phase!
            nums = zeros(Float64, Hs ÷ 2, size(param_range1, 1), size(param_range2, 1))
            num0 = zeros(Float64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        end

        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)
    end

    if plot == true && Hs % 2 == 0
        nums_half = sum(@view(nums[1:Hs÷2, :, :]), dims=1)[1, :, :]
        plot2D(nums_half, param_range1, param_range2) # half-filling case
    end

    nums
end