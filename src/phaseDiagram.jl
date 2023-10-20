function update1D!(nums, num0, H, alg!, range1::T, p::Params) where {T<:AbstractVector}
    for i in eachindex(range1)
        Ham0(k) = H(k, range1[i])
        p.Hamiltonian = Ham0
        alg!(num0, p)
        nums[:, i] .= num0
    end
end

function update2D!(nums, num0, H, alg!, range1::T1, range2::T2, p::Params) where {T1<:AbstractVector,T2<:AbstractVector}
    for i in eachindex(range1), j in eachindex(range2)
        param = (range1[i], range2[j])
        Ham0(k) = H(k, param)
        p.Hamiltonian = Ham0
        alg!(num0, p)
        nums[:, i, j] .= num0
    end
end

function calc_data1D(H, param_range, alg, p::Params)
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, size(param_range, 1))
    num0 = zeros(Float64, Hs)

    # if rounds == true
    #     num0 = zeros(Int64, Hs)

    #     if alg == "BerryPhase"
    #         algorithm! = BerryPhase_round!
    #     elseif alg == "Z2"
    #         algorithm! = Z2Phase_round!
    #         nums = zeros(Float64, Hs ÷ 2, size(param_range, 1))
    #         num0 = zeros(Int64, Hs ÷ 2)
    #     elseif alg == "Chern"
    #         algorithm! = ChernPhase!
    #     else
    #         throw(ArgumentError("Unknown algorithm $alg"))
    #     end

    #     update1D!(nums, num0, H, algorithm!, param_range, p)
    #     # nums = Int.(transpose(nums))
    #     return Int.(transpose(nums))
    # elseif rounds == false
    #     num0 = zeros(Float64, Hs)

    #     if alg == "BerryPhase"
    #         algorithm! = BerryPhase!
    #     elseif alg == "Z2"
    #         algorithm! = Z2Phase!
    #         nums = zeros(Float64, Hs ÷ 2, size(param_range, 1))
    #         num0 = zeros(Float64, Hs ÷ 2)
    #     elseif alg == "Chern"
    #         algorithm! = ChernPhase!
    #     else
    #         throw(ArgumentError("Unknown algorithm $alg"))
    #     end

    #     update1D!(nums, num0, H, algorithm!, param_range, p)
    #     # nums = transpose(nums)
    #     return transpose(nums)
    # end
    
    if rounds == true

        if alg == "BerryPhase"
            algorithm! = BerryPhase!
        elseif alg == "Z2"
            algorithm! = Z2Phase!
            nums = zeros(Float64, Hs ÷ 2, size(param_range, 1))
            num0 = zeros(Float64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        else
            throw(ArgumentError("Unknown algorithm $alg"))
        end

        update1D!(nums, num0, H, algorithm!, param_range, p)
        nums = round.(Int, transpose(nums))
    elseif rounds == false

        if alg == "BerryPhase"
            algorithm! = BerryPhase!
        elseif alg == "Z2"
            algorithm! = Z2Phase!
            nums = zeros(Float64, Hs ÷ 2, size(param_range, 1))
            num0 = zeros(Float64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        else
            throw(ArgumentError("Unknown algorithm $alg"))
        end

        update1D!(nums, num0, H, algorithm!, param_range, p)
        nums = transpose(nums)
    end
end

function calc_data2D(H, param_range1, param_range2, alg, p::Params)
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, size(param_range1, 1), size(param_range2, 1))
    num0 = zeros(Float64, Hs)

    # if rounds == true
    #     num0 = zeros(Int64, Hs)

    #     if alg == "BerryPhase"
    #         algorithm! = BerryPhase_round!
    #     elseif alg == "Z2"
    #         algorithm! = Z2Phase_round!
    #         nums = zeros(Float64, Hs ÷ 2, size(param_range1, 1), size(param_range2, 1))
    #         num0 = zeros(Int64, Hs ÷ 2)
    #     elseif alg == "Chern"
    #         algorithm! = ChernPhase!
    #     else
    #         throw(ArgumentError("Unknown algorithm $alg"))
    #     end

    #     update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)
    #     # nums = Int.(nums)
    #     return Int.(nums)
    # elseif rounds == false
    #     num0 = zeros(Float64, Hs)

    #     if alg == "BerryPhase"
    #         algorithm! = BerryPhase!
    #     elseif alg == "Z2"
    #         algorithm! = Z2Phase!
    #         nums = zeros(Float64, Hs ÷ 2, size(param_range1, 1), size(param_range2, 1))
    #         num0 = zeros(Float64, Hs ÷ 2)
    #     elseif alg == "Chern"
    #         algorithm! = ChernPhase!
    #     else
    #         throw(ArgumentError("Unknown algorithm $alg"))
    #     end

    #     update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)
    #     return nums
    # end
    
    if rounds == true

        if alg == "BerryPhase"
            algorithm! = BerryPhase!
        elseif alg == "Z2"
            algorithm! = Z2Phase!
            nums = zeros(Float64, Hs ÷ 2, size(param_range1, 1), size(param_range2, 1))
            num0 = zeros(Float64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        else
            throw(ArgumentError("Unknown algorithm $alg"))
        end

        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)
        nums = round.(Int, nums)
    elseif rounds == false

        if alg == "BerryPhase"
            algorithm! = BerryPhase!
        elseif alg == "Z2"
            algorithm! = Z2Phase!
            nums = zeros(Float64, Hs ÷ 2, size(param_range1, 1), size(param_range2, 1))
            num0 = zeros(Float64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        else
            throw(ArgumentError("Unknown algorithm $alg"))
        end

        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)
    end
end

@doc raw"""
    calcPhaseDiagram(H::Function, param_range::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}
"""
function calcPhaseDiagram(H::Function, param_range::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}

    dim = Hs = 0
    Hamiltonian(k) = H(k, 0.0)

    try
        Hs = size(Hamiltonian(0.0), 1)
        dim = 1
    catch
        Hs = size(Hamiltonian(zeros(2)), 1)
        dim = 2
    end

    p = Params(; Hamiltonian, dim, N, gapless, rounds, Hs)

    # nums = zeros(Float64, Hs, size(param_range, 1))
    # num0 = zeros(Float64, Hs)

    nums = calc_data1D(H, param_range, alg, p)

    # if rounds == true

    #     if alg == "BerryPhase"
    #         algorithm! = BerryPhase!
    #     elseif alg == "Z2"
    #         algorithm! = Z2Phase!
    #         nums = zeros(Float64, Hs ÷ 2, size(param_range, 1))
    #         num0 = zeros(Float64, Hs ÷ 2)
    #     elseif alg == "Chern"
    #         algorithm! = ChernPhase!
    #     else
    #         throw(ArgumentError("Unknown algorithm $alg"))
    #     end

    #     update1D!(nums, num0, H, algorithm!, param_range, p)
    #     nums = round.(Int, transpose(nums))
    # elseif rounds == false

    #     if alg == "BerryPhase"
    #         algorithm! = BerryPhase!
    #     elseif alg == "Z2"
    #         algorithm! = Z2Phase!
    #         nums = zeros(Float64, Hs ÷ 2, size(param_range, 1))
    #         num0 = zeros(Float64, Hs ÷ 2)
    #     elseif alg == "Chern"
    #         algorithm! = ChernPhase!
    #     else
    #         throw(ArgumentError("Unknown algorithm $alg"))
    #     end

    #     update1D!(nums, num0, H, algorithm!, param_range, p)
    #     nums = transpose(nums)
    # end

    if plot == true
        plot1D(nums, param_range)
    end

    (; param=param_range, nums)
end



@doc raw"""
    calcPhaseDiagram(H::Function, param_range1::T, param_range2::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}
"""
function calcPhaseDiagram(H::Function, param_range1::T1, param_range2::T2, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector}

    dim = Hs = 0
    Hamiltonian(k) = H(k, (0.0, 0.0))

    try
        Hs = size(Hamiltonian(0.0), 1)
        dim = 1
    catch
        Hs = size(Hamiltonian(zeros(2)), 1)
        dim = 2
    end

    p = Params(; Hamiltonian, dim, N, gapless, rounds, Hs)

    # nums = zeros(Float64, Hs, size(param_range1, 1), size(param_range2, 1))
    # num0 = zeros(Float64, Hs)

    nums = calc_data2D(H, param_range1, param_range2, alg, p)

    if rounds == true

        if alg == "BerryPhase"
            algorithm! = BerryPhase!
        elseif alg == "Z2"
            algorithm! = Z2Phase!
            nums = zeros(Float64, Hs ÷ 2, size(param_range1, 1), size(param_range2, 1))
            num0 = zeros(Float64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        else
            throw(ArgumentError("Unknown algorithm $alg"))
        end

        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)
        nums = round.(Int, nums)
    elseif rounds == false

        if alg == "BerryPhase"
            algorithm! = BerryPhase!
        elseif alg == "Z2"
            algorithm! = Z2Phase!
            nums = zeros(Float64, Hs ÷ 2, size(param_range1, 1), size(param_range2, 1))
            num0 = zeros(Float64, Hs ÷ 2)
        elseif alg == "Chern"
            algorithm! = ChernPhase!
        else
            throw(ArgumentError("Unknown algorithm $alg"))
        end

        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)
    end

    if plot == true && Hs % 2 == 0
        nums_half = sum(@view(nums[1:end÷2, :, :]), dims=1)[1, :, :]
        plot2D(nums_half, param_range1, param_range2) # half-filling case
    end

    (; param1=param_range1, param2=param_range2, nums)
end