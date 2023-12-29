function update1D!(::T, nums, num0, H, alg!, range1, p::Params) where {T<:SecondChernAlgorithms}
    for i in eachindex(range1)
        Ham0(k) = H(k, range1[i])
        p.Hamiltonian = Ham0

        v = setParams(p) # Set parameters
        setBasis!(v, p) # Set the basis

        alg!(v)
        nums[i] = v.sys.chern
    end
end

function update1D!(::T, nums, num0, H, alg!, range1, idxs::ProgressBar, p::Params) where {T<:SecondChernAlgorithms}
    for i in idxs
        Ham0(k) = H(k, range1[i])
        p.Hamiltonian = Ham0

        v = setParams(p) # Set parameters
        setBasis!(v, p) # Set the basis

        alg!(v)
        nums[i] = v.sys.chern
    end
end

function update2D!(::T, nums, num0, H, alg!, range1, range2, p::Params) where {T<:SecondChernAlgorithms}
    for i in eachindex(range1), j in eachindex(range2)
        param = (range1[i], range2[j])
        Ham0(k) = H(k, param)
        p.Hamiltonian = Ham0

        v = setParams(p) # Set parameters
        setBasis!(v, p) # Set the basis

        alg!(v)
        nums[i, j] = v.sys.chern
    end
end

function update2D!(::T, nums, num0, H, alg!, range1, range2, idxs::ProgressBar, p::Params) where {T<:SecondChernAlgorithms}
    for i in idxs, j in eachindex(range2)
        param = (range1[i], range2[j])
        Ham0(k) = H(k, param)
        p.Hamiltonian = Ham0

        v = setParams(p) # Set parameters
        setBasis!(v, p) # Set the basis

        alg!(v)
        nums[i, j] = v.sys.chern
    end
end

function update1D!(nums, num0, H, alg!, range1, p::Params)
    for i in eachindex(range1)
        Ham0(k) = H(k, range1[i])
        p.Hamiltonian = Ham0
        alg!(num0, p)
        nums[:, i] .= num0
    end
end

function update1D!(nums, num0, H, alg!, range1, idxs::ProgressBar, p::Params)
    for i in idxs
        Ham0(k) = H(k, range1[i])
        p.Hamiltonian = Ham0
        alg!(num0, p)
        nums[:, i] .= num0
    end
end

function update2D!(nums, num0, H, alg!, range1, range2, p::Params)
    for i in eachindex(range1), j in eachindex(range2)
        param = (range1[i], range2[j])
        Ham0(k) = H(k, param)
        p.Hamiltonian = Ham0
        alg!(num0, p)
        nums[:, i, j] .= num0
    end
end

function update2D!(nums, num0, H, alg!, range1, range2, idxs::ProgressBar, p::Params)
    for i in idxs, j in eachindex(range2)
        param = (range1[i], range2[j])
        Ham0(k) = H(k, param)
        p.Hamiltonian = Ham0
        alg!(num0, p)
        nums[:, i, j] .= num0
    end
end

function calc_data1D(H, param_range, alg::T, p::Params) where {T<:TopologicalNumbersAlgorithms}
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range))
    num0 = zeros(Float64, Hs)

    if alg == SecondChern_FHS()
        algorithm! = SecondChernPhase!
        type0 = ifelse(rounds, Float64, ComplexF64)
        nums = zeros(ComplexF64, length(param_range))
        num0 = zero(type0)
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update1D!(alg, nums, num0, H, algorithm!, param_range, p)

    if alg == SecondChern_FHS()
        if rounds == true
            real.(nums)
        elseif rounds == false
            nums
        end
    else
        if rounds == true
            round.(Int, transpose(nums))
        elseif rounds == false
            transpose(nums)
        end
    end
end

function calc_data1D(H, param_range, idxs::ProgressBar, alg::T, p::Params) where {T<:TopologicalNumbersAlgorithms}
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range))
    num0 = zeros(Float64, Hs)

    if alg == SecondChern_FHS()
        algorithm! = SecondChernPhase!
        type0 = ifelse(rounds, Float64, ComplexF64)
        nums = zeros(ComplexF64, length(param_range))
        num0 = zero(type0)
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update1D!(alg, nums, num0, H, algorithm!, param_range, idxs, p)

    if alg == SecondChern_FHS()
        if rounds == true
            real.(nums)
        elseif rounds == false
            nums
        end
    else
        if rounds == true
            round.(Int, transpose(nums))
        elseif rounds == false
            transpose(nums)
        end
    end
end

function calc_data2D(H, param_range1, param_range2, alg::T, p::Params) where {T<:TopologicalNumbersAlgorithms}
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range1), length(param_range2))
    num0 = zeros(Float64, Hs)

    if alg == SecondChern_FHS()
        algorithm! = SecondChernPhase!
        type0 = ifelse(rounds, Float64, ComplexF64)
        nums = zeros(ComplexF64, length(param_range1), length(param_range2))
        num0 = zero(type0)
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update2D!(alg, nums, num0, H, algorithm!, param_range1, param_range2, p)

    if alg == SecondChern_FHS()
        if rounds == true
            real.(nums)
        elseif rounds == false
            nums
        end
    else
        if rounds == true
            round.(Int, nums)
        elseif rounds == false
            nums
        end
    end
end

function calc_data2D(H, param_range1, param_range2, idxs::ProgressBar, alg::T, p::Params) where {T<:TopologicalNumbersAlgorithms}
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range1), length(param_range2))
    num0 = zeros(Float64, Hs)

    if alg == SecondChern_FHS()
        algorithm! = SecondChernPhase!
        type0 = ifelse(rounds, Float64, ComplexF64)
        nums = zeros(ComplexF64, length(param_range1), length(param_range2))
        num0 = zero(type0)
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update2D!(alg, nums, num0, H, algorithm!, param_range1, param_range2, idxs, p)

    if alg == SecondChern_FHS()
        if rounds == true
            real.(nums)
        elseif rounds == false
            nums
        end
    else
        if rounds == true
            round.(Int, nums)
        elseif rounds == false
            nums
        end
    end
end

# Old method
function calc_data1D(H, param_range, alg::String, p::Params)
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range))
    num0 = zeros(Float64, Hs)

    if alg == "BerryPhase"
        algorithm! = BerryPhase!
    elseif alg == "Z2"
        algorithm! = Z2Phase!
        nums = zeros(Float64, Hs ÷ 2, length(param_range))
        num0 = zeros(Float64, Hs ÷ 2)
    elseif alg == "Chern"
        algorithm! = ChernPhase!
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update1D!(nums, num0, H, algorithm!, param_range, p)

    if rounds == true
        round.(Int, transpose(nums))
    elseif rounds == false
        transpose(nums)
    end
end

# Old method
function calc_data1D(H, param_range, idxs::ProgressBar, alg::String, p::Params)
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range))
    num0 = zeros(Float64, Hs)

    if alg == "BerryPhase"
        algorithm! = BerryPhase!
    elseif alg == "Z2"
        algorithm! = Z2Phase!
        nums = zeros(Float64, Hs ÷ 2, length(param_range))
        num0 = zeros(Float64, Hs ÷ 2)
    elseif alg == "Chern"
        algorithm! = ChernPhase!
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update1D!(nums, num0, H, algorithm!, param_range, idxs, p)

    if rounds == true
        round.(Int, transpose(nums))
    elseif rounds == false
        transpose(nums)
    end
end

# Old method
function calc_data2D(H, param_range1, param_range2, alg::String, p::Params)
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range1), length(param_range2))
    num0 = zeros(Float64, Hs)

    if alg == "BerryPhase"
        algorithm! = BerryPhase!
    elseif alg == "Z2"
        algorithm! = Z2Phase!
        nums = zeros(Float64, Hs ÷ 2, length(param_range1), length(param_range2))
        num0 = zeros(Float64, Hs ÷ 2)
    elseif alg == "Chern"
        algorithm! = ChernPhase!
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update2D!(nums, num0, H, algorithm!, param_range1, param_range2, p)

    if rounds == true
        round.(Int, nums)
    elseif rounds == false
        nums
    end
end

# Old method
function calc_data2D(H, param_range1, param_range2, idxs::ProgressBar, alg::String, p::Params)
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range1), length(param_range2))
    num0 = zeros(Float64, Hs)

    if alg == "BerryPhase"
        algorithm! = BerryPhase!
    elseif alg == "Z2"
        algorithm! = Z2Phase!
        nums = zeros(Float64, Hs ÷ 2, length(param_range1), length(param_range2))
        num0 = zeros(Float64, Hs ÷ 2)
    elseif alg == "Chern"
        algorithm! = ChernPhase!
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update2D!(nums, num0, H, algorithm!, param_range1, param_range2, idxs, p)

    if rounds == true
        round.(Int, nums)
    elseif rounds == false
        nums
    end
end

@doc raw"""
    calcPhaseDiagram(H::Function, param_range::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}
"""
function calcPhaseDiagram(H::Function, param_range::T1, alg::T2; N::T3=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:Union{String,TopologicalNumbersAlgorithms},T3<:Union{Int,Tuple,AbstractVector}}

    dim = Hs = 0
    Hamiltonian(k) = H(k, 0.0)

    try
        Hs = size(Hamiltonian(0.0), 1)
        dim = 1
    catch
        try
            Hs = size(Hamiltonian(zeros(2)), 1)
            dim = 2
        catch
            try
                Hs = size(Hamiltonian(zeros(3)), 1)
                dim = 3
            catch
                Hs = size(Hamiltonian(zeros(4)), 1)
                dim = 4
            end
        end
    end

    p = Params(; Hamiltonian, dim, N, gapless, rounds, Hs)
    if alg == SecondChern_FHS()
        @reset p.N = (N, N, N, N)
        @reset p.Nfill = Hs ÷ 2
    end

    nums = if progress == true
        calc_data1D(H, param_range, ProgressBar(eachindex(param_range)), alg, p)
    else
        calc_data1D(H, param_range, alg, p)
    end

    if plot == true
        if alg == SecondChern_FHS()
            plot1D(alg, nums, param_range)
        else
            plot1D(nums, param_range)
        end
    end

    (; param=param_range, nums)
end



@doc raw"""
    calcPhaseDiagram(H::Function, param_range1::T, param_range2::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}
"""
function calcPhaseDiagram(H::Function, param_range1::T1, param_range2::T2, alg::T3; N::T4=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:Union{String,TopologicalNumbersAlgorithms},T4<:Union{Int,Tuple,AbstractVector}}

    dim = Hs = 0
    Hamiltonian(k) = H(k, (0.0, 0.0))

    try
        Hs = size(Hamiltonian(0.0), 1)
        dim = 1
    catch
        try
            Hs = size(Hamiltonian(zeros(2)), 1)
            dim = 2
        catch
            try
                Hs = size(Hamiltonian(zeros(3)), 1)
                dim = 3
            catch
                Hs = size(Hamiltonian(zeros(4)), 1)
                dim = 4
            end
        end
    end

    p = Params(; Hamiltonian, dim, N, gapless, rounds, Hs)
    if alg == SecondChern_FHS()
        @reset p.N = (N, N, N, N)
        @reset p.Nfill = Hs ÷ 2
    end

    nums = if progress == true
        calc_data2D(H, param_range1, param_range2, ProgressBar(eachindex(param_range1)), alg, p)
    else
        calc_data2D(H, param_range1, param_range2, alg, p)
    end

    if alg == SecondChern_FHS()
        if plot == true
            plot2D(transpose(nums_half), param_range1, param_range2)
        end
    else
        if plot == true && Hs % 2 == 0
            nums_half = sum(@view(nums[1:end÷2, :, :]), dims=1)[1, :, :]
            plot2D(transpose(nums_half), param_range1, param_range2) # half-filling case
        end
    end

    (; param1=param_range1, param2=param_range2, nums)
end
