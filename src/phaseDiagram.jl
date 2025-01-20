function update1Din!(
    ::T, i, nums, num0, H, alg!, range1, p::Params
) where {T<:SecondChernAlgorithms}
    Ham0(k) = H(k, range1[i])
    @reset p.Ham = Ham0

    v = setParams(p) # Set parameters
    setBasis!(v, p) # Set the basis

    alg!(v)
    return nums[i] = v.sys.chern
end

function update1Din!(
    ::T, i, nums, num0, H, alg!, range1, p::Params
) where {T<:Union{BerryPhaseAlgorithms,FirstChernAlgorithms}}
    Ham0(k) = H(k, range1[i])
    @reset p.Ham = Ham0
    alg!(num0, p)
    return nums[:, i] .= num0
end

function update1Din!(
    ::T, i, nums, num0, H, alg!, range1, p::Params
) where {T<:Union{Z2Algorithms}}
    Ham0(k) = H(k, range1[i])
    @reset p.Ham = Ham0

    v = setTemporalZ2(p)
    alg!(v, p)
    return nums[:, i] .= v.num
end

# Old method
function update1Din!(i, nums, num0, H, alg!, range1, p::Params)
    Ham0(k) = H(k, range1[i])
    @reset p.Ham = Ham0
    alg!(num0, p)
    return nums[:, i] .= num0
end

function update2Din!(
    ::T, i, j, nums, num0, H, alg!, range1, range2, p::Params
) where {T<:SecondChernAlgorithms}
    param = (range1[i], range2[j])
    Ham0(k) = H(k, param)
    @reset p.Ham = Ham0

    v = setParams(p) # Set parameters
    setBasis!(v, p) # Set the basis

    alg!(v)
    return nums[i, j] = v.sys.chern
end

function update2Din!(
    ::T, i, j, nums, num0, H, alg!, range1, range2, p::Params
) where {T<:Union{BerryPhaseAlgorithms,FirstChernAlgorithms}}
    param = (range1[i], range2[j])
    Ham0(k) = H(k, param)
    @reset p.Ham = Ham0
    alg!(num0, p)
    return nums[:, i, j] .= num0
end

function update2Din!(
    ::T, i, j, nums, num0, H, alg!, range1, range2, p::Params
) where {T<:Union{Z2Algorithms}}
    param = (range1[i], range2[j])
    Ham0(k) = H(k, param)
    @reset p.Ham = Ham0

    v = setTemporalZ2(p)
    alg!(v, p)
    return nums[:, i, j] .= v.num
end

# Old method
function update2Din!(i, j, nums, num0, H, alg!, range1, range2, p::Params)
    param = (range1[i], range2[j])
    Ham0(k) = H(k, param)
    @reset p.Ham = Ham0
    alg!(num0, p)
    return nums[:, i, j] .= num0
end

function update1D!(
    ::T, nums, num0, H, alg!, range1, ::UseSingleThread, p::Params
) where {T<:TopologicalNumbersAlgorithms}
    for i in eachindex(range1)
        update1Din!(T(), i, nums, num0, H, alg!, range1, p)
    end
end

function update1D!(
    ::T, nums, num0, H, alg!, range1, idxs::ProgressBar, ::UseSingleThread, p::Params
) where {T<:TopologicalNumbersAlgorithms}
    for i in idxs
        update1Din!(T(), i, nums, num0, H, alg!, range1, p)
    end
end

function update2D!(
    ::T, nums, num0, H, alg!, range1, range2, ::UseSingleThread, p::Params
) where {T<:TopologicalNumbersAlgorithms}
    for i in eachindex(range1), j in eachindex(range2)
        update2Din!(T(), i, j, nums, num0, H, alg!, range1, range2, p)
    end
end

function update2D!(
    ::T,
    nums,
    num0,
    H,
    alg!,
    range1,
    range2,
    idxs::ProgressBar,
    ::UseSingleThread,
    p::Params,
) where {T<:TopologicalNumbersAlgorithms}
    for i in idxs, j in eachindex(range2)
        update2Din!(T(), i, j, nums, num0, H, alg!, range1, range2, p)
    end
end

# Old method
function update1D!(nums, num0, H, alg!, range1, ::UseSingleThread, p::Params)
    for i in eachindex(range1)
        update1Din!(i, nums, num0, H, alg!, range1, p)
    end
end

# Old method
function update1D!(
    nums, num0, H, alg!, range1, idxs::ProgressBar, ::UseSingleThread, p::Params
)
    for i in idxs
        update1Din!(i, nums, num0, H, alg!, range1, p)
    end
end

# Old method
function update2D!(nums, num0, H, alg!, range1, range2, ::UseSingleThread, p::Params)
    for i in eachindex(range1), j in eachindex(range2)
        update2Din!(i, j, nums, num0, H, alg!, range1, range2, p)
    end
end

# Old method
function update2D!(
    nums, num0, H, alg!, range1, range2, idxs::ProgressBar, ::UseSingleThread, p::Params
)
    for i in idxs, j in eachindex(range2)
        update2Din!(i, j, nums, num0, H, alg!, range1, range2, p)
    end
end

# function update1D!(::T, nums, num0, H, alg!, range1, ::UseThreads, p::Params) where {T<:SecondChernAlgorithms}
#     Base.Threads.@threads for i in eachindex(range1)
#         update1Din!(T(), i, nums, num0, H, alg!, range1, p)
#     end
# end

# function update1D!(::T, nums, num0, H, alg!, range1, idxs::ProgressBar, ::UseThreads, p::Params) where {T<:SecondChernAlgorithms}
#     Base.Threads.@threads for i in idxs
#         update1Din!(T(), i, nums, num0, H, alg!, range1, p)
#     end
# end

# function update2D!(::T, nums, num0, H, alg!, range1, range2, ::UseThreads, p::Params) where {T<:SecondChernAlgorithms}
#     Base.Threads.@threads for i in eachindex(range1)
#         for j in eachindex(range2)
#             update2Din!(T(), i, j, nums, num0, H, alg!, range1, range2, p)
#         end
#     end
# end

# function update2D!(::T, nums, num0, H, alg!, range1, range2, idxs::ProgressBar, ::UseThreads, p::Params) where {T<:SecondChernAlgorithms}
#     Base.Threads.@threads for i in idxs
#         for j in eachindex(range2)
#             update2Din!(T(), i, j, nums, num0, H, alg!, range1, range2, p)
#         end
#     end
# end

# function update1D!(nums, num0, H, alg!, range1, ::UseThreads, p::Params)
#     Base.Threads.@threads for i in eachindex(range1)
#         update1Din!(i, nums, num0, H, alg!, range1, p)
#     end
# end

# function update1D!(nums, num0, H, alg!, range1, idxs::ProgressBar, ::UseThreads, p::Params)
#     Base.Threads.@threads for i in idxs
#         update1Din!(i, nums, num0, H, alg!, range1, p)
#     end
# end

# function update2D!(nums, num0, H, alg!, range1, range2, ::UseThreads, p::Params)
#     Base.Threads.@threads for i in eachindex(range1)
#         for j in eachindex(range2)
#             update2Din!(i, j, nums, num0, H, alg!, range1, range2, p)
#         end
#     end
# end

# function update2D!(nums, num0, H, alg!, range1, range2, idxs::ProgressBar, ::UseThreads, p::Params)
#     Base.Threads.@threads for i in idxs
#         for j in eachindex(range2)
#             update2Din!(i, j, nums, num0, H, alg!, range1, range2, p)
#         end
#     end
# end

function update1D!(
    ::T, nums, num0, H, alg!, range1, mod::UseMPI, p::Params
) where {T<:TopologicalNumbersAlgorithms}
    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(range1), nprocs)[myrank + 1]

    for i in idxs
        update1Din!(T(), i, nums, num0, H, alg!, range1, p)
    end

    nums .= mod.MPI.Allreduce!(nums, mod.MPI.SUM, comm)

    return mod.MPI.Barrier(comm)
end

function update1D!(
    ::T, nums, num0, H, alg!, range1, idxs::ProgressBar, mod::UseMPI, p::Params
) where {T<:TopologicalNumbersAlgorithms}
    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(range1), nprocs)[myrank + 1]
    if myrank == 0
        idxs = ProgressBar(idxs)
    end

    for i in idxs
        update1Din!(T(), i, nums, num0, H, alg!, range1, p)
    end

    nums .= mod.MPI.Allreduce!(nums, mod.MPI.SUM, comm)

    return mod.MPI.Barrier(comm)
end

function update2D!(
    ::T, nums, num0, H, alg!, range1, range2, mod::UseMPI, p::Params
) where {T<:TopologicalNumbersAlgorithms}
    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(range1), nprocs)[myrank + 1]

    for i in idxs
        for j in eachindex(range2)
            update2Din!(T(), i, j, nums, num0, H, alg!, range1, range2, p)
        end
    end

    nums .= mod.MPI.Allreduce!(nums, mod.MPI.SUM, comm)

    return mod.MPI.Barrier(comm)
end

function update2D!(
    ::T, nums, num0, H, alg!, range1, range2, idxs::ProgressBar, mod::UseMPI, p::Params
) where {T<:TopologicalNumbersAlgorithms}
    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(range1), nprocs)[myrank + 1]
    if myrank == 0
        idxs = ProgressBar(idxs)
    end

    for i in idxs
        for j in eachindex(range2)
            update2Din!(T(), i, j, nums, num0, H, alg!, range1, range2, p)
        end
    end

    nums .= mod.MPI.Allreduce!(nums, mod.MPI.SUM, comm)

    return mod.MPI.Barrier(comm)
end

# Old method
function update1D!(nums, num0, H, alg!, range1, mod::UseMPI, p::Params)
    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(range1), nprocs)[myrank + 1]

    for i in idxs
        update1Din!(i, nums, num0, H, alg!, range1, p)
    end

    nums .= mod.MPI.Allreduce!(nums, mod.MPI.SUM, comm)

    return mod.MPI.Barrier(comm)
end

# Old method
function update1D!(nums, num0, H, alg!, range1, idxs::ProgressBar, mod::UseMPI, p::Params)
    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(range1), nprocs)[myrank + 1]
    if myrank == 0
        idxs = ProgressBar(idxs)
    end

    for i in idxs
        update1Din!(i, nums, num0, H, alg!, range1, p)
    end

    nums .= mod.MPI.Allreduce!(nums, mod.MPI.SUM, comm)

    return mod.MPI.Barrier(comm)
end

# Old method
function update2D!(nums, num0, H, alg!, range1, range2, mod::UseMPI, p::Params)
    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(range1), nprocs)[myrank + 1]

    for i in idxs
        for j in eachindex(range2)
            update2Din!(i, j, nums, num0, H, alg!, range1, range2, p)
        end
    end

    nums .= mod.MPI.Allreduce!(nums, mod.MPI.SUM, comm)

    return mod.MPI.Barrier(comm)
end

# Old method
function update2D!(
    nums, num0, H, alg!, range1, range2, idxs::ProgressBar, mod::UseMPI, p::Params
)
    mod.MPI.Init()
    comm = mod.MPI.COMM_WORLD
    myrank = mod.MPI.Comm_rank(comm)
    nprocs = mod.MPI.Comm_size(comm)

    idxs = Distributed.splitrange(1, length(range1), nprocs)[myrank + 1]
    if myrank == 0
        idxs = ProgressBar(idxs)
    end

    for i in idxs
        for j in eachindex(range2)
            update2Din!(i, j, nums, num0, H, alg!, range1, range2, p)
        end
    end

    nums .= mod.MPI.Allreduce!(nums, mod.MPI.SUM, comm)

    return mod.MPI.Barrier(comm)
end

function calc_data1D(
    H, param_range, alg::T1, parallel::T2, p::Params
) where {T1<:TopologicalNumbersAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack rounds, Hs, returnRealValue = p

    nums = zeros(Float64, Hs, length(param_range))
    num0 = zeros(Float64, Hs)

    if alg isa FHS2
        algorithm! = SecondChernPhase!
        type0 = ifelse(returnRealValue, Float64, ComplexF64)
        nums = zeros(ComplexF64, length(param_range))
        num0 = zero(type0)
    elseif alg isa BP
        algorithm! = BerryPhase!
    elseif alg isa Shio
        algorithm! = Z2Phase!
        nums = zeros(Float64, 2, length(param_range))
        num0 = zeros(Float64, 2)
    elseif alg isa FHS
        algorithm! = ChernPhase!
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update1D!(alg, nums, num0, H, algorithm!, param_range, parallel, p)

    if alg isa FHS2
        if returnRealValue == true
            real.(nums)
        elseif returnRealValue == false
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

function calc_data1D(
    H, param_range, idxs::ProgressBar, alg::T1, parallel::T2, p::Params
) where {T1<:TopologicalNumbersAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack rounds, Hs, returnRealValue = p

    nums = zeros(Float64, Hs, length(param_range))
    num0 = zeros(Float64, Hs)

    if alg isa FHS2
        algorithm! = SecondChernPhase!
        type0 = ifelse(returnRealValue, Float64, ComplexF64)
        nums = zeros(ComplexF64, length(param_range))
        num0 = zero(type0)
    elseif alg isa BP
        algorithm! = BerryPhase!
    elseif alg isa Shio
        algorithm! = Z2Phase!
        nums = zeros(Float64, 2, length(param_range))
        num0 = zeros(Float64, 2)
    elseif alg isa FHS
        algorithm! = ChernPhase!
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update1D!(alg, nums, num0, H, algorithm!, param_range, idxs, parallel, p)

    if alg isa FHS2
        if returnRealValue == true
            real.(nums)
        elseif returnRealValue == false
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

function calc_data2D(
    H, param_range1, param_range2, alg::T1, parallel::T2, p::Params
) where {T1<:TopologicalNumbersAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack rounds, Hs, returnRealValue = p

    nums = zeros(Float64, Hs, length(param_range1), length(param_range2))
    num0 = zeros(Float64, Hs)

    if alg isa FHS2
        algorithm! = SecondChernPhase!
        type0 = ifelse(returnRealValue, Float64, ComplexF64)
        nums = zeros(ComplexF64, length(param_range1), length(param_range2))
        num0 = zero(type0)
    elseif alg isa BP
        algorithm! = BerryPhase!
    elseif alg isa Shio
        algorithm! = Z2Phase!
        nums = zeros(Float64, 2, length(param_range1), length(param_range2))
        num0 = zeros(Float64, 2)
    elseif alg isa FHS
        algorithm! = ChernPhase!
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update2D!(alg, nums, num0, H, algorithm!, param_range1, param_range2, parallel, p)

    if alg isa FHS2
        if returnRealValue == true
            real.(nums)
        elseif returnRealValue == false
            nums
        end
    else
        if rounds == true
            if all(!isnan, nums)
                nums = round.(Int, nums)
            else
                for i in eachindex(nums)
                    if nums[i] !== NaN
                        nums[i] = round(Int, nums[i])
                    end
                end
            end
            nums
        elseif rounds == false
            nums
        end
    end
end

function calc_data2D(
    H, param_range1, param_range2, idxs::ProgressBar, alg::T1, parallel::T2, p::Params
) where {T1<:TopologicalNumbersAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack rounds, Hs, returnRealValue = p

    nums = zeros(Float64, Hs, length(param_range1), length(param_range2))
    num0 = zeros(Float64, Hs)

    if alg isa FHS2
        algorithm! = SecondChernPhase!
        type0 = ifelse(returnRealValue, Float64, ComplexF64)
        nums = zeros(ComplexF64, length(param_range1), length(param_range2))
        num0 = zero(type0)
    elseif alg isa BP
        algorithm! = BerryPhase!
    elseif alg isa Shio
        algorithm! = Z2Phase!
        nums = zeros(Float64, 2, length(param_range1), length(param_range2))
        num0 = zeros(Float64, 2)
    elseif alg isa FHS
        algorithm! = ChernPhase!
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end
    update2D!(alg, nums, num0, H, algorithm!, param_range1, param_range2, idxs, parallel, p)

    if alg isa FHS2
        if returnRealValue == true
            real.(nums)
        elseif returnRealValue == false
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
function calc_data1D(
    H, param_range, alg::String, parallel::T, p::Params
) where {T<:TopologicalNumbersParallel}
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range))
    num0 = zeros(Float64, Hs)

    if alg == "BerryPhase"
        algorithm! = BerryPhase!
        update1D!(nums, num0, H, algorithm!, param_range, parallel, p)
    elseif alg == "Z2"
        algorithm! = Z2Phase!
        nums = zeros(Float64, 2, length(param_range))
        num0 = zeros(Float64, 2)
        @reset p.Nfill = Hs ÷ 2
        update1D!(Shio(), nums, num0, H, algorithm!, param_range, parallel, p)
    elseif alg == "Chern"
        algorithm! = ChernPhase!
        update1D!(nums, num0, H, algorithm!, param_range, parallel, p)
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end

    if rounds == true
        round.(Int, transpose(nums))
    elseif rounds == false
        transpose(nums)
    end
end

# Old method
function calc_data1D(
    H, param_range, idxs::ProgressBar, alg::String, parallel::T, p::Params
) where {T<:TopologicalNumbersParallel}
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range))
    num0 = zeros(Float64, Hs)

    if alg == "BerryPhase"
        algorithm! = BerryPhase!
        update1D!(nums, num0, H, algorithm!, param_range, idxs, parallel, p)
    elseif alg == "Z2"
        algorithm! = Z2Phase!
        nums = zeros(Float64, 2, length(param_range))
        num0 = zeros(Float64, 2)
        @reset p.Nfill = Hs ÷ 2
        update1D!(Shio(), nums, num0, H, algorithm!, param_range, idxs, parallel, p)
    elseif alg == "Chern"
        algorithm! = ChernPhase!
        update1D!(nums, num0, H, algorithm!, param_range, idxs, parallel, p)
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end

    if rounds == true
        round.(Int, transpose(nums))
    elseif rounds == false
        transpose(nums)
    end
end

# Old method
function calc_data2D(
    H, param_range1, param_range2, alg::String, parallel::T, p::Params
) where {T<:TopologicalNumbersParallel}
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range1), length(param_range2))
    num0 = zeros(Float64, Hs)

    if alg == "BerryPhase"
        algorithm! = BerryPhase!
        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, parallel, p)
    elseif alg == "Z2"
        algorithm! = Z2Phase!
        nums = zeros(Float64, 2, length(param_range1), length(param_range2))
        num0 = zeros(Float64, 2)
        @reset p.Nfill = Hs ÷ 2
        update2D!(
            Shio(), nums, num0, H, algorithm!, param_range1, param_range2, parallel, p
        )
    elseif alg == "Chern"
        algorithm! = ChernPhase!
        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, parallel, p)
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end

    if rounds == true
        round.(Int, nums)
    elseif rounds == false
        nums
    end
end

# Old method
function calc_data2D(
    H, param_range1, param_range2, idxs::ProgressBar, alg::String, parallel::T, p::Params
) where {T<:TopologicalNumbersParallel}
    @unpack rounds, Hs = p

    nums = zeros(Float64, Hs, length(param_range1), length(param_range2))
    num0 = zeros(Float64, Hs)

    if alg == "BerryPhase"
        algorithm! = BerryPhase!
        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, idxs, parallel, p)
    elseif alg == "Z2"
        algorithm! = Z2Phase!
        nums = zeros(Float64, 2, length(param_range1), length(param_range2))
        num0 = zeros(Float64, 2)
        @reset p.Nfill = Hs ÷ 2
        update2D!(
            Shio(), nums, num0, H, algorithm!, param_range1, param_range2, idxs, parallel, p
        )
    elseif alg == "Chern"
        algorithm! = ChernPhase!
        update2D!(nums, num0, H, algorithm!, param_range1, param_range2, idxs, parallel, p)
    else
        throw(ArgumentError("Unknown algorithm $alg"))
    end

    if rounds == true
        round.(Int, nums)
    elseif rounds == false
        nums
    end
end

# Old method
@doc raw"""
    calcPhaseDiagram(H::Function, param_range::T1, alg::T2; N::T3=51, parallel::T4=UseSingleThread(), gapless::Real=0.0, rounds::Bool=true, plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:Union{String,TopologicalNumbersAlgorithms},T3<:Union{Int,Tuple,AbstractVector},T4<:TopologicalNumbersParallel}"""
function calcPhaseDiagram(
    H::Function,
    param_range::T1,
    alg::T2;
    N::T3=51,
    parallel::T4=UseSingleThread(),
    gapless::Real=0.0,
    rounds::Bool=true,
    plot::Bool=false,
    progress::Bool=false,
) where {
    T1<:AbstractVector,
    T2<:Union{String,TopologicalNumbersAlgorithms},
    T3<:Union{Int,Tuple,AbstractVector},
    T4<:TopologicalNumbersParallel,
}
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

    p = Params(; Ham=Hamiltonian, dim, N, gapless, rounds, Hs)
    if alg == FHS2()
        @reset p.N = (N, N, N, N)
        @reset p.Nfill = Hs ÷ 2
    end

    nums = if progress == true
        calc_data1D(H, param_range, ProgressBar(eachindex(param_range)), alg, parallel, p)
    else
        calc_data1D(H, param_range, alg, parallel, p)
    end

    if plot == true
        if alg == FHS2()
            plot1D(alg, nums, param_range)
        else
            plot1D(nums, param_range)
        end
    end

    return (; param=param_range, nums)
end

# Old method
@doc raw"""
    calcPhaseDiagram(H::Function, param_range1::T1, param_range2::T2, alg::T3; N::T4=51, parallel::T5=UseSingleThread(), gapless::Real=0.0, rounds::Bool=true, plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:Union{String,TopologicalNumbersAlgorithms},T4<:Union{Int,Tuple,AbstractVector},T5<:TopologicalNumbersParallel}
"""
function calcPhaseDiagram(
    H::Function,
    param_range1::T1,
    param_range2::T2,
    alg::T3;
    N::T4=51,
    parallel::T5=UseSingleThread(),
    gapless::Real=0.0,
    rounds::Bool=true,
    plot::Bool=false,
    progress::Bool=false,
) where {
    T1<:AbstractVector,
    T2<:AbstractVector,
    T3<:Union{String,TopologicalNumbersAlgorithms},
    T4<:Union{Int,Tuple,AbstractVector},
    T5<:TopologicalNumbersParallel,
}
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

    p = Params(; Ham=Hamiltonian, dim, N, gapless, rounds, Hs)
    if alg == FHS2()
        @reset p.N = (N, N, N, N)
        @reset p.Nfill = Hs ÷ 2
    end

    nums = if progress == true
        calc_data2D(
            H,
            param_range1,
            param_range2,
            ProgressBar(eachindex(param_range1)),
            alg,
            parallel,
            p,
        )
    else
        calc_data2D(H, param_range1, param_range2, alg, parallel, p)
    end

    if alg == FHS2()
        if plot == true
            plot2D(transpose(nums_half), param_range1, param_range2)
        end
    else
        if plot == true && Hs % 2 == 0
            nums_half = sum(@view(nums[1:(end ÷ 2), :, :]); dims=1)[1, :, :]
            plot2D(transpose(nums_half), param_range1, param_range2) # half-filling case
        end
    end

    return (; param1=param_range1, param2=param_range2, nums)
end

function calcPhaseDiagram1D_core(H, param_range, alg, p; parallel, plot, progress)
    nums = if progress == true
        calc_data1D(H, param_range, ProgressBar(eachindex(param_range)), alg, parallel, p)
    else
        calc_data1D(H, param_range, alg, parallel, p)
    end

    if plot == true
        plot1D(nums, param_range)
    end
    return nums
end

function calcPhaseDiagram2D_core(
    H, param_range1, param_range2, alg, p; parallel, plot, progress
)
    nums = if progress == true
        calc_data2D(
            H,
            param_range1,
            param_range2,
            ProgressBar(eachindex(param_range1)),
            alg,
            parallel,
            p,
        )
    else
        calc_data2D(H, param_range1, param_range2, alg, parallel, p)
    end

    if alg isa FHS2
        if plot == true
            plot2D(transpose(nums_half), param_range1, param_range2)
        end
    else
        if plot == true && p.Hs % 2 == 0
            nums_half = sum(@view(nums[1:(end ÷ 2), :, :]); dims=1)[1, :, :]
            plot2D(transpose(nums_half), param_range1, param_range2) # half-filling case
        end
    end
    return nums
end

# Berry phase
@doc raw"""
Calculate the one-dimensional phase diagram for the Berry phase over a specified parameter range.

    calcPhaseDiagram(prob::BPProblem, param_range::T1, alg::T2=BP(); parallel::T3=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:BerryPhaseAlgorithms,T3<:TopologicalNumbersParallel}

# Description
`calcPhaseDiagram` calculates the Berry phase for a one-dimensional system described by `BPProblem` over a range of parameters (`param_range`), constructing a phase diagram. By changing the `alg` or `parallel` options, you can flexibly switch between different calculation methods and parallelization strategies. Moreover, setting `plot=true` generates a simple plot of the results upon completion of the calculation.

# Arguments
- `prob::BPProblem`: 
  A structure that contains the problem settings for Berry phase calculations, including the Hamiltonian, mesh size, gap information, and other parameters.

- `param_range::AbstractVector`: 
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the Berry phase.

- `alg::BerryPhaseAlgorithms=BP()`: 
  The algorithm used to compute the Berry phase. The default is `BP()`. You can implement and specify other algorithms as needed.

- `parallel::TopologicalNumbersParallel=UseSingleThread()`: 
  The parallelization strategy. The default is single-threaded (`UseSingleThread()`). If you wish to run in a distributed environment, you can specify it accordingly. (Note that, at present, the Berry phase calculation does not support multithreading.)

- `plot::Bool=false`: 
  If set to `true`, a simple plot of the calculated results is shown after the computation finishes. This is helpful for visually confirming topological phase transitions.

- `progress::Bool=false`: 
  If set to `true`, progress is displayed during loop calculations, which can be useful for large-scale computations.

# Returns
- `NamedTuple{(:param, :nums)}`: 
  - `param::AbstractVector`: A copy of the input `param_range`.
  - `nums::AbstractMatrix`: A matrix containing the computed Berry phases for each parameter. The size of this matrix is (**number of bands** × **number of parameter values**).

# Examples

```julia
julia> using TopologicalNumbers
julia> # Define the Hamiltonian
       H₀(k, p) = SSH(k, p)
julia> # Set up the Berry phase problem
       prob = BPProblem(H₀)
julia> # Define a parameter range from -2.0 to 2.0 with 9 divisions
       param_range = range(-2.0, 2.0, length=9)
julia> # Calculate the phase diagram (using the default BP algorithm, single-threaded,
       # with plotting enabled and progress display)
       result = calcPhaseDiagram(prob, param_range; plot=true, progress=true)
(param = -2.0:0.5:2.0, nums = [1 1; 1 1; 0 0; 0 0; 0 0; 0 0; 0 0; 1 1; 1 1])

julia> result.param
-2.0:0.5:2.0

julia> result.nums
9×2 Matrix{Int64}:
 1  1
 1  1
 0  0
 0  0
 0  0
 0  0
 0  0
 1  1
 1  1
```

"""
function calcPhaseDiagram(
    prob::BPProblem,
    param_range::T1,
    alg::T2=BP();
    parallel::T3=UseSingleThread(),
    plot::Bool=false,
    progress::Bool=false,
) where {T1<:AbstractVector,T2<:BerryPhaseAlgorithms,T3<:TopologicalNumbersParallel}
    @unpack H, N, gapless, rounds = prob

    dim = 1
    Ham(k) = H(k, 0.0)
    Hs = size(Ham(zero(Float64)), 1)

    p = Params(; Ham, dim, N, gapless, rounds, Hs)

    nums = calcPhaseDiagram1D_core(H, param_range, alg, p; parallel, plot, progress)

    return (; param=param_range, nums)
end

@doc raw"""
Calculate the two-dimensional phase diagram for the Berry phase over a specified parameter range.

    calcPhaseDiagram(prob::BPProblem, param_range1::T1, param_range2::T2, alg::T3=BP(); parallel::T4=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:BerryPhaseAlgorithms,T4<:TopologicalNumbersParallel}

# Description
`calcPhaseDiagram` calculates the Berry phase for a two-dimensional system described by `BPProblem` over a range of parameters (`param_range`), constructing a phase diagram. By changing the `alg` or `parallel` options, you can flexibly switch between different calculation methods and parallelization strategies. Moreover, setting `plot=true` generates a simple plot of the results upon completion of the calculation.

# Arguments
- `prob::BPProblem`: 
  A structure that contains the problem settings for Berry phase calculations, including the Hamiltonian, mesh size, gap information, and other parameters.

- `param_range1::AbstractVector`: 
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the Berry phase.

- `param_range2::AbstractVector`:
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the Berry phase.

- `alg::BerryPhaseAlgorithms=BP()`: 
  The algorithm used to compute the Berry phase. The default is `BP()`. You can implement and specify other algorithms as needed.

- `parallel::TopologicalNumbersParallel=UseSingleThread()`: 
  The parallelization strategy. The default is single-threaded (`UseSingleThread()`). If you wish to run in a distributed environment, you can specify it accordingly. (Note that, at present, the Berry phase calculation does not support multithreading.)

- `plot::Bool=false`: 
  If set to `true`, a simple plot of the calculated results is shown after the computation finishes. This is helpful for visually confirming topological phase transitions.

- `progress::Bool=false`: 
  If set to `true`, progress is displayed during loop calculations, which can be useful for large-scale computations.

# Returns
- `NamedTuple{(:param, :nums)}`: 
  - `param1::AbstractVector`: A copy of the input `param_range1`.
  - `param2::AbstractVector`: A copy of the input `param_range2`.
  - `nums::AbstractArray`: An array containing the computed Berry phases for each parameter. The size of this array is (**number of bands** × **number of param1** × **number of param2**).

# Examples

```julia
julia> using TopologicalNumbers
julia> # Define the Hamiltonian
       H₀(k, p) = KitaevChain(k, p)
julia> # Set up the Berry phase problem
       prob = BPProblem(H₀)
julia> # Define a parameter range
       param1 = range(-3.0, 3.0, length=4)
       param2 = range(-1.0, 1.0, length=4)
julia> # Calculate the phase diagram (using the default BP algorithm, single-threaded,
       # with plotting enabled and progress display)
       result = calcPhaseDiagram(prob, param1, param2; plot=true, progress=true)
(param1 = -3.0:2.0:3.0, param2 = -1.0:0.6666666666666666:1.0, nums = [0 1 1 0; 0 1 1 0;;; 0 1 1 0; 0 1 1 0;;; 0 1 1 0; 0 1 1 0;;; 0 1 1 0; 0 1 1 0])

julia> result.param1
-3.0:2.0:3.0

julia> result.param2
-1.0:0.6666666666666666:1.0

julia> result.nums
2×4×4 Array{Int64, 3}:
[:, :, 1] =
 0  1  1  0
 0  1  1  0

[:, :, 2] =
 0  1  1  0
 0  1  1  0

[:, :, 3] =
 0  1  1  0
 0  1  1  0

[:, :, 4] =
 0  1  1  0
 0  1  1  0
```

"""
function calcPhaseDiagram(
    prob::BPProblem,
    param_range1::T1,
    param_range2::T2,
    alg::T3=BP();
    parallel::T4=UseSingleThread(),
    plot::Bool=false,
    progress::Bool=false,
) where {
    T1<:AbstractVector,
    T2<:AbstractVector,
    T3<:BerryPhaseAlgorithms,
    T4<:TopologicalNumbersParallel,
}
    @unpack H, N, gapless, rounds = prob

    dim = 1
    Ham(k) = H(k, (0.0, 0.0))
    Hs = size(Ham(zero(Float64)), 1)

    p = Params(; Ham, dim, N, gapless, rounds, Hs)

    nums = calcPhaseDiagram2D_core(
        H, param_range1, param_range2, alg, p; parallel, plot, progress
    )

    return (; param1=param_range1, param2=param_range2, nums)
end

# First Chern number
@doc raw"""
Calculate the one-dimensional phase diagram for the first Chern number over a specified parameter range.

    calcPhaseDiagram(prob::FCProblem, param_range::T1, alg::T2=FHS(); parallel::T3=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:FirstChernAlgorithms,T3<:TopologicalNumbersParallel}

# Description
`calcPhaseDiagram` calculates the first Chern number for a one-dimensional system described by `FCProblem` over a range of parameters (`param_range`), constructing a phase diagram. By changing the `alg` or `parallel` options, you can flexibly switch between different calculation methods and parallelization strategies. Moreover, setting `plot=true` generates a simple plot of the results upon completion of the calculation.

# Arguments
- `prob::FCProblem`: 
  A structure that contains the problem settings for first Chern number calculations, including the Hamiltonian, mesh size, gap information, and other parameters.

- `param_range::AbstractVector`: 
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the first Chern number.

- `alg::FirstChernAlgorithms=FHS()`: 
  The algorithm used to compute the first Chern number. The default is `FHS()`. You can implement and specify other algorithms as needed.

- `parallel::TopologicalNumbersParallel=UseSingleThread()`: 
  The parallelization strategy. The default is single-threaded (`UseSingleThread()`). If you wish to run in a distributed environment, you can specify it accordingly. (Note that, at present, the first Chern number calculation does not support multithreading.)

- `plot::Bool=false`: 
  If set to `true`, a simple plot of the calculated results is shown after the computation finishes. This is helpful for visually confirming topological phase transitions.

- `progress::Bool=false`: 
  If set to `true`, progress is displayed during loop calculations, which can be useful for large-scale computations.

# Returns
- `NamedTuple{(:param, :nums)}`: 
  - `param::AbstractVector`: A copy of the input `param_range`.
  - `nums::AbstractMatrix`: A matrix containing the computed first Chern numbers for each parameter. The size of this matrix is (**number of bands** × **number of parameter values**).

# Examples

```julia
julia> using TopologicalNumbers
julia> # Define the Hamiltonian
       H₀(k, p) = Haldane(k, p)
       H(k, p) = H₀(k, (1, p, 2.5))
julia> # Set up the first Chern number problem
       prob = FCProblem(H)
julia> # Define a parameter range
       param = range(-π, π, length=10);
julia> # Calculate the phase diagram (using the default FHS algorithm, single-threaded,
       # with plotting enabled and progress display)
       result = calcPhaseDiagram(prob, param; plot=true, progress=true)
(param = -3.141592653589793:0.6981317007977318:3.141592653589793, nums = [0 0; -1 1; -1 1; -1 1; 0 0; 0 0; 1 -1; 1 -1; 1 -1; 0 0])

julia> result.param
-3.141592653589793:0.6981317007977318:3.141592653589793

julia> result.nums
10×2 Matrix{Int64}:
  0   0
 -1   1
 -1   1
 -1   1
  0   0
  0   0
  1  -1
  1  -1
  1  -1
  0   0
```

"""
function calcPhaseDiagram(
    prob::FCProblem,
    param_range::T1,
    alg::T2=FHS();
    parallel::T3=UseSingleThread(),
    plot::Bool=false,
    progress::Bool=false,
) where {T1<:AbstractVector,T2<:FirstChernAlgorithms,T3<:TopologicalNumbersParallel}
    @unpack H, N, gapless, rounds = prob

    dim = 2
    Ham(k) = H(k, 0.0)
    Hs = size(Ham(zeros(2)), 1)

    p = Params(; Ham, dim, N, gapless, rounds, Hs)

    nums = calcPhaseDiagram1D_core(H, param_range, alg, p; parallel, plot, progress)

    return (; param=param_range, nums)
end

@doc raw"""
Calculate the two-dimensional phase diagram for the first Chern number over a specified parameter range.

    calcPhaseDiagram(prob::FCProblem, param_range1::T1, param_range2::T2, alg::T3=FHS(); parallel::T4=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:FirstChernAlgorithms,T4<:TopologicalNumbersParallel}

# Description
`calcPhaseDiagram` calculates the first Chern number for a two-dimensional system described by `FCProblem` over a range of parameters (`param_range`), constructing a phase diagram. By changing the `alg` or `parallel` options, you can flexibly switch between different calculation methods and parallelization strategies. Moreover, setting `plot=true` generates a simple plot of the results upon completion of the calculation.

# Arguments
- `prob::FCProblem`: 
  A structure that contains the problem settings for first Chern number calculations, including the Hamiltonian, mesh size, gap information, and other parameters.

- `param_range1::AbstractVector`: 
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the first Chern number.

- `param_range2::AbstractVector`:
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the first Chern number.

- `alg::FirstChernAlgorithms=FHS()`: 
  The algorithm used to compute the first Chern number. The default is `FHS()`. You can implement and specify other algorithms as needed.

- `parallel::TopologicalNumbersParallel=UseSingleThread()`: 
  The parallelization strategy. The default is single-threaded (`UseSingleThread()`). If you wish to run in a distributed environment, you can specify it accordingly. (Note that, at present, the first Chern number calculation does not support multithreading.)

- `plot::Bool=false`: 
  If set to `true`, a simple plot of the calculated results is shown after the computation finishes. This is helpful for visually confirming topological phase transitions.

- `progress::Bool=false`: 
  If set to `true`, progress is displayed during loop calculations, which can be useful for large-scale computations.

# Returns
- `NamedTuple{(:param, :nums)}`: 
  - `param1::AbstractVector`: A copy of the input `param_range1`.
  - `param2::AbstractVector`: A copy of the input `param_range2`.
  - `nums::AbstractArray`: An array containing the computed first Chern numbers for each parameter. The size of this array is (**number of bands** × **number of param1** × **number of param2**).

# Examples

```julia
julia> using TopologicalNumbers
julia> # Define the Hamiltonian
       H₀(k, p) = Haldane(k, p)
       H(k, p) = H₀(k, (1, p[1], p[2]))
julia> # Set up the first Chern number problem
       prob = FCProblem(H₀)
julia> # Define a parameter range
       param1 = range(-π, π, length=4)
       param2 = range(-6.0, 6.0, length=4)
julia> # Calculate the phase diagram (using the default FHS algorithm, single-threaded,
       # with plotting enabled and progress display)
       result = calcPhaseDiagram(prob, param1, param2; plot=true, progress=true)
(param1 = -3.141592653589793:2.0943951023931953:3.141592653589793, param2 = -6.0:4.0:6.0, nums = [0 0 0 0; 0 0 0 0;;; 0 -1 1 0; 0 1 -1 0;;; 0 -1 1 0; 0 1 -1 0;;; 0 0 0 0; 0 0 0 0])

julia> result.param1
-3.141592653589793:2.0943951023931953:3.141592653589793

julia> result.param2
-6.0:4.0:6.0

julia> result.nums
2×4×4 Array{Int64, 3}:
[:, :, 1] =
 0  0  0  0
 0  0  0  0

[:, :, 2] =
 0  -1   1  0
 0   1  -1  0

[:, :, 3] =
 0  -1   1  0
 0   1  -1  0

[:, :, 4] =
 0  0  0  0
 0  0  0  0
```

"""
function calcPhaseDiagram(
    prob::FCProblem,
    param_range1::T1,
    param_range2::T2,
    alg::T3=FHS();
    parallel::T4=UseSingleThread(),
    plot::Bool=false,
    progress::Bool=false,
) where {
    T1<:AbstractVector,
    T2<:AbstractVector,
    T3<:FirstChernAlgorithms,
    T4<:TopologicalNumbersParallel,
}
    @unpack H, N, gapless, rounds = prob

    dim = 2
    Ham(k) = H(k, (0.0, 0.0))
    Hs = size(Ham(zeros(2)), 1)

    p = Params(; Ham, dim, N, gapless, rounds, Hs)

    nums = calcPhaseDiagram2D_core(
        H, param_range1, param_range2, alg, p; parallel, plot, progress
    )

    return (; param1=param_range1, param2=param_range2, nums)
end

# Z2
@doc raw"""
Calculate the one-dimensional phase diagram for the $\mathbb{Z}_2$ number over a specified parameter range.

    calcPhaseDiagram(prob::Z2Problem, param_range::T1, alg::T2=Shio(); parallel::T3=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:Z2Algorithms,T3<:TopologicalNumbersParallel}

# Description
`calcPhaseDiagram` calculates the $\mathbb{Z}_2$ number for a one-dimensional system described by `Z2Problem` over a range of parameters (`param_range`), constructing a phase diagram. By changing the `alg` or `parallel` options, you can flexibly switch between different calculation methods and parallelization strategies. Moreover, setting `plot=true` generates a simple plot of the results upon completion of the calculation.

# Arguments
- `prob::Z2Problem`: 
  A structure that contains the problem settings for $\mathbb{Z}_2$ number calculations, including the Hamiltonian, mesh size, gap information, and other parameters.

- `param_range::AbstractVector`: 
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the $\mathbb{Z}_2$ number.

- `alg::Z2Algorithms=Shio()`: 
  The algorithm used to compute the $\mathbb{Z}_2$ number. The default is `Shio()`. You can implement and specify other algorithms as needed.

- `parallel::TopologicalNumbersParallel=UseSingleThread()`: 
  The parallelization strategy. The default is single-threaded (`UseSingleThread()`). If you wish to run in a distributed environment, you can specify it accordingly. (Note that, at present, the $\mathbb{Z}_2$ number calculation does not support multithreading.)

- `plot::Bool=false`: 
  If set to `true`, a simple plot of the calculated results is shown after the computation finishes. This is helpful for visually confirming topological phase transitions.

- `progress::Bool=false`: 
  If set to `true`, progress is displayed during loop calculations, which can be useful for large-scale computations.

# Returns
- `NamedTuple{(:param, :nums)}`: 
  - `param::AbstractVector`: A copy of the input `param_range`.
  - `nums::AbstractMatrix`: A matrix containing the computed $\mathbb{Z}_2$ numbers for each parameter. The size of this matrix is (**number of bands** × **number of parameter values**).

# Examples

```julia
julia> using TopologicalNumbers
julia> # Define the Hamiltonian
       H₀(k, p) = BHZ(k, p)
       H(k, p) = H₀(k, (p, 0.25))
julia> # Set up the $\mathbb{Z}_2$ number problem
       prob = Z2Problem(H)
julia> # Define a parameter range
       param = range(-2, 2, length=10)
julia> # Calculate the phase diagram (using the default Shio algorithm, single-threaded,
       # with plotting enabled and progress display)
       result = calcPhaseDiagram(prob, param; plot=true, progress=true)
(param = -2.0:0.4444444444444444:2.0, nums = [0 0; 0 0; 0 0; 1 1; 1 1; 1 1; 1 1; 0 0; 0 0; 0 0])

julia> result.param
-2.0:0.4444444444444444:2.0

julia> result.nums
10×2 Matrix{Int64}:
 0  0
 0  0
 0  0
 1  1
 1  1
 1  1
 1  1
 0  0
 0  0
 0  0
```

"""
function calcPhaseDiagram(
    prob::Z2Problem,
    param_range::T1,
    alg::T2=Shio();
    parallel::T3=UseSingleThread(),
    plot::Bool=false,
    progress::Bool=false,
) where {T1<:AbstractVector,T2<:Z2Algorithms,T3<:TopologicalNumbersParallel}
    @unpack H, Nfill, N, rounds = prob

    dim = 2
    Ham(k) = H(k, 0.0)
    Hs = size(Ham(zeros(2)), 1)

    if isnothing(Nfill)
        Nfill = Hs ÷ 2
    end

    p = Params(; Ham, dim, Nfill, N, rounds, Hs)

    nums = calcPhaseDiagram1D_core(H, param_range, alg, p; parallel, plot, progress)

    return (; param=param_range, nums)
end

@doc raw"""
Calculate the two-dimensional phase diagram for the $\mathbb{Z}_2$ number over a specified parameter range.

    calcPhaseDiagram(prob::Z2Problem, param_range1::T1, param_range2::T2, alg::T3=Shio(); parallel::T4=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:Z2Algorithms,T4<:TopologicalNumbersParallel}

# Description
`calcPhaseDiagram` calculates the $\mathbb{Z}_2$ number for a two-dimensional system described by `Z2Problem` over a range of parameters (`param_range`), constructing a phase diagram. By changing the `alg` or `parallel` options, you can flexibly switch between different calculation methods and parallelization strategies. Moreover, setting `plot=true` generates a simple plot of the results upon completion of the calculation.

# Arguments
- `prob::Z2Problem`: 
  A structure that contains the problem settings for $\mathbb{Z}_2$ number calculations, including the Hamiltonian, mesh size, gap information, and other parameters.

- `param_range1::AbstractVector`: 
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the $\mathbb{Z}_2$ number.

- `param_range2::AbstractVector`:
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the $\mathbb{Z}_2$ number.

- `alg::Z2Algorithms=Shio()`: 
  The algorithm used to compute the $\mathbb{Z}_2$ number. The default is `Shio()`. You can implement and specify other algorithms as needed.

- `parallel::TopologicalNumbersParallel=UseSingleThread()`: 
  The parallelization strategy. The default is single-threaded (`UseSingleThread()`). If you wish to run in a distributed environment, you can specify it accordingly. (Note that, at present, the $\mathbb{Z}_2$ number calculation does not support multithreading.)

- `plot::Bool=false`: 
  If set to `true`, a simple plot of the calculated results is shown after the computation finishes. This is helpful for visually confirming topological phase transitions.

- `progress::Bool=false`: 
  If set to `true`, progress is displayed during loop calculations, which can be useful for large-scale computations.

# Returns
- `NamedTuple{(:param, :nums)}`: 
  - `param1::AbstractVector`: A copy of the input `param_range1`.
  - `param2::AbstractVector`: A copy of the input `param_range2`.
  - `nums::AbstractArray`: An array containing the computed $\mathbb{Z}_2$ numbers for each parameter. The size of this array is (**number of bands** × **number of param1** × **number of param2**).

# Examples

```julia
julia> using TopologicalNumbers
julia> # Define the Hamiltonian
       H₀(k, p) = BHZ(k, p)
julia> # Set up the $\mathbb{Z}_2$ number problem
       prob = Z2Problem(H₀)
julia> # Define a parameter range
       param1 = range(-2, 2, length=4);
       param2 = range(-0.4, 0.4, length=4)
julia> # Calculate the phase diagram (using the default Shio algorithm, single-threaded,
       # with plotting enabled and progress display)
       result = calcPhaseDiagram(prob, param1, param2; plot=true, progress=true)
(param1 = -2.0:1.3333333333333333:2.0, param2 = -0.4:0.26666666666666666:0.4, nums = [0 1 1 0; 0 1 1 0;;; 0 0 0 0; 0 0 0 0;;; 0 0 0 0; 0 0 0 0;;; 0 1 1 0; 0 1 1 0])

julia> result.param1
-2.0:1.3333333333333333:2.0

julia> result.param2
-0.4:0.26666666666666666:0.4

julia> result.nums
2×4×4 Array{Int64, 3}:
[:, :, 1] =
 0  1  1  0
 0  1  1  0

[:, :, 2] =
 0  0  0  0
 0  0  0  0

[:, :, 3] =
 0  0  0  0
 0  0  0  0

[:, :, 4] =
 0  1  1  0
 0  1  1  0
```

"""
function calcPhaseDiagram(
    prob::Z2Problem,
    param_range1::T1,
    param_range2::T2,
    alg::T3=Shio();
    parallel::T4=UseSingleThread(),
    plot::Bool=false,
    progress::Bool=false,
) where {
    T1<:AbstractVector,T2<:AbstractVector,T3<:Z2Algorithms,T4<:TopologicalNumbersParallel
}
    @unpack H, Nfill, N, rounds = prob

    dim = 2
    Ham(k) = H(k, (0.0, 0.0))
    Hs = size(Ham(zeros(2)), 1)

    if isnothing(Nfill)
        Nfill = Hs ÷ 2
    end

    p = Params(; Ham, dim, Nfill, N, rounds, Hs)

    nums = calcPhaseDiagram2D_core(
        H, param_range1, param_range2, alg, p; parallel, plot, progress
    )

    return (; param1=param_range1, param2=param_range2, nums)
end

# Second Chern number
@doc raw"""
Calculate the one-dimensional phase diagram for the second Chern number over a specified parameter range.

    calcPhaseDiagram(prob::SCProblem, param_range::T1, alg::T2=FHS2(); parallel::T3=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:SecondChernAlgorithms,T3<:TopologicalNumbersParallel}

# Description
`calcPhaseDiagram` calculates the second Chern number for a one-dimensional system described by `SCProblem` over a range of parameters (`param_range`), constructing a phase diagram. By changing the `alg` or `parallel` options, you can flexibly switch between different calculation methods and parallelization strategies. Moreover, setting `plot=true` generates a simple plot of the results upon completion of the calculation.

# Arguments
- `prob::SCProblem`: 
  A structure that contains the problem settings for second Chern number calculations, including the Hamiltonian, mesh size, number of filled bands, gap information, and other parameters. The default is a half-filled band case.

- `param_range::AbstractVector`: 
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the second Chern number.

- `alg::SecondChernAlgorithms=FHS2()`: 
  The algorithm used to compute the second Chern number. The default is `FHS2()`. You can implement and specify other algorithms as needed.

- `parallel::TopologicalNumbersParallel=UseSingleThread()`: 
  The parallelization strategy. The default is single-threaded (`UseSingleThread()`). If you wish to run in a distributed environment, you can specify it accordingly.

- `plot::Bool=false`: 
  If set to `true`, a simple plot of the calculated results is shown after the computation finishes. This is helpful for visually confirming topological phase transitions.

- `progress::Bool=false`: 
  If set to `true`, progress is displayed during loop calculations, which can be useful for large-scale computations.

# Returns
- `NamedTuple{(:param, :nums)}`: 
  - `param::AbstractVector`: A copy of the input `param_range`.
  - `nums::AbstractVector`: A vector containing the computed second Chern numbers for each parameter.

# Examples

```julia
julia> using TopologicalNumbers
julia> # Define the Hamiltonian
       H₀(k, p) = LatticeDirac(k, p)
julia> N = 10; # number of k-points in each direction
julia> # Set up the second Chern number problem
       prob = SCProblem(H₀, N)
julia> # Define a parameter range
       param = range(-4.9, 4.9, length=4)
julia> # Calculate the phase diagram (using the default FHS2 algorithm, single-threaded,
       # with plotting enabled and progress display)
       result = calcPhaseDiagram(prob, param; plot=true, progress=true)
(param = -4.9:3.2666666666666666:4.9, nums = [0.0010237313095167136, -2.0667333080974726, 2.1572606447321445, -0.0009805850180973081])

julia> result.param
-4.9:3.2666666666666666:4.9

julia> result.nums
4-element Vector{Float64}:
  0.0010237313095167136
 -2.0667333080974726
  2.1572606447321445
 -0.0009805850180973081
```


"""
function calcPhaseDiagram(
    prob::SCProblem,
    param_range::T1,
    alg::T2=FHS2();
    parallel::T3=UseSingleThread(),
    plot::Bool=false,
    progress::Bool=false,
) where {T1<:AbstractVector,T2<:SecondChernAlgorithms,T3<:TopologicalNumbersParallel}
    @unpack H, N, Nfill, RV = prob

    dim = 4
    Ham(k) = H(k, 0.0)
    Hs = size(Ham(zeros(4)), 1)

    if N isa Int
        N = (N, N, N, N)
    end
    if isnothing(Nfill)
        Nfill = Hs ÷ 2
    end

    p = Params(; Ham, dim, N, Nfill, rounds=false, returnRealValue=RV, Hs)

    nums = calcPhaseDiagram1D_core(H, param_range, alg, p; parallel, plot, progress)

    return (; param=param_range, nums)
end

@doc raw"""
Calculate the two-dimensional phase diagram for the second Chern number over a specified parameter range.

    calcPhaseDiagram(prob::SCProblem, param_range1::T1, param_range2::T2, alg::T3=FHS2(); parallel::T4=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:SecondChernAlgorithms,T4<:TopologicalNumbersParallel}

# Description
`calcPhaseDiagram` calculates the second Chern number for a two-dimensional system described by `SCProblem` over a range of parameters (`param_range`), constructing a phase diagram. By changing the `alg` or `parallel` options, you can flexibly switch between different calculation methods and parallelization strategies. Moreover, setting `plot=true` generates a simple plot of the results upon completion of the calculation.

# Arguments
- `prob::SCProblem`: 
  A structure that contains the problem settings for second Chern number calculations, including the Hamiltonian, mesh size, gap information, and other parameters.

- `param_range1::AbstractVector`: 
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the second Chern number.

- `param_range2::AbstractVector`:
  A list or range of parameters (e.g., `range(-2.0, 2.0, length=50)`) over which to scan. The Hamiltonian is evaluated at each parameter value in order to compute the second Chern number.

- `alg::SecondChernAlgorithms=FHS2()`: 
  The algorithm used to compute the second Chern number. The default is `FHS2()`. You can implement and specify other algorithms as needed.

- `parallel::TopologicalNumbersParallel=UseSingleThread()`: 
  The parallelization strategy. The default is single-threaded (`UseSingleThread()`). If you wish to run in a distributed environment, you can specify it accordingly. (Note that, at present, the second Chern number calculation does not support multithreading.)

- `plot::Bool=false`: 
  If set to `true`, a simple plot of the calculated results is shown after the computation finishes. This is helpful for visually confirming topological phase transitions.

- `progress::Bool=false`: 
  If set to `true`, progress is displayed during loop calculations, which can be useful for large-scale computations.

# Returns
- `NamedTuple{(:param, :nums)}`: 
  - `param1::AbstractVector`: A copy of the input `param_range1`.
  - `param2::AbstractVector`: A copy of the input `param_range2`.
  - `nums::AbstractArray`: An array containing the computed second Chern numbers for each parameter. The size of this array is (**number of bands** × **number of param1** × **number of param2**).

"""
function calcPhaseDiagram(
    prob::SCProblem,
    param_range1::T1,
    param_range2::T2,
    alg::T3=FHS2();
    parallel::T4=UseSingleThread(),
    plot::Bool=false,
    progress::Bool=false,
) where {
    T1<:AbstractVector,
    T2<:AbstractVector,
    T3<:SecondChernAlgorithms,
    T4<:TopologicalNumbersParallel,
}
    @unpack H, N, Nfill, RV = prob

    dim = 4
    Ham(k) = H(k, (0.0, 0.0))
    Hs = size(Ham(zeros(4)), 1)

    if N isa Int
        N = (N, N, N, N)
    end
    if isnothing(Nfill)
        Nfill = Hs ÷ 2
    end

    p = Params(; Ham, dim, N, Nfill, rounds=false, returnRealValue=RV, Hs)

    nums = calcPhaseDiagram2D_core(
        H, param_range1, param_range2, alg, p; parallel, plot, progress
    )

    return (; param1=param_range1, param2=param_range2, nums)
end
