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
    calcPhaseDiagram(prob::BPProblem, param_range::T1, alg::T2=BP(); parallel::T3=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:BerryPhaseAlgorithms,T3<:TopologicalNumbersParallel}
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
    calcPhaseDiagram(prob::BPProblem, param_range1::T1, param_range2::T2, alg::T3=BP(); parallel::T4=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:BerryPhaseAlgorithms,T4<:TopologicalNumbersParallel}
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
    calcPhaseDiagram(prob::FCProblem, param_range::T1, alg::T2=FHS(); parallel::T3=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:FirstChernAlgorithms,T3<:TopologicalNumbersParallel}
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
    calcPhaseDiagram(prob::FCProblem, param_range1::T1, param_range2::T2, alg::T3=FHS(); parallel::T4=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:FirstChernAlgorithms,T4<:TopologicalNumbersParallel}
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
    calcPhaseDiagram(prob::Z2Problem, param_range::T1, alg::T2=Shio(); parallel::T3=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:Z2Algorithms,T3<:TopologicalNumbersParallel}
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
    calcPhaseDiagram(prob::Z2Problem, param_range1::T1, param_range2::T2, alg::T3=Shio(); parallel::T4=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:Z2Algorithms,T4<:TopologicalNumbersParallel}
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
    calcPhaseDiagram(prob::SCProblem, param_range::T1, alg::T2=FHS2(); parallel::T3=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:SecondChernAlgorithms,T3<:TopologicalNumbersParallel}
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
    calcPhaseDiagram(prob::SCProblem, param_range1::T1, param_range2::T2, alg::T3=FHS2(); parallel::T4=UseSingleThread(), plot::Bool=false, progress::Bool=false) where {T1<:AbstractVector,T2<:AbstractVector,T3<:SecondChernAlgorithms,T4<:TopologicalNumbersParallel}
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
