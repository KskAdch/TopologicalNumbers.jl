function make_k0list!(E, k0list, N, Hs, gapless)
    Evec = zeros(Hs)

    for n1 in 0:N
        for n2 in 0:N
            for n3 in 0:N
                k = 2pi / N .* [n1, n2, n3]
                Evec .= E(k)

                for b in 1:Hs-1
                    dE = Evec[b+1] .- Evec[b]
                    if dE < gapless
                        append!(k0list[b], [[n1, n2, n3]])
                        append!(k0list[b+1], [[n1, n2, n3]])
                    end
                    unique!(k0list[b])
                end
            end
        end
    end
end

function update_k0list!(E, k0list, N, Ni, Hs, gapless)
    Evec = zeros(2)

    for b in 1:Hs-1

        if b == 1
            k0list_b = copy(k0list[b])
            empty!(k0list[b])
        end

        k0list_b1 = copy(k0list[b+1])
        empty!(k0list[b+1])

        for i in 1:length(k0list_b)

            for n1 in (-N+1):(N-1)

                if N * k0list_b[i][1] .+ n1 .< 0 || N * k0list_b[i][1] .+ n1 .> Ni
                    @goto n1_loop
                end

                for n2 in (-N+1):(N-1)

                    if N * k0list_b[i][2] .+ n2 .< 0 || N * k0list_b[i][2] .+ n2 .> Ni
                        @goto n2_loop
                    end

                    for n3 in (-N+1):(N-1)

                        if N * k0list_b[i][3] .+ n3 .< 0 || N * k0list_b[i][3] .+ n3 .> Ni
                            @goto n3_loop
                        end

                        n = N .* k0list_b[i] .+ [n1, n2, n3]
                        Evec .= E(n .* 2pi / Ni)[b:b+1]
                        dE = Evec[2] .- Evec[1]

                        if dE < gapless
                            append!(k0list[b], [n])
                            append!(k0list[b+1], [n])
                        end

                        @label n3_loop
                    end
                    @label n2_loop
                end
                @label n1_loop
            end
        end
        unique!(k0list[b])
        k0list_b = copy(k0list_b1)
    end
end

@doc raw"""
"""
function weylpoint!(Hamiltonian, k0list, Nodes, Ni, Hs, rounds)

    H(k) = Hamiltonian([k[1] - pi / Ni, k[2] - pi / Ni, k[3] - pi / Ni])

    for b in 1:Hs

        for i in 1:length(k0list[b])
            k0list[b][i][k0list[b][i].==Ni] .= 0
        end

        unique!(k0list[b])

        j = 0
        for i in 1:length(k0list[b])

            node = calcWeylNode(H, k0list[b][i-j]; N=Ni, rounds=rounds).TopologicalNumber[b]
            if abs(node) > 1e-10
                append!(Nodes[b], node)
            else
                deleteat!(k0list[b], i - j)
                j += 1
            end
        end
    end
end

@doc raw"""
    findWeylPoint(Hamiltonian::Function; N::Int=10, gapless::T=[1e-1, 1e-2, 1e-3, 1e-4], rounds::Bool=true) where {T<:AbstractVector{Float64}}

 Arguments
 - Hamiltionian::Function: The Hamiltonian matrix with three-dimensional wavenumber `k` as an argument.
 - N::Int=10: The number of meshes when discretizing the Brillouin Zone. The $n$th iteration divides the Brillouin zone into $1/N^n$.
 - gapless<:AbstractVector{Float64}: The threshold that determines the state to be degenerate. The $n$th iteration adopts the threshold value of the $n$th value of the vector. The number of iterations can be varied by the length of the vector.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.
   
"""
function findWeylPoint(
    Hamiltonian::Function;
    N::Int=10,
    gapless::T=[1e-1, 1e-2, 1e-3, 1e-4],
    rounds::Bool=true
) where {T<:AbstractVector{Float64}}

    Hs = size(Hamiltonian(zeros(3)), 1)
    k0list = [Vector{Int64}[] for i in 1:Hs]

    if rounds == true
        Nodes = [Int64[] for i in 1:Hs]
    else
        Nodes = [Float64[] for i in :Hs]
    end

    E(k) = eigvals(Hamiltonian(k))

    make_k0list!(E, k0list, N, Hs, gapless[1])

    Ni = N
    for iter in 1:(length(gapless)-1)
        Ni *= N
        update_k0list!(E, k0list, N, Ni, Hs, gapless[iter])
    end

    weylpoint!(Hamiltonian, k0list, Nodes, Ni, Hs, rounds)

    (; WeylPoint=k0list, N=Ni, Nodes)
end


@doc raw"""
    solve(prob::WPProblem, alg::T1=Evar(); parallel::T2=UseSingleThread()) where {T1<:WeylPointsAlgorithms,T2<:TopologicalNumbersParallel}

 Arguments
 - Hamiltionian::Function: The Hamiltonian matrix with three-dimensional wavenumber `k` as an argument.
 - N::Int=10: The number of meshes when discretizing the Brillouin Zone. The $n$th iteration divides the Brillouin zone into $1/N^n$.
 - gapless<:AbstractVector{Float64}: The threshold that determines the state to be degenerate. The $n$th iteration adopts the threshold value of the $n$th value of the vector. The number of iterations can be varied by the length of the vector.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.
   
"""
function solve(prob::WPProblem,
    alg::T1=Evar();
    parallel::T2=UseSingleThread()
) where {T1<:WeylPointsAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack H, N, gapless, rounds = prob

    Hs = size(H(zeros(3)), 1)
    k0list = [Vector{Int64}[] for i in 1:Hs]

    if rounds == true
        Nodes = [Int64[] for i in 1:Hs]
    else
        Nodes = [Float64[] for i in :Hs]
    end

    E(k) = eigvals(H(k))

    make_k0list!(E, k0list, N, Hs, gapless[1])

    Ni = N
    for iter in 1:(length(gapless)-1)
        Ni *= N
        update_k0list!(E, k0list, N, Ni, Hs, gapless[iter])
    end

    weylpoint!(H, k0list, Nodes, Ni, Hs, rounds)

    WPSolution(; WeylPoint=k0list, N=Ni, Nodes)
end
