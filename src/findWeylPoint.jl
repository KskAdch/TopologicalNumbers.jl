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
Calculate the Weyl points in the three-dimensional case using energy variational method.

    solve(prob::WPProblem, alg::T1=Evar(); parallel::T2=UseSingleThread()) where {T1<:WeylPointsAlgorithms,T2<:TopologicalNumbersParallel}

# Arguments
- `prob::WPProblem`: The WPProblem struct that contains the Hamiltonian matrix function in the wave number space and other parameters.
- `alg::T1=Evar()`: The algorithm to use for calculating the Weyl points. Default is `Evar` algorithm.
- `parallel::T2=UseSingleThread()`: The parallelization strategy to use. Default is to use a single thread.

# Returns
- `WPSolution`: A struct that contains the calculated Weyl points.

# Examples

```julia
julia> function H₀(k, p) # Weyl
    k1, k2, k3 = k
    t1, t2, t3, m, k0 = p

    h0 = 0
    hx = 2t1*(cos(k1) - cos(k0)) + m*(2 - cos(k2) - cos(k3))
    hy = 2t2*sin(k2)
    hz = 2t3*sin(k3)

    s0 = [1 0; 0 1]
    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    sz = [1 0; 0 -1]

    h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
end
julia> p0 = (1, 1, 1, 2, 2pi*2/5);
julia> H(k) = H₀(k, p0);
julia> prob = WPProblem(H)
julia> result = solve(prob)
WPSolution{Vector{Vector{Vector{Int64}}}, Int64, Vector{Vector{Int64}}}([[[4000, 0, 0], [6000, 0, 0]], [[4000, 0, 0], [6000, 0, 0]]], 10000, [[1, -1], [-1, 1]])
julia> 2pi*result.WeylPoint[1] / result.N .- pi*[ones(3), ones(3)]
2-element Vector{Vector{Float64}}:
 [-0.6283185307179586, -3.141592653589793, -3.141592653589793]
 [0.6283185307179586, -3.141592653589793, -3.141592653589793]
```

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
