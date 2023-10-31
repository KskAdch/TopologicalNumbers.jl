function make_k0list!(E, k0list, N, Hs, gapless)
    Evec = zeros(Hs)

    for n1 in 0:N
        for n2 in 0:N
            for n3 in 0:N
                k = 2pi/N .* [n1, n2, n3]
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
                        Evec .= E(n .* 2pi/Ni)[b:b+1]
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

function weylpoint!(Hamiltonian, k0list, Nodes, Ni, Hs, rounds)

    H(k) = Hamiltonian([k[1] - pi/Ni, k[2] - pi/Ni, k[3] - pi/Ni])

    for b in 1:Hs

        for i in 1:length(k0list[b])
            k0list[b][i][k0list[b][i] .== Ni] .= 0
        end

        unique!(k0list[b])

        j = 0
        for i in 1:length(k0list[b])

            node = calcWeylNode(H, k0list[b][i-j]; N = Ni, rounds = rounds).TopologicalNumber[b]
            if abs(node) > 1e-10
                append!(Nodes[b], node)
            else
                deleteat!(k0list[b], i-j)
                j += 1
            end
        end
    end
end

@doc raw"""
    findWeylPoint(Hamiltonian::Function; N::Int=10, gapless::Vector{Float64}=[1e-1, 1e-2, 1e-3, 1e-4], rounds::Bool=true)
"""
function findWeylPoint(Hamiltonian::Function; N::Int=10, gapless::Vector{Float64}=[1e-1, 1e-2, 1e-3, 1e-4], rounds::Bool=true)

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

    (; WeylPoint = k0list, N = Ni, Nodes)
end