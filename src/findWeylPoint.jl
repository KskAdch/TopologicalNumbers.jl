@doc raw"""
    findWeylPoint(Hamiltonian::Function; N::Int=51)
"""
function findWeylPoint(Hamiltonian::Function; N::Int=51, itteration::Int=3)

    Hs = size(Hamiltonian(zeros(3)), 1)
    
    E(k) = eigvals(Hamiltonian(k))
    jacobE(k) = jacobian(E, [k[1], k[2], k[3]])[1]

    Ene = zeros(Hs, 3)
    k = zeros(3)
    k0_list = [Float64[]]
    for b in 2:Hs
        k0_list = append!(k0_list, [Float64[]])
    end

    for i in 1:N
        k[1] = 2pi*(i-1) / N
        for j in 1:N
            k[2] = 2pi*(j-1) / N
            for l in 1:N
                k[3] = 2pi*(l-1) / N
                Ene .= jacobE(k)
                for b in 1:Hs
                    if abs(Ene[b, 1]) < 5e-1 && abs(Ene[b, 2]) < 5e-1 && abs(Ene[b, 3]) < 5e-1
                        k0_list[b] = append!(k0_list[b], k)
                    end
                end
            end
        end
    end

    # for n in 1:(length(k0_list) ÷ 3)
    #     k1_range = range(k0_list[3*n-2]-2pi/N, k0_list[3*n-2]+2pi/N, length = N)
    #     k2_range = range(k0_list[3*n-1]-2pi/N, k0_list[3*n-1]+2pi/N, length = N)
    #     k3_range = range(k0_list[3*n]-2pi/N, k0_list[3*n]+2pi/N, length = N)
    #     for i in eachindex(k1_range)
    #         k[1] = k1_range[i]
    #         for j in eachindex(k2_range)
    #             k[2] = k1_range[j]
    #             for l in eachindex(k3_range)
    #                 k[3] = k1_range[l]
    #                 Ene .= jacobE(k)
    #                 for b in 1:Hs
    #                     if abs(Ene[b, 1]) < 1e-1 && abs(Ene[b, 2]) < 1e-1 && abs(Ene[b, 3]) < 1e-1
    #                         k0_list = append!(k, k0_list)
    #                     end
    #                 end
    #             end
    #         end
    #     end
    # end

    WeylPoint = [reshape(k0_list[1], 3, :)]
    for b in 2:Hs
        WeylPoint = append!(WeylPoint, [reshape(k0_list[b], 3, :)])
    end
    
    (; WeylPoint)
end