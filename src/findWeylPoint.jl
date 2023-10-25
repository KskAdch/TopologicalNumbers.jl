@doc raw"""
    findWeylPoint(Hamiltonian::Function; N::Int=51)
"""
function findWeylPoint(Hamiltonian::Function; N::Int=51)

    E(k) = eigvals(Hamiltonian(k))
    jacobE(k) = jacobian(E, [k[1], k[2], k[3]])[1]

    Ene = zeros(Hs, 3)
    k = zeros(3)
    k_list = 
    n = 0
    for i in 1:N
        k[1] = 2pi*(i-1) / N
        for j in 1:N
            k[2] = 2pi*(j-1) / N
            for l in 1:N, 
                k[3] = 2pi*(l-1) / N
                Ene .= jacobE(k)
                if abs(Ene[1, 1]) < 1e-1 && abs(Ene[1, 2]) < 1e-1 && abs(Ene[1, 3]) < 1e-1 && n == 0
                    k0_list = k
                    n = 1
                elseif abs(Ene[1, 1]) < 1e-1 && abs(Ene[1, 2]) < 1e-1 && abs(Ene[1, 3]) < 1e-1
                    k0_list = [k0_list k]
                end
            end
        end
    end
end