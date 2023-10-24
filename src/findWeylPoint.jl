@doc raw"""
    findWeylPoint(Hamiltonian::Function; N::Int=51)
"""
function findWeylPoint(Hamiltonian::Function; N::Int=51)

    E(k) = eigvals(Hamiltonian(k))
    jacobE(k) = jacobian(E, [k[1], k[2], k[3]])[1]

    Ene = zeros(N, N, N, Hs)
    ene = zeros(Hs, 3)
    k = zeros(3)
    for l in 1:N
        k3 = 2pi*(l-1) / N
        for j in 1:N
            k2 = 2pi*(j-1) / N
            for i in 1:N, 
                k1 = 2pi*(i-1) / N
                k .= [k1, k2, k3]
                ene .= jacobE(k)
                Ene[i, j, l, :] .= [norm(ene[i, :], 2) for i in 1:Hs]
            end
        end
    end
end