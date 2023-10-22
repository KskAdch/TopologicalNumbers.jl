@doc raw"""
    findWeylPoint(Hamiltonian::Function)
"""
function findWeylPoint(Hamiltonian::Function)

    E(k) = eigvals(Hamiltonian(k))
    jacobE(k) = jacobian(E, [k[1], k[2], k[3]])[1]

end