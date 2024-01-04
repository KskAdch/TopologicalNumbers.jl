function psimat_cube!(n, psimat, Evec, p::Params) # wave function □
    @unpack Ham, N = p

    n100 = n .+ [1, 0, 0]
    n110 = n .+ [1, 1, 0]
    n010 = n .+ [0, 1, 0]
    n001 = n .+ [0, 0, 1]
    n101 = n .+ [1, 0, 1]
    n111 = n .+ [1, 1, 1]
    n011 = n .+ [0, 1, 1]

    n100 .= [mod(n100[i], N) for i in 1:3]
    n110 .= [mod(n110[i], N) for i in 1:3]
    n010 .= [mod(n010[i], N) for i in 1:3]
    n001 .= [mod(n001[i], N) for i in 1:3]
    n101 .= [mod(n101[i], N) for i in 1:3]
    n111 .= [mod(n111[i], N) for i in 1:3]
    n011 .= [mod(n011[i], N) for i in 1:3]

    k1 = 2pi * n / N
    k2 = 2pi * n100 / N
    k3 = 2pi * n110 / N
    k4 = 2pi * n010 / N
    k5 = 2pi * n001 / N
    k6 = 2pi * n101 / N
    k7 = 2pi * n111 / N
    k8 = 2pi * n011 / N

    eigens = eigen!(Ham(k1))
    psimat[1, :, :] .= eigens.vectors
    Evec[:] .= eigens.values

    psimat[2, :, :] .= eigen!(Ham(k2)).vectors
    psimat[3, :, :] .= eigen!(Ham(k3)).vectors
    psimat[4, :, :] .= eigen!(Ham(k4)).vectors
    psimat[5, :, :] .= eigen!(Ham(k5)).vectors
    psimat[6, :, :] .= eigen!(Ham(k6)).vectors
    psimat[7, :, :] .= eigen!(Ham(k7)).vectors
    psimat[8, :, :] .= eigen!(Ham(k8)).vectors
end

@views function Linkmat_cube!(psimat, Evec, Linkmat, p::Params)
    @unpack gapless, Hs = p

    l = 1
    while l <= Hs
        l0 = Hs - count(Evec .> (gapless + Evec[l]))

        if l == l0
            Linkmat[1, l:l0] .= dot(psimat[1, :, l:l0], psimat[2, :, l:l0])
            Linkmat[2, l:l0] .= dot(psimat[2, :, l:l0], psimat[3, :, l:l0])
            Linkmat[3, l:l0] .= dot(psimat[4, :, l:l0], psimat[3, :, l:l0])
            Linkmat[4, l:l0] .= dot(psimat[1, :, l:l0], psimat[4, :, l:l0])

            Linkmat[5, l:l0] .= dot(psimat[1, :, l:l0], psimat[5, :, l:l0])
            Linkmat[6, l:l0] .= dot(psimat[2, :, l:l0], psimat[6, :, l:l0])
            Linkmat[7, l:l0] .= dot(psimat[3, :, l:l0], psimat[7, :, l:l0])
            Linkmat[8, l:l0] .= dot(psimat[4, :, l:l0], psimat[8, :, l:l0])

            Linkmat[9, l:l0] .= dot(psimat[5, :, l:l0], psimat[6, :, l:l0])
            Linkmat[10, l:l0] .= dot(psimat[6, :, l:l0], psimat[7, :, l:l0])
            Linkmat[11, l:l0] .= dot(psimat[8, :, l:l0], psimat[7, :, l:l0])
            Linkmat[12, l:l0] .= dot(psimat[5, :, l:l0], psimat[8, :, l:l0])
        else
            Linkmat[1, l:l0] .= det(psimat[1, :, l:l0]' * psimat[2, :, l:l0])
            Linkmat[2, l:l0] .= det(psimat[2, :, l:l0]' * psimat[3, :, l:l0])
            Linkmat[3, l:l0] .= det(psimat[4, :, l:l0]' * psimat[3, :, l:l0])
            Linkmat[4, l:l0] .= det(psimat[1, :, l:l0]' * psimat[4, :, l:l0])

            Linkmat[5, l:l0] .= det(psimat[1, :, l:l0]' * psimat[5, :, l:l0])
            Linkmat[6, l:l0] .= det(psimat[2, :, l:l0]' * psimat[6, :, l:l0])
            Linkmat[7, l:l0] .= det(psimat[3, :, l:l0]' * psimat[7, :, l:l0])
            Linkmat[8, l:l0] .= det(psimat[4, :, l:l0]' * psimat[8, :, l:l0])

            Linkmat[9, l:l0] .= det(psimat[5, :, l:l0]' * psimat[6, :, l:l0])
            Linkmat[10, l:l0] .= det(psimat[6, :, l:l0]' * psimat[7, :, l:l0])
            Linkmat[11, l:l0] .= det(psimat[8, :, l:l0]' * psimat[7, :, l:l0])
            Linkmat[12, l:l0] .= det(psimat[5, :, l:l0]' * psimat[8, :, l:l0])
        end

        l = 1 + l0
    end
end

function F!(Linkmat, phi, TopologicalNumber, p::Params)
    @unpack rounds, Hs = p

    dphi = zeros(6, Hs)

    phi[1, :] = [angle(conj(Linkmat[4, l]) * Linkmat[5, l] * Linkmat[12, l] * conj(Linkmat[8, l])) for l in 1:Hs]
    dphi[1, :] = [-angle(Linkmat[4, l]) + angle(Linkmat[5, l]) + angle(Linkmat[12, l]) - angle(Linkmat[8, l]) for l in 1:Hs]

    phi[2, :] = [angle(Linkmat[2, l] * Linkmat[7, l] * conj(Linkmat[10, l]) * conj(Linkmat[6, l])) for l in 1:Hs]
    dphi[2, :] = [angle(Linkmat[2, l]) + angle(Linkmat[7, l]) - angle(Linkmat[10, l]) - angle(Linkmat[6, l]) for l in 1:Hs]

    phi[3, :] = [angle(Linkmat[1, l] * Linkmat[6, l] * conj(Linkmat[9, l]) * conj(Linkmat[5, l])) for l in 1:Hs]
    dphi[3, :] = [angle(Linkmat[1, l]) + angle(Linkmat[6, l]) - angle(Linkmat[9, l]) - angle(Linkmat[5, l]) for l in 1:Hs]

    phi[4, :] = [angle(conj(Linkmat[3, l]) * Linkmat[8, l] * Linkmat[11, l] * conj(Linkmat[7, l])) for l in 1:Hs]
    dphi[4, :] = [-angle(Linkmat[3, l]) + angle(Linkmat[8, l]) + angle(Linkmat[11, l]) - angle(Linkmat[7, l]) for l in 1:Hs]

    phi[5, :] = [angle(Linkmat[4, l] * Linkmat[3, l] * conj(Linkmat[2, l]) * conj(Linkmat[1, l])) for l in 1:Hs]
    dphi[5, :] = [angle(Linkmat[4, l]) + angle(Linkmat[3, l]) - angle(Linkmat[2, l]) - angle(Linkmat[1, l]) for l in 1:Hs]

    phi[6, :] = [angle(Linkmat[9, l] * Linkmat[10, l] * conj(Linkmat[11, l]) * conj(Linkmat[12, l])) for l in 1:Hs]
    dphi[6, :] = [angle(Linkmat[9, l]) + angle(Linkmat[10, l]) - angle(Linkmat[11, l]) - angle(Linkmat[12, l]) for l in 1:Hs]

    for j in 1:6
        if rounds == true
            phi[j, :] = [round(Int, (phi[j, i] - dphi[j, i]) / 2pi) for i in 1:Hs]
        else
            phi[j, :] .= (phi[j, :] - dphi[j, :]) ./ 2pi
        end

        TopologicalNumber .+= phi[j, :]
    end
end

@doc raw"""

 Calculate the Weyl node in the three-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    calcWeylNode(Hamiltonian::Function, n::Vector{Int64}; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

 Arguments
 - Hamiltionian::Function: the Hamiltonian matrix with three-dimensional wavenumber `k` as an argument.
 - n::Vector{Int64}: The wavenumber($2\pi\bm{n}/N$) when calculating Weyl node.
 - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.

"""
function calcWeylNode(
    Hamiltonian::Function,
    n::T;
    N::Int=51,
    gapless::Real=0.0,
    rounds::Bool=true
) where {T<:AbstractVector{Int64}}

    Hs = size(Hamiltonian(n), 1)
    p = Params(; Ham=Hamiltonian, N, gapless, rounds, Hs, dim=3)

    n .= [mod(n[i], N) for i in 1:3]

    psimat = zeros(ComplexF64, 8, Hs, Hs)
    Evec = zeros(Hs)
    Linkmat = zeros(ComplexF64, 12, Hs)
    phi = zeros(6, Hs)
    TopologicalNumber = zeros(Hs)

    psimat_cube!(n, psimat, Evec, p)
    Linkmat_cube!(psimat, Evec, Linkmat, p)
    F!(Linkmat, phi, TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
    end

    (; TopologicalNumber, n, N)
end



@doc raw"""

 Calculate the Weyl node in the three-dimensional case with reference to Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite).

    calcWeylNode(Hamiltonian::Function, n::Vector{Int64}; N::Int=51, gapless::Real=0.0, rounds::Bool=true)

 Arguments
 - Hamiltionian::Function: the Hamiltonian matrix with three-dimensional wavenumber `k` as an argument.
 - n::Vector{Int64}: The wavenumber($2\pi\bm{n}/N$) when calculating Weyl node.
 - N::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - gapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(`N`) but increasing `gapless` will increase the accuracy of the calculation.
 - rounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.

"""
function solve(prob::WNProblem,
    alg::T1=FHSlocal3();
    parallel::T2=UseSingleThread()
) where {T1<:WeylPointsAlgorithms,T2<:TopologicalNumbersParallel}
    @unpack H, n, N, gapless, rounds = prob

    Hs = size(H(n), 1)
    p = Params(; Ham=H, N, gapless, rounds, Hs, dim=3)

    n .= [mod(n[i], N) for i in 1:3]

    psimat = zeros(ComplexF64, 8, Hs, Hs)
    Evec = zeros(Hs)
    Linkmat = zeros(ComplexF64, 12, Hs)
    phi = zeros(6, Hs)
    TopologicalNumber = zeros(Hs)

    psimat_cube!(n, psimat, Evec, p)
    Linkmat_cube!(psimat, Evec, Linkmat, p)
    F!(Linkmat, phi, TopologicalNumber, p)

    if rounds == true
        TopologicalNumber = round.(Int, TopologicalNumber)
    end

    WNSolution(; TopologicalNumber, n, N)
end
