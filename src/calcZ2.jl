function psi_j!(j, v::TemporalZ2, p::Params) # wave function
    @unpack Ham, N = p
    for i in 1:N
        v.k[1] = (i - 1) * 2pi / N
        v.k[2] = (j - 1) * 2pi / N
        v.psi_1[i, :, :] .= eigen!(Ham(v.k)).vectors
    end
end

# @views function Link!(v::TemporalZ2)
@views function Link(v::TemporalZ2, p::Params)
    @unpack Nfill = p
    # for l in 1:v.Hshalf
    #     v.link10[l] = det(v.psi00[:, 2l-1:2l]' * v.psi10[:, 2l-1:2l])
    #     v.link01[l] = det(v.psi00[:, 2l-1:2l]' * v.psi01[:, 2l-1:2l])
    # end
    [
        det(v.psi00[:, 1:Nfill]' * v.psi10[:, 1:Nfill]) det(v.psi00[:, 1:Nfill]' * v.psi01[:, 1:Nfill])
        det(v.psi00[:, Nfill+1:end]' * v.psi10[:, Nfill+1:end]) det(v.psi00[:, Nfill+1:end]' * v.psi01[:, Nfill+1:end])
    ]
end

@views function U!(j, v::TemporalZ2, p::Params) # link variable
    @unpack Nfill, N = p

    if j != 1
        v.Link1 .= v.Link2
    end

    if j != N
        if j == 1

            psi_j!(1, v, p)
            v.psi_0 .= v.psi_1
            v.psi_N .= v.psi_1
            psi_j!(2, v, p)

            for i in 1:N

                v.psi00[:, :] = v.psi_0[i, :, :]
                if i == N
                    v.psi10[:, :] = v.psi_0[1, :, :]
                    v.psi01[:, :] = v.psi_1[N, :, :]
                else
                    v.psi10[:, :] = v.psi_0[i+1, :, :]
                    v.psi01[:, :] = v.psi_1[i, :, :]
                end

                # Link!(v)

                # v.Link1[:, 1, i] .= v.link10
                # v.Link1[:, 2, i] .= v.link01
                v.Link1[:, :, i] .= Link(v, p)
                v.LinkN1 .= v.Link1
            end

            v.w00 .= v.psi_0[1, :, :]' * v.T * conj(v.psi_0[1, :, :])
            if isapprox(round.(v.w00[1:Nfill, Nfill+1:end], digits = 5), zero(v.w00[1:Nfill, Nfill+1:end])) == false
                v.T .= kron(v.sy, v.s0)
            end
            v.w00 .= v.psi_0[1, :, :]' * v.T * conj(v.psi_0[1, :, :])
            v.wp0 .= v.psi_0[v.Nhalf, :, :]' * v.T * conj(v.psi_0[v.Nhalf, :, :])
        end

        v.psi_0 .= v.psi_1

        if j != N - 1
            psi_j!(j + 2, v, p)
        end

        if j == v.Nhalf - 1
            v.w0p .= v.psi_0[1, :, :]' * v.T * conj(v.psi_0[1, :, :])
            v.wpp .= v.psi_0[v.Nhalf, :, :]' * v.T * conj(v.psi_0[v.Nhalf, :, :])
        end

        for i in 1:N

            v.psi00[:, :] = v.psi_0[i, :, :]
            if i == N && j == N - 1
                v.psi10[:, :] = v.psi_0[1, :, :]
                v.psi01[:, :] = v.psi_N[N, :, :]
            elseif i == N
                v.psi10[:, :] = v.psi_0[1, :, :]
                v.psi01[:, :] = v.psi_1[N, :, :]
            elseif j == N - 1
                v.psi10[:, :] = v.psi_0[i+1, :, :]
                v.psi01[:, :] = v.psi_N[i, :, :]
            else
                v.psi10[:, :] = v.psi_0[i+1, :, :]
                v.psi01[:, :] = v.psi_1[i, :, :]
            end

            # Link!(v)

            # v.Link2[:, 1, i] .= v.link10
            # v.Link2[:, 2, i] .= v.link01
            v.Link2[:, :, i] .= Link(v, p)
        end
    end
end

# @views function F!(phi, Px0, Pxp, i, j, v::TemporalZ2, p::Params) # lattice field strength
@views function F!(i, j, v::TemporalZ2, p::Params) # lattice field strength
    @unpack N = p

    if j == 1 && i < v.Nhalf
        # for l in 1:v.Hshalf
        for l in 1:2
            # Px0[l] -= angle(v.Link1[l, 1, i])
            v.Px0[l] -= angle(v.Link1[l, 1, i])
        end
    elseif j == v.Nhalf && i < v.Nhalf
        # for l in 1:v.Hshalf
        for l in 1:2
            # Pxp[l] -= angle(v.Link1[l, 1, i])
            v.Pxp[l] -= angle(v.Link1[l, 1, i])
        end
    end

    if i != N && j != N
        # for l in 1:v.Hshalf
        for l in 1:2
        #     phi[l] = angle(v.Link1[l, 1, i] * v.Link1[l, 2, i+1] * conj(v.Link2[l, 1, i]) * conj(v.Link1[l, 2, i]))
            v.phi[l] = angle(v.Link1[l, 1, i] * v.Link1[l, 2, i+1] * conj(v.Link2[l, 1, i]) * conj(v.Link1[l, 2, i]))
        end
        # phi = angle(v.Link1[1, i] * v.Link1[2, i+1] * conj(v.Link2[1, i]) * conj(v.Link1[2, i]))
    elseif i == N && j != N
        # for l in 1:v.Hshalf
        for l in 1:2
        #     phi[l] = angle(v.Link1[l, 1, N] * v.Link1[l, 2, 1] * conj(v.Link2[l, 1, N]) * conj(v.Link1[l, 2, N]))
            v.phi[l] = angle(v.Link1[l, 1, N] * v.Link1[l, 2, 1] * conj(v.Link2[l, 1, N]) * conj(v.Link1[l, 2, N]))
        end
        # phi = angle(v.Link1[1, N] * v.Link1[2, 1] * conj(v.Link2[1, N]) * conj(v.Link1[2, N]))
    elseif i != N && j == N
        # for l in 1:v.Hshalf
        for l in 1:2
        #     phi[l] = angle(v.Link1[l, 1, i] * v.Link1[l, 2, i+1] * conj(v.LinkN1[l, 1, i]) * conj(v.Link1[l, 2, i]))
            v.phi[l] = angle(v.Link1[l, 1, i] * v.Link1[l, 2, i+1] * conj(v.LinkN1[l, 1, i]) * conj(v.Link1[l, 2, i]))
        end
        # phi = angle(v.Link1[1, i] * v.Link1[2, i+1] * conj(v.LinkN1[1, i]) * conj(v.Link1[2, i]))
    elseif i == N && j == N
        # for l in 1:v.Hshalf
        for l in 1:2
        #     phi[l] = angle(v.Link1[l, 1, N] * v.Link1[l, 2, 1] * conj(v.LinkN1[l, 1, N]) * conj(v.Link1[l, 2, N]))
            v.phi[l] = angle(v.Link1[l, 1, N] * v.Link1[l, 2, 1] * conj(v.LinkN1[l, 1, N]) * conj(v.Link1[l, 2, N]))
        end
        # phi = angle(v.Link1[1, N] * v.Link1[2, 1] * conj(v.LinkN1[1, N]) * conj(v.Link1[2, N]))
    end
end

@doc raw"""
"""
# @views function Z2Phase!(TopologicalNumber, p::Params) # chern number
@views function Z2Phase!(v::TemporalZ2, p::Params) # chern number
    @unpack Nfill, N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1

    if size(v.num) == (2,)
        for j in 1:Nhalf
            U!(j, v, p)
            for i in 1:N
                # F!(phi, Px0, Pxp, i, j, v, p)
                F!(i, j, v, p)
                if j < Nhalf
                    # TN[:] .+= phi[:]
                    v.num[:] .+= v.phi[:]
                end
            end
        end
    elseif size(v.num) == (2, 2)
        for j in 1:N
            U!(j, v, p)
            for i in 1:N
                # F!(phi, Px0, Pxp, i, j, v, p)
                F!(i, j, v, p)
                if j < Nhalf
                    # TN[:, 1] .+= phi[:]
                    v.num[:, 1] .+= v.phi[:]
                else
                    # TN[:, 2] .+= phi[:]
                    v.num[:, 2] .+= v.phi[:]
                end
            end
        end
    end

    v.Px0[1] += angle(pyconvert(ComplexF64, pf.pfaffian(np.matrix(v.w00[1:Nfill, 1:Nfill])) / pf.pfaffian(np.matrix(v.wp0[1:Nfill, 1:Nfill]))))
    v.Pxp[1] += angle(pyconvert(ComplexF64, pf.pfaffian(np.matrix(v.w0p[1:Nfill, 1:Nfill])) / pf.pfaffian(np.matrix(v.wpp[1:Nfill, 1:Nfill]))))

    v.Px0[2] += angle(pyconvert(ComplexF64, pf.pfaffian(np.matrix(v.w00[Nfill+1:end, Nfill+1:end])) / pf.pfaffian(np.matrix(v.wp0[Nfill+1:end, Nfill+1:end]))))
    v.Pxp[2] += angle(pyconvert(ComplexF64, pf.pfaffian(np.matrix(v.w0p[Nfill+1:end, Nfill+1:end])) / pf.pfaffian(np.matrix(v.wpp[Nfill+1:end, Nfill+1:end]))))

    # for l in 1:Hshalf
    #     Px0[l] += angle((v.w00[2l-1, 2l]) / (v.wp0[2l-1, 2l]))
    #     Pxp[l] += angle((v.w0p[2l-1, 2l]) / (v.wpp[2l-1, 2l]))

    #     if TN[l] - 2Px0[l] + 2Pxp[l] !== NaN
    #         TopologicalNumber[l] = 1 - abs(1 - rem(abs(TN[l] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
    #     end
    # end
    for l in 1:2
        v.num[l, :] = 1 .- abs.(1 .- rem.(abs.(v.num[l, :] .- 2v.Px0[l] .+ 2v.Pxp[l]) ./ 2pi, 2))
    end
    # v.num .= 1 - abs.(1 - rem.(abs.(v.num .- 2v.Px0 .+ 2v.Pxp) ./ 2pi, 2))
end

# @doc raw"""
# """
# # @views function Z2Phase!(TopologicalNumber, TRTopologicalNumber, p::Params) # chern number
# @views function TRZ2Phase(p::Params) # chern number
#     @unpack Nfill, N, Hs, rounds = p
#     Nhalf = N ÷ 2 + 1
#     Hshalf = Hs ÷ 2

#     k = zeros(2)

#     # Link1 = zeros(ComplexF64, Hshalf, 2, N)
#     # Link2 = zeros(ComplexF64, Hshalf, 2, N)
#     # LinkN1 = zeros(ComplexF64, Hshalf, 2, N)
#     Link1 = zeros(ComplexF64, 2, N)
#     Link2 = zeros(ComplexF64, 2, N)
#     LinkN1 = zeros(ComplexF64, 2, N)
#     # link10 = zeros(ComplexF64, Hshalf)
#     # link01 = zeros(ComplexF64, Hshalf)

#     psi_0 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
#     psi_1 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
#     psi_N = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)

#     psi00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     psi10 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     psi01 = zeros(ComplexF64, 2Hshalf, 2Hshalf)

#     # phi = zeros(Hshalf)
#     # Px0 = zeros(Hshalf)
#     # Pxp = zeros(Hshalf)
#     Px = zeros(2)

#     s0 = Matrix{ComplexF64}(I, Hshalf, Hshalf)
#     sy = [0 -1; 1 0] # imaginary???
#     T = kron(s0, sy)

#     w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
#     wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

#     # v = TemporalZ2(k, T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, link10, link01, psi_0, psi_1, psi_N, psi00, psi10, psi01, Nhalf, Hshalf)
#     v = TemporalZ2(k, T, w00, w0p, wp0, wpp, Link1, Link2, LinkN1, psi_0, psi_1, psi_N, psi00, psi10, psi01, Nhalf, Hshalf)

#     # TN = zeros(Hshalf, 2)
#     TN = zeros(2)

#     for j in 1:N
#         U!(j, v, p)
#         for i in 1:N
#             # F!(phi, Px0, Pxp, i, j, v, p)
#             Pol!(i, j, v)
#             if j < Nhalf
#                 # TN[:, 1] .+= phi[:]
#                 TN[1] += F(i, j, v, p)
#             else
#                 # TN[:, 2] .+= phi[:]
#                 TN[2] += F(i, j, v, p)
#             end
#         end
#     end

#     Px[1] += angle(pyconvert(ComplexF64, pf.pfaffian(np.matrix(v.w00[1:Nfill, 1:Nfill])) / pf.pfaffian(np.matrix(v.wp0[1:Nfill, 1:Nfill]))))
#     Px[2] += angle(pyconvert(ComplexF64, pf.pfaffian(np.matrix(v.w0p[1:Nfill, 1:Nfill])) / pf.pfaffian(np.matrix(v.wpp[1:Nfill, 1:Nfill]))))

#     # for l in 1:Hshalf
#     #     Px0[l] += angle((v.w00[2l-1, 2l]) / (v.wp0[2l-1, 2l]))
#     #     Pxp[l] += angle((v.w0p[2l-1, 2l]) / (v.wpp[2l-1, 2l]))

#     #     if TN[l] - 2Px0[l] + 2Pxp[l] !== NaN
#     #         TopologicalNumber[l] = 1 - abs(1 - rem(abs(TN[l, 1] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
#     #         TRTopologicalNumber[l] = 1 - abs(1 - rem(abs(TN[l, 2] - 2Px0[l] + 2Pxp[l]) / 2pi, 2))
#     #     end
#     # end
#     1 - abs.(1 - rem.(abs.(TN .- 2Px[1] .+ 2Px[2]) / 2pi, 2))
# end

function setTemporalZ2(p::Params)
    @unpack Nfill, N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1
    Hshalf = Hs ÷ 2

    k = zeros(2)

    s0 = Matrix{ComplexF64}(I, Hshalf, Hshalf)
    sy = [0 -1; 1 0] # imaginary???
    T = kron(s0, sy)

    w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    Px0 = zeros(2)
    Pxp = zeros(2)
    phi = zeros(2)

    Link1 = zeros(ComplexF64, 2, 2, N)
    Link2 = zeros(ComplexF64, 2, 2, N)
    LinkN1 = zeros(ComplexF64, 2, 2, N)
    
    psi_0 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    psi_1 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    psi_N = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    
    psi00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    psi10 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    psi01 = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    num = zeros(2)
    
    TemporalZ2(num, k, sy, s0, T, w00, w0p, wp0, wpp, Px0, Pxp, phi, Link1, Link2, LinkN1, psi_0, psi_1, psi_N, psi00, psi10, psi01, Nhalf, Hshalf)
end

function setTemporalZ2TR(p::Params)
    @unpack Nfill, N, Hs, rounds = p
    Nhalf = N ÷ 2 + 1
    Hshalf = Hs ÷ 2

    k = zeros(2)

    s0 = Matrix{ComplexF64}(I, Hshalf, Hshalf)
    sy = [0 -1; 1 0] # imaginary???
    T = kron(s0, sy)

    w00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    w0p = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wp0 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    wpp = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    Px0 = zeros(2)
    Pxp = zeros(2)
    phi = zeros(2)

    Link1 = zeros(ComplexF64, 2, 2, N)
    Link2 = zeros(ComplexF64, 2, 2, N)
    LinkN1 = zeros(ComplexF64, 2, 2, N)
    
    psi_0 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    psi_1 = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    psi_N = zeros(ComplexF64, N, 2Hshalf, 2Hshalf)
    
    psi00 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    psi10 = zeros(ComplexF64, 2Hshalf, 2Hshalf)
    psi01 = zeros(ComplexF64, 2Hshalf, 2Hshalf)

    num = zeros(2, 2)
    
    TemporalZ2(num, k, sy, s0, T, w00, w0p, wp0, wpp, Px0, Pxp, phi, Link1, Link2, LinkN1, psi_0, psi_1, psi_N, psi00, psi10, psi01, Nhalf, Hshalf)
end

@doc raw"""

 Calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case with reference to Shiozaki method [Fukui2007Quantum,Shiozaki2023discrete](@cite).

    calcZ2(Hamiltonian::Function; N::Int=50, rounds::Bool=true, TR::Bool=false)

 Arguments
 - `Hamiltonian::Function` is a matrix with one-dimensional wavenumber `k` as an argument.
 - `N::Int` is the number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - `rounds::Bool` is an option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.

# Definition
 The $\mathbb{Z}_{2}$ number of the $2n$th (and $2n-1$th) band $\nu_{n}$ is defined by
```math
\nu_{n}=F_{n}-\left(P_{n}(0)-P_{n}(\pi)\right)
```
 $F_{n}$ is the Berry flux of the $n$th band in the $\mathrm{BZ}'$. The range $\mathrm{BZ}'$ is $\bm{k}\in[0,2\pi]\times[0,\pi]$ half of BZ(Brillouin Zone).
```math
F_{n}=\frac{1}{2\pi}\sum_{\bm{k}\in\mathrm{BZ}'}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,1}^{*}(\bm{k})\right]\right]
```
 $P_{n}(k_{2})$ is the time-reversal polarization at wavenumber $k_{2}$.
```math
P_{n}(k_{2})=\frac{1}{2\pi}\frac{\mathrm{Pf}[\omega(0,k_{2})]}{\mathrm{Pf}[\omega(\pi,k_{2})]}\sum_{k_{1}=0}^{\pi-e_{1}}U_{n,1}(\bm{k})
```
 $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $2n$th (and $2n-1$th) band. $\omega(\bm{k})$ is the unitary matrix given by
```math
\omega(\bm{k})=\bra{\Psi(-\bm{k})}T\ket{\Psi(\bm{k})}
```
 $T$ is the time-reversal operator.
"""
# function calcZ2(Hamiltonian::Function; N::Int=50, rounds::Bool=true, TR::Bool=false)
function calcZ2(Hamiltonian::Function; Nfill::T1=nothing, N::Int=50, rounds::Bool=true, TR::Bool=false) where {T1<:Union{Int,Nothing}}

    Hs = size(Hamiltonian(zeros(2)), 1)
    Hshalf = Hs ÷ 2
    if isodd(N)
        throw(ArgumentError("N should be an even number"))
    else
        if isnothing(Nfill)
            Nfill = Hshalf
        elseif isodd(Nfill)
            throw(ArgumentError("Nfill should be an even number"))
        end
    end
    p = Params(; Ham=Hamiltonian, Nfill, N, Hs, gapless=0.0, rounds, dim=2)

    if TR == false

        v = setTemporalZ2(p)
        Z2Phase!(v, p)

        if rounds == true
            TopologicalNumber = round.(Int, v.num)
            Total = mod(sum(TopologicalNumber), 2)
        elseif rounds == false
            TopologicalNumber = v.num
            Total = sum(TopologicalNumber)
            Total = mod(Total, 2*sign(Total-2))
        end

        (; TopologicalNumber, Total)
    elseif TR == true
        
        v = setTemporalZ2TR(p)
        Z2Phase!(v, p)

        if rounds == true
            TopologicalNumber = round.(Int, v.num[:, 1])
            TRTopologicalNumber = round.(Int, v.num[:, 2])

            Total = mod(sum(v.num[:, 1]), 2)
        elseif rounds == false
            TopologicalNumber = v.num[:, 1]
            TRTopologicalNumber = v.num[:, 2]
            Total = sum(TopologicalNumber)
            Total = mod(Total, 2*sign(Total-2))
        end

        (; TopologicalNumber, TRTopologicalNumber, Total)
    end

    # TopologicalNumber = zeros(Hshalf)
    # if TR == false
    #     # Z2Phase!(TopologicalNumber, p)
    #     TopologicalNumber = Z2Phase(p)

    #     if rounds == true
    #         TopologicalNumber = round.(Int, TopologicalNumber)
    #         # Total = rem(sum(TopologicalNumber), 2)
    #     else
    #         # Total = abs(sum(TopologicalNumber))
    #         # Total = abs(rem(Total, 2))
    #     end

    #     # (; TopologicalNumber, Total)
    #     (; TopologicalNumber)
    # else
    #     # TRTopologicalNumber = zeros(Hshalf)
    #     # Z2Phase!(TopologicalNumber, TRTopologicalNumber, p)
    #     TopologicalNumber, TRTopologicalNumber = TRZ2Phase(p)

    #     if rounds == true
    #         TopologicalNumber = round.(Int, TopologicalNumber)
    #         TRTopologicalNumber = round.(Int, TRTopologicalNumber)
    #         # Total = rem(sum(TopologicalNumber), 2)
    #     else
    #         # Total = abs(sum(TopologicalNumber))
    #         # Total = rem(1 - abs(1 - Total), 2)
    #     end

    #     # (; TopologicalNumber, TRTopologicalNumber, Total)
    #     (; TopologicalNumber, TRTopologicalNumber)
    # end
end


@doc raw"""

 Calculate the $\mathbb{Z}_2$ numbers in the two-dimensional case with reference to Shiozaki method [Fukui2007Quantum,Shiozaki2023discrete](@cite).

    solve(prob::Z2Problem, alg::T1=Shio(); parallel::T2=UseSingleThread()) where {T1<:Z2Algorithms,T2<:TopologicalNumbersParallel}

 Arguments
 - `Hamiltonian::Function` is a matrix with one-dimensional wavenumber `k` as an argument.
 - `N::Int` is the number of meshes when discretizing the Brillouin Zone. It is preferable for `N` to be an odd number to increase the accuracy of the calculation.
 - `rounds::Bool` is an option to round the value of the topological number to an integer value. The topological number returns a value of type `Int` when `true`, and a value of type `Float` when `false`.

# Definition
 The $\mathbb{Z}_{2}$ number of the $2n$th (and $2n-1$th) band $\nu_{n}$ is defined by
```math
\nu_{n}=F_{n}-\left(P_{n}(0)-P_{n}(\pi)\right)
```
 $F_{n}$ is the Berry flux of the $n$th band in the $\mathrm{BZ}'$. The range $\mathrm{BZ}'$ is $\bm{k}\in[0,2\pi]\times[0,\pi]$ half of BZ(Brillouin Zone).
```math
F_{n}=\frac{1}{2\pi}\sum_{\bm{k}\in\mathrm{BZ}'}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,1}^{*}(\bm{k})\right]\right]
```
 $P_{n}(k_{2})$ is the time-reversal polarization at wavenumber $k_{2}$.
```math
P_{n}(k_{2})=\frac{1}{2\pi}\frac{\mathrm{Pf}[\omega(0,k_{2})]}{\mathrm{Pf}[\omega(\pi,k_{2})]}\sum_{k_{1}=0}^{\pi-e_{1}}U_{n,1}(\bm{k})
```
 $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $2n$th (and $2n-1$th) band. $\omega(\bm{k})$ is the unitary matrix given by
```math
\omega(\bm{k})=\bra{\Psi(-\bm{k})}T\ket{\Psi(\bm{k})}
```
 $T$ is the time-reversal operator.
"""
function solve(
    prob::Z2Problem,
    alg::T1=Shio();
    parallel::T2=UseSingleThread()
) where {T1<:Z2Algorithms,T2<:TopologicalNumbersParallel}
    @unpack H, N, rounds, TR = prob

    Hs = size(H(zeros(2)), 1)
    Hshalf = Hs ÷ 2
    p = Params(; Ham=H, N, Hs, rounds, dim=2)

    TopologicalNumber = zeros(Hshalf)
    if TR == false
        Z2Phase!(TopologicalNumber, p)

        if rounds == true
            TopologicalNumber = round.(Int, TopologicalNumber)
            Total = rem(sum(TopologicalNumber), 2)
        else
            Total = abs(sum(TopologicalNumber))
            Total = abs(rem(Total, 2))
        end

        Z2Solution(; TopologicalNumber, Total)
    else
        TRTopologicalNumber = zeros(Hshalf)
        Z2Phase!(TopologicalNumber, TRTopologicalNumber, p)

        if rounds == true
            TopologicalNumber = round.(Int, TopologicalNumber)
            TRTopologicalNumber = round.(Int, TRTopologicalNumber)
            Total = rem(sum(TopologicalNumber), 2)
        else
            Total = abs(sum(TopologicalNumber))
            Total = rem(1 - abs(1 - Total), 2)
        end

        Z2Solution(; TopologicalNumber, TRTopologicalNumber, Total)
    end
end