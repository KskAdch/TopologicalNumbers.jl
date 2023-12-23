using LinearAlgebra
using PythonPlot
using LaTeXStrings
using StaticArrays
using Accessors
using Parameters

# Kronecker delta function
KroneckerDelta(i, j) = i == j ? 1 : 0

# Hamiltonian
function H(kx, ky, kz, kw, p)
    return sin(kx) * p.G2 + sin(ky) * p.G3 + sin(kz) * p.G4 + sin(kw) * p.G5 + (p.m0 + (cos(kx) + cos(ky) + cos(kz) + cos(kw))) * p.G1
end

# Compute eigenstates
function compute_eigenstates(p)
    vecs = Array{Complex{Float64},6}(undef, p.pts + 1, p.pts + 1, p.pts + 1, p.pts + 1, 4, p.Nocc)
    for kx in range(-π, stop=π, length=p.pts + 1)
        for ky in range(-π, stop=π, length=p.pts + 1)
            for kz in range(-π, stop=π, length=p.pts + 1)
                for kw in range(-π, stop=π, length=p.pts + 1)
                    Hk = H(kx, ky, kz, kw, p)
                    vals, vecs_k = eigen(Hermitian(Hk))
                    sorted_indices = sortperm(vals)
                    vecs_k = vecs_k[:, sorted_indices[1:p.Nocc]]
                    vecs[round(Int, (kx + π) / (2π / p.pts))+1, round(Int, (ky + π) / (2π / p.pts))+1, round(Int, (kz + π) / (2π / p.pts))+1, round(Int, (kw + π) / (2π / p.pts))+1, :, :] = vecs_k
                end
            end
        end
    end
    return vecs
end

# Projector
function projector(vecs, a, b, c, d)
    y = vecs[a, b, c, d, :, :]
    x = y * y'
    return x
end

# Compute all projectors
function compute_all_projectors(vecs, p)
    allProjs = Array{Complex{Float64},6}(undef, p.pts + 1, p.pts + 1, p.pts + 1, p.pts + 1, 4, 4)
    for a in 1:(p.pts+1)
        for b in 1:(p.pts+1)
            for c in 1:(p.pts+1)
                for d in 1:(p.pts+1)
                    allProjs[a, b, c, d, :, :] = projector(vecs, a, b, c, d)
                end
            end
        end
    end
    return allProjs
end

# Loop function (F)
function proj_loop(allProjs, a, b, c, d, e, f, vecs, p)
    P1 = allProjs[a, b, c, d, :, :]
    P2 = allProjs[a+KroneckerDelta(e, 1),
        b+KroneckerDelta(e, 2),
        c+KroneckerDelta(e, 3),
        d+KroneckerDelta(e, 4), :, :]
    P3 = allProjs[a+KroneckerDelta(e, 1)+KroneckerDelta(f, 1),
        b+KroneckerDelta(e, 2)+KroneckerDelta(f, 2),
        c+KroneckerDelta(e, 3)+KroneckerDelta(f, 3),
        d+KroneckerDelta(e, 4)+KroneckerDelta(f, 4), :, :]
    P4 = allProjs[a+KroneckerDelta(f, 1),
        b+KroneckerDelta(f, 2),
        c+KroneckerDelta(f, 3),
        d+KroneckerDelta(f, 4), :, :]

    v = vecs[a, b, c, d, :, :]
    y = v' * P4 * P3 * P2 * v
    z = y * p.Nocc / tr(y) # U
    x = log(z)
    return x
end

# Compute second Chern number
function compute_second_chern_number(allProjs, vecs, p)
    FwedgeF = zeros(Complex{Float64}, p.pts, p.pts, p.pts, p.pts)
    for a in 1:p.pts, b in 1:p.pts, c in 1:p.pts, d in 1:p.pts
        FwedgeF[a, b, c, d] = tr(proj_loop(allProjs, a, b, c, d, 1, 2, vecs, p) * proj_loop(allProjs, a, b, c, d, 3, 4, vecs, p)) +
                              tr(proj_loop(allProjs, a, b, c, d, 2, 3, vecs, p) * proj_loop(allProjs, a, b, c, d, 1, 4, vecs, p)) +
                              tr(proj_loop(allProjs, a, b, c, d, 3, 1, vecs, p) * proj_loop(allProjs, a, b, c, d, 2, 4, vecs, p))
    end
    C2 = sum(FwedgeF) / (4 * π^2)
    return C2
end

function makeFigure(mList, ChernList)

    fig = figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.plot(mList, ChernList, "o--")
    ax.set_xlabel(L"m")
    ax.set_ylabel("Second Chern number")
    ax.grid(true)
    tight_layout()
    savefig("secondChern_GPT.png")
    plotshow()
end

# Main computation
function main()

    # Parameters
    mList = range(-5.0, 5.0, length=10)
    # mList = [-3.0]
    ChernList = zeros(ComplexF64, length(mList))
    Nocc, m0, pts = 2, -3, 30

    # Pauli matrices and Dirac matrices
    sx, sy, sz = [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]
    s0 = Matrix{Complex{Float64}}(I, 2, 2)
    G1, G2, G3, G4, G5 = kron(sx, s0), kron(sy, s0), kron(sz, sx), kron(sz, sy), kron(sz, sz)

    p = (; Nocc, m0, pts, sx, sy, sz, s0, G1, G2, G3, G4, G5)

    for i in eachindex(mList)
        @reset p.m0 = mList[i]

        vecs = compute_eigenstates(p)
        allProjs = compute_all_projectors(vecs, p)
        C2 = compute_second_chern_number(allProjs, vecs, p)
        println("m: ", p.m0, ", Second Chern number: ", C2)
        ChernList[i] = C2
    end

    # makeFigure(mList, ChernList)


end

@time main()