@testset "Pfaffian Tests" begin

    N = 100
    Ar = rand(N, N)
    Ar .= Ar .- transpose(Ar)

    Ac = rand(ComplexF64, N, N)
    Ac .= Ac .- transpose(Ac)

    @testset "householder_real" begin
        T = [Vector{Float64}, Int64, Float64]

        i = 1
        A0 = Ar[i+1:end, i]

        B1 = pf.householder_real(np.array(A0))
        B1 = [pyconvert(T[i], B1[i-1]) for i in eachindex(T)]
        B2 = TopologicalNumbers.householder_real(A0)

        for i in eachindex(T)
            @test B1[i] ≈ B2[i]
        end

        t = zeros(eltype(A0), size(A0))

        B2 = TopologicalNumbers.householder_real!(t, A0)
        @test B1[1] ≈ t

        for i in eachindex(T)[2:3]
            @test B1[i] ≈ B2[i-1]
        end

    end

    @testset "householder_complex" begin
        T = [Vector{ComplexF64}, Int64, ComplexF64]

        i = 1
        A0 = Ac[i+1:end, i]

        B1 = pf.householder_complex(np.array(A0))
        B1 = [pyconvert(T[i], B1[i-1]) for i in eachindex(T)]
        B2 = TopologicalNumbers.householder_complex(A0)

        for i in eachindex(T)
            @test B1[i] ≈ B2[i]
        end

        t = zeros(eltype(A0), size(A0))

        B2 = TopologicalNumbers.householder_complex!(t, A0)
        @test B1[1] ≈ t

        for i in eachindex(T)[2:3]
            @test B1[i] ≈ B2[i-1]
        end

    end

    @testset "skew_tridiagonalize" begin

        T = [Matrix{ComplexF64}, Matrix{ComplexF64}]
        A = Ac

        B1 = pf.skew_tridiagonalize(np.array(A))
        B1 = [pyconvert(T[i], B1[i-1]) for i in eachindex(T)]
        B2 = skew_tridiagonalize(A)

        for i in eachindex(T)
            @test B1[i] ≈ B2[i]
        end


        T = [Matrix{Float64}, Matrix{Float64}]
        A = Ar

        B1 = pf.skew_tridiagonalize(np.array(A))
        B1 = [pyconvert(T[i], B1[i-1]) for i in eachindex(T)]
        B2 = skew_tridiagonalize(A)

        for i in eachindex(T)
            @test B1[i] ≈ B2[i]
        end

    end

    @testset "skew_LTL" begin

        T = [Matrix{ComplexF64}, Matrix{ComplexF64}, Matrix{Float64}]
        A = Ac

        B1 = pf.skew_LTL(np.array(A))
        B1 = [[pyconvert(T[i], B1[i-1]) for i in eachindex(T)[1:end-1]]..., pyconvert(T[end], B1[2].toarray())]
        B2 = skew_LTL(A)
        B2 = (B2[1], B2[2], Matrix(B2[3]))

        for i in eachindex(T)
            @test B1[i] ≈ B2[i]
        end


        T = [Matrix{Float64}, Matrix{Float64}, Matrix{Float64}]
        A = Ar

        B1 = pf.skew_LTL(np.array(A))
        B1 = [[pyconvert(T[i], B1[i-1]) for i in eachindex(T)[1:end-1]]..., pyconvert(T[end], B1[2].toarray())]
        B2 = skew_LTL(A)
        B2[3] = Matrix(B2[3])

        for i in eachindex(T)
            @test B1[i] ≈ B2[i]
        end

    end

    @testset "pfaffian_LTL" begin

        T = ComplexF64
        A = Ac

        B1 = pyconvert(T, pf.pfaffian_LTL(np.array(A)))
        B2 = pfaffian_LTL(A)

        @test B1 ≈ B2


        T = Float64
        A = Ar

        B1 = pyconvert(T, pf.pfaffian_LTL(np.array(A)))
        B2 = pfaffian_LTL(A)

        @test B1 ≈ B2

    end

    @testset "pfaffian_householder" begin

        T = ComplexF64
        A = Ac

        B1 = pyconvert(T, pf.pfaffian_householder(np.array(A)))
        B2 = pfaffian_householder(A)

        @test B1 ≈ B2


        T = Float64
        A = Ar

        B1 = pyconvert(T, pf.pfaffian_householder(np.array(A)))
        B2 = pfaffian_householder(A)

        @test B1 ≈ B2

    end

    @testset "pfaffian_schur" begin

        T = Float64
        A = Ar

        B1 = pyconvert(T, pf.pfaffian_schur(np.array(A)))
        B2 = pfaffian_schur(A)

        @test B1 ≈ B2

    end

    @testset "pfaffian_P" begin

        T = ComplexF64
        A = Ac

        B1 = pyconvert(T, pf.pfaffian(np.array(A), "P"))
        B2 = pfaffian(A, method="P")

        @test B1 ≈ B2

        B3 = det(A)
        @test B1^2 ≈ B3
        @test B2^2 ≈ B3


        T = Float64
        A = Ar

        B1 = pyconvert(T, pf.pfaffian(np.array(A), "P"))
        B2 = pfaffian(A, method="P")

        @test B1 ≈ B2

        B3 = det(A)
        @test B1^2 ≈ B3
        @test B2^2 ≈ B3

    end

    @testset "pfaffian_H" begin

        T = ComplexF64
        A = Ac

        B1 = pyconvert(T, pf.pfaffian(np.array(A), "H"))
        B2 = pfaffian(A, method="H")

        @test B1 ≈ B2

        B3 = det(A)
        @test B1^2 ≈ B3
        @test B2^2 ≈ B3


        T = Float64
        A = Ar

        B1 = pyconvert(T, pf.pfaffian(np.array(A), "H"))
        B2 = pfaffian(A, method="H")

        @test B1 ≈ B2

        B3 = det(A)
        @test B1^2 ≈ B3
        @test B2^2 ≈ B3

    end


    # From PFAPACK source https://github.com/basnijholt/pfapack
    # Migration from PFAPACK(Python) 
    EPS = 1e-12

    @testset "Pfaffian Tests" begin
        # First test with real matrices
        A = Ar

        pfa1 = pfaffian(A)
        pfa2 = pfaffian(A, method="H")
        pfa3 = pfaffian_schur(A)

        deta = det(A)

        @test abs((pfa1 - pfa2) / pfa1) < EPS
        @test abs((pfa1 - pfa3) / pfa1) < EPS
        @test abs((pfa1^2 - deta) / deta) < EPS

        # Then test with complex matrices
        A = Ac

        pfa1 = pfaffian(A)
        pfa2 = pfaffian(A, method="H")

        deta = det(A)

        @test abs((pfa1 - pfa2) / pfa1) < EPS
        @test abs((pfa1^2 - deta) / deta) < EPS
    end


    @testset "Decomposition Tests" begin
        # First test with real matrices
        A = Ar

        T, L, P = skew_LTL(A)

        @test norm(P * A * transpose(P) - L * T * transpose(L)) / norm(A) < EPS

        T, Q = skew_tridiagonalize(A)

        @test norm(A - Q * T * transpose(Q)) / norm(A) < EPS

        # Then test with complex matrices
        A = Ac

        T, L, P = skew_LTL(A)

        @test norm(P * A * transpose(P) - L * T * transpose(L)) / norm(A) < EPS

        T, Q = skew_tridiagonalize(A)

        @test norm(A - Q * T * transpose(Q)) / norm(A) < EPS
    end

end