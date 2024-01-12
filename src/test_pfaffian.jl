# From PFAPACK source https://github.com/basnijholt/pfapack
# Migration to Julia by: 


using LinearAlgebra
using SparseArrays
using Random
using Test
include("pfaffian.jl")

EPS = 1e-12


const EPS = 1e-12

@testset "Pfaffian Tests" begin
    # First test with real matrices
    A = rand(100, 100) - rand(100, 100)'
    pfa1 = pfaffian(A)
    pfa2 = pfaffian(A, method="H")
    pfa3 = pfaffian_schur(A)

    deta = det(A)

    @test abs((pfa1 - pfa2) / pfa1) < EPS
    @test abs((pfa1 - pfa3) / pfa1) < EPS
    @test abs((pfa1^2 - deta) / deta) < EPS

    # Then test with complex matrices
    A = rand(Complex{Float64}, 100, 100) - rand(Complex{Float64}, 100, 100)'
    pfa1 = pfaffian(A)
    pfa2 = pfaffian(A, method="H")

    deta = det(A)

    @test abs((pfa1 - pfa2) / pfa1) < EPS
    @test abs((pfa1^2 - deta) / deta) < EPS
end


@testset "Decomposition Tests" begin
    # First test with real matrices
    A = rand(100, 100) - rand(100, 100)'
    T, L, P = skew_LTL(A)

    @test norm(P * A * P' - L * T * L') / norm(A) < EPS

    T, Q = skew_tridiagonalize(A)

    @test norm(A - Q * T * Q') / norm(A) < EPS

    # Then test with complex matrices
    A = rand(Complex{Float64}, 100, 100) - rand(Complex{Float64}, 100, 100)'
    T, L, P = skew_LTL(A)

    @test norm(P * A * P' - L * T * L') / norm(A) < EPS

    T, Q = skew_tridiagonalize(A)

    @test norm(A - Q * T * Q') / norm(A) < EPS
end


# @pytest.mark.skipif(not with_ctypes, reason="the libs might not be installed")
# def test_ctypes():
#     for method in ("P", "H"):
#         # first real matrices
#         A = numpy.matlib.rand(100, 100)
#         A = A - A.T
#         pf_a = cpfaffian(A, uplo="L", method=method)
#         pf_a2 = cpfaffian(A, uplo="L", avoid_overflow=True, method=method)

#         np.testing.assert_almost_equal(pf_a / pf_a2, 1)
#         np.testing.assert_almost_equal(pf_a / pf.pfaffian(A), 1)

#         # then complex matrices
#         A = numpy.matlib.rand(100, 100) + 1.0j * numpy.matlib.rand(100, 100)
#         A = A - A.T
#         pf_a = cpfaffian(A, uplo="L", method=method)
#         pf_a2 = cpfaffian(A, uplo="L", avoid_overflow=True, method=method)

#         np.testing.assert_almost_equal(pf_a / pf_a2, 1)
#         np.testing.assert_almost_equal(pf_a / pf.pfaffian(A), 1)