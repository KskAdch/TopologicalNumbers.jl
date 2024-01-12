# From PFAPACK source https://github.com/basnijholt/pfapack
# Make Julia package:

"""A package for computing Pfaffians"""

using LinearAlgebra
using SparseArrays

"""
    (v, tau, alpha) = householder_real(x)

Compute a Householder transformation such that
(1-tau v v^T) x = alpha e_1
where x and v a real vectors, tau is 0 or 2, and
alpha a real number (e_1 is the first unit vector)
"""
function householder_real(x)

    assert(length(x) > 0)

    sigma = dot(x[2:end], x[2:end])

    if sigma == 0
        return zeros(length(x)), 0, x[1]
    else
        norm_x = sqrt(x[1]^2 + sigma)
        v = copy(x)

        if x[1] <= 0
            v[1] -= norm_x
            alpha = norm_x
        else
            v[1] += norm_x
            alpha = -norm_x
        end

        v /= norm(v)
        return v, 2, alpha
    end
end


"""
    (v, tau, alpha) = householder_real(x)

Compute a Householder transformation such that
(1-tau v v^T) x = alpha e_1
where x and v a complex vectors, tau is 0 or 2, and
alpha a complex number (e_1 is the first unit vector)
"""
function householder_complex(x)
    assert(length(x) > 0)

    sigma = dot(conj.(x[2:end]), x[2:end])

    if sigma == 0
        return zeros(eltype(x), length(x)), 0, x[1]
    else
        norm_x = sqrt(conj(x[1]) * x[1] + sigma)
        v = copy(x)

        phase = exp(im * atan2(imag(x[1]), real(x[1])))

        v[1] += phase * norm_x
        v /= norm(v)
    end

    return (v, 2, -phase * norm_x)
end

"""
    T, Q = skew_tridiagonalize(A, overwrite_a, calc_q=True)

or

    T = skew_tridiagonalize(A, overwrite_a, calc_q=False)

Bring a real or complex skew-symmetric matrix (A=-A^T) into
tridiagonal form T (with zero diagonal) with a orthogonal
(real case) or unitary (complex case) matrix U such that
A = Q T Q^T
(Note that Q^T and *not* Q^dagger also in the complex case)

A is overwritten if overwrite_a=True (default: False), and
Q only calculated if calc_q=True (default: True)
"""
function skew_tridiagonalize(A; overwrite_a=false, calc_q=true)
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A + A')) < 1e-14

    # Check if we have a complex data type
    if eltype(A) <: Complex
        householder = householder_complex
    elseif !(eltype(A) <: Number)
        throw(TypeError(skew_tridiagonalize, "skew_tridiagonalize can only work on numeric input", Number))
    else
        householder = householder_real
    end

    if !overwrite_a
        A = copy(A)
    end

    if calc_q
        Q = I(size(A, 1))
    end

    for i in 1:size(A, 1)-2
        # Find a Householder vector to eliminate the i-th column
        v, tau, alpha = householder(A[i+1:end, i])
        A[i+1, i] = alpha
        A[i, i+1] = -alpha
        A[i+2:end, i] = zeros(eltype(A), size(A, 1) - i - 1)
        A[i, i+2:end] = zeros(eltype(A), size(A, 1) - i - 1)

        # Update the matrix block A(i+1:N,i+1:N)
        w = tau * (A[i+1:end, i+1:end] * conj(v))
        A[i+1:end, i+1:end] += v * w' - w * v'

        if calc_q
            # Accumulate the individual Householder reflections
            y = tau * (Q[:, i+1:end] * v)
            Q[:, i+1:end] -= y * v'
        end
    end

    if calc_q
        return A, Q
    else
        return A
    end
end

"""
    T, L, P = skew_LTL(A, overwrite_a, calc_q=True)

Bring a real or complex skew-symmetric matrix (A=-A^T) into
tridiagonal form T (with zero diagonal) with a lower unit
triangular matrix L such that
P A P^T= L T L^T

A is overwritten if overwrite_a=True (default: False),
L and P only calculated if calc_L=True or calc_P=True,
respectively (default: True).
"""
function skew_LTL(A; overwrite_a=false, calc_L=true, calc_P=true)
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A + A')) < 1e-14

    n = size(A, 1)
    if !overwrite_a
        A = copy(A)
    end

    if calc_L
        L = I(n)
    end

    if calc_P
        Pv = collect(1:n)
    end

    for k in 1:n-2
        # First, find the largest entry in A[k+1:end,k] and permute it to A[k+1,k]
        kp = k + 1 + argmax(abs.(A[k+1:end, k]))

        # Check if we need to pivot
        if kp != k + 1
            # Interchange rows k+1 and kp
            A[[k + 1, kp], k:end] = A[[kp, k + 1], k:end]

            # Then interchange columns k+1 and kp
            A[k:end, [k + 1, kp]] = A[k:end, [kp, k + 1]]

            if calc_L
                # Permute L accordingly
                L[[k + 1, kp], 1:k+1] = L[[kp, k + 1], 1:k+1]
            end

            if calc_P
                # Accumulate the permutation matrix
                Pv[k+1], Pv[kp] = Pv[kp], Pv[k+1]
            end
        end

        # Now form the Gauss vector
        if A[k+1, k] != 0.0
            tau = A[k+2:end, k] / A[k+1, k]

            # Clear eliminated row and column
            A[k+2:end, k] .= 0.0
            A[k, k+2:end] .= 0.0

            # Update the matrix block A[k+2:end,k+2:end]
            A[k+2:end, k+2:end] .+= tau * A[k+2:end, k+1]'
            A[k+2:end, k+2:end] .-= A[k+2:end, k+1] * tau'

            if calc_L
                L[k+2:end, k+1] = tau
            end
        end
    end

    if calc_P
        # Form the permutation matrix as a sparse matrix
        P = spzeros(n, n)
        for i in 1:n
            P[i, Pv[i]] = 1
        end
    end

    if calc_L
        if calc_P
            return A, L, P
        else
            return A, L
        end
    else
        if calc_P
            return A, P
        else
            return A
        end
    end
end

"""
    pfaffian(A, overwrite_a=False, method='P')

Compute the Pfaffian of a real or complex skew-symmetric
matrix A (A=-A^T). If overwrite_a=True, the matrix A
is overwritten in the process. This function uses
either the Parlett-Reid algorithm (method='P', default),
or the Householder tridiagonalization (method='H')
"""
function pfaffian(A; overwrite_a=false, method="P")
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A + A')) < 1e-14
    # Check that the method variable is appropriately set
    @assert method in ["P", "H"]

    if method == "P"
        return pfaffian_LTL(A, overwrite_a=overwrite_a)
    else
        return pfaffian_householder(A, overwrite_a=overwrite_a)
    end
end


"""
    pfaffian_LTL(A, overwrite_a=False)

Compute the Pfaffian of a real or complex skew-symmetric
matrix A (A=-A^T). If overwrite_a=True, the matrix A
is overwritten in the process. This function uses
the Parlett-Reid algorithm.
"""
function pfaffian_LTL(A; overwrite_a=false)
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A + A')) < 1e-14

    n, m = size(A)
    if eltype(A) != Complex{Float64}
        A = convert(Array{Float64,2}, A)
    end

    # Quick return if possible
    if n % 2 == 1
        return 0.0
    end

    if !overwrite_a
        A = copy(A)
    end

    pfaffian_val = 1.0

    for k in 1:2:n-1
        # First, find the largest entry in A[k+1:end, k] and permute it to A[k+1, k]
        kp = k + 1 + argmax(abs.(A[k+1:end, k]))

        # Check if we need to pivot
        if kp != k + 1
            # Interchange rows k+1 and kp
            A[[k + 1, kp], k:end] = A[[kp, k + 1], k:end]

            # Then interchange columns k+1 and kp
            A[k:end, [k + 1, kp]] = A[k:end, [kp, k + 1]]

            # Every interchange corresponds to a "-" in det(P)
            pfaffian_val *= -1
        end

        # Now form the Gauss vector
        if A[k+1, k] != 0.0
            tau = A[k, k+2:end] / A[k, k+1]

            pfaffian_val *= A[k, k+1]

            if k + 2 <= n
                # Update the matrix block A[k+2:end, k+2:end]
                A[k+2:end, k+2:end] .+= tau * A[k+2:end, k+1]'
                A[k+2:end, k+2:end] .-= A[k+2:end, k+1] * tau'
            end
        else
            # If we encounter a zero on the super/subdiagonal, the Pfaffian is 0
            return 0.0
        end
    end

    return pfaffian_val
end

"""
    pfaffian(A, overwrite_a=False)

Compute the Pfaffian of a real or complex skew-symmetric
matrix A (A=-A^T). If overwrite_a=True, the matrix A
is overwritten in the process. This function uses the
Householder tridiagonalization.

Note that the function pfaffian_schur() can also be used in the
real case. That function does not make use of the skew-symmetry
and is only slightly slower than pfaffian_householder().
"""
function pfaffian_householder(A; overwrite_a=false)
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A + A')) < 1e-14

    n = size(A, 1)

    if eltype(A) != Complex{Float64}
        A = convert(Array{Float64,2}, A)
    end

    # Quick return if possible
    if n % 2 == 1
        return 0.0
    end

    # Determine the appropriate Householder transformation function
    if eltype(A) <: Complex
        householder = householder_complex
    else
        householder = householder_real
    end

    if !overwrite_a
        A = copy(A)
    end

    pfaffian_val = 1.0

    for i in 1:n-2
        # Find a Householder vector to eliminate the i-th column
        v, tau, alpha = householder(A[i+1:end, i])
        A[i+1, i] = alpha
        A[i, i+1] = -alpha
        A[i+2:end, i] .= 0
        A[i, i+2:end] .= 0

        # Update the matrix block A[i+1:end, i+1:end]
        w = tau * (A[i+1:end, i+1:end] * conj(v))
        A[i+1:end, i+1:end] .+= v * w' - w * v'

        if tau != 0
            pfaffian_val *= 1 - tau
        end
        if i % 2 == 0
            pfaffian_val *= -alpha
        end
    end

    pfaffian_val *= A[n-2, n-1]

    return pfaffian_val
end

"""
Calculate Pfaffian of a real antisymmetric matrix using
the Schur decomposition. (Hessenberg would in principle be faster,
but scipy-0.8 messed up the performance for scipy.linalg.hessenberg()).

This function does not make use of the skew-symmetry of the matrix A,
but uses a LAPACK routine that is coded in FORTRAN and hence faster
than python. As a consequence, pfaffian_schur is only slightly slower
than pfaffian().
"""
function pfaffian_schur(A; overwrite_a=false)
    @assert eltype(A) <: Real
    @assert size(A, 1) == size(A, 2) > 0
    @assert maximum(abs.(A + A')) < 1e-14

    # Quick return if possible
    if size(A, 1) % 2 == 1
        return 0.0
    end

    if !overwrite_a
        A = copy(A)
    end

    T, Z = schur(A)
    l = diag(T, 1)  # Get the superdiagonal of T
    return prod(l[1:2:end]) * det(Z)
end