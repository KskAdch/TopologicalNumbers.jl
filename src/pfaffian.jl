# From PFAPACK source https://github.com/basnijholt/pfapack
# Migration from PFAPACK(Python)

"""A package for computing Pfaffians"""

@doc raw"""
    (v, tau, alpha) = householder_real(x)

Compute a Householder transformation such that
(1-tau v v^T) x = alpha e_1
where x and v a real vectors, tau is 0 or 2, and
alpha a real number (e_1 is the first unit vector)
"""
function householder_real(x)

    @assert length(x) > 0

    sigma = sum(abs2, @view(x[2:end]))

    if sigma == 0
        return zeros(length(x)), 0, x[1]
    else
        norm_x = sqrt(abs2(x[1]) + sigma)
        v = copy(x)
        # depending on whether x[0] is positive or negatvie
        # choose the sign
        if x[1] <= 0
            v[1] -= norm_x
            alpha = norm_x
        else
            v[1] += norm_x
            alpha = -norm_x
        end

        v ./= norm(v)
        return v, 2, alpha
    end
end

@doc raw"""
    (tau, alpha) = householder_real!(v, x)

Compute a Householder transformation such that
(1-tau v v^T) x = alpha e_1
where x and v a real vectors, tau is 0 or 2, and
alpha a real number (e_1 is the first unit vector)
"""
function householder_real!(v, x)

    @assert length(x) > 0

    sigma = sum(abs2, @view(x[2:end]))

    if sigma == 0
        v .= zero(eltype(x))
        return 0, x[1]
    else
        v .= x
        norm_x = sqrt(abs2(x[1]) + sigma)
        # depending on whether x[0] is positive or negatvie
        # choose the sign
        if x[1] <= 0
            v[1] -= norm_x
            alpha = norm_x
        else
            v[1] += norm_x
            alpha = -norm_x
        end

        v ./= norm(v)
        return 2, alpha
    end
end


@doc raw"""
    (v, tau, alpha) = householder_complex(x)

Compute a Householder transformation such that
(1-tau v v^T) x = alpha e_1
where x and v a complex vectors, tau is 0 or 2, and
alpha a complex number (e_1 is the first unit vector)
"""
function householder_complex(x)
    @assert length(x) > 0

    sigma = sum(abs2, @view(x[2:end]))

    if sigma == 0
        return zeros(eltype(x), length(x)), 0, x[1]
    else
        norm_x = sqrt(abs2(x[1]) + sigma)
        v = copy(x)

        phase = exp(im * atan(imag(x[1]), real(x[1])))

        v[1] += phase * norm_x
        v ./= norm(v)
    end

    return v, 2, -phase * norm_x
end

@doc raw"""
    (tau, alpha) = householder_complex!(v, x)

Compute a Householder transformation such that
(1-tau v v^T) x = alpha e_1
where x and v a complex vectors, tau is 0 or 2, and
alpha a complex number (e_1 is the first unit vector)
"""
function householder_complex!(v, x)
    @assert length(x) > 0

    sigma = sum(abs2, @view(x[2:end]))

    if sigma == 0
        v .= zero(eltype(x))
        return 0, x[1]
    else
        v .= x
        norm_x = sqrt(abs2(x[1]) + sigma)

        phase = exp(im * atan(imag(x[1]), real(x[1])))

        v[1] += phase * norm_x
        v ./= norm(v)
    end

    return 2, -phase * norm_x
end

@doc raw"""
    T, Q = skew_tridiagonalize(A::AbstractMatrix{T}; overwrite_a=false, calc_q=true) where {T<:Number}

or

    T = skew_tridiagonalize(A::AbstractMatrix{T}; overwrite_a=false, calc_q=false) where {T<:Number}

Bring a real or complex skew-symmetric matrix (A=-A^T) into
tridiagonal form T (with zero diagonal) with a orthogonal
(real case) or unitary (complex case) matrix U such that
A = Q T Q^T
(Note that Q^T and *not* Q^dagger also in the complex case)

A is overwritten if overwrite_a=true (default: false), and
Q only calculated if calc_q=true (default: true)
"""
function skew_tridiagonalize(A::AbstractMatrix{T}; overwrite_a=false, calc_q=true) where {T<:Number}
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A + transpose(A))) < 1e-14

    # Check if we have a complex data type
    if T <: Complex
        householder = householder_complex
    else
        householder = householder_real
    end

    if !overwrite_a
        A = copy(A)
    end

    if calc_q
        Q = Matrix{T}(I, size(A))
    end

    @inbounds for i in 1:size(A, 1)-2
        # Find a Householder vector to eliminate the i-th column
        v, tau, alpha = householder(@view(A[i+1:end, i]))
        A[i+1, i] = alpha
        A[i, i+1] = -alpha
        A[i+2:end, i] .= zero(T)
        A[i, i+2:end] .= zero(T)

        # Update the matrix block A(i+1:N,i+1:N)
        w = tau .* (@view(A[i+1:end, i+1:end]) * conj(v))
        A[i+1:end, i+1:end] .+= v * transpose(w) .- w * transpose(v)

        if calc_q
            # Accumulate the individual Householder reflections
            # Accumulate them in the form P_1*P_2*..., which is
            # (..*P_2*P_1)^dagger
            y = tau .* (@view(Q[:, i+1:end]) * v)
            Q[:, i+1:end] .-= y * v'
        end
    end

    if calc_q
        return A, Q
    else
        return A
    end
end

@doc raw"""
    T, L, P = skew_LTL(A::AbstractMatrix{T}; overwrite_a=false, calc_L=true, calc_P=true) where {T<:Number}

Bring a real or complex skew-symmetric matrix (A=-A^T) into
tridiagonal form T (with zero diagonal) with a lower unit
triangular matrix L such that
P A P^T= L T L^T

A is overwritten if overwrite_a=true (default: false),
L and P only calculated if calc_L=true or calc_P=true,
respectively (default: true).
"""
function skew_LTL(A::AbstractMatrix{T}; overwrite_a=false, calc_L=true, calc_P=true) where {T<:Number}
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A + transpose(A))) < 1e-14

    n = size(A, 1)
    if !overwrite_a
        A = copy(A)
    end

    if calc_L
        L = Matrix{eltype(A)}(I, n, n)
    end

    if calc_P
        Pv = collect(1:n)
    end

    for k in 1:n-2
        # First, find the largest entry in A[k+1:end,k] and permute it to A[k+1,k]
        kp = k + findmax(abs.(@view(A[k+1:end, k])))[2]

        # Check if we need to pivot
        if kp != k + 1
            # Interchange rows k+1 and kp
            temp = copy(@view(A[k+1, k:end]))
            A[k+1, k:end] .= @view(A[kp, k:end])
            A[kp, k:end] .= temp

            # Then interchange columns k+1 and kp
            temp = copy(A[k:end, k+1])
            A[k:end, k+1] .= A[k:end, kp]
            A[k:end, kp] .= temp

            if calc_L
                # Permute L accordingly
                temp = copy(L[k+1, 1:k])
                L[k+1, 1:k] = L[kp, 1:k]
                L[kp, 1:k] = temp
            end

            if calc_P
                # Accumulate the permutation matrix
                Pv[k+1], Pv[kp] = Pv[kp], Pv[k+1]
            end
        end

        # Now form the Gauss vector
        if A[k+1, k] != 0.0
            tau = A[k+2:end, k] ./ A[k+1, k]

            # Clear eliminated row and column
            A[k+2:end, k] .= 0.0
            A[k, k+2:end] .= 0.0

            # Update the matrix block A[k+2:end,k+2:end]
            A[k+2:end, k+2:end] .+= tau * transpose(A[k+2:end, k+1])
            A[k+2:end, k+2:end] .-= A[k+2:end, k+1] * transpose(tau)

            if calc_L
                L[k+2:end, k+1] .= tau
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

@doc raw"""
    pfaffian(A::AbstractMatrix{T}, overwrite_a=false, method='P') where {T<:Number}

Compute the Pfaffian of a real or complex skew-symmetric
matrix A (A=-A^T). If overwrite_a=true, the matrix A
is overwritten in the process. This function uses
either the Parlett-Reid algorithm (method='P', default),
or the Householder tridiagonalization (method='H')
"""
function pfaffian(A::AbstractMatrix{T}; overwrite_a=false, method="P") where {T<:Number}
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A .+ transpose(A))) < 1e-14
    # Check that the method variable is appropriately set
    @assert method in ("P", "H")

    if method == "P"
        return pfaffian_LTL(A, overwrite_a=overwrite_a)
    else
        return pfaffian_householder(A, overwrite_a=overwrite_a)
    end
end


@doc raw"""
    pfaffian_LTL(A::AbstractMatrix{T}; overwrite_a=false) where {T<:Number}

Compute the Pfaffian of a real or complex skew-symmetric
matrix A (A=-A^T). If overwrite_a=true, the matrix A
is overwritten in the process. This function uses
the Parlett-Reid algorithm.
"""
function pfaffian_LTL(A::AbstractMatrix{T}; overwrite_a=false) where {T<:Number}
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A .+ transpose(A))) < 1e-14

    n, m = size(A)
    # if !(eltype(A) <: Complex)
    #     A = convert(Array{Float64,2}, A)
    # end

    # Quick return if possible
    if n % 2 == 1
        return zero(T)
    end

    if !overwrite_a
        A = copy(A)
    end

    tau = Array{T}(undef, n - 2)

    pfaffian_val = one(T)

    @inbounds for k in 1:2:n-1
        tau0 = @view tau[k:end]

        # First, find the largest entry in A[k+1:end, k] and permute it to A[k+1, k]
        @views kp = k + findmax(abs.(A[k+1:end, k]))[2]

        # Check if we need to pivot
        if kp != k + 1
            # Interchange rows k+1 and kp
            @inbounds @simd for l in k:n
                t = A[k+1, l]
                A[k+1, l] = A[kp, l]
                A[kp, l] = t
            end

            # Then interchange columns k+1 and kp
            @inbounds @simd for l in k:n
                t = A[l, k+1]
                A[l, k+1] = A[l, kp]
                A[l, kp] = t
            end

            # Every interchange corresponds to a "-" in det(P)
            pfaffian_val *= -1
        end

        # Now form the Gauss vector
        @inbounds if A[k+1, k] != zero(T)
            @inbounds @views tau0 .= A[k, k+2:end] ./ A[k, k+1]

            pfaffian_val *= @inbounds A[k, k+1]

            if k + 2 <= n
                # Update the matrix block A[k+2:end, k+2:end]
                @inbounds for l1 in eachindex(tau0)
                    @simd for l2 in eachindex(tau0)
                        @fastmath A[k+1+l2, k+1+l1] += tau0[l2] * A[k+1+l1, k+1] - tau0[l1] * A[k+1+l2, k+1]
                    end
                end
            end
        else
            # If we encounter a zero on the super/subdiagonal, the Pfaffian is 0
            return zero(T)
        end
    end

    return pfaffian_val
end

@doc raw"""
    pfaffian_householder(A::AbstractMatrix{T}; overwrite_a=false) where {T<:Complex}

Compute the Pfaffian of a real or complex skew-symmetric
matrix A (A=-A^T). If overwrite_a=true, the matrix A
is overwritten in the process. This function uses the
Householder tridiagonalization.

Note that the function pfaffian_schur() can also be used in the
real case. That function does not make use of the skew-symmetry
and is only slightly slower than pfaffian_householder().
"""
function pfaffian_householder(A::AbstractMatrix{T}; overwrite_a=false) where {T<:Complex}
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A .+ transpose(A))) < 1e-14

    n, m = size(A)

    # if !(T <: Complex)
    #     A = convert(Array{Float64,2}, A)
    # end

    # Quick return if possible
    if n % 2 == 1
        return zero(T)
    end

    # Determine the appropriate Householder transformation function
    householder! = householder_complex!

    if !overwrite_a
        A = copy(A)
    end

    # Z = Array{T}(undef, n - 1, m - 1)
    v = Array{T}(undef, n - 1)
    w = copy(v)

    pfaffian_val = one(T)

    @inbounds for i in 1:(n-2)
        # Z0 = @view Z[i:end, i:end]
        v0 = @view v[i:end]
        w0 = @view w[i:end]

        # Find a Householder vector to eliminate the i-th column
        @views tau, alpha = householder!(v0, A[i+1:end, i])
        A[i+1, i] = alpha
        A[i, i+1] = -alpha
        @views A[i+2:end, i] .= zero(T)
        @views A[i, i+2:end] .= zero(T)

        # Update the matrix block A[i+1:end, i+1:end]
        # A0 = @view A[i+1:end, i+1:end]
        # @views mul!(w0, A0, conj.(v0))
        @views mul!(w0, A[i+1:end, i+1:end], conj.(v0))
        w0 .*= tau
        # mul!(Z0, v0, transpose(w0))
        # @views A0 .+= Z0
        # @inbounds @simd for j in eachindex(A0)
        #     A0[j] = A0[j] + Z0[j]
        # end
        # mul!(Z0, w0, transpose(v0))
        # @views A0 .-= Z0
        # @inbounds @simd for j in eachindex(A0)
        #     A0[j] = A0[j] - Z0[j]
        # end
        @inbounds for j in eachindex(w0)
            @simd for k in eachindex(v0)
                # A0[k, j] += v0[k] * w0[j] - v0[j] * w0[k]
                @fastmath A[i+k, i+j] += v0[k] * w0[j] - v0[j] * w0[k]
            end
        end

        if tau != 0
            pfaffian_val *= 1 - tau
        end
        if (i - 1) % 2 == 0
            pfaffian_val *= -alpha
        end
    end

    pfaffian_val *= @inbounds A[n-1, n]

    return pfaffian_val
end

@doc raw"""
    pfaffian_householder(A::AbstractMatrix{T}; overwrite_a=false) where {T<:Real}

Compute the Pfaffian of a real or complex skew-symmetric
matrix A (A=-A^T). If overwrite_a=true, the matrix A
is overwritten in the process. This function uses the
Householder tridiagonalization.

Note that the function pfaffian_schur() can also be used in the
real case. That function does not make use of the skew-symmetry
and is only slightly slower than pfaffian_householder().
"""
function pfaffian_householder(A::AbstractMatrix{T}; overwrite_a=false) where {T<:Real}
    # Check if matrix is square
    @assert size(A, 1) == size(A, 2) > 0
    # Check if it's skew-symmetric
    @assert maximum(abs.(A .+ transpose(A))) < 1e-14

    n, m = size(A)

    # if !(T <: Complex)
    #     A = convert(Array{Float64,2}, A)
    # end

    # Quick return if possible
    if n % 2 == 1
        return zero(T)
    end

    # Determine the appropriate Householder transformation function
    householder! = householder_real!

    if !overwrite_a
        A = copy(A)
    end

    # Z = Array{T}(undef, n - 1, m - 1)
    v = Array{T}(undef, n - 1)
    w = copy(v)

    pfaffian_val = one(T)

    @inbounds for i in 1:(n-2)
        # Z0 = @view Z[i:end, i:end]
        v0 = @view v[i:end]
        w0 = @view w[i:end]

        # Find a Householder vector to eliminate the i-th column
        @views tau, alpha = householder!(v0, A[i+1:end, i])
        A[i+1, i] = alpha
        A[i, i+1] = -alpha
        @views A[i+2:end, i] .= zero(T)
        @views A[i, i+2:end] .= zero(T)

        # Update the matrix block A[i+1:end, i+1:end]
        # A0 = @view A[i+1:end, i+1:end]
        # @views mul!(w0, A0, conj.(v0))
        @views mul!(w0, A[i+1:end, i+1:end], v0)
        w0 .*= tau
        # mul!(Z0, v0, transpose(w0))
        # @views A0 .+= Z0
        # @inbounds @simd for j in eachindex(A0)
        #     A0[j] = A0[j] + Z0[j]
        # end
        # mul!(Z0, w0, transpose(v0))
        # @views A0 .-= Z0
        # @inbounds @simd for j in eachindex(A0)
        #     A0[j] = A0[j] - Z0[j]
        # end
        @inbounds for j in eachindex(w0)
            @simd for k in eachindex(v0)
                # A0[k, j] += v0[k] * w0[j] - v0[j] * w0[k]
                @fastmath A[i+k, i+j] += v0[k] * w0[j] - v0[j] * w0[k]
            end
        end

        if tau != 0
            pfaffian_val *= 1 - tau
        end
        if (i - 1) % 2 == 0
            pfaffian_val *= -alpha
        end
    end

    pfaffian_val *= @inbounds A[n-1, n]

    return pfaffian_val
end

@doc raw"""
    pfaffian_schur(A::AbstractMatrix{T}; overwrite_a=false) where {T<:Real}

Calculate Pfaffian of a real antisymmetric matrix using
the Schur decomposition. (Hessenberg would in principle be faster,
but scipy-0.8 messed up the performance for scipy.linalg.hessenberg()).

This function does not make use of the skew-symmetry of the matrix A,
but uses a LAPACK routine that is coded in FORTRAN and hence faster
than python. As a consequence, pfaffian_schur is only slightly slower
than pfaffian().
"""
function pfaffian_schur(A::AbstractMatrix{TY}; overwrite_a=false) where {TY<:Real}
    # @assert eltype(A) <: Real
    @assert size(A, 1) == size(A, 2) > 0
    @assert maximum(abs.(A .+ transpose(A))) < 1e-14

    # Quick return if possible
    if size(A, 1) % 2 == 1
        return 0.0
    end

    if !overwrite_a
        A = copy(A)
    end

    T, Z = schur(A)
    l = diag(T, 1)  # Get the superdiagonal of T
    return prod(@view(l[1:2:end])) * det(Z)
end