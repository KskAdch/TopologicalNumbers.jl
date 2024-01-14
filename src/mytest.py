import cmath
import math
import sys

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from pfapack import pfaffian as pf

def test_pfaffian():
    A = np.array([
        [0.3+0.2j, 0.6+0.8j, 0.1+0.7j, 0.9+0.4j],
        [0.4+0.2j, 0.2+0.3j, 0.5+0.1j, 0.7+0.6j],
        [0.9+0.2j, 0.1+0.3j, 0.4+0.6j, 0.8+0.2j],
        [0.6+0.6j, 0.7+0.1j, 0.3+0.8j, 0.1+0.5j]])
    # A1 = A0.transpose()
    # A1
    A = A - A.transpose()

    # Check if matrix is square
    assert A.shape[0] == A.shape[1] > 0
    # Check if it's skew-symmetric
    assert abs((A + A.T).max()) < 1e-14

    n = A.shape[0]

    # Check if we have a complex data type
    if np.issubdtype(A.dtype, np.complexfloating):
        householder = pf.householder_complex
    elif not np.issubdtype(A.dtype, np.number):
        raise TypeError("pfaffian() can only work on numeric input")
    else:
        householder = pf.householder_real

    A = np.asarray(A)  # the slice views work only properly for arrays

    pfaffian_val = 1.0

    for i in range(A.shape[0] - 2):
        # Find a Householder vector to eliminate the i-th column
        v, tau, alpha = householder(A[i + 1 :, i])
        # print(v, tau, alpha, file=sys.stderr)
        A[i + 1, i] = alpha
        A[i, i + 1] = -alpha
        A[i + 2 :, i] = 0
        A[i, i + 2 :] = 0
        # print(A, file=sys.stderr)

        # Update the matrix block A(i+1:N,i+1:N)
        w = tau * np.dot(A[i + 1 :, i + 1 :], v.conj())
        A[i + 1 :, i + 1 :] = A[i + 1 :, i + 1 :] + np.outer(v, w) - np.outer(w, v)
        # print(A, file=sys.stderr)
        # print(w, file=sys.stderr)
        # print(np.outer(v, w), file=sys.stderr)
        # print(np.outer(w, v), file=sys.stderr)
        # print(np.outer(v, w) - np.outer(w, v), file=sys.stderr)

        if tau != 0:
            pfaffian_val *= 1 - tau
        if i % 2 == 0:
            pfaffian_val *= -alpha

        print(A[n - 2, n - 1], file=sys.stderr)

    pfaffian_val *= A[n - 2, n - 1]
    print(A[n - 2, n - 1], file=sys.stderr)
    print(A, file=sys.stderr)

    return pfaffian_val, A