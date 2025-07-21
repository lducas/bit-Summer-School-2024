import numpy as np
import matplotlib.pyplot as plt

############
# Exercise 0 — Arrays
############

def gen_zeros(n, k):
    """Create a zero matrix of size n x k."""
    return np.zeros((n, k))


def gen_identity(n):
    """Create the identity matrix of size n x n."""
    return np.eye(n)


def copy(X):
    """Return a deep copy of the input matrix."""
    return np.array(X, copy=True)


def transpose(X):
    """Return the transpose of the input matrix."""
    return np.transpose(X)

def inverse(A):
    """Return the iverse of the input matrix."""

    assert A.shape[0] == A.shape[1], "Matrix must be square"
    
    B = np.linalg.inv(A)
    return B


def sum_l_rows(n, k, a, l):
    """
    Generate a random matrix of shape (n, k) with entries in [0, a),
    then return the sum of the first l rows.
    """
    assert l <= n, "Parameter l must be <= n."
    A = gen_rand(n, k, a)
    return np.sum(A[:l], axis=0, keepdims=True)


def gen_struct_mat(n, k):
    """
    Create a matrix where each row follows the pattern:

    [n-1      n+1      n+3    ...  n+2k-2
     ...                        
     3        5       7     ... 2k+1
     2        4       6     ... 2k
     1        3       5     ... 2k-1]

    So values increase by 2 across columns and decrease in base as row index increases.
    """
    A = np.zeros((n, k), dtype=int)
    for i in reversed(range(n)):
        for j in range(k):
            A[i][j] = n - i + 2 * j
    return A

############
# Exercise 1 — Random generation
############

def gen_rand(n, k, q):
    """Generate a random matrix with values in [0, q)."""
    return q * np.random.rand(n, k)


def gen_rand_centered(n, k, q):
    """Generate a random matrix with values in [-q, q)."""
    return q * (np.random.rand(n, k) - 0.5) * 2

############
# Exercise 2 — Linear Algebra
############

def inner_product(n, k, q):
    """
    Create a random vector of shape (1, n) and matrix of shape (n, k),
    both with values in [-q, q), and return their matrix product (1 x k).
    """
    x = gen_rand_centered(1, n, q)
    A = gen_rand_centered(n, k, q)
    return x @ A


def solv_lineq(n, q):
    """
    Solve a linear system Ax = b, where A is a random n x n matrix
    and b is a random n x 1 vector with entries in [-q, q).
    """
    A = gen_rand_centered(n, n, q)
    b = gen_rand_centered(n, 1, q)
    return np.linalg.solve(A, b)


def get_rref(A):
    """
    Return a matrix B such that B @ A = RREF(A), using only NumPy.
    """
    A = A.astype(float)  # Ensure float division
    n, k = A.shape

    B = np.eye(n)        # Initialize transformation matrix

    lead = 0
    for r in range(n):
        if lead >= k:
            break
        i = r
        while abs(A[i, lead]) < 1e-12:
            i += 1
            if i == n:
                i = r
                lead += 1
                if k == lead:
                    break
        if lead >= k:
            break

        # Swap rows in both A and B
        A[[r, i]] = A[[i, r]]
        B[[r, i]] = B[[i, r]]

        # Normalize pivot row
        lv = A[r, lead]
        A[r] = A[r] / lv
        B[r] = B[r] / lv

        # Eliminate column in other rows
        for i in range(n):
            if i != r:
                lv = A[i, lead]
                A[i] = A[i] - lv * A[r]
                B[i] = B[i] - lv * B[r]

        lead += 1

    print(A)
    return B, A  # A is now in RREF form

