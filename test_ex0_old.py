import numpy as np
import pytest
from itertools import product
from numpy import array
from sol0 import (
    gen_zeros, gen_identity, copy, transpose, inverse, sum_l_rows,
    gen_struct_mat, gen_rand, inner_product, solv_lineq, get_rref
)

def is_rref(A, tol=1e-8):
    """
    Check if matrix A is in reduced row echelon form (RREF).
    - Leading 1s in each non-zero row.
    - Each leading 1 is the only non-zero entry in its column.
    - Leading 1s move strictly to the right as you go down the rows.
    - Zero rows are at the bottom.
    """
    m, n = A.shape
    lead_col = -1

    for i in range(m):
        row = A[i]
        nonzero = np.where(np.abs(row) > tol)[0]
        if len(nonzero) == 0:
            continue  # Zero row — allowed at bottom
        lead = nonzero[0]

        if lead <= lead_col:
            return False  # Pivot not strictly to the right
        lead_col = lead

        if not np.isclose(row[lead], 1.0, atol=tol):
            return False  # Leading entry must be 1

        for j in range(m):
            if j != i and np.abs(A[j, lead]) > tol:
                return False  # Non-pivot rows must have 0 in pivot col

    return True


# ========== Predefined Values ==========

n_values = [3, 5, 10]
k_values = [3, 4, 6]
q_values = [2, 3, 5]
a_values = [0.3, 0.5]
b_values = [0.8, 1]
l_values = [2, 3]

# ========== Generate All Parameter Combinations ==========

test_cases_nk = list(product(n_values, k_values))
test_cases_nkab = list(product(n_values, k_values, a_values, b_values))
test_cases_nab = list(product(n_values, a_values, b_values))
test_cases_nkabl = list(product(n_values, k_values, a_values, b_values, l_values))

# ========== Random Matrices ==========

np.random.seed(42)  # For reproducibility

def gen_random_matrix(n, k, low=-10.0, high=10.0):
    return np.random.uniform(low, high, size=(n, k))

# Random float matrices for tests
matrices_copy = [gen_random_matrix(n, k) for n, k in test_cases_nk[:3]]
matrices_transpose = [gen_random_matrix(n, k) for n, k in test_cases_nk[3:6]]

# For inverse tests: square float matrices, ensure full rank
matrices_inverse = []
for n in n_values:
    if n > 5:
        continue  # keep inverse matrices small
    A = gen_random_matrix(n, n)
    if np.linalg.matrix_rank(A) == n:
        matrices_inverse.append(A)


# For RREF tests: rectangular float matrices (n ≤ k), ensure nonzero rank
rref_test_matrices = []
for n in n_values:
    for k in k_values:
        if n > k:
            continue  # Only use wide or square matrices
        A = gen_random_matrix(n, k)
        if np.linalg.matrix_rank(A) > 0:
            rref_test_matrices.append(A)

# ========== Tests ==========

@pytest.mark.parametrize("n,k", test_cases_nk)
def test_gen_zeros(n, k):
    Z = gen_zeros(n, k)
    assert Z.shape == (n, k)
    assert np.all(Z == 0)


@pytest.mark.parametrize("n", n_values)
def test_gen_identity(n):
    I = gen_identity(n)
    assert I.shape == (n, n)
    assert np.allclose(I @ I, I)
    assert np.allclose(I.T, I)


@pytest.mark.parametrize("X", matrices_copy)
def test_copy(X):
    Y = copy(X)
    assert np.allclose(Y, X)
    assert Y is not X


@pytest.mark.parametrize("X", matrices_transpose)
def test_transpose(X):
    T = transpose(X)
    assert np.allclose(T, X.T)
    assert T.shape == X.T.shape


@pytest.mark.parametrize("A", matrices_inverse)
def test_inverse(A):
    A_inv = inverse(A)
    I = np.eye(A.shape[0])
    assert np.allclose(A @ A_inv, I, atol=1e-5)
    assert np.allclose(A_inv @ A, I, atol=1e-5)

@pytest.mark.parametrize("n,k,a,b", test_cases_nkab)
def test_gen_rand(n, k, a, b):
    R = gen_rand(n, k, a, b)
    assert R.shape == (n, k)
    assert np.all(R >= a)
    assert np.all(R < b)


@pytest.mark.parametrize("n,k,a,b,l", test_cases_nkabl)
def test_sum_l_rows(n, k, a, b, l):
    result = sum_l_rows(n, k, a, b, l)
    assert result.shape == (1, k)
    assert np.all(result > a)
    assert np.all(result < l * b)


@pytest.mark.parametrize("n,k", test_cases_nk)
def test_gen_struct_mat(n, k):
    A = gen_struct_mat(n, k)
    assert A.shape == (n, k)
    for i in range(n):
        for j in range(k - 1):
            assert np.isclose(A[i, j+1] - A[i, j], 2)
            if i < n - 1:
                assert np.isclose(A[i, j] - A[i + 1, j], 1)


@pytest.mark.parametrize("n,k,a,b", test_cases_nkab)
def test_inner_product(n, k, a, b):
    result = inner_product(n, k, a, b)
    assert result.shape == (1, k)
    assert isinstance(result, np.ndarray)


@pytest.mark.parametrize("n,a,b", test_cases_nab)
def test_solv_lineq(n, a, b):
    x = solv_lineq(n, a, b)
    assert x.shape == (n, 1)
    assert isinstance(x, np.ndarray)


@pytest.mark.parametrize("A", rref_test_matrices)
def test_get_rref(A):
    A_orig = A.copy()
    B, A_rref = get_rref(A.copy())

    # 1. Check transformation: B @ A = RREF(A)
    assert np.allclose(B @ A_orig, A_rref, atol=1e-6), "B @ A does not equal RREF(A)"

    # 2. Check that A_rref is actually in RREF
    assert is_rref(A_rref), "Matrix is not in RREF form"