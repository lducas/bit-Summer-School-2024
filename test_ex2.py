import matplotlib.pyplot
import numpy as np
from numpy import array
from math import sqrt
import pytest
import matplotlib

# Set random seed for reproducibility
np.random.seed(42)
np.random.seed(42)

from sol2 import lagrange_reduce, size_reduce, LLL, anim_LLL
from sol1 import Gram_Schmidt_orth, in_lattice

# -------------------- Test Data --------------------

def iarray(x):
    return array(x, dtype=int)

bases = [
    array([[1]]),
    array([[2]]),
    array([[200]]),
    iarray([[4, 0], [3, 1]]),
    iarray([[50, 0], [0, 1]]),
    iarray([[50, 0], [40, 1]]),
    iarray([[50, -30, 14], [0, 200, -4], [350, 600, -12]]),
    iarray([[50, 33, -30, 14], [200, 40, -45, -4], [1, 35, 0, -1200], [-15, 3, 8, -7]])
]

def gen_same_lat(B, C):
    return all(in_lattice(C, b) for b in B) and all(in_lattice(B, c) for c in C)

# -------------------- Exercise 1 --------------------

@pytest.mark.parametrize("B", [
    iarray([[999, 0], [np.random.randint(0, 999), 1]])
    for _ in range(30)
])
def test_lagrange_reduce(B):
    B_orig = np.copy(B)
    U = lagrange_reduce(B)

    assert U is not None, "Returned matrix U is None"
    assert isinstance(U[0, 0], np.int64), "U is not integer matrix"
    assert np.round(abs(np.linalg.det(U))) == 1, "U is not unimodular"
    assert np.all(B == U @ B_orig), "B_out != U * B_in"
    assert np.linalg.norm(B[1]) >= np.linalg.norm(B[0]), "Not Lagrange-reduced"
    
    mu = abs(B[0].dot(B[1]) / B[0].dot(B[0]))
    assert mu <= 0.500001, f"mu = {mu} > 0.5 — Not Lagrange-reduced"

# -------------------- Exercise 2 --------------------

@pytest.mark.parametrize("B", bases)
def test_size_reduce(B):
    
    Bs = Gram_Schmidt_orth(B)

    B_reduced = np.copy(B)
    size_reduce(B_reduced, Bs)

    assert gen_same_lat(B_reduced, B), "Output basis does not generate the same lattice"

    n, _ = B.shape
    for i in range(n):
        for j in range(i):
            mu = abs(B_reduced[i] @ Bs[j]) / (Bs[j] @ Bs[j])
            assert mu <= 0.500001, f"mu = {mu} > 0.5 — Not size-reduced"

# -------------------- Exercise 3 --------------------

@pytest.mark.parametrize("n", list(range(2, 20)))  # reduce upper limit for speed
def test_lll(n):
    q = 999999
    B_ = np.identity(n, dtype=int)
    B_[0, 0] = q
    for i in range(1, n):
        B_[i, 0] = np.random.randint(0, q)

    B = np.copy(B_)
    epsilon = 0.01
    gamma_2 = sqrt(4 / 3)

    try:
        steps = list(LLL(B, epsilon=epsilon, animate=False))
    except RuntimeError:
        pytest.fail("LLL did not converge")

    assert gen_same_lat(B, B_), "LLL output does not generate same lattice"

    Bs = Gram_Schmidt_orth(B)

    for i in range(n - 1):
        assert np.linalg.norm(Bs[i]) <= (gamma_2 + epsilon) * np.linalg.norm(Bs[i + 1]), "Not weak-LLL reduced"

    for i in range(n):
        for j in range(i):
            mu = abs(B[i] @ Bs[j]) / (Bs[j] @ Bs[j])
            assert mu <= 0.500001, f"mu = {mu} > 0.5 — Not size-reduced"

# -------------------- Exercise 4 (Animation) --------------------

def test_anim_lll_runs():
    try:
        # This should not raise an error even in headless environments
        ani = anim_LLL(n=20, q=999999)
        matplotlib.pyplot.show()

    except Exception as e:
        pytest.fail(f"anim_LLL raised an unexpected exception: {e}")
