import numpy as np
import pytest
from numpy import array
from ex1 import (
    iarray, in_lattice, simple_rounding,
    orth_proj, Gram_Schmidt_orth, nearest_plane,
    compare_norm_distrib
)
from test_bases import B2, B4, B24

# All test bases
bases = [
iarray([[1],]),
iarray([[2],]),
iarray([[200],]),
iarray([[1,0],[0,1]]),
iarray([[50,0],[0,1]]),
iarray([[50,0],[20,1]]),
iarray([[50, -30, 14], [0, 20, -4], [35, 0, -12]]),
iarray([[50, 33, -30, 14], [0, 20, 4, -4], [1, 35, 0, -12], [-15, 3, 8, -7]])
]

@pytest.mark.parametrize("B", bases)
def test_ex1_simple_rounding(B):
    n, _ = B.shape
    for _ in range(20):
        x = 50 * (np.random.rand(n) - np.random.rand(n))
        xr = np.round(x)

        # Skip borderline cases due to numerical tolerance
        if max(abs(x - xr)) > 0.499:
            continue

        t = x.dot(B)
        res = simple_rounding(B, t)

        assert res is not None, "Returned None"
        assert np.allclose(res, xr.dot(B)), "Rounded point doesn't match"

@pytest.mark.parametrize("n", [2, 3, 5, 10, 50])
def test_ex2a_orth_proj(n):
    for _ in range(20):
        x = 50 * (np.random.rand(n) - np.random.rand(n))
        y = 50 * (np.random.rand(n) - np.random.rand(n))

        z = orth_proj(x, y)

         # Test 1: z is orthogonal to y
        assert np.allclose(z.dot(y), 0), "z is not orthogonal to y"

        # Test 2: inner product with x equals norm squared
        assert np.allclose(z.dot(x), z.dot(z)), "z·x != z·z"

@pytest.mark.parametrize("B", bases)
def test_ex2b_gram_schmidt(B):
    Bs = Gram_Schmidt_orth(B)
    n, _ = B.shape

    # Determinants should match (or very close)
    assert np.allclose(np.linalg.det(B), np.linalg.det(Bs)), "Determinants mismatch"

    D = Bs.dot(Bs.transpose())
    T = B.dot(Bs.transpose())

    assert np.allclose(np.diagonal(D), np.diagonal(T)), "Diagonal mismatch"

    for i in range(n):
        for j in range(i):
            assert abs(D[j, i]) < 1e-8
            assert abs(D[i, j]) < 1e-8
            assert abs(T[j, i]) < 1e-8

@pytest.mark.parametrize("B", bases)
def test_ex3_nearest_plane(B):
    n, _ = B.shape
    Bs = Gram_Schmidt_orth(B)

    for _ in range(20):
        x = 50 * (np.random.rand(n) - np.random.rand(n))
        t = x.dot(B)
        v = nearest_plane(B, Bs, np.copy(t))

        assert in_lattice(B, v), "Returned point is not in the lattice"

        for bs in Bs:
            lhs = abs((t - v).dot(bs))
            rhs = 0.501 * bs.dot(bs)
            assert lhs <= rhs, "Point not nearest in lattice direction"

plot_bases = [B2, B4, B24]
@pytest.mark.parametrize("B", plot_bases)
def test_ex4_compare_norm_distrb(B):
    compare_norm_distrib(B, 50000)