import numpy as np
import pytest
from numpy import array
from sol1 import ( in_lattice, simple_rounding,
    orth_proj, Gram_Schmidt_orth, nearest_plane,
    compare_norm_distrib
)
from test_bases import B2, B4, B24

def iarray(x):
	return array(x, dtype=int)

def in_span(A, v):
    x = np.linalg.solve(A.T, v)
    return np.allclose(x @ A, v)

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
def test_in_lattice(B):
    n, _ = B.shape

    # Test known lattice points
    for _ in range(20):
        x = np.random.randint(-50, 50, size=n)
        v = x @ B
        assert in_lattice(B, v), f"Known lattice point {v} not recognized"

    # Test non-lattice points
    for _ in range(20):
        x = np.random.randint(-50, 50, size=n)
        v = x @ B + np.random.uniform(0.1, 0.9, size=B.shape[1])  # Add small non-integer shift
        assert not in_lattice(B, v), f"Non-lattice point {v} incorrectly accepted"

    # Test borderline numerical tolerance
    for _ in range(20):
        x = np.random.randint(-50, 50, size=n)
        v = x @ B + 1e-9 * np.random.randn(B.shape[1])  # Add tiny numerical noise
        assert in_lattice(B, v), f"Point {v} within tolerance was not accepted"

@pytest.mark.parametrize("B", bases)
def test_simple_rounding(B):
    n, _ = B.shape

    for _ in range(20):
        x = np.random.uniform(-50, 50, size=n)
        xr = np.round(x)

        if not np.allclose(x, xr, atol=0.499): # Avoid numerical borderline cases
            continue
        
        
        t = x @ B
        res = simple_rounding(B, t)

        assert res is not None, "Result is NaN value."
        assert in_lattice(B, res), "Result vector is not in the lattice."
        assert np.allclose(res, xr @ B)


@pytest.mark.parametrize("n", [2, 3, 5, 10, 50])
def test_orth_proj(n):
    for _ in range(20):
        x = np.random.uniform(-50, 50, size=n)
        y = np.random.uniform(-50, 50, size=n)

        z = orth_proj(x, y)

        assert np.allclose(z @ y, 0) # Check if z is orthogonal to y
        assert np.allclose(z @ z, z @ x) # Check inner products 

@pytest.mark.parametrize("B", bases)
def test_gram_schmidt(B):
    Bs = Gram_Schmidt_orth(B)
    n, d = B.shape

    assert Bs.shape == B.shape, f"Expected shape {B.shape}, got {Bs.shape}"
    assert np.allclose(B[0], Bs[0])

    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                assert np.allclose(Bs[i] @ Bs[j], 0), f"Bs{i} and Bs{j} are not orthogonal."

    for i in range(1, n):
        for j in range(1, i):
            assert np.allclose(Bs[i] @ B[j], 0), f"B{i} and Bs{j} are not orthogonal."
    
    for b in B:
        assert in_span(Bs, b), "Spans of first i basis vectors are not equal."

    for bs in Bs:
        assert in_span(B, bs), "Spans of first i basis vectors are not equal."
        

# @pytest.mark.parametrize("B", bases)
# def test_gram_schmidt(B):
#     Bs = Gram_Schmidt_orth(B)
#     n, _ = B.shape

#     # Determinants should match (or very close)
#     assert np.allclose(np.linalg.det(B), np.linalg.det(Bs)), "Determinants mismatch"

#     D = Bs.dot(Bs.transpose())
#     T = B.dot(Bs.transpose())

#     assert np.allclose(np.diagonal(D), np.diagonal(T)), "Diagonal mismatch"

#     for i in range(n):
#         for j in range(i):
#             assert abs(D[j, i]) < 1e-8
#             assert abs(D[i, j]) < 1e-8
#             assert abs(T[j, i]) < 1e-8

@pytest.mark.parametrize("B", bases)
def test_nearest_plane(B):
    n, _ = B.shape
    Bs = Gram_Schmidt_orth(B)

    for _ in range(20):
        x = np.random.uniform(-50, 50, size=n)
        t = x @ B
        v = nearest_plane(B, Bs, np.copy(t))

        assert in_lattice(B, v), "Result vector is not in the lattice."

        e = t - v
        for j in range(n):
            shadow = (e @ Bs[j]) / (Bs[j] @ Bs[j])
            assert shadow < 0.501 or shadow > 0.501

# @pytest.mark.parametrize("B", bases)
# def test_nearest_plane(B):
#     n, _ = B.shape
#     Bs = Gram_Schmidt_orth(B)

#     for _ in range(20):
#         x = 50 * (np.random.rand(n) - np.random.rand(n))
#         t = x.dot(B)
#         v = nearest_plane(B, Bs, np.copy(t))

#         assert in_lattice(B, v), "Returned point is not in the lattice"

#         for bs in Bs:
#             lhs = abs((t - v).dot(bs))
#             rhs = 0.501 * bs.dot(bs)
#             assert lhs <= rhs, "Point not nearest in lattice direction"

plot_bases = [B2, B4, B24]

@pytest.mark.parametrize("B", plot_bases)
def test_ex4_compare_norm_distrb(B):
    compare_norm_distrib(B, 50000)