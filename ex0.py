import numpy as np
import matplotlib.pyplot as plt

# Replace 'pass' with your implementation in each exercise.

############
# Exercise 0 — Arrays
############

def gen_zeros(n, k):
    """Create a zero matrix of size n x k."""
    pass

def gen_identity(n):
    """Create the identity matrix of size n x n."""
    pass


def copy(X):
    """Return a deep copy of the input matrix."""
    pass


def transpose(X):
    """Return the transpose of the input matrix."""
    pass

def inverse(A):
    """Return the iverse of the input matrix."""

    assert A.shape[0] == A.shape[1], "Matrix must be square"
    
    pass

def sum_l_rows(n, k, a, l):
    """
    Generate a random matrix of shape (n, k) with entries in [0, a),
    then return the sum of the first l rows.
    """
    assert l <= n, "Parameter l must be <= n."
    
    pass


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
    pass

############
# Exercise 1 — Random generation
############

def gen_rand(n, k, q):
    """Generate a random matrix with values in [0, q)."""
    pass

def gen_rand_centered(n, k, q):
    """Generate a random matrix with values in [-q, q)."""
    pass

############
# Exercise 2 — Linear Algebra
############

def inner_product(n, k, q):
    """
    Create a random vector of shape (1, n) and matrix of shape (n, k),
    both with values in [-q, q), and return their matrix product (1 x k).
    """
    pass


def solv_lineq(n, q):
    """
    Solve a linear system Ax = b, where A is a random n x n matrix
    and b is a random n x 1 vector with entries in [-q, q).
    """
    pass


def get_rref(A):
    """
    Return a matrix B such that B @ A = RREF(A), using only NumPy.
    """
    A = A.astype(float)  # Ensure float division
    n, k = A.shape

    B = np.eye(n)        # Initialize transformation matrix

    pass

