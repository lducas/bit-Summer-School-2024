import numpy as np
import matplotlib.pyplot as plt

# --- Exercises using loops (patterns) ---

def gen_checkerboard(n, k):
    """Return an n x k checkerboard of 0s and 1s."""
    A = np.zeros((n, k), dtype=int)
    for i in range(n):
        for j in range(k):
            A[i, j] = (i + j) % 2
    return A


def gen_triangle_mat(n):
    """Return an n x n lower-triangular matrix of 1s."""
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1):
            A[i, j] = 1
    return A

# --- Exercises using slicing / indexing ---

def reverse_rows(A):
    """Return a matrix with rows reversed."""
    return A[::-1, :]


def extract_diag(A):
    """Return the main diagonal of a square matrix as a vector."""
    return np.diag(A)

# --- Exercises using random generation ---

def gen_rand_int(n, k, low, high):
    """Return an n x k matrix of random integers in [low, high)."""
    return np.random.randint(low, high, size=(n, k))


def sum_first_rows(n, k, l):
    """Return the sum of the first l rows of a random matrix."""
    A = np.random.rand(n, k)
    row_sum = np.sum(A[:l, :], axis=0)
    print("Matrix:\n", A)
    print("Sum of first", l, "rows:", row_sum)
    return row_sum

# --- Exercises using algebra (dot product, trace, etc.) ---
def mul_by_vector(n, k):
    """Return the product of an n x k matrix and a k x 1 vector of ones."""
    A = np.arange(1, n*k+1).reshape(n, k)
    v = np.ones((k, 1))
    return A @ v


def trace_of_square(n):
    """Generate an n x n matrix and return its trace (sum of diagonal)."""
    A = np.arange(1, n*n+1).reshape(n, n)
    return np.trace(A)

# --- More interesting exercises ---
def demo_identity_effect(A):
    """Show that multiplying by identity leaves a matrix unchanged."""
    I = np.identity(A.shape[0])
    result = A @ I
    print("Original matrix:\n", A)
    print("A * I = \n", result)
    return result


def check_orthogonality(u, v):
    """Check if two vectors are orthogonal using inner product."""
    ip = np.dot(u, v)
    print(f"Inner product = {ip}")
    if np.isclose(ip, 0):
        print("Vectors are orthogonal!")
    else:
        print("Vectors are not orthogonal.")
    return ip


def demo_transpose_properties(A):
    """Show that transpose of transpose is the original, and test for symmetry."""
    T = A.T
    TT = T.T
    print("Matrix A:\n", A)
    print("Transpose A^T:\n", T)
    print("Double transpose (A^T)^T:\n", TT)
    if np.allclose(A, T):
        print("Matrix is symmetric!")
    return T


def demo_inverse(A):
    """Show that multiplying a matrix by its inverse yields identity."""
    try:
        invA = np.linalg.inv(A)
        I_check = A @ invA
        print("Matrix A:\n", A)
        print("Inverse A^-1:\n", invA)
        print("A * A^-1 = \n", I_check)
        return invA
    except np.linalg.LinAlgError:
        print("Matrix is not invertible.")
        return None

# --- Plotting helpers ---
def plot_matrix(A, title="Matrix", cmap="viridis"):
    """Plot a matrix using imshow with a colorbar."""
    plt.imshow(A, cmap=cmap, interpolation='nearest')
    plt.colorbar()
    plt.title(title)
    plt.show()


def plot_checkerboard(A):
    """Plot a checkerboard as a chessboard (black & white)."""
    plt.imshow(A, cmap="gray", interpolation='nearest')
    plt.title("Checkerboard / Chessboard")
    plt.show()

# --- Worksheet Example ---
if __name__ == "__main__":
    print("\n1. Checkerboard:")
    cb = gen_checkerboard(8, 8)
    print(cb)
    plot_checkerboard(cb)

    print("\n2. Lower-triangular:")
    tri = gen_triangle_mat(5)
    print(tri)
    plot_matrix(tri, title="Lower-triangular")

    print("\n3. Reverse rows:")
    mat3 = np.arange(1, 10).reshape(3, 3)
    rev = reverse_rows(mat3)
    print(rev)
    plot_matrix(rev, title="Rows reversed")

    print("\n4. Extract diagonal:")
    mat4 = np.arange(1, 17).reshape(4, 4)
    diag = extract_diag(mat4)
    print(diag)

    print("\n5. Random integers:")
    rand_ints = gen_rand_int(3, 4, 0, 10)
    print(rand_ints)
    plot_matrix(rand_ints, title="Random integers")

    print("\n6. Sum of first rows:")
    sum_first_rows(5, 3, 2)

    print("\n7. Multiply by vector of ones:")
    print(mul_by_vector(3, 4))

    print("\n8. Trace of square:")
    print(trace_of_square(4))

    print("\n9. Identity effect:")
    demo_identity_effect(np.array([[2, 1], [0, 3]]))

    print("\n10. Orthogonality check:")
    check_orthogonality(np.array([1, 0, 0]), np.array([0, 1, 0]))

    print("\n11. Transpose properties:")
    demo_transpose_properties(np.array([[1, 2], [2, 1]]))

    print("\n12. Inverse demo:")
    demo_inverse(np.array([[4, 7], [2, 6]]))
